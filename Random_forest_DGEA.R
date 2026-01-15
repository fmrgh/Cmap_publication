# analyze_zebrafish_expression.R
# Purpose: Analyze zebrafish gene expression using `lfc_all.csv` and `MoA_classification.csv`.
# - Handles semicolon-delimited CSVs; decimal comma in lfc_all
# - Normalizes substance names for robust joins to MoA
# - Filtering:
#     (1) Keep genes with any padj < 0.05
#     (2) Keep genes significant (padj <= 0.05) for ALL substances of at least one MoA
#     (3) Keep genes with consistent LFC sign within such an MoA (all >=0 or all <=0)
# - Summaries, volcano plots, top-variable heatmap (PNG+PDF)
# - Random Forest (stratified, robust) to pick Top-100 discriminating genes + heatmap (PNG+PDF)
#   *Only genes that are significantly regulated (padj <= 0.05) within at least one MoA and
#    with consistent LFC sign across its substances are kept for RF; MoA-wise masking
#    sets other (gene,MoA) combos to 0 before RF.*
# - Ordered RF heatmap (fixed substance order, 180° x-labels, SYMBOL on right, MoA colors Set2/fill_map)
# - Diagnostics: list unmapped substances -> outputs/substances_missing_moa.csv
# Author: (auto-generated)

options(stringsAsFactors = FALSE)
OUTPUT_DIR <- "outputs"
INPUT_LFC_FILE <- "lfc_all.csv"
INPUT_MOA_FILE <- "MoA_classification.csv"

# ----------------------------
# Helper: install/load packages
# ----------------------------
ensure_pkg <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p, repos = "https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}
ensure_pkg(c(
  "readr","dplyr","tidyr","stringr","ggplot2","pheatmap","scales",
  "purrr","tibble","randomForest","forcats","ComplexHeatmap",
  "circlize","grid","RColorBrewer"
))

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
stopifnot(file.exists(INPUT_LFC_FILE), file.exists(INPUT_MOA_FILE))

# ----------------------------
# Config & helpers
# ----------------------------
FILL_UNKNOWN_IN_PLOTS <- TRUE   # if FALSE: drop columns without MoA from heatmaps
TOP_VAR_N <- 50                 # top-N variable genes for the general heatmap
RF_N_ITER <- 50                 # RF bootstrap iters
RF_NTREE  <- 800                # trees per RF

# Normalize substance names to maximize matching between lfc_all and MoA table
norm_subst <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

# ----------------------------
# Shared MoA color helper (Set2 / fill_map-compatible)
# ----------------------------
get_moa_color_map <- function(moa_levels) {
  moa_levels <- unique(as.character(moa_levels))
  # 1) If fill_map exists (from GSEA script), reuse/extend it
  if (exists("fill_map")) {
    col_map <- fill_map
    missing_lvls <- setdiff(moa_levels, names(col_map))
    if (length(missing_lvls) == 0) {
      return(col_map[moa_levels])
    } else {
      n_miss <- length(missing_lvls)
      base_cols <- RColorBrewer::brewer.pal(
        min(8, max(3, n_miss)),
        "Set2"
      )
      extra_cols <- if (n_miss > length(base_cols)) {
        grDevices::colorRampPalette(base_cols)(n_miss)
      } else {
        base_cols[seq_len(n_miss)]
      }
      names(extra_cols) <- missing_lvls
      col_map <- c(col_map, extra_cols)
      return(col_map[moa_levels])
    }
  }
  
  # 2) Otherwise, create a fresh Set2-based palette
  n <- length(moa_levels)
  base_cols <- RColorBrewer::brewer.pal(
    min(8, max(3, n)),
    "Set2"
  )
  cols <- if (n > length(base_cols)) {
    grDevices::colorRampPalette(base_cols)(n)
  } else {
    base_cols[seq_len(n)]
  }
  names(cols) <- moa_levels
  cols
}

# ----------------------------
# Read data with correct delimiters/locales
# ----------------------------
# lfc_all.csv: semicolon delimiter, decimal comma
lfc_raw <- readr::read_delim(
  INPUT_LFC_FILE,
  delim = ";",
  locale = readr::locale(decimal_mark = ",", grouping_mark = ".", encoding = "UTF-8"),
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA")
)

# MoA_classification.csv: semicolon delimiter
moa_raw <- readr::read_delim(
  INPUT_MOA_FILE,
  delim = ";",
  locale = readr::locale(encoding = "UTF-8"),
  show_col_types = FALSE,
  progress = FALSE,
  na = c("", "NA")
)

# ----------------------------
# Normalize/clean MoA table (explicit columns from file)
# ----------------------------
if (!all(c("Substance","Mode_of_action") %in% names(moa_raw))) {
  stop("MoA_classification.csv must contain columns 'Substance' and 'Mode_of_action'.")
}
moa_tbl <- moa_raw %>%
  dplyr::select(Substance, Mode_of_action) %>%
  dplyr::rename(substance = Substance) %>%
  dplyr::mutate(
    substance = as.character(substance),
    Mode_of_action = as.character(Mode_of_action)
  ) %>%
  dplyr::mutate(substance_norm = norm_subst(substance))

# ----------------------------
# Identify ID and metric columns
# ----------------------------
id_candidates <- tolower(names(lfc_raw))
ensembl_col <- if ("ensembl" %in% id_candidates) names(lfc_raw)[which(id_candidates == "ensembl")[1]] else names(lfc_raw)[1]
symbol_col  <- if ("symbol"  %in% id_candidates) names(lfc_raw)[which(id_candidates == "symbol")[1]]  else names(lfc_raw)[2]

lfc_cols  <- names(lfc_raw)[stringr::str_detect(names(lfc_raw), "_lfc$")]
padj_cols <- names(lfc_raw)[stringr::str_detect(names(lfc_raw), "_padj$")]
if (length(lfc_cols) == 0 || length(padj_cols) == 0) stop("No *_lfc / *_padj columns detected.")

lfc_subs  <- stringr::str_replace(lfc_cols, "_lfc$", "")
padj_subs <- stringr::str_replace(padj_cols, "_padj$", "")
common_subs <- intersect(lfc_subs, padj_subs)
if (length(common_subs) == 0) stop("No matching substances for *_lfc and *_padj columns.")

# Diagnostics: unmapped substances after normalization
expr_norm_map <- tibble::tibble(substance = common_subs, substance_norm = norm_subst(common_subs))
unmapped <- dplyr::anti_join(expr_norm_map, moa_tbl %>% dplyr::select(substance_norm), by = "substance_norm")
if (nrow(unmapped) > 0) {
  readr::write_csv(unmapped, file.path(OUTPUT_DIR, "substances_missing_moa.csv"))
  message("Diagnostics: wrote 'substances_missing_moa.csv' listing substances without MoA mapping (after normalization).")
}

# ----------------------------
# FILTER STEP 1: Keep genes with any padj < 0.05
# ----------------------------
padj_mat <- as.data.frame(lfc_raw[, paste0(common_subs, "_padj"), drop = FALSE])
keep_any_signif <- apply(padj_mat, 1, function(x) any(as.numeric(x) < 0.05, na.rm = TRUE))
lfc_raw <- lfc_raw[keep_any_signif, , drop = FALSE]
message("Step 1: removed ", sum(!keep_any_signif), " genes with no significant padj across any substance.")

# ----------------------------
# FILTER STEPS 2 & 3: MoA-consistency of significance and LFC sign
# ----------------------------
# Build moa_map using normalized keys
moa_map <- dplyr::inner_join(
  tibble::tibble(substance = common_subs, substance_norm = norm_subst(common_subs)),
  moa_tbl %>% dplyr::select(substance_norm, Mode_of_action),
  by = "substance_norm"
) %>% dplyr::select(substance, Mode_of_action)

subs_by_moa <- split(moa_map$substance, moa_map$Mode_of_action)
padj_cols_by_moa <- lapply(subs_by_moa, function(subs) paste0(subs, "_padj"))
lfc_cols_by_moa  <- lapply(subs_by_moa, function(subs) paste0(subs, "_lfc"))

if (length(subs_by_moa) == 0) {
  warning("No overlapping MoA mapping; skipping MoA-consistency filters (steps 2 & 3).")
} else {
  # STEP 2: significance in ALL substances within at least one MoA
  keep_moa_consistent <- vapply(seq_len(nrow(lfc_raw)), function(i) {
    any(vapply(names(padj_cols_by_moa), function(moa_name) {
      cols <- intersect(padj_cols_by_moa[[moa_name]], names(lfc_raw))
      if (length(cols) == 0) return(FALSE)
      vals <- suppressWarnings(as.numeric(lfc_raw[i, cols, drop = TRUE]))
      if (any(is.na(vals))) return(FALSE)
      all(vals <= 0.05)
    }, logical(1)))
  }, logical(1))
  message("Step 2: removed ", sum(!keep_moa_consistent), " genes not consistently significant within any MoA group.")
  lfc_raw <- lfc_raw[keep_moa_consistent, , drop = FALSE]
  
  # STEP 3: sign consistency within such an MoA (all >=0 or all <=0)
  keep_moa_sign <- vapply(seq_len(nrow(lfc_raw)), function(i) {
    any(vapply(names(padj_cols_by_moa), function(moa_name) {
      padj_cols_moa <- intersect(padj_cols_by_moa[[moa_name]], names(lfc_raw))
      lfc_cols_moa  <- intersect(lfc_cols_by_moa[[moa_name]],  names(lfc_raw))
      if (length(padj_cols_moa) == 0 || length(lfc_cols_moa) == 0) return(FALSE)
      padj_vals <- suppressWarnings(as.numeric(lfc_raw[i, padj_cols_moa, drop = TRUE]))
      if (any(is.na(padj_vals)) || !all(padj_vals <= 0.05)) return(FALSE)
      lfc_vals <- suppressWarnings(as.numeric(lfc_raw[i, lfc_cols_moa, drop = TRUE]))
      if (any(is.na(lfc_vals))) return(FALSE)
      all(lfc_vals >= 0) || all(lfc_vals <= 0)
    }, logical(1)))
  }, logical(1))
  message("Step 3: removed ", sum(!keep_moa_sign), " genes with mixed LFC signs within MoA groups (despite significance).")
  lfc_raw <- lfc_raw[keep_moa_sign, , drop = FALSE]
}

# ----------------------------
# Tidy long format and join MoA (normalized join)
# ----------------------------
tidy_long <- lfc_raw %>%
  dplyr::select(all_of(c(ensembl_col, symbol_col, lfc_cols, padj_cols))) %>%
  dplyr::rename(ENSEMBL = all_of(ensembl_col), SYMBOL = all_of(symbol_col)) %>%
  tidyr::pivot_longer(cols = -c(ENSEMBL, SYMBOL),
                      names_to = c("substance", "metric"),
                      names_pattern = "^(.*)_(lfc|padj)$",
                      values_to = "value") %>%
  tidyr::pivot_wider(names_from = metric, values_from = value) %>%
  dplyr::mutate(substance = as.character(substance),
                substance_norm = norm_subst(substance)) %>%
  dplyr::left_join(moa_tbl %>% dplyr::select(substance_norm, Mode_of_action),
                   by = "substance_norm") %>%
  dplyr::mutate(
    padj = as.numeric(padj),
    lfc  = as.numeric(lfc),
    is_signif = !is.na(padj) & (padj < 0.05)
  )

readr::write_csv(tidy_long, file.path(OUTPUT_DIR, "tidy_long_table.csv"))

# ----------------------------
# Summaries
# ----------------------------
sig_counts_substance <- tidy_long %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::group_by(substance) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    signif = sum(padj < 0.05, na.rm = TRUE),
    up = sum(padj < 0.05 & lfc > 0, na.rm = TRUE),
    down = sum(padj < 0.05 & lfc < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(signif))
readr::write_csv(sig_counts_substance, file.path(OUTPUT_DIR, "sig_counts_per_substance.csv"))

sig_counts_moa <- tidy_long %>%
  dplyr::filter(!is.na(Mode_of_action), !is.na(padj)) %>%
  dplyr::group_by(Mode_of_action) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    signif = sum(padj < 0.05, na.rm = TRUE),
    up = sum(padj < 0.05 & lfc > 0, na.rm = TRUE),
    down = sum(padj < 0.05 & lfc < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(signif))
readr::write_csv(sig_counts_moa, file.path(OUTPUT_DIR, "sig_counts_per_moa.csv"))

# ----------------------------
# Volcano plots
# ----------------------------
volcano_data <- tidy_long %>%
  dplyr::mutate(
    padj_clipped = pmin(pmax(padj, 1e-300), 1),
    neglog10_padj = -log10(padj_clipped),
    signif_label = dplyr::case_when(
      is.na(padj) ~ "NA",
      padj < 0.05 & lfc > 0 ~ "Up (padj<0.05)",
      padj < 0.05 & lfc < 0 ~ "Down (padj<0.05)",
      TRUE ~ "NS"
    )
  )

p_volcano <- ggplot2::ggplot(volcano_data, ggplot2::aes(x = lfc, y = neglog10_padj)) +
  ggplot2::geom_point(ggplot2::aes(shape = signif_label), alpha = 0.5, size = 0.7) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggplot2::facet_wrap(~ substance, scales = "free_x") +
  ggplot2::labs(title = "Volcano plots per substance",
                x = "log2 fold change (lfc)",
                y = expression(-log[10](padj)),
                shape = "Significance") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(file.path(OUTPUT_DIR, "volcano_plots_faceted.png"), p_volcano, width = 14, height = 9, dpi = 150)

# ----------------------------
# Heatmap: top variable genes across substances (LFC) + PDF
# ----------------------------
lfc_mat <- lfc_raw %>%
  dplyr::select(all_of(c(ensembl_col, symbol_col, paste0(common_subs, "_lfc")))) %>%
  dplyr::rename(ENSEMBL = all_of(ensembl_col), SYMBOL = all_of(symbol_col))

padj_mat_full <- lfc_raw %>%
  dplyr::select(all_of(c(ensembl_col, symbol_col, paste0(common_subs, "_padj")))) %>%
  dplyr::rename(ENSEMBL = all_of(ensembl_col), SYMBOL = all_of(symbol_col))

lfc_numeric <- as.matrix(lfc_mat[, -c(1, 2)])
rownames(lfc_numeric) <- ifelse(is.na(lfc_mat$SYMBOL) | lfc_mat$SYMBOL == "", lfc_mat$ENSEMBL, lfc_mat$SYMBOL)

# robust numeric conversion that preserves dimensions & names
padj_numeric <- padj_mat_full[, -c(1, 2)]
padj_numeric <- as.data.frame(lapply(padj_numeric, function(x) as.numeric(x)))
padj_numeric <- as.matrix(padj_numeric)
rownames(padj_numeric) <- rownames(lfc_numeric)

# Gene variance is computed on the original LFC values
gene_var <- apply(lfc_numeric, 1, function(x) stats::var(x, na.rm = TRUE))
keep_idx <- order(gene_var, decreasing = TRUE)[seq_len(min(TOP_VAR_N, length(gene_var)))]
lfc_top  <- lfc_numeric[keep_idx, , drop = FALSE]
padj_top <- padj_numeric[keep_idx, , drop = FALSE]

# Set LFC to 0 for non-significant cells (padj > 0.05 OR NA)
lfc_top[is.na(padj_top) | padj_top > 0.05] <- 0

# Build robust column annotation using normalized join
anno_col <- tibble::tibble(substance = colnames(lfc_top),
                           substance_norm = norm_subst(colnames(lfc_top))) %>%
  dplyr::left_join(moa_tbl %>% dplyr::select(substance_norm, Mode_of_action), by = "substance_norm") %>%
  dplyr::select(substance, Mode_of_action) %>%
  tibble::column_to_rownames("substance")

if (!FILL_UNKNOWN_IN_PLOTS) {
  keep_cols <- rownames(anno_col)[!is.na(anno_col$Mode_of_action) & anno_col$Mode_of_action != ""]
  lfc_top <- lfc_top[, keep_cols, drop = FALSE]
  anno_col <- anno_col[keep_cols, , drop = FALSE]
}
if (FILL_UNKNOWN_IN_PLOTS) {
  anno_col$Mode_of_action[is.na(anno_col$Mode_of_action) | anno_col$Mode_of_action == ""] <- "Unknown"
}

if (nrow(anno_col) > 0) {
  moa_levels <- unique(as.character(anno_col$Mode_of_action))
  anno_col$Mode_of_action <- factor(anno_col$Mode_of_action, levels = moa_levels)
  palette_moa <- get_moa_color_map(moa_levels)
  ann_colors <- list(Mode_of_action = palette_moa)
} else {
  ann_colors <- NULL
}

png(filename = file.path(OUTPUT_DIR, "heatmap_top_variable_genes.png"), width = 1400, height = 1000, res = 150)
pheatmap::pheatmap(lfc_top,
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   clustering_method = "complete",
                   scale = "row",
                   annotation_col = anno_col,
                   annotation_colors = ann_colors,
                   main = "Top variable genes (by LFC variance; non-significant LFC set to 0)")
dev.off()

pdf(file = file.path(OUTPUT_DIR, "heatmap_top_variable_genes.pdf"), width = 11, height = 8.5)
pheatmap::pheatmap(lfc_top,
                   clustering_distance_rows = "euclidean",
                   clustering_distance_cols = "euclidean",
                   clustering_method = "complete",
                   scale = "row",
                   annotation_col = anno_col,
                   annotation_colors = ann_colors,
                   main = "Top variable genes (by LFC variance; non-significant LFC set to 0)")
dev.off()

# ----------------------------
# RF gene eligibility (gene-MoA specific):
# For each (gene, MoA): all padj <= 0.05 AND all LFC >=0 OR all <=0  -> pass
# We'll use:
#  - keep_genes_rf: genes passing in at least one MoA (eligible for RF)
#  - MoA-wise masking: for MoAs where pass == FALSE, LFCs are set to 0 before RF.
# ----------------------------
rf_gene_elig <- tidy_long %>%
  dplyr::filter(!is.na(Mode_of_action), !is.na(padj), !is.na(lfc)) %>%
  dplyr::group_by(ENSEMBL, Mode_of_action) %>%
  dplyr::summarise(
    all_sig  = all(padj <= 0.05),
    all_up   = all(lfc >= 0),
    all_down = all(lfc <= 0),
    .groups  = "drop"
  ) %>%
  dplyr::mutate(pass = all_sig & (all_up | all_down))

keep_genes_rf <- rf_gene_elig %>%
  dplyr::group_by(ENSEMBL) %>%
  dplyr::summarise(pass_any = any(pass), .groups = "drop") %>%
  dplyr::filter(pass_any) %>%
  dplyr::pull(ENSEMBL)

if (length(keep_genes_rf) == 0) {
  warning("No genes fulfill RF significance+sign criteria (padj <= 0.05 within a MoA, same LFC sign). Random Forest will be skipped.")
}

# ----------------------------
# Random-Forest-based Top-100 discriminative genes (robust sampling)
# ----------------------------
message("Running Random Forest to identify MoA-discriminating gene signatures (Top 100)...")

if (length(keep_genes_rf) > 0) {
  
  # MoA-wise masking: join rf_gene_elig and set lfc_masked = lfc if pass, else 0
  lfc_for_rf <- tidy_long %>%
    dplyr::filter(!is.na(Mode_of_action),
                  ENSEMBL %in% keep_genes_rf) %>%
    dplyr::left_join(
      rf_gene_elig %>%
        dplyr::select(ENSEMBL, Mode_of_action, pass),
      by = c("ENSEMBL", "Mode_of_action")
    ) %>%
    dplyr::mutate(
      GENE_ID = ENSEMBL,
      GENE_LABEL = dplyr::if_else(is.na(SYMBOL) | SYMBOL == "", ENSEMBL, SYMBOL),
      lfc_masked = dplyr::if_else(!is.na(pass) & pass, lfc, 0)
    ) %>%
    dplyr::group_by(substance, Mode_of_action, GENE_ID) %>%
    dplyr::summarise(
      lfc = mean(lfc_masked, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = GENE_ID, values_from = lfc,
      values_fn = list(lfc = function(x) mean(as.numeric(x), na.rm = TRUE)),
      values_fill = list(lfc = NA_real_)
    )
  
  if (!("Mode_of_action" %in% names(lfc_for_rf))) {
    warning("No Mode_of_action labels found after RF gene filtering; skipping Random Forest step.")
  } else {
    labels <- as.factor(lfc_for_rf$Mode_of_action)
    feature_df <- lfc_for_rf %>% dplyr::select(-substance, -Mode_of_action)
    
    # Remove all-NA or zero-variance features
    nzv <- vapply(feature_df, function(col) {
      v <- as.numeric(col); v <- v[!is.na(v)]
      if (length(v) <= 1) return(FALSE)
      stats::var(v) > 0
    }, logical(1))
    feature_df <- feature_df[, nzv, drop = FALSE]
    
    # Impute NAs by gene-wise median
    for (j in seq_len(ncol(feature_df))) {
      col <- as.numeric(feature_df[[j]])
      med <- stats::median(col, na.rm = TRUE)
      col[is.na(col)] <- med
      feature_df[[j]] <- col
    }
    
    # Drop classes with <2 samples
    class_tab <- table(labels)
    keep_classes <- names(class_tab)[class_tab >= 2]
    if (length(keep_classes) < 2) {
      warning("Random Forest skipped: fewer than two MoA classes have >=2 samples.")
    } else {
      keep_idx_rows <- labels %in% keep_classes
      feature_df <- feature_df[keep_idx_rows, , drop = FALSE]
      labels <- droplevels(labels[keep_idx_rows])
      
      # Balanced sampling vector
      class_tab <- table(labels)
      sampsize_vec <- rep(min(class_tab), length(class_tab))
      names(sampsize_vec) <- names(class_tab)
      
      set.seed(123)
      n_iter <- RF_N_ITER
      ntree <- RF_NTREE
      imp_mat <- matrix(0, nrow = ncol(feature_df), ncol = n_iter,
                        dimnames = list(colnames(feature_df), NULL))
      
      for (i in seq_len(n_iter)) {
        rf_fit <- randomForest::randomForest(
          x = feature_df,
          y = labels,
          ntree = ntree,
          importance = TRUE,
          replace = TRUE,
          strata = labels,
          sampsize = sampsize_vec
        )
        imp <- randomForest::importance(rf_fit, type = 2) # MeanDecreaseGini
        imp_vec <- rep(0, ncol(feature_df)); names(imp_vec) <- colnames(feature_df)
        common <- intersect(names(imp_vec), rownames(imp))
        imp_vec[common] <- imp[common, "MeanDecreaseGini"]
        imp_mat[, i] <- imp_vec
      }
      
      imp_mean <- rowMeans(imp_mat, na.rm = TRUE)
      imp_sd   <- apply(imp_mat, 1, sd, na.rm = TRUE)
      imp_tbl <- tibble::tibble(GENE_ID = names(imp_mean),
                                rf_importance_mean = as.numeric(imp_mean),
                                rf_importance_sd   = as.numeric(imp_sd)) %>%
        dplyr::arrange(dplyr::desc(rf_importance_mean))
      readr::write_csv(imp_tbl, file.path(OUTPUT_DIR, "RF_feature_importance_all_genes.csv"))
      
      TOP_K <- min(100, nrow(imp_tbl))
      top100 <- imp_tbl %>% dplyr::slice_head(n = TOP_K)
      readr::write_csv(top100, file.path(OUTPUT_DIR, "RF_top100_discriminative_genes.csv"))
      
      # Build heatmap matrix for Top-100 (rows = genes (ENSEMBL), cols = substances)
      # Use the same MoA-wise masking logic as above.
      heat_long <- tidy_long %>%
        dplyr::filter(!is.na(Mode_of_action),
                      ENSEMBL %in% keep_genes_rf) %>%
        dplyr::left_join(
          rf_gene_elig %>%
            dplyr::select(ENSEMBL, Mode_of_action, pass),
          by = c("ENSEMBL", "Mode_of_action")
        ) %>%
        dplyr::mutate(
          GENE_ID = ENSEMBL,
          GENE_LABEL = dplyr::if_else(is.na(SYMBOL) | SYMBOL == "", ENSEMBL, SYMBOL),
          lfc_masked = dplyr::if_else(!is.na(pass) & pass, lfc, 0)
        ) %>%
        dplyr::filter(GENE_ID %in% top100$GENE_ID) %>%
        dplyr::group_by(GENE_ID, substance) %>%
        dplyr::summarise(
          lfc  = mean(lfc_masked,  na.rm = TRUE),
          padj = if (all(is.na(padj))) NA_real_ else min(padj, na.rm = TRUE),
          .groups = "drop"
        )
      
      # Wide LFC and padj matrices
      heat_lfc_wide <- heat_long %>%
        dplyr::select(GENE_ID, substance, lfc) %>%
        tidyr::pivot_wider(
          names_from = substance,
          values_from = lfc
        )
      
      heat_padj_wide <- heat_long %>%
        dplyr::select(GENE_ID, substance, padj) %>%
        tidyr::pivot_wider(
          names_from = substance,
          values_from = padj
        )
      
      heat_wide <- heat_lfc_wide %>%
        tibble::column_to_rownames("GENE_ID") %>%
        as.matrix()
      
      padj_wide <- heat_padj_wide %>%
        tibble::column_to_rownames("GENE_ID") %>%
        as.matrix()
      
      # Set LFC to 0 for non-significant cells (padj > 0.05 OR NA)
      heat_wide[is.na(padj_wide) | padj_wide > 0.05] <- 0
      
      # Create SYMBOL lookup for rows
      gene_map <- tidy_long %>% dplyr::distinct(ENSEMBL, SYMBOL)
      
      if (nrow(heat_wide) > 1 && ncol(heat_wide) > 1) {
        # Impute remaining NAs by gene-wise median for plotting
        for (i in seq_len(nrow(heat_wide))) {
          r <- heat_wide[i, ]
          med <- stats::median(r, na.rm = TRUE)
          r[is.na(r)] <- med
          heat_wide[i, ] <- r
        }
        
        # Robust annotations via normalized join
        anno_col_rf <- tibble::tibble(substance = colnames(heat_wide)) %>%
          dplyr::mutate(substance_norm = norm_subst(substance)) %>%
          dplyr::left_join(moa_tbl %>% dplyr::select(substance_norm, Mode_of_action), by = "substance_norm") %>%
          dplyr::select(substance, Mode_of_action) %>%
          tibble::column_to_rownames("substance")
        
        if (!FILL_UNKNOWN_IN_PLOTS) {
          keep_cols_rf <- rownames(anno_col_rf)[!is.na(anno_col_rf$Mode_of_action) & anno_col_rf$Mode_of_action != ""]
          heat_wide <- heat_wide[, keep_cols_rf, drop = FALSE]
          anno_col_rf <- anno_col_rf[keep_cols_rf, , drop = FALSE]
        }
        if (FILL_UNKNOWN_IN_PLOTS) {
          anno_col_rf$Mode_of_action[is.na(anno_col_rf$Mode_of_action) | anno_col_rf$Mode_of_action == ""] <- "Unknown"
        }
        
        moa_levels_rf <- unique(as.character(anno_col_rf$Mode_of_action))
        anno_col_rf$Mode_of_action <- factor(anno_col_rf$Mode_of_action, levels = moa_levels_rf)
        palette_moa_rf <- get_moa_color_map(moa_levels_rf)
        ann_colors_rf <- list(Mode_of_action = palette_moa_rf)
        
        # Default clustered RF heatmap
        png(filename = file.path(OUTPUT_DIR, "heatmap_top100_rf_genes.png"), width = 1400, height = 1100, res = 150)
        pheatmap::pheatmap(heat_wide, scale = "row",
                           clustering_distance_rows = "euclidean",
                           clustering_distance_cols = "euclidean",
                           clustering_method = "complete",
                           annotation_col = anno_col_rf,
                           annotation_colors = ann_colors_rf,
                           main = "Top 100 MoA-discriminating genes (RF; MoA-masked, non-significant LFC set to 0)")
        dev.off()
        
        pdf(file = file.path(OUTPUT_DIR, "heatmap_top100_rf_genes.pdf"), width = 11, height = 8.5)
        pheatmap::pheatmap(heat_wide, scale = "row",
                           clustering_distance_rows = "euclidean",
                           clustering_distance_cols = "euclidean",
                           clustering_method = "complete",
                           annotation_col = anno_col_rf,
                           annotation_colors = ann_colors_rf,
                           main = "Top 100 MoA-discriminating genes (RF; MoA-masked, non-significant LFC set to 0)")
        dev.off()
        
        # ----------------------------
        # Ordered RF heatmap (ComplexHeatmap) per requested order & styling
        # ----------------------------
        # User-specified substance order (left to right)
        desired_order <- c(
          "6-PTU", "Methimazole", "Iopanoic_acid", "T3",
          "Methyltestosterone", "Trenbolone", "Flutamid",
          "Estradiol", "Ethinylestradiol", "Tamoxifen", 
          "Tebuconazole", "Difenoconazole", "Epoxiconazole", 
          "DNOC", "Fenbutatin_oxide", "Fenazaquin", "Boscalid", "Metalaxyl",
          "Carbaryl", "Chlorpyrifos", "Imidacloprid", "Nicotine",
          "Methoxychlor", "Abamectin", "Endosulfan", "Fipronil"
        )
        cols_present <- intersect(desired_order, colnames(heat_wide))
        if (length(cols_present) >= 2) {
          heat_ord <- heat_wide[, cols_present, drop = FALSE]
          # Clip values to [-2, 2] to strengthen coloring and avoid outlier domination
          heat_ord_clipped <- pmin(pmax(heat_ord, -2), 2)
          
          # Row labels: SYMBOL (fallback to ENSEMBL)
          row_symbols <- gene_map$SYMBOL[match(rownames(heat_ord_clipped), gene_map$ENSEMBL)]
          row_symbols[is.na(row_symbols) | row_symbols == ""] <- rownames(heat_ord_clipped)
          
          # MoA annotation for ordered columns
          ann_rf_ord <- tibble::tibble(substance = colnames(heat_ord_clipped)) %>%
            dplyr::mutate(substance_norm = norm_subst(substance)) %>%
            dplyr::left_join(moa_tbl %>% dplyr::select(substance_norm, Mode_of_action), by = "substance_norm") %>%
            dplyr::select(substance, Mode_of_action)
          
          if (!FILL_UNKNOWN_IN_PLOTS) {
            keep_cols_ord <- ann_rf_ord$substance[!is.na(ann_rf_ord$Mode_of_action) & ann_rf_ord$Mode_of_action != ""]
            heat_ord_clipped <- heat_ord_clipped[, keep_cols_ord, drop = FALSE]
            ann_rf_ord <- ann_rf_ord[ann_rf_ord$substance %in% keep_cols_ord, , drop = FALSE]
          } else {
            ann_rf_ord$Mode_of_action[is.na(ann_rf_ord$Mode_of_action) | ann_rf_ord$Mode_of_action == ""] <- "Unknown"
          }
          
          moa_vec <- ann_rf_ord$Mode_of_action
          names(moa_vec) <- ann_rf_ord$substance
          
          moa_levels_ord <- unique(as.character(moa_vec))
          col_map <- get_moa_color_map(moa_levels_ord)
          
          # ---- NEW: compute separator positions between MoA groups ----
          moa_groups <- rle(as.character(moa_vec))
          block_lengths <- moa_groups$lengths
          block_ends <- cumsum(block_lengths)
          # boundaries between MoA groups (ignore last group end)
          separator_positions <- block_ends[-length(block_ends)]
          
          # MoA annotation
          top_ann <- ComplexHeatmap::HeatmapAnnotation(
            Mode_of_action = factor(moa_vec, levels = moa_levels_ord),
            col = list(Mode_of_action = col_map),
            annotation_legend_param = list(title = "Mode of action")
          )
          
          # Diverging color scale with fixed limits at ±2 for LFC
          col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#3B4CC0", "#FFFFFF", "#B40426"))
          
          # Draw heatmap: fixed columns, cluster rows; show SYMBOLs on right; rotate x-labels 90° clockwise (-90)
          ht <- ComplexHeatmap::Heatmap(
            heat_ord_clipped,
            name = "lfc",
            col = col_fun,
            cluster_columns = FALSE,
            cluster_rows = TRUE,
            top_annotation = top_ann,
            show_row_names = TRUE,
            row_names_side = "right",
            row_names_gp = grid::gpar(fontsize = 8),  # slightly larger for readability
            row_labels = row_symbols,
            column_names_rot = -90,
            column_names_gp = grid::gpar(fontsize = 8),
            heatmap_legend_param = list(title = "LFC", legend_direction = "vertical")
          )
          
          # PNG output with increased height
          png(
            filename = file.path(OUTPUT_DIR, "heatmap_top100_rf_genes_ORDERED.png"),
            width = 1400,
            height = 1300,  # increased vertical space
            res = 150
          )
          
          ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
          
          # Add vertical black lines between MoA groups *after* drawing
          ComplexHeatmap::decorate_heatmap_body("lfc", {
            n_cols <- ncol(heat_ord_clipped)
            for (b in separator_positions) {
              grid::grid.lines(
                x = grid::unit(rep(b / n_cols, 2), "npc"),
                y = grid::unit(c(0, 1), "npc"),
                gp = grid::gpar(col = "black", lwd = 1.2)
              )
            }
          })
          
          dev.off()
          
          # PDF output with increased height
          pdf(
            file = file.path(OUTPUT_DIR, "heatmap_top100_rf_genes_ORDERED.pdf"),
            width = 11,
            height = 10    # increased vertical space
          )
          
          ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
          
          ComplexHeatmap::decorate_heatmap_body("lfc", {
            n_cols <- ncol(heat_ord_clipped)
            for (b in separator_positions) {
              grid::grid.lines(
                x = grid::unit(rep(b / n_cols, 2), "npc"),
                y = grid::unit(c(0, 1), "npc"),
                gp = grid::gpar(col = "black", lwd = 1.2)
              )
            }
          })
          
          dev.off()
        } else {
          message("Ordered RF heatmap skipped: fewer than 2 of the requested substances are present in data.")
        }
        
        # Top-30 importance barplot
        top30 <- top100 %>% dplyr::slice_head(n = min(30, nrow(top100)))
        bp <- ggplot2::ggplot(top30, ggplot2::aes(x = stats::reorder(GENE_ID, rf_importance_mean), y = rf_importance_mean)) +
          ggplot2::geom_col() + ggplot2::coord_flip() +
          ggplot2::labs(title = "Random Forest feature importance (Top 30)", x = "Gene", y = "MeanDecreaseGini") +
          ggplot2::theme_minimal(base_size = 11)
        ggplot2::ggsave(file.path(OUTPUT_DIR, "RF_top30_importance_barplot.png"), bp, width = 8, height = 10, dpi = 150)
      } else {
        warning("Heatmap for Top-100 genes skipped: insufficient matrix size after filtering.")
      }
    }
  }
}

# ----------------------------
# Session info & README
# ----------------------------
sink(file.path(OUTPUT_DIR, "sessionInfo.txt")); print(sessionInfo()); sink()

readme_text <- c(
  "# Outputs generated by analyze_zebrafish_expression.R",
  "",
  "## Filters applied",
  "1) Keep genes with any padj < 0.05 across substances.",
  "2) Keep genes that are padj <= 0.05 for **all** substances of at least one MoA.",
  "3) Keep genes that, within such an MoA, have consistent LFC sign across its substances.",
  "",
  "## Random Forest",
  "- Gene eligibility: only genes that are significantly regulated (padj <= 0.05)",
  "  within at least one MoA and with consistent LFC sign (all >=0 or all <=0) are kept.",
  "- MoA-wise masking: for (gene, MoA) pairs that do NOT meet this criterion, all LFCs",
  "  of that gene for substances in that MoA are set to 0 before RF.",
  "",
  "## Tables",
  "- tidy_long_table.csv",
  "- sig_counts_per_substance.csv",
  "- sig_counts_per_moa.csv",
  "- RF_feature_importance_all_genes.csv",
  "- RF_top100_discriminative_genes.csv",
  "- substances_missing_moa.csv (diagnostics)",
  "",
  "## Plots",
  "- volcano_plots_faceted.png",
  "- heatmap_top_variable_genes.png",
  "- heatmap_top_variable_genes.pdf",
  "- heatmap_top100_rf_genes.png",
  "- heatmap_top100_rf_genes.pdf",
  "- heatmap_top100_rf_genes_ORDERED.png",
  "- heatmap_top100_rf_genes_ORDERED.pdf",
  "- RF_top30_importance_barplot.png",
  "",
  "## Notes",
  "- Parsing set to ';' delimiter. Decimal comma handled for lfc_all via readr::locale(decimal_mark=\",\").",
  "- Normalized join by 'substance' to improve matching to MoA. Toggle FILL_UNKNOWN_IN_PLOTS to drop or fill unmapped substances.",
  "- Ordered RF heatmap uses Set2/fill_map-based colors for Mode_of_action (matching the GSEA script).",
  "- In all heatmaps, LFC values for padj > 0.05 or NA are set to 0 before plotting.",
  "- For RF, additional MoA-wise masking sets LFCs to 0 where the (gene, MoA) pair fails the significance+sign criterion."
)
writeLines(readme_text, con = file.path(OUTPUT_DIR, "README_outputs.txt"))

message("Done. Outputs are in the 'outputs' directory.")
