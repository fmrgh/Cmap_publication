# =========================
# Packages
# =========================
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(forcats)
library(randomForest)
library(RColorBrewer)
library(stringr)

# =========================
# Parameters (easily adjustable)
# =========================
infile   <- "gsea_results_he.csv"
out_topN <- 50                      # Top-N GO terms based on RF importance
pdf_w1   <- 16; pdf_h1 <- 9         # Size of Plot 1 (Substance on x-axis)
pdf_w2   <- 14; pdf_h2 <- 9         # Size of Plot 2 (Mode_of_action on x-axis)

# =========================
# 0) Robust CSV reading (DE/EN formats)
#    Tries semicolon with comma-decimal and semicolon with dot-decimal.
# =========================
.read_gsea <- function(path) {
  try1 <- try(
    read_delim(path, delim = ";",
               locale = locale(decimal_mark = ","),
               col_types = cols(.default = col_character())),
    silent = TRUE
  )
  if (!inherits(try1, "try-error")) return(try1)
  
  try2 <- try(
    read_delim(path, delim = ";",
               locale = locale(decimal_mark = "."),
               col_types = cols(.default = col_character())),
    silent = TRUE
  )
  if (!inherits(try2, "try-error")) return(try2)
  
  stop("Could not read file robustly. Please check file format.")
}

# =========================
# 1) Read & clean data
# =========================
df_filtered <- .read_gsea(infile)
colnames(df_filtered) <- c(
  "Substance", "Mode_of_action", "Group", "GO_term",
  "GO_term_description", "NES", "p_adjust", "rank"
)

df_filtered <- df_filtered %>%
  mutate(
    GO_term = na_if(GO_term, "NA"),
    GO_term_description = na_if(GO_term_description, "NA")
  ) %>%
  filter(!is.na(GO_term), !is.na(GO_term_description))

# Numeric parser that tolerates comma or dot decimals
parse_num <- function(x) suppressWarnings(as.numeric(str_replace_all(x, ",", ".")))

# =========================
# 2) Presence/absence matrix (RF input)
# =========================
mat <- df_filtered %>%
  distinct(Substance, Mode_of_action, Group, GO_term) %>%
  mutate(present = 1, GO_col_name = make.names(GO_term)) %>%
  filter(!is.na(GO_col_name), GO_col_name != "NA") %>%
  pivot_wider(names_from = GO_col_name, values_from = present, values_fill = 0) %>%
  mutate(
    Mode_of_action = as.factor(Mode_of_action),
    Group = as.factor(Group)
  )

# =========================
# 3) Random Forest importance on Mode_of_action
# =========================
set.seed(42)
rf <- randomForest(Mode_of_action ~ . - Substance - Group,
                   data = mat, ntree = 100, importance = TRUE)
imp <- importance(rf)
rf_table <- as.data.frame(imp)
rf_table$GO_col_name <- rownames(rf_table)

# Top-N GO terms by MeanDecreaseGini
top_rf_terms <- rf_table %>%
  filter(GO_col_name != "NA") %>%
  arrange(desc(MeanDecreaseGini)) %>%
  slice_head(n = out_topN) %>%
  pull(GO_col_name)

# GO descriptions for Top-N
top_rf_go <- df_filtered %>%
  mutate(GO_col_name = make.names(GO_term)) %>%
  filter(GO_col_name %in% top_rf_terms, !is.na(GO_term_description)) %>%
  distinct(GO_col_name, GO_term_description)

# =========================
# 4) Plot data (restricted to Top-N GO terms)
# =========================
plot_data <- df_filtered %>%
  mutate(
    GO_col_name = make.names(GO_term),
    p_adjust_num = parse_num(p_adjust),
    NES_num = parse_num(NES),
    neg_log_padj = -log10(p_adjust_num)
  ) %>%
  filter(GO_col_name %in% top_rf_terms, !is.na(GO_term_description))

# Keep only GO terms that appear in ALL substances for a given Mode_of_action
sub_per_moa <- df_filtered %>%
  group_by(Mode_of_action) %>%
  summarise(total_sub = n_distinct(Substance), .groups = "drop")

shared_go_terms <- df_filtered %>%
  group_by(Mode_of_action, GO_term) %>%
  summarise(substance_count = n_distinct(Substance), .groups = "drop") %>%
  inner_join(sub_per_moa, by = "Mode_of_action") %>%
  filter(substance_count == total_sub) %>%
  select(Mode_of_action, GO_term)

plot_data <- plot_data %>%
  inner_join(shared_go_terms, by = c("Mode_of_action", "GO_term"))

# =========================
# 5) Best-discriminating Mode_of_action per GO term (for clustering)
# =========================
mode_gini <- mat %>%
  select(Mode_of_action, all_of(top_rf_terms)) %>%
  group_by(Mode_of_action) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  pivot_longer(-Mode_of_action, names_to = "GO_col_name", values_to = "count") %>%
  group_by(GO_col_name) %>%
  slice_max(order_by = count, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(GO_col_name, best_Mode_of_action = Mode_of_action)

plot_data <- plot_data %>%
  left_join(mode_gini, by = "GO_col_name")

# =========================
# 6) Axis orders
# =========================
# y-axis order: by best_Mode_of_action then by significance
go_order <- plot_data %>%
  arrange(best_Mode_of_action, desc(neg_log_padj)) %>%
  pull(GO_term_description) %>%
  unique()

plot_data <- plot_data %>%
  mutate(GO_term_description = factor(GO_term_description, levels = go_order))

# Manual sorting order for Mode_of_action (as requested)
custom_moa_order <- c(
  "acetylcholine_esterase_inhibitor",
  "acetylcholine_receptor_agonist",
  "sodium_channel_modulator",
  "glutamate_chloride_channel_blocker",
  "gaba_channel_blocker",
  "TPO_inhibitor",
  "DIO_inhibitor",
  "thyroid_hormone_receptor_agonist",
  "androgen_receptor_agonist",
  "androgen_receptor_antagonist",
  "estrogen_receptor_agonist",
  "estrogen_receptor_antagonist",
  "sterol_biosynthesis_inhibitor",
  "oxidative_phosphorylation_uncoupler",
  "oxidative_phosphorylation_inhibitor",
  "SDH_inhibitor",
  "RNA_metabolism_inhibitor"
)

# x-axis substance order within Group > Mode_of_action > Substance
plot_data <- plot_data %>%
  mutate(Mode_of_action = factor(Mode_of_action, levels = custom_moa_order)) %>%
  arrange(Group, Mode_of_action, Substance) %>%     # now uses your order
  mutate(Substance = factor(Substance, levels = unique(Substance)))

# =========================
# 7) Background panels for Mode_of_action groups on x-axis
# =========================
substance_levels <- levels(plot_data$Substance)

boundaries <- plot_data %>%
  select(Substance, Mode_of_action) %>%
  distinct() %>%
  mutate(Substance_num = as.numeric(factor(Substance, levels = substance_levels))) %>%
  arrange(Substance_num) %>%
  group_by(Mode_of_action) %>%
  summarise(
    xmin = min(Substance_num) - 0.5,
    xmax = max(Substance_num) + 0.5,
    .groups = "drop"
  )

num_moa <- nrow(boundaries)
base_cols <- brewer.pal(min(8, max(3, num_moa)), "Set2")
fills <- if (num_moa > length(base_cols)) colorRampPalette(base_cols)(num_moa) else base_cols[seq_len(num_moa)]
boundaries$fill <- fills
fill_map <- setNames(boundaries$fill, boundaries$Mode_of_action)

# =========================
# 8) Dotplot per Substance
# =========================
pdf("Dotplot_GOterm_by_Modeofaction_topN_grouped.pdf", width = pdf_w1, height = pdf_h1)
ggplot(plot_data, aes(x = Substance, y = GO_term_description)) +
  geom_rect(
    data = boundaries,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Mode_of_action),
    inherit.aes = FALSE, alpha = 0.45
  ) +
  geom_point(aes(size = neg_log_padj, color = NES_num)) +
  scale_fill_manual(values = fill_map, guide = guide_legend(title = "Mode of action")) +
  scale_color_viridis_c() +
  scale_size_continuous(range = c(1, 12)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20),
    axis.text.y = element_text(size = 20),
    legend.position = "right"
  ) +
  labs(
    x = "Test compounds (sorted by Group & Mode_of_action)",
    y = "GO term description (clustered by best-differentiating Mode_of_action)",
    size = "-log10(p_adjust)",
    color = "NES",
    title = paste0(
      "Dotplot: Top ", out_topN,
      " RF-selected GO terms per compound (grouped by Group, Mode_of_action)"
    )
  )
dev.off()

# =========================
# 9) Aggregated Dotplot by Mode_of_action
# =========================
plot_data_ma <- plot_data %>%
  group_by(Mode_of_action, GO_term_description) %>%
  summarise(
    neg_log_padj = max(neg_log_padj, na.rm = TRUE),
    NES_num = mean(NES_num, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Mode_of_action = factor(Mode_of_action, levels = custom_moa_order),
    GO_term_description = factor(GO_term_description, levels = go_order)
  )

pdf("Dotplot_GOterm_by_Modeofaction_x_MA.pdf", width = pdf_w2, height = pdf_h2)
ggplot(plot_data_ma, aes(x = Mode_of_action, y = GO_term_description)) +
  geom_point(aes(size = neg_log_padj, color = NES_num)) +
  scale_color_viridis_c() +
  scale_size_continuous(range = c(1, 12)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    x = "Mode of action",
    y = "GO term description (clustered)",
    size = "-log10(p_adjust)",
    color = "NES",
    title = paste0("Dotplot: Top ", out_topN, " RF-selected GO terms grouped by Mode_of_action")
  )
dev.off()

# =========================
# 10) Useful exports
# =========================
# 1) Presence matrix used as RF input
write_csv(mat, "rf_presence_matrix.csv")

# 2) RF importance table (sorted)
rf_table %>% arrange(desc(MeanDecreaseGini)) %>% write_csv("rf_importance_sorted.csv")

# 3) Mapping GO_col_name -> description (for Top-N)
top_rf_go %>% arrange(GO_col_name) %>% write_csv("topN_go_descriptions.csv")
save.image("new_analysis.RData")
