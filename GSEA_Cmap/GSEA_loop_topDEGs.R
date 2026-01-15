######### Load necessary libraries ###########
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(Hmisc)
library(fgsea)
library(stats)

######### Load data ###########
#load a list of HE, LE or intersect dataframes <- subsetted_dataframes

# Set this as the main input for the rest of the analysis
subsetted_dataframes <- le_list

######### Functions ###########

prepare_ranked_list <- function(df) {
  df <- df %>% filter(!is.na(stat))
  ranks <- df$stat
  names(ranks) <- df$gene
  sort(ranks, decreasing = TRUE)
}

create_gene_sets <- function(df, padj_cutoff = 0.05, lfc_cutoff = 0) {
  up_genes <- df %>% filter(padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>% pull(gene)
  down_genes <- df %>% filter(padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>% pull(gene)
  list(up = up_genes, down = down_genes)
}

run_connectivity <- function(rank_list, gene_sets) {
  res_up <- fgsea(pathways = list(up = gene_sets$up), stats = rank_list)
  res_down <- fgsea(pathways = list(down = gene_sets$down), stats = rank_list)
  
  ES_up <- ifelse(length(res_up$ES) == 0, 0, res_up$ES[1])
  ES_down <- ifelse(length(res_down$ES) == 0, 0, res_down$ES[1])
  NES_up <- ifelse(length(res_up$NES) == 0, 0, res_up$NES[1])
  NES_down <- ifelse(length(res_down$NES) == 0, 0, res_down$NES[1])
  p_up <- ifelse(length(res_up$pval) == 0, 1, res_up$pval[1])
  p_down <- ifelse(length(res_down$pval) == 0, 1, res_down$pval[1])
  
  raw_score <- ES_up - ES_down
  norm_score <- NES_up - NES_down
  combined_pval <- max(p_up, p_down)
  
  return(list(
    raw_score = raw_score,
    norm_score = norm_score,
    pval = combined_pval,
    NES_up = NES_up,
    NES_down = NES_down
  ))
}

######### Main loop over top N ###########
top_n_values <- c(20, 50, 100)

for (n_top in top_n_values) {
  cat("=== Running analysis for top", n_top, "DEGs ===\n")
  
  # Top N selection functions
  get_top_n_padj <- function(df) {
    df[order(df$padj), ][1:min(n_top, nrow(df)), ]
  }
  get_top_n_log2fold <- function(df) {
    df[order(abs(df$log2FoldChange), decreasing = TRUE), ][1:min(n_top, nrow(df)), ]
  }
  get_top_n_stat <- function(df) {
    df <- df %>%
      mutate(stat = log2FoldChange / lfcSE)
    df[order(abs(df$stat), decreasing = TRUE), ][1:min(n_top, nrow(df)), ]
  }
  
  # Apply selection
  top_padj_list <- lapply(subsetted_dataframes, get_top_n_padj) %>% lapply(na.omit)
  top_log2fold_list <- lapply(subsetted_dataframes, get_top_n_log2fold) %>% lapply(na.omit)
  top_stat_list <- lapply(subsetted_dataframes, get_top_n_stat) %>% lapply(na.omit)
  
  gene_vectors <- list(
    padj = unique(unlist(lapply(top_padj_list, rownames))),
    log2fold = unique(unlist(lapply(top_log2fold_list, rownames))),
    stat = unique(unlist(lapply(top_stat_list, rownames)))
  )
  
  for (vector_name in names(gene_vectors)) {
    cat("  → Processing", vector_name, "genes...\n")
    selected_genes <- gene_vectors[[vector_name]]
    
    # Subset HE dataframes
    filtered_data_list <- lapply(subsetted_dataframes, function(df) {
      df[rownames(df) %in% selected_genes, , drop = FALSE]
    })
    
    # Add stat column
    df_list_with_stat <- lapply(filtered_data_list, function(df) {
      df$stat <- df$log2FoldChange / df$lfcSE
      df
    })
    
    # Keep selected columns (2nd, 8th, 5th)
    columns_to_keep_indices <- c(2, 8, 5)
    df_list_subset <- lapply(df_list_with_stat, function(df) {
      df <- df[, columns_to_keep_indices]
      df$gene <- rownames(df)
      rownames(df) <- NULL
      return(df)
    })
    
    test_list <- df_list_subset
    ref_list <- df_list_subset
    
    # Run fgsea connectivity
    results_all <- data.frame()
    
    for (ref_name in names(ref_list)) {
      ref_ranked <- prepare_ranked_list(ref_list[[ref_name]])
      
      for (test_name in names(test_list)) {
        test_gene_sets <- create_gene_sets(test_list[[test_name]])
        res <- run_connectivity(ref_ranked, test_gene_sets)
        
        results_all <- rbind(results_all,
                             data.frame(
                               reference = ref_name,
                               test = test_name,
                               raw_score = res$raw_score,
                               norm_score = res$norm_score,
                               NES_up = res$NES_up,
                               NES_down = res$NES_down,
                               pval = res$pval
                             ))
      }
    }
    
    # Save results
    output_filename <- paste0("connectivity_results_top", n_top, "_", vector_name, ".RData")
    save(results_all, file = output_filename)
    cat("    ✔ Saved:", output_filename, "\n\n")
  }
}
