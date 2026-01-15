######### Load necessary libraries###########
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(Hmisc)



load("y:/manuscripts/manuscript 3/1- data_reanalysis/desec2_data_updated.RData")
#load("~/../Desktop/man3/figures/fig1/test_last_meeting/degpcit.RData")

list_of_dfs<- degPcut.ls
df_names <- names(list_of_dfs)

# Find the names of the substances
substances <- unique(sub(".\\w+$", "", df_names))

# Initialize an empty list to store the final dataframes
final_dfs <- list()

# Loop through each substance
for (substance in substances) {
  # Find the conditions for this substance
  conditions <- sub("^.*\\.", "", df_names[grep(substance, df_names)])
  
  # If there's a Mid condition, keep only High and Mid
  if ("Mid" %in% conditions) {
    final_dfs <- c(final_dfs, list_of_dfs[paste0(substance, ".High")], list_of_dfs[paste0(substance, ".Mid")])
  } else {
    # If not, keep all conditions
    final_dfs <- c(final_dfs, list_of_dfs[grep(substance, df_names)])
  }
}

# The final_dfs list now contains the filtered dataframes
#####get the core DEGs for each substance#####
get_intersection <- function(df1, df2) {
  intersect(rownames(df1), rownames(df2))
}

# Initialize a list to store intersections
intersections <- list()
degsel.ls<- final_dfs
# Loop through the list in pairs
for (i in seq(1, length(degsel.ls), by=2)) {
  intersections[[length(intersections) + 1]] <- get_intersection(degsel.ls[[i]], degsel.ls[[i + 1]])
}
################################
###for androstendion and spirotetramat selection will be based on HE due to 
#low nomber of DEGs at LE --> that also indicates that HE responses still specific
#intersections[[19]]<- rownames(final_dfs[["Androstendion.High"]])
#intersections[[22]]<- rownames(final_dfs[["Spirotetramat.High"]])


#########################################################
####subset HE dataframes####
#####keep only HE dataframes in the list#####
degsel.ls  <- degsel.ls[grep("\\.High$", names(degsel.ls))]
####subset the HE dataframes for the defined core DEGs######
subsetted_dataframes <- lapply(seq_along(degsel.ls), function(i) {
  degsel.ls[[i]][rownames(degsel.ls[[i]]) %in% intersections[[i]], , drop = FALSE]
})
names(subsetted_dataframes) <- names(degsel.ls)
subsetted_dataframes<- subsetted_dataframes[-c(14,3,4,29,33,35,37,39)]

#####################################################################
#####select for each top 20 DEGs#####
# Function to select top 20 rows by pvalue
get_top_20 <- function(df) {
  df[order(df$padj), ][1:20, ]
}
get_top_log2fold <- function(df) {
  # Order by absolute values of log2FoldChange and get the top 2
  top_genes <- df[order(abs(df$log2FoldChange), decreasing = TRUE), ][1:20, ]
  
  return(top_genes)
}
get_top_stat <- function(df, n = 20) {
  # Calculate Wald statistic (log2FoldChange / SE)
  df <- df %>%
    dplyr::mutate(stat = log2FoldChange / lfcSE)
  
  # Order by absolute value of the statistic and return top n
  top_genes <- df %>%
    dplyr::arrange(desc(abs(stat))) %>%
    dplyr::slice(1:n)
  
  return(top_genes)
}

# Initialize a list to store top 20 dataframes
top_20_list <- lapply(subsetted_dataframes, get_top_20)
top_20_list<- lapply(top_20_list, na.omit)
top_20_log2fold <- lapply(subsetted_dataframes, get_top_log2fold)
top_20_log2fold<- lapply(top_20_log2fold, na.omit)
top_20_stat <- lapply(subsetted_dataframes, get_top_stat)
top_20_stat<- lapply(top_20_stat, na.omit)

###################################################

# For a list of dataframes
rownames_vector_padj<- unique(unlist(lapply(top_20_list, rownames)))
rownames_vector_log2fold<- unique(unlist(lapply(top_20_log2fold, rownames)))
rownames_vector_stat<- unique(unlist(lapply(top_20_stat, rownames)))
rownames_vector_all<- unique(unlist(lapply(subsetted_dataframes, rownames)))

#######################
filtered_data_list <- lapply(subsetted_dataframes, function(df) df[rownames(df) %in% rownames_vector_all, ])
modified_data_list <- list()

# Iterate through each data frame in the filtered list
for (i in seq_along(filtered_data_list)) {
  df <- filtered_data_list[[i]]
  fungicide_name <- names(filtered_data_list)[i]  # Get the name of the data frame
  
  # Select column 1 and column 3
  df$x<- rownames(df)
  rownames(df)<- NULL
  modified_df <- df[, c(8,2)]  # Using base R to select columns
  
  # Rename column 3 to the name of the data frame
  colnames(modified_df)[2] <- fungicide_name
  #rownames(modified_df)<- modified_df$X
  #modified_df$X<- NULL
  
  # Store the modified data frame in the new list
  modified_data_list[[fungicide_name]] <- modified_df
}



merged_final_data <- Reduce(function(x, y) full_join(x, y, by = names(x)[1]), modified_data_list)
# Replace NAs with 0 or any other value if necessary
merged_final_data[is.na(merged_final_data)] <- 0
rownames(merged_final_data)<- merged_final_data$x
merged_final_data$x<- NULL
colnames(merged_final_data)<- sub("\\..*$", "", colnames(merged_final_data))

#########calculate the correlation matrix######
#correlation_matrix <- cor(merged_final_data)
#matrix <- correlation_matrix[,c(1,2,3,4,5,6)]
#correlation_matrix[is.na(cor_matrix)] <- 0
correlation_matrix <- rcorr(as.matrix(merged_final_data), type = "pearson")
##############################
cor_matrix <- correlation_matrix$r
p_matrix <- correlation_matrix$P

# Set correlation values to zero where p-value < 0.05
cor_matrix[p_matrix > 0.05] <- 0



#corr_mat<- cor_matrix[-c(19),-c(19)] 
######plot heatmap of correlation matrix####
cor_matrix[is.na(cor_matrix)] <- 0

mat<- as.matrix(cor_matrix)
#conn_mat<- mat
# --------- Step 0: Load Data and Define MoA Mapping ---------
#load("~/../Desktop/man3/mat.RData")

moa_groups <- list(
  acetylcholine_modulator = c("Imidacloprid.High", "Nicotine.High"),
  acetylcholine_blocker = c("Cartaphydrochloride.High"),
  Acetylcholinesterase_inhibitor = c("Chlorpyrifos.High", "Carbaryl.High"),
  sodium_channel_modulators = c("Methoxychlor.High"),
  gaba_channel_blockers = c("Fipronil.High", "Endosulfan.High"),
  glutamate_chloride_channel_blocker = c("Abamectin.High"),
  antiandrogen = c("Flutamid.High"),
  acetylcoa_carboxylase_inhibitor = c("Spirotetramat.High"),
  oxidative_phosphorylation_uncoupler = c("DNOC.High"),
  SDH_inhibitor = c("Boscalid.High"),
  mitochondrial_inhibitor = c("Fenbutatin Oxide.High", "Fenazaquin.High"),
  RNA_metabolism_inhibitor = c("Metalaxyl.High"),
  androgenic = c("Methyltestosterone.High", "Androstendion.High", "Trenbolone.High"),
  antiestrogenic = c("Fulvestrant.High", "Tamoxifen.High"),
  estrogenic = c("Estradiol.High", "Ethinylestradiol.High", "BisphenolA.High"),
  thyroid_TPO_inhibitor = c("Methimazole.High", "6-PTU.High"),
  thyroid_DIO_inhibitor = c("Iopanoic_acid.High"),
  thyroid_hormone = c("T3.High"),
  sterol_biosynthesis_inhibitor = c("Prochloraz.High", "Epoxiconazole.High", "Tebuconazole.High", "Difenoconazole.High")
)

moa_map <- do.call(rbind, lapply(names(moa_groups), function(moa) {
  data.frame(
    compound = moa_groups[[moa]],
    moa = moa,
    stringsAsFactors = FALSE
  )
}))
moa_map$compound <- sub("\\..*$", "", moa_map$compound)


# --------- Step 1: Prepare the Connectivity Matrix ---------

conn_mat <- mat  # Use loaded matrix

# Simulate instead, if needed:
# set.seed(123)
# compounds <- moa_map$compound
# conn_mat <- matrix(runif(length(compounds)^2, min = -1, max = 1), 
#                    nrow = length(compounds), ncol = length(compounds),
#                    dimnames = list(compounds, compounds))
# diag(conn_mat) <- 1

# --------- Step 2: Create Label Matrix ---------

compound_names <- rownames(conn_mat)
moa_list <- setNames(moa_map$moa, moa_map$compound)[compound_names]

n <- length(compound_names)
label_mat <- matrix(0, n, n, dimnames = list(compound_names, compound_names))

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j && moa_list[i] == moa_list[j]) {
      label_mat[i, j] <- 1
    }
  }
}
diag(label_mat) <- NA  # Ignore self-comparisons

# --------- Step 3: Predict Shared MoA Based on Threshold ---------

pred_mat <- matrix(as.integer(conn_mat > 0), n, n)
diag(pred_mat) <- NA

# --------- Step 4: Compute Standard Classification Metrics ---------

labels <- label_mat[upper.tri(label_mat)]
preds <- pred_mat[upper.tri(pred_mat)]

TP <- sum(preds == 1 & labels == 1)
FP <- sum(preds == 1 & labels == 0)
TN <- sum(preds == 0 & labels == 0)
FN <- sum(preds == 0 & labels == 1)

# ✅ Standard definitions
TPR <- TP / (TP + FN)  # Sensitivity / Recall
FPR <- FP / (FP + TN)  # False Positive Rate
TNR <- TN / (TN + FP)  # Specificity
FNR <- FN / (FN + TP)  # Miss Rate
accuracy <- (TP + TN) / (TP + TN + FP + FN)
error_rate <- 1 - accuracy

cat(sprintf("Threshold > 0:\nTPR = %.3f, FPR = %.3f, TNR = %.3f, FNR = %.3f\n",
            TPR, FPR, TNR, FNR))
cat(sprintf("Accuracy = %.3f, Error Rate = %.3f\n", accuracy, error_rate))

