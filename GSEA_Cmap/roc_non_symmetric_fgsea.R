###########################
library(fgsea)
library(dplyr)
library(stats)



#results_all$reference <- sub("\\..*$", ".ref", results_all$reference)
#results_all$test <- sub("\\..*$", ".test", results_all$test)


############################

#results_all$padj <- p.adjust(results_all$pval, method = "BH")
########################
library(tidyr)

wide_raw_scores <- results_all %>%
  select(reference, test, raw_score) %>%
  pivot_wider(names_from = test, values_from = raw_score)

wide_norm_scores <- results_all %>%
  select(reference, test, norm_score) %>%
  pivot_wider(names_from = test, values_from = norm_score)

wide_pvals <- results_all %>%
  select(reference, test, pval) %>%
  pivot_wider(names_from = test, values_from = pval)

#wide_padj <- results_all %>%
 # select(reference, test, padj) %>%
 # pivot_wider(names_from = test, values_from = padj)
####################
# Load necessary package
library(tibble)


# Convert first column to row names
df <- as.data.frame(wide_norm_scores)                     # Convert tibble to data.frame
df<- df[-c(33),]
rownames(df) <- df[[1]]                     # Set row names from first column
df <- df[,-c(1,33)]
colnames(df) <- sub("\\..*$", "", colnames(df))
rownames(df) <- sub("\\..*$", "", rownames(df))




#########
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

conn_mat <- df  # Use loaded matrix

# Simulate instead, if needed:
# set.seed(123)
# compounds <- moa_map$compound
# conn_mat <- matrix(runif(length(compounds)^2, min = -1, max = 1), 
#                    nrow = length(compounds), ncol = length(compounds),
#                    dimnames = list(compounds, compounds))
# diag(conn_mat) <- 1

# --------- Step 2: Create Label Matrix ---------

# --------- Step 1: Create label matrix ---------
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

# --------- Step 2: Predict Shared MoA Based on Threshold ---------

threshold <- 0  # You can change this threshold
pred_mat <- matrix(as.integer(conn_mat > threshold), n, n)
diag(pred_mat) <- NA

# --------- Step 3: Compute Classification Metrics ---------

# Use all off-diagonal entries (excluding NAs)
valid_idx <- !is.na(label_mat) & !is.na(pred_mat)
labels <- label_mat[valid_idx]
preds  <- pred_mat[valid_idx]

TP <- sum(preds == 1 & labels == 1)
FP <- sum(preds == 1 & labels == 0)
TN <- sum(preds == 0 & labels == 0)
FN <- sum(preds == 0 & labels == 1)

# Avoid division by zero
safe_div <- function(x, y) ifelse((x + y) == 0, NA, x / (x + y))

TPR <- safe_div(TP, FN)           # Sensitivity / Recall
FPR <- safe_div(FP, TN)           # False Positive Rate
TNR <- safe_div(TN, FP)           # Specificity
FNR <- safe_div(FN, TP)           # Miss Rate
accuracy <- safe_div(TP + TN, TP + TN + FP + FN)
error_rate <- 1 - accuracy

# --------- Output ---------
cat(sprintf("Threshold > 0:\nTPR = %.3f, FPR = %.3f, TNR = %.3f, FNR = %.3f\n, Accuracy = %.3f",
            TPR, FPR, TNR, FNR, accuracy))

rm(list = ls())

