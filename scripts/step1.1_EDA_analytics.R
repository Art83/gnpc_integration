library(dplyr)
library(ggplot2)
library(uwot)
library(ggrepel)
library(ranger)

merged_proteomics_subset <- readRDS(gzcon(file("output/proc/dpgen.rds.gz", open = "rb"))) # datatset is not in github due to user agreement restrictions 

ignore_vars <- c("contributor_code","sample_id", "sex", "race", "plate_id", 
                 "outcome", "non_na_protein_count", "inferred_assay")


# ---------- Step 1: Filter to 7K only ----------
proteomics_7k <- merged_proteomics_subset %>%
  filter(inferred_assay == "7K")

# ---------- Step 2: Impute ----------
threshold <- 0.1
proteomics_7k <- proteomics_7k %>%
  mutate(across(!all_of(ignore_vars), ~ na_if(., -1))) %>%
  {
    protein_cols <- setdiff(names(.), ignore_vars)
    na_prop <- colMeans(is.na(select(., all_of(protein_cols))))
    good_cols <- names(na_prop[na_prop < threshold])
    select(., all_of(c(good_cols, ignore_vars))) %>%
      group_by(outcome) %>%
      mutate(across(all_of(good_cols), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))
  }

# ---------- Step 4: Log2 transform ----------
proteomics_7k[, !colnames(proteomics_7k) %in% ignore_vars] <-
  log2(proteomics_7k[, !colnames(proteomics_7k) %in% ignore_vars] + 1)


# ---------- Step 5: PCA & UMAP ----------
# Extract protein data
protein_cols <- setdiff(colnames(proteomics_7k), ignore_vars)
protein_matrix <- proteomics_7k[, protein_cols]
protein_matrix <- as.matrix(protein_matrix)

# Scale the protein data
protein_scaled <- scale(protein_matrix)

# Metadata
metadata <- proteomics_7k[, c("sample_id", "contributor_code", "plate_id", "outcome")]

# PCA
pca <- prcomp(protein_scaled, center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x[, 1:3])
pca_df <- cbind(pca_df, metadata)

ggplot(pca_df, aes(x = PC1, y = PC2, color = contributor_code)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "PCA of 7K Assay Samples (Colored by Site)")

ggplot(pca_df, aes(x = PC1, y = PC2, color = outcome)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "PCA of 7K Assay Samples (Colored by Diagnosis)")

# UMAP
set.seed(42)
umap_out <- umap(protein_scaled, n_neighbors = 15, min_dist = 0.1)
umap_df <- as.data.frame(umap_out)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df <- cbind(umap_df, metadata)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = contributor_code)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "UMAP of 7K Assay Samples (Colored by Site)")

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = outcome)) +
  geom_point(alpha = 0.6) +
  theme_minimal() +
  labs(title = "UMAP of 7K Assay Samples (Colored by Diagnosis)")



# ---------- Step 6: Evaluation of Corrections ----------
# Separate features and metadata
features <- proteomics_7k[, !colnames(proteomics_7k) %in% ignore_vars]
meta <- proteomics_7k[,colnames(proteomics_7k) %in% ignore_vars]
outcome <- as.factor(meta$outcome)
site <- as.factor(meta$contributor_code)

batch_corrected_data <- list()
batch_corrected_data$raw <- as.data.frame(correct_data("raw", features, meta))
batch_corrected_data$lm <- as.data.frame(correct_data("lm", features, meta))
batch_corrected_data$zscore <- as.data.frame(correct_data("zscore", features, meta))
batch_corrected_data$limma <- as.data.frame(correct_data("limma", features, meta))
batch_corrected_data$combat <- as.data.frame(correct_data("combat", features, meta))
batch_confound <- c("B", "K", "M", "P", "T", "S")
idx_confound <- which(!meta$contributor_code %in% batch_confound)
batch_corrected_data$combat_prec <- as.data.frame(correct_data("combat", features[idx_confound, ], meta[idx_confound,]))
gc()


metrics_list <- vector("list", length(batch_corrected_data))
names(metrics_list) <- names(batch_corrected_data)


library(ranger)
# Loop over each correction method
for (m in names(batch_corrected_data)) {
  cat("Now working on:", m, "\n")
  if (m == "combat_prec") {
    md  <- meta[idx_confound, ]
  } else {
    md  <- meta
  }
  mat <- batch_corrected_data[[m]]
  # Compute metrics
  cat("Now working on:", m, "sv", "\n")
  sv   <- compute_site_variance_fast(as.matrix(mat), md$contributor_code)
  cat("Now working on:", m, "pca", "\n")
  pc_vars  <- compute_pca_var_top2(mat)
  cat("Now working on:", m, "ml", "\n")
  auc  <- compute_outcome_auc(mat, md$contributor_code)
  cat("Now working on:", m, "sil", "\n")
  sil  <- compute_silhouette_fast(mat, md$contributor_code)
  
  # Store in a one-row data.frame
  metrics_list[[m]] <- data.frame(
    Method = m,
    Median_Site_Variance = sv,
    PC1_Variance <- pc_vars["PC1"],
    PC1_Variance <- pc_vars["PC2"],
    Outcome_AUC = auc,
    Silhouette_Score = sil,
    stringsAsFactors     = FALSE
  )
}

# Combine all rows into a single data.frame
evaluation_results <- do.call(rbind, metrics_list)
saveRDS(evaluation_results, file = gzfile(paste0(here(), "/output/tables/eval_ad_ctr.rds.gz")))



# Variance
methods <- names(batch_corrected_data)
results_var <- data.frame(
  Method            = methods,
  site_variance     = NA_real_,
  outcome_variance  = NA_real_,
  stringsAsFactors  = FALSE
)

vp_list <- list()
for (i in seq_along(methods)) {
  m   <- methods[i]
  mat <- batch_corrected_data[[m]]
  md  <- meta
  
  # Transpose: samples = rows
  expr <- as.data.frame(mat)
  expr$sample_id <- rownames(expr)
  
  # Attach sample_id to metadata for alignment
  md2 <- md %>% select(sample_id, outcome, contributor_code) 
  df_vp <- variancePartition::fitExtractVarPartModel(
    t(mat),            # genes × samples
    ~ (1|outcome) + (1|contributor_code),
    md2
  )
  f_vp$Method <- m
  vp_list[[m]] <- df_vp
  
  # Summarize
  results_var$site_variance[i]    <- median(df_vp$contributor_code, na.rm = TRUE)
  results_var$outcome_variance[i] <- median(df_vp$outcome,         na.rm = TRUE)
}

saveRDS(vp_list, file = gzfile(paste0(here(), "output/tables/vp_list.rds.gz")))





