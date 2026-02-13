library(dplyr)
library(stringr)

# Base data directory
data_dir <- "GNPC_HDS_general"

# Columns to keep as metadata
ignore_vars <- c(
  "sample_id", "contributor_code", "sex", "race",
  "age_at_visit", "plate_id", "outcome",
  "non_na_protein_count", "inferred_assay"
)

# ===============================
# 1. Main preprocessing pipeline
# ===============================

# 2.1 Load metadata and annotation
clinical_df <- load_gnpc_file("HDS V1.3MS APOE/ClinicalV1_3ms_apoe.csv")
mapping_df  <- load_gnpc_file("HDS V1.3MS APOE/PersonMappingV1_3ms_apoe.csv")
sample_meta <- load_gnpc_file("HDS V1.3MS APOE/SomalogicMetaV1_3ms_apoe.csv")
annotation_uniprot <- load_gnpc_file("HDS V1.3MS APOE/SomalogicAnalyteInfoV1_3ms_apoe.csv")

# 2.2 Define outcomes
df_outcome <- define_outcomes(clinical_df, mapping_df, sample_meta)

# 2.3 Determine plasma sample IDs
plasma_ids <- get_plasma_ids(sample_meta)

# 2.4 Load & merge all 12 plates
plate_list <- lapply(1:12, load_plate, valid_ids = plasma_ids)
merged_proteomics <- Reduce(
  function(x, y) inner_join(x, y, by = "sample_id"),
  plate_list
)

# 2.5 Join clinical outcomes & metadata
df <- merged_proteomics %>%
  inner_join(df_outcome, by = "sample_id")

# 2.6 Compute non-missing aptamer counts & infer assay type
protein_cols <- setdiff(colnames(df), ignore_vars)
df$non_na_protein_count <- rowSums(df[protein_cols] != -1, na.rm = TRUE)
df$inferred_assay <- case_when(
  df$non_na_protein_count >= 7500                        ~ "7K",
  df$non_na_protein_count >= 5000                        ~ "5K+",
  df$non_na_protein_count >= 4000                        ~ "5K",
  df$non_na_protein_count < 2000                         ~ "Restricted/Other",
  TRUE                                                   ~ "Unknown"
)


# 2.7 Annotate aptamers with UniProt IDs
df <- annotate_aptamers(df, annotation_uniprot)

# 2.8 Remove any accidental empty/duplicated columns
df <- df[,!is.na(colnames(df)) &
           !duplicated(colnames(df)) &
           colnames(df) != ""]

# 2.9 Save for downstream analyses
saveRDS(df7K, gzfile("AS_CF/proc/dpgen.rds.gz"))


