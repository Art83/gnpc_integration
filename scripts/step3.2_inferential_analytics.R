# Final data analysis for AD, for other dosorders the script is the same

library(dplyr)
library(here)
source(paste0(here(),"/scripts/utils.R"))

merged_proteomics_subset <- readRDS(gzcon(file("output/proc/dpgen.rds.gz", open = "rb"))) # datatset is not in github due to user agreement restrictions 

ignore_vars <- c("contributor_code","sample_id", "sex", "race", "plate_id", 
                 "outcome", "non_na_protein_count", "inferred_assay")

proteomics_7k <- merged_proteomics_subset %>%
  filter(inferred_assay == "7K") %>%
  filter(outcome %in% c("AD", "HC"))


df_subset <- proteomics_7k


threshold <- 0.1
df_subset  <- df_subset %>%
  mutate(across(!all_of(ignore_vars), ~ na_if(., -1))) %>%
  {
    protein_cols <- setdiff(names(.), ignore_vars)
    na_prop <- colMeans(is.na(select(., all_of(protein_cols))))
    good_cols <- names(na_prop[na_prop < threshold])
    select(., all_of(c(good_cols, ignore_vars))) %>%
      group_by(outcome) %>%
      mutate(across(all_of(good_cols), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))
  }

df_subset[, !colnames(df_subset) %in% ignore_vars] <-
  log2(df_subset[, !colnames(df_subset) %in% ignore_vars] + 1)



matched_AD <- readRDS("/proc/output/AD.rds.gz") # datatset is not in github due to user agreement restrictions 
ad_cases <- matched_AD$matched_cases$sample_id
hc_cases <- matched_AD$matched_controls$sample_id

# Separate features and metadata
features <- df_subset[df_subset$sample_id %in% c(ad_cases, hc_cases), !colnames(df_subset) %in% ignore_vars]
meta <- df_subset[df_subset$sample_id %in% c(ad_cases, hc_cases),colnames(df_subset) %in% ignore_vars]
outcome <- as.factor(meta$outcome)
meta$outcome <- relevel(factor(meta$outcome), ref = "HC")
meta$batch <- factor(paste(meta$contributor_code, meta$plate_id, sep="_"))


batch_corrected_data <- list()
batch_corrected_data$raw <- as.data.frame(correct_data("raw", features, meta))
batch_corrected_data$limma <- as.data.frame(correct_data("limma", features, meta))

sites <- LETTERS[1:18]
mx_sites <- sites[!sites %in% c("B", "K", "M", "P", "O", "H")]
idx_mixed_sites <- which(meta$contributor_code %in% mx_sites)

# 1) site-wise limma
list_res <- list()
for (i in seq_along(mx_sites) ){
  s = mx_sites[i]
  idx <- meta$outcome %in% c("HC","AD") & meta$contributor_code == s
  X   <- t(batch_corrected_data[["raw"]][idx,])
  m   <- meta[idx, ]
  design <- model.matrix(~ outcome + batch, m)
  fit    <- limma::lmFit(X, design) |> limma::eBayes()
  tt     <- limma::topTable(fit, coef = "outcomeAD", n = Inf) |>
    tibble::rownames_to_column("protein") |>
    mutate(site = s,
           se   = abs(logFC / t))           # limma SE from t-stat
  list_res[[s]] <- tt[, c("protein","site","logFC","se", "t")]
}

site_tbl <- bind_rows(list_res)

# 2) inverse-variance meta per protein
meta_df <- site_tbl |>
  group_by(protein) |>
  group_modify(~{
    yi <- .x$logFC; sei <- .x$se
    fe <- metafor::rma.uni(yi = yi, sei = sei, method = "FE")
    data.frame(logFC = fe$b[1], se_meta = fe$se, p_meta = fe$pval,
               k = nrow(.x), I2 = fe$I2)
  }) |>
  ungroup() |>
  mutate(P.Value = p_meta,
         adj.P.Val = p.adjust(p_meta, "BH")) %>%
  arrange(adj.P.Val) %>%
  select(protein, logFC, P.Value, adj.P.Val) %>%
  mutate(method = "meta")


#pooled limma
design <- model.matrix(~ outcome+plate_id, data = meta)
  
# Run limma
fit <- limma::lmFit(t(batch_corrected_data[["limma"]]), design) |> limma::eBayes()
  
pooled <- limma::topTable(fit, coef = "outcomeAD", n = Inf) %>%
  tibble::rownames_to_column("protein") %>%
  mutate(method = "limma") %>%
  select(protein, logFC_pooled = logFC,
         P.Value_pooled = P.Value,
         adjP_pooled = adj.P.Val)

site_tbl <- site_tbl %>%
  transmute(
    protein = as.character(protein),
    site    = as.character(site),
    logFC   = as.numeric(logFC),
    t       = as.numeric(t),
    se      = abs(logFC / t)
  ) %>%
  filter(is.finite(logFC), is.finite(se), se > 0)

pooled_tbl <- pooled %>%
  transmute(
    protein      = as.character(protein),
    logFC_pooled = as.numeric(logFC_pooled),
    P_pooled     = as.numeric(P.Value_pooled),
    adjP_pooled  = as.numeric(adjP_pooled)
  ) %>%
  distinct(protein, .keep_all = TRUE)


## ---------------------------------------------------
## 1) Per-protein meta diagnostics + robust flags
##    - FE meta (primary estimate)
##    - RE prediction interval for direction guard
##    - CI-based directional discord (neutral sites allowed)
##    - Weighted discord relative to non-neutral weight
##    - LOO fragility only if FE meta is significant
## ---------------------------------------------------

meta_diag <- site_tbl %>%
  group_by(protein) %>%
  group_modify(~ meta_diag_one(.x)) %>%
  ungroup() %>%
  mutate(adjP_meta = p.adjust(p_meta, "BH"))

## --------------------------------------------
## 2) Gate candidates & assign final tiers
## --------------------------------------------
opp_cut   <- 0.30  # 30% of non-neutral weight on the opposite side = discord
dom_cut   <- 0.60  # single site dominates
het_low   <- 30    # for Supported (underpowered)
require_same_sign_for_core <- TRUE

calls <- pooled_tbl %>%
  inner_join(meta_diag, by = "protein") %>%
  mutate(
    is_candidate = (adjP_pooled < 0.05 | adjP_meta < 0.10),
    same_sign    = sign(logFC_pooled) == sign(logFC_meta),
    RE_PI_zero   = mapply(pi_covers_zero, re_pi_lb, re_pi_ub),
    high_conf = (adjP_meta < 0.05 & adjP_pooled < 0.05) &
      (opp_wt_rel < opp_cut) &
      !lo_flip &
      !(wmax > dom_cut & lo_p_cross) &
      (if (require_same_sign_for_core) same_sign else TRUE),
    ## Supported (underpowered): pooled FDR, coherent direction, nominal meta p
    supported = (adjP_pooled < 0.05) & same_sign & (p_meta < 0.05) & (I2 <= het_low) &
      (opp_wt_rel < opp_cut) & !RE_PI_zero,
    
    ## Likely site-driven: MUST have directional problems (weighted discord) OR fragility.
    ## RE PI crossing 0 alone is NOT enough to exclude.
    site_driven = (opp_wt_rel >= opp_cut) | lo_flip | (wmax > dom_cut & lo_p_cross),
    
    tier = dplyr::case_when(
      is_candidate & high_conf   ~ "High-confidence (meta, concordant)",
      is_candidate & supported   ~ "Supported (underpowered)",
      is_candidate & site_driven ~ "Likely site-driven",
      is_candidate               ~ "Ambiguous",
      TRUE                       ~ "Background (not evaluated)"
    ),
    
    ## Optional annotation columns for reporting, not filtering:
    note_hetero = dplyr::case_when(
      I2 > 50 & !site_driven ~ "high I2 (magnitude)",
      RE_PI_zero & !site_driven ~ "RE PI spans 0 (magnitude)",
      TRUE ~ NA_character_
    )
  )

tier_counts <- calls %>% count(tier, name = "n") %>% arrange(desc(n))
print(tier_counts)

write.csv(calls, paste0(here(), "/output/tables/dt_res_AD.csv"), row.names = FALSE)