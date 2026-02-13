library(dplyr)
library(here)
source(paste0(here(),"/scripts/utils.R"))


merged_proteomics_subset <- readRDS(gzcon(file("output/proc/dpgen.rds.gz", open = "rb"))) # dataset is not in github due to user agreement restrictions 


ignore_vars <- c("contributor_code","sample_id", "sex", "race", "plate_id", 
                 "outcome", "non_na_protein_count", "inferred_assay")


proteomics_7k <- merged_proteomics_subset %>%
  filter(inferred_assay == "7K") %>%
  filter(!(contributor_code == "K" & outcome == "HC" | (contributor_code %in% c("N","S") & outcome == "ALS")) )
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

clinical_df <- read.csv("output/proc/ClinicalV1_3ms_apoe.csv") # dataset is not in github due to user agreement restrictions 

df <- clinical_df %>%
  select(sample_id, age_at_visit, sex) %>%
  inner_join(proteomics_7k[,c("sample_id", "outcome", "contributor_code")])


# Ensure correct factor levels
df <- df %>% 
  mutate(outcome = factor(outcome))

# Pairwise comparison: each disorder vs Control
disorders <- levels(df$outcome)[levels(df$outcome) != "HC"]

age_comparisons <- lapply(disorders, function(dis) {
  if(dis=="ALS"){
    tmp <- df %>% filter(contributor_code=="M", 
                         outcome %in% c("HC", dis))
    tmp$outcome <- factor(tmp$outcome, levels=c("HC", dis))
  } else {
    tmp <- df %>% filter(outcome %in% c("HC", dis))
    tmp$outcome <- factor(tmp$outcome, levels=c("HC", dis))
  }
  
  eff <- rstatix::wilcox_effsize(tmp, age_at_visit ~ outcome, ref.group = "HC", ci = TRUE)
  eff$comparison <- paste(dis, "vs Control")
  eff
}) %>% bind_rows()





matched_AD <- downsample_HC_to_match_age_sex(df, group_case = "AD")
matched_HC_AD <- matched_AD$matched_HC
matched_case_AD <- matched_AD$case_group
plot_age_sex_balance(df, matched_AD, "AD")


matched_ALS <- downsample_HC_to_match_age_sex(df[df$contributor_code=="M",], group_case = "HC", group_control = "ALS")
matched_HC_ALS <- matched_ALS$matched_HC
matched_case_ALS <- matched_ALS$case_group
plot_age_sex_balance(df, matched_ALS, "ALS")

matched_PDD <- downsample_HC_to_match_age_sex(df, group_case = "PDD")
matched_HC_PDD <- matched_PDD$matched_HC
matched_case_PDD <- matched_PDD$case_group
plot_age_sex_balance(df, matched_PDD, "PDD")

matched_FTD <- downsample_HC_to_match_age_sex(df, group_case = "FTD")
matched_HC_FTD <- matched_FTD$matched_HC
matched_case_FTD <- matched_FTD$case_group
plot_age_sex_balance(df, matched_FTD, "FTD")

matched_PD <- downsample_HC_to_match_age_sex(df, group_case = "PD")
matched_HC_PD <- matched_PD$matched_HC
matched_case_PD <- matched_PD$case_group
plot_age_sex_balance(df, matched_PD, "PD")



matchedAD <- match_controls_by_age_sex(meta_df = df, group_case = "AD")
matchedFTD <- match_controls_by_age_sex(meta_df = df, group_case = "FTD")
matchedPD <- match_controls_by_age_sex(meta_df = df, group_case = "PD")
matchedALS <- match_controls_by_age_sex(meta_df = df, group_case = "ALS")
matchedPDD <- match_controls_by_age_sex(meta_df = df, group_case = "PDD")


















plot_age_distributions <- function(case, ctrl, case_label = "AD") {
  ggplot(bind_rows(
    mutate(case, group = case_label),
    mutate(ctrl, group = "HC (matched)")
  ), aes(x = age_at_visit, fill = group)) +
    geom_density(alpha = 0.5) +
    theme_minimal() +
    ggtitle(paste("Age distribution matching:", case_label, "vs HC")) +
    xlab("Age") + ylab("Density")
}

plot_age_distributions(matched_case_AD, matched_HC_AD, "AD")
plot_age_distributions(df[df$outcome=="AD",], df[df$outcome=="HC",], "AD")



plot_age_distributions(matched_case_ALS, matched_HC_ALS, "ALS")
plot_age_distributions(df[df$outcome=="ALS",], df[df$outcome=="HC",], "ALS")


plot_age_distributions(matched_case_PDD, matched_HC_PDD, "PDD")
plot_age_distributions(df[df$outcome=="PDD",], df[df$outcome=="HC",], "PDD")

saveRDS(matched_AD, file=gzfile("~/files/AS_CF/proc/be_man/proc/matched_AD.rds.gz"))
saveRDS(matched_PDD, file=gzfile("~/files/AS_CF/proc/be_man/proc/matched_PDD.rds.gz"))
saveRDS(matched_FTD, file=gzfile("~/files/AS_CF/proc/be_man/proc/matched_FTD.rds.gz"))
saveRDS(matched_PD, file=gzfile("~/files/AS_CF/proc/be_man/proc/matched_PD.rds.gz"))
saveRDS(matched_ALS, file=gzfile("~/files/AS_CF/proc/be_man/proc/matched_ALS.rds.gz"))
