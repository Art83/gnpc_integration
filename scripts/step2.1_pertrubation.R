library(dplyr)

merged_proteomics_subset <- readRDS(gzcon(file("output/proc/dpgen.rds.gz", open = "rb"))) # datatset is not in github due to user agreement restrictions 


ignore_vars <- c("contributor_code","sample_id", "sex", "race", "plate_id", 
                 "outcome", "non_na_protein_count", "inferred_assay")


# ---------- Step 2: Filter to 7K only ----------
proteomics_7k <- merged_proteomics_subset %>%
  filter(inferred_assay == "7K") %>%
  filter(outcome %in% c("AD", "HC"))


rm(merged_proteomics_subset)

threshold <- 0.1
df_subset  <- proteomics_7k %>%
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
df_subset[, !colnames(df_subset) %in% ignore_vars] <-
  log2(df_subset[, !colnames(df_subset) %in% ignore_vars] + 1) %>%
  ungroup()

rm(proteomics_7k)
gc()

N_SAMPLES <- 2000
N_PROTEINS <- 1000
REPETITIONS <- 100


effect_sizes <- c(0, 0.5, 1.0, 2.0)
confounding_levels <- c("none", "moderate", "severe")
heteroskedastic_flags <- c(FALSE, TRUE)
nonlinear_batch_flags <- c(FALSE, TRUE)
signal_sparsity_levels <- c(0, 0.05, 0.1) # 0 = null signal





grid <- expand.grid(
  rep            = 1:REPETITIONS,
  effect_size    = effect_sizes,
  confounding    = confounding_levels,
  heteroskedastic = heteroskedastic_flags,
  nonlinear      = nonlinear_batch_flags,
  signal_sparsity= signal_sparsity_levels,
  stringsAsFactors = FALSE
)


meta <- df_subset[, c("sample_id","contributor_code","outcome","plate_id","sex","race", "non_na_protein_count", "inferred_assay")]

future::plan(future::multicore,   workers = parallel::detectCores() - 2)
options(future.globals.maxSize = 2 * 1024^3)

results_list <- future.apply::future_lapply(
    seq_len(nrow(grid)),           
    function(i) run_scenario(grid[i, ]),
    future.seed=TRUE
  )
  grid_rows <- split(grid, seq_len(nrow(grid)))
  results_with_meta <- Map(
    function(res_dt, meta_row) {
      cbind(data.table::as.data.table(meta_row), res_dt)
    },
    results_list,
    grid_rows
  )
final_results <- data.table::rbindlist(results_with_meta, use.names = TRUE, fill = TRUE)
saveRDS(final_results, gzfile(paste0(here(), "/output/tables/simulation_results.rds.gz")))
