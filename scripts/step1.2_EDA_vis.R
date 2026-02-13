library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(here)

# 1. evaluation_table: assume summary_table has columns
# Method, Median_Site_Variance, site_variance, outcome_variance,
# PC1_variance, PC2_variance, Silhouette_Score

summary_table <- readRDS(gzcon(file(paste0(here(), "/output/tables/eval_ad_ctr.rds.gz"), open = "rb")))
vp_list <- readRDS(gzfile(paste0(here(), "/output/tables/vp_list.rds.gz")))
vp_long <- do.call(rbind, lapply(names(vp_list), function(m) {
  df  <- as.data.frame(vp_list[[m]], stringsAsFactors=FALSE)
  df$Method <- m
  df
}))

summary_table <- summary_table %>%
  mutate(Method = ifelse(Method == "combat_prec", "combat(bal)", Method) ) %>%
  mutate(Method = factor(Method, levels = c("raw", "zscore", "lm", "limma", "combat", "combat(bal)")))

vp_long <- vp_long %>%
  mutate(Method = ifelse(Method == "combat_prec", "combat(bal)", Method) ) %>%
  mutate(Method = factor(Method, levels = c("raw", "zscore", "lm", "limma", "combat", "combat(bal)")))

###### General Methods ######
# Figure 1. PCA (PC1, PC2)
pc_colors <- c("PC1" = "grey40", "PC2" = "grey70")

method_colors <- c(
  "raw" = "#7F7F7F",
  "zscore" = "#1F78B4",
  "lm" = "#FF7F00",
  "limma" = "#33A02C",
  "combat" = "#6A3D9A",
  "combat(bal)" = "#E31A1C"
)
df3 <- summary_table %>%
  select(Method, PC1_variance, PC2_variance) %>%
  tidyr::pivot_longer(-Method, names_to="PC", values_to="Variance") %>%
  mutate(PC = recode(PC,
                     PC1_variance = "PC1 Variance (%)",
                     PC2_variance = "PC2 Variance (%)"
  ))

tiff(paste0(here(), "/output/pics/PCs.tiff"), width = 6, height = 4, units = "in", res = 300)
ggplot(df3, aes(x = Method, y = Variance, fill = Method)) +
  geom_col() +
  facet_wrap(~PC) +
  scale_fill_manual(values = method_colors)+
  theme_minimal() +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14, angle = 45, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16))
dev.off()

ggplot(df3, aes(x = Method, y = Variance, fill = PC)) +
  geom_col(position = "dodge") +
  labs(
    title = "Variance Explained by PC1 and PC2",
    x = "Method", y = "Variance (%)", fill = ""
  ) +
  theme_minimal()

# Figure 2 Silhouette Scores
tiff(paste0(here(), "/output/pics/SC.tiff"), width = 6, height = 4, units = "in", res = 300)
summary_table %>%
  ggplot(aes(x = reorder(Method, Silhouette_Score), y = Silhouette_Score, fill = Method)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = method_colors)+
  labs(x = "Method", y = "Silhouette Score")+
  theme_minimal()+
    theme(axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          legend.text = element_text(size=14),
          legend.title = element_text(size=16))
dev.off()

############Retention of Biological Signal #################
raw_outcome <- summary_table$outcome_variance[summary_table$Method=="raw"]

retention_df <- summary_table %>%
  mutate(
    Retention = outcome_variance / raw_outcome * 100
  )

# Figure 3 Retention of outcome
tiff(paste0(here(), "/output/pics/retention_bio.tiff"), width = 6, height = 4, units = "in", res = 300)
ggplot(retention_df, aes(x = reorder(Method, -Retention), y = Retention, fill = Method)) +
  geom_col() +
  coord_flip(ylim = c(0,100)) +
  labs(
    x     = "Method",
    y     = "Median Outcome Variance Retained (%)"
  ) +
  scale_fill_manual(values = method_colors)+
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16))
dev.off()

# Figure 4. AUC of ranger in predicting AD vs Ctr
tiff(paste0(here(), "/output/pics/auc.tiff"), width = 6, height = 4, units = "in", res = 300)
ggplot(summary_table, aes(x = reorder(Method, -Outcome_AUC), y = Outcome_AUC, fill = Method)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "AD vs HC Classification AUC by Correction Method",
    x     = "Method",
    y     = "AUC"
  ) +
  theme_minimal() +
  theme(legend.position="none")
dev.off()

############Variance#################
library(scales)

# define the breaks
my_breaks <- c(1e-10, 1e-4, 1e-3, 1e-2, 1e-1, 1, 5, 10)

# and human‐friendly labels
my_labels <- c("0", "0.0001%", "0.001%", "0.01%", "0.1%", "1%", "5%", "10%")

vp_tidy <- vp_long %>%
  select(Method, site = contributor_code, outcome) %>%
  tidyr::pivot_longer(
    cols = c(site, outcome),
    names_to  = "Component",
    values_to = "Variance"
  ) %>%
  mutate(Variance = Variance * 100)


raw_med <- median(vp_tidy$Variance[vp_tidy$Method == "raw" & vp_tidy$Component=="outcome" ], na.rm = TRUE)
tiff(paste0(here(), "/output/pics/outcome_sim.tiff"), width = 6, height = 5, units = "in", res = 300)
ggplot(vp_tidy[vp_tidy$Component=="outcome",], aes(x = Method, y = Variance, fill = Method)) +
  geom_violin(trim = TRUE, scale = "width", na.rm = TRUE) +
  stat_summary(fun = median, geom = "point", color = "white", size = 2, na.rm = TRUE) +
  scale_y_continuous(
    trans  = pseudo_log_trans(sigma = 1e-4, base = 10),
    breaks = my_breaks,
    labels = my_labels
  ) +
  geom_hline(                                 
    yintercept = raw_med,
    linetype   = "dashed",
    color      = "black",
    size       = 0.5
  )+
  labs(
    x     = "Method",
    y     = "Variance (%) (pseudo-log scale)"
  ) +
  scale_fill_manual(values = method_colors)+
  theme_minimal() +
  theme(legend.position = "none") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    legend.position  = "none"
  )+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16))

dev.off()


tiff(paste0(here(), "/output/pics/site_sim.tiff"), width = 6, height = 6, units = "in", res = 300)
ggplot(vp_tidy[vp_tidy$Component=="site",], aes(x = Method, y = Variance, fill = Method)) +
  geom_violin(trim = TRUE, scale = "width", na.rm = TRUE) +
  stat_summary(fun = median, geom = "point", color = "white", size = 2, na.rm = TRUE) +
  scale_y_continuous(
    trans  = pseudo_log_trans(sigma = 1e-4, base = 10),
    breaks = my_breaks,
    labels = my_labels
  ) +
  labs(
    title = "Per-Protein VariancePartition Distributions",
    x     = "Method",
    y     = "Variance (%) (pseudo-log scale)"
  ) +
  scale_fill_manual(values = method_colors)+
  theme_minimal() +
  theme(legend.position = "none") +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey95"),
    legend.position  = "none"
  )+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16))
dev.off()



