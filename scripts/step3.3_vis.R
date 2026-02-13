library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(here)


ad_df <- read.csv(paste0(here(), "/output/tables/dt_res_AD.csv"))

volc_df <- ad_df %>%
  mutate(neglog10 = -log10(pmax(p_meta, 1e-300)),
         tier = case_when(tier == "High-confidence (meta, concordant)" ~ "High-confidence",
                          tier == "Likely site-driven" ~ "Site-driven",
                          tier == "Supported (underpowered)" ~ "Underpowered",
                          tier == "Ambiguous" ~ "Ambiguous",
                          tier == "Background (not evaluated)" ~ "Background"),
        tier = factor(tier, level = c("High-confidence",
                                      "Site-driven",
                                      "Ambiguous",
                                      "Underpowered",
                                      "Background") )  ) %>%
  rename(Tier = tier)


tiff(paste0(here(), "/output/pics/volcano_good.tiff"), width = 6, height = 5, units = "in", res = 300)
ggplot(volc_df, aes(x = logFC_meta, y = neglog10, color = Tier) ) + 
  geom_point(alpha = 0.6, size = 1) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_colour_manual(
    values = c(
      "High-confidence" = "red", 
      "Site-driven" = "black",  
      "Underpowered"  = "blue",
      "Ambiguous" = "purple",
      "Background" = "grey"
    )
  ) +
  labs(x = "Meta logFC", y = "-log10(p-value)") + 
  theme_classic()+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        legend.text = element_text(size=14),
        legend.title = element_text(size=16))
dev.off()

tiff(paste0(here(), "/output/pics/volcano_bad.tiff"), width = 6, height = 5, units = "in", res = 300)
ggplot(volc_df[volc_df$Tier != "High-confidence",], aes(x = logFC_meta, y = neglog10, color = Tier) ) + 
  geom_point(alpha = 0.6, size = 1) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_colour_manual(
    values = c(
      "Site-driven" = "black",  
      "Underpowered"= "blue",
      "Ambiguous" = "purple",
      "Background" = "grey"
    )) +
  labs(x = "Meta logFC", y = "-log10(Meta p-value)") + 
  theme_classic()+
  theme(axis.title.x = element_text(size=16),
                       axis.title.y = element_text(size=16),
                       axis.text.x = element_text(size=14),
                       axis.text.y = element_text(size=14),
                       legend.text = element_text(size=14),
                       legend.title = element_text(size=16))
dev.off()
