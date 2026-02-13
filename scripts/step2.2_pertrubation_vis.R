library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

sim_result <- readRDS(paste0(here(),"/output/tables/simulation_results.rds"))

final_results <- sim_result %>%
  group_by(rep, effect_size, confounding, heteroskedastic, 
           nonlinear, signal_sparsity) %>%
  mutate(scenario_id = cur_group_id()) %>%
  ungroup() %>%
  mutate(signal_fraction = recovery_rate * (varpart_outcome / (varpart_outcome + site_var + varpart_site + 1e-30)))





m_site <- lme4::lmer(site_var ~ method
                     + effect_size 
                     + confounding 
                     + heteroskedastic 
                     + nonlinear 
                     + signal_sparsity
                     + (1|scenario_id),
                     data = final_results)





m_out  <- lme4::lmer(varpart_outcome ~ method
                     + effect_size 
                     + confounding 
                     + heteroskedastic 
                     + nonlinear 
                     + signal_sparsity
                     + (1|scenario_id),
                     data = final_results)


emm_conf <- emmeans::emmeans(m_out, ~ method | confounding, pbkrtest.limit = 76553)
emm_hetero <- emmeans::emmeans(m_out, ~ method |  heteroskedastic, pbkrtest.limit = 76553)
emm_nonlinear <- emmeans::emmeans(m_out, ~ method | nonlinear, pbkrtest.limit = 76553)

emm_conf <- emmeans::emmeans(m_site, ~ method | confounding, pbkrtest.limit = 76553)
emm_hetero <- emmeans::emmeans(m_site, ~ method |  heteroskedastic, pbkrtest.limit = 76553)
emm_nonlinear <- emmeans::emmeans(m_site, ~ method | nonlinear, pbkrtest.limit = 76553)




emm_conf2 <- as.data.frame(emm_conf) %>%
  rename(level = confounding) %>%
  mutate(grouping = "Confounding",
         level=as.character(level))

emm_hetero2 <- as.data.frame(emm_hetero) %>%
  rename(level = heteroskedastic) %>%
  mutate(grouping = "Heteroskedastic",
         level=as.character(level))

emm_nonlin2 <- as.data.frame(emm_nonlinear) %>%
  rename(level = nonlinear) %>%
  mutate(grouping = "Nonlinear",
         level=as.character(level))

emm_all <- bind_rows(emm_conf2, emm_hetero2, emm_nonlin2) %>%
  filter(!is.na(emmean)) %>%           # drop the non-estimable rows
  mutate(
    grouping = factor(grouping,      # order your strips top→bottom if you like
                      c("Confounding","Heteroskedastic","Nonlinear")),
    level    = as.character(level)   # so ggplot doesn’t try to re-combine levels
  )


emm_all <- emm_all %>%
  mutate(
    level = if_else(
      grouping == "Confounding",
      forcats::fct_relevel(level, "none", "moderate", "severe"),
      level
    )
  ) %>%
  mutate(method = ifelse(as.character(method) == "combat_prec", "combat(bal)", as.character(method)))

tiff(paste0(here(), "/output/tables/emm_site.tiff"), width = 7, height = 4, units = "in", res = 300)

emm_all %>%
  mutate(
    x_disp = case_when(
      grouping == "Confounding" ~ factor(level,
                                         levels = c("none", "moderate", "severe")
      ),
      TRUE ~ factor(level)  # default factor for other facets ("FALSE","TRUE")
    )
  ) %>%
  ggplot(aes(x = x_disp, y = emmean, color = method)) +
  geom_pointrange(
    aes(ymin = lower.CL, ymax = upper.CL),
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~ grouping, scales = "free_x", nrow = 1) +
  #scale_y_continuous(limits = c(0.01, 0.1)) +
  scale_color_manual(values = method_colors) +
  labs(
    x = NULL,
    y = "EMM site",
    color = "Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    legend.position  = "bottom"
  )

dev.off()




