### Libraries
library(dplyr)
library(readr)
library(readxl)
library(writexl)
library(ggplot2)
library(tibble)
library(ggrepel)
library(tidyr)
library(cowplot)
library(palmerpenguins)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(skimr)
library(ggExtra)
library(RColorBrewer)
library(gridExtra)
library(Matrix)
library(lme4)
library(lmerTest)
library(pillar)
library(DT)
library(forcats)
library(emmeans)
library(openxlsx)
library(FSA)
library(car)
library(coin)
library(glmmTMB)
library(ggeffects)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(svglite)

# ###  Statistic and Visualization p75ECD

# Raw data, without SYMP and NA, just the 4 groups to analyze 
Urine_data <-  read_excel("results/clinico_molecular_without_outliers_final.xlsx") 

#### Adjusting for sex, age and center using LM which will "normalize" the residuals.
#### Typically, after adjustment it is recommended to use ANCOVA, which is an ANOVA with covariates
#### However, I checked the normality of the data of the residuals of the LM, and the data is not normal..

####### UPDATE:
#### Ana (24/04): i was searching better and perhaps there is a better way to deal with the covariate adjustment 
#### not using ANCOVA, since the data is not normal. The problem is that with KW the interpretation becomes
#### harder as we are comparing median ranks and not really means of the data.
#### After searching a lot, I have understood that we can use a linear model based on Gamma distribution
#### this allows for non-normal data, and can be applied for data which is always positive and skewed, 
#### which is our case, as the urine data is always positive and skewed to the right distribution
#### I hava modified the code to use this instead 


### p75ECD (adjustment for age, sex and country)


### STEP 1: Violin plot and statistic to compare ALS, PGMC, CTR and mimics


# Function to adjust for age, sex, and country, here include type, for adjusted values to plot.
adjust_for_covariates <- function(data, visit_name) {
  data_visit <- data %>% filter(Visit == visit_name)
  
  model <- lm(`ng p75ECD/mg creatinine` ~ type + age + sex + Country, data = data_visit)
  
  data_visit <- data_visit %>%
    mutate(score_adjusted = fitted(model))
  
  list(model = model, data = data_visit)
}


# Function to plot
# Urine_data = merged_table_final_group # Not used? can we remove this line?

plot_data <- Urine_data %>%
  mutate(
    `ng p75ECD/mg creatinine` = as.numeric(`ng p75ECD/mg creatinine`),
    type = factor(type,levels = c("ALS", "PGMC", "CTR", "mimic")),
    sex = factor(sex),
    age = as.numeric(age),
    Country = factor(Country)  
  )

# Residuals, new adjusted values used for plotting and statistic corrected for sex age and country  

adjusted_V0 <- adjust_for_covariates(plot_data, "V0")
adjusted_V1 <- adjust_for_covariates(plot_data, "V1")

# check normality of residuals, separatly for V0 and V1? As here we will do two plots. 
shapiro.test(residuals(adjusted_V0$model))
# p-value = 3.054e-09 -> not normaly distributed
qqnorm(residuals(adjusted_V0$model)); qqline(residuals(adjusted_V0$model))
# verification data is right skewed
ggplot(adjusted_V0$data, aes(x = `ng p75ECD/mg creatinine`)) +
  geom_histogram(bins = 30)

shapiro.test(residuals(adjusted_V1$model))
# p-value = 2.097e-06 -> not normaly distributed
qqnorm(residuals(adjusted_V1$model)); qqline(residuals(adjusted_V1$model))
# verification data is right skewed
ggplot(adjusted_V1$data, aes(x = `ng p75ECD/mg creatinine`)) +
  geom_histogram(bins = 30)

### Summary: not normalised data based on residuals, and data always positive and skewed to right
### possible: to use a glm model with gamma distribution

run_gamma_analysis <- function(data, visit_name) {
  
  data_visit <- data %>% filter(Visit == visit_name)
  
  # fit a Gamma GLM (non normaly distributed data, lmr is for normaly distributed data)
  model_gamma <- glm(`ng p75ECD/mg creatinine` ~ type + age + sex + Country, 
                     family = Gamma(link = "log"), 
                     data = data_visit)
  
  # stats per variable
  stats_table <- car::Anova(model_gamma, type = "III")
  
  # stats overall
  anova_pvalue = as.data.frame(stats_table) 
  anova_pvalue = anova_pvalue$`Pr(>Chisq)`[1]
  
  # post-hoc with type == "response" to transform log data to original units
  adj_comparisons <- emmeans(model_gamma, pairwise ~ type, type = "response")
  
  # get means adjusted for covariates
  adj_means_df <- as.data.frame(adj_comparisons$emmeans)
  
  # remove effect of covariates for the visualisation
  model_cov = glm(`ng p75ECD/mg creatinine` ~ age + sex + Country,
                  family = Gamma(link="log"),
                  data = data_visit)
  predict_cov = predict(model_cov,type = "link")
  
  data_visit$score_adj = exp(log(data_visit$`ng p75ECD/mg creatinine`) - predict_cov + mean(predict_cov))
  
  list(model = model_gamma,
    anova = stats_table,
    anova_pvalue = anova_pvalue,
    comparisons = adj_comparisons$contrasts,
    plot_data = adj_means_df,
    data_adj = data_visit)
}

# get results for each visit
results_V0 <- run_gamma_analysis(plot_data, "V0")
print(results_V0$anova)
print(results_V0$comparisons)

# contrast     ratio     SE  df null t.ratio p.value
# ALS / PGMC   1.384 0.1300 177    1   3.449  0.0039
# ALS / CTR    1.356 0.1250 177    1   3.294  0.0065

results_V1 <- run_gamma_analysis(plot_data, "V1")
print(results_V1$anova)
print(results_V1$comparisons)

# contrast     ratio    SE  df null t.ratio p.value
# ALS / PGMC   1.327 0.141 101    1   2.666  0.0436
# ALS / CTR    1.552 0.150 101    1   4.544 <0.0001

# function for violin plots
create_violin_gamma_plot <- function(gamma_results,visit_name) {
  
  adjusted_data = gamma_results$data_adj
  
  group_colors <- c(
    "ALS" = "#F8766D",    # pink
    "PGMC" = "#7CAE00",   # greenish-blue
    "CTR" = "#00BFC4",    # turquoise
    "mimic" = "#C77CFF",   # purple
    "C9orf72" = "#cfb624",
    "SOD1" = "#45ba78",
    "TARDBP" = "#9b5ccb",
    "Other" = "#936c77"
  )
  
  contrast_df <- as.data.frame(gamma_results$comparisons) %>%
    tidyr::separate(contrast, into = c("group1", "group2"), sep = " / ") %>%
    arrange(p.value) %>%
    mutate(
      y.position = max(adjusted_data$score_adj, na.rm = TRUE) * 1.2 +
        seq(0, by = 0.5, length.out = n()),
      p.signif = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        p.value < 0.05 ~ "*",
        p.value < 0.1  ~ paste0("p = ", signif(p.value, 2)),
        TRUE ~ "ns"
      )
    )
  
  p <- ggplot(adjusted_data, aes(x = type, y = score_adj, fill = type)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.8) +
    stat_summary(aes(shape = "Median"),
      fun = median,
      geom = "point",
      size = 3,
      color = "black",
      fill = "grey70") +
    geom_point(
      data = gamma_results$plot_data,
      aes(x = type, y = response, shape = "Mean"),
      size = 2,
      fill = "white",
      color = "black",
      stroke = 1.2,
      inherit.aes = FALSE
    ) +
    ggpubr::stat_pvalue_manual(
      contrast_df,
      label = "p.signif",
      hide.ns = TRUE,
      tip.length = 0.01
    ) +
    annotate(
      "text",
      x = 2.5,
      y = max(contrast_df %>% 
                filter(p.signif != "ns") %>%
                pull(y.position)) + 1,
      label = paste0(
        "ANOVA p = ",
        formatC(gamma_results$anova_pvalue, format = "e", digits = 2)
      ),
      size = 4
    ) +
    scale_shape_manual(
      name = "",
      values = c("Median" = 23, "Mean" = 21)
      #values = c("Mean" = 23)
    ) +
    guides(
      fill = "none",  
      shape = guide_legend()
    ) +
    labs(
      title = paste("p75ECD at Visit", visit_name),
      x = "Patient group",
      y = "ng p75ECD / mg creatinine (adjusted)"
    ) +
    scale_fill_manual(values = group_colors) + 
    theme_minimal() +
    theme(
      legend.position = "top",
      legend.box = "horizontal",
      legend.text = element_text(size = 11),
      plot.title = element_text(face = "bold")
    )
  
  return(p)
}

# Generate plots
results_V0$data_adj$type = factor(results_V0$data_adj$type,
                                  levels = c("ALS","PGMC","CTR","mimic"))
plot_V0_adj <- create_violin_gamma_plot(results_V0,visit_name = "V0")
results_V1$data_adj$type = factor(results_V1$data_adj$type,
                                  levels = c("ALS","PGMC","CTR","mimic"))
plot_V1_adj <- create_violin_gamma_plot(results_V1,visit_name = "V1")

combined_plot_adj <- plot_grid(plot_V0_adj, plot_V1_adj, labels = c("V0", "V1"), ncol = 2)
combined_plot_adj

# Save
ggsave(
  filename = "plots/boxplots_urine/p75ECD_V0V1_fulladjusted.svg",
  plot = combined_plot_adj,
  width = 10,
  height = 5
)

#### PGMCs
plot_data_PGMC <- Urine_data %>%
  filter(type %in% c("PGMC", "CTR", "ALS")) %>%
  mutate(
    type = trimws(type),
    
    # Adjust mutation column for PGMC
    mutation_adjusted = ifelse(mutation %in% c("C9orf72", "SOD1", "TARDBP"), mutation, "Other"),
    
    # Create unified plotting group
    plot_group = case_when(
      type == "PGMC" ~ mutation_adjusted,  # split PGMC by mutation
      TRUE ~ type                          # ALS or CTR remain
    ) %>% factor(levels = c("ALS", "C9orf72", "SOD1", "TARDBP", "Other", "CTR")),
    
    `ng p75ECD/mg creatinine` = as.numeric(`ng p75ECD/mg creatinine`),
    sex = factor(sex),
    age = as.numeric(age),
    Country = factor(Country)
  ) %>%
  filter(!is.na(plot_group))   # remove NAs

plot_data_PGMC$type = plot_data_PGMC$plot_group

# get results for each visit
results_PGMC_V0 <- run_gamma_analysis(plot_data_PGMC, "V0")
print(results_PGMC_V0$anova)
print(results_PGMC_V0$comparisons)

results_PGMC_V1 <- run_gamma_analysis(plot_data_PGMC, "V1")
print(results_PGMC_V1$anova)
print(results_PGMC_V1$comparisons)

# Generate plots
plot_V0_adj_PGMC <- create_violin_gamma_plot(results_PGMC_V0,visit_name = "V0")
plot_V1_adj_PGMC <- create_violin_gamma_plot(results_PGMC_V1,visit_name = "V1")

# On this graph, for ENCALS poster, I want only significant comparison to CTR, not all comparions -> I can d later manually not to loose any info here.
combined_plot_adj_PGMC <- plot_grid(plot_V0_adj_PGMC, plot_V1_adj_PGMC, labels = c("V0", "V1"), ncol = 2)
combined_plot_adj_PGMC

# Save
ggsave(
  filename = "plots/boxplots_urine/p75ECD_V0V1_PGMCtype_ALS_CTR_fulladjusted.svg",
  plot = combined_plot_adj_PGMC,
  width = 13,
  height = 6
)

### STEP 3: Longitudinal evolution of p75ECD in ALS, PGMC, CTR and mimic

# Combine V0 and V1 adjusted data

# ANA: NOT USE ADJUSTED DATA TO PERFORM AGAIN REGRESSION, here raw data
# Here Visit -> visit and type -> status
data_model <- plot_data %>% 
  filter(Visit %in% c("V0","V1")) %>%
  mutate(
    visit = factor(Visit),
    status = fct_relevel(type, "CTR"),  # CTR as reference
    PatientID = factor(PatientID)
  )

# ANA: changed to raw data and adjust it again, and applied gamma distribution
model_glmm <- glmmTMB(`ng p75ECD/mg creatinine` ~ status * visit + age + sex + Country + (1 | PatientID),
                      family = Gamma(link = "log"),
                      data = data_model)

summary(model_glmm)
# Test if effect significantly contribute to the model. Here age significant.
# And values differe between visit differently depending on status. Country effect? 

Anova(model_glmm, type = "III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: ng p75ECD/mg creatinine
# Chisq Df Pr(>Chisq)    
# (Intercept)  13.7039  1  0.0002140 ***
#   status       16.7019  3  0.0008139 ***
#   visit         1.1975  1  0.2738184    
# age          22.3927  1  2.222e-06 ***
#   sex           0.0327  1  0.8564601    
# Country      45.4491  4  3.207e-09 ***
#   status:visit  3.9965  3  0.2618388    

# Calculate the adjusted mean per group and visit. 
emm <- emmeans(model_glmm, ~ status | visit, type = "response")  # groups at each visit
pairs(emm, adjust = "tukey")

# compare visits (V0-V1) for each group, delta p75ECD significant in ALS and PGMC
emm_visits <-emmeans(model_glmm, ~ visit | status, type = "response")
pairs(emm_visits, adjust = "tukey", reverse = TRUE) # reverse = TRUE to calculate the delta V1-V0 and not V0-V1 
# to consider, here it is ratios, not difference.
# interpretation: ALS V1/V0 = 1.27, P=0.0004, exp(log(1.27))-exp(log(1))=(1.27-1)x100 = 27% increase from V0 to V1
# PGMC V1/V0 = 1.22, P=0.0019, 22% increase from V0 to V1

# status = CTR:
# contrast ratio     SE  df null z.ratio p.value
# V1 / V0   1.07 0.0681 Inf    1   1.094  0.2738
# 
# status = ALS:
# contrast ratio     SE  df null z.ratio p.value
# V1 / V0   1.27 0.0865 Inf    1   3.572  0.0004
# 
# status = PGMC:
# contrast ratio     SE  df null z.ratio p.value
# V1 / V0   1.22 0.0790 Inf    1   3.110  0.0019
# 
# status = mimic:
# contrast ratio     SE  df null z.ratio p.value
# V1 / V0   1.13 0.1400 Inf    1   0.974  0.3302



# Can we here look at the slop difference between groups? Ex: ALS vs CTRL, PGMC CTRL, mimic vs CTRL..
# ANA: You have this in the results statusALS:visitV1, statusPGMC:visitV1, ..., these are the diff from V0 to V1 
# between ALS vs CTR, PGMC vs CTR, etc. Here we have almost significace between ALS vs CTR only

# visible in summary

summary(model_glmm) ### stayed here

#Is the statusCTRL:visitV1 the same as visitV1 here, (CRL is intercept no)? 

# statusALS:visitV1    0.172977   0.092896   1.862 0.062596 .  
# statusPGMC:visitV1   0.131502   0.090366   1.455 0.145610    
# statusmimic:visitV1  0.051430   0.139324   0.369 0.712025 

# age                   0.010001   0.002113   4.732 2.22e-06 ***


# data frame for Error bars (95% CI) per visit
emm_df <- data.frame(emm_visits)

# add slope info
slope_info = emm_df %>%
  dplyr::select(visit, status, response) %>%
  tidyr::pivot_wider(names_from = visit, values_from = response) %>%
  dplyr::mutate(
    slope = V1 / V0)

p_df <- as.data.frame(pairs(emm_visits, reverse = TRUE)) %>%
  dplyr::mutate(status = c("CTR","ALS","PGMC","mimic")) %>%
  select(status, p.value)

slope_info <- slope_info %>%
  dplyr::left_join(p_df, by = "status") %>%
  dplyr::mutate(
    stars = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "" ),
    label = ifelse(p.value > 0.05,
                   paste0("Fold-change (V1 vs V0) = ", round(slope, 2)),
                   paste0("Fold-change (V1 vs V0) = ", round(slope, 2)," (", stars,")")))

# Plot results
create_longitudinal_plot <- function(data, emm_df,group_name) {
  
  data_group <- data %>% filter(status == group_name)
  data_group$visit = factor(data_group$visit, levels = c("V0","V1"))
  emm_group = emm_df %>% filter(status == group_name)
  
  # Define colors to match violin plots
  group_colors <- c(
    "ALS" = "#F8766D",    # pink
    "PGMC" = "#7CAE00",   # greenish-blue
    "CTR" = "#00BFC4",    # turquoise
    "mimic" = "#C77CFF",   # purple
    "C9orf72" = "#cfb624",
    "SOD1" = "#45ba78",
    "TARDBP" = "#9b5ccb",
    "Other" = "#936c77"
  )
  
  
  ggplot() +
    
    # raw data
    geom_line(
      data = data_group,
      aes(x = visit, y = `ng p75ECD/mg creatinine`, group = PatientID),
      alpha = 0.3,
      color = "grey50") +
    
    geom_point(
      data = data_group,
      aes(x = visit, y = `ng p75ECD/mg creatinine`),
      alpha = 0.3,
      size = 2,
      color = "grey50") +
    
    # model means with 95% error bars
    geom_line(
      data = emm_group,
      aes(x = visit, y = response, group = status),
      color = group_colors[group_name],
      size = 1) +
    
    geom_point(
      data = emm_group,
      aes(x = visit, y = response),
      color = group_colors[group_name],
      size = 3) +
    
    geom_errorbar(
      data = emm_group,
      aes(x = visit, ymin = asymp.LCL, ymax = asymp.UCL),
      color = group_colors[group_name],
      width = 0.1,
      size = 1) +
    
    # slope info
    geom_text(
      data = slope_info %>% filter(status == group_name), 
      aes(
        x = 1.5,
        y = 15, # y = max(data_group$`ng p75ECD/mg creatinine`, na.rm = TRUE) * 1.05,
        label = label
      ),
      inherit.aes = FALSE,
      color = group_colors[group_name],
      fontface = "bold",
      size = 3.5
    ) +
    
    labs(
      x = "Visit",
      y = "ng p75ECD / mg creatinine"
    ) +
    
    coord_cartesian(ylim = c(0, 15)) + #coord_cartesian(ylim = c(0, max(data_group$`ng p75ECD/mg creatinine`, na.rm = TRUE) * 1.1))
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_blank()
    )
}

groups <- c("ALS", "PGMC", "CTR", "mimic")

plots_list <- lapply(groups, function(grp) create_longitudinal_plot(data_model,emm_df, grp))

combined_plot <- plot_grid(
  plotlist = plots_list,
  labels = groups,  # labels for each subplot
  ncol = 2          # adjust columns as needed
)

final_plot <- ggdraw() +
  draw_label(
    "Urinary p75ECD", 
    fontface = 'bold', 
    x = 0.5, y = 0.97, hjust = 0.5
  ) +
  draw_plot(combined_plot, y = 0, height = 0.95)

final_plot

# Save
ggsave(
  filename = "plots/boxplots_urine/p75ECD_longitudinal_adjusted.svg",
  plot = final_plot,
  width = 10,
  height = 8
)

### STEP 4: Longitudinal evolution of p75ECD in PGMC subtype, ALS and CTR

data_model_PGMC <- plot_data %>% 
  filter(Visit %in% c("V0","V1") & type %in% c("PGMC","ALS","CTR")) %>%
  mutate(
    type = trimws(type),
    # split PGMC by mutation
    mutation_adjusted = ifelse(mutation %in% c("C9orf72","SOD1","TARDBP"), mutation, "Other"),
    status = case_when(
      type == "PGMC" ~ mutation_adjusted,
      TRUE ~ type
    ) %>% factor(levels = c("ALS","C9orf72","SOD1","TARDBP","Other","CTR")),
    visit = factor(Visit),
    PatientID = factor(PatientID)
  ) %>% filter(!is.na(status))


# CTR as reference
data_model_PGMC$status <- relevel(data_model_PGMC$status, "CTR")

# change to gamma 
model_glmm_PGMC <- glmmTMB(`ng p75ECD/mg creatinine` ~ status * visit + age + sex + Country + (1 | PatientID),
                           family = Gamma(link = "log"),
                           data = data_model_PGMC)

summary(model_glmm_PGMC)
Anova(model_glmm_PGMC, type = "III")

# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: ng p75ECD/mg creatinine
# Chisq Df Pr(>Chisq)    
# (Intercept)  17.9792  1  2.233e-05 ***
# status       20.2899  5   0.001102 ** 
# visit         1.1110  1   0.291869    
# age          15.3493  1  8.935e-05 ***
# sex           0.1161  1   0.733334    
# Country      46.5123  4  1.927e-09 ***
# status:visit  7.4067  5   0.192107  

# Calculate the adjusted mean per group and visit.
emm_PGMC <- emmeans(model_glmm_PGMC, ~ status | visit,type = "response")
pairs(emm_PGMC, adjust = "tukey")

# Pairwise comparisons between visits (V1 - V0), delta p75ECD significant in ALS, SOD1 and C9
emm_visits_PGMC <- emmeans(model_glmm_PGMC, ~ visit | status, type = "response")
pairs(emm_visits_PGMC, adjust = "tukey", reverse = TRUE)

# status = CTR:
#   contrast ratio     SE  df null z.ratio p.value
# V1 / V0  1.069 0.0680 Inf    1   1.054  0.2919
# 
# status = ALS:
#   contrast ratio     SE  df null z.ratio p.value
# V1 / V0  1.275 0.0866 Inf    1   3.576  0.0003
# 
# status = C9orf72:
#   contrast ratio     SE  df null z.ratio p.value
# V1 / V0  1.305 0.1630 Inf    1   2.133  0.0329
# 
# status = SOD1:
#   contrast ratio     SE  df null z.ratio p.value
# V1 / V0  1.293 0.1430 Inf    1   2.334  0.0196
# 
# status = TARDBP:
#   contrast ratio     SE  df null z.ratio p.value
# V1 / V0  1.292 0.2080 Inf    1   1.593  0.1112
# 
# status = Other:
#   contrast ratio     SE  df null z.ratio p.value
# V1 / V0  0.972 0.1310 Inf    1  -0.208  0.8350

# data frame for Error bars (95% CI) per visit
emm_df_PGMC <- data.frame(emm_visits_PGMC)

# add slope info
slope_info = emm_df_PGMC %>%
  dplyr::select(visit, status, response) %>%
  tidyr::pivot_wider(names_from = visit, values_from = response) %>%
  dplyr::mutate(
    slope = V1 / V0)

p_df <- as.data.frame(pairs(emm_visits_PGMC, reverse = TRUE)) %>%
  dplyr::mutate(status = c("CTR","ALS","C9orf72","SOD1","TARDBP","Other")) %>%
  select(status, p.value)

slope_info <- slope_info %>%
  dplyr::left_join(p_df, by = "status") %>%
  dplyr::mutate(
    stars = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "" ),
    label = ifelse(p.value > 0.05,
                   paste0("Fold-change (V1 vs V0) = ", round(slope, 2)),
                   paste0("Fold-change (V1 vs V0) = ", round(slope, 2)," (", stars,")")))


# Plot results
# Plot longitudinal trajectories per group
groups <- c("ALS","C9orf72","SOD1","TARDBP","Other","CTR")

plots_list <- lapply(groups, function(grp) create_longitudinal_plot(data_model_PGMC, emm_df_PGMC,grp))

combined_plot <- plot_grid(
  plotlist = plots_list,
  labels = groups,
  ncol = 2
)

final_plot <- ggdraw() +
  draw_label(
    "Urinary p75ECD PGMC subtypes", 
    fontface = 'bold', 
    x = 0.5, y = 0.97, hjust = 0.5
  ) +
  draw_plot(combined_plot, y = 0, height = 0.95)

final_plot

# Save
ggsave(
  filename = "plots/boxplots_urine/p75ECD_longitudinal_PGMCsybtype_adjusted.svg",
  plot = final_plot,
  width = 10,
  height = 12
)

####  STEP 5: check analyses of urine data with clinical variables ALS
df_base <- Urine_data %>%
  mutate(
    `ng p75ECD/mg creatinine` = as.numeric(`ng p75ECD/mg creatinine`),
    type   = factor(type, levels = c("ALS", "PGMC", "CTR", "mimic")),
    sex    = factor(sex),
    age    = as.numeric(age),
    Country = factor(Country),
    Visit  = factor(Visit, levels = c("V0", "V1")),
    PatientID = factor(PatientID)
  ) %>%
  filter(type == "ALS")

## -> ALSFRS-R
df_A <- df_base %>%
  mutate(ALSFRSR_Total = as.numeric(`Total ALSFRS-r`)) %>%
  filter(!is.na(ALSFRSR_Total), !is.na(`ng p75ECD/mg creatinine`))

model_A <- glmmTMB(
  `ng p75ECD/mg creatinine` ~
    ALSFRSR_Total + Visit + age + sex + Country + (1 | PatientID),
  family = Gamma(link = "log"),
  data = df_A
)

# p-value and effect size
anova_A  <- car::Anova(model_A, type = "III")
p_A      <- anova_A$`Pr(>Chisq)`[rownames(anova_A) == "ALSFRSR_Total"]
beta_A   <- fixef(model_A)$cond["ALSFRSR_Total"]
pct_A    <- (exp(-beta_A) - 1) * 100  # per 1-point decrease

label_A <- paste0(
  round(pct_A, 1), "% increase in p75ECD per 1-pt ALSFRS-R decrease\n",
  "p = ", signif(p_A, 2)
)

# Marginal prediction (Visit, age, sex, Country held at reference)
pred_A <- data.frame(
  ALSFRSR_Total = seq(min(df_A$ALSFRSR_Total), max(df_A$ALSFRSR_Total), length.out = 300),
  Visit   = levels(df_A$Visit)[1],
  age     = mean(df_A$age, na.rm = TRUE),
  sex     = levels(df_A$sex)[1],
  Country = levels(df_A$Country)[1]
)
pred_A$predicted <- predict(model_A, newdata = pred_A, type = "response", re.form = NA)

pred_link <- predict(model_A, newdata = pred_A, type = "link", re.form = NA, se.fit = TRUE)

pred_A$predicted <- exp(pred_link$fit)
pred_A$lower     <- exp(pred_link$fit - 1.96 * pred_link$se.fit)
pred_A$upper     <- exp(pred_link$fit + 1.96 * pred_link$se.fit)

plot_A <- ggplot(df_A, aes(x = ALSFRSR_Total, y = `ng p75ECD/mg creatinine`)) +
  geom_point(aes(shape = Visit), alpha = 0.6, size = 2) +
  geom_ribbon(
    data = pred_A,
    aes(x = ALSFRSR_Total, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.4) +
  geom_line(data = pred_A, aes(x = ALSFRSR_Total, y = predicted),
            color = "black", size = 1.2, inherit.aes = FALSE,
            linetype = "dashed") +
  annotate("text", x = 2, y = max(df_A$`ng p75ECD/mg creatinine`, na.rm = TRUE) + 1.5,
           label = label_A, hjust = 0, vjust = 1, size = 3.5, fontface = "italic") +
  scale_x_continuous(limits = c(0, 48), expand = c(0, 0)) +
  labs(x = "ALSFRS-R Total Score",
       y = "Urinary p75ECD / mg creatinine",
       title = "Association between urinary p75ECD and ALSFRS-R (all visits)") +
  theme_minimal()

plot_A

ggsave("plots/boxplots_urine/panel_ALSFRS-R.svg", plot_A, width = 5, height = 8)

# --------------------------
## -> Age in ALS only
df_B <- df_base %>%
  filter(!is.na(age), !is.na(`ng p75ECD/mg creatinine`))

model_B <- glmmTMB(
  `ng p75ECD/mg creatinine` ~
    age + Visit + sex + Country + (1 | PatientID),
  family = Gamma(link = "log"),
  data = df_B
)

# p-value and effect size
anova_B <- car::Anova(model_B, type = "III")
p_B     <- anova_B$`Pr(>Chisq)`[rownames(anova_B) == "age"]
beta_B  <- fixef(model_B)$cond["age"]
pct_B   <- (exp(beta_B) - 1) * 100  # per 1-year increase

label_B <- paste0(
  round(abs(pct_B), 1), "% ", ifelse(pct_B > 0, "increase", "decrease"),
  " in p75ECD per 1-yr age increase\n",
  "p = ", signif(p_B, 2)
)

# Marginal prediction (Visit, sex, Country held at reference)
pred_B <- data.frame(
  age     = seq(min(df_B$age, na.rm = TRUE), max(df_B$age, na.rm = TRUE), length.out = 300),
  Visit   = levels(df_B$Visit)[1],
  sex     = levels(df_B$sex)[1],
  Country = levels(df_B$Country)[1])

pred_link_B      <- predict(model_B, newdata = pred_B, type = "link", re.form = NA, se.fit = TRUE)
pred_B$predicted <- exp(pred_link_B$fit)
pred_B$lower     <- exp(pred_link_B$fit - 1.96 * pred_link_B$se.fit)
pred_B$upper     <- exp(pred_link_B$fit + 1.96 * pred_link_B$se.fit)

plot_B <- ggplot(df_B, aes(x = age, y = `ng p75ECD/mg creatinine`)) +
  geom_point(aes(shape = Visit), alpha = 0.6, size = 2) +
  geom_ribbon(
    data = pred_B,
    aes(x = age, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.4
  ) +
  geom_line(
    data = pred_B,
    aes(x = age, y = predicted),
    color = "black",
    size = 1.2,
    inherit.aes = FALSE,
    linetype = "dashed"
  ) +
  annotate("text",
           x = min(df_B$age, na.rm = TRUE),
           y = max(df_B$`ng p75ECD/mg creatinine`, na.rm = TRUE) + 1.5,
           label = label_B,
           hjust = 0, vjust = 1, size = 3.5, fontface = "italic") +
  scale_x_continuous(
    limits = c(floor(min(df_B$age, na.rm = TRUE)), ceiling(max(df_B$age, na.rm = TRUE))),
    expand = c(0, 0)
  ) +
  labs(
    x = "Age (years)",
    y = "Urinary p75ECD / mg creatinine",
    title = "Association between urinary p75ECD and Age (all visits)"
  ) +
  theme_minimal()

plot_B

ggsave("plots/boxplots_urine/panel_age.svg", plot_B, width = 5, height = 8)


##### Age all cohort

df_base_all <- Urine_data %>%
  mutate(
    `ng p75ECD/mg creatinine` = as.numeric(`ng p75ECD/mg creatinine`),
    type   = factor(type, levels = c("ALS", "PGMC", "CTR", "mimic")),
    sex    = factor(sex),
    age    = as.numeric(age),
    Country = factor(Country),
    Visit  = factor(Visit, levels = c("V0", "V1")),
    PatientID = factor(PatientID)
  )

df_B_all <- df_base_all %>%
  filter(!is.na(age), !is.na(`ng p75ECD/mg creatinine`))

model_B <- glmmTMB(
  `ng p75ECD/mg creatinine` ~
    age + Visit + sex + Country + (1 | PatientID),
  family = Gamma(link = "log"),
  data = df_B_all
)

# p-value and effect size
anova_B <- car::Anova(model_B, type = "III")
p_B     <- anova_B$`Pr(>Chisq)`[rownames(anova_B) == "age"]
beta_B  <- fixef(model_B)$cond["age"]
pct_B   <- (exp(beta_B) - 1) * 100  # per 1-year increase

label_B <- paste0(
  round(abs(pct_B), 1), "% ", ifelse(pct_B > 0, "increase", "decrease"),
  " in p75ECD per 1-yr age increase\n",
  "p = ", signif(p_B, 2)
)

# Marginal prediction (Visit, sex, Country held at reference)
pred_B <- data.frame(
  age     = seq(min(df_B_all$age, na.rm = TRUE), max(df_B_all$age, na.rm = TRUE), length.out = 300),
  Visit   = levels(df_B_all$Visit)[1],
  sex     = levels(df_B_all$sex)[1],
  Country = levels(df_B_all$Country)[1])

pred_link_B      <- predict(model_B, newdata = pred_B, type = "link", re.form = NA, se.fit = TRUE)
pred_B$predicted <- exp(pred_link_B$fit)
pred_B$lower     <- exp(pred_link_B$fit - 1.96 * pred_link_B$se.fit)
pred_B$upper     <- exp(pred_link_B$fit + 1.96 * pred_link_B$se.fit)

plot_B <- ggplot(df_B_all, aes(x = age, y = `ng p75ECD/mg creatinine`)) +
  geom_point(aes(shape = Visit), alpha = 0.6, size = 2) +
  geom_ribbon(
    data = pred_B,
    aes(x = age, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.4
  ) +
  geom_line(
    data = pred_B,
    aes(x = age, y = predicted),
    color = "black",
    size = 1.2,
    inherit.aes = FALSE,
    linetype = "dashed"
  ) +
  annotate("text",
           x = min(df_B_all$age, na.rm = TRUE),
           y = max(df_B_all$`ng p75ECD/mg creatinine`, na.rm = TRUE) + 1.5,
           label = label_B,
           hjust = 0, vjust = 1, size = 3.5, fontface = "italic") +
  scale_x_continuous(
    limits = c(floor(min(df_B_all$age, na.rm = TRUE)), ceiling(max(df_B_all$age, na.rm = TRUE))),
    expand = c(0, 0)
  ) +
  labs(
    x = "Age (years)",
    y = "Urinary p75ECD / mg creatinine",
    title = "Association between urinary p75ECD and Age"
  ) +
  theme_minimal()

plot_B

ggsave("plots/boxplots_urine/panel_age_all.svg", plot_B, width = 5, height = 8)

# Age in control only

df_base <- Urine_data %>%
  mutate(
    `ng p75ECD/mg creatinine` = as.numeric(`ng p75ECD/mg creatinine`),
    type   = factor(type, levels = c("ALS", "PGMC", "CTR", "mimic")),
    sex    = factor(sex),
    age    = as.numeric(age),
    Country = factor(Country),
    Visit  = factor(Visit, levels = c("V0", "V1")),
    PatientID = factor(PatientID)
  ) %>%
  filter(type == "CTR")

df_B <- df_base %>%
  filter(!is.na(age), !is.na(`ng p75ECD/mg creatinine`))

model_B <- glmmTMB(
  `ng p75ECD/mg creatinine` ~
    age + Visit + sex + Country + (1 | PatientID),
  family = Gamma(link = "log"),
  data = df_B
)

# p-value and effect size
anova_B <- car::Anova(model_B, type = "III")
p_B     <- anova_B$`Pr(>Chisq)`[rownames(anova_B) == "age"]
beta_B  <- fixef(model_B)$cond["age"]
pct_B   <- (exp(beta_B) - 1) * 100  # per 1-year increase

label_B <- paste0(
  round(abs(pct_B), 1), "% ", ifelse(pct_B > 0, "increase", "decrease"),
  " in p75ECD per 1-yr age increase\n",
  "p = ", signif(p_B, 2)
)

# Marginal prediction (Visit, sex, Country held at reference)
pred_B <- data.frame(
  age     = seq(min(df_B$age, na.rm = TRUE), max(df_B$age, na.rm = TRUE), length.out = 300),
  Visit   = levels(df_B$Visit)[1],
  sex     = levels(df_B$sex)[1],
  Country = levels(df_B$Country)[1])

pred_link_B      <- predict(model_B, newdata = pred_B, type = "link", re.form = NA, se.fit = TRUE)
pred_B$predicted <- exp(pred_link_B$fit)
pred_B$lower     <- exp(pred_link_B$fit - 1.96 * pred_link_B$se.fit)
pred_B$upper     <- exp(pred_link_B$fit + 1.96 * pred_link_B$se.fit)

plot_B <- ggplot(df_B, aes(x = age, y = `ng p75ECD/mg creatinine`)) +
  geom_point(aes(shape = Visit), alpha = 0.6, size = 2) +
  geom_ribbon(
    data = pred_B,
    aes(x = age, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.4
  ) +
  geom_line(
    data = pred_B,
    aes(x = age, y = predicted),
    color = "black",
    size = 1.2,
    inherit.aes = FALSE,
    linetype = "dashed"
  ) +
  annotate("text",
           x = min(df_B$age, na.rm = TRUE),
           y = max(df_B$`ng p75ECD/mg creatinine`, na.rm = TRUE) + 1.5,
           label = label_B,
           hjust = 0, vjust = 1, size = 3.5, fontface = "italic") +
  scale_x_continuous(
    limits = c(floor(min(df_B$age, na.rm = TRUE)), ceiling(max(df_B$age, na.rm = TRUE))),
    expand = c(0, 0)
  ) +
  labs(
    x = "Age (years)",
    y = "Urinary p75ECD / mg creatinine",
    title = "Association between urinary p75ECD and Age (all visits)"
  ) +
  theme_minimal()

plot_B

ggsave("plots/boxplots_urine/panel_age_CTR.svg", plot_B, width = 5, height = 8)



# --------------------------
## -> Correlation with motor neuron score (ALS only) (Laura test -> graph looks not good). (Or does it look better with LMN or UMM signs only?)

df_motor <- Urine_data %>%
  mutate(
    `ng p75ECD/mg creatinine` = as.numeric(`ng p75ECD/mg creatinine`),
    type              = factor(type, levels = c("ALS", "PGMC", "CTR", "mimic")),
    sex               = factor(sex),
    age               = as.numeric(age),
    Country           = factor(Country),
    Visit             = factor(Visit, levels = c("V0", "V1")),
    PatientID         = factor(PatientID),
    TotalMotorSigns   = as.numeric(TotalMotorSigns)
  )

df_E_all <- df_motor %>%
  filter(
    type == "ALS",
    !is.na(age),
    !is.na(`ng p75ECD/mg creatinine`),
    !is.na(TotalMotorSigns)
  )

# Model:
# Association between motor signs and urinary p75ECD,
# adjusted for age, visit, sex, country
model_E <- glmmTMB(
  `ng p75ECD/mg creatinine` ~
    TotalMotorSigns + age + Visit + sex + Country + (1 | PatientID),
  family = Gamma(link = "log"),
  data = df_E_all
)

# p-value and effect size
anova_E <- car::Anova(model_E, type = "III")

p_E <- anova_E$`Pr(>Chisq)`[
  rownames(anova_E) == "TotalMotorSigns"
]

beta_E <- fixef(model_E)$cond["TotalMotorSigns"]

# % change in p75ECD per 1-point increase in motor score
pct_E <- (exp(beta_E) - 1) * 100

label_E <- paste0(
  round(abs(pct_E), 1), "% ",
  ifelse(pct_E > 0, "increase", "decrease"),
  " in p75ECD per 1-point motor score increase\n",
  "p = ", signif(p_E, 2)
)

# Marginal prediction across motor score
# Age, Visit, sex and Country are held constant
pred_E <- data.frame(
  TotalMotorSigns = seq(
    min(df_E_all$TotalMotorSigns, na.rm = TRUE),
    max(df_E_all$TotalMotorSigns, na.rm = TRUE),
    length.out = 300
  ),
  age     = median(df_E_all$age, na.rm = TRUE), # ANA to Laura: changed here to median
  Visit   = levels(df_E_all$Visit)[1],
  sex     = levels(df_E_all$sex)[1],
  Country = levels(df_E_all$Country)[1]
)

pred_link_E <- predict(
  model_E,
  newdata = pred_E,
  type = "link",
  re.form = NA,
  se.fit = TRUE
)

pred_E$predicted <- exp(pred_link_E$fit)
pred_E$lower     <- exp(pred_link_E$fit - 1.96 * pred_link_E$se.fit)
pred_E$upper     <- exp(pred_link_E$fit + 1.96 * pred_link_E$se.fit)

# Plot
plot_E <- ggplot(
  df_E_all,
  aes(
    x = TotalMotorSigns,
    y = `ng p75ECD/mg creatinine`
  )
) +
  geom_point(
    aes(shape = Visit),
    alpha = 0.6,
    size = 2,
    position = position_jitter(width = 0.2, height = 0) # new line for jitter of points
  ) +
  geom_ribbon(
    data = pred_E,
    aes(
      x = TotalMotorSigns,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.4
  ) +
  geom_line(
    data = pred_E,
    aes(
      x = TotalMotorSigns,
      y = predicted
    ),
    color = "black",
    size = 1.2,
    inherit.aes = FALSE,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = min(df_E_all$TotalMotorSigns, na.rm = TRUE),
    y = max(
      df_E_all$`ng p75ECD/mg creatinine`,
      na.rm = TRUE
    ) + 1.5,
    label = label_E,
    hjust = 0,
    vjust = 1,
    size = 3.5,
    fontface = "italic"
  ) +
  scale_x_continuous(
    breaks = 0:8
    #limits = c(0, 8),
    #expand = c(0, 0)
  ) +
  labs(
    x = "Total motor signs score",
    y = "Urinary p75ECD / mg creatinine",
    title = "Association between urinary p75ECD and motor signs"
  ) +
  theme_minimal()

plot_E

ggsave(
  "plots/boxplots_urine/panel_motor_signs_all.svg",
  plot_E,
  width = 5,
  height = 8
)

## Ana: separate the upper and lower motor neuron analyses (08/09)
## -> Correlation with UPPER motor neuron signs (ALS only)

df_motor <- Urine_data %>%
  mutate(
    `ng p75ECD/mg creatinine` = as.numeric(`ng p75ECD/mg creatinine`),
    type              = factor(type, levels = c("ALS", "PGMC", "CTR", "mimic")),
    sex               = factor(sex),
    age               = as.numeric(age),
    Country           = factor(Country),
    Visit             = factor(Visit, levels = c("V0", "V1")),
    PatientID         = factor(PatientID),
    UpperMotonScore   = as.numeric(UpperMotonScore),
    LowerMotonScore   = as.numeric(LowerMotonScore)
  )

df_E_UMN <- df_motor %>%
  filter(
    type == "ALS",
    !is.na(age),
    !is.na(`ng p75ECD/mg creatinine`),
    !is.na(UpperMotonScore)
  )

# Model
model_E_UMN <- glmmTMB(
  `ng p75ECD/mg creatinine` ~
    UpperMotonScore + age + Visit + sex + Country + (1 | PatientID),
  family = Gamma(link = "log"),
  data = df_E_UMN
)

# p-value and effect size
anova_E_UMN <- car::Anova(model_E_UMN, type = "III")

p_E_UMN <- anova_E_UMN$`Pr(>Chisq)`[
  rownames(anova_E_UMN) == "UpperMotonScore"
]

beta_E_UMN <- fixef(model_E_UMN)$cond["UpperMotonScore"]

# % change in p75ECD per 1-point increase in UMN score
pct_E_UMN <- (exp(beta_E_UMN) - 1) * 100

label_E_UMN <- paste0(
  round(abs(pct_E_UMN), 1), "% ",
  ifelse(pct_E_UMN > 0, "increase", "decrease"),
  " in p75ECD per 1-point UMN score increase\n",
  "p = ", signif(p_E_UMN, 2)
)

# Marginal prediction across UMN score
# Age, Visit, sex and Country are held constant
pred_E_UMN <- data.frame(
  UpperMotonScore = seq(
    min(df_E_UMN$UpperMotonScore, na.rm = TRUE),
    max(df_E_UMN$UpperMotonScore, na.rm = TRUE),
    length.out = 300
  ),
  age     = mean(df_E_UMN$age, na.rm = TRUE),
  Visit   = levels(df_E_UMN$Visit)[1],
  sex     = levels(df_E_UMN$sex)[1],
  Country = levels(df_E_UMN$Country)[1]
)

pred_link_E_UMN <- predict(
  model_E_UMN,
  newdata = pred_E_UMN,
  type = "link",
  re.form = NA,
  se.fit = TRUE
)

pred_E_UMN$predicted <- exp(pred_link_E_UMN$fit)
pred_E_UMN$lower     <- exp(pred_link_E_UMN$fit - 1.96 * pred_link_E_UMN$se.fit)
pred_E_UMN$upper     <- exp(pred_link_E_UMN$fit + 1.96 * pred_link_E_UMN$se.fit)

# x-axis padding to accommodate jitter without clipping points
umn_range <- range(df_E_UMN$UpperMotonScore, na.rm = TRUE)
umn_pad   <- 0.3

# Plot
plot_E_UMN <- ggplot(
  df_E_UMN,
  aes(
    x = UpperMotonScore,
    y = `ng p75ECD/mg creatinine`
  )
) +
  geom_point(
    aes(shape = Visit),
    alpha = 0.6,
    size = 2,
    position = position_jitter(width = 0.2, height = 0)
  ) +
  geom_ribbon(
    data = pred_E_UMN,
    aes(
      x = UpperMotonScore,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.4
  ) +
  geom_line(
    data = pred_E_UMN,
    aes(
      x = UpperMotonScore,
      y = predicted
    ),
    color = "black",
    size = 1.2,
    inherit.aes = FALSE,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = umn_range[1],
    y = max(
      df_E_UMN$`ng p75ECD/mg creatinine`,
      na.rm = TRUE
    ) + 1.5,
    label = label_E_UMN,
    hjust = 0,
    vjust = 1,
    size = 3.5,
    fontface = "italic"
  ) +
  scale_x_continuous(
    breaks = seq(umn_range[1], umn_range[2], by = 1),
    limits = c(umn_range[1] - umn_pad, umn_range[2] + umn_pad),
    expand = c(0, 0)
  ) +
  labs(
    x = "Upper motor neuron signs score",
    y = "Urinary p75ECD / mg creatinine",
    title = "Association between urinary p75ECD and UMN signs"
  ) +
  theme_minimal()

plot_E_UMN

ggsave(
  "plots/boxplots_urine/panel_motor_signs_UMN.svg",
  plot_E_UMN,
  width = 5, height = 8
)

## -> Correlation with LOWER motor neuron signs (ALS only)

df_E_LMN <- df_motor %>%
  filter(
    type == "ALS",
    !is.na(age),
    !is.na(`ng p75ECD/mg creatinine`),
    !is.na(LowerMotonScore)
  )

# Model
model_E_LMN <- glmmTMB(
  `ng p75ECD/mg creatinine` ~
    LowerMotonScore + age + Visit + sex + Country + (1 | PatientID),
  family = Gamma(link = "log"),
  data = df_E_LMN
)

# p-value and effect size
anova_E_LMN <- car::Anova(model_E_LMN, type = "III")

p_E_LMN <- anova_E_LMN$`Pr(>Chisq)`[
  rownames(anova_E_LMN) == "LowerMotonScore"
]

beta_E_LMN <- fixef(model_E_LMN)$cond["LowerMotonScore"]

# % change in p75ECD per 1-point increase in LMN score
pct_E_LMN <- (exp(beta_E_LMN) - 1) * 100

label_E_LMN <- paste0(
  round(abs(pct_E_LMN), 1), "% ",
  ifelse(pct_E_LMN > 0, "increase", "decrease"),
  " in p75ECD per 1-point LMN score increase\n",
  "p = ", signif(p_E_LMN, 2)
)

# Marginal prediction across LMN score
# Age, Visit, sex and Country are held constant
pred_E_LMN <- data.frame(
  LowerMotonScore = seq(
    min(df_E_LMN$LowerMotonScore, na.rm = TRUE),
    max(df_E_LMN$LowerMotonScore, na.rm = TRUE),
    length.out = 300
  ),
  age     = mean(df_E_LMN$age, na.rm = TRUE),
  Visit   = levels(df_E_LMN$Visit)[1],
  sex     = levels(df_E_LMN$sex)[1],
  Country = levels(df_E_LMN$Country)[1]
)

pred_link_E_LMN <- predict(
  model_E_LMN,
  newdata = pred_E_LMN,
  type = "link",
  re.form = NA,
  se.fit = TRUE
)

pred_E_LMN$predicted <- exp(pred_link_E_LMN$fit)
pred_E_LMN$lower     <- exp(pred_link_E_LMN$fit - 1.96 * pred_link_E_LMN$se.fit)
pred_E_LMN$upper     <- exp(pred_link_E_LMN$fit + 1.96 * pred_link_E_LMN$se.fit)

# x-axis padding to accommodate jitter without clipping points
lmn_range <- range(df_E_LMN$LowerMotonScore, na.rm = TRUE)
lmn_pad   <- 0.3

# Plot
plot_E_LMN <- ggplot(
  df_E_LMN,
  aes(
    x = LowerMotonScore,
    y = `ng p75ECD/mg creatinine`
  )
) +
  geom_point(
    aes(shape = Visit),
    alpha = 0.6,
    size = 2,
    position = position_jitter(width = 0.2, height = 0)
  ) +
  geom_ribbon(
    data = pred_E_LMN,
    aes(
      x = LowerMotonScore,
      ymin = lower,
      ymax = upper
    ),
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.4
  ) +
  geom_line(
    data = pred_E_LMN,
    aes(
      x = LowerMotonScore,
      y = predicted
    ),
    color = "black",
    size = 1.2,
    inherit.aes = FALSE,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = lmn_range[1],
    y = max(
      df_E_LMN$`ng p75ECD/mg creatinine`,
      na.rm = TRUE
    ) + 1.5,
    label = label_E_LMN,
    hjust = 0,
    vjust = 1,
    size = 3.5,
    fontface = "italic"
  ) +
  scale_x_continuous(
    breaks = seq(lmn_range[1], lmn_range[2], by = 1),
    limits = c(lmn_range[1] - lmn_pad, lmn_range[2] + lmn_pad),
    expand = c(0, 0)
  ) +
  labs(
    x = "Lower motor neuron signs score",
    y = "Urinary p75ECD / mg creatinine",
    title = "Association between urinary p75ECD and LMN signs"
  ) +
  theme_minimal()

plot_E_LMN

ggsave(
  "plots/boxplots_urine/panel_motor_signs_LMN.svg",
  plot_E_LMN,
  width = 5, height = 8
)



# --------------------------
## -> Sex
df_C <- Urine_data %>%
  mutate(
    `ng p75ECD/mg creatinine` = as.numeric(`ng p75ECD/mg creatinine`),
    type             = factor(type, levels = c("ALS", "PGMC", "CTR", "mimic")),
    sex              = factor(sex, levels = c("M", "F")),
    age              = as.numeric(age),
    Country          = factor(Country),
    Visit            = factor(Visit, levels = c("V0", "V1")),
    PatientID        = factor(PatientID),
    SiteDiseaseOnset = factor(SiteDiseaseOnset, levels = c("spinal", "bulbar"))
  ) %>%
  filter(
    type == "ALS",
    sex %in% c("M", "F"),
    !is.na(`ng p75ECD/mg creatinine`)
  )

# Mixed model
model_C <- glmmTMB(
  `ng p75ECD/mg creatinine` ~
    sex + Visit + age  + Country + (1 | PatientID),
  family = Gamma(link = "log"),
  data = df_C
)

# p-value
anova_C <- car::Anova(model_C, type = "III")
p_C     <- anova_C$`Pr(>Chisq)`[rownames(anova_C) == "sex"]
label_C <- paste0("p = ", signif(p_C, 2))

# summary stats for boxplot annotation positions
y_max <- max(df_C$`ng p75ECD/mg creatinine`, na.rm = TRUE)


plot_C <- ggplot(df_C, aes(x = sex,
                           y = `ng p75ECD/mg creatinine`,
                           fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5,
               show.legend = FALSE) +
  geom_jitter(aes(shape = Visit), width = 0.15, alpha = 0.5, size = 1.8,
              show.legend = TRUE) +
  
  # p-value bracket
  annotate("segment",
           x = 1, xend = 2,
           y = y_max + 1, yend = y_max + 1) +
  annotate("segment", x = 1, xend = 1,
           y = y_max + 0.5, yend = y_max + 1) +
  annotate("segment", x = 2, xend = 2,
           y = y_max + 0.5, yend = y_max + 1) +
  annotate("text",
           x = 1.5, y = y_max + 1.5,
           label = label_C,
           hjust = 0.5, size = 4, fontface = "italic") +
  
  scale_fill_manual(values = c("M" = "lightblue", "F" = "pink")) +
  scale_x_discrete(labels = c("M" = "Male", "F" = "Female")) +
  
  labs(
    x = "Sex",
    y = "Urinary p75ECD / mg creatinine",
    title = "Urinary p75ECD by Sex (all visits)"
  ) +
  
  theme_minimal() +
  theme(legend.position = "bottom")

plot_C

ggsave("plots/boxplots_urine/panel_sex.svg", plot_C, width = 5, height = 5)

## -> Sex all cohort (Laura)
df_C_all <- Urine_data %>%
  mutate(
    `ng p75ECD/mg creatinine` = as.numeric(`ng p75ECD/mg creatinine`),
    type             = factor(type, levels = c("ALS", "PGMC", "CTR", "mimic")),
    sex              = factor(sex, levels = c("M", "F")),
    age              = as.numeric(age),
    Country          = factor(Country),
    Visit            = factor(Visit, levels = c("V0", "V1")),
    PatientID        = factor(PatientID),
    SiteDiseaseOnset = factor(SiteDiseaseOnset, levels = c("spinal", "bulbar"))
  ) %>%
  filter(
    sex %in% c("M", "F"),
    !is.na(`ng p75ECD/mg creatinine`)
  )

# Mixed model
model_C <- glmmTMB(
  `ng p75ECD/mg creatinine` ~
    sex + Visit + age  + Country + (1 | PatientID),
  family = Gamma(link = "log"),
  data = df_C_all
)

# p-value
anova_C <- car::Anova(model_C, type = "III")
p_C     <- anova_C$`Pr(>Chisq)`[rownames(anova_C) == "sex"]
label_C <- paste0("p = ", signif(p_C, 2))

# summary stats for boxplot annotation positions
y_max <- max(df_C_all$`ng p75ECD/mg creatinine`, na.rm = TRUE)


plot_C <- ggplot(df_C_all, aes(x = sex,
                           y = `ng p75ECD/mg creatinine`,
                           fill = sex)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5,
               show.legend = FALSE) +
  geom_jitter(aes(shape = Visit), width = 0.15, alpha = 0.5, size = 1.8,
              show.legend = TRUE) +
  
  # p-value bracket
  annotate("segment",
           x = 1, xend = 2,
           y = y_max + 1, yend = y_max + 1) +
  annotate("segment", x = 1, xend = 1,
           y = y_max + 0.5, yend = y_max + 1) +
  annotate("segment", x = 2, xend = 2,
           y = y_max + 0.5, yend = y_max + 1) +
  annotate("text",
           x = 1.5, y = y_max + 1.5,
           label = label_C,
           hjust = 0.5, size = 4, fontface = "italic") +
  
  scale_fill_manual(values = c("M" = "lightblue", "F" = "pink")) +
  scale_x_discrete(labels = c("M" = "Male", "F" = "Female")) +
  
  labs(
    x = "Sex",
    y = "Urinary p75ECD / mg creatinine",
    title = "Urinary p75ECD by Sex (all visits)"
  ) +
  
  theme_minimal() +
  theme(legend.position = "bottom")

plot_C

ggsave("plots/boxplots_urine/panel_sex_all.svg", plot_C, width = 5, height = 5)

# --------------------------
## -> Site of onset 
df_D <- Urine_data %>%
  mutate(
    `ng p75ECD/mg creatinine` = as.numeric(`ng p75ECD/mg creatinine`),
    type             = factor(type, levels = c("ALS", "PGMC", "CTR", "mimic")),
    sex              = factor(sex),
    age              = as.numeric(age),
    Country          = factor(Country),
    Visit            = factor(Visit, levels = c("V0", "V1")),
    PatientID        = factor(PatientID),
    SiteDiseaseOnset = factor(SiteDiseaseOnset, levels = c("spinal", "bulbar"))
  ) %>%
  filter(
    type == "ALS",
    SiteDiseaseOnset %in% c("spinal", "bulbar"),
    !is.na(`ng p75ECD/mg creatinine`)
  )

# Mixed model
model_D <- glmmTMB(
  `ng p75ECD/mg creatinine` ~
    SiteDiseaseOnset + Visit + age + sex + Country + (1 | PatientID),
  family = Gamma(link = "log"),
  data = df_D
)

# p-value
anova_D <- car::Anova(model_D, type = "III")
p_D     <- anova_D$`Pr(>Chisq)`[rownames(anova_D) == "SiteDiseaseOnset"]
label_D <- paste0("p = ", signif(p_D, 2))

p_D

# summary stats for boxplot annotation positions
y_max <- max(df_D$`ng p75ECD/mg creatinine`, na.rm = TRUE)

plot_D <- ggplot(df_D, aes(x = SiteDiseaseOnset,
                           y = `ng p75ECD/mg creatinine`,
                           fill = SiteDiseaseOnset)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5,
               show.legend = FALSE) +
  geom_jitter(aes(shape = Visit), width = 0.15, alpha = 0.5, size = 1.8,
              show.legend = TRUE) +
  
  # p-value bracket
  annotate("segment",
           x = 1, xend = 2,
           y = y_max + 1, yend = y_max + 1) +
  annotate("segment", x = 1, xend = 1,
           y = y_max + 0.5, yend = y_max + 1) +
  annotate("segment", x = 2, xend = 2,
           y = y_max + 0.5, yend = y_max + 1) +
  annotate("text",
           x = 1.5, y = y_max + 1.5,
           label = label_D,
           hjust = 0.5, size = 4, fontface = "italic") +
  
  scale_fill_manual(values = c("spinal" = "lightgrey", "bulbar" = "#3e3e3e")) + # scale_fill_manual(values = c("spinal" = "#4E9BB9", "bulbar" = "#E07B54")) 
  
  labs(
    x = "Site of Disease Onset",
    y = "Urinary p75ECD / mg creatinine",
    title = "Urinary p75ECD by Site of Onset (all visits)"
  ) +
  
  theme_minimal() +
  theme(legend.position = "bottom")

plot_D

ggsave("plots/boxplots_urine/panel_site_of_onset.svg", plot_D, width = 5, height = 5)

###################
# HEATMAP

all_participants_ID_visits <- read_excel("results/all_participants_ID_visits.xlsx")

heatmap_df = Urine_data %>%
  left_join(all_participants_ID_visits %>% 
              select(PatientID,TotalMotorSigns,Visit)) %>%
  mutate(
    ParticipantCode = factor(ParticipantCode),
    Visit = factor(Visit, levels = c("V0", "V1")),
    p75ECD = as.numeric(`ng p75ECD/mg creatinine`),
    MotorSigns = as.numeric(TotalMotorSigns) ,
    type = factor(type,
                  levels = c("ALS", "PGMC", "CTR", "mimic"))
  ) %>%
  filter(
    !is.na(p75ECD),
    !is.na(MotorSigns)
  ) %>%
  distinct()

# changee in motor signs V0 to V1
change_df = heatmap_df %>%
  select(ParticipantCode, Visit, MotorSigns,type) %>%
  pivot_wider(names_from = Visit, values_from = MotorSigns) %>%
  mutate(
    motor_change = case_when(
      V1 > V0 ~ "increase",
      V1 < V0 ~ "decrease",
      TRUE ~ "stable"
    )
  )

create_heatmap <- function(visit_name) {
  
  df_visit = heatmap_df %>% 
    filter(Visit == visit_name) %>%
    arrange(desc(p75ECD))
  
  mat_df = df_visit %>%
    group_by(ParticipantCode) %>%
    summarise(p75ECD = first(p75ECD)) %>%
    arrange(desc(p75ECD)) %>%
    pivot_wider(names_from = ParticipantCode, values_from = p75ECD)
  
  participant_ids = colnames(mat_df)
  mat = as.matrix(mat_df)
  rownames(mat) = "Urinary p75ECD"
  
  # Create annotations 
  annot_df = df_visit %>%
    group_by(ParticipantCode) %>%
    summarise(
      p75ECD_val = first(p75ECD),
      MotorSigns = first(MotorSigns),
      type = first(type)
    ) %>%
    left_join(change_df %>% select(ParticipantCode, motor_change), by = "ParticipantCode") %>%
    slice(match(participant_ids, ParticipantCode))
  
  # Define specific colors for the new status annotation
  status_colors = c(
    "increase" = "red", 
    "decrease" = "blue", 
    "stable"   = "grey90"
  )
  
  type_colors = group_colors <- c(
    "ALS" = "#F8766D",    # pink
    "PGMC" = "#7CAE00",   # greenish-blue
    "CTR" = "#00BFC4",    # turquoise
    "mimic" = "#C77CFF"   # purple
  )
  
  ha = HeatmapAnnotation(
    MotorSigns = annot_df$MotorSigns,
    MotorStatus = annot_df$motor_change, 
    DiseaseType = annot_df$type,
    col = list(
      MotorSigns = colorRamp2(c(0, 4, 8), c("#F7F7F7", "#BDBDBD", "#525252")),
      MotorStatus = status_colors,
      DiseaseType = type_colors
    ),
    annotation_name_side = "left"
  )
  
  # Heatmap
  ht = Heatmap(
    mat,
    name = "p75ECD",
    top_annotation = ha,
    col = colorRamp2(
      c(min(mat, na.rm = TRUE), median(mat, na.rm = TRUE), max(mat, na.rm = TRUE)),
      c("navy", "white", "firebrick")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 8, col = "black"), 
    row_names_gp = gpar(fontsize = 10),
    column_title = paste0("Urinary p75ECD - ", visit_name)
  )
  
  return(ht)
}


ht_V0 = create_heatmap("V0")
ht_V1 = create_heatmap("V1")

ht_V0

ht_V1

svglite("plots/boxplots_urine/heatmap_V0.svg", width = 13, height = 5)
draw(ht_V0,heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()
svglite("plots/boxplots_urine/heatmap_V1.svg", width = 13, height = 4)
draw(ht_V1,heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom")
dev.off()

