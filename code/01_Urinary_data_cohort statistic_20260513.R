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
library(emmeans)
library(coin)
library(fuzzyjoin)


merged_table00 = read_excel("results/clinico_molecular_final_table.xlsx")


### STEP 1: Quality checks, tube, patients, visits number

# how many IDs in the raw data 
length(merged_table00$Tube_number)
#317
# how many unique IDs in the raw data CSF
length(unique(merged_table00$Tube_number))
#317

# how many unique Patient IDs
length(unique(na.omit(merged_table00$ParticipantCode)))
#199

# how many tube numbers and patient IDs per V0 and V1 - including NA and SYMP
merged_table00 %>%
  group_by(Visit) %>%
  summarise(
    total_patients = n_distinct(ParticipantCode, na.rm = TRUE), # amount of patient IDs
    total_tubes = n_distinct(Tube_number, na.rm = TRUE), # amount of tube numbers
    total_patients_tubes_pairs = n_distinct(
      paste(ParticipantCode, Tube_number)[!is.na(ParticipantCode) & !is.na(Tube_number)])) # amount of pairs

# Visit total_patients total_tubes total_patients_tubes_pairs
# <chr>          <int>       <int>                      <int>
# 1 V0               199         199                        199
# 2 V1               118         118                        118

# patients at V0 and V1
patients_V0 <- merged_table00 %>%
  filter(Visit == "V0") %>%
  pull(ParticipantCode) %>%
  unique()
patients_V1 <- merged_table00 %>%
  filter(Visit == "V1") %>%
  pull(ParticipantCode) %>%
  unique() %>% 
  na.omit()

# check if all patients at V1 are at V0
all(patients_V1 %in% patients_V0)
# TRUE

merged_table00 %>%
  count(type, Visit)

# type  Visit     n
# <chr> <chr> <int>
# 1 ALS   V0       70
# 2 ALS   V1       31
# 3 CTR   V0       47
# 4 CTR   V1       37
# 5 NA    V0        1
# 6 NA    V1        1
# 7 PGMC  V0       53
# 8 PGMC  V1       35
# 9 SYMP  V0        8
# 10 SYMP  V1        5
# 11 mimic V0       20
# 12 mimic V1        9


##### STEP 2: Remove outliers, DE320, TR302 and TR310

merged_table <- merged_table00 %>%
  filter(!ParticipantCode %in% c("DE320", "TR302", "TR310"))

##### STEP 3: Check number of participant per visit only for the 4 groups we analyze here

merged_table_final_group <- merged_table %>%
  filter(type %in% c("CTR", "PGMC", "ALS", "mimic"))

# how many IDs in the raw data 
length(merged_table_final_group$Tube_number)
#297
# how many unique IDs in the raw data CSF
length(unique(merged_table_final_group$Tube_number))
#297

# how many unique Patient IDs
length(unique(na.omit(merged_table_final_group$ParticipantCode)))
#187 = 199 - 3 outliers - 8 SYMP - 1 NA, OK

# how many tube numbers and patient IDs per V0 and V1 - including NA and SYMP
merged_table_final_group %>%
  group_by(Visit) %>%
  summarise(
    total_patients = n_distinct(ParticipantCode, na.rm = TRUE), # amount of patient IDs
    total_tubes = n_distinct(Tube_number, na.rm = TRUE), # amount of tube numbers
    total_patients_tubes_pairs = n_distinct(
      paste(ParticipantCode, Tube_number)[!is.na(ParticipantCode) & !is.na(Tube_number)])) # amount of pairs

# Visit total_patients total_tubes total_patients_tubes_pairs
# <chr>          <int>       <int>                      <int>
# 1 V0               187         187                        187
# 2 V1               110         110                        110


# patients at V0 and V1
patients_V0 <- merged_table_final_group %>%
  filter(Visit == "V0") %>%
  pull(ParticipantCode) %>%
  unique()
patients_V1 <- merged_table_final_group %>%
  filter(Visit == "V1") %>%
  pull(ParticipantCode) %>%
  unique() %>% 
  na.omit()

# check if all patients at V1 are at V0
all(patients_V1 %in% patients_V0)
# TRUE

merged_table_final_group %>%
  count(type, Visit)

# type  Visit     n
# <chr> <chr> <int>
# 1 ALS   V0       67
# 2 ALS   V1       29
# 3 CTR   V0       47
# 4 CTR   V1       37
# 5 PGMC  V0       53
# 6 PGMC  V1       35
# 7 mimic V0       20
# 8 mimic V1        9

#### - Clinical statistic

IDs_urine_data <- merged_table_final_group %>%
  distinct(ParticipantCode, Visit, .keep_all = TRUE) %>%
  mutate(
    Country = case_when(
      grepl("TR", ParticipantCode) ~ "Turkey",
      grepl("CH", ParticipantCode) ~ "Switzerland",
      grepl("DE", ParticipantCode) ~ "Germany",
      grepl("SK", ParticipantCode) ~ "Slovakia",
      grepl("FR", ParticipantCode) ~ "France",
      grepl("IL", ParticipantCode) ~ "Israel",
      grepl("CZ", ParticipantCode) ~ "Czechia",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(ParticipantCode)

DT::datatable(IDs_urine_data)

overview_country_urine = IDs_urine_data %>%
  group_by(Visit,Country) %>%
  summarise(n_unique_patients = n_distinct(ParticipantCode),.groups = "drop") %>%
  arrange(Visit,desc(n_unique_patients))
overview_country_urine

# Visit Country     n_unique_patients
# <chr> <chr>                   <int>
# 1 V0    Turkey                     74
# 2 V0    Germany                    43
# 3 V0    Slovakia                   41
# 4 V0    Switzerland                15
# 5 V0    France                     14
# 6 V1    Turkey                     67
# 7 V1    Germany                    19
# 8 V1    Slovakia                   19
# 9 V1    Switzerland                 5


#group per country 
Type_per_country <-
  merged_table_final_group %>%
  group_by(type, Visit, Country) %>%
  summarise(nr_samples = n_distinct(ParticipantCode), .groups = "drop") %>%
  arrange(Visit, Country) %>%
  pivot_wider(
    names_from  = type,
    values_from = nr_samples,
    values_fill = 0
  )
Type_per_country

# Visit Country       ALS   CTR  PGMC mimic
# <chr> <chr>       <int> <int> <int> <int>
# 1 V0    France          4     1     9     0
# 2 V0    Germany-M      23     6     7     7
# 3 V0    Slovakia       24     9     2     6
# 4 V0    Switzerland     1     6     4     4
# 5 V0    Turkey         15    25    31     3
# 6 V1    Germany-M       9     4     4     2
# 7 V1    Slovakia        9     7     0     3
# 8 V1    Switzerland     0     2     2     1
# 9 V1    Turkey         11    24    29     3


# counts by mutation for PGMC only
mutation_PGMC_per_country <-
  merged_table_final_group %>%
  filter(type == "PGMC") %>%
  mutate(Visit = factor(Visit, levels = c("V0", "V1", "V2"))) %>%
  group_by(mutation, Visit, Country) %>%
  summarise(nr_samples = n_distinct(ParticipantCode), .groups = "drop") %>%
  arrange(Visit, Country) %>%
  pivot_wider(
    names_from  = mutation,
    values_from = nr_samples,
    values_fill = 0
  )
mutation_PGMC_per_country

# # A tibble: 8 ? 10
#    Visit Country     C9orf72  SOD1   FUS TARDBP  FIG4 UBQLN2 ANXA11 HSPB1
#    <fct> <chr>         <int> <int> <int>  <int> <int>  <int>  <int> <int>
# 1 V0    France            7     2     0      0     0      0      0     0
# 2 V0    Germany           3     1     2      1     0      0      0     0
# 3 V0    Slovakia          1     0     0      0     1      0      0     0
# 4 V0    Switzerland       2     1     0      0     0      1      0     0
# 5 V0    Turkey            6    13     0      5     1      1      3     2
# 6 V1    Germany           1     0     2      1     0      0      0     0
# 7 V1    Switzerland       2     0     0      0     0      0      0     0
# 8 V1    Turkey            6    12     0      5     1      1      3     1

# counts by mutation for ALS only
mutation_ALS_per_country <-
  merged_table_final_group %>%
  filter(type == "ALS") %>%
  group_by(mutation, Visit, Country) %>%
  summarise(nr_samples = n_distinct(ParticipantCode), .groups = "drop") %>%
  arrange(Visit, Country) %>%
  pivot_wider(
    names_from  = mutation,
    values_from = nr_samples,
    values_fill = 0
  )

mutation_ALS_per_country

# Visit Country     C9orf72  `NA`  SOD1  FIG4  NEK1 ANXA11 TARDBP
# <chr> <chr>         <int> <int> <int> <int> <int>  <int>  <int>
# 1 V0    France            1     3     0     0     0      0      0
# 2 V0    Germany-M         1    21     1     0     0      0      0
# 3 V0    Slovakia          1    21     0     1     1      0      0
# 4 V0    Switzerland       0     1     0     0     0      0      0
# 5 V0    Turkey            0    12     1     0     0      1      1
# 6 V1    Germany-M         1     7     1     0     0      0      0
# 7 V1    Slovakia          1     7     0     1     0      0      0
# 8 V1    Turkey            0     9     1     0     0      1      0

# Age statistic per group (unique urine ID) -> rather extract V0 than taking unique no?
merged_table_final_group_unique <- merged_table_final_group %>%
  distinct(PatientID, .keep_all = TRUE)

merged_table_final_group_unique %>%
  filter(type == "CTR") %>%
  skim(age)
# 47
# ?????? Variable type: numeric ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# skim_variable n_missing complete_rate mean   sd p0 p25 p50  p75 p100 hist 
# 1 age                   0             1 46.7 11.7 23  39  47 54.5   78 ???????????????

merged_table_final_group_unique %>%
  filter(type == "ALS") %>%
  skim(age)
# 67
# ?????? Variable type: numeric ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# skim_variable n_missing complete_rate mean   sd  p0 p25  p50 p75 p100 hist 
# 1 age                   0            1 60.5 12.0 28 56.2  62  67   86 ???????????????

merged_table_final_group_unique %>%
  filter(type == "PGMC") %>%
  skim(age)
# 53
# ?????? Variable type: numeric ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# skim_variable n_missing complete_rate mean   sd p0 p25 p50 p75 p100 hist 
# 1 age                   0             1 48.9 15.2 22  38  47  57   86 ???????????????

merged_table_final_group_unique %>%
  filter(type == "mimic") %>%
  skim(age)
# 20
# ?????? Variable type: numeric ????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
# skim_variable n_missing complete_rate   mean   sd p0  p25 p50  p75 p100 hist 
# 1 age                   0             1 54.8 14.7 23 45.5  58 63.8   77 ???????????????


# correction: test normality per disease group
shapiro.test(merged_table_final_group_unique %>% 
               filter(type == "CTR") %>%
               pull(age))
# p-value = 0.5875
shapiro.test(merged_table_final_group_unique %>% 
               filter(type == "PGMC") %>%
               pull(age))
# p-value = 0.05866
shapiro.test(merged_table_final_group_unique %>% 
               filter(type == "ALS") %>%
               pull(age))
# p-value = 0.01462, not normal
shapiro.test(merged_table_final_group_unique %>% 
               filter(type == "mimic") %>%
               pull(age))
# p-value = 0.2871

# test normality for all data set, not per group
shapiro.test(merged_table_final_group_unique$age)
# follow normality p-value = 0.127

# Test homogeneity of variance per group here then? 
leveneTest(age ~ type, data = merged_table_final_group_unique)
# p = 0.1997, variance are homogeneous

## better to use KW here
kruskal.test(age ~ type, data = merged_table_final_group_unique)
# p-value = 4.028e-08

# # -> ANOVA two sided
# anova_model <- aov(age ~ type, data = merged_table_final_group_unique)
# summary(anova_model)
# #  2.14e-07

# # Post hoc anlysis, Tukey automatically correct for multiple comparison -> adj P value
# TukeyHSD(anova_model)
# 
# # diff        lwr       upr     p adj
# # CTR-ALS    -13.683709 -20.202656 -7.164762 0.0000010  ****
# # mimic-ALS   -5.543284 -14.273370  3.186803 0.3553837
# # PGMC-ALS   -11.475359 -17.773707 -5.177011 0.0000271  ****
# # mimic-CTR    8.140426  -1.006704 17.287555 0.1000686
# # PGMC-CTR     2.208350  -4.656396  9.073096 0.8381869
# # PGMC-mimic  -5.932075 -14.923321  3.059170 0.3211091

# with Dunn's (quite concordant with ANOVA)
dunnTest(age ~ type, data = merged_table_final_group_unique,
         method="bonferroni")
# Comparison         Z      P.unadj        P.adj
# 1    ALS - CTR  5.350809 8.756179e-08 5.253707e-07
# 2  ALS - mimic  1.539417 1.237025e-01 7.422151e-01
# 3  CTR - mimic -2.344167 1.906964e-02 1.144179e-01
# 4   ALS - PGMC  4.853948 1.210276e-06 7.261659e-06
# 5   CTR - PGMC -0.627815 5.301251e-01 1.000000e+00
# 6 mimic - PGMC  1.905476 5.671817e-02 3.403090e-01

## Sex statistic per group (unique urine ID)

merged_table_final_group_unique %>%
  count(type, sex) %>%
  group_by(type) %>%
  mutate(percentage = round(100 * n / sum(n), 1))

#   type  sex       n percentage
#   <chr> <chr> <int>      <dbl>
# 1 ALS   F        26       38.8
# 2 ALS   M        41       61.2
# 3 CTR   F        21       44.7
# 4 CTR   M        26       55.3
# 5 PGMC  F        31       58.5
# 6 PGMC  M        22       41.5
# 7 mimic F         4       20  
# 8 mimic M        16       80  

# Test for more than 2 groups
table(merged_table_final_group_unique$sex, merged_table_final_group_unique$type)
chisq.test(table(merged_table_final_group_unique$sex, merged_table_final_group_unique$type))
#  p-value = 0.01908-> sex distribution significantly differs across the 4 groups


#  Post hoc: Since some expected counts are <5 (e.g., Mimic F = 4), Fisher's exact test is more robust:
pairwise_fisher_test(table(merged_table_final_group_unique$type,merged_table_final_group_unique$sex), p.adjust.method = "holm")

#  <chr>  <chr>  <int>   <dbl>  <dbl> <chr>       
# 1 ALS    CTR      114 0.566   0.566  ns          
# 2 ALS    mimic     87 0.18    0.54   ns          
# 3 ALS    PGMC     120 0.0429  0.214  ns          
# 4 CTR    mimic     67 0.0962  0.385  ns          
# 5 CTR    PGMC     100 0.229   0.54   ns          
# 6 mimic  PGMC      73 0.00399 0.0239 *  


# Mutation in the ALS group (urine cohort)

ALS_mutation_table <- merged_table_final_group_unique %>%
  filter(type == "ALS") %>%
  count(mutation) %>%
  mutate(
    percent = round(100 * n / sum(n), 1),
    summary = paste0(n, " (", percent, "%)")
  )

ALS_mutation_table
# mutation     n percent summary   
# <chr>    <int>   <dbl> <chr>     
# 1 ANXA11       1     1.5 1 (1.5%)  
# 2 C9orf72      3     4.5 3 (4.5%)  
# 3 FIG4         1     1.5 1 (1.5%)  
# 4 NEK1         1     1.5 1 (1.5%)  
# 5 SOD1         2     3   2 (3%)    
# 6 TARDBP       1     1.5 1 (1.5%)  
# 7 NA          58    86.6 58 (86.6%)


other_mutations_ALS <- ALS_mutation_table %>%
  filter(!mutation %in% c("C9orf72", "SOD1", "TARDBP"), !is.na(mutation)) %>%
  summarise(
    n = sum(n, na.rm = TRUE),
    percent = sum(percent, na.rm = TRUE),
    summary = paste0(sum(n, na.rm = TRUE), " (", round(sum(percent, na.rm = TRUE),1), "%)")
  )

other_mutations_ALS
#   
# n percent summary 
# <int>   <dbl> <chr>   
#   3     4.5 3 (4.5%)


#Mutation in tne the PGMC group (urine cohort)

PGMC_mutation_table <- merged_table_final_group_unique %>%
  filter(type == "PGMC") %>%
  count(mutation) %>%
  mutate(
    percent = round(100 * n / sum(n), 1),
    summary = paste0(n, " (", percent, "%)")
  )

PGMC_mutation_table

# mutation     n percent summary   
# <chr>    <int>   <dbl> <chr>     
#   1 ANXA11       3     5.7 3 (5.7%)  
# 2 C9orf72     19    35.8 19 (35.8%)
# 3 FIG4         2     3.8 2 (3.8%)  
# 4 FUS          2     3.8 2 (3.8%)  
# 5 HSPB1        2     3.8 2 (3.8%)  
# 6 SOD1        17    32.1 17 (32.1%)
# 7 TARDBP       6    11.3 6 (11.3%) 
# 8 UBQLN2       2     3.8 2 (3.8%)  


other_mutations_PGMC <- PGMC_mutation_table %>%
  filter(!mutation %in% c("C9orf72", "SOD1", "TARDBP"), !is.na(mutation)) %>%
  summarise(
    n = sum(n, na.rm = TRUE),
    percent = sum(percent, na.rm = TRUE),
    summary = paste0(sum(n, na.rm = TRUE), " (", round(sum(percent, na.rm = TRUE),1), "%)")
  )

other_mutations_PGMC

# n percent summary   
# <int>   <dbl> <chr>     
#   1    11    20.9 11 (20.9%)


# Site at onset statistic per group

ALS_site_onset <- merged_table_final_group_unique %>%
  filter(type == "ALS") %>%
  count(SiteDiseaseOnset) %>%
  mutate(
    percent = round(100 * n / sum(n), 1),
    summary = paste0(n, " (", percent, "%)")
  )

ALS_site_onset

# SiteDiseaseOnset     n percent summary   
# <chr>            <int>   <dbl> <chr>     
#   1 bulbar              18    26.9 18 (26.9%)
# 2 other                1     1.5 1 (1.5%)  
# 3 spinal              48    71.6 48 (71.6%)

Other_siteonset <- merged_table_final_group_unique %>%
  filter(type == "ALS", SiteDiseaseOnset == "other")
# Axial

# ALSFRS-r statistic per group and visit

merged_table_final_group %>%
  group_by(type, Visit) %>%
  summarise(
    mean_ALSFRS = mean(`Total ALSFRS-r`, na.rm = TRUE),
    sd_ALSFRS   = sd(`Total ALSFRS-r`, na.rm = TRUE),
    n = sum(!is.na(`Total ALSFRS-r`)),
    .groups = "drop"
  )

# type  Visit     mean_ALSFRS sd_ALSFRS     n
# <chr> <chr>       <dbl>     <dbl> <int>
# 1 ALS   V0           41.1     5.10     67
# 2 ALS   V1           32.6     9.84     29
# 3 CTR   V0           47.9     0.729    47
# 4 CTR   V1           48       0        37
# 5 PGMC  V0           47.6     1.85     53
# 6 PGMC  V1           46.5     5.19     35
# 7 mimic V0           45.0     4.61     20
# 8 mimic V1           45.3     3.81      9



write.xlsx(
  merged_table_final_group,
  file = "results/Prticipant_without_outliers_final.xlsx",
  sheetName = "without_outliers_final",
  overwrite = TRUE
)



