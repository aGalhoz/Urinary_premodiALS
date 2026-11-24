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



### STEP 1: Directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# create directory for results
dir.create(file.path(getwd(),'results'), showWarnings = FALSE)
# create directory for plots
dir.create(file.path(getwd(),'plots'), showWarnings = FALSE)
dir.create(file.path(getwd(),'plots/boxplots_urine'), showWarnings = FALSE)



## Load tube and patient ID

# Loading data including:
#
# - All templates: information on the tube ID, rack, visit, country, type of fluid, position
# - Shipment of all boxes file: information on the patient ID, visit, position
# - Raw data: lipidomics expression with info on the tube ID for CSF and Plasma


### STEP 2: merge info of all patients

# The **first step** is to *merge* the information of the *"All templates"* and
# *"Shipment of all boxes"* files, and *create 3 datasets* (one for *plasma*, one for
# *CSF*, one for *urine*) with the information of *Patient ID, Position, Biofluid, Visit, Tube number, Country, Rack and Box*.

# Add visit function
make_visit_vec <- function(position, urine_id) {

  # Extract plate row letter (A-H)
  row_letter <- substr(position, 1, 1)

  # Default rule
  visit <- ifelse(row_letter %in% c("A", "C", "E", "G"), "V0", "V1")

  # Special rule for urine_4
  if (urine_id == 4) {
    visit[row_letter %in% c("C","D","E","F","G","H")] <- "V0"
  }

  visit
}

# -> urine
urine_1 = read_excel("Rack_5500275403_Urine Australia_1.xlsx",
                                sheet = "Aliquots") %>%
  dplyr::rename(Position = Position,
                Rack = `Rack-Name`,
                Tube_number = `Barcode`,
                Biofluid = Material) %>%
  select(Position,Biofluid,Tube_number,Rack)

urine_1$Visit <- make_visit_vec(urine_1$Position, urine_id = 1)

urine_2 = read_excel("Rack_5500275416_Urine Australia_2.xlsx",
                     sheet = "Aliquots") %>%
  dplyr::rename(Position = Position,
                Rack = `Rack-Name`,
                Tube_number = `Barcode`,
                Biofluid = Material) %>%
  select(Position,Biofluid,Tube_number,Rack)

urine_2$Visit <- make_visit_vec(urine_2$Position, urine_id = 2)

urine_3 = read_excel("Rack_5500275448_Urine Australia_3.xlsx",
                     sheet = "Aliquots") %>%
  dplyr::rename(Position = Position,
                Rack = `Rack-Name`,
                Tube_number = `Barcode`,
                Biofluid = Material) %>%
  select(Position,Biofluid,Tube_number,Rack)

urine_3$Visit <- make_visit_vec(urine_3$Position, urine_id = 3)

urine_4_1 = read_excel("Rack_5500275461_Urine Australia_4.xlsx",
                     sheet = "Aliquots") %>%
  dplyr::rename(Position = Position,
                Rack = `Rack-Name`,
                Tube_number = `Barcode`,
                Biofluid = Material) %>%
  select(Position,Biofluid,Tube_number,Rack)

urine_4_1$Visit <- make_visit_vec(urine_4_1$Position, urine_id = 4)

urine_4_2 <- read_excel("Rack_5500275461_Urine Australia_4.xlsx",
                        sheet = "Aliquots_4_2") %>%
  dplyr::rename(Position = Position,
                Rack = `Rack-Name`,
                Tube_number = `Barcode`,
                Biofluid = Material) %>%
  select(Position, Biofluid, Tube_number, Rack) %>%
  mutate(
    Tube_number = as.character(Tube_number),
    Tube_number = if_else(
      str_length(Tube_number) < 4,
      str_pad(Tube_number, width = 4, side = "left", pad = "0"),
      Tube_number
    )
  )

urine_4_2$Visit <- "V0"

#Merge 4_1 and 4_2

urine_4_1 <- urine_4_1 %>% mutate(Rack = as.character(Rack))
urine_4_2 <- urine_4_2 %>% mutate(Rack = as.character(Rack))

urine_4 <- bind_rows(urine_4_1, urine_4_2)


### STEP 3: get patient ID per box, position and visit
info_patientID_1 =  read_excel("premodiALS shipment all boxes list for omics analyses-2.xlsx",
                               sheet = "Box1") %>%
  dplyr::rename(Patient = `Patient ID`) %>%
  select(Patient,Position,Visit)

info_patientID_2 =  read_excel("premodiALS shipment all boxes list for omics analyses-2.xlsx",
                               sheet = "Box2") %>%
  dplyr::rename(Patient = `Patient ID`) %>%
  select(Patient,Position,Visit)

info_patientID_3 =  read_excel("premodiALS shipment all boxes list for omics analyses-2.xlsx",
                               sheet = "Box3") %>%
  dplyr::rename(Patient = `Patient ID`) %>%
  select(Patient,Position,Visit)


info_patientID_4_1 <- read_excel(
  "premodiALS shipment all boxes list for omics analyses-2.xlsx",
  sheet = "Box4"
) %>%
  dplyr::rename(Patient = `Patient ID`) %>%
  select(Patient, Position, Visit)
# keep only first 9 rows
info_patientID_4_1 <- info_patientID_4_1[c(1:9, 13:21), ]


info_patientID_4_2 <- read_excel(
  "premodiALS shipment all boxes list for omics analyses-2.xlsx",
  sheet = "unscanned"
) |>
  setNames(c("Patient","Position","Biofluid","Visit",
             "Tube_nr","Country","Rack")) |>
  mutate(
    Tube_nr = str_pad(as.character(Tube_nr), width = 4, pad = "0")
  ) |>
  filter(
    str_trim(str_to_lower(Country)) == "australia"
  ) |>
  select(Patient, Position, Visit)


info_patientID_4 <- bind_rows(info_patientID_4_1, info_patientID_4_2)


# Urine, match info
merge_urine_1 <- urine_1 %>% left_join(info_patientID_1) %>%
  select(Patient,Position,Biofluid,Visit,Tube_number,Rack) %>%
  mutate(Box = "Box 1")

merge_urine_2 <- urine_2 %>% left_join(info_patientID_2) %>%
  select(Patient,Position,Biofluid,Visit,Tube_number,Rack) %>%
  mutate(Box = "Box 2")

merge_urine_3 <- urine_3 %>% left_join(info_patientID_3) %>%
  select(Patient,Position,Biofluid,Visit,Tube_number,Rack) %>%
  mutate(Box = "Box 3")

merge_urine_4 <- urine_4 %>% left_join(info_patientID_4) %>%
  select(Patient,Position,Biofluid,Visit,Tube_number,Rack) %>%
  mutate(Box = "Box 4")


merge_urine_all = do.call("rbind",list(merge_urine_1,
                                            merge_urine_2,
                                            merge_urine_3,
                                            merge_urine_4))

length(unique(merge_urine_all$Patient))
#199




### STEP 4: get raw data
# p75ECD

raw_data_p75ECD =  read_excel("2025 Premodi Germany results summary p75ECD_creatinine_neopterin.xlsx", 
                             sheet = "Data summary - p75ECD Urine") 

raw_data_p75ECD_values <- raw_data_p75ECD %>%
  select(3, 4, 5, 33) %>%
  { .[-2, ] } %>%  # remove row 2 safely
  setNames(c("Tube_number", "Rack", "Position", "ng p75ECD/mg creatinine")) %>%
  filter(
    if_any(everything(), ~ !is.na(.) & . != ""),
    Tube_number != "nothing"
  )



# Neopterin

raw_data_Neopterin = read_excel("2025 Premodi Germany results summary p75ECD_creatinine_neopterin.xlsx", 
                                sheet = "Data summary - Neopterin Urine") 

raw_data_Neopterin_values <- raw_data_Neopterin %>%
  select(3, 4, 5, 34) %>%
  { .[-2, ] } %>%  # remove row 2 safely
  setNames(c("Tube_number", "Rack", "Position", "umol neopterin/ mol creatinine")) %>%
  filter(
    if_any(everything(), ~ !is.na(.) & . != ""),
    Tube_number != "nothing"
  )

# final-table

final_table_urine <- merge_urine_all %>%
  left_join(raw_data_Neopterin_values, by = c("Tube_number", "Rack", "Position")) %>%
  left_join(raw_data_p75ECD_values, by = c("Tube_number", "Rack", "Position")) %>%
  select(
    Patient, Position, Biofluid, Visit, Tube_number, Rack,
    `umol neopterin/ mol creatinine`,
    `ng p75ECD/mg creatinine`
  )

final_table_urine <- final_table_urine %>%
  rename(ParticipantCode = Patient)


##### SOLVING PATIENT ID eror, biofluid ID (only name) CH301 and CH303 were inverted

final_table_urine$ParticipantCode <- final_table_urine$ParticipantCode %>%
  dplyr::recode(
    "CH301" = "TEMP",
    "CH303" = "CH301"
  )

final_table_urine$ParticipantCode[final_table_urine$ParticipantCode  == "TEMP"]  <- "CH303"

write.xlsx(
  final_table_urine,
  file = "results/final_table_urine.xlsx",
  sheetName = "final_table_urine",
  overwrite = TRUE
)

final_table_urine = read_excel("results/final_table_urine.xlsx")


### STEP 5: Quality checks, tube, patients, visits number

# how many IDs in the raw data 
length(final_table_urine$Tube_number)
#317
# how many unique IDs in the raw data CSF
length(unique(final_table_urine$Tube_number))
#317

# how many unique Patient IDs
length(unique(na.omit(final_table_urine$ParticipantCode)))
#199

# how many tube numbers and patient IDs per V0 and V1
final_table_urine %>%
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
patients_V0 <- final_table_urine %>%
  filter(Visit == "V0") %>%
  pull(ParticipantCode) %>%
  unique()
patients_V1 <- final_table_urine %>%
  filter(Visit == "V1") %>%
  pull(ParticipantCode) %>%
  unique() %>% 
  na.omit()

# check if all patients at V1 are at V0
all(patients_V1 %in% patients_V0)
# TRUE



#### Step 6: merge with participant data from eCRF 

# Output file from the "initiation" script 

all_participants_ID_visits = read_excel("results/all_participants_ID_visits.xlsx")


#### - Urine (the patient IDs existing in the raw data)

merged_table <- inner_join(
  final_table_urine,
  all_participants_ID_visits,
  by = c("ParticipantCode", "Visit")
)

# missing_rows <- anti_join(
#   final_table_urine,
#   all_participants_ID_visits,
#   by = c("ParticipantCode", "Visit")
# )
# 
# missing_rows
# 0

write.xlsx(
  merged_table,
  file = "results/clinico_molecular_final_table.xlsx",
  sheetName = "clinico_molecular_final_table",
  overwrite = TRUE
)


