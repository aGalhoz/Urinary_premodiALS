library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(svglite)
library(readxl)

###################
# HEATMAP
# V0 p75ECD
# TotalMotorSigns
###################

# --------------------------------------------------
# Prepare data
# --------------------------------------------------

heatmap_df = Urine_data %>%
  select(
    ParticipantCode,
    Visit,
    TotalMotorSigns,
    `ng p75ECD/mg creatinine`,
    type
  ) %>%
  mutate(
    ParticipantCode = as.character(ParticipantCode),
    Visit = as.character(Visit),
    TotalMotorSigns = as.numeric(TotalMotorSigns),
    p75ECD = as.numeric(`ng p75ECD/mg creatinine`)
  ) %>%
  filter(
    type == "PGMC"
  )


# --------------------------------------------------
# V0 p75ECD
# --------------------------------------------------

p75ECD_V0 = heatmap_df %>%
  filter(
    Visit == "V0",
    !is.na(p75ECD)
  ) %>%
  group_by(ParticipantCode) %>%
  summarise(
    p75ECD = first(p75ECD),
    .groups = "drop"
  ) %>%
  arrange(desc(p75ECD))


# --------------------------------------------------
# V1 p75ECD
# --------------------------------------------------

p75ECD_V1 = heatmap_df %>%
  filter(
    Visit == "V1"
  ) %>%
  group_by(ParticipantCode) %>%
  summarise(
    p75ECD_V1 = first(p75ECD),
    .groups = "drop"
  )

# --------------------------------------------------
# Participant order based on V0 p75ECD
# --------------------------------------------------

participant_ids <- p75ECD_V0 %>%
  pull(ParticipantCode)


# --------------------------------------------------
# Cohort / participant type
# --------------------------------------------------

participant_type <- heatmap_df %>%
  select(ParticipantCode, type) %>%
  distinct() %>%
  filter(ParticipantCode %in% participant_ids) %>%
  right_join(
    data.frame(ParticipantCode = participant_ids),
    by = "ParticipantCode"
  ) %>%
  pull(type)

type_colors <- c(
  "ALS"   = "#F8766D",
  "PGMC"  = "#7CAE00",
  "CTR"   = "#00BFC4",
  "mimic" = "#C77CFF"
)


# Import data

all_participants_ID_visits <- read_excel(
  "results/all_participants_ID_visits.xlsx"
) %>%
  filter(
    type %in% c("ALS", "PGMC", "CTR", "mimic")
  )

# --------------------------------------------------
# Participant group + mutation
# --------------------------------------------------

participant_info <- data.frame(
  ParticipantCode = participant_ids
) %>%
  left_join(
    all_participants_ID_visits %>%
      select(
        ParticipantCode,
        type,
        mutation
      ) %>%
      mutate(
        ParticipantCode = as.character(ParticipantCode),
        mutation = as.character(mutation)
      ) %>%
      filter(!is.na(type)) %>%
      group_by(ParticipantCode) %>%
      summarise(
        type = first(type),
        mutation = ifelse(
          all(is.na(mutation)),
          NA_character_,
          first(na.omit(mutation))
        ),
        .groups = "drop"
      ),
    by = "ParticipantCode"
  )

participant_type <- participant_info$type
participant_mutation <- participant_info$mutation

participant_mutation[is.na(participant_mutation)] <- ""

# --------------------------------------------------
# Add V1 p75ECD in the same V0-based participant order
# --------------------------------------------------

p75ECD_data = data.frame(
  ParticipantCode = participant_ids
) %>%
  left_join(
    p75ECD_V0 %>%
      rename(p75ECD_V0 = p75ECD),
    by = "ParticipantCode"
  ) %>%
  left_join(
    p75ECD_V1,
    by = "ParticipantCode"
  )


# --------------------------------------------------
# Create p75ECD matrix: V1 above V0
# --------------------------------------------------

mat = rbind(
  "Urinary p75ECD - V1" = p75ECD_data$p75ECD_V1,
  "Urinary p75ECD - V0" = p75ECD_data$p75ECD_V0
)

colnames(mat) = participant_ids

# --------------------------------------------------
# Motor signs + phenoconversion data
# --------------------------------------------------

motor_data = all_participants_ID_visits %>%
  select(
    ParticipantCode,
    Visit,
    TotalMotorSigns,
    DiseaseOnset,
    `Disease duration`,
    mutation
  ) %>%
  mutate(
    ParticipantCode = as.character(ParticipantCode),
    Visit = as.character(Visit),
    TotalMotorSigns = as.numeric(TotalMotorSigns),
    `Disease duration` = as.numeric(`Disease duration`),
    mutation = as.character(mutation)
  ) %>%
  filter(
    ParticipantCode %in% participant_ids,
    Visit %in% c("V0", "V1", "V2", "V3")
  )



# --------------------------------------------------
# TotalMotorSigns V0
# --------------------------------------------------

motor_V0 = motor_data %>%
  filter(Visit == "V0") %>%
  group_by(ParticipantCode) %>%
  summarise(
    TotalMotorSigns_V0 = first(TotalMotorSigns),
    .groups = "drop"
  )


# --------------------------------------------------
# TotalMotorSigns V1
# --------------------------------------------------

motor_V1 = motor_data %>%
  filter(Visit == "V1") %>%
  group_by(ParticipantCode) %>%
  summarise(
    TotalMotorSigns_V1 = first(TotalMotorSigns),
    .groups = "drop"
  )


# --------------------------------------------------
# TotalMotorSigns V2
# --------------------------------------------------

motor_V2 = motor_data %>%
  filter(Visit == "V2") %>%
  group_by(ParticipantCode) %>%
  summarise(
    TotalMotorSigns_V2 = first(TotalMotorSigns),
    .groups = "drop"
  )


# --------------------------------------------------
# TotalMotorSigns V3
# --------------------------------------------------

motor_V3 = motor_data %>%
  filter(Visit == "V3") %>%
  group_by(ParticipantCode) %>%
  summarise(
    TotalMotorSigns_V3 = first(TotalMotorSigns),
    .groups = "drop"
  )


# --------------------------------------------------
# Combine TotalMotorSigns
# in the same order as V0 p75ECD
# --------------------------------------------------

annot_df = data.frame(
  ParticipantCode = participant_ids
) %>%
  left_join(
    motor_V0,
    by = "ParticipantCode"
  ) %>%
  left_join(
    motor_V1,
    by = "ParticipantCode"
  ) %>%
  left_join(
    motor_V2,
    by = "ParticipantCode"
  ) %>%
  left_join(
    motor_V3,
    by = "ParticipantCode"
  )


# Check values
print(annot_df)

# --------------------------------------------------
# TotalMotorSigns color scale
# --------------------------------------------------

motor_colors = colorRamp2(
  c(0, 4, 8),
  c("#FFFDEB", "#FEC44F", "#D95F0E")
)


# --------------------------------------------------
# PHENOCONVERSION
# --------------------------------------------------

phenoconversion_data = all_participants_ID_visits %>%
  filter(
    type == "PGMC",
    DiseaseOnset == "ALS",
    `Disease duration` > 0,
    ParticipantCode %in% participant_ids
  ) %>%
  group_by(ParticipantCode) %>%
  slice_min(
    `Disease duration`,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()


# --------------------------------------------------
# Determine phenoconversion visit
# --------------------------------------------------

phenoconversion_visit = phenoconversion_data %>%
  select(
    ParticipantCode,
    Visit
  ) %>%
  mutate(
    Visit = as.character(Visit)
  )

print(phenoconversion_visit)

# --------------------------------------------------
# Participants phenoconverting at each transition
# --------------------------------------------------

phenoconv_V0_V1 <- phenoconversion_visit %>%
  filter(Visit == "V1") %>%
  pull(ParticipantCode)

phenoconv_V1_V2 <- phenoconversion_visit %>%
  filter(Visit == "V2") %>%
  pull(ParticipantCode)

phenoconv_V2_V3 <- phenoconversion_visit %>%
  filter(Visit == "V3") %>%
  pull(ParticipantCode)


# Convert to TRUE/FALSE vectors in participant_ids order

phenoconv_V0_V1 <- participant_ids %in% phenoconv_V0_V1

phenoconv_V1_V2 <- participant_ids %in% phenoconv_V1_V2

phenoconv_V2_V3 <- participant_ids %in% phenoconv_V2_V3


participant_type <- data.frame(
  ParticipantCode = participant_ids
) %>%
  left_join(
    all_participants_ID_visits %>%
      select(ParticipantCode, type) %>%
      filter(!is.na(type)) %>%
      distinct(ParticipantCode, .keep_all = TRUE),
    by = "ParticipantCode"
  ) %>%
  pull(type)

# --------------------------------------------------
# Phenoconversion annotation function
# --------------------------------------------------

ha = HeatmapAnnotation(
  
  Mutation = anno_text(
    participant_mutation,
    gp = gpar(
      fontsize = 8,
      fontface = "plain"
    ),
    rot = 90,
    just = "center",
    height = unit(8, "mm")
  ),
  
  Group = participant_type,
  
  TotalMotorSigns_V3 = annot_df$TotalMotorSigns_V3,
  TotalMotorSigns_V2 = annot_df$TotalMotorSigns_V2,
  TotalMotorSigns_V1 = annot_df$TotalMotorSigns_V1,
  TotalMotorSigns_V0 = annot_df$TotalMotorSigns_V0,
  
  col = list(
    Group = type_colors,
    TotalMotorSigns_V3 = motor_colors,
    TotalMotorSigns_V2 = motor_colors,
    TotalMotorSigns_V1 = motor_colors,
    TotalMotorSigns_V0 = motor_colors
  ),
  
  na_col = "white",
  border = FALSE,
  
  annotation_name_side = "right",
  
  annotation_name_gp = gpar(
    fontsize = 12,
    col = "black"
  ),
  
  show_annotation_name = TRUE,
  
  show_legend = c(
    FALSE,
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    FALSE
  ),
  
  annotation_legend_param = list(
    Group = list(
      title = "Group"
    ),
    TotalMotorSigns_V3 = list(
      title = "TotalMotorSigns"
    )
  )
)

# --------------------------------------------------
# phenoconvertion legend
# --------------------------------------------------

phenoconversion_legend = Legend(
  title = NULL,
  labels = "Phenoconversion",
  type = "points",
  pch = 16,
  size = unit(3, "mm"),
  legend_gp = gpar(col = "black"),
  labels_gp = gpar(fontface = "bold")
)

# --------------------------------------------------
# Heatmap
# --------------------------------------------------

ht_V0 = Heatmap(
  
  mat,
  
  name = "p75ECD",
  
  top_annotation = ha,
  
  column_title = "Urinary p75ECD - PGMC\n",
  
  column_title_gp = gpar(
    fontsize = 14
  ),
  
  column_title_side = "top",
  
  column_title_rot = 0,
  
  col = colorRamp2(
    c(
      min(mat, na.rm = TRUE),
      median(mat, na.rm = TRUE),
      max(mat, na.rm = TRUE)
    ),
    c("navy", "#D9D9D9", "firebrick")
  ),
  
  na_col = "white",
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  column_names_gp = gpar(
    fontsize = 8,
    col = "black"
  ),
  
  row_names_gp = gpar(
    fontsize = 12
  )
)

# --------------------------------------------------
# Display
# --------------------------------------------------

# ht_V0


# --------------------------------------------------
# Save
# --------------------------------------------------

svglite(
  "plots/p75ECD/heatmap_V0_PGMC_TotalMotorSigns_V0_V1.svg",
  width = 13,
  height = 4.5
)

draw(
  ht_V0,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  annotation_legend_list = list(phenoconversion_legend),
  merge_legends = TRUE
)

decorate_annotation("Mutation", {
  grid.text(
    "Mutation",
    x = unit(1, "npc") + unit(9, "mm"),
    y = unit(0.5, "npc"),
    just = "center",
    gp = gpar(
      fontsize = 12,
      col = "black"
    )
  )
})

# --------------------------------------------------
# Add phenoconversion dots on MotorSigns boundaries
# --------------------------------------------------

# V0 -> V1
decorate_annotation("TotalMotorSigns_V0", {
  
  grid.clip("off")
  
  idx <- which(phenoconv_V0_V1)
  
  if (length(idx) > 0) {
    grid.points(
      x = unit(idx, "native"),
      y = unit(
        rep(1, length(idx)),
        "npc"
      ),
      pch = 16,
      size = unit(2.5, "mm"),
      gp = gpar(
        col = "black",
        fill = "black"
      )
    )
  }
})


# V1 -> V2
decorate_annotation("TotalMotorSigns_V1", {
  
  grid.clip("off")
  
  idx <- which(phenoconv_V1_V2)
  
  if (length(idx) > 0) {
    grid.points(
      x = unit(idx, "native"),
      y = unit(
        rep(1, length(idx)),
        "npc"
      ),
      pch = 16,
      size = unit(2.5, "mm"),
      gp = gpar(
        col = "black",
        fill = "black"
      )
    )
  }
})


# V2 -> V3
decorate_annotation("TotalMotorSigns_V2", {
  
  grid.clip("off")
  
  idx <- which(phenoconv_V2_V3)
  
  if (length(idx) > 0) {
    grid.points(
      x = unit(idx, "native"),
      y = unit(
        rep(1, length(idx)),
        "npc"
      ),
      pch = 16,
      size = unit(2.5, "mm"),
      gp = gpar(
        col = "black",
        fill = "black"
      )
    )
  }
})

dev.off()












# Similar graph but ploting and organizing according to delta V1-V0 p75ECD rather than V0 and V1 


###################
# HEATMAP
# Delta p75ECD V1 - V0
# Participants ordered by delta
###################

# --------------------------------------------------
# Calculate delta p75ECD V1 - V0
# Remove participants without V1 p75ECD
# --------------------------------------------------

delta_p75ECD_data = p75ECD_data %>%
  filter(
    !is.na(p75ECD_V0),
    !is.na(p75ECD_V1)
  ) %>%
  mutate(
    delta_p75ECD = p75ECD_V1 - p75ECD_V0
  ) %>%
  arrange(desc(delta_p75ECD))

# --------------------------------------------------
# New participant order based on delta p75ECD
# Higher -> lower
# --------------------------------------------------

participant_ids_delta <- delta_p75ECD_data %>%
  pull(ParticipantCode)



# --------------------------------------------------
# Reorder delta data
# --------------------------------------------------

delta_p75ECD_data = delta_p75ECD_data %>%
  slice(match(
    participant_ids_delta,
    ParticipantCode
  ))


# --------------------------------------------------
# Delta p75ECD matrix
# --------------------------------------------------

delta_mat = matrix(
  delta_p75ECD_data$delta_p75ECD,
  nrow = 1
)

rownames(delta_mat) = "?? p75ECD (V1???V0)"
colnames(delta_mat) = participant_ids_delta


# --------------------------------------------------
# Reorder motor-sign annotations according to
# delta p75ECD participant order
# --------------------------------------------------

annot_df_delta = data.frame(
  ParticipantCode = participant_ids_delta
) %>%
  left_join(
    motor_V0,
    by = "ParticipantCode"
  ) %>%
  left_join(
    motor_V1,
    by = "ParticipantCode"
  ) %>%
  left_join(
    motor_V2,
    by = "ParticipantCode"
  ) %>%
  left_join(
    motor_V3,
    by = "ParticipantCode"
  )


# --------------------------------------------------
# Reorder participant type
# --------------------------------------------------

participant_type_delta = data.frame(
  ParticipantCode = participant_ids_delta
) %>%
  left_join(
    all_participants_ID_visits %>%
      select(ParticipantCode, type) %>%
      filter(!is.na(type)) %>%
      distinct(ParticipantCode, .keep_all = TRUE),
    by = "ParticipantCode"
  ) %>%
  pull(type)

# --------------------------------------------------
# Reorder mutation according to delta p75ECD
# --------------------------------------------------

participant_mutation_delta = data.frame(
  ParticipantCode = participant_ids_delta
) %>%
  left_join(
    participant_info %>%
      select(ParticipantCode, mutation),
    by = "ParticipantCode"
  ) %>%
  pull(mutation)

participant_mutation_delta[is.na(participant_mutation_delta)] <- ""

# --------------------------------------------------
# Reorder phenoconversion vectors
# --------------------------------------------------

phenoconv_V0_V1_delta_ids <- phenoconversion_visit %>%
  filter(Visit == "V1") %>%
  pull(ParticipantCode)

phenoconv_V1_V2_delta_ids <- phenoconversion_visit %>%
  filter(Visit == "V2") %>%
  pull(ParticipantCode)

phenoconv_V2_V3_delta_ids <- phenoconversion_visit %>%
  filter(Visit == "V3") %>%
  pull(ParticipantCode)


phenoconv_V0_V1_delta <- participant_ids_delta %in%
  phenoconv_V0_V1_delta_ids

phenoconv_V1_V2_delta <- participant_ids_delta %in%
  phenoconv_V1_V2_delta_ids

phenoconv_V2_V3_delta <- participant_ids_delta %in%
  phenoconv_V2_V3_delta_ids


# --------------------------------------------------
# Delta p75ECD color scale
# --------------------------------------------------

delta_range <- max(
  abs(delta_mat),
  na.rm = TRUE
)

delta_colors = colorRamp2(
  c(
    -delta_range,
    0,
    delta_range
  ),
  c(
    "navy",
    "#D9D9D9",
    "firebrick"
  )
)


# --------------------------------------------------
# Annotation
# --------------------------------------------------

ha_delta = HeatmapAnnotation(
  
  Mutation = anno_text(
    participant_mutation_delta,
    gp = gpar(
      fontsize = 8,
      fontface = "plain"
    ),
    rot = 90,
    just = "center",
    height = unit(8, "mm")
  ),
  
  Group = participant_type_delta,
  
  TotalMotorSigns_V3 = annot_df_delta$TotalMotorSigns_V3,
  TotalMotorSigns_V2 = annot_df_delta$TotalMotorSigns_V2,
  TotalMotorSigns_V1 = annot_df_delta$TotalMotorSigns_V1,
  TotalMotorSigns_V0 = annot_df_delta$TotalMotorSigns_V0,
  
  col = list(
    
    Group = type_colors,
    
    TotalMotorSigns_V3 = motor_colors,
    TotalMotorSigns_V2 = motor_colors,
    TotalMotorSigns_V1 = motor_colors,
    TotalMotorSigns_V0 = motor_colors
  ),
  
  na_col = "white",
  
  border = FALSE,
  
  annotation_name_side = "right",
  
  annotation_name_gp = gpar(
    fontsize = 12,
    col = "black"
  ),
  
  show_annotation_name = TRUE,
  
  show_legend = c(
    FALSE,
    TRUE,
    TRUE,
    FALSE,
    FALSE,
    FALSE
  ),
  
  annotation_legend_param = list(
    
    Group = list(
      title = "Group"
    ),
    
    TotalMotorSigns_V3 = list(
      title = "TotalMotorSigns"
    )
  )
)

# --------------------------------------------------
# Phenoconversion legend
# --------------------------------------------------

phenoconversion_legend_delta = Legend(
  title = NULL,
  labels = "Phenoconversion",
  type = "points",
  pch = 16,
  size = unit(3, "mm"),
  legend_gp = gpar(col = "black"),
  labels_gp = gpar(fontface = "bold")
)


# --------------------------------------------------
# Heatmap
# --------------------------------------------------

ht_delta = Heatmap(
  
  delta_mat,
  
  name = "?? p75ECD",
  
  top_annotation = ha_delta,
  
  col = delta_colors,
  
  na_col = "white",
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  column_names_gp = gpar(
    fontsize = 8,
    col = "black"
  ),
  
  row_names_gp = gpar(
    fontsize = 12
  ),
  
  column_title = "Urinary p75ECD ?? (V1???V0) - PGMC\n"
)


# --------------------------------------------------
# Save
# --------------------------------------------------

svglite(
  "plots/p75ECD/heatmap_delta_p75ECD_V1_V0_PGMC_TotalMotorSigns.svg",
  width = 13,
  height = 4.5
)


draw(
  ht_delta,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  annotation_legend_list = list(
    phenoconversion_legend_delta
  ),
  merge_legends = TRUE
)

decorate_annotation("Mutation", {
  grid.text(
    "Mutation",
    x = unit(1, "npc") + unit(9, "mm"),
    y = unit(0.5, "npc"),
    just = "center",
    gp = gpar(
      fontsize = 12,
      col = "black"
    )
  )
})


# --------------------------------------------------
# Add phenoconversion dots
# --------------------------------------------------

# V0 -> V1
decorate_annotation("TotalMotorSigns_V0", {
  
  grid.clip("off")
  
  idx <- which(phenoconv_V0_V1_delta)
  
  if (length(idx) > 0) {
    grid.points(
      x = unit(idx, "native"),
      y = unit(
        rep(1, length(idx)),
        "npc"
      ),
      pch = 16,
      size = unit(2.5, "mm"),
      gp = gpar(
        col = "black",
        fill = "black"
      )
    )
  }
})


# V1 -> V2
decorate_annotation("TotalMotorSigns_V1", {
  
  grid.clip("off")
  
  idx <- which(phenoconv_V1_V2_delta)
  
  if (length(idx) > 0) {
    grid.points(
      x = unit(idx, "native"),
      y = unit(
        rep(1, length(idx)),
        "npc"
      ),
      pch = 16,
      size = unit(2.5, "mm"),
      gp = gpar(
        col = "black",
        fill = "black"
      )
    )
  }
})


# V2 -> V3
decorate_annotation("TotalMotorSigns_V2", {
  
  grid.clip("off")
  
  idx <- which(phenoconv_V2_V3_delta)
  
  if (length(idx) > 0) {
    grid.points(
      x = unit(idx, "native"),
      y = unit(
        rep(1, length(idx)),
        "npc"
      ),
      pch = 16,
      size = unit(2.5, "mm"),
      gp = gpar(
        col = "black",
        fill = "black"
      )
    )
  }
})


dev.off()