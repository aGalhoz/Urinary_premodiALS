library(readxl)
library(dplyr)
library(writexl)
library(tidyr)
library(ggplot2)
library(patchwork)


#### Correlation between p75ECD and NEFL NULISA values at V0

# NULISA data
NULISA_data = read_excel("data_NULISA/Updated_P004_NPQ_09042025.xlsx")
sample_ID_info <- read_excel("data_NULISA/Biospecimen Manifest Form_Banner Biomarker Program_premodiALS.xlsx", sheet = "Sample Information") %>% select("Sample ID","Optional Informarion (patient ID)")

# -> change NA2 info to DE107
sample_ID_info$`Optional Informarion (patient ID)` <- ifelse(sample_ID_info$`Optional Informarion (patient ID)` == "NA2",
                                                             "DE107",sample_ID_info$`Optional Informarion (patient ID)`)
# -> change samples of size 4 with TR116 to TR126 (writing mistake)
sample_ID_info$`Optional Informarion (patient ID)` <- ifelse(sample_ID_info$`Sample ID` %in% c("9772","9636","9626"),"TR126",
                                                             sample_ID_info$`Optional Informarion (patient ID)`)

# urinary data at V0 : plot_data in the previous script ?? No can not find it anymore... Urine data
urinary_data_V0 <- plot_data %>%
  filter(Visit == "V0") 


proteins <- c("NEFL", "MAPT", "FABP3", "pTau-217", "TAFA5", "NEFH", "pTau-181", "pTau-231", "IL18", "GFAP")


analyse_NULISA_protein <- function(protein,
                                   NULISA_data,
                                   sample_ID_info,
                                   urinary_data_V0) {
  
  # --------------------------------------------------
  # 1. Select protein
  # --------------------------------------------------
  
  data_protein <- NULISA_data %>%
    filter(
      Target == protein,
      SampleMatrixType %in% c("CSF", "PLASMA", "SERUM")
    )
  
  
  # --------------------------------------------------
  # 2. Add Participant IDs
  # --------------------------------------------------
  
  data_protein_IDs <- data_protein %>%
    left_join(
      sample_ID_info %>%
        rename(
          SampleName = `Sample ID`,
          ParticipantCode = `Optional Informarion (patient ID)`
        )
    )
  
  
  # --------------------------------------------------
  # 3. Create V0 dataset
  # --------------------------------------------------
  
  protein_V0 <- data_protein_IDs %>%
    select(
      ParticipantCode,
      SampleMatrixType,
      NPQ
    ) %>%
    mutate(
      Visit = "V0"
    ) %>%
    select(
      ParticipantCode,
      Visit,
      SampleMatrixType,
      NPQ
    )
  
  
  # --------------------------------------------------
  # 4. Wide format
  # --------------------------------------------------
  
  protein_wide <- protein_V0 %>%
    mutate(
      SampleMatrixType = case_when(
        SampleMatrixType == "PLASMA" ~ paste("Plasma", protein),
        SampleMatrixType == "SERUM"  ~ paste("Serum", protein),
        SampleMatrixType == "CSF"    ~ paste("CSF", protein),
        TRUE ~ SampleMatrixType
      )
    ) %>%
    select(
      ParticipantCode,
      SampleMatrixType,
      NPQ
    ) %>%
    pivot_wider(
      names_from = SampleMatrixType,
      values_from = NPQ
    )
  
  
  # --------------------------------------------------
  # 5. Merge with urinary p75ECD
  # --------------------------------------------------
  
  protein_urinary <- data_protein_IDs %>%
    full_join(
      urinary_data_V0 %>%
        select(
          ParticipantCode,
          `ng p75ECD/mg creatinine`,
          type
        ),
      by = "ParticipantCode"
    )
  
  
  # --------------------------------------------------
  # 6. Keep samples with both values
  # --------------------------------------------------
  
  protein_urinary_filter <- protein_urinary %>%
    filter(
      !is.na(NPQ),
      !is.na(`ng p75ECD/mg creatinine`)
    )
  
  
  # --------------------------------------------------
  # 7. Correlation plotting function
  # --------------------------------------------------
  
  plot_correlation <- function(data, fluid, protein_name, cohort_label) {
    
    df <- data %>%
      filter(SampleMatrixType == fluid)
    
    # Remove missing values
    df <- df %>%
      filter(
        !is.na(NPQ),
        !is.na(`ng p75ECD/mg creatinine`)
      )
    
    # Avoid cor.test error if too few observations
    if (nrow(df) < 3) {
      return(
        ggplot() +
          theme_void() +
          annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = paste0(
              fluid, "\n",
              cohort_label, "\n",
              "n < 3"
            ),
            size = 5
          )
      )
    }
    
    # Spearman correlation
    cor_res <- cor.test(
      df$NPQ,
      df$`ng p75ECD/mg creatinine`,
      method = "spearman"
    )
    
    rho <- unname(cor_res$estimate)
    p   <- cor_res$p.value
    
    label <- paste0(
      "\u03C1 = ", round(rho, 3),
      "\np = ", signif(p, 3),
      "\nn = ", nrow(df)
    )
    
    ggplot(
      df,
      aes(
        x = NPQ,
        y = `ng p75ECD/mg creatinine`
      )
    ) +
      geom_point(
        size = 2,
        alpha = 0.7
      ) +
      geom_smooth(
        method = "lm",
        se = TRUE,
        colour = "black",
        linetype = "dashed"
      ) +
      annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = label,
        hjust = -0.1,
        vjust = 1.1,
        size = 4
      ) +
      labs(
        x = paste0(fluid, " ", protein_name, " (NPQ)"),
        y = "Urinary p75ECD / mg creatinine",
        title = fluid
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        )
      )
  }
  
  
  # --------------------------------------------------
  # 8. Generate ALL-COHORT plots
  # --------------------------------------------------
  
  plots_all <- list(
    
    Plasma = plot_correlation(
      protein_urinary_filter,
      "PLASMA",
      protein,
      "All cohorts"
    ),
    
    Serum = plot_correlation(
      protein_urinary_filter,
      "SERUM",
      protein,
      "All cohorts"
    ),
    
    CSF = plot_correlation(
      protein_urinary_filter,
      "CSF",
      protein,
      "All cohorts"
    )
  )
  
  
  # --------------------------------------------------
  # 9. ALS-only dataset
  # --------------------------------------------------
  
  protein_urinary_filter_ALS <- protein_urinary_filter %>%
    filter(type == "ALS")
  
  
  # --------------------------------------------------
  # 10. Generate ALS-only plots
  # --------------------------------------------------
  
  plots_ALS <- list(
    
    Plasma = plot_correlation(
      protein_urinary_filter_ALS,
      "PLASMA",
      protein,
      "ALS only"
    ),
    
    Serum = plot_correlation(
      protein_urinary_filter_ALS,
      "SERUM",
      protein,
      "ALS only"
    ),
    
    CSF = plot_correlation(
      protein_urinary_filter_ALS,
      "CSF",
      protein,
      "ALS only"
    )
  )
  
  
  # --------------------------------------------------
  # 11. Combine into ONE PAGE
  # --------------------------------------------------
  
  combined_plot <- (
    
    plots_all$Plasma |
      plots_all$Serum |
      plots_all$CSF
    
  ) / (
    
    plots_ALS$Plasma |
      plots_ALS$Serum |
      plots_ALS$CSF
    
  ) +
    
    plot_annotation(
      title = paste0(
        protein,
        " vs urinary p75ECD at V0"
      )
    ) &
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 16
      )
    )
  
  
  # --------------------------------------------------
  # 12. Return everything
  # --------------------------------------------------
  
  return(
    list(
      protein_V0 = protein_V0,
      protein_wide = protein_wide,
      protein_urinary = protein_urinary,
      protein_urinary_filter = protein_urinary_filter,
      protein_urinary_filter_ALS = protein_urinary_filter_ALS,
      plots = plots_all,
      plots_ALS = plots_ALS,
      combined_plot = combined_plot
    )
  )
}



results <- lapply(
  proteins,
  analyse_NULISA_protein,
  NULISA_data = NULISA_data,
  sample_ID_info = sample_ID_info,
  urinary_data_V0 = urinary_data_V0
)

names(results) <- proteins


dir.create("NULISA_p75ECD_plots", showWarnings = FALSE)

for (protein in proteins) {
  
  ggsave(
    filename = paste0(
      "NULISA_p75ECD_plots/",
      protein,
      "_p75ECD_correlations.svg"
    ),
    plot = results[[protein]]$combined_plot,
    width = 12,
    height = 8,
    units = "in"
  )
}

# results$NEFL
# results$MAPT
# results$MAPT


########################  EXPORT TABLE WITH INTERESTING TARGETS

# List of proteins to include
proteins_to_extract <- c(
  "NEFL",
  "NEFH",
  "MAPT",
  "pTau-217",
  "pTau-181",
  "pTau-231",
  "TAFA5",
  "GFAP"
)


# Extract and combine all selected proteins
NULISA_table <- bind_rows(
  lapply(proteins_to_extract, function(protein) {
    
    results[[protein]]$protein_V0 %>%
      mutate(Protein = protein)
    
  })
) %>%
  filter(Visit == "V0") %>%
  mutate(
    SampleMatrixType = recode(
      SampleMatrixType,
      "PLASMA" = "Plasma",
      "SERUM"  = "Serum",
      "CSF"    = "CSF"
    ),
    Column = paste(SampleMatrixType, Protein, sep = "_")
  ) %>%
  select(ParticipantCode, Column, NPQ) %>%
  pivot_wider(
    names_from = Column,
    values_from = NPQ
  ) %>%
  arrange(ParticipantCode)

# View table
View(NULISA_table)

# Export to Excel
write_xlsx(NULISA_table, "NULISA_table.xlsx")
