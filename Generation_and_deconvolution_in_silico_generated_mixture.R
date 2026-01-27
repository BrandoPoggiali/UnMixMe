# Project: Age prediction in DNAm mixture.
# Script Title: Generation of in-silico DNAm profiles with different DNA ratios, deconvolution of them and age prediction.
# Description: The script has the aim of in-silico prove that feasibility of the method for DNAm deconvolution and age prediction.
# The first step will be the creation of in-silico DNAm mixture containing different DNAm ratios.
# Author: Brando Poggiali
# Date: 25-04-2024


#### 1.Setting libraries and paths. ------------------------------------------------
BiocManager::install("")
install.packages("") 


#Load package
library(methylclock)
library(readxl)
library(writexl)
library(tibble)
library(tidyr)
library(tidyverse)
library(ggplot2)

#Set paths
#idat_dir <- "N:/projects/age_epic_array_blood_stains_4103/local_data"
EPIC_bloodstain_path <- "/mnt/ngs/projects/age_epic_array_blood_stains_4103/users/lfw156"
annotation_files_path <- "/home/ri-domain.local/lfw156/EPIC_annotations_files/"
brando_path <- "/mnt/ngs/projects/age_prediction_EPICv2_DNA_mixture/users/lfw156"
results_path <- "/mnt/ngs/projects/age_prediction_EPICv2_DNA_mixture/users/lfw156/Results/4_impact_of_variables/"

#### 2. Upload saved data (If you do not need to upload saved data skip this step)--------------------------------------------
beta_values_autosomal_cgch_noNAs <- readRDS(file.path(EPIC_bloodstain_path, "/Data/beta_values_autosomal_cgch_noNAs_collps_bloodstains_wide_06-02-2024.rds"))
beta_values_autosomal_cgch_noNAs <- as.data.frame(beta_values_autosomal_cgch_noNAs)
cg_ch_probes <- substr(rownames(beta_values_autosomal_cgch_noNAs),1,2) == "cg"
beta_values_autosomal_cg_noNAs <- beta_values_autosomal_cgch_noNAs[cg_ch_probes,]
metadata_bloodstains <- read.table(paste0(EPIC_bloodstain_path, "/Metadata/metadata_bloodstains_tidied_26-01-2024.tsv")) 
predicted_age <- read_xlsx(paste0(EPIC_bloodstain_path,'/Results/4_methylclok_age_prediction/Predicted_age_methylclock_07-02-2024.xlsx'))

betas_autosomal_cpgs_noNAs_2021 <- readRDS("/mnt/ngs/projects/age_prediction_EPICv2_DNA_mixture/users/lfw156/Data/beta_values_longitudinal_study_2021_22-11-2024.rds")
predicted_age_metadata_2021 <- read_xlsx('/mnt/ngs/projects/age_prediction_EPICv2_DNA_mixture/users/lfw156/Data/Metadata_and_predicted_age_longitudinal_study_2021_22-11-2024.xlsx')

df <- readRDS(paste0(results_path, "Age_prediction_STR_ratio_err_simulations_Longitudinal_data_09-01-2025.rds"))    

df_long <- readRDS(paste0(results_path, "Age_prediction_Age_gap_simulations_Longitudinal_data_22-07-2025.rds"))    

#Upload datasets
df <- readRDS(paste0(results_path, "Age_prediction_tecnology_err_10_simulations_Longitudinal_data_22-07-2025.rds"))

#### 3. Generation and Deconvolution of Mixture--------------------------------------------

mixture_generator <- function(beta_victim, beta_offender_1, 
                              beta_offender_2 = NULL, ratio_victim = 1,
                              ratio_offender = 1){
  beta_mixture <-  ((ratio_victim * beta_victim) + (ratio_offender * beta_offender_1)) / (ratio_victim + ratio_offender)
  return(beta_mixture)
}


mixture_deconvolution <- function(beta_mixture, beta_victim, ratio_victim = 1,
                                  ratio_offender = 1){
  beta_offender <- (((ratio_victim + ratio_offender) * beta_mixture) - (ratio_victim * beta_victim)) / ratio_offender
  return(beta_offender)
}


#Deconvolution of different DNA mixture ratios

ratios <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30)

#ratio=1
for (ratio in ratios) {
  #Create mixture
  if (ratio == 1) {
    table_pred_age_main_contr <- data.frame(ratio = character(),
                                                      Horvath = numeric(),
                                                      BLUP = numeric(),
                                                      EN = numeric(),
                                                      skinHorvath = numeric())
    table_pred_age_main_contr <- rbind(table_pred_age_main_contr, 
                                                c("Single_source", predicted_age[predicted_age$Sample_name == "AE_6",
                                                           c("Horvath", "BLUP", "EN", "skinHorvath")]))
    colnames(table_pred_age_main_contr)[1] <- "Ratio" 
  }
  mixture_28_55_ratio <- mixture_generator(beta_values_autosomal_cg_noNAs$AF_0,
                                           beta_values_autosomal_cg_noNAs$AE_0,
                                           ratio_victim = 1, ratio_offender = ratio)

  #Deconvolute mixture, extract beta value offender with known beta value victim
  beta_offender_55_ratio <- mixture_deconvolution(mixture_28_55_ratio, 
                                                  beta_values_autosomal_cg_noNAs$AF_0, 
                                                  ratio_victim = 1, ratio_offender = ratio)

  #sum(beta_offender_55_ratio > 1)
  #sum(beta_offender_55_ratio < 0)
  #Correct beta values of the offender
  beta_offender_55_ratio[beta_offender_55_ratio > 1] <- 1
  beta_offender_55_ratio[beta_offender_55_ratio < 0] <- 0
  #Predict age (we duplicate beta_offender_55_ratio in the dataframe because methylclock package need at least two individuals to infer the age)
  df_for_prediction <- data.frame(rownames(beta_values_autosomal_cg_noNAs),
                                  beta_offender_55_ratio,
                                  beta_offender_55_ratio)
  predicted_age_mixture <- DNAmAge(df_for_prediction, clocks = c("BLUP", "EN", "Horvath", "skinHorvath"))
  colnames(predicted_age_mixture)[1] <- "Ratio"
  predicted_age_mixture$Ratio <- c(paste0("1:", ratio))
  
  table_pred_age_main_contr <- rbind(table_pred_age_main_contr, predicted_age_mixture[1,])
}

write_xlsx(table_pred_age_main_contr, 
           'G:/FAELLES/Dokumenter/BRP/5_Projects/Mixture_deconvolution/Results/In-silico_mixture/Predicted_age_in_silico_deconvolute_mixture_main_contr_offender_22-03-2024.xlsx')


write_xlsx(table_pred_age_main_contr, 
           'G:/FAELLES/Dokumenter/BRP/5_Projects/Mixture_deconvolution/Results/In-silico_mixture/Predicted_age_in_silico_deconvolute_mixture_main_contr_victim_22-03-2024.xlsx')


#Create beta value vector with random error
#beta_with_err <- beta_values_autosomal_cg_noNAs$AF_6 + ifelse(runif(length(beta_values_autosomal_cg_noNAs$AF_6)) > 0.5, 0.1, -0.1)




#### 4. Test impact of technical variability (measure the same sample more times for cohort of 63 individuals from Logitudinal study) -------------------------------
predicted_age_metadata_2021
predicted_age_metadata_2021$AE <- round(abs(predicted_age_metadata_2021$Age - predicted_age_metadata_2021$BLUP), 2)

individuals <- colnames(betas_autosomal_cpgs_noNAs_2021)
shuffled_vector <- sample(individuals)

# Create random pairs
random_pairs <- matrix(individuals, ncol = 2, byrow = TRUE)

#Upload cpg of epigenetic clock to asses reconstruction of DNAm profile
checkClocks(betas_autosomal_cpgs_noNAs_2021)

# Display the pairs
random_pairs
predicted_age_metadata_2021 <- as.data.frame(predicted_age_metadata_2021)

#Initialize dataset:
median_betas_reconstructed_profile <- data.frame(
  Ratio = character(),
  Tecnical_error = numeric(),
  Horvath = numeric(),
  SkinHorvath = numeric(),
  BLUP = numeric(),
  EN = numeric(),
  stringsAsFactors = FALSE
)

f <- 123
for (n in 1:nrow(random_pairs)){
  victim_name <- random_pairs[n,][1] 
  offender_name <- random_pairs[n,][2] 
  
  for (i in c(1:3)){
    ratios <- c("1:1", "1:4", "1:10", "4:1", "10:1")
    #ratio <- "1:10"
    tecn_err <- c(0, 0.01 , 0.025, 0.05, 0.1) 
    #err <- 0
    set.seed(f) 
    for (ratio in ratios) {
      #Create mixture
      ratio_split <- strsplit(ratio, ":")
      ratio_victim <- as.numeric(ratio_split[[1]][1])
      ratio_offender <- as.numeric(ratio_split[[1]][2])
      
      for (err in tecn_err){
        print(paste(i, ratio, err))
        random_signs <- sample(c(-1, 1), nrow(betas_autosomal_cpgs_noNAs_2021), replace = TRUE)
        random_errors <- random_signs * err
        beta_value_victim <- betas_autosomal_cpgs_noNAs_2021[,victim_name] + random_errors
        beta_value_victim[beta_value_victim > 1] <- 1
        beta_value_victim[beta_value_victim < 0] <- 0
        
        if(ratio == "1:1"& err == 0){
          table_pred_age_main_contr <- data.frame(ratio = character(),
                                                  tecn_err = numeric(),
                                                  Horvath = numeric(),
                                                  BLUP = numeric(),
                                                  EN = numeric(),
                                                  skinHorvath = numeric())
          values <- as.numeric(predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == offender_name, 
                                                           c("Horvath", "BLUP", "EN", "skinHorvath")])
          
          # Append a new row to the table_pred_age_main_contr data frame
          table_pred_age_main_contr <- rbind(table_pred_age_main_contr, 
                                             data.frame(ratio = "Single_source", 
                                                        tecn_err = NA, 
                                                        Horvath = values[1], 
                                                        BLUP = values[2], 
                                                        EN = values[3], 
                                                        skinHorvath = values[4]))
          colnames(table_pred_age_main_contr)[c(1,2)] <- c("Ratio", "Tecnical_error") 
        }
        mixture_28_55_ratio <- mixture_generator(beta_victim=betas_autosomal_cpgs_noNAs_2021[,victim_name],
                                                 beta_offender_1=betas_autosomal_cpgs_noNAs_2021[,offender_name],
                                                 ratio_victim = ratio_victim, ratio_offender = ratio_offender)
        
        #Deconvolute mixture, extract beta value offender with known beta value victim
        beta_offender_55_ratio <- mixture_deconvolution(mixture_28_55_ratio, 
                                                        beta_value_victim, 
                                                        ratio_victim = ratio_victim, ratio_offender = ratio_offender)
        
        

        #sum(beta_offender_55_ratio > 1)
        #sum(beta_offender_55_ratio < 0)
        #Correct beta values of the offender
        beta_offender_55_ratio[beta_offender_55_ratio > 1] <- 1
        beta_offender_55_ratio[beta_offender_55_ratio < 0] <- 0
        
        #Assess performance of deconvolution and recreation of DNAm profile of offender  performed 
        absolute_error_mixture <- abs(beta_offender_55_ratio - betas_autosomal_cpgs_noNAs_2021[,offender_name])
        names(absolute_error_mixture) <- rownames(betas_autosomal_cpgs_noNAs_2021)
        median_delta_betas_vector <- data.frame(
          Ratio = ratio,
          Tecnical_error = err,
          Horvath = median(absolute_error_mixture[coefHorvath$CpGmarker[-1]], na.rm = TRUE),
          SkinHorvath = median(absolute_error_mixture[coefSkin$CpGmarker[-1]], na.rm = TRUE),
          BLUP = median(absolute_error_mixture[coefBLUP$CpGmarker[-1]], na.rm = TRUE),
          EN = median(absolute_error_mixture[coefEN$CpGmarker[-1]], na.rm = TRUE)
        )
        
        median_betas_reconstructed_profile <- rbind(median_betas_reconstructed_profile, median_delta_betas_vector)
        
        
        #Predict age (we duplicate beta_offender_55_ratio in the dataframe because methylclock package need at least two individuals to infer the age)
        df_for_prediction <- data.frame(rownames(betas_autosomal_cpgs_noNAs_2021),
                                        beta_offender_55_ratio,
                                        beta_offender_55_ratio)
        predicted_age_mixture <- DNAmAge(df_for_prediction, clocks = c("BLUP", "EN", "Horvath", "skinHorvath"))
        colnames(predicted_age_mixture)[1] <- "Ratio"
        predicted_age_mixture$Ratio <- ratio
        predicted_age_mixture$Tecnical_error <- err


        table_pred_age_main_contr <- rbind(table_pred_age_main_contr, predicted_age_mixture[1,])
      }}
    
    
    table_pred_age_main_contr[, c(3,4,5,6)] <- abs(table_pred_age_main_contr[, c(3,4,5,6)] - 
                                                     predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == offender_name,"Age"])
    
    if (n == 1 & i == 1 & ratio == "10:1" & err == 0.10) {
      df <- table_pred_age_main_contr
    } else {
      df <- rbind(df, table_pred_age_main_contr)
    }
  }
  }


saveRDS(df, paste0(results_path, "Age_prediction_tecnology_err_10_simulations_Longitudinal_data_R_16-01-2026.rds"))
saveRDS(median_betas_reconstructed_profile, paste0(results_path, "Median_betas_reconstructed_profile_R_16-01-2026.rds"))

df <- readRDS(paste0(results_path, "Age_prediction_tecnology_err_10_simulations_Longitudinal_data_11-09-2025.rds"))
median_betas_reconstructed_profile <- readRDS(paste0(results_path, "Median_betas_reconstructed_profile_11-09-2025.rds"))

## Plotting median beta value difference reconstructed profile in EPIC v2.0 precision analysis
median_betas_reconstructed_profile[,c(3:6)] <- round(median_betas_reconstructed_profile[,c(3:6)], 9)

#Average median beta reconstructed profile for the three iteration and 32 individuals 
median_betas_summary <- median_betas_reconstructed_profile %>%
  group_by(Ratio, Tecnical_error) %>%
  summarise(
    Horvath = mean(Horvath, na.rm = TRUE),
    BLUP= mean(BLUP, na.rm = TRUE),
    SkinHorvath = mean(SkinHorvath, na.rm = TRUE),
    EN = mean(EN, na.rm = TRUE),
    .groups = "drop"
  )

#Invert ratio
median_betas_summary$Ratio <- paste0(median_betas_summary$Ratio, "_old")

median_betas_summary$Ratio <- gsub("1:1_old", "1:1", median_betas_summary$Ratio)
median_betas_summary$Ratio <- gsub("1:10_old", "10:1", median_betas_summary$Ratio)
median_betas_summary$Ratio <- gsub("10:1_old", "1:10", median_betas_summary$Ratio)
median_betas_summary$Ratio <- gsub("1:4_old", "4:1", median_betas_summary$Ratio)
median_betas_summary$Ratio <- gsub("4:1_old", "1:4", median_betas_summary$Ratio)
median_betas_summary$Ratio <- factor(median_betas_summary$Ratio, levels = c("10:1", "4:1", "1:1", "1:4", "1:10"))
df_long$Clock <- paste(df_long$Clock, "clock")
median_betas_summary_long <- median_betas_summary %>%
  pivot_longer(cols = -c(Ratio, Tecnical_error),
               names_to = "Clock", 
               values_to = "Errors")

# Factor levels for ratio (legend order)
median_betas_summary_long$Ratio <- factor(
  median_betas_summary_long$Ratio,
  levels = c("10:1", "4:1", "2:1", "1:1", "1:2", "1:4", "1:10")
)

median_betas_summary_long$Clock <- paste(median_betas_summary_long$Clock, "clock")


# Plot
diff_beta_precision_plot <- ggplot(
  median_betas_summary_long,
  aes(
    x = Tecnical_error,
    y = Errors,
    color = Ratio,
    group = interaction(Ratio, Clock)
  )
) +
  geom_line(size = 1) +  # Connect points with lines
  geom_point(size = 2) +  # Optional: Add points at each time point
  facet_wrap(~ Clock) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 0.50)) +
  scale_x_continuous(breaks = c(0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.10),
                     #labels = c("0", "0.005", "0.01", "0.02", "0.03", "0.05", "0.1"),
                     expand = c(0.001, 0.001)) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = -45, vjust = 0.5, hjust = 0),
    legend.position = "top",
    panel.grid.major.x = element_line(color = "grey", size = 0.4),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "lightgrey"),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing.x = unit(1.5, "lines"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    plot.margin = unit(c(1, 1.5, 1, 1), "lines")
  ) + theme(
    plot.margin = unit(c(1, 1, 1, 1.5), "cm")  # top, right, bottom, left
  ) +
  labs(title = "",
       x = "DNAm Technical Noise (|Δβ| Between Replicates)",
       y = "Median |Δβ|",
       color = "Suspect-to-Victim Ratio")

diff_beta_precision_plot

ggsave(paste0(results_path,"/Median_delta_betas_reconstructed_DNAm_profiles_precision_analysis_02-10-2025.png"), 
       diff_beta_precision_plot, width = 11, height = 7, dpi = 600, bg = "white")



## Plotting age prediction accuracy in EPIC v2.0 precision analysis
df[, c(3,4,5,6)] <- df[, c(3,4,5,6)]/(32*3)

df_long <- df %>%
  pivot_longer(cols = -c(Ratio, Tecnical_error),
               names_to = "Clock", 
               values_to = "MAE")


df_long <- df_long[-c(1:4),]
#Invert ratio
df_long$Ratio <- paste0(df_long$Ratio, "_old")

df_long$Ratio <- gsub("1:1_old", "1:1", df_long$Ratio)
df_long$Ratio <- gsub("1:10_old", "10:1", df_long$Ratio)
df_long$Ratio <- gsub("10:1_old", "1:10", df_long$Ratio)
df_long$Ratio <- gsub("1:4_old", "4:1", df_long$Ratio)
df_long$Ratio <- gsub("4:1_old", "1:4", df_long$Ratio)
df_long$Ratio <- factor(df_long$Ratio, levels = c("10:1", "4:1", "1:1", "1:4", "1:10"))
df_long$Clock <- paste(df_long$Clock, "clock")
# df_long <- df_long[df_long$Tecnical_error != 0.1,]

plot_tecn_err <- ggplot(df_long, aes(x = Tecnical_error, y = MAE, color = Ratio, group = interaction(Ratio, Clock))) +
  geom_line(size = 1) +  # Connect points with lines
  geom_point(size = 2) +  # Optional: Add points at each time point
  facet_wrap(~ Clock) +  # Facet by the "Group" variable
  theme_minimal() +  # Use a minimal theme
  labs(title = "",
       x = "DNAm Technical Noise (|Δβ| Between Replicates)",
       y = "MAE (years)",
       color = "Suspect-to-Victim Ratio") +
  coord_cartesian(ylim = c(0, 25)) +
  scale_x_continuous(breaks = c(0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.10),
                     #labels = c("0", "0.005", "0.01", "0.02", "0.03", "0.05", "0.1"),
                     expand = c(0.001, 0.001)) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Customize facet labels
    axis.text = element_text(size = 12),  # Customize axis text
    axis.text.x = element_text(angle = -45,vjust = 0, hjust = 0.1, size = 12),  # Customize axis text
    legend.position = "top",  # Position the legend at the top
    panel.grid.major.x = element_line(color = "grey", size = 0.4),  # Set the x-axis major grid lines
    panel.grid.minor.x = element_blank(),  # Remove minor grid lines on the x-axis
    panel.grid.major.y = element_line(color = "lightgrey"),  # Keep y-axis grid lines
    panel.grid.minor.y = element_blank(),  # Remove y-axis minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Set border
    panel.spacing.x = unit(1.5, "lines"),
    axis.title.x = element_text(size = 14),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 14),   # Increase legend text size
    legend.title = element_text(size = 15),
    plot.margin = unit(c(1, 1.5, 1, 1), "lines")
  ) + theme(
    plot.margin = unit(c(1, 1, 1, 1.5), "cm")  # top, right, bottom, left
  )

plot_tecn_err

ggsave(paste0(results_path,"Plot_tecnical_error_simulation_cohort_3_simulation_02-10-2025.png"), 
       plot_tecn_err, width = 11, height = 7, dpi = 600, bg = "white", limitsize = FALSE)

#### 4. Revision: Test impact of technical variability (bar plots)------------------------------
df <- readRDS(paste0(results_path, "Age_prediction_tecnology_err_10_simulations_Longitudinal_data_R_16-01-2026.rds"))
median_betas_reconstructed_profile <- readRDS(paste0(results_path, "Median_betas_reconstructed_profile_R_16-01-2026.rds"))

## Plotting median beta value difference reconstructed profile in EPIC v2.0 precision analysis
#Invert ratio
median_betas_reconstructed_profile$Ratio <- paste0(median_betas_reconstructed_profile$Ratio, "_old")

median_betas_reconstructed_profile$Ratio <- gsub("1:1_old", "1:1", median_betas_reconstructed_profile$Ratio)
median_betas_reconstructed_profile$Ratio <- gsub("1:10_old", "10:1", median_betas_reconstructed_profile$Ratio)
median_betas_reconstructed_profile$Ratio <- gsub("10:1_old", "1:10", median_betas_reconstructed_profile$Ratio)
median_betas_reconstructed_profile$Ratio <- gsub("1:4_old", "4:1", median_betas_reconstructed_profile$Ratio)
median_betas_reconstructed_profile$Ratio <- gsub("4:1_old", "1:4", median_betas_reconstructed_profile$Ratio)
median_betas_reconstructed_profile$Ratio <- factor(median_betas_reconstructed_profile$Ratio, levels = c("10:1", "4:1", "1:1", "1:4", "1:10"))
median_betas_reconstructed_profile_long <- median_betas_reconstructed_profile %>% 
  pivot_longer(cols = -c(Ratio, Tecnical_error),
               names_to = "Clock", 
               values_to = "Errors")

# Factor levels for ratio (legend order)
median_betas_reconstructed_profile_long$Ratio <- factor(
  median_betas_reconstructed_profile_long$Ratio,
  levels = c("10:1", "4:1", "2:1", "1:1", "1:2", "1:4", "1:10")
)

median_betas_reconstructed_profile_long$Clock <- paste(median_betas_reconstructed_profile_long$Clock, "clock")
median_betas_reconstructed_profile_long$Tecnical_error_f <- factor(
  median_betas_reconstructed_profile_long$Tecnical_error,
  levels = c(0, 0.01, 0.025, 0.05, 0.10)
)

str(median_betas_reconstructed_profile_long)

# Plot
diff_beta_precision_plot <- ggplot(
  median_betas_reconstructed_profile_long,
  aes(
    x = Tecnical_error_f,
    y = Errors,
    color = Ratio,
    group = interaction(Tecnical_error_f, Ratio)
  )
) +
  # Error bars (mean ± SD)
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    position = position_dodge(width = 0.6),
    width = 0.25,
    linewidth = 0.5
  ) +
  # Mean point
  stat_summary(
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.6),
    size = 2
  ) +
  facet_wrap(
    ~ Clock, ncol = 1,
    labeller = labeller(
      Clock = c(
        "Horvath clock" = "Horvath clock (353 CpGs)",
        "SkinHorvath clock"  = "skinHorvath clock (391 CpGs)",
        "EN clock"    = "EN clock (514 CpGs)",
        "BLUP clock"     = "BLUP clock (319,607 CpGs)"
      )
    )) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(
    x = "DNAm Technical Noise (|Δβ| Between Replicates)",
    y = "Median |Δβ|",
    color = "Suspect-to-Victim Ratio"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing.y = unit(1.2, "lines"),
    plot.margin = unit(c(1, 1.5, 1, 1.5), "cm")
  )


ggsave(paste0(results_path,"/Median_delta_betas_reconstructed_DNAm_profiles_precision_analysis_R_16-01-2025.png"), 
       diff_beta_precision_plot, width = 10, height = 8.5, dpi = 600, bg = "white")



## Plotting age prediction accuracy in EPIC v2.0 precision analysis
df <- df[df$Ratio != "Single_source",]

#Invert ratio
df$Ratio <- paste0(df$Ratio, "_old")

df$Ratio <- gsub("1:1_old", "1:1", df$Ratio)
df$Ratio <- gsub("1:10_old", "10:1", df$Ratio)
df$Ratio <- gsub("10:1_old", "1:10", df$Ratio)
df$Ratio <- gsub("1:4_old", "4:1", df$Ratio)
df$Ratio <- gsub("4:1_old", "1:4", df$Ratio)
df$Ratio <- factor(df$Ratio, levels = c("10:1", "4:1", "1:1", "1:4", "1:10"))

df_long <- df %>%
  pivot_longer(cols = -c(Ratio, Tecnical_error),
               names_to = "Clock", 
               values_to = "Errors")

df_long$Clock <- paste(df_long$Clock, "clock")

# Factor levels for ratio (legend order)
df_long$Ratio <- factor(
  df_long$Ratio,
  levels = c("10:1", "4:1", "2:1", "1:1", "1:2", "1:4", "1:10")
)

str(df_long)


df_long$Tecnical_error_f <- factor(
  df_long$Tecnical_error,
  levels = c(0, 0.01, 0.025, 0.05, 0.10)
)

plot_tecn_err <- ggplot(
  df_long,
  aes(
    x = Tecnical_error_f,
    y = Errors,              # or MAE
    fill = Ratio,
    color = Ratio,
    group = interaction(Tecnical_error_f, Ratio)
  )
) +
  geom_boxplot(
    position = position_dodge2(width = 0.85, preserve = "single"),
    width = 0.7,
    color = "black",         # black box outlines
    linewidth = 0.3,
    outlier.alpha = 0.4,
    outlier.size = 0.8
  ) +
  facet_wrap(
    ~ Clock, ncol = 1,
    labeller = labeller(
      Clock = c(
        "Horvath clock" = "Horvath clock (353 CpGs)",
        "skinHorvath clock"  = "skinHorvath clock (391 CpGs)",
        "EN clock"    = "EN clock (514 CpGs)",
        "BLUP clock"     = "BLUP clock (319,607 CpGs)"
      )
    )) +
  coord_cartesian(ylim = c(0, 25)) +
  labs(
    x = "DNAm Technical Noise (|Δβ| Between Replicates)",
    y = "MAE (years)",
    fill = "Suspect-to-Victim Ratio",
    color = "Suspect-to-Victim Ratio"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80"),
    panel.spacing.y = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = unit(c(1, 1.5, 1, 1.5), "cm")
  )


ggsave(paste0(results_path,"Plot_tecnical_error_simulation_cohort_R_16-01-2025.png"), 
       plot_tecn_err, width = 10, height = 8.5, dpi = 600, bg = "white", limitsize = FALSE)




#### 5. Test impact of STR ratio calculation (cohort 63 individuals, Publication) -------------------------------

predicted_age_metadata_2021
predicted_age_metadata_2021$AE <- round(abs(predicted_age_metadata_2021$Age - predicted_age_metadata_2021$BLUP), 2)

individuals <- colnames(betas_autosomal_cpgs_noNAs_2021)
shuffled_vector <- sample(individuals)

# Create random pairs
random_pairs <- matrix(individuals, ncol = 2, byrow = TRUE)

# Display the pairs
random_pairs
predicted_age_metadata_2021 <- as.data.frame(predicted_age_metadata_2021)

#pair <- "AF_0, AD_0"
#ratio <- "10:1"
#err <- 2

#Initialize dataset:
median_betas_reconstructed_profile <- data.frame(
  Ratio = character(),
  Tecnical_error = numeric(),
  Horvath = numeric(),
  SkinHorvath = numeric(),
  BLUP = numeric(),
  EN = numeric(),
  stringsAsFactors = FALSE
)

#checkClocks(betas_autosomal_cpgs_noNAs_2021) #Upload clocks

for (n in 1:nrow(random_pairs)){
  victim_name <- random_pairs[n,][1]
  offender_name <- random_pairs[n,][2]
  #ratios <- c("1:1", "1:4", "1:10", "4:1", "10:1")
  ratios <- c("50:50", "20:80", "9.09090909:90.90909090", "80:20", "90.90909090:9.09090909")
  STR_err <- c(-10, -5, 0, 5, 10)
  
  #random_signs <- sample(c(-1, 1), length(beta_values_autosomal_cg_noNAs$AF_0), replace = TRUE)
  #random_errors <- random_signs * 0.01
  beta_value_victim <- betas_autosomal_cpgs_noNAs_2021[,victim_name] 
  
  for (err in STR_err) {
    #Create mixture
    
    for (ratio in ratios){
      print(paste(n, ratio, err))
      ratio_split <- strsplit(ratio, ":")
      ratio_victim_raw <- as.numeric(ratio_split[[1]][1])
      ratio_offender_raw <- as.numeric(ratio_split[[1]][2])
      
      if (ratio_victim_raw/ratio_offender_raw > 1){
        ratio_victim <- ratio_victim_raw/ratio_offender_raw
        ratio_offender <- 1
      } else {
        ratio_offender <- ratio_offender_raw/ratio_victim_raw
        ratio_victim <- 1
        
      }
      
      
      if (ratio_victim_raw/ratio_offender_raw > 1){
        ratio_victim_err <- (ratio_victim_raw + err)/(ratio_offender_raw - err)
        ratio_offender_err <- 1
        
      } else { if (ratio == "50:50") {
        if (err<0){
          ratio_offender_err <- (ratio_offender_raw + abs(err))/(ratio_victim_raw + err)
          ratio_victim_err <- 1
          
        } else {
          ratio_victim_err <- (ratio_victim_raw + err)/(ratio_offender_raw - err)
          ratio_offender_err <- 1
        }} else {
          ratio_offender_err <- (ratio_offender_raw + err)/(ratio_victim_raw - err)
          ratio_victim_err <- 1
        }}
      
      #print(paste(ratio_victim, ratio_offender))
      
      
      if (ratio == "50:50" & err == -10) { #or -5 substitute base on the analysis
        table_pred_age_main_contr <- data.frame(ratio = character(),
                                                tecn_err = numeric(),
                                                Horvath = numeric(),
                                                BLUP = numeric(),
                                                EN = numeric(),
                                                skinHorvath = numeric())
        table_pred_age_main_contr <- rbind(table_pred_age_main_contr, 
                                           c("Single_source", NA, predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == offender_name,
                                                                                c("Horvath", "BLUP", "EN", "skinHorvath")]))
        colnames(table_pred_age_main_contr)[c(1,2)] <- c("Ratio", "Tecnical_error") 
      }
      mixture_28_55_ratio <- mixture_generator(beta_victim=betas_autosomal_cpgs_noNAs_2021[,victim_name],
                                               beta_offender_1=betas_autosomal_cpgs_noNAs_2021[,offender_name],
                                               ratio_victim = ratio_victim, ratio_offender = ratio_offender)
      
      #Deconvolute mixture, extract beta value offender with known beta value victim
      beta_offender_55_ratio <- mixture_deconvolution(mixture_28_55_ratio, 
                                                      beta_value_victim, 
                                                      ratio_victim = ratio_victim_err, ratio_offender = ratio_offender_err)
      
      #sum(beta_offender_55_ratio > 1)
      #sum(beta_offender_55_ratio < 0)
      #Correct beta values of the offender
      beta_offender_55_ratio[beta_offender_55_ratio > 1] <- 1
      beta_offender_55_ratio[beta_offender_55_ratio < 0] <- 0
      
      #Assess performance of deconvolution and recreation of DNAm profile of offender  performed 
      absolute_error_mixture <- abs(beta_offender_55_ratio - betas_autosomal_cpgs_noNAs_2021[,offender_name])
      names(absolute_error_mixture) <- rownames(betas_autosomal_cpgs_noNAs_2021)
      median_delta_betas_vector <- data.frame(
        Ratio = ratio,
        Tecnical_error = err,
        Horvath = median(absolute_error_mixture[coefHorvath$CpGmarker[-1]], na.rm = TRUE),
        SkinHorvath = median(absolute_error_mixture[coefSkin$CpGmarker[-1]], na.rm = TRUE),
        BLUP = median(absolute_error_mixture[coefBLUP$CpGmarker[-1]], na.rm = TRUE),
        EN = median(absolute_error_mixture[coefEN$CpGmarker[-1]], na.rm = TRUE)
      )
      
      median_betas_reconstructed_profile <- rbind(median_betas_reconstructed_profile, median_delta_betas_vector)
      
      #Predict age (we duplicate beta_offender_55_ratio in the dataframe because methylclock package need at least two individuals to infer the age)
      df_for_prediction <- data.frame(rownames(betas_autosomal_cpgs_noNAs_2021),
                                      beta_offender_55_ratio,
                                      beta_offender_55_ratio)
      predicted_age_mixture <- DNAmAge(df_for_prediction, clocks = c("BLUP", "EN", "Horvath", "skinHorvath"))
      colnames(predicted_age_mixture)[1] <- "Ratio"
      predicted_age_mixture$Ratio <- ratio
      predicted_age_mixture$Tecnical_error <- err
      
      table_pred_age_main_contr <- rbind(table_pred_age_main_contr, predicted_age_mixture[1,])
    }}
  table_pred_age_main_contr[, c(3,4,5,6)] <- abs(table_pred_age_main_contr[, c(3,4,5,6)] - 
                                                   predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == offender_name,"Age"])
  
  if (n == 1) {
    df <- table_pred_age_main_contr
  } else {
    df <- rbind(df, table_pred_age_main_contr)
  }
  
}


#Divide by the total number of DNA mixture pair and save
saveRDS(df, paste0(results_path, "Age_prediction_STR_ratio_err_simulations_Longitudinal_data_R_10_17-01-2026.rds"))
saveRDS(median_betas_reconstructed_profile, paste0(results_path, "Median_betas_reconstructed_profile_STR_analysis_R_10_17-01-2026.rds"))


## Plotting median beta value difference reconstructed profile in EPIC v2.0 precision analysis
#Average median beta reconstructed profile for the three iteration and 32 individuals 
median_betas_summary <- median_betas_reconstructed_profile %>%
  group_by(Ratio, Tecnical_error) %>%
  summarise(
    Horvath = mean(Horvath, na.rm = TRUE),
    BLUP= mean(BLUP, na.rm = TRUE),
    SkinHorvath = mean(SkinHorvath, na.rm = TRUE),
    EN = mean(EN, na.rm = TRUE),
    .groups = "drop"
  )

median_betas_reconstructed_profile[,c(3:6)] <- round(median_betas_reconstructed_profile[,c(3:6)], 9)

median_betas_summary_long <- median_betas_summary %>%
  pivot_longer(cols = -c(Ratio, Tecnical_error),
               names_to = "Clock", 
               values_to = "Errors")

median_betas_summary_long$Ratio[median_betas_summary_long$Ratio == "50:50"] <- "1:1"
median_betas_summary_long$Ratio[median_betas_summary_long$Ratio == "20:80"] <- "1:4"
median_betas_summary_long$Ratio[median_betas_summary_long$Ratio == "80:20"] <- "4:1"
median_betas_summary_long$Ratio[median_betas_summary_long$Ratio == "90.90909090:9.09090909"] <- "10:1"
median_betas_summary_long$Ratio[median_betas_summary_long$Ratio == "9.09090909:90.90909090"] <- "1:10"

median_betas_summary_long$Ratio <- factor(median_betas_summary_long$Ratio, levels = c("10:1", "4:1", "1:1", "1:4", "1:10"))

#Invert ratio
median_betas_summary_long$Ratio <- paste0(median_betas_summary_long$Ratio, "_old")

median_betas_summary_long$Ratio <- gsub("1:1_old", "1:1", median_betas_summary_long$Ratio)
median_betas_summary_long$Ratio <- gsub("1:10_old", "10:1", median_betas_summary_long$Ratio)
median_betas_summary_long$Ratio <- gsub("10:1_old", "1:10", median_betas_summary_long$Ratio)
median_betas_summary_long$Ratio <- gsub("1:4_old", "4:1", median_betas_summary_long$Ratio)
median_betas_summary_long$Ratio <- gsub("4:1_old", "1:4", median_betas_summary_long$Ratio)
median_betas_summary_long$Ratio <- factor(median_betas_summary_long$Ratio, levels = c("10:1", "4:1", "1:1", "1:4", "1:10"))
median_betas_summary_long$Clock <- paste(median_betas_summary_long$Clock, "clock")

beta_err_STR_plot_analysis <- ggplot(median_betas_summary_long, aes(x = Tecnical_error, y = Errors, color = Ratio, group = interaction(Ratio, Clock))) +
  geom_line(size = 1) +  # Connect points with lines
  geom_point(size = 2) +  # Optional: Add points at each time point
  facet_wrap(~ Clock) +  # Facet by the "Group" variable
  theme_minimal() +  # Use a minimal theme
  labs(title = "",
       x = "Error in Victim DNA Proportion (%)",
       y = "Median |Δβ|",
       color = "Suspect-to-Victim Ratio") +
  coord_cartesian(ylim = c(0, 0.1)) +
  scale_x_continuous(breaks = c(-5, -3, -2, -1, 0, 1, 2, 3, 5),
                     #                    #labels = c("0", "0.005", "0.01", "0.02", "0.03", "0.05", "0.1"),
                     expand = c(0.005, 0.005)) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Customize facet labels
    axis.text = element_text(size = 12),  # Customize axis text
    axis.text.x = element_text(angle = 0,vjust = 0, hjust = 0.1, size = 12),  # Customize axis text
    legend.position = "top",  # Position the legend at the top
    panel.grid.major.x = element_line(color = "grey", size = 0.4),  # Set the x-axis major grid lines
    panel.grid.minor.x = element_blank(),  # Remove minor grid lines on the x-axis
    panel.grid.major.y = element_line(color = "lightgrey"),  # Keep y-axis grid lines
    panel.grid.minor.y = element_blank(),  # Remove y-axis minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Set border
    panel.spacing.x = unit(1.5, "lines"),
    axis.title.x = element_text(size = 14, margin = margin(t=10)),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 14),   # Increase legend text size
    legend.title = element_text(size = 15),
    plot.margin = unit(c(1, 1.5, 1, 1), "lines")
  )


ggsave(paste0(results_path,"Median_delta_betas_reconstructed_DNAm_profiles_STR_analysis_15-09-2025.png"), 
       beta_err_STR_plot_analysis, width = 11, height = 7, dpi = 600, bg = "white")



## Plot age prediction accuracy for STR ratio estimation error
df_long <- df %>%
  pivot_longer(cols = -c(Ratio, Tecnical_error),
               names_to = "Clock", 
               values_to = "MAE")

df_long <- df_long[-c(1:4),]

df_long$Ratio[df_long$Ratio == "50:50"] <- "1:1"
df_long$Ratio[df_long$Ratio == "20:80"] <- "1:4"
df_long$Ratio[df_long$Ratio == "80:20"] <- "4:1"
df_long$Ratio[df_long$Ratio == "90.90909090:9.09090909"] <- "10:1"
df_long$Ratio[df_long$Ratio == "9.09090909:90.90909090"] <- "1:10"

df_long$Ratio <- factor(df_long$Ratio, levels = c("10:1", "4:1", "1:1", "1:4", "1:10"))

#Invert ratio
df_long$Ratio <- paste0(df_long$Ratio, "_old")

df_long$Ratio <- gsub("1:1_old", "1:1", df_long$Ratio)
df_long$Ratio <- gsub("1:10_old", "10:1", df_long$Ratio)
df_long$Ratio <- gsub("10:1_old", "1:10", df_long$Ratio)
df_long$Ratio <- gsub("1:4_old", "4:1", df_long$Ratio)
df_long$Ratio <- gsub("4:1_old", "1:4", df_long$Ratio)
df_long$Ratio <- factor(df_long$Ratio, levels = c("10:1", "4:1", "1:1", "1:4", "1:10"))
df_long$Clock <- paste(df_long$Clock, "clock")

STR_plot <- ggplot(df_long, aes(x = Tecnical_error, y = MAE, color = Ratio, group = interaction(Ratio, Clock))) +
  geom_line(size = 1) +  # Connect points with lines
  geom_point(size = 2) +  # Optional: Add points at each time point
  facet_wrap(~ Clock) +  # Facet by the "Group" variable
  theme_minimal() +  # Use a minimal theme
  labs(title = "",
       x = "Error in Victim DNA Proportion (%)",
       y = "MAE (years)",
       color = "Suspect-to-Victim Ratio") +
  coord_cartesian(ylim = c(0, 16)) +
  scale_x_continuous(breaks = c(-5, -3, -2, -1, 0, 1, 2, 3, 5),
                     #                    #labels = c("0", "0.005", "0.01", "0.02", "0.03", "0.05", "0.1"),
                     expand = c(0.005, 0.005)) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),  # Customize facet labels
    axis.text = element_text(size = 12),  # Customize axis text
    axis.text.x = element_text(angle = 0,vjust = 0, hjust = 0.1, size = 12),  # Customize axis text
    legend.position = "top",  # Position the legend at the top
    panel.grid.major.x = element_line(color = "grey", size = 0.4),  # Set the x-axis major grid lines
    panel.grid.minor.x = element_blank(),  # Remove minor grid lines on the x-axis
    panel.grid.major.y = element_line(color = "lightgrey"),  # Keep y-axis grid lines
    panel.grid.minor.y = element_blank(),  # Remove y-axis minor grid lines
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Set border
    panel.spacing.x = unit(1.5, "lines"),
    axis.title.x = element_text(size = 14, margin = margin(t=10)),  # Increase x-axis title size
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 14),   # Increase legend text size
    legend.title = element_text(size = 15),
    plot.margin = unit(c(1, 1.5, 1, 1), "lines")
  )


ggsave(paste0(results_path,"Plot_STR_ratio_err_cohort_15-09-2025.png"), 
       STR_plot, width = 11, height = 7, dpi = 600, bg = "white")


#### 5. Revision: Test impact of technical variability (bar plots)------------------------------
df <- readRDS(paste0(results_path, "Age_prediction_STR_ratio_err_simulations_Longitudinal_data_R_17-01-2026.rds"))
median_betas_reconstructed_profile <- readRDS(paste0(results_path, "Median_betas_reconstructed_profile_STR_analysis_R_17-01-2026.rds"))

## Plotting median beta value difference reconstructed profile in EPIC v2.0 precision analysis

median_betas_reconstructed_profile$Ratio[median_betas_reconstructed_profile$Ratio == "50:50"] <- "1:1"
median_betas_reconstructed_profile$Ratio[median_betas_reconstructed_profile$Ratio == "20:80"] <- "1:4"
median_betas_reconstructed_profile$Ratio[median_betas_reconstructed_profile$Ratio == "80:20"] <- "4:1"
median_betas_reconstructed_profile$Ratio[median_betas_reconstructed_profile$Ratio == "90.90909090:9.09090909"] <- "10:1"
median_betas_reconstructed_profile$Ratio[median_betas_reconstructed_profile$Ratio == "9.09090909:90.90909090"] <- "1:10"

median_betas_reconstructed_profile$Ratio <- factor(median_betas_reconstructed_profile$Ratio, levels = c("10:1", "4:1", "1:1", "1:4", "1:10"))

#Invert ratio
median_betas_reconstructed_profile$Ratio <- paste0(median_betas_reconstructed_profile$Ratio, "_old")

median_betas_reconstructed_profile$Ratio <- gsub("1:1_old", "1:1", median_betas_reconstructed_profile$Ratio)
median_betas_reconstructed_profile$Ratio <- gsub("1:10_old", "10:1", median_betas_reconstructed_profile$Ratio)
median_betas_reconstructed_profile$Ratio <- gsub("10:1_old", "1:10", median_betas_reconstructed_profile$Ratio)
median_betas_reconstructed_profile$Ratio <- gsub("1:4_old", "4:1", median_betas_reconstructed_profile$Ratio)
median_betas_reconstructed_profile$Ratio <- gsub("4:1_old", "1:4", median_betas_reconstructed_profile$Ratio)

median_betas_reconstructed_profile_long <- median_betas_reconstructed_profile %>% 
  pivot_longer(cols = -c(Ratio, Tecnical_error),
               names_to = "Clock", 
               values_to = "Errors")

median_betas_reconstructed_profile_long$Ratio <- factor(median_betas_reconstructed_profile_long$Ratio, levels = c("10:1", "4:1", "1:1", "1:4", "1:10"))
median_betas_reconstructed_profile_long$Clock <- paste(median_betas_reconstructed_profile_long$Clock, "clock")




median_betas_reconstructed_profile_long$Tecnical_error_f <- factor(
  median_betas_reconstructed_profile_long$Tecnical_error,
  levels = c(-5, -2, 0, 2, 5) #c(-10, -5, 0, 5, 10)
)

str(median_betas_reconstructed_profile_long)

# Plot
beta_err_STR_plot_analysis <- ggplot(
  median_betas_reconstructed_profile_long,
  aes(
    x = Tecnical_error_f,
    y = Errors,
    color = Ratio,
    group = interaction(Tecnical_error_f, Ratio)
  )
) +
  # Error bars: mean ± absolute error
  stat_summary(
    fun.data = function(x) {
      m <- mean(x, na.rm = TRUE)
      ae <- mean(abs(x - m), na.rm = TRUE)
      data.frame(y = m, ymin = m - ae, ymax = m + ae)
    },
    geom = "errorbar",
    position = position_dodge(width = 0.6),
    width = 0.15,
    linewidth = 0.4
  ) +
  # Mean point
  stat_summary(
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.6),
    size = 2.5
  ) +
  facet_wrap(
    ~ Clock, ncol = 1,
    labeller = labeller(
      Clock = c(
        "Horvath clock"      = "Horvath clock (353 CpGs)",
        "skinHorvath clock"  = "SkinHorvath clock (391 CpGs)",
        "EN clock"           = "EN clock (514 CpGs)",
        "BLUP clock"         = "BLUP clock (319,607 CpGs)"
      )
    )) +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(
    x = "Error in Victim DNA Proportion (%)",
    y = "Median |Δβ|",
    color = "Suspect-to-Victim Ratio"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14),
    legend.position = "top",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80"),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing.y = unit(1.5, "lines"),
    plot.margin = unit(c(1, 1.5, 1, 1.5), "cm")
  )

beta_err_STR_plot_analysis

ggsave(paste0(results_path,"/Median_delta_betas_reconstructed_DNAm_profiles_STR_analysis_R_10_17-01-2025.png"), 
       beta_err_STR_plot_analysis, width = 10, height = 8.5, dpi = 600, bg = "white")


## Plotting age prediction accuracy in EPIC v2.0 precision analysis
df <- df[df$Ratio != "Single_source",]

df$Ratio[df$Ratio == "50:50"] <- "1:1"
df$Ratio[df$Ratio == "20:80"] <- "1:4"
df$Ratio[df$Ratio == "80:20"] <- "4:1"
df$Ratio[df$Ratio == "90.90909090:9.09090909"] <- "10:1"
df$Ratio[df$Ratio == "9.09090909:90.90909090"] <- "1:10"

df$Ratio <- factor(df$Ratio, levels = c("10:1", "4:1", "1:1", "1:4", "1:10"))

#Invert ratio
df$Ratio <- paste0(df$Ratio, "_old")

df$Ratio <- gsub("1:1_old", "1:1", df$Ratio)
df$Ratio <- gsub("1:10_old", "10:1", df$Ratio)
df$Ratio <- gsub("10:1_old", "1:10", df$Ratio)
df$Ratio <- gsub("1:4_old", "4:1", df$Ratio)
df$Ratio <- gsub("4:1_old", "1:4", df$Ratio)

df_long <- df %>%
  pivot_longer(cols = -c(Ratio, Tecnical_error),
               names_to = "Clock", 
               values_to = "MAE")

df_long$Ratio <- factor(df_long$Ratio, levels = c("10:1", "4:1", "1:1", "1:4", "1:10"))
df_long$Clock <- paste(df_long$Clock, "clock")

df_long$Tecnical_error_f <- factor(
  df_long$Tecnical_error,
  levels = c(-10, -5, 0, 5, 10) #c(-5, -2, 0, 2, 5)
)

#Plot
plot_tecn_err_STR <- ggplot(
  df_long,
  aes(
    x = Tecnical_error_f,
    y = MAE,                     # <-- your column is MAE
    fill = Ratio,
    color = Ratio,
    group = interaction(Tecnical_error_f, Ratio)
  )
) +
  geom_boxplot(
    position = position_dodge2(width = 0.85, preserve = "single"),
    width = 0.7,
    color = "black",             # black box outlines
    linewidth = 0.3,
    outlier.alpha = 0.4,
    outlier.size = 0.8
  ) +
  facet_wrap(
    ~ Clock, ncol = 1,
    labeller = labeller(
      Clock = c(
        "Horvath clock"      = "Horvath clock (353 CpGs)",
        "skinHorvath clock"  = "skinHorvath clock (391 CpGs)",
        "EN clock"           = "EN clock (514 CpGs)",
        "BLUP clock"         = "BLUP clock (319,607 CpGs)"
      )
    )
  ) +
  coord_cartesian(ylim = c(0, 25)) +
  labs(
    x = "Error in Victim DNA Proportion (%)",   # <-- change if you want the old x-label instead
    y = "MAE (years)",
    fill = "Suspect-to-Victim Ratio",
    color = "Suspect-to-Victim Ratio"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80"),
    panel.spacing.y = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = unit(c(1, 1.5, 1, 1.5), "cm")
  )


ggsave(paste0(results_path,"Plot_STR_ratio_err_cohort_R_10_17-01-2025.png"), 
       plot_tecn_err_STR, width = 10, height = 8.5, dpi = 600, bg = "white", limitsize = FALSE)




#### 6. Test impact of Age differences in absolute error of offender from DNA mixture 2021 Longitudinal Cohort (Boxplot, Publication) -------------------

ages_2021 <- predicted_age_metadata_2021[,c(1,6)]

# Create all possible pairs of individuals
pairs <- expand.grid(predicted_age_metadata_2021$Sample_ID, predicted_age_metadata_2021$Sample_ID)
colnames(pairs) <- c("Name1", "Name2")

# Remove self-pairings
pairs <- pairs %>%
  filter(Name1 != Name2)
# Add ages for both individuals in the pair
pairs <- pairs %>%
  left_join(ages_2021, by = c("Name1" = "Sample_ID")) %>%
  rename(Age1 = Age) %>%
  left_join(ages_2021, by = c("Name2" = "Sample_ID")) %>%
  rename(Age2 = Age)
# Calculate the age difference
pairs <- pairs %>%
  mutate(AgeDifference = abs(Age1 - Age2))

# Filter pairs by age difference criteria
pairs_less_5 <- pairs %>%
  filter(AgeDifference < 5)

pairs_more_20 <- pairs %>%
  filter(AgeDifference > 20)


# Function to select unique pairs
select_unique_pairs <- function(pairs, n) {
  selected <- data.frame()
  used_individuals <- c()
  
  for (i in 1:nrow(pairs)) {
    pair <- pairs[i, ]
    if (!(pair$Name1 %in% used_individuals || pair$Name2 %in% used_individuals)) {
      selected <- bind_rows(selected, pair)
      used_individuals <- c(used_individuals, pair$Name1, pair$Name2)
    }
    if (nrow(selected) >= n) break
  }
  
  return(selected)
}


# Select unique pairs for each group
n <- min(nrow(pairs_less_5), nrow(pairs_more_20), 31 / 2) # Balance and limit to 31 pairs

selected_less_5 <- select_unique_pairs(pairs_less_5, n)
selected_more_20 <- select_unique_pairs(pairs_more_20, n)

# Combine selected pairs
random_pairs <- bind_rows(selected_less_5, selected_more_20)

for (n in 1:nrow(random_pairs)){
  ratios <- c("1:1", "10:1")
  victim_name <- random_pairs[n,][[1]] 
  offender_name <- random_pairs[n,][[2]] 
  
  random_signs <- sample(c(-1, 1), nrow(betas_autosomal_cpgs_noNAs_2021), replace = TRUE)
  random_errors <- random_signs * 0.01
  beta_value_victim <- betas_autosomal_cpgs_noNAs_2021[,victim_name] + random_errors
  #Create mixture
  for (ratio in ratios){
    print(paste(victim_name, offender_name, ratio))
    ratio_split <- strsplit(ratio, ":")
    ratio_victim <- as.numeric(ratio_split[[1]][1])
    ratio_offender <- as.numeric(ratio_split[[1]][2])
    
    if(victim_name == random_pairs[1,][[1]] & offender_name == random_pairs[1,][[2]] & ratio == "1:1"){
      table_pred_age_main_contr <- data.frame(Pair = character(),
                                              Ratio = character(),
                                              Age_gap = numeric(),
                                              Horvath = numeric(),
                                              BLUP = numeric(),
                                              EN = numeric(),
                                              skinHorvath = numeric())
    }
    mixture_28_55_ratio <- mixture_generator(beta_victim=betas_autosomal_cpgs_noNAs_2021[,victim_name],
                                             beta_offender_1=betas_autosomal_cpgs_noNAs_2021[,offender_name],
                                             ratio_victim = ratio_victim, ratio_offender = ratio_offender)
    
    #Deconvolute mixture, extract beta value offender with known beta value victim
    beta_offender_55_ratio <- mixture_deconvolution(mixture_28_55_ratio, 
                                                    beta_value_victim, 
                                                    ratio_victim = ratio_victim, ratio_offender = ratio_offender)
    
    #Correct beta values of the offender
    beta_offender_55_ratio[beta_offender_55_ratio > 1] <- 1
    beta_offender_55_ratio[beta_offender_55_ratio < 0] <- 0
    #Predict age (we duplicate beta_offender_55_ratio in the dataframe because methylclock package need at least two individuals to infer the age)
    df_for_prediction <- data.frame(rownames(betas_autosomal_cpgs_noNAs_2021),
                                    beta_offender_55_ratio,
                                    beta_offender_55_ratio)
    predicted_age_mixture <- DNAmAge(df_for_prediction, clocks = c("BLUP", "EN", "Horvath", "skinHorvath"))
    colnames(predicted_age_mixture)[1] <- "Ratio"
    predicted_age_mixture$Pair <- paste0(random_pairs[n,][1], "_", random_pairs[n,][2])
    predicted_age_mixture$Ratio <- ratio
    predicted_age_mixture$Age_gap <- round(abs(predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == victim_name, "Age"][[1]] -
                                                 predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == offender_name, "Age"][[1]]), 2)
    
    table_pred_age_main_contr <- rbind(table_pred_age_main_contr, predicted_age_mixture[1,])
  }
  
  pair <- paste0(random_pairs[n,][1], "_", random_pairs[n,][2])
  table_pred_age_main_contr[table_pred_age_main_contr$Pair == pair, c(2,3,4,5)] <- abs(table_pred_age_main_contr[table_pred_age_main_contr$Pair == pair, c(2,3,4,5)] - 
                                                                                         predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == offender_name,"Age"][[1]])
  
}


df_long <- table_pred_age_main_contr %>%
  pivot_longer(cols = -c(Ratio, Age_gap, Pair),
               names_to = "Clock", 
               values_to = "MAE")


df_long$Ratio <- factor(df_long$Ratio, levels = c("1:1", "10:1"))

df_long$Age_gap <- ifelse(df_long$Age_gap < 5, "<5", ">20")


saveRDS(df_long, paste0(results_path, "Age_prediction_Age_gap_simulations_Longitudinal_data_22-07-2025.rds"))
df_long <- readRDS(paste0(results_path, "Age_prediction_Age_gap_simulations_Longitudinal_data_22-07-2025.rds"))

df_long$Ratio <- gsub("10:1", "1:10", df_long$Ratio)
df_long$Age_gap <- gsub("<", "< ", df_long$Age_gap)
df_long$Age_gap <- gsub(">", "> ", df_long$Age_gap)
df_long$Ratio <- factor(df_long$Ratio, levels = c("1:1", "1:10"))

library(ggpubr)
age_gap_plot_2 <- ggplot(df_long, aes(x = Ratio, y = MAE, fill = Age_gap)) +
  geom_boxplot() +
  facet_wrap(~ Clock) +
  scale_y_continuous(limits = c(0, 20)) +
  xlab('Time (days)') +
  labs(
    title = "",
    x = "Offender-to-Victim Ratio",
    y = "AE (Years)",
    fill = "Chronological Age Difference"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("< 5" = "#009796", "> 20" = "#8A60B0"), limits = c("< 5", "> 20")) +
  stat_compare_means(
    aes(group = Age_gap), # Specify the groups for comparison
    method = "t.test", # You can change this to "t.test" or another test
    label = "p.format", # Display the p-value in a formatted way
    label.y = 19 # Adjust the position of the p-value labels
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),        # Plot title (currently empty)
    axis.title.x = element_text(size = 15),                   # X-axis title
    axis.title.y = element_text(size = 15),                   # Y-axis title
    axis.text.x = element_text(size = 13),                    # X-axis tick labels
    axis.text.y = element_text(size = 13),                    # Y-axis tick labels
    strip.text = element_text(size = 14, face = "bold"),      # Facet panel titles
    legend.title = element_text(size = 14),                   # Legend title
    legend.text = element_text(size = 13),
    legend.position = "top"# Legend items
  )

ggsave(paste0(results_path,"Plot_Age_difference_contributors_boxplot_R_17-01-2026.png"), 
       age_gap_plot_2, width = 10.5, height = 7, dpi = 600, bg = "white")

#### 7. Test impact difference predicted age and chronological age from DNA mixture 2021 Longitudinal Cohort (Boxplot) ---------------------------
##Prediction of age in single sources.
betas_autosomal_cpgs_noNAs_2021 <- data.frame(RowNames = rownames(betas_autosomal_cpgs_noNAs_2021), betas_autosomal_cpgs_noNAs_2021)
predicted_age <- DNAmAge(betas_autosomal_cpgs_noNAs_2021, clocks=c("Horvath", "skinHorvath",  "BLUP" ,"EN" ))
delta_predicted_age <- predicted_age
delta_predicted_age[, -1] <- abs(predicted_age[,-1] - predicted_age_metadata_2021$Age)
sum(delta_predicted_age$BLUP > 4) #14
sum(delta_predicted_age$BLUP < 1) #14


df_delta_pred <- data.frame(matrix(nrow = 28, ncol = 2))
colnames(df_delta_pred) <- c("Name1", "Name2")
df_delta_pred$Name1 <- c(delta_predicted_age$id[delta_predicted_age$BLUP > 4], delta_predicted_age$id[delta_predicted_age$BLUP < 1])
df_delta_pred$Name2 <- predicted_age_metadata_2021$Sample_ID[!predicted_age_metadata_2021$Sample_ID %in% df_delta_pred$Name1][1:28]

df_delta_pred <- df_delta_pred %>%
  left_join(predicted_age[,c(1,4)], by = c("Name1" = "id")) %>%
  rename(Pred_age_BLUP1 = BLUP) %>%
  left_join(ages_2021, by = c("Name1" = "Sample_ID")) %>%
  rename(Age1 = Age)
# Calculate the age difference
df_delta_pred <- df_delta_pred %>%
  mutate(AgeDifference = abs(Pred_age_BLUP1 - Age1))

n <- 123
set.seed(n)
n <- 8
ratio <- "1:1"
for (n in 1:nrow(df_delta_pred)){
  ratios <- c("1:1", "10:1")
  victim_name <- df_delta_pred[n,][[2]] #Change 1 or 2 for changing victim and offender
  offender_name <- df_delta_pred[n,][[1]] #Change 1 or 2 for changing victim and offender 
  
  random_signs <- sample(c(-1, 1), nrow(betas_autosomal_cpgs_noNAs_2021), replace = TRUE)
  random_errors <- random_signs * 0.01
  beta_value_victim <- betas_autosomal_cpgs_noNAs_2021[,victim_name] + random_errors
  #Create mixture
  for (ratio in ratios){
    print(paste(victim_name, offender_name, ratio))
    ratio_split <- strsplit(ratio, ":")
    ratio_victim <- as.numeric(ratio_split[[1]][1])
    ratio_offender <- as.numeric(ratio_split[[1]][2])
    
    if(victim_name == df_delta_pred[1,][[2]] & offender_name == df_delta_pred[1,][[1]] & ratio == "1:1"){
      table_pred_age_main_contr <- data.frame(Pair = character(),
                                              Ratio = character(),
                                              Age_gap = numeric(),
                                              Horvath = numeric(),
                                              BLUP = numeric(),
                                              EN = numeric(),
                                              skinHorvath = numeric())
    }
    mixture_28_55_ratio <- mixture_generator(beta_victim=betas_autosomal_cpgs_noNAs_2021[,victim_name],
                                             beta_offender_1=betas_autosomal_cpgs_noNAs_2021[,offender_name],
                                             ratio_victim = ratio_victim, ratio_offender = ratio_offender)
    
    #Deconvolute mixture, extract beta value offender with known beta value victim
    beta_offender_55_ratio <- mixture_deconvolution(mixture_28_55_ratio, 
                                                    beta_value_victim, 
                                                    ratio_victim = ratio_victim, ratio_offender = ratio_offender)
    
    #Correct beta values of the offender
    beta_offender_55_ratio[beta_offender_55_ratio > 1] <- 1
    beta_offender_55_ratio[beta_offender_55_ratio < 0] <- 0
    #Predict age (we duplicate beta_offender_55_ratio in the dataframe because methylclock package need at least two individuals to infer the age)
    df_for_prediction <- data.frame(rownames(betas_autosomal_cpgs_noNAs_2021),
                                    beta_offender_55_ratio,
                                    beta_offender_55_ratio)
    predicted_age_mixture <- DNAmAge(df_for_prediction, clocks = c("BLUP", "EN", "Horvath", "skinHorvath"))
    colnames(predicted_age_mixture)[1] <- "Ratio"
    predicted_age_mixture$Pair <- paste0(df_delta_pred[n,][2], "_", df_delta_pred[n,][1])
    predicted_age_mixture$Ratio <- ratio
    predicted_age_mixture$Age_gap <- round(df_delta_pred[df_delta_pred$Name1 == offender_name & df_delta_pred$Name2 == victim_name, "AgeDifference"], 2)
    
    table_pred_age_main_contr <- rbind(table_pred_age_main_contr, predicted_age_mixture[1,])
  }
  
  pair <- paste0(df_delta_pred[n,][2], "_", df_delta_pred[n,][1])
  table_pred_age_main_contr[table_pred_age_main_contr$Pair == pair, c(2,3,4,5)] <- abs(table_pred_age_main_contr[table_pred_age_main_contr$Pair == pair, c(2,3,4,5)] - 
                                                                                         predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == offender_name,"Age"][[1]])  #BLUP if you want to substract the predicted age of the single source or #Age for chronological age 
  
}



df_long <- table_pred_age_main_contr %>%
  pivot_longer(cols = -c(Ratio, Age_gap, Pair),
               names_to = "Clock", 
               values_to = "MAE")

df_long$Ratio <- factor(df_long$Ratio, levels = c("1:1", "10:1"))
df_long$Age_gap <- ifelse(df_long$Age_gap < 1, "<1", ">4")
saveRDS(df_long, paste0(results_path, "Age_prediction_err_in_single_source_simulations_Longitudinal_data_22-07-2025.rds"))

#Only BLUP clock
df_long_BLUP <- df_long[df_long$Clock == "BLUP",]

df_to_merge <- df_delta_pred[df_delta_pred$AgeDifference < 1 | df_delta_pred$AgeDifference > 4, c("Name1","AgeDifference")]
df_to_merge <- df_to_merge[df_to_merge$Name1 %in% substr(df_long_BLUP$Pair, 13, 23),]
df_to_merge$Ratio <- "Single_source"
colnames(df_to_merge)[2] <- "MAE"
df_to_merge$Age_gap <- ifelse(df_to_merge$MAE < 1, "<1", ">4")
df_to_merge$Clock <- "BLUP"

df_to_merge <- df_to_merge[match(df_to_merge$Name1, unique(substr(df_long_BLUP$Pair, 13, 23))),]
df_to_merge$Name1 <- unique(df_long_BLUP$Pair)
colnames(df_to_merge)[1] <- "Pair"
df_to_merge <- df_to_merge[,c(3,1,4,5,2)]
df_long_BLUP <- rbind(df_to_merge, df_long_BLUP)


df_long_BLUP$Age_gap <- gsub("<", "< ", df_long_BLUP$Age_gap)
df_long_BLUP$Age_gap <- gsub(">", "> ", df_long_BLUP$Age_gap)
df_long_BLUP$Ratio[df_long_BLUP$Ratio == "Single_source"] <- "Single-Source Sample" 
df_long_BLUP$Ratio[df_long_BLUP$Ratio == "10:1"] <- "1:10" 
df_long_BLUP$Ratio <- factor(df_long_BLUP$Ratio, levels = c("Single-Source Sample", "1:1", "1:10"))


age_gap_BLUP_ss_plot <- ggplot(df_long_BLUP, aes(x = Ratio, y = MAE, color = Age_gap)) +
  geom_point() + 
  geom_line(aes(group = Pair)) +  # Group lines by a unique ID
  scale_y_continuous(limits = c(0, 15)) +
  xlab('Time (days)') +
  labs(
    title = "",
    x = "Offender-to-Victim Ratio",
    y = "AE (Years)",
    color = " Absolute Difference Between Predicted \n Age and Chronological age (Offender)"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),                   # X-axis title
    axis.title.y = element_text(size = 14),                   # Y-axis title
    axis.text.x = element_text(size = 13),                    # X-axis tick labels
    axis.text.y = element_text(size = 13),                    # Y-axis tick labels
    legend.title = element_text(size = 13),                   # Legend title
    legend.text = element_text(size = 13),
  ) +
  scale_color_manual(values = c("< 1" = "#88BDE6", "> 4" = "#FBB258"), limits = c("< 1", "> 4")) +
  stat_compare_means(
    aes(x = Ratio, y = MAE, group = Age_gap),
    method = "t.test",
    label = "p.format",
    comparisons = NULL,
    label.y = 14,
    inherit.aes = FALSE,
    size=4.5
  )


age_gap_BLUP_ss_plot 

ggsave(paste0(results_path,"Plot_pred_chron_age_difference_offender_sing_source_lines_22-07-2025.png"), 
       age_gap_BLUP_ss_plot, width = 8.7, height = 6.4, dpi = 600, bg = "white")


#Use statistics to verify that the two groups produce statistically different results
df_long_BLUP
#Test if data have a normal distribution
shapiro.test(df_long_BLUP[df_long_BLUP$Ratio == "10:1" & df_long_BLUP$Age_gap == "< 1", "MAE"]) #W = 0.90703, p-value = 0.1427
shapiro.test(df_long_BLUP[df_long_BLUP$Ratio == "10:1" & df_long_BLUP$Age_gap == "> 4", "MAE"]) #W = 0.91757, p-value = 0.2029





#### 8. Test impact of age differences on age prediction of mixture with suspect-to-victim 1:10 using Pearson correlation -------------------------
predicted_age_metadata_2021
predicted_age_metadata_2021$AE <- round(abs(predicted_age_metadata_2021$Age - predicted_age_metadata_2021$BLUP), 2)

individuals <- colnames(betas_autosomal_cpgs_noNAs_2021)
shuffled_vector <- sample(individuals)

# Create random pairs
random_pairs <- matrix(individuals, ncol = 2, byrow = TRUE)

#Upload cpg of epigenetic clock to asses reconstruction of DNAm profile
checkClocks(betas_autosomal_cpgs_noNAs_2021)

# Display the pairs
random_pairs
predicted_age_metadata_2021 <- as.data.frame(predicted_age_metadata_2021)

#Create empty dataset
table_pred_age_main_contr <- data.frame(age_gap = character(),
                                        BLUP_ss = numeric(),
                                        EN_ss = numeric(),
                                        Horvath_ss = numeric(),
                                        skinHorvath_ss = numeric(),
                                        BLUP = numeric(),
                                        EN = numeric(),
                                        Horvath = numeric(),
                                        skinHorvath = numeric())

n <- 1
f <- 123
for (n in 1:nrow(random_pairs)){
  victim_name <- random_pairs[n,][1] 
  offender_name <- random_pairs[n,][2] 
  
  for (i in c(1:3)){
    ratio <-  "10:1"
    err <- 0.01 
    set.seed(f) 
      #Create mixture
      ratio_split <- strsplit(ratio, ":")
      ratio_victim <- as.numeric(ratio_split[[1]][1])
      ratio_offender <- as.numeric(ratio_split[[1]][2])

      print(paste(n, i, ratio, err))
      random_signs <- sample(c(-1, 1), nrow(betas_autosomal_cpgs_noNAs_2021), replace = TRUE)
      random_errors <- random_signs * err
      beta_value_victim <- betas_autosomal_cpgs_noNAs_2021[,victim_name] + random_errors
      beta_value_victim[beta_value_victim > 1] <- 1
      beta_value_victim[beta_value_victim < 0] <- 0
      mixture_28_55_ratio <- mixture_generator(beta_victim=betas_autosomal_cpgs_noNAs_2021[,victim_name],
                                               beta_offender_1=betas_autosomal_cpgs_noNAs_2021[,offender_name],
                                               ratio_victim = ratio_victim, ratio_offender = ratio_offender)
      
      #Deconvolute mixture, extract beta value offender with known beta value victim
      beta_offender_55_ratio <- mixture_deconvolution(mixture_28_55_ratio, 
                                                      beta_value_victim, 
                                                      ratio_victim = ratio_victim, ratio_offender = ratio_offender)
      
      #Correct beta values of the offender
      beta_offender_55_ratio[beta_offender_55_ratio > 1] <- 1
      beta_offender_55_ratio[beta_offender_55_ratio < 0] <- 0
      
      #Predict age (we duplicate beta_offender_55_ratio in the dataframe because methylclock package need at least two individuals to infer the age)
      df_for_prediction <- data.frame(rownames(betas_autosomal_cpgs_noNAs_2021),
                                        beta_offender_55_ratio,
                                        beta_offender_55_ratio)
      predicted_age_mixture <- DNAmAge(df_for_prediction, clocks = c("BLUP", "EN", "Horvath", "skinHorvath"))
      ae_predicted_age_mixture <- abs(predicted_age_mixture[1,-1] - predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == offender_name,"Age"])
      if (i ==1 ){
        predicted_age_mixture_mean <- ae_predicted_age_mixture
      } else {
        predicted_age_mixture_mean <- predicted_age_mixture_mean + ae_predicted_age_mixture
      }
  }
        
    #Calculate age gap and age in single-source samples
    age_gap <- abs(predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == victim_name,"Age"] - predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == offender_name,"Age"])
    predicted_age_ss <- abs(predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == offender_name, c("BLUP", "EN", "Horvath", "skinHorvath")] -
      predicted_age_metadata_2021[predicted_age_metadata_2021$Sample_ID == offender_name,"Age"])
    colnames(predicted_age_ss) <- c("BLUP_ss", "EN_ss", "Horvath_ss", "skinHorvath_ss")
    info_pair <- cbind(age_gap, predicted_age_ss, predicted_age_mixture_mean)
    
    table_pred_age_main_contr <- rbind(table_pred_age_main_contr, info_pair)
    
    }
  
table_pred_age_main_contr[, c(6:9)] <- table_pred_age_main_contr[, c(6:9)] / 3

saveRDS(table_pred_age_main_contr, paste0(results_path, "Age_gap_and_ss_for_pearson_corr_simulations_Longitudinal_data_01-10-2025.rds"))

table_pred_age_main_contr <- readRDS(paste0(results_path, "Age_gap_and_ss_for_pearson_corr_simulations_Longitudinal_data_01-10-2025.rds"))

table_pred_age_main_contr_long <- table_pred_age_main_contr[, -c(2,3,4,5)] %>%
  pivot_longer(cols = -c(age_gap),
               names_to = "Clock", 
               values_to = "Errors")

#calculate correlation per clock
cor_df <- table_pred_age_main_contr_long %>%
  group_by(Clock) %>%
  summarise(
    cor_res = list(cor.test(age_gap, Errors, method = "pearson")),
    .groups = "drop"
  ) %>%
  mutate(
    r = sapply(cor_res, function(x) round(x$estimate, 2)),
    p = sapply(cor_res, function(x) signif(x$p.value, 3)),
    label = paste0("r = ", r, ", p = ", p)
  )

# Choose where to place labels (adjust y/x as needed)
cor_df$x <- -0.5    # for example, left side
cor_df$y <- 19   # near the top
cor_df$Clock <- paste(cor_df$Clock, "clock")

# Factor levels for ratio (legend order)
table_pred_age_main_contr_long$Clock <- paste(table_pred_age_main_contr_long$Clock, "clock")

# Plot
library(ggpubr) 
age_gap_plot <- ggplot(
  table_pred_age_main_contr_long,
  aes(x = age_gap, y = Errors)
) +
  geom_point(size = 2) +
  facet_wrap(
    ~ Clock,
    labeller = labeller(
      Clock = c(
        "Horvath clock"      = "Horvath clock (353 CpGs)",
        "skinHorvath clock"  = "SkinHorvath clock (391 CpGs)",
        "EN clock"           = "EN clock (514 CpGs)",
        "BLUP clock"         = "BLUP clock (319,607 CpGs)"
      )
    )) +
  geom_smooth(method = "lm", se = FALSE, color = "#63a7ff", size=1.5) +
  geom_text(
    data = cor_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.9,
    hjust = 0
  ) + 
  theme_minimal() +
  coord_cartesian(ylim = c(0, 20)) +
  scale_x_continuous(expand = c(0.1, 0.1)) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0),
    legend.position = "top",
    panel.grid.major.x = element_line(color = "grey", size = 0.4),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "lightgrey"),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing.x = unit(1.5, "lines"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    plot.margin = unit(c(1, 1.5, 1, 1), "lines")
  ) +
  labs(
    title = "",
    x = "Chronological Age Difference (years)",
    y = "Absolute Error (years)"
  )
age_gap_plot

ggsave(paste0(results_path,"/Age_gap_pearson_corr_simulations_R_17-01-2026.png"), 
       age_gap_plot, width = 10, height = 8, dpi = 600, bg = "white")


#### 9. Test impact of age prediction in single source sample 1:10 using Pearson correlation ----------------
#The dataset used in this section is the same of that one generated in point 8

ss_error_age <- table_pred_age_main_contr[, -1] %>%
  dplyr::rename(
    Age_BLUP = BLUP_ss,
    Age_EN = EN_ss,
    Age_Horvath = Horvath_ss,
    Age_skinHorvath = skinHorvath_ss,
    Error_BLUP = BLUP,
    Error_EN = EN,
    Error_Horvath = Horvath,
    Error_skinHorvath = skinHorvath
  )


ss_error_age <- table_pred_age_main_contr[, -1] %>%
  rename(
    Age_BLUP = BLUP_ss)



ss_error_age_long <- ss_error_age %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c(".value", "Clock"),
    names_sep = "_"
  )

# Factor levels for ratio (legend order)
ss_error_age_long$Clock <- paste(ss_error_age_long$Clock, "clock")

#Calculate correlation
cor_df <- ss_error_age_long %>%
  group_by(Clock) %>%
  summarise(
    cor_res = list(cor.test(Age, Error, method = "pearson")),
    .groups = "drop"
  ) %>%
  mutate(
    r = sapply(cor_res, \(x) round(x$estimate, 2)),
    p = sapply(cor_res, \(x) signif(x$p.value, 3)),
    label = paste0("r = ", r, ", p = ", p)
  )

# Choose where to place labels (adjust y/x as needed)
cor_df$x <- 6   # for example, left side
cor_df$y <- 24   # near the top


age_ss_plot <- ggplot(ss_error_age_long, aes(x = Age, y = Error)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "#63a7ff") +
  facet_wrap(
    ~ Clock,
    labeller = labeller(
      Clock = c(
        "Horvath clock"      = "Horvath clock (353 CpGs)",
        "skinHorvath clock"  = "SkinHorvath clock (391 CpGs)",
        "EN clock"           = "EN clock (514 CpGs)",
        "BLUP clock"         = "BLUP clock (319,607 CpGs)"
      )
    )) +
  geom_text(
    data = cor_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1, vjust = 1, size = 4
  ) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 25)) +
  scale_x_continuous(expand = c(0.1, 0.1)) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0),
    legend.position = "top",
    panel.grid.major.x = element_line(color = "grey", size = 0.4),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "lightgrey"),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing.x = unit(1.5, "lines"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    plot.margin = unit(c(1, 1.5, 1, 1), "lines")
  ) +
  labs(
    title = "",
    x = "Absolute Error in Single-source Sample (years)",
    y = "Absolute Error in Mixture (years)"
  )

ggsave(paste0(results_path,"/Age_single_source_pearson_corr_simulations_R_17-01-2026.png"), 
       age_ss_plot, width = 10, height = 8, dpi = 600, bg = "white")




#### 10. Check missing CpGs EPIC v1.0 ------------------------------------------------
Zhou_probe_annotation_EPIC_v1 <- read_tsv(paste0(annotation_files_path, "EPIC+.hg38.manifest.tsv.gz")) #Two different annotation files, this one contains annotation of probes with probes specification.
missing_df <- data.frame(matrix(1, ncol=4))
colnames(missing_df) <- c("BLUP", "EN", "Horvath", "skinHorvath")

sample <- colnames(beta_values_autosomal_cg_noNAs)[1]
beta_values_autosomal_cg_noNAs <- data.frame(RowNames = rownames(beta_values_autosomal_cg_noNAs), beta_values_autosomal_cg_noNAs)
cpgs_missing_analysis <- cpgs_missing_analysis[c(1,4, 7,8)]

cpgs.missing <- checkClocks(beta_values_autosomal_cg_noNAs)
missing_df[1,] <- c(length(cpgs.missing[[8]]), length(cpgs.missing[[9]]),
                    length(cpgs.missing[[1]]), length(cpgs.missing[[4]]))


rownames(missing_df) <- "N. of missing cpgs" 

missing_df <- data.frame(RowNames = rownames(missing_df), missing_df)
write_xlsx(missing_df, paste0(brando_path, "/Results/2_Age_prediction/Missing_CpGs_clocks_EPICv1.xlsx"))


Zhou_probe_annotation_EPIC_v1$Probe_ID <- substr(Zhou_probe_annotation_EPIC_v1$Probe_ID, 1, 10)

#Missing CpG BLUP
319607 - sum(coefBLUP$CpGmarker %in% Zhou_probe_annotation_EPIC_v1$Probe_ID)
#Missing CpG EN
514 - sum(coefEN$CpGmarker %in% Zhou_probe_annotation_EPIC_v1$Probe_ID)
#Missing CpG Horvath
353 - sum(coefHorvath$CpGmarker %in% Zhou_probe_annotation_EPIC_v1$Probe_ID)
#Missing CpG skinHorvath
391 - sum(coefSkin$CpGmarker %in% Zhou_probe_annotation_EPIC_v1$Probe_ID)

Zhou_probe_annotation_EPIC_v2 <- read_tsv(paste0(annotation_files_path, "EPICv2.hg38.manifest.tsv.gz"))
Zhou_probe_annotation_EPIC_v2$Probe_ID <- substr(Zhou_probe_annotation_EPIC_v2$Probe_ID, 1, 10)

#Missing CpG BLUP
319607 - sum(coefBLUP$CpGmarker %in% Zhou_probe_annotation_EPIC_v2$Probe_ID)
#Missing CpG EN
514 - sum(coefEN$CpGmarker %in% Zhou_probe_annotation_EPIC_v2$Probe_ID)
#Missing CpG Horvath
353 - sum(coefHorvath$CpGmarker %in% Zhou_probe_annotation_EPIC_v2$Probe_ID)
#Missing CpG skinHorvath
391 - sum(coefSkin$CpGmarker %in% Zhou_probe_annotation_EPIC_v2$Probe_ID)

cpgs.missing <- checkClocks(Zhou_probe_annotation_EPIC_v1[,"Probe_ID"])
c(length(cpgs.missing[[8]]), length(cpgs.missing[[9]]),
                    length(cpgs.missing[[1]]), length(cpgs.missing[[4]]))

