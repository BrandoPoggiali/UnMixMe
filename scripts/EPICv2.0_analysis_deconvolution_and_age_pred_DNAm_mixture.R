# Project: Age prediction in DNAm mixture.
# Script Title: Analysis of EPIC v2.0 array data, deconvolution of DNAm mixture profile and age prediction.
# Description: The script has the aim of analysis the EPIC v.20 data related to the DNA mixture created with different ratios of two individuals.
# Then we will try to deconvolute the DNAm mixture and we the deconvolut DNAm profile we will infer DNAm methylation.
# We will try to evaluate the accuracy of the method and in case improve it.
# 
# Author: Brando Poggiali
# Date: 18-07-2024

#### 1.Setting libraries and paths. ------------------------------------------------
install.packages("scrime")
BiocManager::install("minfi")
BiocManager::install("methylclock")
install.packages("gtools") 

#Load package
library(sesame)
library(sesameData)
library(SummarizedExperiment)
library(readxl)
library(tibble)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(factoextra)
library(sva)
library(methylclock)
library(writexl)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra) 
library(reshape)
library(gtools)
library(minfi)

dimsesameDataCache() #In case you install or update SeSAMe

idat_dir <- "G:/FAELLES/Dokumenter/BRP/local_data"
brando_path <- "/mnt/ngs/projects/age_prediction_EPICv2_DNA_mixture/users/lfw156"
results_path <- "/mnt/ngs/projects/age_prediction_EPICv2_DNA_mixture/users/lfw156/Results"
annotation_files_path <- "/home/ri-domain.local/lfw156/EPIC_annotations_files/"

setwd(brando_path)

#### 2. Upload original data & metadata. ----------------------------------------------------------------------------------
Zhou_probe_annotation_EPIC_v2 <- read_tsv(paste0(annotation_files_path, "EPICv2.hg38.manifest.tsv.gz")) #Two different annotation files, this one contains annotation of probes with probes specification.
#Zhou_probe_annotation_EPIC_v2 <- read_tsv(paste0(annotation_files_path, "EPICv2.hg38.manifest.gencode.v41.tsv.gz"))
 
betas <-  openSesame(idat_dir, prep="QCDPB", func=getBetas) # Add collapseToPfx = TRUE if you immediately want to collapse the probes 
metadata <- read_xlsx(path = paste0(brando_path,"/Metadata/20240417_SampleSheet.xlsx"), sheet = "iScan_codes")
metadata_age <- read_xlsx(path = paste0(brando_path,"/Metadata/donor_ages_DNAm_mixture_Alberte_correct.xlsx"))
STR_ratio <- read_xlsx(path = paste0(brando_path,'/Results/0_Ratio_quantification/AgePredictionMixtures/AgePredictionMixtures_STRmixSummary.xlsx'))

#### 3. Tidy up dataset. ----------------------------------------------------------------------------------
metadata$Chip_pos <- paste0(metadata$Sentrix_ID, "_", metadata$Sentrix_position)
metadata$Chip_pos == colnames(betas) #check if order of name is the same

metadata$Individuals <- substring(metadata$Sample_name, 1, 5)
metadata$Individual_1 <- substring(metadata$Sample_name, 1, 2)
metadata$Individual_2 <- substring(metadata$Sample_name, 4, 5)

metadata$Ratio <- substring(metadata$Sample_name, 7, 10)
metadata$Ratio <- ifelse(substring(metadata$Ratio, 4, 4) == "_", substring(metadata$Ratio, 1, 3),metadata$Ratio) # Remove last "_" from ratios element

metadata$Ratio_Individual_1 <- gsub("_", "", substring(metadata$Ratio, 1, 2))
metadata$Ratio_Individual_2 <- gsub("_", "", substring(metadata$Ratio, nchar(metadata$Ratio) - 1, nchar(metadata$Ratio)))

metadata$Individual_1_age <- NA
metadata$Individual_2_age <- NA
ind="AA"

for (ind in unique(metadata$Individual_1)){
  print(ind)
  age <- metadata_age[metadata_age$Donor == ind, "Age"]
  metadata[metadata$Individual_1 == ind, "Individual_1_age"] <- age
  metadata[metadata$Individual_2 == ind, "Individual_2_age"] <- age
}


metadata$Duplicate <- ifelse(nchar(metadata$Sample_name) > 5,
                             substring(metadata$Sample_name, nchar(metadata$Sample_name), nchar(metadata$Sample_name)),
                             NA)

str(metadata)

colnames(betas) <- metadata$Sample_name

STR_ratio <- STR_ratio[, c(1:7)]
STR_ratio$Ratio <- STR_ratio$Component1 / STR_ratio$Component2 


#Save data
setwd(brando_path)
saveRDS(betas, paste0(brando_path, "/beta_values_age_pred_DNAm_mixture_18-07-2024.rds"))
write.table(metadata, file.path(brando_path, "Metadata/metadata_complete_30-08-2025.tsv"), sep="\t")

#### 4. Upload of RSD data after first analysis --------------------------------------------------------------------
betas <- readRDS(paste0(brando_path,"/Data/beta_values_age_pred_DNAm_mixture_18-07-2024.rds")) # no normalize for batch effect
metadata <- read.table(file = paste0(brando_path,'/Metadata/metadata_complete_30-08-2025.tsv'), sep = '\t', header = TRUE)
STR_ratio <- read_xlsx(path = paste0(brando_path,'/Results/0_Ratio_quantification/AgePredictionMixtures/AgePredictionMixtures_STRmixSummary.xlsx'))
STR_ratio <- STR_ratio[, c(1:7)]
STR_ratio$Ratio <- STR_ratio$Component1 / STR_ratio$Component2
betas_cg_autosomal_collps_no_NAs <- readRDS(file.path(brando_path, "/Data/betas_cg_autosomal_collps_30-09-2024.rds"))

#### 5. Quality control. ----------------------------------------------------------------------------------

## 4.1 Check Sesame quality control statistics
#The quality control steps are performed after applying the pre-processing steps "QCDPB" on the sesame package
qcs <- openSesame(idat_dir, prep="QCDPB", func=sesameQC_calcStats) 
qc_df <- do.call(rbind, lapply(qcs, as.data.frame))
#rownames(qc_df) == metadata$Chip_pos #check if order of name is the same
qc_df$Sample_code <- metadata$Sample_name


background_intensities_plot <- ggplot(qc_df, aes(x = mean_oob_grn, y= mean_oob_red, label = Sample_code)) + #Background
  geom_point() + geom_text(hjust = -0.3, vjust = 0.3, position=position_jitter()) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted') +
  xlab('Mean green background intesity') + ylab('Mean red background intensity') + 
  scale_color_manual(values = c("no_degraded" = "darkgrey", "350" = "black", "230" = "#f2d93a", "165" = "#39a2a4", "95" = "#f20000")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,1500)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1500)) + theme_bw() 

ggsave(paste0(results_path,"/1_Quality_control/Green_and_red_background_intensities_18-07-2024.png"), 
       background_intensities_plot, width = 6, height = 4.5)

Mean_intensities_plot <- ggplot(qc_df) + #Mean intensity
  geom_bar(aes(Sample_code, mean_intensity), stat='identity') +
  xlab('Sample Name') + ylab('Mean Intensity') +
  ylim(0,18000) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df$Sample_code) 

ggsave(paste0(results_path,"/1_Quality_control/Mean_intensities_plot_18-07-2024.png"), 
       Mean_intensities_plot, width = 6, height = 4.5)

Ratio_red_to_green_median_int_plot <- ggplot(qc_df) + #Ratio red to green background
  geom_bar(aes(Sample_code, medR / medG), stat='identity') +
  xlab('Sample Name') + ylab('Ratio of Red to Green median Intens.') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df$Sample_code)

ggsave(paste0(results_path,"/1_Quality_control/Ratio_red_to_green_median_int_plot_18-07-2024.png"), 
       Ratio_red_to_green_median_int_plot, width = 6, height = 4.5)

Number_detected_probes_plot <- ggplot(qc_df) + #Number detected Probes
  geom_bar(aes(Sample_code, num_dt), stat='identity') +
  xlab('Sample Name') + ylab('Number of detected probes') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df$Sample_code) +  #ylim(c(0,1000000)) +
  #geom_hline(yintercept = 937690, linetype = "dashed", color = "black", size = 0.4) +
  scale_y_continuous(breaks = c(0, 250000, 500000, 750000, 937690)) +
  labs(title="Number of detected probes")

ggsave(paste0(results_path,"/1_Quality_control/Number_detected_probes_plot_18-07-2024.png"), 
       Number_detected_probes_plot, width = 6, height = 4.5)


Number_of_NAs_plot <- ggplot(qc_df) + #Number of NAs
  geom_bar(aes(Sample_code, num_na_cg), stat='identity') +
  xlab('Sample Name') + ylab('Number of NAs') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df$Sample_code) +
  labs(title="Number of NAs") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0(results_path,"/1_Quality_control/Number_of_NAs_plot_18-07-2024.png"), 
       Number_of_NAs_plot, width = 6, height = 4.5)


bisulfite_conv_efficiency_plot <- ggplot(qc_df) + #Bisulfite conversion efficiency
  geom_bar(aes(Sample_code, frac_unmeth_ch), stat='identity') +
  xlab('Sample Name') + ylab('Fraction of unmethylated ch probes') +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_x_discrete(limits=qc_df$Sample_code)

ggsave(paste0(results_path,"/1_Quality_control/bisulfite_conv_efficiency_plot_18-07-2024.png"), 
       bisulfite_conv_efficiency_plot, width = 6, height = 4.5)

objects_to_remove <- grep("plot", ls(), value = TRUE) # remove plots from environment
rm(list = objects_to_remove)


##4.2 Check concordance of samples gender.
Chromosome_X_and_Y_probes <- Zhou_probe_annotation_EPIC_v2 %>%
  filter(CpG_chrm %in% c("chrX", "chrY")) %>%  # Just select the chromosomes of interests example: c("chrX", "chrY")
  pull(Probe_ID)

beta_values_chr_X_Y_bloodstains_wide <- beta_values_bloodstains_wide[rownames(beta_values_bloodstains_wide) %in%
                                                                       Chromosome_X_and_Y_probes,]
#PCA
PCA_result <- prcomp(t(na.omit(beta_values_chr_X_Y_bloodstains_wide)), center = TRUE, scale. = TRUE) #creates the principal components
PCA_result_scores <- as.data.frame(PCA_result$x)
percentage <- round(PCA_result$sdev^2 / sum(PCA_result$sdev^2) * 100, 2)
percentage <- paste(colnames(PCA_result_scores), "(", paste( as.character(percentage), "%", ")", sep="") )

PCA_result_scores$Sample_name <- as.factor(metadata_bloodstains$Sample_name)
PCA_result_scores$Gender <- as.factor(metadata_bloodstains$Gender)

PCA_chr_X_Y <- ggplot(PCA_result_scores,aes(x=PC1,y=PC2,color = Gender)) +
  geom_point(size=1) + 
  theme_bw() +
  xlab(percentage[1]) + ylab(percentage[2]) +
  geom_text(size= 2.5, aes(label = Sample_name),nudge_y = -1.5, show.legend = FALSE)

ggsave(paste0(brando_path,"/Results/Quality_control/PCA_chromosome_X_and_Y_01-29-2024.png"), 
       PCA_chr_X_Y)


##4.4 Removal of sex chromosomes
Chromosome_X_and_Y_probes <- Zhou_probe_annotation_EPIC_v2 %>%
  filter(CpG_chrm %in% c("chrX", "chrY")) %>%  # Just select the chromosomes of interests example: c("chrX", "chrY")
  pull(Probe_ID)

Chromosome_Y_probes <- Zhou_probe_annotation_EPIC_v2 %>%
  filter(CpG_chrm %in% c("chrY")) %>%  # Just select the chromosomes of interests example: c("chrX", "chrY")
  pull(Probe_ID)

betas_X_Y_chr <- betas[rownames(betas) %in% Chromosome_X_and_Y_probes,]
betas_Y_chr <- betas[rownames(betas) %in% Chromosome_Y_probes,]
betas_autosomal <- betas[!rownames(betas) %in% Chromosome_X_and_Y_probes,]

PCA_result <- prcomp(t(na.omit(betas_X_Y_chr)), center = TRUE, scale. = TRUE) #creates the principal components
PCA_result_scores <- as.data.frame(PCA_result$x)
percentage <- round(PCA_result$sdev^2 / sum(PCA_result$sdev^2) * 100, 2)
percentage <- paste(colnames(PCA_result_scores), "(", paste( as.character(percentage), "%", ")", sep="") )

PCA_result_scores$Ratio <- as.factor(metadata$Ratio)
PCA_result_scores$Individual_1 <- as.factor(metadata$Individual_1)
PCA_result_scores$Individuals <- as.factor(metadata$Individuals)

PCA <- ggplot(PCA_result_scores,aes(x = PC1, y = PC2, color = Individuals, label=Ratio)) +
  geom_point(size=3) + 
  geom_text(hjust = -0.3, vjust = 0.3, position=position_jitter()) +
  theme_bw() +
  xlab(percentage[1]) + ylab(percentage[2]) +
  labs(color = "DNA size", shape = "DNA input", title= "") +
  theme(plot.title = element_text(hjust = 0.5)) #+

ggsave(paste0(results_path,"/1_Quality_control/PCA_on_subset_samples_200K_cg_24-03-2024.png"), 
       PCA, width = 7, height = 5.1)

saveRDS(beta_values_autosomal_bloodstains_wide, file.path(brando_path, "/Data/beta_values_autosomal_bloodstains_wide_30-01-2024.rds"))


##4.5 Removal of probes with NAs
#beta_values_autosomal_noNAs_bloodstains_wide <- na.omit(beta_values_autosomal_bloodstains_wide)
#saveRDS(beta_values_autosomal_noNAs_bloodstains_wide, file.path(brando_path, "/Data/beta_values_autosomal_noNAs_bloodstains_wide_02-02-2024.rds"))

## 4.6 Keeping only cg  sites
cg_ch_probes <- substr(rownames(betas_autosomal),1,2) %in% c("cg")
betas_cg_autosomal <- betas_autosomal[cg_ch_probes,]

##4.7 Check Beta value distribution (This step require an high amount of memory).
beta_vals_autosomal_long <- beta_values_autosomal_cgch_noNAs_bloodstains_wide %>% as.tibble() %>%
  add_column(CpG = rownames(beta_values_autosomal_cgch_noNAs_bloodstains_wide)) %>%
  gather(Sample_name, beta_value, 1:28) 
colnames(beta_vals_autosomal_long)[1] <- "Probe_ID"

beta_vals_autosomal_long <- beta_vals_autosomal_long %>%
  left_join(Zhou_probe_annotation_EPIC_v2[,c("Probe_ID", "type")], by = "Probe_ID") 

samples_distribution <- ggplot(data = beta_vals_autosomal_long) +
  geom_density(aes(x = beta_value, color = Sample_name)) #+
geom_density(aes(x = beta_value))

ggsave(paste0(brando_path,"/Results/1_Quality_control/Samples_beta_distribution_06-02-2024.png"), 
       samples_distribution, width = 8, height = 6)

all_samples_distribution <- ggplot(data = beta_vals_autosomal_long) +
  geom_density(aes(x = beta_value))

ggsave(paste0(brando_path,"/Results/1_Quality_control/All_samples_beta_distribution_06-02-2024.png"), 
       all_samples_distribution, width = 8, height = 6)


## 4.8 Collapse duplicates probes
betas_cg_autosomal_collps <- betasCollapseToPfx(as.matrix(betas_cg_autosomal)) 
betas_cg_autosomal_collps_no_NAs <- na.omit(betas_cg_autosomal_collps)
betas_cg_autosomal_collps_no_NAs <- as.data.frame(betas_cg_autosomal_collps_no_NAs)
saveRDS(betas_cg_autosomal_collps_no_NAs, file.path(brando_path, "/Data/betas_cg_autosomal_collps_30-09-2024.rds"))

##4.9 Check missing and excluded CpG sites
EPICv2_cpgs <-  betasCollapseToPfx(as.matrix(betas))
checkClocks(EPICv2_cpgs)

# Check missing CpGs in the EPIC v2.0 array
length(coefBLUP$CpGmarker[-1]) - sum(coefBLUP$CpGmarker[-1] %in% rownames(EPICv2_cpgs))
length(coefEN$CpGmarker[-1]) - sum(coefEN$CpGmarker[-1] %in% rownames(EPICv2_cpgs))
length(coefHorvath$CpGmarker[-1]) - sum(coefHorvath$CpGmarker[-1] %in% rownames(EPICv2_cpgs))
length(coefSkin$CpGmarker[-1]) - sum(coefSkin$CpGmarker[-1] %in% rownames(EPICv2_cpgs))

##### 6. DNAm Mixture deconvolution and Age prediction offender -----------------------------------------------------
mixture_deconvolution <- function(beta_mixture, beta_victim, proportion_victim = 1,
                                  proportion_offender = 1){
  beta_offender <- (((proportion_victim + proportion_offender) * beta_mixture) - (proportion_victim * beta_victim)) / proportion_offender
  return(beta_offender)
}

mixture_deconvolution_M_values <- function(beta_mixture, beta_victim, proportion_victim = 1, # M-values conversion
                                  proportion_offender = 1){
  M_mixture <- log2(beta_mixture / (1 - beta_mixture))
  M_victim <- log2(beta_victim / (1 - beta_victim))
  M_offender <- (((proportion_victim + proportion_offender) * M_mixture) - (proportion_victim * M_victim)) / proportion_offender
  beta_offender <- 2^M_offender / (2^M_offender + 1)
  return(beta_offender)
}

colnames(betas_cg_autosomal_collps_no_NAs) == metadata$Sample_name

df_predicted_ages <- data.frame(matrix(nrow = ncol(betas_cg_autosomal_collps_no_NAs) - 4, ncol = 5))
df_absolute_error_ages <- data.frame(matrix(nrow = ncol(betas_cg_autosomal_collps_no_NAs) - 4, ncol = 5))
df_reconstructured_DNAm_profile <- data.frame(matrix(nrow = nrow(betas_cg_autosomal_collps_no_NAs), ncol = ncol(betas_cg_autosomal_collps_no_NAs) -4))

#name <- "AE_AF_1_10_A"
n <- 1
for (name in colnames(betas_cg_autosomal_collps_no_NAs)){
  if (nchar(name) < 5){ next }
  betas_mixture <- betas_cg_autosomal_collps_no_NAs[, name]
  name_victim <- substring(name, 4, 5)
  name_offender <- substring(name, 1, 2)
  age_offender <- as.numeric(metadata[metadata$Sample_name == name, "Individual_1_age"])
  betas_victim <- betas_cg_autosomal_collps_no_NAs[, name_victim]
  
  STR_ratio_name <- STR_ratio[STR_ratio$SampleName == substring(name, 1, nchar(name) - 2), ]
  if (STR_ratio_name$Match1 == name_victim) {
    ratio_value_victim <- STR_ratio_name$Ratio
    ratio_value_offender <- 1
  } else {
    ratio_value_offender <- STR_ratio_name$Ratio
    ratio_value_victim <- 1
  }
  
  #ratio_value_victim <- as.numeric(metadata[metadata$Sample_name == name, "Ratio_Individual_2"]) 
  #ratio_value_offender <- as.numeric(metadata[metadata$Sample_name == name, "Ratio_Individual_1"]) 
  
  #Deconvolute profile
  deconv_betas_offender <- mixture_deconvolution(betas_mixture, betas_victim, 
                                         proportion_victim = ratio_value_victim, 
                                         proportion_offender = ratio_value_offender)
  
  
  #Correct beta values of the offender
  deconv_betas_offender[deconv_betas_offender > 1] <- 1
  deconv_betas_offender[deconv_betas_offender < 0] <- 0
  
  #Store reconstructed DNAm profile of suspect
  df_reconstructured_DNAm_profile[,n] <- deconv_betas_offender
  colnames(df_reconstructured_DNAm_profile)[n] <- name
  #Age Prediction
  df_for_prediction <- data.frame(rownames(betas_cg_autosomal_collps_no_NAs), 
                                  deconv_betas_offender, deconv_betas_offender)
  predicted_age_mixture <- as.data.frame(DNAmAge(df_for_prediction, clocks = c("BLUP", "EN", "Horvath", "skinHorvath")))
  predicted_age_mixture <- predicted_age_mixture[1,-1]
  absolute_error_mixture <- abs(predicted_age_mixture - age_offender)
  predicted_age_mixture$Sample_name <- name
  absolute_error_mixture$Sample_name <- name
  df_predicted_ages[n,] <- predicted_age_mixture
  df_absolute_error_ages[n,] <- absolute_error_mixture
  n <- n + 1
  }


colnames(df_predicted_ages) <- c("Horvath", "skinHorvath", "BLUP","EN", "Sample_name")
colnames(df_absolute_error_ages) <- c("AE_Horvath", "AE_skinHorvath", "AE_BLUP","AE_EN", "Sample_name")

metadata_pred_age <- merge(metadata, df_predicted_ages, by="Sample_name")
metadata_AE_pred_age <- merge(metadata, df_absolute_error_ages, by="Sample_name")

write_xlsx(metadata_pred_age, paste0(results_path, "/2_Age_prediction/Predicted_age_DNA_mixtures_02-09-2025.xlsx"))
write_xlsx(metadata_AE_pred_age, paste0(results_path, "/2_Age_prediction/Absolute_errors_predicted_age_DNA_mixtures_02-09-2025.xlsx"))

saveRDS(df_reconstructured_DNAm_profile, paste0(brando_path, "/Reconstructed_beta_values_age_pred_DNAm_mixture_02-09-2025.rds"))

## Result tables creation
write_xlsx(metadata_AE_pred_age, paste0(results_path, "/2_Age_prediction/Absolute_errors_predicted_age_DNA_mixtures_02-09-2025.xlsx"))


AE_predicted_age <- metadata_AE_pred_age[, c(1, 14, 15, 16, 17)]
AE_predicted_age$Mixture_type <- substring(AE_predicted_age$Sample_name, 1, nchar(AE_predicted_age$Sample_name) - 2)

MAE_predicted_age <- AE_predicted_age %>% # Calculate Mean Absolute Error for replicates
  group_by(Mixture_type) %>%
  summarize(across(c("AE_Horvath", "AE_skinHorvath", "AE_BLUP","AE_EN"),
                   mean, .names = "M{.col}"))

MAE_predicted_age$Mixture_type <- factor(MAE_predicted_age$Mixture_type, levels=c("AA_AB_1_1",	"AA_AB_2_1",	"AA_AB_4_1",	"AA_AB_10_1",	"AA_AB_1_2",	"AA_AB_1_4",	"AA_AB_1_10",	"AE_AF_1_1",	
                                   "AE_AF_2_1",	"AE_AF_4_1",	"AE_AF_10_1",	"AE_AF_1_2",	"AE_AF_1_4",	"AE_AF_1_10"))

MAE_predicted_age <- MAE_predicted_age[order(MAE_predicted_age$Mixture_type),]
MAE_predicted_age[,-1] <- round(MAE_predicted_age[,-1], 2)
MAE_predicted_age$Individuals <- substring(as.character(MAE_predicted_age$Mixture_type),1,5)
MAE_predicted_age$Ratio <- substring(MAE_predicted_age$Mixture_type, 7, 10)
MAE_predicted_age$Ratio <- gsub("_", ":", MAE_predicted_age$Ratio)
STR_ratio_add <- STR_ratio[,c(1, 7, 8)]
colnames(STR_ratio_add)[1] <- "Mixture_type"
MAE_predicted_age <- merge(MAE_predicted_age, STR_ratio_add, by="Mixture_type")

for (n in 1:nrow(MAE_predicted_age)){
  if (substring(MAE_predicted_age[n,1], 1, 2) == MAE_predicted_age[n,"Match1"]){
    MAE_predicted_age$STR_ratio[n] <- paste0(round(MAE_predicted_age[n,"Ratio.y"], 1), ":1")
  } else {
    MAE_predicted_age$STR_ratio[n] <- paste0("1:", round(MAE_predicted_age[n,"Ratio.y"], 1))
  }
}

MAE_predicted_age <- MAE_predicted_age[,-c(8,9)]
colnames(MAE_predicted_age)[7] <- "Ratio"

write_xlsx(MAE_predicted_age, paste0(results_path, "/2_Age_prediction/Mean_Absolute_errors_predicted_age_DNA_mixtures_STR_ratio_02-09-2025.xlsx"))


##Check missing CpGs EPIC v2.0
missing_df <- data.frame(matrix(nrow=ncol(betas_cg_autosomal_collps_no_NAs), ncol=4))
colnames(missing_df) <- c("BLUP", "EN", "Horvath", "skinHorvath")

n <- 1
sample <- colnames(betas_cg_autosomal_collps_no_NAs)[1]
betas_cg_autosomal_collps_no_NAs <- data.frame(RowNames = rownames(betas_cg_autosomal_collps_no_NAs), betas_cg_autosomal_collps_no_NAs)
for (sample in colnames(betas_cg_autosomal_collps_no_NAs[-1])){
  betas_cg_collapsed_ind <- betas_cg_autosomal_collps_no_NAs[, c("RowNames",sample)]
  #betas_cg_collapsed_ind <- betas_cg_collapsed_ind[!is.na(betas_cg_collapsed_ind[,2]),]
  
  cpgs.missing <- checkClocks(betas_cg_collapsed_ind)
  missing_df[n,] <- c(length(cpgs.missing[[8]]), length(cpgs.missing[[9]]),
                      length(cpgs.missing[[1]]), length(cpgs.missing[[4]]))
  rownames(missing_df)[n] <- sample 
  n <- n + 1
}

rownames(missing_df) <- gsub("X","",rownames(missing_df))
missing_df <- missing_df[1,]
rownames(missing_df) <- "N. of missing cpgs" 

missing_df <- data.frame(RowNames = rownames(missing_df), missing_df)
write_xlsx(missing_df, paste0(brando_path, "/Results/2_Age_prediction/Missing_CpGs_clocks.xlsx"))

##### 7. Investigation accuracy in reconstructing DNAm profile --------------------------------------------------
df_reconstructured_DNAm_profile <- readRDS(file= paste0(brando_path, "/Reconstructed_beta_values_age_pred_DNAm_mixture_02-09-2025.rds"))

#Upload clocks
checkClocks(df_reconstructured_DNAm_profile)

median_delta_betas_DNAm_profiles <- data.frame(matrix(nrow = ncol(df_reconstructured_DNAm_profile), ncol = 4))

#name <- "AA_AB_1_4_A"
n <- 1
for (name in colnames(df_reconstructured_DNAm_profile)){
  if (nchar(name) < 5){ next }
  print(name)
  name_victim <- substring(name, 4, 5)
  name_offender <- substring(name, 1, 2)
  DNAm_profile_ss <- betas_cg_autosomal_collps_no_NAs[,colnames(betas_cg_autosomal_collps_no_NAs) == name_offender]
  DNAm_profile_sample <- df_reconstructured_DNAm_profile[,colnames(df_reconstructured_DNAm_profile) == name]
  
  absolute_error_mixture <- abs(DNAm_profile_sample - DNAm_profile_ss)
  median_delta_betas_DNAm_profiles[n, 1] <- median(absolute_error_mixture[coefHorvath$CpGmarker[-1]], na.rm = TRUE)
  median_delta_betas_DNAm_profiles[n, 2] <- median(absolute_error_mixture[coefSkin$CpGmarker[-1]], na.rm = TRUE)
  median_delta_betas_DNAm_profiles[n, 3] <- median(absolute_error_mixture[coefBLUP$CpGmarker[-1]], na.rm = TRUE)
  median_delta_betas_DNAm_profiles[n, 4] <- median(absolute_error_mixture[coefEN$CpGmarker[-1]], na.rm = TRUE)
  
  median_delta_betas_DNAm_profiles$Sample_name[n] <- name
  n <- n + 1
}

colnames(median_delta_betas_DNAm_profiles) <- c("Error_Horvath", "Error_skinHorvath", "Error_BLUP", "Error_EN", "Sample_name")

median_delta_betas_DNAm_profiles$Mixture_type <- substring(median_delta_betas_DNAm_profiles$Sample_name, 1, nchar(median_delta_betas_DNAm_profiles$Sample_name) - 2)

median_delta_betas_DNAm_profiles$Mixture_type <- factor(median_delta_betas_DNAm_profiles$Mixture_type, levels=c("AA_AB_1_1",	"AA_AB_2_1",	"AA_AB_4_1",	"AA_AB_10_1",	"AA_AB_1_2",	"AA_AB_1_4",	"AA_AB_1_10",	"AE_AF_1_1",	
                                                                                  "AE_AF_2_1",	"AE_AF_4_1",	"AE_AF_10_1",	"AE_AF_1_2",	"AE_AF_1_4",	"AE_AF_1_10"))
median_delta_betas_DNAm_profiles$Replicate <- substring(median_delta_betas_DNAm_profiles$Sample_name, nchar(median_delta_betas_DNAm_profiles$Sample_name), 
                                                        nchar(median_delta_betas_DNAm_profiles$Sample_name))


median_delta_betas_DNAm_profiles <- median_delta_betas_DNAm_profiles[order(median_delta_betas_DNAm_profiles$Mixture_type),]
median_delta_betas_DNAm_profiles[,-c(5,6,7)] <- round(median_delta_betas_DNAm_profiles[,-c(5,6,7)], 4)
median_delta_betas_DNAm_profiles$Individuals <- substring(as.character(median_delta_betas_DNAm_profiles$Mixture_type),1,5)
median_delta_betas_DNAm_profiles$Ratio <- substring(median_delta_betas_DNAm_profiles$Mixture_type, 7, 10)
median_delta_betas_DNAm_profiles$Ratio <- gsub("_", ":", median_delta_betas_DNAm_profiles$Ratio)
STR_ratio_add <- STR_ratio[,c(1, 7, 8)]
colnames(STR_ratio_add)[1] <- "Mixture_type"
median_delta_betas_DNAm_profiles <- merge(median_delta_betas_DNAm_profiles, STR_ratio_add, by="Mixture_type")

for (n in 1:nrow(median_delta_betas_DNAm_profiles)){
  if (substring(median_delta_betas_DNAm_profiles[n,1], 1, 2) == median_delta_betas_DNAm_profiles[n,"Match1"]){
    median_delta_betas_DNAm_profiles$STR_ratio[n] <- paste0(round(median_delta_betas_DNAm_profiles[n,"Ratio.y"], 1), ":1")
  } else {
    median_delta_betas_DNAm_profiles$STR_ratio[n] <- paste0("1:", round(median_delta_betas_DNAm_profiles[n,"Ratio.y"], 1))
  }
}

median_delta_betas_DNAm_profiles <- median_delta_betas_DNAm_profiles[,-c(10, 11, 12)]
colnames(median_delta_betas_DNAm_profiles)[9] <- "Theoretical_ratio"

write_xlsx(median_delta_betas_DNAm_profiles, paste0(results_path, "/2_Age_prediction/Median_delta_betas_reconstructed_DNAm_profiles_03-09-2025.xlsx"))

#Plot error of reconstructed DNAm profiles
median_delta_betas_DNAm_profiles <- median_delta_betas_DNAm_profiles[,-c(6,8)]

median_delta_betas_DNAm_profiles_long <- median_delta_betas_DNAm_profiles %>%
  pivot_longer(cols = -c(Mixture_type, Theoretical_ratio, Replicate),
               names_to = "Clock", 
               values_to = "Errors")


median_delta_betas_DNAm_profiles_long$Theoretical_ratio <- factor(median_delta_betas_DNAm_profiles_long$Theoretical_ratio, levels = c("10:1", "4:1", "2:1", "1:1","1:2", "1:4", "1:10"))
median_delta_betas_DNAm_profiles_long$Pair <- substr(median_delta_betas_DNAm_profiles_long$Mixture_type, 1, 5)
median_delta_betas_DNAm_profiles_long$Clock <- paste(gsub("Error_","", median_delta_betas_DNAm_profiles_long$Clock), "clock")


#Change AE-AF in AC-AD
median_delta_betas_DNAm_profiles_long$Mixture_type <- gsub("AA_AB", "M1_F1", median_delta_betas_DNAm_profiles_long$Mixture_type)
median_delta_betas_DNAm_profiles_long$Mixture_type <- gsub("AE_AF", "M2_F2", median_delta_betas_DNAm_profiles_long$Mixture_type)

median_delta_betas_DNAm_profiles_long$Pair <- gsub("AA_AB", "M1 F1", median_delta_betas_DNAm_profiles_long$Pair)
median_delta_betas_DNAm_profiles_long$Pair <- gsub("AE_AF", "M2 F2", median_delta_betas_DNAm_profiles_long$Pair)
median_delta_betas_DNAm_profiles_long$Pair <- gsub(" ", "-", median_delta_betas_DNAm_profiles_long$Pair)
#median_delta_betas_DNAm_profiles_long$Clock <- paste(median_delta_betas_DNAm_profiles_long$Clock, "clock")
median_delta_betas_DNAm_profiles_long$Mixture_type <- substr(median_delta_betas_DNAm_profiles_long$Mixture_type, 1, 5)

max(median_delta_betas_DNAm_profiles_long$Errors)
# Define custom colors
pair_colors <- c("M1-F1" = "#028FFB", "M2-F2" = "#FB6E02")
ref_colors <- c("M1" = "#4DB1FE", "M2" = "#FD9341")




MAE_plot <- ggplot(
  median_delta_betas_DNAm_profiles_long,
  aes(
    x = Theoretical_ratio,
    y = Errors,
    color = Pair,
    shape = Replicate,   # distinguish A vs B
    group = interaction(Pair, Replicate)
  )
) +
  geom_point(size = 3, position = position_dodge(width = 0.25)) +
  
  facet_wrap(
    ~ Clock,
    labeller = labeller(
      Clock = c(
        "Horvath clock" = "Horvath clock (353 CpGs)",
        "skinHorvath clock"  = "skinHorvath clock (391 CpGs)",
        "EN clock"    = "EN clock (514 CpGs)",
        "BLUP clock"     = "BLUP clock (319,607 CpGs)"
      )
    )
  ) +
  coord_cartesian(ylim = c(0, 0.18)) +
  theme_minimal() +
  
  # Manual color scale for Pairs
  scale_color_manual(
    values = pair_colors,
    name = "Mixture"
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  ) +
  # Shapes for Replicates
  scale_shape_manual(
    values = c("A" = 16, "B" = 17), # filled circle and triangle
    name = "Replicate"
  ) +
  
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0.1),
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
  xlab("Suspect-to-Victim Ratio") +
  ylab("Median |Δβ|")

ggsave(paste0(results_path,"/2_Age_prediction/Median_delta_betas_reconstructed_DNAm_profiles_R_13-01-2026.png"), 
       MAE_plot, width = 11, height = 7, dpi = 600, bg = "white")

##### 8. Age prediction in single source samples and age prediction from mixture using single-source sample as reference-----------------------------------------------------
betas_single_source <- betas_cg_autosomal_collps_no_NAs[, c("AA","AB", "AE", "AF")]

df_for_prediction <- data.frame(rownames(betas_cg_autosomal_collps_no_NAs), betas_single_source)
predicted_age_single_source <- as.data.frame(DNAmAge(df_for_prediction, clocks = c("BLUP", "EN", "Horvath", "skinHorvath")))

Age_AA <- metadata[metadata$Sample_name == "AA", "Individual_1_age"][[1]]
Age_AB <- metadata[metadata$Sample_name == "AB", "Individual_1_age"][[1]]
Age_AE <- metadata[metadata$Sample_name == "AE", "Individual_1_age"][[1]]
Age_AF <- metadata[metadata$Sample_name == "AF", "Individual_1_age"][[1]]

AE_single_source <- predicted_age_single_source
AE_single_source[1,-1] <- abs(predicted_age_single_source[1,-1] - Age_AA)
AE_single_source[2,-1] <- abs(predicted_age_single_source[2,-1] - Age_AB)
AE_single_source[3,-1] <- abs(predicted_age_single_source[3,-1] - Age_AE)
AE_single_source[4,-1] <- abs(predicted_age_single_source[4,-1] - Age_AF)

write_xlsx(AE_single_source, paste0(results_path, "/2_Age_prediction/Absolute_errors_predicted_age_DNA_single_sources_02-09-2025.xlsx"))

#Age prediction from DNA mixture calculating the MAE using predicted age in single source samples

df_absolute_error_ages_with_single_source <- data.frame(matrix(nrow = ncol(betas_cg_autosomal_collps_no_NAs) - 4, ncol = 4))
predicted_age_ss <- c(AA=Age_AA, AB=Age_AB, AE=Age_AE, AF=Age_AF)
#name <- "AA_AB_1_4_A"
n <- 1
for (name in df_predicted_ages$Sample_name){
  if (nchar(name) < 5){ next }
  print(name)
  name_victim <- substring(name, 4, 5)
  name_offender <- substring(name, 1, 2)
  predicted_age_offender_ss <- predicted_age_single_source[predicted_age_single_source$id == name_offender,][,-1]
  predicted_age_mixture_ss <- df_predicted_ages[df_predicted_ages$Sample_name == name, -5]
  
  absolute_error_mixture <- abs(predicted_age_mixture_ss - predicted_age_offender_ss)
  df_absolute_error_ages_with_single_source[n,] <- absolute_error_mixture
  df_absolute_error_ages_with_single_source$Sample_name[n] <- name
  n <- n + 1
}
colnames(df_absolute_error_ages_with_single_source) <- c("AE_Horvath", "AE_skinHorvath", "AE_BLUP","AE_EN", "Sample_name")

metadata_AE_pred_age_ss <- merge(metadata, df_absolute_error_ages_with_single_source, by="Sample_name")

write_xlsx(metadata_AE_pred_age_ss, paste0(results_path, "/2_Age_prediction/Predicted_age_DNA_mixtures_ss_02-09-2025.xlsx"))

## Result tables creation for MAE calculated on single source sample
AE_predicted_age <- metadata_AE_pred_age_ss[, c(1, 14, 15, 16, 17)]
AE_predicted_age$Mixture_type <- substring(AE_predicted_age$Sample_name, 1, nchar(AE_predicted_age$Sample_name) - 2)

MAE_predicted_age <- AE_predicted_age %>% # Calculate Mean Absolute Error for replicates
  group_by(Mixture_type) %>%
  summarize(across(c("AE_Horvath", "AE_skinHorvath", "AE_BLUP","AE_EN"),
                   mean, .names = "M{.col}"))

MAE_predicted_age$Mixture_type <- factor(MAE_predicted_age$Mixture_type, levels=c("AA_AB_1_1",	"AA_AB_2_1",	"AA_AB_4_1",	"AA_AB_10_1",	"AA_AB_1_2",	"AA_AB_1_4",	"AA_AB_1_10",	"AE_AF_1_1",	
                                                                                  "AE_AF_2_1",	"AE_AF_4_1",	"AE_AF_10_1",	"AE_AF_1_2",	"AE_AF_1_4",	"AE_AF_1_10"))

MAE_predicted_age <- MAE_predicted_age[order(MAE_predicted_age$Mixture_type),]
MAE_predicted_age[,-1] <- round(MAE_predicted_age[,-1], 2)
MAE_predicted_age$Individuals <- substring(as.character(MAE_predicted_age$Mixture_type),1,5)
MAE_predicted_age$Ratio <- substring(MAE_predicted_age$Mixture_type, 7, 10)
MAE_predicted_age$Ratio <- gsub("_", ":", MAE_predicted_age$Ratio)
STR_ratio_add <- STR_ratio[,c(1, 7, 8)]
colnames(STR_ratio_add)[1] <- "Mixture_type"
MAE_predicted_age <- merge(MAE_predicted_age, STR_ratio_add, by="Mixture_type")

for (n in 1:nrow(MAE_predicted_age)){
  if (substring(MAE_predicted_age[n,1], 1, 2) == MAE_predicted_age[n,"Match1"]){
    MAE_predicted_age$STR_ratio[n] <- paste0(round(MAE_predicted_age[n,"Ratio.y"], 1), ":1")
  } else {
    MAE_predicted_age$STR_ratio[n] <- paste0("1:", round(MAE_predicted_age[n,"Ratio.y"], 1))
  }
}

MAE_predicted_age <- MAE_predicted_age[,-c(8,9)]
colnames(MAE_predicted_age)[7] <- "Ratio"

write_xlsx(MAE_predicted_age, paste0(results_path, "/2_Age_prediction/Mean_Absolute_errors_predicted_age_DNA_mixtures_STR_ratio_ss_02-09-2025.xlsx"))


#### 9. Plotting MAE for different DNA ratio (Plot for publication) ------------------------------------------
#This code can be used for Figure 3 and Supplementary Figure 2, you just need to cahnge the input

Predicted_age <- read_xlsx(path = paste0(brando_path,'/Results/2_Age_prediction/Mean_Absolute_errors_predicted_age_DNA_mixtures_STR_ratio_02-09-2025.xlsx')) #ermove or add ss based on the type of plot you want do make
Predicted_age_single_source <- read_xlsx(path = paste0(brando_path,'/Results/2_Age_prediction/Absolute_errors_predicted_age_DNA_single_sources_02-09-2025.xlsx'))
colnames(Predicted_age)[7] <- "Theoretical_ratio"
Predicted_age <- Predicted_age[,-6]

Predicted_age_long <- Predicted_age %>%
  pivot_longer(cols = -c(Mixture_type, Theoretical_ratio, STR_ratio),
               names_to = "Clock", 
               values_to = "MAE")


Predicted_age_long$Theoretical_ratio <- factor(Predicted_age_long$Theoretical_ratio, levels = c("10:1", "4:1", "2:1", "1:1","1:2", "1:4", "1:10"))
Predicted_age_long$Pair <- substr(Predicted_age_long$Mixture_type, 1, 5)
Predicted_age_long$Clock <- substr(Predicted_age_long$Clock, 5, length(Predicted_age_long$Clock))

#df_long$Age_gap <- sprintf("%.1f", df_long$Age_gap)

Predicted_age_single_source
Predicted_age_single_source_long <- Predicted_age_single_source %>%
  pivot_longer(cols = -c(id),
               names_to = "Clock", 
               values_to = "AE")
Predicted_age_single_source_long <- Predicted_age_single_source_long[Predicted_age_single_source_long$id %in% c("AA","AE"),]
colnames(Predicted_age_single_source_long)[1] <- "color"

Predicted_age_long$Pair <- gsub("_", " ", Predicted_age_long$Pair)

#Change AE-AF in AC-AD
Predicted_age_long$Mixture_type <- gsub("AA_AB", "M1_F1", Predicted_age_long$Mixture_type)
Predicted_age_long$Mixture_type <- gsub("AC_AD", "M2_F2", Predicted_age_long$Mixture_type)

Predicted_age_long$Pair <- gsub("AA AB", "M1 F1", Predicted_age_long$Pair)
Predicted_age_long$Pair <- gsub("AE AF", "M2 F2", Predicted_age_long$Pair)
Predicted_age_long$Pair <- gsub(" ", "-", Predicted_age_long$Pair)
Predicted_age_long$Clock <- paste(Predicted_age_long$Clock, "clock")

Predicted_age_single_source_long$color <- gsub("AA", "M1", Predicted_age_single_source_long$color)
Predicted_age_single_source_long$color <- gsub("AE", "M2", Predicted_age_single_source_long$color)
Predicted_age_single_source_long$Clock <- paste(Predicted_age_single_source_long$Clock, "clock")

# Prepare horizontal line data
hline_data <- Predicted_age_single_source_long[Predicted_age_single_source_long$color %in% c("M1", "M2"), ]
hline_data$Reference <- hline_data$color  # Will be "AA" and "AC"

# Define custom colors
pair_colors <- c("M1-F1" = "#028FFB", "M2-F2" = "#FB6E02")
ref_colors <- c("M1" = "#4DB1FE", "M2" = "#FD9341")


## To add the absolute error we need to add all AEs
Predicted_age_abs <- read_xlsx(path = paste0(brando_path,'/Results/2_Age_prediction/Absolute_errors_predicted_age_DNA_mixtures_02-09-2025.xlsx')) #ermove or add ss based on the type of plot you want do make
#Predicted_age_abs <- read_xlsx(path = paste0(brando_path,'/Results/2_Age_prediction/Predicted_age_DNA_mixtures_02-09-2025.xlsx')) #remove or add ss based on the type of plot you want do make
Predicted_age_abs <- Predicted_age_abs[, c(1, 5, 14, 15, 16, 17)]

Predicted_age_abs_long <- Predicted_age_abs %>%
  pivot_longer(
    cols = -c(Sample_name, Individuals),
    names_to = "ClockRep",
    values_to = "AE"
  ) 

Predicted_age_abs_long$Mixture_type <- substring(Predicted_age_abs_long$Sample_name, 1, nchar(Predicted_age_abs_long$Sample_name) - 2)
Predicted_age_abs_long$Mixture_type <- factor(Predicted_age_abs_long$Mixture_type, levels=c("AA_AB_1_1",	"AA_AB_2_1",	"AA_AB_4_1",	"AA_AB_10_1",	"AA_AB_1_2",	"AA_AB_1_4",	"AA_AB_1_10",	"AE_AF_1_1",	
                                                                                  "AE_AF_2_1",	"AE_AF_4_1",	"AE_AF_10_1",	"AE_AF_1_2",	"AE_AF_1_4",	"AE_AF_1_10"))

Predicted_age_abs_long$Replicate <- substr(Predicted_age_abs_long$Sample_name, 
                                      nchar(Predicted_age_abs_long$Sample_name), nchar(Predicted_age_abs_long$Sample_name))

Predicted_age_abs_long$Theoretical_ratio <- substr(Predicted_age_abs_long$Sample_name, 7, nchar(Predicted_age_abs_long$Sample_name)-2)

Predicted_age_abs_long$Theoretical_ratio <- gsub("_", ":", Predicted_age_abs_long$Theoretical_ratio)

Predicted_age_abs_long$Theoretical_ratio <- factor(
  Predicted_age_abs_long$Theoretical_ratio,
  levels = c("10:1", "4:1", "2:1", "1:1","1:2", "1:4", "1:10")
)

Predicted_age_abs_long$Pair <- substr(Predicted_age_abs_long$Mixture_type, 1, 5)
Predicted_age_abs_long$Pair <- gsub("_", " ", Predicted_age_abs_long$Pair)

# Apply the same recoding you used
Predicted_age_abs_long$Mixture_type <- gsub("AA_AB", "M1_F1", Predicted_age_abs_long$Mixture_type)
Predicted_age_abs_long$Mixture_type <- gsub("AC_AD", "M2_F2", Predicted_age_abs_long$Mixture_type)

Predicted_age_abs_long$Pair <- gsub("AA AB", "M1 F1", Predicted_age_abs_long$Pair)
Predicted_age_abs_long$Pair <- gsub("AE AF", "M2 F2", Predicted_age_abs_long$Pair)
Predicted_age_abs_long$Pair <- gsub(" ", "-", Predicted_age_abs_long$Pair)

colnames(Predicted_age_abs_long)[3] <- "Clock"
#add the code below if using comprison with ss
#Predicted_age_abs_long$Clock <- substr(Predicted_age_abs_long$Clock, 4, nchar(Predicted_age_abs_long$Clock))
Predicted_age_abs_long$Clock <- gsub("AE_", "", Predicted_age_abs_long$Clock)
Predicted_age_abs_long$Clock <- paste(Predicted_age_abs_long$Clock, "clock")


AE_bar <- Predicted_age_abs_long %>%
  group_by(Pair, Theoretical_ratio, Clock) %>%
  summarise(
    ymin = AE[1],
    ymax = AE[2],
    .groups = "drop"
  )


#Create the plot using both dataset with MAE and AE
MAE_plot <- ggplot() +
  # 1) Background reference lines (draw FIRST so they stay behind)
  geom_hline( #remove this part if you do not want the AE lines for single source samples
    data = hline_data,
    aes(yintercept = AE, linetype = Reference),
    color = NA,              # keep this to avoid ggplot mapping color
    linewidth = 0.55,
    alpha = 0.6,            # transparency for the background lines
    show.legend = TRUE
  ) +
  geom_hline(
    data = hline_data,
    aes(yintercept = AE, linetype = Reference),
    color = ref_colors[hline_data$Reference],
    linewidth = 0.55,
    alpha = 0.6,            # transparency for the colored lines
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  
  # Bar connecting the two absolute-error measurements
  geom_errorbar(
    data = AE_bar,
    aes(
      x = Theoretical_ratio,
      ymin = ymin,
      ymax = ymax,
      color = Pair
    ),
    alpha = 0.55,
    width = 0.18,
    linewidth = 0.7,
    position = position_dodge(width = 0.35)
  ) +
  
  # Mean point
  geom_point(
    data = Predicted_age_long,
    aes(
      x = Theoretical_ratio,
      y = MAE,
      color = Pair
    ),
    size = 2.5,
    position = position_dodge(width = 0.35)
  ) +
  # Linetype for Reference — colors added in override.aes
  scale_linetype_manual(
    values = c("M1" = "dashed", "M2" = "dashed"),
    name = "Single-Source Sample"
  ) +
  facet_wrap(
    ~ Clock,
    labeller = labeller(
      Clock = c(
        "Horvath clock" = "Horvath clock (353 CpGs)",
        "skinHorvath clock"  = "skinHorvath clock (391 CpGs)",
        "EN clock"    = "EN clock (514 CpGs)",
        "BLUP clock"     = "BLUP clock (319,607 CpGs)"
      )
    )
  ) +
  coord_cartesian(ylim = c(0, 20)) +
  theme_minimal() +
  guides( #remove this part if you do not want the AE lines for single source samples
    color = guide_legend(order = 1),
    linetype = guide_legend(
      order = 2,
      override.aes = list(
        color = unname(ref_colors),  # Inject correct color into legend
        size = 0.6
      )
    )
  ) +
  scale_color_manual(values = pair_colors, name = "Mixture") +
  
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = -45, vjust = 0, hjust = 0.1),
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 14)
  ) +
  xlab("Suspect-to-Victim Ratio") +
  ylab("MAE (years)")




ggsave(paste0(results_path,"/2_Age_prediction/Plot_MAE_offender_in_DNA_mixtures_R_13-01-2025.png"), 
       MAE_plot, width = 10.5, height = 7, dpi = 600, bg = "white")


#### 10. Calculate proportion of Type I and II probes in epienetic clocks -------------
Zhou_probe_annotation_EPIC_v2 <- read_tsv(paste0(annotation_files_path, "EPICv2.hg38.manifest.tsv.gz")) #Two different annotation files, this one contains annotation of probes with probes specification.
Zhou_probe_annotation_EPIC_v2$Probe_ID <- substr(Zhou_probe_annotation_EPIC_v2$Probe_ID, 1, 10) 
probes_type <- Zhou_probe_annotation_EPIC_v2[,c("Probe_ID", "type")]
probes_type <- probes_type[!duplicated(probes_type$Probe_ID), ]

chk <- checkClocks(betas_cg_autosomal_collps_no_NAs)

#Create proportion of Type I and II for each plot
clocks <- list(
  BLUP      = coefBLUP$CpGmarker,
  EN        = coefEN$CpGmarker,
  Horvath = coefHorvath$CpGmarker,
  SkinHorvath     = coefSkin$CpGmarker
)

prop_table <- lapply(names(clocks), function(clock_name) {
  
  cpgs <- clocks[[clock_name]]
  
  prop <- prop.table(
    table(
      probes_type$type[
        probes_type$Probe_ID %in% cpgs
      ]
    )
  )
  
  data.frame(
    Clock     = clock_name,
    ProbeType = names(prop),
    Proportion = as.numeric(prop),
    row.names = NULL
  )
}) |> bind_rows()

#Plot in barplot
type_of_probe_proportion_long <- prop_table %>%
  mutate(
    PercentLabel = paste0(round(Proportion * 100, 1), "%")
  )

type_of_probe_proportion_plot <- ggplot(
  type_of_probe_proportion_long,
  aes(
    x    = Clock,
    y    = Proportion,
    fill = ProbeType
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill"
  ) +
  geom_text(
    aes(label = PercentLabel),
    position = position_fill(vjust = 0.5),
    size = 4,
    color = "black"
  ) +
  theme_bw() +
  xlab("Epigenetic clock") +
  ylab("Relative proportion of CpG probes (%)") +
  labs(fill = "Type of probe") +
  scale_fill_manual(values = c("#7963bc", "#99B24D")) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = c(0, 0)
  ) +
  theme(
    axis.text.x = element_text(size= 11.5, angle = -45, vjust = 1, hjust = 0)
  )

ggsave(paste0(brando_path,"/Results/1_Quality_control/Type_of_probe_proportion_plot_13-01-2026.png"), 
       type_of_probe_proportion_plot, width = 7, height = 5)








