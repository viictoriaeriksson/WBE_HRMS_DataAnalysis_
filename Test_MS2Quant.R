devtools::install_github("kruvelab/MS2Quant",
                         ref="main",
                         INSTALL_opts="--no-multiarch")


library(MS2Quant)

########### test that the package was installed correctly ######################
path_dataframe_calibrants_suspects <- system.file("example_data", "quantification_example.csv", package = "MS2Quant")
path_eluent_file <- system.file("example_data", "eluent.csv", package = "MS2Quant")
path_suspects_sirius_project_folder <- system.file("example_data", "SIRIUS_results", package = "MS2Quant")

MS2Quant_quantification_results <- MS2Quant_quantify(path_dataframe_calibrants_suspects,
                                                     path_eluent_file ,
                                                     organic_modifier = "MeCN",
                                                     pH_aq = 2.7,
                                                     path_suspects_sirius_project_folder)

# Separate calibration plots for each calibrant
MS2Quant_quantification_results$calibrants_separate_plots

# Calibration plot between experimental logRF and logIE
MS2Quant_quantification_results$logIE_logRF_calibration_plot

# Summary of the linear model between logRF and logIE
MS2Quant_quantification_results$calibration_linear_model_summary

# All suspect concentrations
MS2Quant_quantification_results$suspects_concentrations

# Date when the quantification was done
MS2Quant_quantification_results$date




############## Drugmix ###################################################
path_dataframe_calibrants_suspects <- read_delim("C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\03_Library_Building_DrugStD\\Individual_Drug_StD_240815\\quantification.csv", 
                                                 delim = ";", escape_double = FALSE, trim_ws = TRUE)
path_eluent_file <- read_csv("C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\03_Library_Building_DrugStD\\Individual_Drug_StD_240815\\eluent.csv")
path_suspects_sirius_project_folder <- "C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\03_Library_Building_DrugStD\\Individual_Drug_StD_240815\\sirius_results"

MS2Quant_quantification_results <- MS2Quant_quantify(path_dataframe_calibrants_suspects,
                                                     path_eluent_file ,
                                                     organic_modifier = "MeOH",
                                                     pH_aq = 2.87,
                                                     path_suspects_sirius_project_folder)

# Separate calibration plots for each calibrant
MS2Quant_quantification_results$calibrants_separate_plots

# Calibration plot between experimental logRF and logIE
MS2Quant_quantification_results$logIE_logRF_calibration_plot

# Summary of the linear model between logRF and logIE
MS2Quant_quantification_results$calibration_linear_model_summary

# All suspect concentrations
MS2Quant_quantification_results$suspects_concentrations

# Date when the quantification was done
MS2Quant_quantification_results$date



##################### NMCT StD B calibration curve ###############################
path_dataframe_calibrants_suspects <- read_delim("C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\02_NM-CT-stdAB_Collaboration_CCS_MEASURMENT_2024\\NMCTAB_2024_03_21_MSconvert_MSDIAL\\ms2quant_stdB_pos_allconc - Copy.csv", 
                                                 delim = ";", escape_double = FALSE, trim_ws = TRUE)

path_eluent_file <- read_csv("C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\02_NM-CT-stdAB_Collaboration_CCS_MEASURMENT_2024\\NMCTAB_2024_03_21_MSconvert_MSDIAL\\eluent.csv")
path_suspects_sirius_project_folder <- "C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\02_NM-CT-stdAB_Collaboration_CCS_MEASURMENT_2024\\NMCTAB_2024_03_21_MSconvert_MSDIAL\\sirius_results"


MS2Quant_quantification_results <- MS2Quant_quantify(path_dataframe_calibrants_suspects,
                                                     path_eluent_file ,
                                                     organic_modifier = "MeOH",
                                                     pH_aq = 2.87,
                                                     path_suspects_sirius_project_folder)

# Separate calibration plots for each calibrant
MS2Quant_quantification_results$calibrants_separate_plots

# Calibration plot between experimental logRF and logIE
MS2Quant_quantification_results$logIE_logRF_calibration_plot

# Summary of the linear model between logRF and logIE
MS2Quant_quantification_results$calibration_linear_model_summary

# All suspect concentrations
MS2Quant_quantification_results$suspects_concentrations

# Date when the quantification was done
MS2Quant_quantification_results$date

##################### NMCT StD A calibration curve ###############################
path_dataframe_calibrants_suspects <- read_delim("C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\02_NM-CT-stdAB_Collaboration_CCS_MEASURMENT_2024\\NMCTAB_2024_03_21_MSconvert_MSDIAL\\ms2quant_stdA_pos_allconc - Copy.csv", 
                                                 delim = ";", escape_double = FALSE, trim_ws = TRUE)

path_eluent_file <- read_csv("C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\02_NM-CT-stdAB_Collaboration_CCS_MEASURMENT_2024\\NMCTAB_2024_03_21_MSconvert_MSDIAL\\eluent.csv")

MS2Quant_quantification_results <- MS2Quant_quantify(path_dataframe_calibrants_suspects,
                                                     path_eluent_file ,
                                                     organic_modifier = "MeOH",
                                                     pH_aq = 2.87)

# Separate calibration plots for each calibrant
MS2Quant_quantification_results$calibrants_separate_plots

# Calibration plot between experimental logRF and logIE
MS2Quant_quantification_results$logIE_logRF_calibration_plot

# Summary of the linear model between logRF and logIE
MS2Quant_quantification_results$calibration_linear_model_summary

# All suspect concentrations
MS2Quant_quantification_results$suspects_concentrations

# Date when the quantification was done
MS2Quant_quantification_results$date


##################### NMCT StD A and B  calibration curve ###############################
path_dataframe_calibrants_suspects <- read_delim("C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\02_NM-CT-stdAB_Collaboration_CCS_MEASURMENT_2024\\NMCTAB_2024_03_21_MSconvert_MSDIAL\\ms2quant_stdA_stdB_pos_allconc - Copy.csv", 
                                                 delim = ";", escape_double = FALSE, trim_ws = TRUE)

path_eluent_file <- read_csv("C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\02_NM-CT-stdAB_Collaboration_CCS_MEASURMENT_2024\\NMCTAB_2024_03_21_MSconvert_MSDIAL\\eluent.csv")

MS2Quant_quantification_results <- MS2Quant_quantify(path_dataframe_calibrants_suspects,
                                                     path_eluent_file ,
                                                     organic_modifier = "MeOH",
                                                     pH_aq = 2.87)

# Separate calibration plots for each calibrant
MS2Quant_quantification_results$calibrants_separate_plots

# Calibration plot between experimental logRF and logIE
MS2Quant_quantification_results$logIE_logRF_calibration_plot

# Summary of the linear model between logRF and logIE
MS2Quant_quantification_results$calibration_linear_model_summary

# All suspect concentrations
MS2Quant_quantification_results$suspects_concentrations

# Date when the quantification was done
MS2Quant_quantification_results$date