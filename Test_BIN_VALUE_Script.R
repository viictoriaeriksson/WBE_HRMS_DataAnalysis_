
library(mzR)
library(Spectra)
library(tidyverse)
library(viridis)
library(patchwork)
library(crayon)

source("C:/Users/vicer06/OneDrive - Linköpings universitet/Documents/GitHub/WBE_HRMS_DataAnalysis_/BIN_Values_Functions.R")

data_import("C:/Users/vicer06/OneDrive - Linköpings universitet/Documents/01_Projects/03_Library_Building_DrugStD/Individual_Drug_StD_240527/Drug_StD_240527/20240527_drugmix_200ppb_woadv_pos_2_C,1_1.mzML",
            "C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\03_Library_Building_DrugStD\\Individual_Drug_StD_240527\\Drug_StD_240527\\test_arcMSCCS_script_drugmix_woadv_MZmine.csv", type="spectra_library")

extraction_of_data(sp_ims)

find_targets_in_raw(target_compounds, combine_raw_dataframe, tolerance_rt = 0.01)

old_subset_targets <- subset_targets


subset_targets |> filter(compound_name=="Temazepam")
check_dt_spectra(subset_targets, compound_name="Temazepam") 

filter_found_compounds(subset_targets)

filter_found_compounds(subset_targets, output_direction="C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\03_Library_Building_DrugStD\\Individual_Drug_StD\\Drug_StD_240527", file_name= "input_arcMS_CCS")





source("C:/Users/vicer06/OneDrive - Linköpings universitet/Documents/GitHub/WBE_HRMS_DataAnalysis_/BIN_Values_Functions.R")

data_import("C:/Users/vicer06/OneDrive - Linköpings universitet/Documents/01_Projects/03_Library_Building_DrugStD/Individual_Drug_StD_240815/240816_ESI_hdmse_newsettings/20240816_THC_200ppb_woadv_ESI_pos_Newsettings_1_B,7_1.mzML",
            "C:/Users/vicer06/OneDrive - Linköpings universitet/Documents/01_Projects/03_Library_Building_DrugStD/Individual_Drug_StD_240815/thc_aligned_all.csv", type="spectra_library")

extraction_of_data(sp_ims)

find_targets_in_raw(target_compounds, combine_raw_dataframe)


subset_targets |> filter(compound_name=="Temazepam")
check_dt_spectra(subset_targets, compound_name="Temazepam") 

filter_found_compounds(subset_targets)



######################## Drugmix Newsettings

source("C:/Users/vicer06/OneDrive - Linköpings universitet/Documents/GitHub/WBE_HRMS_DataAnalysis_/BIN_Values_Functions.R")

data_import("C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\03_Library_Building_DrugStD\\Individual_Drug_StD_240815\\240816_ESI_hdmse_newsettings\\20240816_drugmix_200ppb_woadv_ESI_pos_Newsettings_1_C,1_1.mzML",
            "C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\03_Library_Building_DrugStD\\Individual_Drug_StD_240815\\for_binscript_drugmix_newsettings.csv", type="spectra_library")

extraction_of_data(sp_ims)

find_targets_in_raw(target_compounds, combine_raw_dataframe)

check_dt_spectra(subset_targets) 

filter_found_compounds(subset_targets, output_direction="C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\03_Library_Building_DrugStD\\Individual_Drug_StD_240815", file_name= "input_arcMS_CCS_drugmixNewsettings")



