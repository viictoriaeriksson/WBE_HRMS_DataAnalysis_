
source("C:/Users/vicer06/OneDrive - Linköpings universitet/Documents/02_Data Processing/01_WBE_HRMS_DataAnalysis_/BIN_Values_Functions.R")

data_import("C:/Users/vicer06/OneDrive - Linköpings universitet/Documents/01_Projects/03_Library_Building_DrugStD/Individual_Drug_StD_240527/Drug_StD_240527/20240527_drugmix_200ppb_woadv_pos_2_C,1_1.mzML",
            "C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\03_Library_Building_DrugStD\\Individual_Drug_StD_240527\\Drug_StD_240527\\test_arcMSCCS_script_drugmix_woadv_MZmine.csv", type="spectra_library")

extraction_of_data(sp_ims)

find_targets_in_raw(target_compounds, combine_raw_dataframe)


subset_targets |> filter(compound_name=="Temazepam")
check_dt_spectra(subset_targets, compound_name="Temazepam") 

filter_found_compounds(subset_targets)

filter_found_compounds(subset_targets, output_direction="C:\\Users\\vicer06\\OneDrive - Linköpings universitet\\Documents\\01_Projects\\03_Library_Building_DrugStD\\Individual_Drug_StD\\Drug_StD_240527", file_name= "input_arcMS_CCS")