##| output: false

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# 
# BiocManager::install("tidyverse")
# BiocManager::install("factoextra")
# BiocManager::install("msdata")
# BiocManager::install("mzR")
# BiocManager::install("rhdf5")
# BiocManager::install("rpx")
# BiocManager::install("MsCoreUtils")
# BiocManager::install("QFeatures")
# BiocManager::install("Spectra")
# BiocManager::install("ProtGenerics")
# BiocManager::install("PSMatch")
# BiocManager::install("pheatmap")
# BiocManager::install("limma")
# BiocManager::install("MSnID")
# BiocManager::install("RforMassSpectrometry/SpectraVis")
# 
# BiocManager::install("synapter")

library(mzR)
library(Spectra)
library(tidyverse)
library(viridis)
library(patchwork)
library(crayon)

############################################################################
#' Imports raw IMS-HRMS data (.mzML) and target compounds (.csv) to be found in the raw data
#'
#' This functions takes a mzML-file and/or a csv-file (output from MZmine) and import 
#' them to the global environment. 
#' For the outtput from MZmine, depending on this method used for finding compounds
#' the type should be put to spectra_library or precursor_search.

#' @param direction_raw Direction to mzML file, should be IMS-HRMS data
#' @param direction_target Direction to csv-file, contining target compounds in output format from MZmine
#' @param type should etiher be spectra_library or precursor_search
#' @return Two data frames, one for the raw data (sp_ims) and one for target compounds (target_compounds), found in global environment.
#' @examples
#'   data_import("path/raw_data", "path/target_data", "spectra_library")
#'   data_import("path/target_data", "spectra_library")
#'   data_import("path/raw_data")

data_import <- function(direction_raw = NULL, direction_target = NULL, type = NULL) {
  if (!is.null(direction_raw)) {
    raw_spectra <- Spectra(direction_raw)
    assign("sp_ims", raw_spectra, envir=globalenv())
  }
  
  if (!is.null(direction_target)) {
    if (missing(type) || !(type %in% c("spectra_library", "precursor_search"))) {
      stop(
        paste0(
          italic("type"),
          " needs ",
          bold('"spectra_library"'),
          " or ",
          bold('"precursor_search"')
        )
      )
    }
    if (type =="spectra_library") {
      #import the output from MZmine, when spectra library search have been used
      target_import <- read_csv(direction_target) |> 
        select("spectral_db_matches:compound_name", "mz", "rt", "ion_mobility") |> 
        rename(compound_name="spectral_db_matches:compound_name")
    } else if (type =="precursor_search") {
      #import the output from MZmine, when precursor search have been used
      target_import <- read_csv(direction_target) |> 
        select("compound_db_identity:compound_name", "mz", "rt","ion_mobility") |>
        rename(compound_name="compound_db_identity:compound_name")
    }
    assign("target_compounds", target_import, envir=globalenv())
  }
}

###############################################################################
#' Extraction of MS1 Level parameters from a Spectra objects
#'
#' This function takes an S4 object and filters it for MS1 and from these 
#' extract useful parameter (mz, retetion time, intensity, drift time) as wel as 
#' creating a extract parameters called bin. The parameters are then combined to a
#' data frame exported to the global environment, for later use.

#' @param spectra_object An S4 object containing mass spectrometry data.
#' @return A dataframe (combine_raw_dataframe) with extracted MS1 parameters, found in global environment.

extraction_of_data <- function(spectra_object) {
  if (!isS4(spectra_object)) {
    stop(
      paste0(
        italic("spectra_object"),
        " needs to be a ",
        bold('"S4"'),
        " object."
      )
    )
  }
  # Filter for MS1 
  sp_ims_ms1 <- filterMsLevel(spectra_object, msLevel. = 1)
  
  # Extract mz, rtime, intensity, and scan index. Futhermore, bin values are created from scan index
  mz_ms1 <- mz(sp_ims_ms1)
  rtime_ms1 <- rtime(sp_ims_ms1)
  intensity_ms1 <- intensity(sp_ims_ms1)
  scanindex_ms1 <- scanIndex(sp_ims_ms1)
  bin_values_ms1 <- rep(1:200, length.out = length(scanindex_ms1))
  dtime_ms1 <- sp_ims_ms1[["ionMobilityDriftTime"]]
  
  #Unlist mz and intensity values so it becomes a vector instead of a list with lists
  mz_ms1_unlisted <- unlist(mz_ms1)
  intensity_ms1_unlisted <- unlist(intensity_ms1)
  
  # Function that repeates the values in some of the objects to be same length as mz values
  repeat_ms1 <- function(A, B){
    unlist(mapply(function(a, b) rep(b, length(a)), A, B))
  }
  
  #Use function to extend rtime, scandindex, and bin values so they are repated to the length of mz
  rtime_ms1_repated <- repeat_ms1(mz_ms1, rtime_ms1)
  scanindex_ms1_repated <- repeat_ms1(mz_ms1, scanindex_ms1)
  bin_values_ms1_repated <- repeat_ms1(mz_ms1, bin_values_ms1)
  dtime_ms1_repated <- repeat_ms1(mz_ms1, dtime_ms1)
  
  
  #Combine the new objects into a dataframe to make it easier to wokr with later
  combine_dataframe <- data.frame(scanindex_ms1 = scanindex_ms1_repated,
                                  bin_values_ms1 = bin_values_ms1_repated,
                                  mz_ms1 = mz_ms1_unlisted,
                                  rtime_ms1 = rtime_ms1_repated/60,
                                  intensity_ms1 = intensity_ms1_unlisted,
                                  dtime_ms1 = dtime_ms1_repated
  )
  #Remove the unlisted and repated objects to same memory
  #rm(rtime_ms1_repated, scanindex_ms1_repated, bin_values_ms1_repated, mz_ms1_unlisted, intensity_ms1_unlisted, dtime_ms1_repated)
  assign("combine_raw_dataframe", combine_dataframe, envir=globalenv())
}

###############################################################################

#' Finds the target compounds in the raw IMS-HRMS MS1 data
#'
#' This functions takes a data frame containing target compunds (should be the
#'  output from MZmine) to match with a data frame containing raw IMS-HRMS MS1
#'  parameters (with at least mz and rt values)

#' @param target_file Direction to mzML file, should be IMS-HRMS data
#' @param raw_MS1parameters_dataframe Direction to csv-file, contining target compounds in output format from MZmine
#' @param tolerance_mz tolerance to be used to match precursor m/z, default 0.003
#' @param tolerance_rt tolerance to be used to match retention time, default 0.001
#' @return Data frames with the found target compounds in the raw data (subset_targets), found in global environment.

find_targets_in_raw <- function(target_file=NULL, raw_MS1parameters_dataframe=NULL, tolerance_mz=0.003, tolerance_rt=0.001) {
  if (is.null(target_file) || is.null(raw_MS1parameters_dataframe)) {
    stop(
      paste0("Need the inputs of ",
             bold("target_file"),
             " and ",
             bold("raw_MS1parameters_dataframe"),
             " \n to be able to find the compounds in the raw data.")
    )
  }
  
  # Adds mz and rt range to the target compounds
  target_file$mz_min <- target_file$mz-tolerance_mz 
  target_file$mz_max <- target_file$mz+tolerance_mz
  target_file$rt_min <- target_file$rt-tolerance_rt
  target_file$rt_max <- target_file$rt+tolerance_rt
  # Finds the compounds from the target compounds in the raw data by using the join function
  subset_targets_in_raw <- target_file |> 
    inner_join(raw_MS1parameters_dataframe, join_by(overlaps("mz_min","mz_max","mz_ms1", "mz_ms1"),overlaps("rt_min","rt_max","rtime_ms1", "rtime_ms1")))
  
  assign("subset_targets", subset_targets_in_raw, envir=globalenv())
}
################################################################################

#' Plots drift time spectra for found compounds
#'
#' This function uses a data frame with found target compounds in raw IMS-HRMS MS1
#'  data and plots the drift time spectra for each compound and retention time as
#'  bin values vs intensities

#' @param found_targets An data frame with found target compounds
#' @return Plots with bin values vs intensities for each compounds retention time

check_dt_spectra <- function(found_targets, compound_name=NULL) {
  
  # Check if compound_name is provided
  if (!is.null(compound_name)) {
    found_targets <- found_targets |> filter(compound_name == !!compound_name)
  }
  
  # Plot of intensities versus bins (related to drift times) for each compound to check the drift time peak forms
  p <- found_targets |> #filter(compound_name == "Ecgonine methyl ester") |> 
    ggplot(aes(x=bin_values_ms1,y=intensity_ms1))+
    geom_line(linewidth=0.4)+
    geom_point(aes(colour=compound_name), size=1)+
    
    #facet_wrap(vars(compound_name), scales="free", nrow=3)+ #separates the plots for each compound and scales separeately for each plot
    facet_wrap(vars(compound_name, round(rtime_ms1, 4)), scales="free_y", nrow=3)+ #separates the plots for each compound and rtime
    
    theme_minimal()+
    theme(
      #panel.border = element_blank(), axis.line = element_line(),
      panel.background = element_rect(fill = "white"),  # Set white background
      panel.grid = element_blank(),  # Remove grid lines
      axis.text = element_text(size = 9),  # Increase font size of axis text
      axis.title.y = element_text(size = 9,face = "bold"),  # Make y-axis label bold and italic
      axis.text.y = element_text(color = "black", size = 9),  # Make x-axis text bold
      axis.title.x = element_text(size = 10,face = "bold"),
      axis.text.x = element_text(color = "black", size = 9),  # Make x-axis text bold
      axis.ticks.x = element_line(color = "black"),  # Add ticks to x-axis in black color
      axis.ticks.y = element_line(color = "black"),  # Add ticks to y-axis in black color
      #aspect.ratio = 1,  # Set aspect ratio
      #plot.margin = margin(0.5, 0.5, 0.4, 0.5, "cm"),  # Set plot margins
      legend.text=element_text(size=9),
      legend.position = "none",
      legend.background = element_rect(fill="white",
                                       size=0.5, linetype="solid", 
                                       colour ="black"),
      legend.title = element_blank())
  return(p)
}

###############################################################################

#' Filter the found target compounds data frame
#'
#' This function uses a data frame with found target compounds in raw IMS-HRMS MS1
#'  data and filters it intensity and deviation between targets compounds drift 
#'  times and raw IMS-HRMS MS1 drift times. Create an output .csv file to be used 
#'  with arcMS packaged to obtain CCS values from UNIFI for target compounds.

#' @param found_targets An data frame with found target compounds
#' @param intensity_cutoff Cutoff intensity to filter compounds, with default 300
#' @param dt_tolerance Tolerance between targets and raw datas drift times, default is 1e-8, if on wants no deviation
#' @param output_direction Optional, Direction for the output .csv file with filtered data
#' @param file_name Optional, Name for the output .csv file to be used in arcMS
#' @return Data frame with the filtered compounds (input_arcMS_CCS), found in global environment and optional .csv file with compounds.
#' @examples
#'   filter_found_compounds(data_frame, "path/output_file", "name_output_file")

filter_found_compounds <- function(found_targets, intensity_cutoff=300, dt_tolerance=1e-8, output_direction=NULL, file_name=NULL) {
  #Filter by intensity and to match the DT from MZmine and Raw data
  tolerance <-dt_tolerance
  filtered_compounds <- found_targets |> group_by(compound_name, rtime_ms1) |> 
    filter(max(intensity_ms1) > intensity_cutoff) |>
    filter(abs(ion_mobility - dtime_ms1) < tolerance)
  
  assign("input_arcMS_CCS", filtered_compounds, envir=globalenv())
  
  if(is.null(output_direction) && !is.null(file_name)){
    stop(
      paste0("Need the ",
             bold("output_direction"),
             ".")
    )
  } 
  if (is.null(file_name) && !is.null(output_direction)){
    stop(
      paste0("Need the ",
             bold("file_name"),
             ".")
    )
  }
  
  if (!is.null(output_direction) && !is.null(file_name)) {
    # Construct the full file path
    output_path <- file.path(output_direction, paste0(file_name, ".csv"))
    #save to a .csv-file for input with arcMS package
    write.csv(filtered_compounds, file = output_path, row.names = FALSE)
  } 
}

###############################################################################

