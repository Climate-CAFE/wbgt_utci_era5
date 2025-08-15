# Created by: Zach Popp
# Date Created: 04/01/2025
# Version Number: v2
# Date Modified: 04/29/2025
# Modifications:
#   Switched queue download to separate script
#
# Overview:
#     This code uses the ecmwfr R package to download ERA5 temperature measures
#     from The Copernicus Climate Data Store. More details on the ERA5 API
#     access using the ecmwfr R package are accessible at:
#         https://github.com/bluegreen-labs/ecmwfr
#     Please explore the vignettes provided and note the provided instructions
#     for how to sign up for a Copernicus account and to access the relevant
#     user and key inputs below required to query data.
#     
#     This is an a companion script to 01A which submits requests in batch.
#     This script reads in the text script output with the url of submitted
#     requests. You can monitor jobs at https://cds.climate.copernicus.eu/requests?tab=all
#     Once all jobs are complete (or in the interim when some subset are 
#     complete), this will download programmatically all files
#
# Load required packages
#
library("ecmwfr")
library("keyring") # Note: the keyring package may not be accessible in your computing
# environment. The package is used to provide a password that
# is otherwise requested directly during an interactive R 
# R Session. If the package cannot be used, this script can 
# be run in your session and the password provided directly.

################### User-Define Parameters #####################################

# Set directory. Establishing this at the start of the script is useful in case
# future adjustments are made to the path. 
#
ecmw_dir <- "RawData/ERA5_Hourly/" # Directory where rasters will be output to (add x to download then move)
home_dir <- "0_codedir/keydir/"                    # Directory where API credential are stored
trac_dir <- "0_codedir/trackdir"       # Directory where request syntax will be saved for download

# Set county This code is developed for a single US county. These could be updated to reflect
# any subdivision of your area of interest, as needed to divide processing into
# more computationally efficient steps.
#
county_in <- 13121 # Example county is Fulton County, GA

# Read in key from file. 
# NOTE: The file storing your API key should never be shared publicly
#
api_key <- scan(paste0(home_dir, "api_key.txt"), what = "", nmax = 1, quiet = TRUE)

# Set key (commented out as this is run without submitting in terminal)
# 
# keyring_pass <- scan(paste0(home_dir, "keyring.txt"), what = "", nmax = 1, quiet = TRUE)

# Set years to download
#
minyear <- 2024
maxyear <- 2024

# NOTE: To run the below, add a directory 'x' to ecmw_dir so download can occur
# The files can then be moved out to the ecmw_dir and the directory deleted.
# The jobs are submitted with an unavailable directory so that they enter the 
# queue without delay waiting for R to download each.

################## Access Filled Requests ######################################
# Set variables list
#
all_vars <- c("2m_temperature",
              "2m_dewpoint_temperature",
              "surface_net_solar_radiation",
              "surface_net_thermal_radiation",
              "surface_solar_radiation_downwards",
              "surface_thermal_radiation_downwards",
              "10m_u_component_of_wind",
              "10m_v_component_of_wind",
              "surface_pressure")

# Loop through years to download processed data
#
for (year_in in c(minyear:maxyear)) {
  
  # Loop through variables
  #
  for (var_in in all_vars) {
    
    # Track
    #
    cat("Processing ", year_in, "_", var_in, "\n")
    
    # Read in log output to extract text for complete jobs download
    #
    api_out <- readLines(paste0(trac_dir, "/console_test_", var_in , year_in,"_", county_in, "_var.txt"))
    
    # Try to extract wf_transfer lines
    #
    api_trans <- api_out[grepl("wf_trans|url|path|filename", api_out) &
                           !grepl("staging", api_out) &
                           !grepl("moved", api_out) &
                           !grepl("Delete", api_out)]
    
    # Transfer lines alone to get n
    #
    api_trans_only <- api_trans[grepl("filename", api_trans)]
    n_trans <- length(api_trans_only) - 1
    
    # Loop through string to build and run function
    #
    for (i in c(0:11)) {
      
      cat(i, "\n")
      
      # Identify start index
      #
      n_start <- 1 + i*5
      n_end <- n_start + 3
      
      # Cut to request, clean and run!
      #
      transfer_lang <- paste(unlist(api_trans[n_start:n_end]), collapse = "")
      transfer_lang <- paste(transfer_lang, ")")
      transfer_lang <- gsub(paste0("filename = '", ecmw_dir, "/x/"),
                            "filename = '", transfer_lang)
      
      e <- simpleError("test error")
      
      # Try to download
      #
      tryCatch(
        eval(parse(text = transfer_lang)),
        error = function(e) e
      ) 
    } 
  }
}


