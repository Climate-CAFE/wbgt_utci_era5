# Created by: Zach Popp
# Date Created: 05/06/2025
# Version Number: v1
# Date Modified: 
# Modifications:
#
# Overview:
#     Note: This script is essentially the same as script 1A but downloads just
#     one varibale for the ERA 25km model that has total_sky_direct_solar_radiation_at_surface
#     which is then downscaled and used in the WBGT and UTCI processing.
#
#     This code uses the ecmwfr R package to download ERA5 temperature measures
#     from The Copernicus Climate Data Store. More details on the ERA5 API
#     access using the ecmwfr R package are accessible at:
#         https://github.com/bluegreen-labs/ecmwfr
#     Please explore the vignettes provided and note the provided instructions
#     for how to sign up for a Copernicus account and to access the relevant
#     user and key inputs below required to query data.
#
#     The Copernicus Climate Data Store API restricts the number of requests
#     actively running at a given time. Requests will sit in a queue, which
#     can be monitored at https://cds.climate.copernicus.eu/requests?tab=all
#
#     Smaller requests move more quickly through the queue. In order to 
#     maximize the request efficiency this request sets up separate requests
#     for the region of interest for each month/year and variable. Based 
#     on information shared through the CDS web forum, the number of hours
#     and number of variables will affect queue time, but the size of the area
#     typically will not. Depending on your time series and area of interest,
#     you may be required to divide your request due to limits on the maximum
#     download size for a single request. 
#
#     This script uses the wf_request_batch function from the ecmwfr package
#     to submit all months of data for a single year and variable at once.
#     This will lead to all requests being added to the queue in short sequence.
#     The batches are submitted in a loop to allow for the addition of more
#     requests to the queue at a given time without waiting for requests to be
#     completed. 
#
#     This script is set up to run via a shell script so that it does not need
#     to be actively monitored in R Studio. A text file with output is written
#     to a Text_Output directory, and this text file is then used in script
#     1C to identify the links where submitted/completed jobs can be 
#     programmatically queried from. Because the ecmwfr requires a user key
#     to be provided, this script is set up to read encrypted text files where
#     API credentials are stored from a directory. It is not recommended to save
#     your API credentials in scripts given that they are sensitive information.
#
#     This example uses the Northeast United States as a test case. The shapefile
#     for the extent is downloaded using the tigris package
#
# Load required packages
#
library("ecmwfr")
library("sf")
library("dplyr")
library("keyring") # Note: the keyring package may not be accessible in your computing
# environment. The package is used to provide a password that
# is otherwise requested directly during an interactive R 
# R Session. If the package cannot be used, this script can 
# be run in your session and the password provided directly.
library("tigris")

# Check package version numbers
#
if (packageVersion("ecmwfr") < "1.5.0"   | packageVersion("sf") < "1.0.16" ) {
  cat("WARNING: packages are outdated and may result in errors.") }

################### User-Define Parameters #####################################
# Set directory. Establishing this at the start of the script is useful in case
# future adjustments are made to the path. 
#
setwd("C:/Users/zpopp/OneDrive - Boston University/Desktop/CAFE/ERA5_WBGT")

ecmw_dir <- "RawData/ERA5_25km/" # Directory where rasters will be output to
home_dir <- "0_codedir/keydir/"                    # Directory where API credential are stored
trac_dir <- "0_codedir/trackdir/"       # Directory where request syntax will be saved for download

# Set county This code is developed for a single US county. These could be updated to reflect
# any subdivision of your area of interest, as needed to divide processing into
# more computationally efficient steps.
#
county_in <- 13121 # Example county is Fulton County, GA

# NOTE: The file storing your API key should never be shared publicly!
# These will need to be set as text files - api_key will have your CDS API
# and keyring will have a password that is otherwise set and then requested by
# RStudio when attempting to unlock use of the CDS API. If not running in a
# batch script, these approaches may not be required.
#
api_key <- scan(paste0(home_dir, "api_key.txt"), what = "", nmax = 1, quiet = TRUE)

# Set key (commented out as this is run without submitting in terminal)
# 
# keyring_pass <- scan(paste0(home_dir, "keyring.txt"), what = "", nmax = 1, quiet = TRUE)

# Identify extent for download.
# LOAD Shapefile. This approach involves a US application for the Northeast US,
# downloaded using TIGRIS. If you are conducting a global analysis or have 
# an existing shapefile, the below can be replaced to just set the shapefile_cut
# input:
#     shapefile_cut <- st_read("shapefile_path")
#
shapefile_cut <- tigris::counties(year = 2020,
                                  state = substr(county_in, 1, 2))

# Subset to Northeast region
#
shapefile_cut <- shapefile_cut[shapefile_cut$GEOID == county_in, ]

# Set years to download
#
minyear <- 2024
maxyear <- 2024

################### Build Requests #############################################
# Set key (commented out as this is run without submitting in terminal)
# 
# keyring_unlock(keyring = "ecmwfr", password = keyring_pass)

# Assess bounding box. The bounding box represents the coordinates of the 
# extent of the shapefile, and will be used to specify the area we would like
# to query from Copernicus Climate Data Store. The API will allow any bounding 
# parameters; however, values that deviate from the original model grid scale
# will be interpolated onto a new grid. Therefore, it’s recommended that for 
# ERA5-Land (which is 0.1˚ resolution) the bounding coordinates be divisible by 
# 0.1 (e.g., 49.5˚N, -66.8˚E, etc.), and that coordinates for ERA5 be divisible
# by 0.25 (e.g., 49.25˚N, -66.75˚E, etc.)
#
input_bbox <- st_bbox(shapefile_cut)

# Add a small buffer around the bounding box to ensure the whole region 
# is queried, and round the parameters to a 0.1 resolution. A 0.1 resolution
# is applied because the resolution of netCDF ERA5 data is .25x.25
# https://confluence.ecmwf.int/display/CKB/ERA5%3A+What+is+the+spatial+reference
#
# To round to 0.25, use the below
#
round_to_025 <- function(x) {
  round(x / 0.25) * 0.25
}

input_bbox$xmin <- round_to_025(input_bbox$xmin[[1]]) - 0.25
input_bbox$ymin <- round_to_025(input_bbox$ymin[[1]]) - 0.25
input_bbox$xmax <- round_to_025(input_bbox$xmax[[1]]) + 0.25
input_bbox$ymax <- round_to_025(input_bbox$ymax[[1]]) + 0.25

# The set of inputs below specify the range of years to request, and sets
# the series of month state/end dates to query. There is 
# a limit on the data size that can be downloaded in a given request and smaller
# requests move more quickly through the CDS API queue, so we set up month-
# and variable-level requests in an effort to increase the efficiency
#
query_years <- c(minyear:maxyear)

# Adjust the input years based on the time period for which you want
# to query data
#
query_starts <- c("01-01", "02-01", "03-01", "04-01", "05-01", "06-01", "07-01", 
                  "08-01", "09-01", "10-01", "11-01", "12-01")
query_ends <- c("01-31", "02-29", "03-31", "04-30", "05-31", "06-30", "07-31", 
                "08-31", "09-30", "10-31", "11-30", "12-31")

# Set variables list
#
all_vars <- c("total_sky_direct_solar_radiation_at_surface")

# Initialize empty lists for each variable
#
combined_request_list <- list()

# To build the requests, we will loop through each year, each variable, and
# each month. We initialize the empty list above and then created nested lists
# within for each year and variable. Separate month-level requests are nested
# in each year-variable sublist. The wf_request_batch function used below
# will pass each request in the list to the API to be filled.
#
for (yr in query_years) {
  
  # Set up list for year
  #
  combined_request_list[[as.character(yr)]] <- list() 
  
  # We have just one variable. 
  # For each year-variable, a monthly request is setup. If a request is
  # too large, it will not be accepted by the CDS servers, so this division
  # of requests is useful and will expedite the process.
  #
  for (i in 1:12) {
    
    # Track progress
    #
    cat("Now processing quarter ", i, "\n")
    
    # Extract inputs for start and end based on list of dates at begin and 
    # end of months around quarter
    #
    query_1 <- query_starts[i]
    query_2 <- query_ends[i]
    
    # Establish query for date periods. This formats the date inputs
    # as they need to be formatted.
    #
    query_dates <- paste0(yr, "-", query_1, "/", yr, "-", query_2)
    
    # Track progress
    #
    cat("Now processing year ", yr, "\n")
    
    # The below is the formatted API request language. All of the inputs
    # specified below in proper formatting can be identified by forming a 
    # request using the Copernicus CDS point-and-click interface for data
    # requests. https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form
    # Select the variables, timing, and netcdf as the output format, and then 
    # select "Show API Request" at the bottom of the screen. 
    #
    # Note that the target is the filename that will be exported to the path
    # specified in the next part of the script. If using a loop, ensure that
    # that the unique features of each request are noted in the output! Here
    # we have each of the year, variable, and months (all our loop parameters)
    # in the filename so we won't accidentally overwrite.
    #
    request_era <- list(
      dataset_short_name = "reanalysis-era5-single-levels",
      product_type = "reanalysis",
      variable = c(all_vars),
      date = query_dates,
      time = c('00:00', '01:00', '02:00',
               '03:00', '04:00', '05:00',
               '06:00', '07:00', '08:00',
               '09:00', '10:00', '11:00',
               '12:00', '13:00', '14:00',
               '15:00', '16:00', '17:00',
               '18:00', '19:00', '20:00',
               '21:00', '22:00', '23:00'),
      data_format = "netcdf",
      download_format = "unarchived",
      area = c(input_bbox$ymax, input_bbox$xmin, input_bbox$ymin, input_bbox$xmax),
      target = paste0("era5-25km-country-",yr , "_", query_1, "_", query_2,"_", region_in, ".nc")
    )
    
    # Add request to batch list
    #
    combined_request_list[[as.character(yr)]] <- c(combined_request_list[[as.character(yr)]], list(request_era))   
    
  }
}


# If you run summary(combined_request_list), you will see that you now
# should have a series of lists of length 12 (1 request for each monht)
# and denoted by a variable name concatenated against a year
# The loop below will pass each list to the wf_request_batch function.
# 
# Within the loop we use 'sink' to ensure all warnings and output are written
# to a text file. These text files will include text that can be uses to reference
# submitted and completed jobs, so that we can programmatically download the data.
# This loop could be used to directly download, but in testing it seems that
# including the download process will be inefficient as there will be a delay in 
# adding requests to the queue while the download is pending. In order to get
# around the download, we have a path set up here that does not exist: 
#             path = paste0(ecmw_dir, "/ERA5_Hourly/x/")
# We add the /x/ so that the script will not download as the path does not
# exist. We can then add the /x/ directory to the ERA5_Hourly directly before
# running script 01C, so that the download is successful
#
# Loop through years to submit batch script
#
for (year in c(minyear:maxyear)) {
  
  # Set up text connection to track processing. This will save the warnings
  # and output to a file on your device, which can then be read in to loop
  # and download.
  #
  log_conn <- textConnection("log_output", "w", local = TRUE)
  sink(log_conn, type = "output")
  sink(log_conn, type = "message")
  
  # We run the submission within a tryCatch function as an error will be 
  # generated by the fake file path. This will only come up after all
  # requests are in
  #
  tryCatch(
    
    # Run the API request
    batch_submit <- wf_request_batch(combined_request_list[[paste0(year)]], user = "ecmwfr", 
                                     path = paste0(ecmw_dir, "/x/"), workers = 12, time_out = 100000)
    
    ,
    
    # If it fails, keep rolling
    error = function(e) e
  ) 
  
  # Stop capturing output to file
  #
  sink(NULL, type = "message")
  sink(NULL, type = "output")
  close(log_conn)
  
  # Update
  #
  cat(year, " done \n")
  
  # Save text warnings to file
  #
  writeLines(log_output, paste0(trac_dir, "console_era_25_", year, "_", county_in, ".txt"))
  
}

