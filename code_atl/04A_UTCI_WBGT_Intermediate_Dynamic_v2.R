# Created by: Zach Popp
# Date Created: 04/01/2025
# Version Number: v2
# Date Modified: 04/29/2025
# Modifications:
#   Switched queue download to separate script
#
# Overview:
#     This code uses downloaded data from script 1A/1C ERA5-Land data, 
#     downloaded FDIR from ERA5 (25km) from script 2A/2C and 
#     static inputs (urbanicity, lat/lon rasters) to estimate two metrics of
#     heat exposure - wet bulb globe temperature (WBGT) and universal
#     thermal climate index (UTCI). 
#
#     These estimations involve a series of intermediary steps to prepare ERA5
#     raw data to be in the correct units, time scale, and to prepare additional
#     measures based on ERA5 inputs. The different steps are described in brief
#     below:
#
#         Resample ERA5 FDIR
#         Solar ERA5 Disaggregation and Unit Conversion
#         Estimation of cosine of zenith of solar angle
#         Calculation of mean radiant temperature (and differential from air temp)
#         Estimation of vapor pressure and relative humidity
#         Estimation of wind speed
#         Application of intermediates to UTCI and WBGT calculation
#
#   These steps are all conducted using the heatmetrics package available on 
#   Figshare (https://figshare.com/articles/software/heatmetrics_R_Package/19739965?file=36838296)
#
#   Please note the Disclaimers_References_and_Changes.Rmd in the downloaded
#   package. The reference notes are copied below:
#
#   When using this package, please cite the accompanying paper:
#       - Spangler, K.R., S. Liang, G.A. Wellenius. "Wet-Bulb Globe Temperature, 
#         Universal Thermal Climate Index, and Other Heat Metrics for US 
#         Counties, 2000-2020." Scientific Data (2022).
#   Additionally, when using data derived with the wbgt() function or any 
#   functions on which wbgt() relies, please also cite the paper published by 
#   the original algorithm writer:
#       - Liljegren, J. C., Carhart, R. A., Lawday, P., Tschopp, S. & Sharp, R. 
#         Modeling the Wet Bulb Globe Temperature Using Standard Meteorological 
#         Measurements. J. Occup. Environ. Hyg. 5, 645-655 (2008). 
#         https://doi.org/10.1080/15459620802310770
#   When using data derived with the utci() function or any functions on which 
#   utci() relies, please also cite the following papers:
#       - C. Brimicombe, C. Di Napoli, T. Quintino, F. Pappenberger, R. 
#         Cornforth, and H.L. Cloke. "Thermofeel: A python thermal comfort 
#         indices library." SoftwareX (2022). 
#         https://doi.org/10.1016/j.softx.2022.101005
#       - Di Napoli, C., Hogan, R.J. & Pappenberger, F. Mean radiant temperature 
#         from global-scale numerical weather prediction models. Int J 
#         Biometeorol 64, 1233–1245 (2020). 
#         https://doi.org/10.1007/s00484-020-01900-5


##### Install / Load Packages / Prep Directories ###############################

library("terra")  # For raster data
library("sf")     # For vector data
library("plyr")   # For data management
library("doBy")   # For aggregation of data across groups
library("tidyverse") # For data management
library("lwgeom")
library("weathermetrics")
library("lubridate")

# Check package version numbers
#
if (packageVersion("terra") < "1.7.78"   | packageVersion("sf") < "1.0.7" | 
    packageVersion("plyr")  < "1.8.7"    |
    packageVersion("doBy")  < "4.6.19"   | packageVersion("lwgeom") < "0.2.8") {
  cat("WARNING: packages are outdated and may result in errors.") }

# Set region. This code is developed for US counties with the region representing
# the four US regions in the contiguous US. These could be updated to reflect
# any subdivision of your area of interest, as needed to divide processing into
# more computationally efficient steps.
#
county_in <- 13121 # Example county is Fulton County, GA

########################## User-Defined Parameters #############################
# Set up directories to read in and output data. Having subfolders by region
# can help to maintain organization 
# 
era_dir <- paste0("RawData/ERA5_Hourly/")        # Should be where all ERA5-Land were downloaded
era_25_dir <- paste0("RawData/ERA5_25km/")# Should be where ERA5 (FDIR) was downloaded
era_interdir <- paste0("InterDir/")   # Where static urbanicity/lat/lon were output
                                        # AND where ERA5 UTCI and WBGT outputs will be saved
# Set years to process 
#
minyear <- 2024
maxyear <- 2024

# LOAD Shapefile. This approach involves a US application for the Northeast US,
# downloaded using TIGRIS. If you are conducting a global analysis or have 
# an existing shapefile, the below can be replaced to just set the shapefile_cut
# input:
#     shapefile_cut <- st_read("shapefile_path")
#
shapefile <- tigris::counties(year = 2020,
                              state = substr(county_in, 1, 2))

# Subset to Northeast region
#
shapefile <- shapefile[shapefile$GEOID == county_in, ]


##### Introduce Batch Arguments (Replace With Manual Index As Needed) ##########
# Get batch input
#
#args <- commandArgs(trailingOnly = TRUE)  # These are arguments passed by the bash script
#b <- as.numeric(args[1])        
# Because we are running for one county, we will use a loop to process different
# months instead of using the batch
#

for (b in c(4:13)) {
  
  
  # Source functions from heatmetrics package:
  #
  # The package can be downloaded from figshare. Update the paths below to where
  # you place the unzipped heatmetrics directory.
  #
  # Reference: K.R. Spangler, S. Liang, and G.A. Wellenius. "Wet-Bulb Globe 
  # Temperature, Universal Thermal Climate Index, and Other Heat Metrics for US 
  # Counties, 2000-2020." Scientific Data (2022). doi: 10.1038/s41597-022-01405-3
  #
  # With the below we can source all functions
  #
  all_functions <- list.files("C:/Users/zpopp/OneDrive - Boston University/Desktop/CAFE/ERA5_WBGT/heatmetrics/R/",
                              full.names = TRUE)
  
  lapply(all_functions, source)
  
  # Track
  # 
  cat(b, "/n")
  
  # This is the index for year/month combinations. We run one at a time because 
  # the processing  is cumbersome with so many variable inputs and intermediates. 
  # The argument "b" will index against a list running from 2000_01 to 2020_12 to 
  # indicate what month/year of data (for all variables) should be processed.
  # If you are processing a small area, or short time series, modifying the script
  # to run without required indexing may be useful. This can be done by removing 
  # these lines and separately noting the rasters which should be read in for
  # processing
  
  # Set series of year_months to index files and bring in one month at a time 
  # for processing
  #
  years <- c(minyear:maxyear)
  months <- sprintf("%02d", 1:12)  # format months as two digits
  year_months <- as.vector(outer(years, months, paste, sep = "_"))
  year_months <- sort(year_months)
  
  # Add one month in advance as the time transformation for disaggregation of 
  # cumulative daily measures and time zone conversion may require added dates
  #
  year_months <- c(paste0(minyear-1, "_12"), year_months)
  
  # Assess input from vector list based on batch input
  #
  year_pre <- year_months[b - 1]
  year_month <- year_months[b]
  year_post <- year_months[b + 1]
  
  # If year_pre is empty,
  
  # Build query subset
  #
  years <- paste0(year_pre, "|", year_month, "|", year_post)
  
  
  ############# Derivation of ERA5 Intermediates ######################
  
  # INITIAL SETUP 
  # Set up ERA5 files
  #
  era_files <- list.files(paste0(era_dir, "/"), pattern=paste0('.*.nc'), full.names = F)
  
  # Read in for subset year month specified from batch input
  #
  era_files <- era_files[grepl(years, era_files)]
  
  # Read in files for the year month specified from batch
  #
  era5_rast <- rast(paste0(era_dir, "/", era_files))
  
  # For at least one date (3/31/2001), the downloaded ERA5 data includes
  # duplicate layers. This causes errors in the code, so to address this
  # we will restrict to unique layers. We can do this based on the layer names
  #
  # Get layer names
  #
  layer_names <- names(era5_rast)
  
  # Identify unique names
  #
  unique_names <- !duplicated(layer_names)
  
  # Subset the raster
  #
  era5_rast <- era5_rast[[unique_names]]
  
  # Add time dimension (NA when downloaded in current API configuration)
  #
  terra::time(era5_rast) <- as_datetime(as.numeric(sub(".*=", "", names(era5_rast))))
  terra::time(era5_rast) <- as.POSIXct(terra::time(era5_rast), format = "%Y-%m-%d %Z", tz = "UTC")
  
  # Build parameters of hours before and after input year
  #
  rpre <- max(time(era5_rast)[grepl(paste0(substr(year_pre, 1, 4), "-", substr(year_pre, 6, 7)), time(era5_rast))])
  rpost <- min(time(era5_rast)[grepl(paste0(substr(year_post, 1, 4), "-", substr(year_post, 6, 7)), time(era5_rast))])
  
  # Subset to maintain last hour of previous month and first of the next month
  #
  era5_rast <- subset(era5_rast, time(era5_rast) >= rpre)
  era5_rast <- subset(era5_rast, time(era5_rast) <= rpost)
  
  # Confirm all variables are available:
  #
  table(substr(names(era5_rast), 1, 4))
  table(substr(time(era5_rast), 6, 10))
  # The data from ERA5 if the previous call was effective should include the names
  # below:
  #
  # d2m - 2-meter dew point temperature
  # sp - surface pressure
  # ssr - surface_net_solar_radiation
  # ssrd - surface_solar_radiation_downwards
  # str - surface_net_thermal_radiation
  # strd - surface_thermal_radiation_downwards
  # t2m - 2-meter air temperature
  # u10 - 10m_u_component_of_wind
  # v10 - 10m_v_component_of_wind
  
  ######################## Resample 25km FDIR ####################################
  
  # FDIR is required for the UTCI and WBGT algorithms. This is only available in 
  # the 25km ERA5 product (as total_sky_direct_solar_radiation_at_surface).
  # To use with other inputs at 9km resolution, we resample the raw data to the 
  # ERA5-Land (9km) resolution
  #
  # List ERA 25km rasters
  #
  era_25_files <- list.files(paste0(era_25_dir, "/"), pattern=paste0('.*.nc'), full.names = F)
  
  # Read in for subset year month specified from batch input
  #
  era_25_files <- era_25_files[grepl(year_month, era_25_files)]
  
  # Read in files for the year month specified from batch
  #
  fdir_25km <- rast(paste0(era_25_dir, "/", era_25_files))
  
  # Update time based on file names
  #
  terra::time(fdir_25km) <- as_datetime(as.numeric(substr(names(fdir_25km), 17, 29)))
  terra::time(fdir_25km) <- as.POSIXct(terra::time(fdir_25km), format = "%Y-%m-%d %Z", tz = "UTC")
  
  # Subset times
  #
  fdir_25km <- subset(fdir_25km, terra::time(fdir_25km) > rpre & terra::time(fdir_25km) < rpost)
  
  # Resample to era5 temp using nearest neighbor. This is in line with Spangler 
  # et al. 
  #
  # "We also obtained total sky direct solar radiation at surface from ERA528,29, 
  # the reanalysis from which ERA5-Land is derived. We interpolated this from the 
  # 0.25-degree ERA5 grid to the 0.1-degree ERA5-Land grid using nearest-neighbor 
  # interpolation, following the approach of Yan et al.30."
  #
  # Yan, Y., Xu, Y. & Yue, S. A high-spatial-resolution dataset of human thermal 
  # stress indices over South and East Asia. Scientific data 8, 229, 
  # https://doi.org/10.1038/s41597-021-01010-w (2021).
  #
  fdir_25km_proj <- terra::project(fdir_25km, crs(era5_rast))
  fdir_25km_resamp <- resample(fdir_25km_proj, era5_rast, method = "near")
  
  # We want FDIR in Watts, so we convert from the hourly to seconds by dividing
  # by 3600 (60 minutes x 60 seconds)
  #
  fdir_25km_resamp_W <- fdir_25km_resamp / 3600
  
  ######################## Solar Estimates ##############################
  
  # Several solar data elements are needed to get mean radian temperature (Tmrt), 
  # which is required for UTCI. Solar estimates are also taken as inputs for 
  # WBGT. These functions require transformed solar estimates 
  #
  # The solar data is cumulative over time, so we need to refit to get the hourly
  # estimates. We will loop through all days included in the data and subtract the
  # previous hour's estimate from the 'current' hours. 
  # For the Tmrt function, we also want to convert from Joules/m2
  # to Watts. That is also done here after converting to hourly measures, by
  # dividing the input by 3600 (60 mins x 60 sec)
  #
  # Get solar rasters
  #
  era5_rast_ssrd_in <- subset(era5_rast, grepl("ssrd", names(era5_rast))) 
  era5_rast_ssr_in <- subset(era5_rast, grepl("ssr_", names(era5_rast))) 
  era5_rast_strd_in <- subset(era5_rast, grepl("strd", names(era5_rast))) 
  era5_rast_str_in <- subset(era5_rast, grepl("str_", names(era5_rast)))
  
  # Set up rasters for output (initialized separately due to required iteration
  # through previous measures)
  #
  era5_rast_ssrd_W <- era5_rast_ssrd_in
  era5_rast_ssr_W <- era5_rast_ssr_in
  era5_rast_strd_W <- era5_rast_strd_in
  era5_rast_str_W <- era5_rast_str_in
  
  # Extract raster times into list
  #
  times_rast <- time(era5_rast_ssrd_in)
  
  # Isolate unique times
  #
  dates_rast <- unique(substr(times_rast, 1, 10))[2:(length(unique(substr(times_rast, 1, 10)))-1) ]
  hours_rast <- c(1:24)
  hours_rast <- formatC(hours_rast, digits = 1, flag = "0")
  
  # For each hour in each day, loop through and recalculate solar estimate
  #
  for (date_in in c(dates_rast)) {
    
    # Track progress
    cat(date_in, "\n")
    
    # Extract index aligning with dates
    index_date <- which(grepl(date_in, time(era5_rast_ssrd_in)))
    index_date <- c(min(index_date) - 1, index_date[1:24])
    
    # Subset rasters to date specified for assessing hourly solar values
    era5_rast_ssrd_date <- era5_rast_ssrd_in[[index_date]]
    era5_rast_ssr_date <- era5_rast_ssr_in[[index_date]]
    era5_rast_strd_date <- era5_rast_strd_in[[index_date]]
    era5_rast_st_date <- era5_rast_str_in[[index_date]]
    
    # Initialize rasters to update
    era5_rast_ssrd_date2 <- era5_rast_ssrd_date
    era5_rast_ssr_date2 <- era5_rast_ssr_date
    era5_rast_strd_date2 <- era5_rast_strd_date
    era5_rast_st_date2 <- era5_rast_st_date
    
    # Loop through times and revise to be subtract previous hour
    #
    for (time_index in c(2:25)) {
      
      # For 1:00AM, no subtraction is needed. Just divide and moveon
      #
      if (time_index == 3) {
        era5_rast_ssrd_out <- subset(era5_rast_ssrd_date, time_index)/3600
        era5_rast_ssr_out <- subset(era5_rast_ssr_date, time_index)/3600
        era5_rast_strd_out <- subset(era5_rast_strd_date, time_index)/3600
        era5_rast_st_out <- subset(era5_rast_st_date, time_index)/3600
        
        # Add to raster, We save results to a new layer because the above code relies 
        # on evaluating against the past solar estimates in their native format
        #
        era5_rast_ssrd_date2[[time_index]] <- era5_rast_ssrd_out
        era5_rast_ssr_date2[[time_index]] <- era5_rast_ssr_out
        era5_rast_strd_date2[[time_index]] <- era5_rast_strd_out
        era5_rast_st_date2[[time_index]] <- era5_rast_st_out
        
        next
      }
      
      # Calculate hourly solar radiation by subtracting the past hour from each
      # estimate and then convert to Watts by dividing by 3600 (60 mins x 60 sec)
      #
      era5_rast_ssrd_out <- (subset(era5_rast_ssrd_date, time_index) - subset(era5_rast_ssrd_date, time_index - 1))/3600
      era5_rast_ssr_out <- (subset(era5_rast_ssr_date, time_index) - subset(era5_rast_ssr_date, time_index - 1))/3600
      era5_rast_strd_out <- (subset(era5_rast_strd_date, time_index) - subset(era5_rast_strd_date, time_index - 1))/3600
      era5_rast_st_out <- (subset(era5_rast_st_date, time_index) - subset(era5_rast_st_date, time_index - 1))/3600
      
      # Add to raster, We save results to a new layer because the above code relies 
      # on evaluating against the past solar estimates in their native format
      #
      era5_rast_ssrd_date2[[time_index]] <- era5_rast_ssrd_out
      era5_rast_ssr_date2[[time_index]] <- era5_rast_ssr_out
      era5_rast_strd_date2[[time_index]] <- era5_rast_strd_out
      era5_rast_st_date2[[time_index]] <- era5_rast_st_out
      
    }
    
    # Write to updated raster output. We have to overwrite on a new raster because
    #
    #
    era5_rast_ssrd_W[[index_date[2:25]]] <- era5_rast_ssrd_date2[[2:25]]
    era5_rast_ssr_W[[index_date[2:25]]] <- era5_rast_ssr_date2[[2:25]]
    era5_rast_strd_W[[index_date[2:25]]] <- era5_rast_strd_date2[[2:25]]
    era5_rast_str_W[[index_date[2:25]]] <- era5_rast_st_date2[[2:25]]
    
  }
  
  # Remove first and last layer
  #
  era5_rast_ssrd_W <- subset(era5_rast_ssrd_W, time(era5_rast_ssrd_W) > rpre & time(era5_rast_ssrd_W) < rpost)
  era5_rast_ssr_W <- subset(era5_rast_ssr_W, time(era5_rast_ssr_W) > rpre & time(era5_rast_ssr_W) < rpost)
  era5_rast_strd_W <- subset(era5_rast_strd_W, time(era5_rast_strd_W) > rpre & time(era5_rast_strd_W) < rpost)
  era5_rast_str_W <- subset(era5_rast_str_W, time(era5_rast_str_W) > rpre & time(era5_rast_str_W) < rpost)
  
  # We only need added bounds, redefine time bounds for rest of processing
  #
  era5_rast <- subset(era5_rast, time(era5_rast) > rpre)
  era5_rast <- subset(era5_rast, time(era5_rast) < rpost)
  
  # Subset to a single varibale. We will use this as a 'template' raster since
  # it will have the time dimension, number of layers, and extent that is repeated
  # across inputs
  #
  era5_rast_t2m <- subset(era5_rast, grepl("t2m", names(era5_rast)))
  
  ######################## Estimate CZA ############################################
  
  # To estimate UTCI, we need to get mean radiant temperature (Tmrt). To get to
  # Tmrt, we need to get CZA - the cosine of the zenith angle of the sun
  #     (see calc_cza_int function for details)
  # Inputs required: lat, lon, time (inc. hours)
  
  # We will loop through layers to grab the year/month/day/hour/lat/lon
  # and get the CZA for future use
  #
  output_list <- vector("list", nlyr(era5_rast_t2m))
  
  for (i in 1:nlyr(era5_rast_t2m)) {
    
    # Set up initial raster for time and number of layers
    r_i <- era5_rast_t2m[[i]]
    
    # Get time extraction
    time_string <- time(era5_rast_t2m)[i]
    y <- substr(time_string, 1, 4)
    m <- substr(time_string, 6, 7)
    d <- substr(time_string, 9, 10)
    h <- as.numeric(substr(time_string, 12, 13))
    
    # For 00:00, h is being dropped. Add based on timestring lenght
    if (nchar(as.character(time_string)) < 15) {
      h <- 0
    }
    
    # Get values of raster and identify non-NA areas
    vals <- values(r_i)
    valid_idx <- which(!is.na(vals))
    coords <- crds(r_i, df = TRUE)
    coords$i <- valid_idx
    
    # Along all cells set as NA, then full with CZA estimation outputs
    result_vals <- rep(NA_real_, ncell(r_i))
    result_vals[valid_idx] <- calc_cza_int(
      lat = coords[coords$i == valid_idx, "y"],
      lon = coords[coords$i == valid_idx, "x"],
      y = y, m = m, d = d, h = h
    )
    
    # List all CZA outputs for all layers
    output_list[[i]] <- setValues(r_i, result_vals)
  }
  
  # Stack all layers into a final output raster
  #
  output_cza_int <- rast(output_list)
  
  ######################## Mean radiant temp ##########################
  
  # "To calculate the mean radiant temperature following the approach
  # used by Di Napoli (2020), https://doi.org/10.1007/s00484-020-01900-5"
  # Tmrt function
  #
  # NOTE: ERROR WAS COMING THROUGH FROM Tmrt.R function - to correct I swapped
  # the ifelse function for terra::if_el. It seems to have made the issue resolve.
  # Mean radiant temperature takes in six different input rasters
  #
  Tmrt_out <- Tmrt(era5_rast_ssrd_W, era5_rast_ssr_W, fdir_25km_resamp_W,
                   era5_rast_strd_W, era5_rast_str_W, output_cza_int)
  
  # Set time
  #
  time(Tmrt_out) <- time(era5_rast_t2m)
  
  # We want to estimate the difference between radiant temperature and air temp
  # for UTCI. This is Tmrt minus ambient temp (t2m). 
  #
  DTmrt <- Tmrt_out - era5_rast_t2m
  
  ######################## Vapor pressures ################################
  
  # We need vapor pressure inputs for UTCI!
  # e = vapor pressure 
  # Source: Lawrence, M. G. The relationship between relative humidity and the 
  # dewpoint temperature in moist air - A simple conversion and applications. 
  # B. Am. Meteorol. Soc. 86, 225–233, https://doi.org/10.1175/Bams-86-2-225 (2005).
  #
  # 2-meter dew point temperature is needed in Celsius - convert from Kelvin
  #
  era5_rast_d2m <- subset(era5_rast, grepl("d2m_", names(era5_rast)))
  era5_rast_d2m_C <- era5_rast_d2m - 273.15
  
  # Estimate vapor pressure
  #
  e <- app(era5_rast_d2m_C, function(x) {
    
    610.94 * exp( (17.625*x) / (243.04+x) )
    
  })
  
  # 2-meter air temperature is also needed in Celsius - convert from Kelvin
  #
  era5_rast_t2m_C <- era5_rast_t2m - 273.15
  
  # es = saturation pressure 
  # Source: Lawrence, M. G. The relationship between relative humidity and the dewpoint temperature in moist air - A simple conversion and applications. B. Am. Meteorol. Soc. 86, 225–233, https://doi.org/10.1175/Bams-86-2-225 (2005).
  #
  e_s <- app(era5_rast_t2m_C, function(x) {
    
    610.94 * exp( (17.625*x) / (243.04+x) )
    
  })
  
  # Convert both to kilo pascals
  #
  e <- e / 1000
  e_s <- e_s / 1000
  
  # Estimate relative humidity
  #
  relhum <- 100 * (e / e_s)
  
  ######################## Wind Speed ###########################################
  
  # Get u speed and v speed
  #
  era5_rast_u10 <- subset(era5_rast, grepl("u10", names(era5_rast)))
  era5_rast_v10 <- subset(era5_rast, grepl("v10", names(era5_rast)))
  
  # To get UTCI, we also need wind speed!
  # Following function from Spangler et al., the below will compute ws_init
  # as the square root of the sum of the u and v components squared
  #
  u2 <- era5_rast_u10 ^ 2
  v2 <- era5_rast_v10 ^ 2
  
  # Get sum
  #
  u2v2 <- u2+ v2
  
  # Get windspeed
  #
  ws_init <- sqrt(u2v2)
  
  ######################## UTCI ##################################################
  
  # UTCI REFERENCE AND SOURCE
  # C. Brimicombe, C. Di Napoli, T. Quintino, F. Pappenberger, R. Cornforth, and 
  # H.L. Cloke. "Thermofeel: A python thermal comfort indices library." SoftwareX 
  # (2022). https://doi.org/10.1016/j.softx.2022.101005
  #
  # Di Napoli, C., Hogan, R.J. & Pappenberger, F. Mean radiant temperature from 
  # global-scale numerical weather prediction models. Int J Biometeorol 64, 
  # 1233–1245 (2020). https://doi.org/10.1007/s00484-020-01900-5
  #
  # To be able to run the utci function across the 5 raster layer inputs, we will
  # "Vectorize" the function as below. This vectorization allows for use of the 
  # lapp terra function which will apply the function across all cells of a layer
  # with the five different inputs
  #
  # Wrap function to be vectorized over rasters
  #
  utci_vec <- Vectorize(utci, vectorize.args = c("Tair", "e", "es", "ws", "D_Tmrt"))
  
  # Initialize raster (we will overwrite this with the UTCI output)
  #
  utci <- era5_rast_t2m_C
  
  # Loop through to apply for all layers
  #
  for (i in c(1:nlyr(era5_rast_t2m_C))) {
    
    # Track progress 
    cat(i, "\n")
    
    # Apply function to layer
    utci_out <- lapp(c(era5_rast_t2m_C[[i]], e[[i]], e_s[[i]], ws_init[[i]], DTmrt[[i]]), fun = utci_vec)
    
    # Add WBGT layer
    names(utci_out) <- "utci"
    
    # Save layer
    utci[[i]] <- utci_out
  }
  
  # Re-apply time
  #
  time(utci) <- time(era5_rast_t2m_C)
  
  # Save universal thermal climate index
  #
  writeRaster(utci, paste0(era_interdir, "utci_final_", year_month, "_", county_in, ".tif"), overwrite = TRUE)
  
  
  ######################## WBGT ################################################
  
  # WET BULB GLOBE REFERENCE AND SOURCE
  # Liljegren, J. C., Carhart, R. A., Lawday, P., Tschopp, S. & Sharp, R. Modeling 
  # the Wet Bulb Globe Temperature Using Standard Meteorological Measurements. J. 
  # Occup. Environ. Hyg. 5, 645-655 (2008). 
  # https://doi.org/10.1080/15459620802310770
  #
  # Get surface pressure
  #
  era5_rast_pressure <- subset(era5_rast, grepl("sp_", names(era5_rast)))
  era5_rast_pressure <- era5_rast_pressure / 100
  
  # Read in urbanicity, lat, and lon rasters (time invariant)
  #
  nlcd_urban <- rast(paste0(era_interdir, "nlcd_urbanicity_", county_in, ".tif"))
  era5_lat <- rast(paste0(era_interdir, "era5_lat_rast.tif"))
  era5_lon <- rast(paste0(era_interdir, "era5_lon_rast.tif"))
  
  # Set empty raster (as above, create template to be over-written)
  #
  output_wbgt <- era5_rast_t2m
  
  # As with UTCI, we will "Vectorize" the function to be able to use lapp and
  # to use raster and numeric inputs
  #
  wbgt_vec <- Vectorize(wbgt, vectorize.args = c("year", "month", "dday", "lat", "lon",
                                                 "solar", "cza", "fdir", "pres",
                                                 "Tair", "relhum", "speed", "zspeed",
                                                 "dT", "urban"))
  
  
  # Estimate fraction of SSR that is direct
  #
  frac_dir <- ifel(era5_rast_ssrd_W > 0, fdir_25km_resamp_W/era5_rast_ssrd_W, 0)
  
  # Run loop to account for time and assess WBGT
  #
  for (i in 1:nlyr(era5_rast_t2m)) {
    
    # Track N layers
    cat(i,"\n")
    
    # Get input layer
    r_i <- era5_rast_t2m[[i]]
    
    # Get times out to pass to function
    time_string <- time(era5_rast_t2m)[i]
    y <- substr(time_string, 1, 4)
    m <- substr(time_string, 6, 7)
    
    # Get decimal day
    d <- as.numeric(substr(time_string, 9, 10))
    h <- as.numeric(substr(time_string, 12, 13))
    
    # For 00:00, h is being dropped. Add based on timestring lenghth
    if (nchar(as.character(time_string)) < 15) {
      h <- 0
    }
    
    # Calculate "decimal day"
    h_f <- h / 24
    dec_d <- d + h_f
    
    # Get valid indices
    vals <- values(r_i)
    valid_idx <- which(!is.na(vals))
    
    # Stack all
    #
    rstack <- c(era5_lat, 
                era5_lon, 
                era5_rast_ssrd_W[[i]], 
                output_cza_int[[i]], 
                frac_dir[[i]], 
                era5_rast_pressure[[i]], 
                era5_rast_t2m_C[[i]],
                relhum[[i]], 
                ws_init[[i]], 
                nlcd_urban)
    
    # Use wrapper function to have numeric and raster inputs
    wbgt_wrapper <- function(x) {
      wbgt(
        year = as.numeric(y),    # your fixed year
        month = as.numeric(m),      # your fixed month
        dday = dec_d,      # your fixed day
        lat = x[1],
        lon = x[2],
        solar = x[3],
        cza = x[4],
        fdir = x[5],
        pres = x[6],
        Tair = x[7],
        relhum = x[8],
        speed = x[9],
        zspeed = 10,
        dT = -0.052,
        urban = x[10]
      )
    }
    
    # Apply wrapper across raster values for the layer
    result <- app(rstack, wbgt_wrapper)
    
    # Add WBGT layer
    names(result) <- "wbgt"
    
    # Add result to base
    output_wbgt[[i]] <- result
    
  }
  
  # Add time back (lost in translation)
  #
  time(output_wbgt) <- time(era5_rast_t2m_C)
  
  # Save wet bulb globe
  #
  writeRaster(output_wbgt, paste0(era_interdir, "wbgt_final_", year_month, "_", county_in, ".tif"), overwrite = TRUE)
  
  
  # Remove raster
  #
  rm(utci, era5_rast, era5_rast_d2m,
     era5_rast_d2m_C, era5_rast_pressure, 
     era5_rast_ssr_date, era5_rast_ssr_date2,
     era5_rast_v10, era5_rast_u10,
     era5_rast_t2m_C, era5_rast_t2m,
     output_wbgt)
  gc()
}

