# Created by: Zach Popp
# Date Created: 04/01/2025
# Version Number: v2
# Date Modified: 04/29/2025
# Modifications:
#   Switched queue download to separate script
#
# Overview:
#     This code uses downloaded data from script 1A/1C ERA5-Land data, and 
#     National Land Cover Database to get time-invariant (static) inputs for the 
#     wet bulb globe temperature estimation in script 04A. As these are time
#     invariant, they are processed separately since the next script is set up
#     to run separately (as a batch script for each year)
#
#     The measures include ERA5 pixel level average land cover, which is then
#     dichotomoized to urban/nonurban (urban = >33% of pixel is urban land use),
#     and rasterized latitude and longitude for the spatial extent of the process
#
# Install or update packages as needed
#
library("terra")  # For raster data
library("sf")     # For vector data
library("plyr")   # For data management
library("doBy")   # For aggregation of data across groups
library("tidyverse") # For data management
library("lwgeom")
library("lubridate")

# Check package version numbers
#
if (packageVersion("terra") < "1.7.78"   | packageVersion("sf") < "1.0.7" | 
    packageVersion("plyr")  < "1.8.7"    |
    packageVersion("doBy")  < "4.6.19"   | packageVersion("lwgeom") < "0.2.8") {
  cat("WARNING: packages are outdated and may result in errors.") }

########################## User-Defined Parameters #############################
# Set up directories to read in and output data
#
era_dir <- "RawData/ERA5_Hourly/"     # ERA5-Land Rasters
era_interdir <- "InterDir/"           # Intermediate where static rasters output
nlcd_dir <- "RawData/NLCD/"           # Where NLCD data is stored
water_dir <- "RawData/JRC_PermWater/" # Where JRC Permanent water is stored

# Set region. This code is developed for US counties with the region representing
# the four US regions in the contiguous US. These could be updated to reflect
# any subdivision of your area of interest, as needed to divide processing into
# more computationally efficient steps.
#
county_in <- 13121 # Example county is Fulton County, GA

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

############# Read in Inputs ###################################################
# Set up path for ERA5 files. The time invariant is only based on the resolution
# and lat/lon which are consistent across all times/variables. So we list
# all files but only read in 1
#
era_files <- list.files(paste0(era_dir, "/"), pattern=paste0('.*.nc'), full.names = F)

# Read in for subset year month specified from batch input
#
era_files <- era_files[1]

# Read in file
#
era5_rast <- rast(paste0(era_dir, "/", era_files))

# We only need one raster
#
era5_rast <- subset(era5_rast, 1)

######################## Urbanicity #######################################
# For WBGT, we need urbanicity. This comes from NLCD. The NLCD
# can be downloaded from https://www.sciencebase.gov/catalog/item/5dfc2280e4b0b207a9fe8235
#
# The file for use is: nlcd_2001_landcover_2011_edition_2014_10_10.zip
# 
# Following the approach of Spangler et al., the 2011 NLCD was used for this
# processing
#
nlcd_files <- list.files(path = nlcd_dir, full.names = TRUE) 
nlcd <- rast(nlcd_files[grepl("img", nlcd_files)])

# Crop NLCD to extent of ERA5 data
#
shape_proj <- project(vect(shapefile), crs(nlcd))

# Crop to larger extent so it will overlap with ERA
# The pivot requirements will vary based on the area of your extent. The pivot
# of the NLCD through projection does required a fairly substanital buffer
# to ensure that the overlap is reasonable. We crop NLCD before projecting
# to the CRS of the ERA data because the NLCD data is massive (30m res)
#
extent_plus <- ext(shape_proj)
extent_plus[1] <- extent_plus[1] - 100000
extent_plus[2] <- extent_plus[2] + 100000
extent_plus[3] <- extent_plus[3] - 100000
extent_plus[4] <- extent_plus[4] + 100000

nlcd_crop <- crop(nlcd, extent_plus)

# Apply restriction based on Keith paper
#
nlcd_crop_urb <- terra::ifel(nlcd_crop == "Developed, Low Intensity" |
                               nlcd_crop == "Developed, Medium Intensity" |
                               nlcd_crop == "Developed, High Intensity", 1, 0)

# Summarize across ERA grid points
#
nlcd_era <- project(nlcd_crop_urb, crs(era5_rast))

# Recrop to smaller area
#
nlcd_crop <- crop(nlcd_era, ext(era5_rast), snap = "near")

# Aggregate to align
#
nlcd_agg <- aggregate(nlcd_crop, fact = 322.8439, fun = "mean")

# Resample to align
#
nlcd_res <- resample(nlcd_agg, era5_rast)
compareGeom(nlcd_res, era5_rast)

# Make it yes no urban
#
nlcd_urban <- ifel(nlcd_res > 0.333, 1, 0)

# Save study area urbanicity (this can be applied for all layers)
#
writeRaster(nlcd_urban, paste0(era_interdir, "nlcd_urbanicity_", county_in, ".tif"), overwrite = TRUE)


######################## Permanent Water #######################################
# For county aggregation, we want to exclude grids that are >50% comprised of 
# water (ie: lakes). This comes the Joint Research Centre Data Catalogue and
# can be downloaded from https://data.jrc.ec.europa.eu/dataset/jrc-floods-floodmapgl_permwb-tif#dataaccess
# 
# Following the approach of Spangler et al., the 2011 NLCD was used for this
# processing
#
jrc_perm <- rast(paste0(water_dir, "floodMapGL_permWB.tif"))

# Crop NLCD to extent of ERA5 data
#
shape_proj <- project(vect(shapefile), crs(jrc_perm))

# Crop to larger extent so it will overlap with ERA
# The pivot requirements will vary based on the area of your extent. The pivot
# of the NLCD through projection does required a fairly substanital buffer
# to ensure that the overlap is reasonable. We crop NLCD before projecting
# to the CRS of the ERA data because the NLCD data is massive (30m res)
#
extent_plus <- ext(shape_proj)
extent_plus[1] <- extent_plus[1] - 1
extent_plus[2] <- extent_plus[2] + 1
extent_plus[3] <- extent_plus[3] - 1
extent_plus[4] <- extent_plus[4] + 1

water_crop <- crop(jrc_perm, extent_plus)

# Project
#
water_era <- project(water_crop, crs(era5_rast))

# Recrop to smaller area
#
water_era <- crop(water_era, ext(era5_rast), snap = "near")

# Aggregate to align
#
water_era_agg <- disagg(water_era, fact = 12)

# Resample to align
#
water_res <- resample(water_era_agg, era5_rast)
compareGeom(water_res, era5_rast)

# Make it yes no urban
#
water_sum <- ifel(water_res > 0.5, 1, 0)

# Save study area urbanicity (this can be applied for all layers)
#
writeRaster(water_sum, paste0(era_interdir, "jrc_perm_era5.tif"), overwrite = TRUE)

######################## Lat/Lon ###################################################

# Create raster with latitude and longitude that can be repeated for all layers
# and used in WBGT function
#
# Get raster template (grid cells where lon and lat will be estimated)
r_i <- era5_rast

# Isolate rasters where there is some data
#
vals <- values(r_i)
valid_idx <- which(!is.na(vals))

# Extract coordinates and add whether ID is valid (not NA) to coordinate data
#
coords <- crds(r_i, df = TRUE)
coords$i <- valid_idx

# Build empty layers for lat and lon
#
result_lat <- rep(NA_real_, ncell(r_i))
result_lon <- rep(NA_real_, ncell(r_i))

# Make lat and lon x and y from coords for valid IDs
#
result_lat[valid_idx] = coords[coords$i == valid_idx, "y"]
result_lon[valid_idx] = coords[coords$i == valid_idx, "x"]

# Set into lists that can be collapsed to a raster layer
#
output_lat <- setValues(r_i, result_lat)
output_lon <- setValues(r_i, result_lon)

# Save study area lat/lon rasters (these can be used across time points)
#
writeRaster(output_lat, paste0(era_interdir, "era5_lat_rast.tif"), overwrite = TRUE)
writeRaster(output_lon, paste0(era_interdir, "era5_lon_rast.tif"), overwrite = TRUE)
