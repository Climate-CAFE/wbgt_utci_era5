# Area-Level Estimation of Wet Bulb Globe Temperature and Universal Thermal Climate Index Pipeline 

## Project Overview
R code is provided for all steps from data download to raster-level summation of hourly wet bulb globe temperature and universal thermal climate index (UTCI). R scripts are accompanied by batch scripts which can be applied to more efficiently and 'passively' download ERA5 data from the Copernicus Climate Data Store, aggregate 10+ input meteorological variables to WBGT and UTCI using the heatmetrics R package, and spatiotemporally aggregate the WBGT and UTCI to daily, county-level metrics for use in epidemiologic analyses. Comments throughout provide guidance on setting up required accounts, downloading additional data, and understanding potential pitfalls in the process.

The download of ERA5 data in this repository advances on the approach described in the [ERA5 Daily Aggregation tutorial](https://github.com/Climate-CAFE/era5-daily-heat-aggregation), and again uses the ecmwfr package to efficiently query hourly ERA5-Land and ERA5 hourly data on single pressure levels data. The heatmetrics R package is then used for the processing of the meteorological and land use variables to relevant metrics of heat stress. The heatmetrics package must be dowloaded for use of script 4 in this repository, it can be downloaded from [Figshare](https://figshare.com/articles/software/heatmetrics_R_Package/19739965?file=36838296)

### Why Use ERA5?
When selecting a dataset for analysis of the relationship between climate change and health, there are several primary considerations. They include:
- **Temporal Resolution**: Depending on your question of interest, hourly, daily, or annual summary metrics could be most appropriate. The ERA5 product we use for this analysis is hourly. There are a range of different temporal resolutions that can be downloaded from the [Climate Data Store](https://cds.climate.copernicus.eu/datasets?q=era5).
- **Time Span**: Selecting a data resource is also informed by the time span available. Some research questions require an extensive historical data record spanning far into the past, and other times the most recent data is more pertinent. ERA5-Land provides data spanning back to 1950, and data for the recent past as well (through a week prior to download).
- **Spatial Resolution**: Re-analysis data such as ERA5 use a combination of various inputs to derive estimates of temperature at a finer scale and across greater areas than individual weather stations could provide. Different data sources provide differing resolution, with some datasets sharing temperature and other metrics at 1m x 1m grids, and others at much larger scales. ERA5-Land provides data at approximately 9km x 9km grids. This larger grid size can be advantageous, in that the storage and processing of less spatially resolved data will be less intensive.
- **Spatial Extent**: One of the key advantages for ERA5 is its global spatial extent. As opposed to country-specific datasets, like the Parameter-elevation Regressions on Independent Slopes Model (PRISM) data for the US, ERA5 is available across all countries.
- **Reliability**: Whatever your research question is, you should ensure the dataset you opt to use is reliable. ERA products have been [widely used in global climate assessments](https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.3803) by organizations including WHO and the IPCC. Previous work has been conducted to compare the ERA5 outputs to available global monitoring data across different [climate regions](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=9658540), [across Europe](https://www.mdpi.com/2073-4441/14/4/543), in [East Africa](https://www.mdpi.com/2073-4433/11/9/996), as well as other regions across the globe. Further assessing whether ERA5 has been tested for reliability in your study region of interest can help inform possible limitations with respect to specific regions, measures, or other data elements.

### Why These Heat Metrics?
This repository provides code for analyzing WBGT and UTCI. These metrics . Code for processing a series of additional metrics is available from a separate [ERA5 Daily Aggregation tutorial](https://github.com/Climate-CAFE/era5-daily-heat-aggregation) which provides a tutorial on estimating three distinct metrics from ERA5 (2-m temperature, 2-m dew point temperature, skin temperature) and two derived metrics (heat index and humidex). Measures like heat index account for temperature and humidity to better approximate heat exposure than temperature alone. WBGT, and more recently UTCI, have been developed as more physiologically relevant alternatives - accounting for solar radiation, radiant temperature, and wind speed to further account for the experience of heat. The measures have been fine-tuned based on epidemiological and physiological models of human heat exposure. For more details, the background of the [Spangler et al. paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9206009/) from which this repository was developed includes several relevant references. A brief description is provided below.
- **Wet Bulb Globe Temperature (WBGT)**: WBGT was developed in the 1950s to better estimate epidemiologically relevant heat exposure thresholds for the prevention of heat-related illness in military settings, and brings together thermal, solar, and convective heat transfers from ambient temperature, humidity, solar radiation, and wind speed. In accounting for the impact of the sun and wind, it is well-positioned to account for exposure for those spending time outside like outdoor workers and athletes. Limitations of WBGT include underestimation of heat impacts when sweating is restricted and variability based on clothing and activity. 
- **Universal Thermal Climate Index (UTCI)**: UTCI has been developed to account for WBGT limitations, and is derived from human energy balance models. The goal of UTCI is broader applicability across contexts, and accounting adaptively for clothing. [Categories](https://climate-adapt.eea.europa.eu/en/metadata/indicators/thermal-comfort-indices-universal-thermal-climate-index-1979-2019) of UTCI have been developed to classify heat stress risk.


## Usage
This repository provides the building blocks for the query of any data from the Copernicus CDS. Code for spatial aggregation of ERA5-Land hourly netCDF files from their native raster format to daily measures averaged across administrative boundaries are also included, using the terra and sf packages. These spatial aggregation methods can be used to derive measures from ERA5-Land in analysis of sociodemographic disparities or epidemiologic studies of temperature and health outcomes. An additional script is provided with some example applications of R package ggplot2 to visualize temperature spatially and over time.

Comments are included within the code repository with descriptions of how to manipulate the ERA5 API language to query additional variables, time frames, or spatial extents. Please note that users will need to set up an account with the Copernicus CDS for the API query to function effectively. For more information on the ecmwfr package and details about how to set up an account and access the necessary user ID and API key inputs for the API query please see: https://github.com/bluegreen-labs/ecmwfr

### Notes on Computation

## Data Sources
- [ERA5-Land hourly data from 1950 to present](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land)
- [ERA5 hourly data on single levels from 1940 to present](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land)
- [National Land Cover Database](https://www.sciencebase.gov/catalog/item/5dfc2280e4b0b207a9fe8235)
  
## Workflow
A series of R and bash scripts are included in the code repository. Bash scripts are used to divide processing into monthly or annual increments, dramatically shortening the time of processing. The scripts included were created for use on a university cloud computing network, and may very for your use case. Learn more about high performance computing here: (https://github.com/Climate-CAFE/hpc_batch_jobs_micro_tutorial)[https://github.com/Climate-CAFE/hpc_batch_jobs_micro_tutorial]

1) 01A_Query_ERA5_Land_CDSAPI_Batch_v2.R applies the ecmwfr package to query hourly ERA5-Land data from the Copernicus CDS using the ecmwfr package. Eight different variables are queried. To increase the efficiency, batch requests are submitted by month and variable for each year to be processed. The script is accompanied by a batch script that can be used to submit to a computing cluster, to allow for processing without active engagement. An additional R script is included to download from the Copernicus queue. Script 1A adds jobs to the queue for a given month, then moves to the next without waiting to download. 
2) 02A_Query_ERA5_25km_CDSAPI_Batch_v1.R applies the ecmwfr package to query hourly ERA5 hourly data on single levels from the Copernicus CDS using the ecmwfr package. A similar process as script 1A is used, although only one variable from the ERA5 25km product is used. The same sequence of scripts is also provided for this download process.
3) 03A_UTCI_WBGT_Intermediate_Static_v1.R generates several time invariant inputs that are then used in the time varying estimation of hourly UTCI and WBGT. These include rasters with latitude and longitude, and a raster-layer with estimates of urbanicity at the ERA5-Land grid scale
4) 04A_UTCI_WBGT_Intermediate_Dynamic_v2.R generates the time varying inputs and carries through processing to the hourly WBGT and UTCI measures. The meteorological variables including solar radiation, wind, temperature, and vapor pressure are converted as needed for use and then applied to calculate WBGT and UTCI using the Spangler et al. heatmetrics R package. The script is set to run as a batch job by month.
5) 05A_ERA5_Extraction_Pts_v1.R is a pre-processing step for the aggregation of the raster data to county summaries. A grid of points from the ERA5 grid cells is created and linked with the polygons for aggregation (US counties).
6) 06A_Aggregate_ERA5_UTCI_UseFishnet_v2.R uses the hourly UTCI and WBGT outputs, and the extraction points to generate the county-level summaries. The conversion and aggregation from hourly UTCI and WBGT at the UTC time zone to daily summaries based on local time are also conducted. The output of this script is daily minimum, mean, and maximum UTCI and WBGT spatially aggregated to the county-level.
   
## Dependencies
Packages used in this repository include:
- library("ecmwfr")         for era5 data query
- library("terra")          for raster data
- library("sf")             for vector data
- library("plyr")           for data management
- library("doBy")           for aggregation of data across groups
- library("tidyverse")      for data management
- library("lwgeom")         for spatial data management
  
## References
Please note the Disclaimers_References_and_Changes.Rmd in the downloaded heatmetrics R package. The reference notes are copied below:
 When using this package, please cite the accompanying paper:
  - Spangler, K.R., S. Liang, G.A. Wellenius. "Wet-Bulb Globe Temperature, 
    Universal Thermal Climate Index, and Other Heat Metrics for US 
    Counties, 2000-2020." Scientific Data (2022).
Additionally, when using data derived with the wbgt() function or any functions on which wbgt() relies, please also cite the paper published by the original algorithm writer:
  - Liljegren, J. C., Carhart, R. A., Lawday, P., Tschopp, S. & Sharp, R. 
    Modeling the Wet Bulb Globe Temperature Using Standard Meteorological 
    Measurements. J. Occup. Environ. Hyg. 5, 645-655 (2008). 
    https://doi.org/10.1080/15459620802310770
When using data derived with the utci() function or any functions on which utci() relies, please also cite the following papers:
  - C. Brimicombe, C. Di Napoli, T. Quintino, F. Pappenberger, R. 
    Cornforth, and H.L. Cloke. "Thermofeel: A python thermal comfort 
    indices library." SoftwareX (2022). 
    https://doi.org/10.1016/j.softx.2022.101005
  - Di Napoli, C., Hogan, R.J. & Pappenberger, F. Mean radiant temperature 
    from global-scale numerical weather prediction models. Int J 
    Biometeorol 64, 1233–1245 (2020). 
    https://doi.org/10.1007/s00484-020-01900-5

ECMWFR package reference:
- Hufkens, K., R. Stauffer, & E. Campitelli. (2019). ecmwfr: Programmatic interface to the two European Centre for Medium-Range Weather Forecasts API services. Zenodo. http://doi.org/10.5281/zenodo.2647531.
