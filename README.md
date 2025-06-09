# Pipeline for Area-Level Estimation of Wet Bulb Globe Temperature and Universal Thermal Climate Index

## Project Overview
R code is provided for all steps from data download to raster-level summation of hourly wet bulb globe temperature and universal thermal climate index (UTCI). R scripts are accompanied by batch scripts which can be applied to more efficiently and 'passively' download ERA5 data from the Copernicus Climate Data Store, aggregate 10+ input meteorological variables to WBGT and UTCI using the heatmetrics R package, and spatiotemporally aggregate the WBGT and UTCI to daily, county-level metrics for use in epidemiologic analyses. Comments throughout provide guidance on setting up required accounts, downloading additional data, and understanding potential pitfalls in the process.

The download of ERA5 data in this repository advances on the approach described in the ERA5 Daily Heat Aggregation repository (https://github.com/Climate-CAFE/era5-daily-heat-aggregation), and again uses the ecmwfr package to efficiently query hourly ERA5-Land and ERA5 hourly data on single pressure levels data. The heatmetrics R package is then used for the processing of the meteorological and land use variables to relevant metrics of heat stress. The heatmetrics package must be dowloaded for use of script 4 in this repository, it can be downloaded from Figshare (https://figshare.com/articles/software/heatmetrics_R_Package/19739965?file=36838296)

## References
#   Please note the Disclaimers_References_and_Changes.Rmd in the downloaded heatmetrics R package. The reference notes are copied below:
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

#   ECMWFR package reference:
- Hufkens, K., R. Stauffer, & E. Campitelli. (2019). ecmwfr: Programmatic interface to the two European Centre for Medium-Range Weather Forecasts API services. Zenodo. http://doi.org/10.5281/zenodo.2647531.
