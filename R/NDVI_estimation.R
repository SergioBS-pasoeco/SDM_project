#NDVI execution

#!/usr/bin/env Rscript
install.packages(c("rgee","sf","terra","dplyr","lubridate"))

getwd()
setwd("C:/Users/santamaria/Documents/Spatial_analysis/GBIF/Plants")

library(rgee)
ee_install()
ee_Initialize()

source("R/ndvi_mean_sd.R") # Llamando las funciones

compute_ndvi_mean_sd_gee(
  aoi_path = "../../Vectors/CO/Casanare.shp",              # or .gpkg
  out_dir = "../SDM_project/results/ndvi",
  source = "sentinel2",                           # or "landsat"
  start_date = Sys.Date() - lubridate::years(20),
  end_date = Sys.Date(),
  scale = 10,                                     # 10 m for S2; use 30 for Landsat
  via = "drive"                                   # "drive" or "gcs"
)
