#!/usr/bin/env Rscript

# NDVI mean/sd generator using Google Earth Engine (Sentinel-2 or Landsat)
#
# Requirements:
# - R packages: rgee, sf, terra, dplyr, lubridate
# - Google Earth Engine account and authentication (once):
#     library(rgee); ee_Initialize()
# - A polygon AOI file (e.g., shapefile or GeoPackage) defining the area of interest
#
# Outputs:
# - GeoTIFF with two bands: ndvi_mean, ndvi_sd
#
# Usage example (interactive):
#   source('R/ndvi_mean_sd.R')
#   compute_ndvi_mean_sd_gee(
#     aoi_path = 'path/to/your/aoi.shp',
#     out_dir = 'data/processed/ndvi',
#     source = 'sentinel2',  # or 'landsat'
#     start_date = Sys.Date() - lubridate::years(20),
#     end_date = Sys.Date(),
#     scale = 10              # 10 m for S2, 30 m for Landsat
#   )

suppressPackageStartupMessages({
  library(rgee)
  library(sf)
  library(terra)
  library(dplyr)
  library(lubridate)
})

initialize_gee <- function() {
  # Initialize GEE (will prompt auth on first run)
  rgee::ee_Initialize()
}

read_aoi <- function(aoi_path) {
  aoi_sf <- sf::read_sf(aoi_path)
  if (!inherits(aoi_sf, "sf")) stop("AOI could not be read as sf")
  # Dissolve multipart into single geometry
  aoi_sf <- aoi_sf %>% st_make_valid() %>% summarise(geometry = st_union(geometry))
  aoi_sf <- st_transform(aoi_sf, 4326)
  aoi_ee <- rgee::sf_as_ee(aoi_sf)
  list(sf = aoi_sf, ee = aoi_ee$geometry())
}

# ---------------------- Sentinel-2 helpers ----------------------

mask_s2_sr <- function(image) {
  # Mask clouds, cirrus, cloud shadows, snow using SCL band
  scl <- image$select("SCL")
  valid <- scl$eq(4)$Or(scl$eq(5))$Or(scl$eq(6))$Or(scl$eq(7))
  image$updateMask(valid)
}

add_s2_ndvi <- function(image) {
  ndvi <- image$expression(
    "(NIR - RED) / (NIR + RED)",
    list(
      NIR = image$select("B8"),
      RED = image$select("B4")
    )
  )$rename("NDVI")$toFloat()
  image$addBands(ndvi)
}

get_s2_ndvi_ic <- function(aoi, start_date, end_date) {
  ee$ImageCollection("COPERNICUS/S2_SR_HARMONIZED")$
    filterBounds(aoi)$
    filterDate(as.character(start_date), as.character(end_date))$
    map(mask_s2_sr)$
    map(add_s2_ndvi)$
    select("NDVI")
}

# ---------------------- Landsat helpers ----------------------

mask_landsat_qa <- function(image) {
  # QA_PIXEL bits: 3 = cloud shadow, 4 = snow, 5 = cloud (C2 L2)
  qa <- image$select("QA_PIXEL")
  cloud_shadow <- qa$bitwiseAnd(bitwShiftL(1, 3))$eq(0)
  snow <- qa$bitwiseAnd(bitwShiftL(1, 4))$eq(0)
  cloud <- qa$bitwiseAnd(bitwShiftL(1, 5))$eq(0)
  image$updateMask(cloud_shadow$And(snow)$And(cloud))
}

scale_landsat_sr <- function(image) {
  # Apply scale and offset to SR bands
  scaled <- image$select("SR_B.*")$multiply(0.0000275)$add(-0.2)
  image$addBands(scaled, NULL, TRUE)
}

add_ndvi_l57 <- function(image) {
  # Landsat 5/7: RED = SR_B3, NIR = SR_B4
  ndvi <- image$expression(
    "(NIR - RED) / (NIR + RED)",
    list(
      NIR = image$select("SR_B4"),
      RED = image$select("SR_B3")
    )
  )$rename("NDVI")$toFloat()
  image$addBands(ndvi)
}

add_ndvi_l89 <- function(image) {
  # Landsat 8/9: RED = SR_B4, NIR = SR_B5
  ndvi <- image$expression(
    "(NIR - RED) / (NIR + RED)",
    list(
      NIR = image$select("SR_B5"),
      RED = image$select("SR_B4")
    )
  )$rename("NDVI")$toFloat()
  image$addBands(ndvi)
}

get_landsat_ndvi_ic <- function(aoi, start_date, end_date) {
  l5 <- ee$ImageCollection("LANDSAT/LT05/C02/T1_L2")$
    filterBounds(aoi)$
    filterDate(as.character(start_date), as.character(end_date))$
    map(mask_landsat_qa)$
    map(scale_landsat_sr)$
    map(add_ndvi_l57)$
    select("NDVI")

  l7 <- ee$ImageCollection("LANDSAT/LE07/C02/T1_L2")$
    filterBounds(aoi)$
    filterDate(as.character(start_date), as.character(end_date))$
    map(mask_landsat_qa)$
    map(scale_landsat_sr)$
    map(add_ndvi_l57)$
    select("NDVI")

  l8 <- ee$ImageCollection("LANDSAT/LC08/C02/T1_L2")$
    filterBounds(aoi)$
    filterDate(as.character(start_date), as.character(end_date))$
    map(mask_landsat_qa)$
    map(scale_landsat_sr)$
    map(add_ndvi_l89)$
    select("NDVI")

  l9 <- ee$ImageCollection("LANDSAT/LC09/C02/T1_L2")$
    filterBounds(aoi)$
    filterDate(as.character(start_date), as.character(end_date))$
    map(mask_landsat_qa)$
    map(scale_landsat_sr)$
    map(add_ndvi_l89)$
    select("NDVI")

  l5$merge(l7)$merge(l8)$merge(l9)
}

# ---------------------- Main compute function ----------------------

compute_ndvi_mean_sd_gee <- function(
  aoi_path,
  out_dir,
  source = c("sentinel2", "landsat"),
  start_date = Sys.Date() - years(20),
  end_date = Sys.Date(),
  scale = NULL,
  file_prefix = NULL,
  via = c("drive", "gcs")
) {
  source <- match.arg(source)
  via <- match.arg(via)

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  initialize_gee()
  aoi <- read_aoi(aoi_path)

  if (is.null(file_prefix)) {
    file_prefix <- paste0(source, "_", format(start_date, "%Y%m%d"), "_", format(end_date, "%Y%m%d"))
  }

  if (source == "sentinel2") {
    ic <- get_s2_ndvi_ic(aoi$ee, start_date, end_date)
    if (is.null(scale)) scale <- 10
  } else {
    ic <- get_landsat_ndvi_ic(aoi$ee, start_date, end_date)
    if (is.null(scale)) scale <- 30
  }

  # Reduce to per-pixel mean and stdDev
  ndvi_mean <- ic$mean()$rename("ndvi_mean")
  ndvi_sd <- ic$reduce(ee$Reducer$stdDev())$rename("ndvi_sd")
  out_img <- ndvi_mean$addBands(ndvi_sd)$clip(aoi$ee)

  out_path <- file.path(out_dir, paste0(file_prefix, "_NDVI_mean_sd.tif"))

  message("Exporting to ", out_path)
  rgee::ee_as_raster(
    image = out_img,
    region = aoi$ee,
    dsn = out_path,
    scale = scale,
    via = via,
    container = "rgee"
  )

  invisible(out_path)
}

# ---------------------- CLI runner (optional) ----------------------

run_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) return(invisible(NULL))

  # Simple CLI: aoi_path out_dir source start_date end_date scale
  aoi_path <- args[[1]]
  out_dir <- args[[2]]
  source <- ifelse(length(args) >= 3, args[[3]], "sentinel2")
  start_date <- ifelse(length(args) >= 4, args[[4]], as.character(Sys.Date() - years(20)))
  end_date <- ifelse(length(args) >= 5, args[[5]], as.character(Sys.Date()))
  scale <- ifelse(length(args) >= 6, as.numeric(args[[6]]), NA_real_)

  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  if (is.na(scale)) scale <- NULL

  compute_ndvi_mean_sd_gee(
    aoi_path = aoi_path,
    out_dir = out_dir,
    source = source,
    start_date = start_date,
    end_date = end_date,
    scale = scale
  )
}

if (identical(environment(), globalenv())) {
  try(run_cli(), silent = TRUE)
}


