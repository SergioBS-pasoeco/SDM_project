### SDM key species for conservation

library(sf)
library(terra)

setwd("C:/Users/santamaria/Documents/Spatial_analysis/GBIF/Plants")

plants_AR_PA = read_sf("Plants_AR/Plants_ARG_PA.shp")

plants_AR_oPA = read_sf("Plants_AR/Plants_ARG_outPA.shp")


plants_CO_PA = read_sf("Plants_CO/Plants_COL_PA.shp")

plants_CO_oPA = read_sf("Plants_CO/Plants_COL_outPA.shp")


# aggregate by year -------------------------------------------------------
dat_all_AR = rbind(plants_AR_PA, plants_AR_oPA)
dat_all_AR$type = c(rep("Reserve",nrow(plants_AR_PA)), rep( "No Reserve", nrow(plants_AR_oPA)))
dat_all_AR$country = rep("Argentina",nrow(dat_all_AR))

library(tidyverse)

dat_all_CO = bind_rows(plants_CO_PA, plants_CO_oPA)
dat_all_CO$type = c(rep("Reserve",nrow(plants_CO_PA)), rep( "No Reserve", nrow(plants_CO_oPA)))
dat_all_CO$country = rep("Colombia",nrow(dat_all_CO))


length(unique(dat_all_AR[which(dat_all_AR$type== "Reserve"),]$species)) #Number of species in Reserves in Argentina 1697
(1697/2900)*100 # 58.5%

length(unique(dat_all_AR[which(dat_all_AR$type== "No Reserve"),]$species)) #Number of species reported outside reserves in Argentina 2677
(2677/2900)*100 # 92.3%

((2900-2677)/2900)*100 # 7.68% of plants species have been only reported inside PAs

#In total there are 2900 species reported for Buenos Aires province

length(unique(dat_all_CO[which(dat_all_CO$type== "Reserve"),]$species)) #Number of plant species in Reserves in Colombia 1010
(1010/3250)*100 # 31.0%

length(unique(dat_all_CO[which(dat_all_CO$type== "No Reserve"),]$species)) #Number of species reported outside reserves in Colombia 3094
(3094/3250)*100 # 95.2%

#In total there are 3250 species reported for Casanare
((3250-3094)/3250)*100 # 4.8 % of the species is only reported for protected areas


#Unique species Reserves
sp_reserves = unique(dat_all_AR[which(dat_all_AR$type== "Reserve"),]$species)

sp_reserves_nr = unique(dat_all_AR[which(dat_all_AR$type== "No Reserve"),]$species)

target_species = setdiff(sp_reserves, sp_reserves_nr)


# Species distribution model ----------------------------------------------

# install.packages(c("sf", "terra", "dismo", "spThin", "geodata", "dplyr", "purrr", "maxnet"))


library(sf)
library(terra)
library(dismo)
library(geodata)
library(dplyr)
library(purrr)
library(maxnet)

# Example Inputs
# sf object with GADM polygons (gadm_sf)
# data.frame with species occurrence: occurrences_df with species, longitude, latitude
# target_species <- unique(dat_all_AR$species)

run_sdm_loop <- function(gadm_sf, occurrences_df, target_species, resolution = 2.5, sdm_method = "maxent") {

  for (sp in target_species) {
    message("Modeling: ", sp)

    # 1. Filter species occurrences
    sp_occ <- occurrences_df %>%
      filter(species == sp) %>%
      distinct(decimalLon, decimalLat, .keep_all = TRUE)
    message("Number of occurrences: ", sp, "_ ",nrow(sp_occ))


    library(ggplot2)
    ggplot() +
      geom_sf(data = gadm_sf, fill = "lightblue") +
      geom_sf(data = st_as_sf(sp_occ, coords = c("decimalLon", "decimalLat"), crs = 4326), color = "red", size = 3) +
      ggtitle(sp)

    if (nrow(sp_occ) < 10) {
      warning("Skipping ", sp, ": insufficient data (<10 points)")
      next
    }

    # 2. Convert to spatial object
    sp_sf <- st_as_sf(sp_occ, coords = c("decimalLon", "decimalLat"), crs = 4326)
    crs(sp_sf)

    gadm_sf <- st_transform(gadm_sf, crs =  crs(sp_sf))
    crs(gadm_sf)

    study_area <- gadm_sf
    bbox <- st_bbox(study_area)

    # 4. Download WorldClim climatic data
    env <- geodata::worldclim_global(var = "bio", res = resolution, path = tempdir())
    env_crop <- crop(env, ext(c(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"])))
    env_masked <- mask(env_crop, vect(study_area))

    # 5. Prepare data for modeling
    occ_coords <- st_coordinates(sp_sf)
    occ_points <- as.data.frame(occ_coords)
    colnames(occ_points) <- c("lon", "lat", "ID")

    # 6. Generate background points (pseudo-absences) within 5 km buffer around occurrences
    # Transform to a projected CRS for buffering (Web Mercator)
    sp_proj <- st_transform(sp_sf, 3857)

    # Create 5 km buffer and union to a single polygon
    buffer_5km <- st_buffer(sp_proj, dist = 5000)
    buffer_union <- st_union(buffer_5km)

    # Transform buffer back to geographic CRS (WGS84)
    buffer_union_ll <- st_transform(buffer_union, 4326)

    # Convert buffer to terra SpatVector and sample background points within it
    buffer_vect <- vect(buffer_union_ll)
    set.seed(2025)
    bg_points <- spatSample(buffer_vect, size = 200, method = "random")

    bg_coords <- terra::crds(bg_points)
    colnames(bg_coords) <- c("lon", "lat")

    bg_sf <- st_as_sf(as.data.frame(bg_coords), coords = c("lon", "lat"), crs = 4326)
    occ_sf <- st_as_sf(as.data.frame(occ_coords), coords = c("X", "Y"), crs = 4326)


    ggplot() +
      geom_sf(data = gadm_sf, fill = "lightblue") +
      geom_sf(data = bg_sf, color = "black", size = 3)  +
      geom_sf(data = occ_sf, color = "red", size = 2.5) +
      ggtitle(sp)

    # 7. Run MaxEnt model
    if (sdm_method == "maxent") {
      env_vals_occ <- extract(env_masked, occ_sf)
      env_vals_bg <- extract(env_masked, bg_sf)

      pres <- rep(1, nrow(env_vals_occ))
      abs <- rep(0, nrow(env_vals_bg))

      sdm_data <- rbind(env_vals_occ, env_vals_bg)
      sdm_data$presence <- c(pres, abs)

      sdm_data <- na.omit(sdm_data)
      sdm_data = sdm_data[,-1]

      model <- maxnet(p = sdm_data$presence, data = sdm_data[, -ncol(sdm_data)], f = maxnet.formula(p = sdm_data$presence, data = sdm_data[, -ncol(sdm_data)], classes = "default"))

      # 8. Predict to raster
      prediction <- predict(env_masked, model, type = "cloglog", na.rm = TRUE)

      # 9. Save output
      out_path <- file.path("C:/Users/santamaria/Documents/Spatial_analysis/GBIF/SDM_project/models/SDM_output_plants/AR", paste0(gsub("[ /]", "_", sp), "_SDM.tif"))

      dir.create(dirname(out_path), showWarnings = FALSE)
      writeRaster(prediction, out_path, overwrite = TRUE)

      message("Saved: ", out_path)
    }
  }
}


# Running the loop --------------------------------------------------------
gadm_BA = read_sf("../../Vectors/AR/Buenos_Aires_Province.shp")

target_species_test = unique(dat_all_AR$species)
target_species_test[1:20]

# Run the loop Argentina
run_sdm_loop(
  gadm_sf = gadm_BA,
  occurrences_df = dat_all_AR,
  target_species = target_species_test[1:10],
  resolution = 2.5,  # WorldClim resolution (10 = 10 arc-min)
  sdm_method = "maxent"
)



gadm_Cas = read_sf("../../Vectors/CO/Casanare.shp")
crs(gadm_Cas)


target_species_test = unique(dat_all_CO$species)

# Run the loop Colombia
run_sdm_loop(
  gadm_sf = gadm_Cas,
  occurrences_df = dat_all_CO,
  target_species = target_species_test,
  resolution = 2.5,  # WorldClim resolution (10 = 10 arc-min)
  sdm_method = "maxent"
)
