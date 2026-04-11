########################################################
# Global Ecometric Analysis: Modern + Fossil + Sensitivity
# Author: Maria A. Hurtado-Materon
# Date: [Date]
# Description: Full workflow for summarizing traits,
# building ecometric models, fossil environment reconstruction,
# and model sensitivity testing (quantitative and qualitative).
########################################################

# ====== Load Required Packages ======
library(dplyr)
library(purrr)
library(sf)
library(leaflet)
library(raster)
library(ggplot2)

# ====== Load Required Functions ======
sapply(list.files("scripts/funcs", pattern = "\\.R$", full.names = TRUE), source)

# ====== Load Data ======
points <- read.csv("inputs/S_L_allpointsdata2.csv")

# continents <- raster::shapefile("inputs/continents.shp")

traits <- read.csv("inputs/RBL_global_Dec2023.csv", header = TRUE, na = ".", stringsAsFactors = FALSE)
traits$TaxonName <- gsub("_", " ", traits$TaxonName)

geometry <- sf::st_read("inputs/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
geography_carnivora <- geometry[geometry$order_ == "CARNIVORA", ]
rm(geometry)

# ====== Summarize Traits at Sampling Points ======
res <- summarize_traits_by_point(
  points_df = points,
  trait_df = traits,
  species_polygons = geography_carnivora,
  species_name_col = "sci_name",
  trait_column = "RBL",
  continent_shp = FALSE,
  lon_col = "Longitude",
  lat_col = "Latitude"
)

# ====== Determine number of bins based on  Scott’s Rule ======

optimal_bins_scott(res$points$mean_trait)
optimal_bins_scott(res$points$sd_trait)

# ====== Inspect Point Species Composition ======
inspect_point_species(
  res_list = res,
  min_species_valid = 5
)

inspect_point_species(
  res_list = res,
  min_species_valid = 5,
  env_var = "BIO12"
)


# ====== Quantitative Ecometric Model: Precipitation ======

EMprec <- ecometric_model(
  points_df = res$points,
  env_var = "BIO12",
  transform_fun = function(x) log(x + 1),
  inv_transform_fun = function(x) exp(x) - 1,
  min_species = 3
)

# ====== Quantitative Fossil Projection ======

# Load Fossil Data
fossil <- read.csv("inputs/Koobi_fora_full.csv", header = TRUE)

# Aggregate fossil community data
fossil_summary <- fossil[1:66, 1:9] %>%
  dplyr::group_by(Time) %>%
  dplyr::summarise(
    Mean = mean(RBL, na.rm = TRUE),
    SD   = sd(RBL, na.rm = TRUE),
    Lat  = dplyr::first(Lat),
    Long = dplyr::first(Long),
    Age_approx_oldest = dplyr::first(Age_approx_oldest),
    Age_approx_new    = dplyr::first(Age_approx_new),
    Paleo_envir       = dplyr::first(Paleo_envir),
    .groups = "drop"
  )

# Project fossils using precipitation model
fossil_results_prec <- reconstruct_env(
  fossildata = fossil_summary,
  model_out = EMprec,
  match_nearest = TRUE,
  fossil_lon = "Long",
  fossil_lat = "Lat",
  modern_id = "GlobalID",
  modern_lon = "Longitude",
  modern_lat = "Latitude"
)

# ====== Plot Quantitative Ecometric Space with Fossils ======
ecometric_space(
  model_out = EMprec,
  env_name = "Precipitation (mm)",
  fossil_data = fossil_results_prec
  
)

# ====== Sensitivity Analysis: Quantitative ======
sensitivity_results_quant <- sensitivity_analysis(
  points_df = res$points,
  env_var = "BIO12",
  iterations = 20,
  test_split = 0.2,
  transform_fun = function(x) log(x + 1),
  parallel = FALSE
)

# ====== Qualitative Ecometric Model: Vegetation Cover ======

EMveg <- ecometric_model_qualitative(
  points_df = res$points,
  category_col = "DOM_DESC",
  min_species = 3
)

# ====== Qualitative Fossil Projection ======

fossil_results_veg <- reconstruct_env_qualitative(
  fossildata = fossil_summary,
  model_out = EMveg,
  match_nearest = TRUE,
  fossil_lon = "Long",
  fossil_lat = "Lat",
  modern_id = "GlobalID",
  modern_lon = "Longitude",
  modern_lat = "Latitude"
)

# ====== Plot Qualitative Ecometric Space with Fossils ======
ecometric_space_qualitative(
  model_out = EMveg,
  category_labels = NULL, # or provide named vector if desired
  palette = NULL,         # or set your custom color palette
  fossil_data = fossil_results_veg
)

# ====== Sensitivity Analysis: Qualitative ======
sensitivity_results_qual <- sensitivity_analysis_qual(
  points_df = res$points,
  category_col = "DOM_DESC",
  sample_sizes = seq(100, 10000, 1000),
  iterations = 20,
  test_split = 0.2,
  parallel = TRUE
)

########################################################
# End of Full Workflow
########################################################