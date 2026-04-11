library(sf)
library(dplyr)
library(commecometrics)
library(ggplot2)

# Data ----
points <- read.csv("inputs/S_L_allpointsdata2.csv") %>% 
  filter(!is.na(VegSimple)) %>% 
  mutate(VegSimple = recode(VegSimple,
                          "1" = "Evergreen",
                          "2" = "Deciduous",
                          "3" = "Desert",
                          "4" = "Arctic",
                          "5" = "Grassland"))

RBL_data <- read.csv("inputs/RBL_global_Dec2023.csv")

geometry <- sf::st_read("inputs/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")

geography_carnivora <- geometry %>% 
  filter(order_ == "CARNIVORA") %>% 
  filter(origin != 3)
geography_carnivora$sci_name <- gsub(" ", "_", geography_carnivora$sci_name)

# Community trait values ----

traitsByPoint <- summarize_traits_by_point(
  points_df = points,
  trait_df = RBL_data,
  species_polygons = geography_carnivora,
  trait_column = "RBL",
  species_name_col = "sci_name",
  continent = TRUE,
  n_cores = 4)

points_africa <- traitsByPoint$points %>% filter(continent == "Africa")

plot(points_africa$Longitude, points_africa$Latitude)

inspect_point_species(traitsByPoint)

# Quantitative environmental variables ----

precip_global <- ecometric_model(
  points_df = traitsByPoint$points,
  env_var = "BIO12",
  transform_fun = function(x) log(x + 1),
  inv_transform_fun = function(x) exp(x) - 1,
  min_species = 3,
  grid_bins_1 = 25,
  grid_bins_2 = 25
)

precip_africa <- ecometric_model(
  points_df = points_africa,
  env_var = "BIO12",
  transform_fun = function(x) log(x + 1),
  inv_transform_fun = function(x) exp(x) - 1,
  min_species = 3,
  grid_bins_1 = 25,
  grid_bins_2 = 25
)

temp_global <- ecometric_model(
  points_df = traitsByPoint$points,
  env_var = "BIO1",
  transform_fun = function(x) x/10,
  inv_transform_fun = function(x) x*10,
  min_species = 3,
  grid_bins_1 = 25,
  grid_bins_2 = 25
)

temp_africa <- ecometric_model(
  points_df = points_africa,
  env_var = "BIO1",
  transform_fun = function(x) x/10,
  inv_transform_fun = function(x) x*10,
  min_species = 3,
  grid_bins_1 = 25,
  grid_bins_2 = 25
)

# Qualitative environmental variables ----

veg_global <- ecometric_model_qual(
  points_df = traitsByPoint$points,
  category_col = "VegSimple",
  min_species = 3,
  grid_bins_1 = 25,
  grid_bins_2 = 25
)

veg_africa <- ecometric_model_qual(
  points_df = points_africa,
  category_col = "VegSimple",
  min_species = 3,
  grid_bins_1 = 25,
  grid_bins_2 = 25
)

ecometricModels <- list(precip_global, precip_africa, temp_global, 
                        temp_africa, veg_global, veg_africa)

names(ecometricModels) <- c("precip_global", "precip_africa", "temp_global", 
                            "temp_africa", "veg_global", "veg_africa")

saveRDS(ecometricModels, "outputs/ecometricmodels.RDS")
