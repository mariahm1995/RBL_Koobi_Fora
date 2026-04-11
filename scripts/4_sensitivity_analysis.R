library(commecometrics)
library(dplyr)

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

png(file="figures/sensitivity_global_prec.png", width = 8, height = 6, 
    units = "in", res = 300)
sensitivity_prec <- sensitivity_analysis(
  points_df = traitsByPoint$points,
  env_var = "BIO12",
  sample_sizes = seq(100, 10000, 1000),
  iterations = 20,
  test_split = 0.2,
  grid_bins_1 = 25,
  grid_bins_2 = 25,
  transform_fun = function(x) log(x + 1),
  parallel = TRUE,
  n_cores = 6
)
dev.off()

png(file="figures/sensitivity_africa_prec.png", width = 8, height = 6,
    units = "in", res = 300)
sensitivity_temp <- sensitivity_analysis(
  points_df = points_africa,
  env_var = "BIO12",
  sample_sizes = seq(100, 10000, 1000),
  iterations = 20,
  test_split = 0.2,
  grid_bins_1 = 25,
  grid_bins_2 = 25,
  transform_fun = function(x) log(x + 1),
  parallel = TRUE,
  n_cores = 6
)
dev.off()

png(file="figures/sensitivity_global_temp.png", width = 8, height = 6, 
    units = "in", res = 300)
sensitivity_prec <- sensitivity_analysis(
  points_df = traitsByPoint$points,
  env_var = "BIO1",
  sample_sizes = seq(100, 10000, 1000),
  iterations = 20,
  test_split = 0.2,
  grid_bins_1 = 25,
  grid_bins_2 = 25,
  transform_fun = function(x) x/10,
  parallel = TRUE,
  n_cores = 8
)
dev.off()

png(file="figures/sensitivity_africa_temp.png", width = 8, height = 6, 
    units = "in", res = 300)
sensitivity_prec <- sensitivity_analysis(
  points_df = points_africa,
  env_var = "BIO1",
  sample_sizes = seq(100, 10000, 1000),
  iterations = 20,
  test_split = 0.2,
  grid_bins_1 = 25,
  grid_bins_2 = 25,
  transform_fun = function(x) x/10,
  parallel = TRUE,
  n_cores = 8
)
dev.off()

png(file="figures/sensitivity_global_veg.png", width = 8, height = 6, 
    units = "in", res = 300)
sensitivity_prec <- sensitivity_analysis_qual(
  points_df = traitsByPoint$points,
  category_col = "VegSimple",
  sample_sizes = seq(100, 10000, 1000),
  iterations = 20,
  test_split = 0.2,
  grid_bins_1 = 25,
  grid_bins_2 = 25,
  parallel = TRUE,
  n_cores = 8
)
dev.off()

png(file="figures/sensitivity_africa_veg.png", width = 8, height = 6, 
    units = "in", res = 300)
sensitivity_prec <- sensitivity_analysis_qual(
  points_df = points_africa,
  category_col = "VegSimple",
  sample_sizes = seq(100, 10000, 1000),
  iterations = 20,
  test_split = 0.2,
  grid_bins_1 = 25,
  grid_bins_2 = 25,
  parallel = TRUE,
  n_cores = 8
)
dev.off()
