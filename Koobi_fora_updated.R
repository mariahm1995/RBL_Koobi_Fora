########################################################
# Dual Ecometric Models: Precipitation and Temperature
# Author: Maria A. Hurtado-Materon
# Date: [Date]
########################################################

# ====== Load Required Packages ======
library(dplyr)
library(purrr)
library(sf)
library(leaflet)
library(raster)
library(ggplot2)
library(commecometrics)


# ====== Load Data ======
points <- read.csv("inputs/S_L_allpointsdata2.csv")

continents <- raster::shapefile("inputs/continents.shp")
traits <- read.csv("inputs/RBL_global_Dec2023.csv", header = TRUE, na = ".", stringsAsFactors = FALSE)
traits$TaxonName <- gsub("_", " ", traits$TaxonName)

geometry <- sf::st_read("inputs/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp")
geography_carnivora <- geometry[geometry$order_ == "CARNIVORA", ]
rm(geometry)

# ====== GLOBAL ======
## ====== Summarize Traits at Sampling Points ======
res <- summarize_traits_by_point(
  points_df = points,
  trait_df = traits,
  species_polygons = geography_carnivora,
  species_name_col = "sci_name",
  trait_column = "RBL",
  continent_shp = continents,
  lon_col = "Longitude",
  lat_col = "Latitude"
)

## ====== Inspect Point Species Composition ======
inspect_point_species(
  res_list = res,
  min_species_valid = 3
)

## ====== Fossil Data ======
fossil <- read.csv("inputs/Koobi_fora_full.csv", header = TRUE)
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

## ====== Ecometric Model 1: Precipitation (BIO12) ======

EM_prec <- ecometric_model(
  points_df = res$points,
  env_var = "BIO12",
  transform_fun = function(x) log(x + 1),
  inv_transform_fun = function(x) exp(x) - 1,
  grid_bins = 25,
  min_species = 3
)

## Project fossils into precipitation ecometric space
fossil_results_prec <- reconstruct_env(
  fossildata = fossil_summary,
  model_out = EM_prec,
  match_nearest = TRUE,
  fossil_lon = "Long",
  fossil_lat = "Lat",
  modern_id = "GlobalID",
  modern_lon = "Longitude",
  modern_lat = "Latitude"
)

## Plot Precipitation Ecometric Space
EM_prec_plot <- ecometric_space(
  model_out = EM_prec,
  modern_color = "#4361ee",
  env_name = "",
  fossil_data = fossil_results_prec
)

text <- tibble(point = c("A","C","D","F"),
               x = fossil_results_prec$fossil_mbc[-c(2,5)] - 0.5,
               y = fossil_results_prec$fossil_sdc[-c(2,5)] - 0.5)

text2 <- tibble(point = c("B","E"),
               x = fossil_results_prec$fossil_mbc[c(2,5)] + c(-1, 0),
               y = fossil_results_prec$fossil_sdc[c(2,5)] + c(0.5, 1.5))

plot_prec_global <- EM_prec_plot + geom_text(data = text,
                                             aes(x = x, y = y, label = point),
                                             inherit.aes = FALSE, size = 3) + 
  geom_text(data = text2,
            aes(x = x, y = y, label = point),
            inherit.aes = FALSE, size = 3)

ggsave("figures/plot_prec_global.png", plot_prec_global, width = 4.8, height = 3.6)

## ====== Ecometric Model 2: Temperature (BIO1) ======

EM_temp <- ecometric_model(
  points_df = res$points,
  env_var = "BIO1",
  transform_fun = function(x) x / 10,
  inv_transform_fun = function(x) x * 10,
  grid_bins = 25,
  min_species = 3
)

# Project fossils into temperature ecometric space
fossil_results_temp <- reconstruct_env(
  fossildata = fossil_summary,
  model_out = EM_temp,
  match_nearest = TRUE,
  fossil_lon = "Long",
  fossil_lat = "Lat",
  modern_id = "GlobalID",
  modern_lon = "Longitude",
  modern_lat = "Latitude"
)

# Plot Temperature Ecometric Space
EM_temp_plot <- ecometric_space(
  model_out = EM_temp,
  env_name = "",
  fossil_data = fossil_results_temp,
  modern_color = "#4361ee",
  palette = c("#0081a7","#fdfcdc", "#f07167")
)

text <- tibble(point = c("A","C","D","F"),
               x = fossil_results_temp$fossil_mbc[-c(2,5)] - 0.5,
               y = fossil_results_temp$fossil_sdc[-c(2,5)] - 0.5)

text2 <- tibble(point = c("B","E"),
                x = fossil_results_temp$fossil_mbc[c(2,5)] + c(-1, 0),
                y = fossil_results_temp$fossil_sdc[c(2,5)] + c(0.5, 1.5))

plot_temp_global <- EM_temp_plot + geom_text(data = text,
                                             aes(x = x, y = y, label = point),
                                             inherit.aes = FALSE, size = 3) + 
  geom_text(data = text2,
            aes(x = x, y = y, label = point),
            inherit.aes = FALSE, size = 3)

ggsave("figures/plot_temp_global.png", plot_temp_global, width = 4.8, height = 3.6)

## ====== Maps ======

## ====== Sensitivity Analysis for Both Models ======
png(file="figures/sensitivity_global_prec.png", width = 8, height = 6, 
    units = "in", res = 300)
sensitivity_prec <- sensitivity_analysis(
  points_df = res$points,
  env_var = "BIO12",
  sample_sizes = seq(100, 10000, 1000),
  iterations = 20,
  test_split = 0.2,
  bins = 25,
  transform_fun = function(x) log(x + 1),
  parallel = TRUE
)
dev.off()

png(file="figures/sensitivity_global_temp.png", width = 8, height = 6,
    units = "in", res = 300)
sensitivity_temp <- sensitivity_analysis(
  points_df = res$points,
  env_var = "BIO1",
  sample_sizes = seq(100, 10000, 1000),
  iterations = 20,
  test_split = 0.2,
  bins = 25,
  transform_fun = function(x) x / 10,
  parallel = TRUE
)
dev.off()

# ====== AFRICA ======

res_africa <- res$points %>% filter(continent == "Africa")

## ====== Ecometric Model 1: Precipitation (BIO12) ======

EM_prec_africa <- ecometric_model(
  points_df = res_africa,
  env_var = "BIO12",
  transform_fun = function(x) log(x + 1),
  inv_transform_fun = function(x) exp(x) - 1,
  grid_bins = 25,
  min_species = 3
)

## Project fossils into precipitation ecometric space
fossil_results_prec_africa <- reconstruct_env(
  fossildata = fossil_summary,
  model_out = EM_prec_africa,
  match_nearest = TRUE,
  fossil_lon = "Long",
  fossil_lat = "Lat",
  modern_id = "GlobalID",
  modern_lon = "Longitude",
  modern_lat = "Latitude"
)

# Plot Precipitation Ecometric Space
EM_prec_plot_africa <- ecometric_space(
  model_out = EM_prec_africa,
  modern_color = "#4361ee",
  env_name = "",
  fossil_data = fossil_results_prec_africa
)

text <- tibble(point = c("A","B","C","E","F"),
               x = fossil_results_prec_africa$fossil_mbc[-c(4)] - 0.5,
               y = fossil_results_prec_africa$fossil_sdc[-c(4)] - 0.5)

text2 <- tibble(point = c("D"),
                x = fossil_results_prec_africa$fossil_mbc[c(4)]- 0.5,
                y = fossil_results_prec_africa$fossil_sdc[c(4)] + 0.5)

plot_prec_africa <- EM_prec_plot_africa + geom_text(data = text,
                                             aes(x = x, y = y, label = point),
                                             inherit.aes = FALSE, size = 3) + 
  geom_text(data = text2,
            aes(x = x, y = y, label = point),
            inherit.aes = FALSE, size = 3)

ggsave("figures/plot_prec_africa.png", plot_prec_africa, width = 4.8, height = 3.6)

## ====== Ecometric Model 2: Temperature (BIO1) ======

EM_temp_africa <- ecometric_model(
  points_df = res_africa,
  env_var = "BIO1",
  transform_fun = function(x) x / 10,
  inv_transform_fun = function(x) x * 10,
  grid_bins = 25,
  min_species = 3
)

# Project fossils into temperature ecometric space
fossil_results_temp_africa <- reconstruct_env(
  fossildata = fossil_summary,
  model_out = EM_temp_africa,
  match_nearest = TRUE,
  fossil_lon = "Long",
  fossil_lat = "Lat",
  modern_id = "GlobalID",
  modern_lon = "Longitude",
  modern_lat = "Latitude"
)

# Plot Temperature Ecometric Space
EM_temp_plot_africa <- ecometric_space(
  model_out = EM_temp_africa,
  env_name = "",
  fossil_data = fossil_results_temp_africa,
  modern_color = "#4361ee",
  palette = c("#0081a7","#fdfcdc", "#f07167")
)

text <- tibble(point = c("A","B","C","E","F"),
               x = fossil_results_temp_africa$fossil_mbc[-c(4)] - 0.5,
               y = fossil_results_temp_africa$fossil_sdc[-c(4)] - 0.5)

text2 <- tibble(point = c("D"),
                x = fossil_results_temp_africa$fossil_mbc[c(4)]- 0.5,
                y = fossil_results_temp_africa$fossil_sdc[c(4)] + 0.5)

plot_temp_africa <- EM_temp_plot_africa + geom_text(data = text,
                                                    aes(x = x, y = y, label = point),
                                                    inherit.aes = FALSE, size = 3) + 
  geom_text(data = text2,
            aes(x = x, y = y, label = point),
            inherit.aes = FALSE, size = 3)

ggsave("figures/plot_temp_africa.png", plot_temp_africa, width = 4.8, height = 3.6)

## ====== Maps ======

act_precip_fossil <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black") +
  theme_void() +
  geom_sf(data = EM_temp_africa$points_df, aes(color = env_trans), size = 2, pch = 16) +
  scale_color_gradientn(
    limits = c(-2, 10), 
    colors = rev(color), 
    name = "Observed precipitation Log(mm)", 
    na.value = "transparent", 
    breaks = c(0, 2, 4, 6, 8, 10), 
    labels = c(0, 2, 4, 6, 8, 10)
  ) +
  labs(color = "LogPrecip") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), 
               fill = NA, color = "black") +
  theme(
    legend.position = "top",  # ← move to top
    legend.direction = "horizontal",  # ← make horizontal
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", colour = "white")
  )

Time_A <- points.rest3 %>% filter(mbc3 == fossil.sf$fossilmbc1[1]) %>% filter(sdbc3 == fossil.sf$fossilsdbc1[1])
Time_BDE <- points.rest3 %>% filter(mbc3 == fossil.sf$fossilmbc1[2]) %>% filter(sdbc3 == fossil.sf$fossilsdbc1[2])
Time_C <- points.rest3 %>% filter(mbc3 == fossil.sf$fossilmbc1[3]) %>% filter(sdbc3 == fossil.sf$fossilsdbc1[3])

Time_A_points <- SpatialPointsDataFrame(Time_A[,2:3], Time_A)
sf_A_modern_points <- st_as_sf(Time_A_points)

Time_BDE_points <- SpatialPointsDataFrame(Time_BDE[,2:3], Time_BDE)
sf_BDE_modern_points <- st_as_sf(Time_BDE_points)

Time_C_points <- SpatialPointsDataFrame(Time_C[,2:3], Time_C)
sf_C_modern_points <- st_as_sf(Time_C_points)


Time_A_fig <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme_void()  +
  geom_sf(data = sf_A_modern_points, aes(color = "#c44536"), size = 2, pch = 16) +
  # scale_color_gradientn(limits=as.numeric(c(-2, 10)), colors = rev(color), name = "Observed \nPrecipitation \n Log (mm)", na.value = "transparent", breaks = c(0, 2, 4, 6, 8, 10), labels = c(0, 2, 4, 6, 8, 10)) +
  labs(title = "Time A") + theme(legend.position="none") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")
# theme(legend.position = c(0.1, 0.3), legend.title = element_text(size = 12), 
#       legend.text = element_text(size = 10), legend.background = element_rect(fill = "white", colour = "white"))

Time_BDE_fig <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme_void()  +
  geom_sf(data = sf_BDE_modern_points, aes(color = "#c44536"), size = 2, pch = 16) +
  labs(title = "Times B D E") + theme(legend.position="none") +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")

Time_C_fig <- ggplot() +
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  theme_void()  +
  geom_sf(data = sf_C_modern_points, aes(color = "#c44536"), size = 2, pch = 16) +
  labs(title = "Time C") + theme(legend.position="none")+
  geom_polygon(data = continents, aes(x = long, y = lat, group = group), fill = NA, color = "black")

combined_plot <- act_precip_fossil + Time_A_fig + Time_BDE_fig+ Time_C_fig

ggsave("figures/combined_plot_prec.png", plot = combined_plot, width = 12, height = 8, dpi = 300)


## ====== Sensitivity Analysis for Both Models ======
png(file="figures/sensitivity_africa_prec.png", width = 8, height = 6, 
    units = "in", res = 300)
sensitivity_prec <- sensitivity_analysis(
  points_df = res_africa,
  env_var = "BIO12",
  sample_sizes = seq(100, 10000, 1000),
  iterations = 20,
  test_split = 0.2,
  bins = 25,
  transform_fun = function(x) log(x + 1),
  parallel = TRUE
)
dev.off()

png(file="figures/sensitivity_africa_temp.png", width = 8, height = 6,
    units = "in", res = 300)
sensitivity_temp <- sensitivity_analysis(
  points_df = res_africa,
  env_var = "BIO1",
  sample_sizes = seq(100, 10000, 1000),
  iterations = 20,
  test_split = 0.2,
  bins = 25,
  transform_fun = function(x) x / 10,
  parallel = TRUE
)
dev.off()



