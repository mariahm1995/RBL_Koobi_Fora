library(dplyr)
library(sf)
library(ggplot2)
library(sp)
library(rnaturalearth)

sf_use_s2(FALSE)
continents <- ne_coastline(returnclass = "sf") %>%
  st_crop(ymin = -60, ymax = 90, xmin = -180, xmax = 180)

africa_shp <- rnaturalearth::ne_countries(scale = "medium", 
                                             returnclass = "sf") %>% 
  filter(continent == "Africa")

ecometricModels <- readRDS("outputs/ecometricModels.RDS")

koobi_fora <- read.csv("outputs/koobi_paleo_summary_nov.csv")

points <- read.csv("inputs/S_L_allpointsdata2.csv") %>% 
  select(GlobalID, Longitude, Latitude)

# Figure 1
koobi_fora_sf <- st_as_sf(
  data.frame(lon = 36.2, lat = 3.95),
  coords = c("lon", "lat"),
  crs = 4326
)

africa <- ggplot() +
  geom_sf(data = africa_shp, fill = "#AAAAAA", color = "white") +  # base map
  geom_sf(data = subset(africa_shp, sovereignt == "Kenya"),
    fill = "#e56b6f", color = "white") + theme_void()

ggsave("figures/Africa.png",width = 15, height = 13, units = "cm")

Kenya <- ggplot() +
  geom_sf(data = subset(africa_shp, sovereignt == "Kenya"), fill = "#e56b6f", color = "white") + 
  geom_sf(data = koobi_fora_sf, color = "black", size = 5) + theme_void()

ggsave("figures/Kenya.png",width = 15, height = 13, units = "cm")

# Figure 5

Time_A <- ecometricModels[["precip_global"]]$points_df %>%
  filter(bin_1 == koobi_fora$global_fossil_bin_1[1]) %>% 
  filter(bin_2 == koobi_fora$global_fossil_bin_2[1]) %>% 
  mutate(time = "A")
  
Time_B <- ecometricModels[["precip_global"]]$points_df %>%
  filter(bin_1 == koobi_fora$global_fossil_bin_1[2]) %>% 
  filter(bin_2 == koobi_fora$global_fossil_bin_2[2]) %>% 
  mutate(time = "B")

Time_C <- ecometricModels[["precip_global"]]$points_df %>%
  filter(bin_1 == koobi_fora$global_fossil_bin_1[3]) %>% 
  filter(bin_2 == koobi_fora$global_fossil_bin_2[3])  %>% 
  mutate(time = "C")

Time_D <- ecometricModels[["precip_global"]]$points_df %>%
  filter(bin_1 == koobi_fora$global_fossil_bin_1[4]) %>% 
  filter(bin_2 == koobi_fora$global_fossil_bin_2[4])  %>% 
  mutate(time = "D")

Time_E <- ecometricModels[["precip_global"]]$points_df %>%
  filter(bin_1 == koobi_fora$global_fossil_bin_1[5]) %>% 
  filter(bin_2 == koobi_fora$global_fossil_bin_2[5])  %>% 
  mutate(time = "E")

times_df <- rbind(Time_A, Time_B, Time_C, Time_D, Time_E)

times_points <- SpatialPointsDataFrame(times_df[,29:30], times_df) %>% st_as_sf()
st_crs(times_points)<- st_crs(continents)

times_col <- c("A" = "#0c4767",
               "B" = "#71b340",
               "C" = "#F59E0B",
               "D" = "#EF4444",
               "E" = "#8B5CF6")

times_map <- ggplot() +
  geom_sf(data = times_points, aes(color = time), size = 0.8, pch = 16, alpha = 0.6) +
  geom_sf(data = continents, fill = NA, color = "black", size = 0.2) +
  scale_color_manual(values = times_col, breaks = names(times_col), name = "Times") +
  theme_void() + labs(color = "Times") +
  theme(legend.position = c(0.1, 0.1),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = grid::unit(6, "mm")) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1, shape = 16)))

ggsave("figures/times_map.png", plot = times_map, width = 6, height = 3, dpi = 300)

Time_B_africa <- ecometricModels[["precip_africa"]]$points_df %>%
  filter(bin_1 == koobi_fora$africa_fossil_bin_1[2]) %>% 
  filter(bin_2 == koobi_fora$africa_fossil_bin_2[2]) %>% 
  mutate(time = "B")

Time_C_africa <- ecometricModels[["precip_africa"]]$points_df %>%
  filter(bin_1 == koobi_fora$africa_fossil_bin_1[3]) %>% 
  filter(bin_2 == koobi_fora$africa_fossil_bin_2[3])  %>% 
  mutate(time = "C")

Time_D_africa <- ecometricModels[["precip_africa"]]$points_df %>%
  filter(bin_1 == koobi_fora$africa_fossil_bin_1[4]) %>% 
  filter(bin_2 == koobi_fora$africa_fossil_bin_2[4])  %>% 
  mutate(time = "D")

Time_E_africa <- ecometricModels[["precip_africa"]]$points_df %>%
  filter(bin_1 == koobi_fora$africa_fossil_bin_1[5]) %>% 
  filter(bin_2 == koobi_fora$africa_fossil_bin_2[5])  %>% 
  mutate(time = "E")

times_df_africa <- rbind(Time_B_africa, Time_C_africa, Time_D_africa, Time_E_africa)

times_points_africa <- SpatialPointsDataFrame(times_df_africa[,29:30], times_df_africa) %>% st_as_sf()
st_crs(times_points_africa) <- st_crs(continents)

times_col_africa <- c("B" = "#71b340",
                      "C" = "#F59E0B",
                      "D" = "#EF4444",
                      "E" = "#8B5CF6")

times_map_africa <- ggplot() +
  geom_sf(data = times_points_africa, aes(color = time), size = 0.8, pch = 16, alpha = 0.6) +
  geom_sf(data = africa_shp, fill = NA, color = "black", size = 0.2) +
  scale_color_manual(values = times_col_africa, breaks = names(times_col_africa), name = "Times") +
  theme_void() + labs(color = "Times") +
  theme(legend.position = c(0.12, 0.1),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = grid::unit(6, "mm")) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1, shape = 16)))

ggsave("figures/times_map_africa.png", plot = times_map_africa, width = 6, height = 3, dpi = 300)

# Anomaly maps
precip_globa_sf <- SpatialPointsDataFrame(ecometricModels[["precip_global"]]$points_df[,29:30],
                                       ecometricModels[["precip_global"]]$points_df) %>% st_as_sf()
temp_globa_sf <- SpatialPointsDataFrame(ecometricModels[["temp_global"]]$points_df[,29:30],
                                          ecometricModels[["temp_global"]]$points_df) %>% st_as_sf()
veg_globa_sf <- SpatialPointsDataFrame(ecometricModels[["veg_global"]]$points_df[,29:30],
                                          ecometricModels[["veg_global"]]$points_df) %>% st_as_sf()

precip_africa_sf <- SpatialPointsDataFrame(ecometricModels[["precip_africa"]]$points_df[,29:30],
                                          ecometricModels[["precip_africa"]]$points_df) %>% st_as_sf()
temp_africa_sf <- SpatialPointsDataFrame(ecometricModels[["temp_africa"]]$points_df[,29:30],
                                        ecometricModels[["temp_africa"]]$points_df) %>% st_as_sf()
veg_africa_sf <- SpatialPointsDataFrame(ecometricModels[["veg_africa"]]$points_df[,29:30],
                                       ecometricModels[["veg_africa"]]$points_df) %>% st_as_sf()


st_crs(temp_africa_sf) <- st_crs(veg_africa_sf) <- st_crs(precip_africa_sf) <- st_crs(temp_globa_sf) <- st_crs(veg_globa_sf) <- st_crs(precip_globa_sf) <- st_crs(continents)

# Define shared color palettes
anomaly_col <- c("#1dd3b0", "#fdf0d5", "#858ae3")
veg_col <- c("Yes" = "#fdf0d5", "No" = "#858ae3")

# ---------------- GLOBAL MAPS ---------------- #

plot_precip_global <- ggplot() +
  geom_sf(data = precip_globa_sf, aes(color = env_anom),
          size = 0.5, pch = 16) +
  geom_sf(data = continents, fill = NA, color = "black", size = 0.1) +
  scale_color_gradientn(
    limits = c(-8, 7),
    colors = anomaly_col,
    name = "Anomaly",
    na.value = "transparent",
    breaks = c(-8, 0, 7),
    labels = c(-8, 0, 7)
  ) + theme_void() +theme(legend.position = "none")


ggsave("figures/anomaly_precip_global.png",
       plot = plot_precip_global,
       width = 10, height = 8, units = "cm")


plot_temp_global <- ggplot() +
  geom_sf(data = temp_globa_sf, aes(color = env_anom),
          size = 0.5, pch = 16) +
  geom_sf(data = continents, fill = NA, color = "black", size = 0.1) +
  scale_color_gradientn(
    limits = c(-44, 40),
    colors = anomaly_col,
    name = "Anomaly",
    na.value = "transparent",
    breaks = c(-44, 0, 40),
    labels = c(-44, 0, 40)
  ) + theme_void() +theme(legend.position = "none")

ggsave("figures/anomaly_temp_global.png",
       plot = plot_temp_global,
       width = 10, height = 8, units = "cm")


plot_veg_global <- ggplot() +
  geom_sf(data = veg_globa_sf, aes(color = factor(correct_prediction)),
          size = 0.5, pch = 16) +
  geom_sf(data = continents, fill = NA, color = "black", size = 0.1) +
  scale_color_manual(values = veg_col, name = "Correct prediction") +
  theme_void() +theme(legend.position = "none")

ggsave("figures/anomaly_veg_global.png",
       plot = plot_veg_global,
       width = 10, height = 8, units = "cm")


# ---------------- AFRICA MAPS ---------------- #

plot_precip_africa <- ggplot() +
  geom_sf(data = precip_africa_sf, aes(color = env_anom),
          size = 0.8, pch = 16) +
  geom_sf(data = africa_shp, fill = NA, color = "black", size = 0.2) +
  scale_color_gradientn(
    limits = c(-8, 7),
    colors = anomaly_col,
    name = "Anomaly",
    na.value = "transparent",
    breaks = c(-8, 0, 7),
    labels = c(-8, 0, 7)
  ) +
  theme_void()

ggsave("figures/anomaly_precip_africa.png",
       plot = plot_precip_africa,
       width = 10, height = 8, units = "cm")


plot_temp_africa <- ggplot() +
  geom_sf(data = temp_africa_sf, aes(color = env_anom),
          size = 0.8, pch = 16) +
  geom_sf(data = africa_shp, fill = NA, color = "black", size = 0.2) +
  scale_color_gradientn(
    limits = c(-44, 40),
    colors = anomaly_col,
    name = "Anomaly",
    na.value = "transparent",
    breaks = c(-44, 0, 40),
    labels = c(-44, 0, 40)
  ) +
  theme_void()

ggsave("figures/anomaly_temp_africa.png",
       plot = plot_temp_africa,
       width = 10, height = 8, units = "cm")


plot_veg_africa <- ggplot() +
  geom_sf(data = veg_africa_sf, aes(color = factor(correct_prediction)),
          size = 0.8, pch = 16) +
  geom_sf(data = africa_shp, fill = NA, color = "black", size = 0.2) +
  scale_color_manual(values = veg_col, name = "Correct prediction") +
  theme_void()

ggsave("figures/anomaly_veg_africa.png",
       plot = plot_veg_africa,
       width = 10, height = 8, units = "cm")

