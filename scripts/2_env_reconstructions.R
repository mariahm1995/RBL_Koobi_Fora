library(sf)
library(dplyr)
library(commecometrics)
library(ggplot2)

ecometricModels <- readRDS("outputs/ecometricModels.RDS")

koobi_fora <- read.csv("inputs/Koobi_fora_10_21_25.csv") %>% 
  group_by(Time, Lat, Long) %>%
  summarise(fossil_summ_trait_1 = mean(RBL, na.rm = TRUE),
            fossil_summ_trait_2 = sd(RBL, na.rm = TRUE))

# --- Paleo reconstructions at Koobi Fora ------------------------------
# Helper that keeps only the columns we need (and avoids name clashes later)
pick_quant_cols <- function(df) {
  df %>%
    dplyr::select(
      Time, Lat, Long,
      fossil_env_est_UN, fossil_minlimit_UN, fossil_maxlimit_UN,
      bin_1, bin_2, fossil_bin_1, fossil_bin_2
    )
}

pick_qual_cols <- function(df) {
  df %>%
    dplyr::select(
      Time, Lat, Long,
      fossil_env_est,
      bin_1, bin_2, fossil_bin_1, fossil_bin_2
    )
}

# Precipitation
koobi_global_precip <- reconstruct_env(
  fossildata  = koobi_fora,
  model_out   = ecometricModels[["precip_global"]],
  fossil_lon  = "Long",
  fossil_lat  = "Lat",
  modern_id   = "GlobalID",
  modern_lat  = "Latitude",
  modern_lon  = "Longitude"
)

koobi_africa_precip <- reconstruct_env(
  fossildata  = koobi_fora,
  model_out   = ecometricModels[["precip_africa"]],
  fossil_lon  = "Long",
  fossil_lat  = "Lat",
  modern_id   = "GlobalID",
  modern_lat  = "Latitude",
  modern_lon  = "Longitude"
)

# Temperature
koobi_global_temp <- reconstruct_env(
  fossildata  = koobi_fora,
  model_out   = ecometricModels[["temp_global"]],
  fossil_lon  = "Long",
  fossil_lat  = "Lat",
  modern_id   = "GlobalID",
  modern_lat  = "Latitude",
  modern_lon  = "Longitude"
) 

koobi_africa_temp <- reconstruct_env(
  fossildata  = koobi_fora,
  model_out   = ecometricModels[["temp_africa"]],
  fossil_lon  = "Long",
  fossil_lat  = "Lat",
  modern_id   = "GlobalID",
  modern_lat  = "Latitude",
  modern_lon  = "Longitude"
)

# Vegetation
koobi_global_veg <- reconstruct_env_qual(
  fossildata  = koobi_fora,
  model_out   = ecometricModels[["veg_global"]],
  fossil_lon  = "Long",
  fossil_lat  = "Lat",
  modern_id   = "GlobalID",
  modern_lat  = "Latitude",
  modern_lon  = "Longitude"
) 

koobi_africa_veg <- reconstruct_env_qual(
  fossildata  = koobi_fora,
  model_out   = ecometricModels[["veg_africa"]],
  fossil_lon  = "Long",
  fossil_lat  = "Lat",
  modern_id   = "GlobalID",
  modern_lat  = "Latitude",
  modern_lon  = "Longitude"
) 

# --- Build one-row-per-Time tables (no summarising needed) ------------

# Quantitative: just select/rename the three columns we need
prep_quant <- function(df, prefix) {
  df %>%
    dplyr::select(
      Time,
      !!paste0(prefix, "_est") := fossil_env_est_UN,
      !!paste0(prefix, "_min") := fossil_minlimit_UN,
      !!paste0(prefix, "_max") := fossil_maxlimit_UN,
      !!paste0(prefix, "_bin_1") := bin_1,
      !!paste0(prefix, "_bin_2") := bin_2,
      !!paste0(prefix, "_fossil_bin_1") := fossil_bin_1,
      !!paste0(prefix, "_fossil_bin_2") := fossil_bin_2
    )
}

# Qualitative: take the estimated category as-is
prep_qual <- function(df, prefix) {
  df %>%
    dplyr::select(
      Time,
      !!paste0(prefix, "_veg") := fossil_env_est,
      !!paste0(prefix, "_bin_1") := bin_1,
      !!paste0(prefix, "_bin_2") := bin_2,
      !!paste0(prefix, "_fossil_bin_1") := fossil_bin_1,
      !!paste0(prefix, "_fossil_bin_2") := fossil_bin_2
    )
}

temp_global_tbl  <- prep_quant(pick_quant_cols(koobi_global_temp),  "temp_global")
temp_africa_tbl  <- prep_quant(pick_quant_cols(koobi_africa_temp),  "temp_africa")
prec_global_tbl  <- prep_quant(pick_quant_cols(koobi_global_precip),"precip_global")
prec_africa_tbl  <- prep_quant(pick_quant_cols(koobi_africa_precip),"precip_africa")
veg_global_tbl   <- prep_qual(pick_qual_cols(koobi_global_veg), "veg_global")
veg_africa_tbl   <- prep_qual(pick_qual_cols(koobi_africa_veg),   "veg_africa")

# --- Final 6-row summary table ----------------------------------------
koobi_paleo_summary <- temp_global_tbl %>%
  dplyr::full_join(temp_africa_tbl, by = "Time") %>%
  dplyr::full_join(prec_global_tbl, by = "Time") %>%
  dplyr::full_join(prec_africa_tbl, by = "Time") %>%
  dplyr::full_join(veg_global_tbl,  by = "Time") %>%
  dplyr::full_join(veg_africa_tbl,  by = "Time") %>%
  dplyr::arrange(Time) %>%
  dplyr::select(
    Time,
    temp_global_est, temp_global_min, temp_global_max,
    temp_africa_est, temp_africa_min, temp_africa_max,
    precip_global_est, precip_global_min, precip_global_max,
    precip_africa_est, precip_africa_min, precip_africa_max,
    veg_global_veg, veg_africa_veg, temp_africa_bin_1,
    temp_africa_bin_2, temp_africa_fossil_bin_1, temp_africa_fossil_bin_2,
    temp_global_bin_1, temp_global_bin_2, temp_global_fossil_bin_1, temp_global_fossil_bin_2
  )

print(koobi_paleo_summary)

write.csv(koobi_paleo_summary, "outputs/koobi_paleo_summary.csv", row.names = FALSE)

# Ecometric spaces ----

text_global_veg <- tibble(point = c("A","B","C", "D", "E", "F"),
                          x = koobi_global_veg$fossil_bin_1 - 0.5,
                          y = koobi_global_veg$fossil_bin_2 - 0.5)

text_africa <- tibble(point = c("A","B","C", "D", "E", "F"),
                      x = koobi_africa_temp$fossil_bin_1 - 0.5,
                      y = koobi_africa_temp$fossil_bin_2 - 0.5)


pal_biomes <- c("#56B4E9", "#f4a261", "#99582a", "#a7c957", "#CC79A7")
pal_biomes_afri <- c("#f4a261", "#99582a", "#a7c957", "#CC79A7")

plot_veg_global <- ecometric_space_qual(
  model_out = ecometricModels[["veg_global"]],
  palette = pal_biomes, 
  fossil_data = koobi_global_veg,
  modern_color = "white",
  x_label = NULL,
  y_label = NULL)$ecometric_space_plot + 
  geom_text(data = text_global_veg, 
            aes(x = x, y = y, label = point),
            inherit.aes = FALSE, size = 4, 
            color = "black") +
  theme(axis.text = element_text(size = 18))

ggsave("figures/ecospace_veg_global.png",width = 15, height = 12, units = "cm")

plot_veg_africa <- ecometric_space_qual(
  model_out = ecometricModels[["veg_africa"]],
  palette = pal_biomes_afri,
  fossil_data = koobi_africa_veg,
  modern_color = "white",
  x_label = NULL,
  y_label = NULL)$ecometric_space_plot + 
  geom_text(data = text_africa, 
            aes(x = x, y = y, label = point),
            inherit.aes = FALSE, size = 4, 
            color = "black") +
  theme(axis.text = element_text(size = 18))

ggsave("figures/ecospace_veg_africa.png",width = 15, height = 12, units = "cm")

pal_precp <- c("#99521E", "#D38B5D", "#EFEBCE", "#739E82", "#2C5530")

plot_prec_africa <- ecometric_space(
  model_out   = ecometricModels[["precip_africa"]],
  palette = pal_precp,
  env_name    = "Precipitation",
  fossil_data = koobi_africa_precip,
  modern_color = "white",
  x_label     = NULL,
  y_label     = NULL) + 
  geom_text(data = text_africa, 
              aes(x = x, y = y, label = point),
              inherit.aes = FALSE, size = 4, 
            color = "black") +
  theme(axis.text = element_text(size = 18))

ggsave("figures/ecospace_precip_africa.png",width = 15, height = 12, units = "cm")

plot_prec_global <- ecometric_space(
  model_out   = ecometricModels[["precip_global"]],
  palette = pal_precp,
  env_name    = "Precipitation",
  fossil_data = koobi_global_precip,
  modern_color = "white",
  x_label     = NULL,
  y_label     = NULL) + 
  geom_text(data = text_global_veg, 
              aes(x = x, y = y, label = point),
              inherit.aes = FALSE, size = 4, 
            color = "black") +
  theme(axis.text = element_text(size = 18))

ggsave("figures/ecospace_precip_global.png",width = 15, height = 12, units = "cm")

# --- Temperature -------------------------------------------------------
pal_temp <- c("#005f73","#0a9396","#94d2bd","#e9d8a6","#ca6702","#bb3e03","#ae2012")

plot_temp_africa <- ecometric_space(
  model_out   = ecometricModels[["temp_africa"]],
  palette = pal_temp,
  env_name    = "Temperature",
  fossil_data = koobi_africa_temp,
  modern_color = "white",
  x_label     = NULL,
  y_label     = NULL) + 
  geom_text(data = text_africa, 
              aes(x = x, y = y, label = point),
              inherit.aes = FALSE, size = 4, 
            color = "black") +
  theme(axis.text = element_text(size = 18))

ggsave("figures/ecospace_temp_africa.png",width = 15, height = 12, units = "cm")

plot_temp_global <- ecometric_space(
  model_out   = ecometricModels[["temp_global"]],
  palette = pal_temp,
  env_name    = "Temperature",
  fossil_data = koobi_global_temp,
  modern_color = "white",
  x_label     = NULL,
  y_label     = NULL) +   
  geom_text(data = text_global_veg, 
                aes(x = x, y = y, label = point),
                inherit.aes = FALSE, size = 4, 
            color = "black") +
  theme(axis.text = element_text(size = 18))

ggsave("figures/ecospace_temp_global.png",width = 15, height = 12, units = "cm")

lapply(names(ecometricModels)[1:2], function(i){
  ecometricModel <- ecometricModels[[i]]
  ecometric_space(
    model_out   = ecometricModel,
    palette = pal_temp,
    env_name    = "Precipitation",
    fossil_data = koobi_global_temp,
    modern_color = "white",
    x_label     = NULL,
    y_label     = NULL) +   
    geom_text(data = text_global_veg, 
              aes(x = x, y = y, label = point),
              inherit.aes = FALSE, size = 4, 
              color = "black") +
    theme(axis.text = element_text(size = 18))
  
  ggsave(paste0("figures/ecospace_", names(ecometricModels[i]), ".png"),
         width = 15, height = 12, units = "cm")
})

lapply(names(ecometricModels)[3:4], function(i){
  ecometricModel <- ecometricModels[[i]]
  ecometric_space(
    model_out   = ecometricModel,
    palette = pal_temp,
    env_name    = "Temperature",
    fossil_data = koobi_global_temp,
    modern_color = "white",
    x_label     = NULL,
    y_label     = NULL) +   
    geom_text(data = text_global_veg, 
              aes(x = x, y = y, label = point),
              inherit.aes = FALSE, size = 4, 
              color = "black") +
    theme(axis.text = element_text(size = 18))
  
  ggsave(paste0("figures/ecospace_", names(ecometricModels[i]), ".png"),
         width = 15, height = 12, units = "cm")
})

# lapply(names(ecometricModels)[1:2], function(i){
#   ecometricModel <- ecometricModels[[i]]
#   ecometric_space(
#     model_out   = ecometricModel,
#     palette = pal_temp,
#     env_name    = "Precipitation",
#     fossil_data = koobi_global_temp,
#     modern_color = "white",
#     x_label     = NULL,
#     y_label     = NULL) +   
#     geom_text(data = text_global_veg, 
#               aes(x = x, y = y, label = point),
#               inherit.aes = FALSE, size = 4, 
#               color = "black") +
#     theme(axis.text = element_text(size = 18))
#   
#   ggsave(paste0("figures/ecospace_", names(ecometricModels[i]), ".png"),
#          width = 15, height = 12, units = "cm")
# })