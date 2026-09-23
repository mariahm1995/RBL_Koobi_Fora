library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(commecometrics)
library(sf)
library(stringr)
source("scripts/funcs.R")

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

ecometricModels <- readRDS("outputs/ecometricModels.RDS")

koobi_fora <- read.csv("inputs/Koobi_fora_10_21_25.csv") %>% 
  group_by(Time, Lat, Long)

koobi_fora |> count(Time)

richness_levels <- c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)
n_iterations <- 1000

rarefaction_results <- map_dfr(
  richness_levels,
  function(target_richness) {
    # target_richness <- richness_levels[5]
    
   koobi_fora %>%
      ungroup() %>%
      filter(!is.na(RBL)) %>%
      group_by(
        Time,
        Time_Name,
        Lat,
        Long
      ) %>%
      
      # Only retain sites that can support this richness
      filter(n() >= target_richness) %>%
      
      nest() %>%
      
      mutate(
        results = map(
          data,
          function(site_data) {
            
            map_dfr(
              seq_len(n_iterations),
              function(iter) {
                
                sampled <- site_data %>%
                  slice_sample(n = target_richness)
                
                tibble(
                  iteration = iter,
                  richness = target_richness,
                  
                  fossil_summ_trait_1 =
                    mean(sampled$RBL, na.rm = TRUE),
                  
                  fossil_summ_trait_2 =
                    sd(sampled$RBL, na.rm = TRUE)
                )
              }
            )
          }
        )
      ) %>%
      
      select(-data) %>%
      unnest(results)
  }
)

rarefaction_results %>%
  distinct(Time, richness) %>%
  arrange(Time, richness)

# ------------------------------------------------------------
# Reconstruct models
# ------------------------------------------------------------

model_names <- c("africa" = "precip_africa", "global" = "precip_global")

reconstructions <- model_names %>%
  map(~ reconstruct_env_fast(
    fossildata = rarefaction_results,
    model_out = ecometricModels[[.x]]
  ))

all_rare <- bind_rows(reconstructions, .id = "region") %>%
  mutate(region = str_to_title(region)) %>%
  ungroup()

analog_summary_all <- all_rare %>%
  group_by(region, Time, Time_Name, richness) %>%
  summarise(
    iterations = n(),
    valid = sum(!is.na(fossil_env_est_UN)),
    missing = sum(is.na(fossil_env_est_UN)),
    analog_percent = 100 * valid / iterations,
    .groups = "drop"
  )

plot_analog_curve <- function(data) {
  data %>%
    ggplot(aes(x = richness, y = analog_percent)) +
    geom_line() +
    geom_point(size = 1.5) +
    facet_grid(region ~ Time) +
    scale_x_continuous(
      labels = scales::label_number(accuracy = 1)
    ) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
    labs(x = "Species richness", y = "Fossil communities with a modern analog (%)") +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(hjust = 1, size = 8)
    )
}
plot_analog <- plot_analog_curve(analog_summary_all)
plot_analog

ggsave("figures/analog_curve.png",
       plot = plot_analog,
       width = 17, height = 12, units = "cm")

