library(ggplot2)
library(dplyr)
library(tidyr)

all_data <- read.csv("inputs/comparasion_data.csv")

delta_car_AP_afri <- diff(all_data$car_AP_afri)
delta_car_AP_global <- diff(all_data$car_AP_glob)
delta_her_AP <- diff(all_data$herb_AP)
delta_her_AP_K <- diff(all_data$herb_AP_K)

plot_df <- data.frame(
  # interval = 1:(nrow(all_data)-1),
  interval = c("1-2", "2-3", "3-4", "4-5", "5-6"),
  carnivore <- diff(all_data$car_AP_glob),
  herbivore <- diff(all_data$herb_AP)
)
colnames(plot_df) <- c("interval", "carnivore", "herbivore")
# plot_df <- plot_df %>%
#   mutate(match = sign(delta_car_AP_global) == sign(delta_her_AP))

# Convert to long format

plot_long <- plot_df %>%
  pivot_longer(cols = c(carnivore, herbivore),
               names_to = "model",
               values_to = "delta")

# Plot
MAT <- ggplot(plot_long, aes(x = factor(interval), y = delta,
                      fill = model, group = model)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.6,
           color = "black",
           aes(linetype = model),
           size = 0.8) +
  scale_fill_manual(values = c("grey40", "grey80")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(x = "Time Interval",
       y = expression(Delta~AP),
       fill = "Reconstruction",
       linetype = "Reconstruction") +
  theme_minimal() 

ggsave("figures/AP_direction.png", MAT, width = 6, height = 5)
