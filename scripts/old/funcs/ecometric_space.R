#' Plot Ecometric Trait Space with Optional Fossil Overlay
#'
#' Visualizes the binned trait–environment ecometric space and optionally overlays fossil trait bins.
#'
#' @param model_out Output from \code{run_ecometric_model()}, containing \code{raster_df} and \code{diagnostics}.
#' @param env_name Character string for the environmental variable label (e.g., "Precipitation (mm)").
#' @param fossil_data Optional. A data frame from \code{reconstruct_past_environment()}, including \code{fossil_mbc} and \code{fossil_sdc}.
#' @param fossil_color Color for fossil bin overlay (default = "#c44536").
#' @param palette Color palette for the continuous gradient (default = c("#bc6c25", "#fefae0", "#606c38")).
#'
#' @return A ggplot2 object showing the ecometric space.
#'
#' @import ggplot2
#' @export
ecometric_space <- function(model_out,
                            env_name = "Environment",
                            fossil_data = NULL,
                            fossil_color = "black",
                            modern_color = "#bc4749",
                            palette = c("#bc6c25", "#fefae0", "#606c38")) {
  
  # Extract model outputs
  raster_df <- model_out$eco_space
  mbreaks <- model_out$diagnostics$mbrks
  sd_breaks <- model_out$diagnostics$sdbrks
  grid_bins <- length(mbreaks) - 1
  
  # Axes formatting
  middle_idx <- if (grid_bins %% 2 == 0) { (grid_bins / 2) + 1 } else { ceiling(grid_bins / 2) }
  
  x_breaks <- c(mbreaks[1], mbreaks[middle_idx], mbreaks[grid_bins - 1])
  x_labels <- round(x_breaks, 2)
  x_pos <- c(0.5, middle_idx - 0.5, grid_bins - 0.5)
  
  y_breaks <- c(sd_breaks[1], sd_breaks[middle_idx], sd_breaks[grid_bins - 1])
  y_labels <- round(y_breaks, 2)
  y_pos <- c(0.5, middle_idx - 0.5, grid_bins - 0.5)
  
  # Base plot
  ecospace <- ggplot(raster_df, aes(x = x, y = y, fill = layer)) +
    geom_raster() +
    scale_fill_gradientn(colors = palette, name = env_name, na.value = "transparent") +
    scale_x_continuous(name = "Mean", breaks = x_pos, labels = x_labels, expand = c(0, 0), limits = c(0, grid_bins)) +
    scale_y_continuous(name = "SD", breaks = y_pos, labels = y_labels, expand = c(0, 0), limits = c(0, grid_bins)) +
    coord_fixed() +
    theme_bw()
  
  # Optional fossil overlay
  if (!is.null(fossil_data)) {
    ecospace <- ecospace +
      geom_rect(data = fossil_data,
                aes(xmin = as.numeric(fossil_mbc) - 1,
                    xmax = as.numeric(fossil_mbc),
                    ymin = as.numeric(fossil_sdc) - 1,
                    ymax = as.numeric(fossil_sdc)),
                inherit.aes = FALSE,
                colour = fossil_color, alpha = 0, size = 1) +
      geom_rect(data = fossil_data,
                aes(xmin = as.numeric(mbc) - 1,
                    xmax = as.numeric(mbc),
                    ymin = as.numeric(sdc) - 1,
                    ymax = as.numeric(sdc)),
                inherit.aes = FALSE,
                colour = modern_color, alpha = 0, size = 1)
  }
  
  return(ecospace)
}
