reconstruct_env_fast <- function(
    fossildata,
    model_out,
    inv_transform = NULL,
    ci = 0.05
) {
  
  # ------------------------------------------------------------
  # Check required columns
  # ------------------------------------------------------------
  
  required_cols <- c(
    "fossil_summ_trait_1",
    "fossil_summ_trait_2"
  )
  
  missing_cols <- setdiff(
    required_cols,
    names(fossildata)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in fossildata: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  
  # ------------------------------------------------------------
  # Extract model information
  # ------------------------------------------------------------
  
  mbrks  <- model_out$diagnostics$brks_1
  sdbrks <- model_out$diagnostics$brks_2
  
  modern_points <- model_out$points_df
  
  env_col <- if (model_out$settings$transformed) {
    "env_trans"
  } else {
    model_out$settings$env_var
  }
  
  
  # ------------------------------------------------------------
  # Inverse transformation
  # ------------------------------------------------------------
  
  if (missing(inv_transform) || is.null(inv_transform)) {
    
    inv_transform <- model_out$settings$inv_transform_fun
    
  }
  
  
  # ------------------------------------------------------------
  # Assign fossil observations to ecometric bins
  # ------------------------------------------------------------
  
  fossildata <- fossildata %>%
    mutate(
      fossil_bin_1 = .bincode(
        fossil_summ_trait_1,
        breaks = mbrks
      ),
      
      fossil_bin_2 = .bincode(
        fossil_summ_trait_2,
        breaks = sdbrks
      )
    )
  
  
  # ------------------------------------------------------------
  # Only calculate bins actually used by fossils
  # ------------------------------------------------------------
  
  needed_bins <- fossildata %>%
    ungroup() %>%
    distinct(
      fossil_bin_1,
      fossil_bin_2
    ) %>%
    filter(
      !is.na(fossil_bin_1),
      !is.na(fossil_bin_2)
    )
  
  # ------------------------------------------------------------
  # Compute environmental reconstruction ONCE per bin
  # ------------------------------------------------------------
  
  bin_lookup <- needed_bins %>%
    rowwise() %>%
    mutate(
      
      reconstruction = list({
        
        dat <- modern_points %>%
          filter(
            bin_1 == fossil_bin_1,
            bin_2 == fossil_bin_2
          ) %>%
          pull(all_of(env_col))
        
        dat <- dat[!is.na(dat)]
        
        if (length(dat) < 2) {
          
          tibble(
            fossil_env_est = NA_real_,
            fossil_minlimit = NA_real_,
            fossil_maxlimit = NA_real_
          )
          
        } else {
          
          dens <- density(
            dat,
            bw = 1
          )
          
          mode_idx <- which.max(dens$y)
          
          bound_n <- floor(
            length(dens$x) * ci
          )
          
          lower_idx <- max(
            1,
            mode_idx - bound_n
          )
          
          upper_idx <- min(
            length(dens$x),
            mode_idx + bound_n
          )
          
          tibble(
            fossil_env_est =
              dens$x[mode_idx],
            
            fossil_minlimit =
              dens$x[lower_idx],
            
            fossil_maxlimit =
              dens$x[upper_idx]
          )
        }
        
      })
    ) %>%
    tidyr::unnest(reconstruction) %>%
    ungroup()
  
  
  # ------------------------------------------------------------
  # Join reconstruction to every fossil simulation
  # ------------------------------------------------------------
  
  fossildata <- fossildata %>%
    left_join(
      bin_lookup,
      by = c(
        "fossil_bin_1",
        "fossil_bin_2"
      )
    )
  
  
  # ------------------------------------------------------------
  # Back transform
  # ------------------------------------------------------------
  
  if (!is.null(inv_transform)) {
    
    fossildata <- fossildata %>%
      mutate(
        fossil_env_est_UN =
          inv_transform(fossil_env_est),
        
        fossil_minlimit_UN =
          inv_transform(fossil_minlimit),
        
        fossil_maxlimit_UN =
          inv_transform(fossil_maxlimit)
      )
  }
  
  
  return(fossildata)
}
