one_spatial_rel_abs_plot_limited_ensemble <- function(katz_ensemble_gam_biovars_predictions,
                                                      base_ensemble_gam_biovars_predictions,
                                                      euforest,
                                                      species,
                                                      display_species,
                                                      palette,
                                                      terraclimate_sdm) {
  katz <- rast(lapply(
    1:length(display_species),
    function(i) {
      mask(
        rast(katz_ensemble_gam_biovars_predictions[[i]]),
        vect(euforest[[which(species == display_species[[i]])]]$ca)
      )
    }
  ))
  names(katz) <- display_species
  base <- rast(lapply(
    1:length(display_species),
    function(i) {
      mask(
        rast(base_ensemble_gam_biovars_predictions[[i]]),
        vect(euforest[[which(species == display_species[[i]])]]$ca)
      )
    }
  ))

  pal_abs <- wes_palette("Zissou1")
  low_ind_abs <- 5
  high_ind_abs <- 1
  pal_rel <- wes_palette("Darjeeling1")
  low_ind_rel <- 4
  high_ind_rel <- 2

  eu <- rast(terraclimate_sdm[[1]]$file)
  eu[!is.na(eu)] <- 1
  eu_pol <- as.polygons(eu)

  base_lim <- rep(NA, 4)
  sp_limits <- rep(NA, 4)

  plots_abs <- lapply(1:length(display_species), function(i) {
    ca <- vect(euforest[[which(species == display_species[[i]])]]$ca)
    pol <- crop(ca, eu_pol)

    diff <- katz[[i]] - base[[i]]

    background <- base[[i]]
    base_lim[i] <- global(background, quantile, probs = 0.95, na.rm = T)[[1]]
    background[background > base_lim[i]] <- base_lim[i]

    sp_limits[i] <- base_lim[i] * 0.66
    diff[diff > sp_limits[i]] <- sp_limits[i]
    diff[diff < -1 * sp_limits[i]] <- -1 * sp_limits[i]

    ggplot() +
      geom_spatraster(data = background) +
      scale_fill_gradient(
        low = "white",
        high = "grey25",
        na.value = NA,
        limits = c(0, base_lim[i])
      ) +
      new_scale_fill() +
      geom_spatraster(data = diff) +
      scale_fill_gradient2(
        low = pal_abs[low_ind_abs],
        mid = NA,
        high = pal_abs[high_ind_abs],
        midpoint = 0,
        limits = c(-1, 1) * sp_limits[i],
        na.value = NA
      ) +
      geom_spatvector(data = pol, fill = NA) +
      labs(
        x = display_species[i],
      ) +
      theme_pubr() +
      theme(
        axis.title.x = element_text(
          face = "bold.italic",
          size = 20
        ),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
  })

  plots_rel <- lapply(1:length(display_species), function(i) {
    ca <- vect(euforest[[which(species == display_species[[i]])]]$ca)
    pol <- crop(ca, eu_pol)

    diff <- katz[[i]] - base[[i]]
    diff[diff < 0] <- diff / (base[[i]])
    diff[diff > 0] <- diff / (1 - base[[i]])

    background <- base[[i]]
    base_lim[i] <- global(background, quantile, probs = 0.95, na.rm = T)[[1]]
    background[background > base_lim[i]] <- base_lim[i]

    sp_limits[i] <- global(abs(diff), fun = "max", na.rm = T)$max
    diff[diff > sp_limits[i]] <- sp_limits[i]
    diff[diff < -1 * sp_limits[i]] <- -1 * sp_limits[i]

    ggplot() +
      geom_spatraster(data = background) +
      scale_fill_gradient(
        low = "white",
        high = "grey25",
        na.value = NA,
        limits = c(0, base_lim[i])
      ) +
      new_scale_fill() +
      geom_spatraster(data = diff) +
      scale_fill_gradient2(
        low = pal_rel[low_ind_rel],
        mid = NA,
        high = pal_rel[high_ind_rel],
        midpoint = 0,
        limits = c(-1, 1) * sp_limits[i],
        na.value = NA
      ) +
      geom_spatvector(data = pol, fill = NA) +
      labs(
        x = display_species[i],
      ) +
      theme_pubr() +
      theme(
        axis.title.x = element_text(
          face = "bold.italic",
          size = 20
        ),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
  })

  heights <- sapply(1:length(display_species), function(i) {
    ca <- vect(euforest[[which(species == display_species[[i]])]]$ca)
    ext(ca)[4] - ext(ca)[3]
  })

  ## leg <- get_legend(plots_rel[[1]])
  ## ggsave("rel-legend.svg", leg)

  one <- ggarrange(
    plotlist = c(plots_abs, plots_rel), ncol = 4, nrow = 2,
    common.legend = T
  )

  ggsave("./plot/one_spatial/limit_rel_abs_one_spatial.svg", one,
    scale = 1, height = 297, width = 210 * 2, units = "mm"
  ) # A4

  return(ggplotGrob(one))
}

one_spatial_rel_abs_plot_other_ensemble <- function(katz_ensemble_gam_biovars_predictions,
                                                    base_ensemble_gam_biovars_predictions,
                                                    euforest,
                                                    species,
                                                    display_species,
                                                    palette,
                                                    terraclimate_sdm) {
  katz <- rast(lapply(
    1:length(display_species),
    function(i) {
      mask(
        rast(katz_ensemble_gam_biovars_predictions[[i]]),
        vect(euforest[[which(species == display_species[[i]])]]$ca)
      )
    }
  ))
  names(katz) <- display_species
  base <- rast(lapply(
    1:length(display_species),
    function(i) {
      mask(
        rast(base_ensemble_gam_biovars_predictions[[i]]),
        vect(euforest[[which(species == display_species[[i]])]]$ca)
      )
    }
  ))

  pal_abs <- wes_palette("Zissou1")
  low_ind_abs <- 5
  high_ind_abs <- 1
  pal_rel <- wes_palette("Darjeeling1")
  low_ind_rel <- 4
  high_ind_rel <- 2

  eu <- rast(terraclimate_sdm[[1]]$file)
  eu[!is.na(eu)] <- 1
  eu_pol <- as.polygons(eu)

  base_lim <- rep(NA, 4)
  sp_limits <- rep(NA, 4)

  plots_abs <- lapply(1:length(display_species), function(i) {
    ca <- vect(euforest[[which(species == display_species[[i]])]]$ca)
    pol <- crop(ca, eu_pol)

    diff <- katz[[i]] - base[[i]]

    background <- base[[i]]
    base_lim[i] <- global(background, quantile, probs = 0.95, na.rm = T)[[1]]
    background[background > base_lim[i]] <- base_lim[i]

    sp_limits[i] <- base_lim[i] * 0.66
    diff[diff > sp_limits[i]] <- sp_limits[i]
    diff[diff < -1 * sp_limits[i]] <- -1 * sp_limits[i]

    ggplot() +
      geom_spatraster(data = background) +
      scale_fill_gradient(
        low = "white",
        high = "grey25",
        na.value = NA,
        limits = c(0, base_lim[i])
      ) +
      new_scale_fill() +
      geom_spatraster(data = diff) +
      scale_fill_gradient2(
        low = pal_abs[low_ind_abs],
        mid = NA,
        high = pal_abs[high_ind_abs],
        midpoint = 0,
        limits = c(-1, 1) * sp_limits[i],
        na.value = NA
      ) +
      geom_spatvector(data = pol, fill = NA) +
      labs(
        x = display_species[i],
      ) +
      theme_pubr() +
      theme(
        axis.title.x = element_text(
          face = "bold.italic",
          size = 20
        ),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
  })

  plots_rel <- lapply(1:length(display_species), function(i) {
    ca <- vect(euforest[[which(species == display_species[[i]])]]$ca)
    pol <- crop(ca, eu_pol)

    diff <- katz[[i]] - base[[i]]
    diff[diff < 0] <- diff / (base[[i]])
    diff[diff > 0] <- diff / (1 - base[[i]])

    background <- base[[i]]
    base_lim[i] <- global(background, quantile, probs = 0.95, na.rm = T)[[1]]
    background[background > base_lim[i]] <- base_lim[i]

    sp_limits[i] <- global(abs(diff), fun = "max", na.rm = T)$max
    diff[diff > sp_limits[i]] <- sp_limits[i]
    diff[diff < -1 * sp_limits[i]] <- -1 * sp_limits[i]

    ggplot() +
      geom_spatraster(data = background) +
      scale_fill_gradient(
        low = "white",
        high = "grey25",
        na.value = NA,
        limits = c(0, base_lim[i])
      ) +
      new_scale_fill() +
      geom_spatraster(data = diff) +
      scale_fill_gradient2(
        low = pal_rel[low_ind_rel],
        mid = NA,
        high = pal_rel[high_ind_rel],
        midpoint = 0,
        limits = c(-1, 1) * sp_limits[i],
        na.value = NA
      ) +
      geom_spatvector(data = pol, fill = NA) +
      labs(
        x = display_species[i],
      ) +
      theme_pubr() +
      theme(
        axis.title.x = element_text(
          face = "bold.italic",
          size = 20
        ),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
      )
  })

  one <- ggarrange(
    plotlist = c(plots_abs, plots_rel), nrow = 2, ncol = 24,
    common.legend = T
  )
  ggsave("./plot/one_spatial/other_rel_abs_one_spatial.svg", one,
    scale = 1, height = 297, width = 210 * 12, units = "mm", limitsize = F
  )

  return(ggplotGrob(one))
}
