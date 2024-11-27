create_palette <- function() {
  source_palette <- palette.colors(palette = "Okabe-Ito")
  palette <- c(
    "B" = source_palette[8],
    "E" = source_palette[2],
    "N" = source_palette[3],
    "G" = source_palette[4]
  )
  names(palette) <- c("B", "E", "N", "G")

  return(palette)
}

name_models <- function(x) {
  model_names <- ordered(
    case_when(
      x == "base" ~ "B",
      x == "zimmerman" ~ "N",
      x == "stewart10" ~ "E10",
      x == "stewart15" ~ "E15",
      x == "stewart20" ~ "E20",
      x == "stewart25" ~ "E25",
      x == "katz" ~ "G",
      x == "katz2" ~ "G2"
    ),
    levels = c("B", "E10", "E15", "E20", "E25", "N", "G2", "G")
  )
}

performance_difference_to_B_comparison_boxplot <- function(metrics, edge_perf, palette, species_properties) {
  metrics <- metrics %>%
    tibble() %>%
    left_join(edge_perf %>%
      rename(species = sp, model = "mod"), by = c("species", "model", "rep"))

  variables <- c("AIC", "AUC", "edge_AUC")

  best_E <- metrics %>%
    tibble() %>%
    mutate(model = name_models(model)) %>%
    filter(model %in% c("E10", "E15", "E20", "E25")) %>%
    group_by(model, species) %>%
    summarise(mean_AUC = mean(AUC, na.rm = T)) %>%
    group_by(species) %>%
    filter(mean_AUC == max(mean_AUC, na.rm = T)) %>%
    mutate(comb = paste0(species, model))

  datas <- lapply(variables, function(variable) {
    data <- metrics %>%
      tibble() %>%
      mutate(id = factor(paste0(rep))) %>%
      rename(var = all_of(variable)) %>%
      mutate(model = name_models(model)) %>%
      filter(model %in% c("G", "N", "B") | paste0(species, model) %in% best_E$comb) %>%
      select(rep, species, var, id, model) %>%
      mutate(model = case_when(
        model %in% c("G", "N", "B") ~ model,
        T ~ "E"
      )) %>%
      pivot_wider(names_from = model, values_from = var) %>%
      mutate(
        diff_NB = N - B,
        diff_GB = G - B,
        diff_EB = E - B,
        diff_BG = B - G,
        diff_NG = N - G,
        diff_EG = E - G
      ) %>% # /get(comparison[2]))
      mutate(metric = variable)
  })

  summary_table <-
    rbind(
      datas[[1]] %>%
        mutate(metric = "AIC"),
      datas[[2]] %>%
        mutate(metric = "AUC"),
      datas[[2]] %>%
        mutate(metric = "edge_AUC")
    ) %>%
    group_by(species, metric) %>%
    summarise(
      B = median(B),
      G = median(G),
      E = median(E),
      N = median(N),
      diff_NB = median(diff_NB),
      diff_GB = median(diff_GB),
      diff_EB = median(diff_EB),
      diff_BG = median(diff_BG),
      diff_NG = median(diff_NG),
      diff_EG = median(diff_EG)
    )

  write.csv(summary_table, "./output/summary/summary-table.csv")

  tests <- lapply(variables, function(variable) {
    test <- metrics %>%
      tibble() %>%
      mutate(id = factor(paste0(rep))) %>%
      rename(var = all_of(variable)) %>%
      mutate(model = name_models(model)) %>%
      filter(model %in% c("G", "N", "B") | paste0(species, model) %in% best_E$comb) %>%
      select(rep, species, var, id, model) %>%
      mutate(model = case_when(
        model %in% c("G", "N", "B") ~ model,
        T ~ "E"
      )) %>%
      group_by(species) %>%
      mutate(model = ordered(model, levels = c("B", "G", "N", "E"))) %>%
      arrange(species, rep, model) %>%
      mutate(
        fr_p_diff_B = friedmanTest(var,
          groups = model, blocks = rep,
          alternative = "two.sided"
        )$p.value,
        fr_ny_p_diff_B = rep(c(
          NA,
          frdManyOneNemenyiTest(var,
            groups = model,
            blocks = rep,
            alternative = "two.sided"
          )$p.value
        ), length(unique(rep))),
        n = length(unique(rep))
      ) %>%
      mutate(model = ordered(model, levels = c("G", "B", "E", "N"))) %>%
      arrange(species, rep, model) %>%
      mutate(
        fr_p_diff_G = friedmanTest(var,
          groups = model, blocks = rep,
          alternative = "two.sided"
        )$p.value,
        fr_ny_p_diff_G = rep(c(
          NA,
          frdManyOneNemenyiTest(var,
            groups = model,
            blocks = rep,
            alternative = "two.sided"
          )$p.value
        ), length(unique(rep))),
        n = length(unique(rep))
      ) %>%
      select(-rep, -var, -id) %>%
      distinct() %>%
      mutate(
        sign_B = fr_p_diff_B <= 0.05 & fr_ny_p_diff_B <= 0.05,
        sign_G = fr_p_diff_G <= 0.05 & fr_ny_p_diff_G <= 0.05
      ) %>%
      select(species, model, sign_B, sign_G) %>%
      rename(comparison = model) %>%
      mutate(metric = variable)
  })

  species_order <- datas[[2]] %>%
    group_by(species) %>%
    summarise(median = median(diff_GB)) %>%
    left_join(species_properties, by = "species") %>%
    arrange(Leaf.type, -median) %>%
    pull(species)

  plot_data <- do.call(rbind, datas) %>%
    left_join(species_properties, by = "species") %>%
    mutate(Leaf.type = case_when(
      Leaf.type == "B" ~ "Broadleaf",
      Leaf.type == "N" ~ "Needle"
    )) %>%
    pivot_longer(c("diff_NB", "diff_GB", "diff_EB"),
      names_to = "comparison",
      values_to = "diff"
    ) %>%
    mutate(comparison = case_when(
      comparison == "diff_NB" ~ "N",
      comparison == "diff_GB" ~ "G",
      comparison == "diff_EB" ~ "E"
    )) %>%
    select(species, Leaf.type, diff, comparison, metric) %>%
    left_join(do.call(rbind, tests) %>% select(species, comparison, metric, sign_B),
      by = c("species", "comparison", "metric")
    ) %>%
    mutate(metric = factor(metric,
      levels = c("AIC", "AUC", "edge_AUC"),
      labels = c(
        expression(Delta * " " * AIC[model]),
        expression(Delta * " AUC"),
        expression(Delta * " " * AUC[edge])
      )
    )) %>%
    mutate(
      comparison = ordered(comparison,
        levels = c("G", "N", "E")
      ),
      sign_B = factor(sign_B, levels = c(TRUE, FALSE))
    )

  plot <- ggplot(
    data = plot_data %>%
      mutate(species = ordered(species, levels = species_order)),
    aes(x = interaction(species, Leaf.type), y = diff, col = comparison, fill = sign_B)
  ) +
    geom_hline(yintercept = 0) +
    geom_boxplot() +
    labs(x = NULL) +
    theme_pubr() +
    facet_wrap(~metric,
      labeller = label_parsed,
      strip.position = "left", scales = "free_y", ncol = 1
    ) +
    ylab(NULL) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(
      values = palette,
      labels = c("GEV - Mean", "Normal - Mean", "Quantile - Mean")
    ) +
    scale_fill_manual(
      values = c("grey75", "white"),
      labels = c("Sign.", "Not sign."),
      drop = T
    ) +
    labs(col = NULL, fill = NULL) +
    scale_x_discrete(guide = "axis_nested") +
    theme(
      axis.text.x = element_text(face = "italic"),
      panel.background = element_rect(colour = "black", size = 0.25)
    ) +
    rotate_x_text(angle = 90)

  return(list(
    plot = plot,
    summary_table = summary_table
  ))
}

performance_difference_to_G_comparison_boxplot <- function(metrics, edge_perf, palette, species_properties) {
  metrics <- metrics %>%
    tibble() %>%
    left_join(edge_perf %>%
      rename(species = sp, model = "mod"), by = c("species", "model", "rep"))

  variables <- c("AIC", "AUC", "edge_AUC")

  best_E <- metrics %>%
    tibble() %>%
    mutate(model = name_models(model)) %>%
    filter(model %in% c("E10", "E15", "E20", "E25")) %>%
    group_by(model, species) %>%
    summarise(mean_AUC = mean(AUC, na.rm = T)) %>%
    group_by(species) %>%
    filter(mean_AUC == max(mean_AUC, na.rm = T)) %>%
    mutate(comb = paste0(species, model))

  datas <- lapply(variables, function(variable) {
    data <- metrics %>%
      tibble() %>%
      mutate(id = factor(paste0(rep))) %>%
      rename(var = all_of(variable)) %>%
      mutate(model = name_models(model)) %>%
      filter(model %in% c("G", "N", "B") | paste0(species, model) %in% best_E$comb) %>%
      select(rep, species, var, id, model) %>%
      mutate(model = case_when(
        model %in% c("G", "N", "B") ~ model,
        T ~ "E"
      )) %>%
      pivot_wider(names_from = model, values_from = var) %>%
      mutate(
        diff_BG = B - G,
        diff_NG = N - G,
        diff_EG = E - G
      ) %>% # /get(comparison[2]))
      mutate(metric = variable)
  })

  tests <- lapply(variables, function(variable) {
    test <- metrics %>%
      tibble() %>%
      mutate(id = factor(paste0(rep))) %>%
      rename(var = all_of(variable)) %>%
      mutate(model = name_models(model)) %>%
      filter(model %in% c("G", "N", "B") | paste0(species, model) %in% best_E$comb) %>%
      select(rep, species, var, id, model) %>%
      mutate(model = case_when(
        model %in% c("G", "N", "B") ~ model,
        T ~ "E"
      )) %>%
      group_by(species) %>%
      mutate(model = ordered(model, levels = c("B", "G", "N", "E"))) %>%
      arrange(species, rep, model) %>%
      mutate(
        fr_p_diff_B = friedmanTest(var,
          groups = model, blocks = rep,
          alternative = "two.sided"
        )$p.value,
        fr_ny_p_diff_B = rep(c(
          NA,
          frdManyOneNemenyiTest(var,
            groups = model,
            blocks = rep,
            alternative = "two.sided"
          )$p.value
        ), length(unique(rep))),
        n = length(unique(rep))
      ) %>%
      mutate(model = ordered(model, levels = c("G", "B", "E", "N"))) %>%
      arrange(species, rep, model) %>%
      mutate(
        fr_p_diff_G = friedmanTest(var,
          groups = model, blocks = rep,
          alternative = "two.sided"
        )$p.value,
        fr_ny_p_diff_G = rep(c(
          NA,
          frdManyOneNemenyiTest(var,
            groups = model,
            blocks = rep,
            alternative = "two.sided"
          )$p.value
        ), length(unique(rep))),
        n = length(unique(rep))
      ) %>%
      select(-rep, -var, -id) %>%
      distinct() %>%
      mutate(
        sign_B = fr_p_diff_B <= 0.05 & fr_ny_p_diff_B <= 0.05,
        sign_G = fr_p_diff_G <= 0.05 & fr_ny_p_diff_G <= 0.05
      ) %>%
      select(species, model, sign_B, sign_G) %>%
      rename(comparison = model) %>%
      mutate(metric = variable)
  })

  species_order <- datas[[2]] %>%
    group_by(species) %>%
    summarise(median = median(diff_BG)) %>%
    left_join(species_properties, by = "species") %>%
    arrange(Leaf.type, median) %>%
    pull(species)

  plot_data <- do.call(rbind, datas) %>%
    left_join(species_properties, by = "species") %>%
    mutate(Leaf.type = case_when(
      Leaf.type == "B" ~ "Broadleaf",
      Leaf.type == "N" ~ "Needle"
    )) %>%
    pivot_longer(c("diff_BG", "diff_NG", "diff_EG"),
      names_to = "comparison",
      values_to = "diff"
    ) %>%
    mutate(
      comparison =
        case_when(
          comparison == "diff_BG" ~ "B",
          comparison == "diff_NG" ~ "N",
          comparison == "diff_EG" ~ "E"
        )
    ) %>%
    select(species, Leaf.type, diff, comparison, metric) %>%
    left_join(do.call(rbind, tests) %>% select(species, comparison, metric, sign_G),
      by = c("species", "comparison", "metric")
    ) %>%
    mutate(metric = factor(metric,
      levels = c("AIC", "AUC", "edge_AUC"),
      labels = c(
        expression(Delta * " " * AIC[model]),
        expression(Delta * " AUC"),
        expression(Delta * " " * AUC[edge])
      )
    )) %>%
    mutate(
      comparison = ordered(comparison,
        levels = c("B", "E", "N")
      ),
      sign_G = factor(sign_G, levels = c(TRUE, FALSE))
    )

  plot <- ggplot(
    data = plot_data %>%
      mutate(species = ordered(species, levels = species_order)),
    aes(x = interaction(species, Leaf.type), y = diff, col = comparison, fill = sign_G)
  ) +
    geom_hline(yintercept = 0) +
    geom_boxplot() +
    labs(x = NULL) +
    theme_pubr() +
    facet_wrap(~metric,
      labeller = label_parsed,
      strip.position = "left", scales = "free_y", ncol = 1
    ) +
    ylab(NULL) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(face = "bold", size = 12)
    ) +
    scale_color_manual(
      values = palette,
      labels = c("Mean - GEV", "Quantile - GEV", "Normal - GEV")
    ) +
    scale_fill_manual(
      values = c("grey75", "white"),
      labels = c("Sign.", "Not sign."),
      drop = T
    ) +
    labs(col = NULL, fill = NULL) +
    scale_x_discrete(guide = "axis_nested") +
    theme(
      axis.text.x = element_text(face = "italic"),
      panel.background = element_rect(colour = "black", size = 0.25)
    ) +
    rotate_x_text(angle = 90)

  return(list(
    plot = plot
  ))
}

get_cor_plot <- function(sdm_file, biovars_file) {
  r <- c(rast(sdm_file), crop(rast(biovars_file), rast(sdm_file)))
  plots <- lapply(c("tmaxhm", "tmincm", "pptqdq"), function(var) {
    cor <- layerCor(subset(r, grepl(paste0(var, "|bio"), names(r))),
      "pearson",
      na.rm = T, maxcell = 100000
    )$pearson

    nms <- c(
      "Mean",
      "Std. Dev",
      "10 year anomaly",
      "15 year anomaly",
      "20 year anomaly",
      "25 year anomaly",
      "Location",
      "Scale",
      "Shape",
      paste("Bio", 1:19)
    )

    colnames(cor) <- nms
    rownames(cor) <- nms
    g <- ggcorrplot(cor,
      lab_size = 1,
      ## method = "circle",
      title = case_when(
        var == "tmaxhm" ~ "Max. temperature hottest month",
        var == "tmincm" ~ "Min. temperature coldest month",
        var == "pptqdq" ~ "Precipitation driest quarter"
      )
    )
    ggsave(
      paste0("./plot/", var, "corplot.svg"),
      g,
    )
  })

  return(plots)
  ## ggarrange(plotlist = plots, ncol = 3, nrow = 1, common.legend = T)
}

per_pixel_delta_aic_plot <- function(terraclimate_delta_aic) {
  r <- rast(terraclimate_delta_aic)

  pal <- wes_palette("Zissou1")

  tmax <- ggplot() +
    geom_spatraster(data = r[[1]]) +
    theme_pubr() +
    scale_fill_gradient2(
      high = pal[1],
      mid = "grey90",
      low = pal[5],
      na.value = NA,
      limits = c(-1, 1) * as.numeric(global(abs(r[[1]]), fun = "max", na.rm = T)[1])
    ) +
    labs(title = "TxHm", fill = expression(Delta * "AIC"[distr.])) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )

  tmin <- ggplot() +
    geom_spatraster(data = r[[2]]) +
    theme_pubr() +
    scale_fill_gradient2(
      high = pal[1],
      mid = "grey90",
      low = pal[5],
      na.value = NA,
      limits = c(-1, 1) * as.numeric(global(abs(r[[2]]), fun = "max", na.rm = T)[1])
    ) +
    labs(title = "TnCm", fill = expression(Delta * "AIC"[distr.])) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )

  ppt <- ggplot() +
    geom_spatraster(data = r[[3]]) +
    theme_pubr() +
    scale_fill_gradient2(
      high = pal[1],
      mid = "grey90",
      low = pal[5],
      na.value = NA,
      limits = c(-1, 1) * as.numeric(global(abs(r[[3]]), fun = "max", na.rm = T)[1]),
      breaks = c(-150, 0, 150)
    ) +
    labs(title = "PDq", fill = expression(Delta * "AIC"[distr.])) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )

  p_delta_aic <- ggarrange(
    tmax, tmin, ppt,
    ncol = 3,
    labels = "auto"
  )

  return(ggplotGrob(p_delta_aic))
}
