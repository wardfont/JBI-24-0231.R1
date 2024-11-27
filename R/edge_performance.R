edge_performance <- function(species, repetitions,
                             biovars_models) {
  order <- expand.grid(
    model = biovars_models,
    repetition = repetitions,
    species = species
  ) %>%
    mutate(id = as.numeric(rownames(.))) %>%
    as.data.frame() %>%
    tibble()

  species_tables <- lapply(species, function(sp) {
    repetition_tables <- lapply(repetitions, function(rep) {
      inds <- order %>% filter(
        species == sp,
        repetition == rep
      )
      base <- tar_read(euforest_gam_biovars,
        branches = inds %>%
          filter(model == "base") %>%
          pull(id)
      )[[1]]

      ## base_edge <- base$gam$data_ens %>%
      ##   filter(pred > 0.025 & pred < edge_lims[2])

      base_edge <- base$gam$data_ens %>%
        filter(pred > 0.025 & pred < base$gam$performance$thr_value)

      edge_AUCs <- sapply(biovars_models, function(mod) {
        ind <- order %>%
          filter(
            species == sp,
            repetition == rep,
            model == mod
          ) %>%
          pull(id)
        mod_gam <- tar_read(euforest_gam_biovars,
          branches = ind
        )[[1]]

        edge_preds <- mod_gam$gam$data_ens %>%
          filter(rnames %in% base_edge$rnames)

        edge_AUC <- auc(
          response = edge_preds$pr_ab,
          predictor = edge_preds$pred
        )

        return(edge_AUC)
      })

      return(tibble(
        sp = sp, rep = rep,
        mod = biovars_models,
        edge_AUC = edge_AUCs
      ))
    })
    table <- do.call(rbind, repetition_tables)
  })

  edge_AUC_data <- do.call(rbind, species_tables)

  ## ggplot(edge_AUC_data %>%
  ##   mutate(
  ##     id = paste0(sp, rep),
  ##     mod = ordered(mod, levels = c("base", "stewart10", "stewart15", "stewart20", "stewart25", "zimmerman", "katz"))
  ##   )) +
  ##   geom_line(aes(x = mod, y = edge_AUC, group = id))
}
