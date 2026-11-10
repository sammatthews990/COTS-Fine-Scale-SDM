# Report figures from the frozen calibration summaries made by script 22.
build_metacommunity_figures <- function(cfg) {
  read_mc_table <- function(n) utils::read.csv(file.path(cfg$dirs$tables, n), stringsAsFactors = FALSE)
  tr <- read_mc_table("metacommunity_trajectories.csv")
  obs <- readRDS(file.path(cfg$dirs$cache, "metacommunity_observations.rds"))
  reefs <- c("Lizard Island Reef", "MacGillivray Reef", "North Direction Reef", "Eyrie Reef")
  tr$reef_name <- factor(tr$reef_name, reefs)
  obs$reef_name <- factor(obs$reef_name, reefs)
  tr$candidate <- factor(tr$candidate, c("W", "C3a", "C6"))
  cols <- c(W = layer_colours[["legacy_clean_year2026"]],
            C3a = layer_colours[["legacy_rg_year2026"]], C6 = layer_colours[["hmsc_joint"]])
  p <- ggplot2::ggplot(tr, ggplot2::aes(.data$year, .data$simulated_cpue_mean, colour = .data$candidate, linetype = .data$candidate)) +
    ggplot2::geom_hline(yintercept = 0.22, linetype = "dotted", colour = ink$muted) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(0, .data$simulated_cpue_mean - .data$simulated_cpue_sd),
      ymax = .data$simulated_cpue_mean + .data$simulated_cpue_sd, fill = .data$candidate), alpha = 0.14, colour = NA) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(data = obs, ggplot2::aes(.data$year, .data$cotsptow), inherit.aes = FALSE,
      colour = ink$primary, size = 1.2) +
    ggplot2::facet_wrap(~reef_name, ncol = 2, scales = "free_y") +
    ggplot2::scale_colour_manual(values = cols, name = "Candidate") +
    ggplot2::scale_fill_manual(values = cols, name = "Candidate") +
    ggplot2::scale_linetype_manual(values = c(W = "solid", C3a = "longdash", C6 = "dotdash"), name = "Candidate") +
    ggplot2::labs(x = "LTMP reporting year", y = "COTS per manta tow",
      title = "Calibration candidates retain mismatches in outbreak height",
      subtitle = "Lines: five-seed means; ribbons: +/- one seed SD; points: LTMP manta observations",
      caption = "Linear conversion: 0.015 COTS/tow per adult/ha. Dotted line: 0.22 COTS/tow. Calibration, not held-out validation.") +
    theme_report() + ggplot2::theme(legend.position = "bottom")
  save_figure(p, cfg, "fig18_metacommunity_candidates.png", width = 7.5, height = 5.4)

  late <- read_mc_table("metacommunity_2026_trajectories.csv")
  late$reef_name <- factor(late$reef_name, reefs)
  late$coral <- factor(late$coral, c("earlier thermal settings", "revised thermal settings"))
  p <- ggplot2::ggplot(late, ggplot2::aes(.data$year, .data$adults_ha, colour = .data$dataset, linetype = .data$dataset)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_grid(reef_name ~ coral, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c(original = cols[["W"]], corrected = cols[["C3a"]]),
      labels = c(original = "Original outbreak boundary", corrected = "Corrected outbreak boundary"), name = NULL) +
    ggplot2::scale_linetype_manual(values = c(original = "solid", corrected = "longdash"),
      labels = c(original = "Original outbreak boundary", corrected = "Corrected outbreak boundary"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = c(2018, 2022, 2026)) +
    ggplot2::labs(x = "Model year", y = "Adults per hectare of site area",
      title = "The 2020s response depends on both coral forcing and the outbreak boundary",
      subtitle = "V0.4 to 2026; one seed pair; chlorophyll modifiers off; culling off",
      caption = "Corrected boundary uses cull/SALAD evidence also used for assessment: this is a diagnostic, not independent validation.") +
    theme_report() + ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8))
  save_figure(p, cfg, "fig19_metacommunity_2026_diagnostic.png", width = 7.5, height = 7.8)
}
