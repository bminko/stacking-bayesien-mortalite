# Figures fusionnees destinees au memoire.
#
# Ce script lit uniquement des resultats deja valides. Il ne relance aucun
# ajustement Stan, aucune simulation et aucun calcul predictif.
#
# Contrat visuel 1
# - Question : la recuperation des poids et la performance predictive changent-
#   elles selon le scenario simule ?
# - Forme : matrice 3 x 4 de boxplots (RMSE, LogS et CRPS).
# - Comparaison : 30 repetitions, trois methodes, quatre scenarios.
#
# Contrat visuel 2
# - Question : ou se situent l'heterogeneite du meilleur modele individuel et
#   les ecarts entre poids hierarchiques et non regularises ?
# - Forme : heatmap LFO et cinq heatmaps d'ecart de poids avec echelle commune.
#
# Contrat visuel 3
# - Question : comment les scores predictifs varient-ils selon l'age et
#   l'horizon ?
# - Forme : matrice 3 x 2 de courbes avec une legende commune.
#
# Palette commune des methodes : bleu (global), orange (contextuel non
# regularise), vert (hierarchique). Les types de ligne et symboles fournissent
# une distinction supplementaire independante de la couleur.

options(stringsAsFactors = FALSE, scipen = 999)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "figures_memoire")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Fichier source absent : ", path, call. = FALSE)
  }
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

method_order <- c(
  "stacking_global",
  "stacking_contextual",
  "stacking_hierarchical"
)
method_axis_labels <- c(
  stacking_global = "Global",
  stacking_contextual = "Non rég.",
  stacking_hierarchical = "Hiér."
)
method_labels <- c(
  stacking_global = "Stacking global",
  stacking_contextual = "Contextuel non régularisé",
  stacking_hierarchical = "Contextuel hiérarchique"
)
method_colors <- c(
  stacking_global = "#4477AA",
  stacking_contextual = "#EE7733",
  stacking_hierarchical = "#228833"
)
method_lines <- c(
  stacking_global = 2,
  stacking_contextual = 3,
  stacking_hierarchical = 1
)
method_points <- c(
  stacking_global = 1,
  stacking_contextual = 17,
  stacking_hierarchical = 16
)

scenario_order <- c(
  "constant_weights",
  "horizon_only",
  "age_horizon",
  "low_information"
)
scenario_labels <- c(
  constant_weights = "Poids constants",
  horizon_only = "Poids selon l'horizon",
  age_horizon = "Poids selon l'âge et l'horizon",
  low_information = "Faible information"
)

model_order <- c("LC", "RH", "APC", "CBD", "M6")
model_colors <- c(
  LC = "#4477AA",
  RH = "#CC6677",
  APC = "#DDCC77",
  CBD = "#117733",
  M6 = "#AA4499"
)

render_pair <- function(filename, width, height, draw) {
  pdf_path <- file.path(output_dir, paste0(filename, ".pdf"))
  png_path <- file.path(output_dir, paste0(filename, ".png"))

  grDevices::cairo_pdf(
    pdf_path,
    width = width,
    height = height,
    family = "sans",
    bg = "white"
  )
  draw()
  grDevices::dev.off()

  grDevices::png(
    png_path,
    width = width,
    height = height,
    units = "in",
    res = 300,
    type = "cairo",
    family = "sans",
    bg = "white"
  )
  draw()
  grDevices::dev.off()

  invisible(c(pdf = pdf_path, png = png_path))
}

draw_quiet_grid <- function(horizontal = TRUE, vertical = FALSE) {
  graphics::grid(
    nx = if (vertical) NULL else NA,
    ny = if (horizontal) NULL else NA,
    col = "#E3E7EA",
    lty = 1
  )
}

draw_box_panel <- function(
  data, metric, scenario, row_index, column_index, total_rows
) {
  part <- data[
    data$scenario == scenario & data$method %in% method_order,
    ,
    drop = FALSE
  ]
  part$method <- factor(part$method, levels = method_order)

  bottom_margin <- if (row_index == total_rows) 4.0 else 1.1
  left_margin <- if (column_index == 1L) 4.2 else 2.0
  top_margin <- if (row_index == 1L) 2.5 else 1.0
  graphics::par(
    mar = c(bottom_margin, left_margin, top_margin, 0.7),
    las = 1,
    bty = "l",
    col.axis = "#30343B",
    col.lab = "#30343B",
    col.main = "#20242A",
    fg = "#30343B"
  )

  split_values <- split(part[[metric]], part$method)
  values <- unlist(split_values, use.names = FALSE)
  padding <- diff(range(values)) * 0.08
  if (!is.finite(padding) || padding == 0) {
    padding <- max(abs(values), 1) * 0.05
  }
  y_limits <- range(values) + c(-padding, padding)
  graphics::boxplot(
    split_values,
    names = unname(method_axis_labels[method_order]),
    col = unname(method_colors[method_order]),
    border = "#30343B",
    medcol = "#20242A",
    medlwd = 2,
    whisklty = 2,
    staplewex = 0.6,
    outline = FALSE,
    ylim = y_limits,
    axes = FALSE
  )
  draw_quiet_grid(horizontal = TRUE)
  graphics::axis(
    1,
    at = seq_along(method_order),
    labels = unname(method_axis_labels[method_order]),
    cex.axis = 0.78
  )

  ticks <- pretty(y_limits, n = 4)
  tick_labels <- if (metric == "mean_crps") {
    formatC(ticks, format = "f", digits = 5)
  } else if (metric == "mean_logs") {
    formatC(ticks, format = "f", digits = 2)
  } else {
    formatC(ticks, format = "f", digits = 2)
  }
  graphics::axis(2, at = ticks, labels = tick_labels, cex.axis = 0.76)
  graphics::box(bty = "l")

  if (row_index == 1L) {
    graphics::title(
      main = scenario_labels[[scenario]],
      font.main = 2,
      cex.main = 0.96,
      line = 0.7
    )
  }
  if (column_index == 1L) {
    row_label <- c(
      weight_rmse = "RMSE des poids",
      mean_logs = "LogS moyen",
      mean_crps = "CRPS moyen"
    )[[metric]]
    graphics::mtext(
      row_label,
      side = 2,
      line = 2.8,
      las = 0,
      font = 2,
      cex = 0.88
    )
  }
  if (row_index == total_rows) {
    graphics::mtext("Méthode", side = 1, line = 2.6, cex = 0.82)
  }
}

load_simulation_results <- function() {
  recovery <- read_csv(file.path(
    root, "results", "full", "simulation", "main",
    "weight_recovery_by_repetition.csv"
  ))
  performance <- read_csv(file.path(
    root, "results", "full", "simulation", "main",
    "performance_by_repetition.csv"
  ))

  merge(
    recovery[c("scenario", "repetition", "method", "weight_rmse")],
    performance[c("scenario", "repetition", "method", "mean_logs", "mean_crps")],
    by = c("scenario", "repetition", "method"),
    all = FALSE
  )
}

draw_simulation_recovery <- function() {
  plot_data <- load_simulation_results()
  layout_matrix <- rbind(
    1:4
  )
  graphics::layout(layout_matrix)
  graphics::par(
    oma = c(0.4, 0.5, 4.7, 0.5),
    family = "sans",
    xaxs = "r",
    yaxs = "r"
  )

  for (column_index in seq_along(scenario_order)) {
    draw_box_panel(
      plot_data,
      "weight_rmse",
      scenario_order[[column_index]],
      row_index = 1L,
      column_index = column_index,
      total_rows = 1L
    )
  }

  graphics::mtext(
    "Récupération des surfaces de poids selon le scénario simulé",
    side = 3,
    outer = TRUE,
    line = 3.0,
    font = 2,
    cex = 1.22
  )
  graphics::mtext(
    "RMSE calculée sur 30 répétitions indépendantes par scénario ; une valeur plus faible est meilleure.",
    side = 3,
    outer = TRUE,
    line = 1.55,
    cex = 0.82,
    col = "#4A4F55"
  )
}

draw_simulation_performance <- function() {
  plot_data <- load_simulation_results()
  graphics::layout(rbind(1:4, 5:8))
  graphics::par(
    oma = c(0.4, 0.5, 4.7, 0.5),
    family = "sans",
    xaxs = "r",
    yaxs = "r"
  )

  metrics <- c("mean_logs", "mean_crps")
  for (row_index in seq_along(metrics)) {
    for (column_index in seq_along(scenario_order)) {
      draw_box_panel(
        plot_data,
        metrics[[row_index]],
        scenario_order[[column_index]],
        row_index,
        column_index,
        total_rows = length(metrics)
      )
    }
  }

  graphics::mtext(
    "Performance prédictive selon le scénario simulé",
    side = 3,
    outer = TRUE,
    line = 3.0,
    font = 2,
    cex = 1.22
  )
  graphics::mtext(
    "LogS et CRPS sur 30 répétitions indépendantes par scénario ; une valeur plus faible est meilleure.",
    side = 3,
    outer = TRUE,
    line = 1.55,
    cex = 0.82,
    col = "#4A4F55"
  )
}

matrix_from_grid <- function(data, value, x_values, y_values) {
  keys <- paste(
    rep(x_values, times = length(y_values)),
    rep(y_values, each = length(x_values))
  )
  source_keys <- paste(data$horizon, data$age)
  matrix(
    data[[value]][match(keys, source_keys)],
    nrow = length(x_values),
    ncol = length(y_values)
  )
}

draw_lfo_heatmap <- function(
  best_lfo,
  title = "Meilleur modèle individuel dans la validation LFO"
) {
  graphics::par(
    mar = c(3.7, 4.4, 2.7, 0.8),
    las = 1,
    bty = "l",
    col.axis = "#30343B",
    col.lab = "#30343B",
    col.main = "#20242A",
    fg = "#30343B"
  )
  best_lfo$model_index <- match(best_lfo$best_model, model_order)
  matrix_value <- matrix_from_grid(
    best_lfo, "model_index", 1:10, 50:90
  )
  graphics::image(
    1:10,
    50:90,
    matrix_value,
    breaks = seq(0.5, 5.5, by = 1),
    col = unname(model_colors[model_order]),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  graphics::axis(1, at = c(1, 2, 4, 6, 8, 10))
  graphics::axis(2, at = c(50, 60, 70, 80, 90))
  graphics::mtext("Horizon", side = 1, line = 2.2)
  graphics::mtext("Âge", side = 2, line = 2.7, las = 0)
  if (nzchar(title)) {
    graphics::title(
      main = title,
      cex.main = 1.0,
      line = 0.8
    )
  }
  graphics::box()
}

draw_model_legend <- function() {
  graphics::par(mar = c(0.2, 0.2, 0.2, 0.2))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = model_order,
    fill = unname(model_colors[model_order]),
    border = "#30343B",
    bty = "n",
    cex = 0.9,
    title = "Modèle"
  )
}

draw_difference_heatmap <- function(
  data, model, limit, palette, row_index, column_index
) {
  part <- data[data$model == model, , drop = FALSE]
  matrix_value <- matrix_from_grid(part, "difference", 1:10, 50:90)
  graphics::par(
    mar = c(
      if (row_index == 2L) 3.4 else 1.2,
      if (column_index == 1L) 3.7 else 1.4,
      2.4,
      0.6
    ),
    las = 1,
    bty = "l",
    col.axis = "#30343B",
    col.lab = "#30343B",
    col.main = "#20242A",
    fg = "#30343B"
  )
  graphics::image(
    1:10,
    50:90,
    matrix_value,
    zlim = c(-limit, limit),
    col = palette,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  if (row_index == 2L) {
    graphics::axis(1, at = c(1, 4, 7, 10))
    graphics::mtext("Horizon", side = 1, line = 2.0, cex = 0.82)
  }
  if (column_index == 1L) {
    graphics::axis(2, at = c(50, 70, 90))
    graphics::mtext("Âge", side = 2, line = 2.3, las = 0, cex = 0.82)
  }
  graphics::title(main = model, cex.main = 0.95, line = 0.7)
  graphics::box()
}

draw_difference_scale <- function(limit, palette) {
  graphics::par(mar = rep(0, 4), las = 1)
  boundaries <- seq(-limit, limit, length.out = length(palette) + 1L)
  graphics::plot.new()
  padding <- max(0.12 * limit, 0.01)
  graphics::plot.window(
    xlim = c(-limit - padding, limit + padding),
    ylim = c(-0.45, 1.35),
    xaxs = "i",
    yaxs = "i"
  )
  graphics::rect(
    boundaries[-length(boundaries)],
    0.35,
    boundaries[-1L],
    0.72,
    col = palette,
    border = NA
  )
  graphics::rect(-limit, 0.35, limit, 0.72, border = "#30343B")
  tick_values <- c(-limit, 0, limit)
  graphics::segments(
    tick_values, 0.35, tick_values, 0.26,
    col = "#30343B"
  )
  graphics::text(
    tick_values,
    0.12,
    labels = formatC(tick_values, format = "f", digits = 2),
    cex = 0.72,
    col = "#30343B"
  )
  graphics::text(
    0,
    1.02,
    labels = "Échelle commune",
    font = 2,
    cex = 0.88,
    col = "#20242A"
  )
  graphics::text(
    0,
    -0.18,
    labels = "Hiérarchique - non régularisé",
    cex = 0.72,
    col = "#30343B"
  )
}

load_weight_differences <- function() {
  chapter_dir <- file.path(root, "results", "full", "chapter3_results")
  nonregularized <- read_csv(file.path(
    chapter_dir, "05_nonregularized_weights.csv"
  ))
  hierarchical <- read_csv(file.path(
    chapter_dir, "06_hierarchical_weights.csv"
  ))
  nonregularized$model <- toupper(nonregularized$model)
  hierarchical$model <- toupper(hierarchical$model)

  differences <- merge(
    hierarchical[c("age", "horizon", "model", "mean")],
    nonregularized[c("age", "horizon", "model", "weight")],
    by = c("age", "horizon", "model"),
    all = FALSE
  )
  differences$difference <- differences$mean - differences$weight
  differences
}

draw_empirical_heterogeneity <- function() {
  chapter_dir <- file.path(root, "results", "full", "chapter3_results")
  best_lfo <- read_csv(file.path(
    chapter_dir, "lfo_best_model_age_horizon.csv"
  ))

  graphics::layout(matrix(c(1, 1, 1, 1, 1, 2), nrow = 1))
  graphics::par(
    oma = c(0.4, 0.5, 4.7, 0.5),
    family = "sans"
  )

  draw_lfo_heatmap(best_lfo, title = "")
  draw_model_legend()

  graphics::mtext(
    "Hétérogénéité des performances individuelles selon l'âge et l'horizon",
    side = 3,
    outer = TRUE,
    line = 3.0,
    font = 2,
    cex = 1.22
  )
  graphics::mtext(
    "Meilleur modèle individuel dans la validation LFO, pour les âges 50-90 et les horizons 1-10.",
    side = 3,
    outer = TRUE,
    line = 1.55,
    cex = 0.82,
    col = "#4A4F55"
  )
}

draw_empirical_weights <- function() {
  differences <- load_weight_differences()
  limit <- max(abs(differences$difference))
  palette <- grDevices::hcl.colors(61, "Blue-Red 3", rev = TRUE)

  graphics::layout(rbind(1:3, 4:6))
  graphics::par(
    oma = c(0.4, 0.5, 4.7, 0.5),
    family = "sans"
  )
  draw_difference_heatmap(
    differences, "LC", limit, palette, 1L, 1L
  )
  draw_difference_heatmap(
    differences, "RH", limit, palette, 1L, 2L
  )
  draw_difference_heatmap(
    differences, "APC", limit, palette, 1L, 3L
  )
  draw_difference_heatmap(
    differences, "CBD", limit, palette, 2L, 1L
  )
  draw_difference_heatmap(
    differences, "M6", limit, palette, 2L, 2L
  )
  draw_difference_scale(limit, palette)

  graphics::mtext(
    "Poids de stacking selon l'âge et l'horizon",
    side = 3,
    outer = TRUE,
    line = 3.0,
    font = 2,
    cex = 1.22
  )
  graphics::mtext(
    "Écart hiérarchique - contextuel non régularisé, calculé cellule par cellule sur les âges 50-90 et les horizons 1-10.",
    side = 3,
    outer = TRUE,
    line = 1.55,
    cex = 0.83,
    col = "#4A4F55"
  )
}

metric_specs <- list(
  list(field = "logs", label = "LogS moyen", digits = 2),
  list(field = "crps", label = "CRPS moyen", digits = 4),
  list(field = "mae_deaths", label = "MAE des décès", digits = 0)
)

draw_line_panel <- function(
  data, x_field, metric_spec, row_index, column_index
) {
  part <- data[data$method %in% method_order, , drop = FALSE]
  x_values <- part[[x_field]]
  y_values <- part[[metric_spec$field]]
  y_padding <- diff(range(y_values)) * 0.08
  if (!is.finite(y_padding) || y_padding == 0) {
    y_padding <- max(abs(y_values), 1) * 0.05
  }
  y_limits <- range(y_values) + c(-y_padding, y_padding)

  graphics::par(
    mar = c(
      if (row_index == 3L) 3.6 else 1.2,
      if (column_index == 1L) 4.7 else 3.7,
      if (row_index == 1L) 2.7 else 1.0,
      0.8
    ),
    las = 1,
    bty = "l",
    col.axis = "#30343B",
    col.lab = "#30343B",
    col.main = "#20242A",
    fg = "#30343B"
  )
  graphics::plot(
    range(x_values),
    y_limits,
    type = "n",
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  draw_quiet_grid(horizontal = TRUE)

  for (method in method_order) {
    values <- part[part$method == method, , drop = FALSE]
    values <- values[order(values[[x_field]]), , drop = FALSE]
    graphics::lines(
      values[[x_field]],
      values[[metric_spec$field]],
      type = "o",
      pch = method_points[[method]],
      lty = method_lines[[method]],
      col = method_colors[[method]],
      lwd = 1.8,
      cex = if (x_field == "age") 0.42 else 0.62
    )
  }

  x_ticks <- if (x_field == "age") {
    c(50, 60, 70, 80, 90)
  } else {
    c(1, 2, 4, 6, 8, 10)
  }
  graphics::axis(1, at = x_ticks, labels = x_ticks, cex.axis = 0.82)
  y_ticks <- pretty(y_limits, n = 5)
  y_labels <- formatC(
    y_ticks,
    format = "f",
    digits = metric_spec$digits
  )
  graphics::axis(2, at = y_ticks, labels = y_labels, cex.axis = 0.78)
  graphics::box(bty = "l")

  if (row_index == 1L) {
    graphics::title(
      main = if (x_field == "age") "Selon l'âge" else "Selon l'horizon",
      cex.main = 1.02,
      line = 0.8
    )
  }
  if (column_index == 1L) {
    graphics::mtext(
      metric_spec$label,
      side = 2,
      line = 3.3,
      las = 0,
      font = 2,
      cex = 0.9
    )
  }
  if (row_index == 3L) {
    graphics::mtext(
      if (x_field == "age") "Âge" else "Horizon",
      side = 1,
      line = 2.3,
      cex = 0.86
    )
  }
}

draw_empirical_performance_fusion <- function() {
  chapter_dir <- file.path(root, "results", "full", "chapter3_results")
  metrics_age <- read_csv(file.path(chapter_dir, "09_metrics_by_age.csv"))
  metrics_horizon <- read_csv(file.path(
    chapter_dir, "10_metrics_by_horizon.csv"
  ))

  layout_matrix <- rbind(
    c(1, 2),
    c(3, 4),
    c(5, 6),
    c(7, 7)
  )
  graphics::layout(layout_matrix, heights = c(1, 1, 1, 0.18))
  graphics::par(oma = c(0.4, 0.5, 4.8, 0.5), family = "sans")

  panel_id <- 1L
  for (row_index in seq_along(metric_specs)) {
    draw_line_panel(
      metrics_age, "age", metric_specs[[row_index]], row_index, 1L
    )
    panel_id <- panel_id + 1L
    draw_line_panel(
      metrics_horizon,
      "horizon",
      metric_specs[[row_index]],
      row_index,
      2L
    )
    panel_id <- panel_id + 1L
  }

  graphics::par(mar = rep(0, 4))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = unname(method_labels[method_order]),
    col = unname(method_colors[method_order]),
    lty = unname(method_lines[method_order]),
    pch = unname(method_points[method_order]),
    pt.bg = unname(method_colors[method_order]),
    lwd = 1.8,
    horiz = TRUE,
    bty = "n",
    cex = 0.88,
    x.intersp = 0.8
  )

  graphics::mtext(
    "Performance prédictive des méthodes de stacking selon l'âge et l'horizon",
    side = 3,
    outer = TRUE,
    line = 3.0,
    font = 2,
    cex = 1.28
  )
  graphics::mtext(
    "Évaluation empirique sur la période de test ; une valeur plus faible est meilleure.",
    side = 3,
    outer = TRUE,
    line = 1.55,
    cex = 0.84,
    col = "#4A4F55"
  )
}

render_pair(
  "figure_02_recuperation_surfaces_poids",
  width = 14.5,
  height = 4.8,
  draw = draw_simulation_recovery
)
render_pair(
  "figure_12_performance_predictive_scenarios",
  width = 14.5,
  height = 7.2,
  draw = draw_simulation_performance
)
render_pair(
  "figure_10_heterogeneite_performances_individuelles",
  width = 12.5,
  height = 5.6,
  draw = draw_empirical_heterogeneity
)
render_pair(
  "figure_14_poids_stacking_age_horizon",
  width = 12.5,
  height = 8.4,
  draw = draw_empirical_weights
)
render_pair(
  "figures_15_17_performance_stacking_age_horizon",
  width = 12.5,
  height = 9.8,
  draw = draw_empirical_performance_fusion
)

message("Figures du mémoire générées dans : ", output_dir)
