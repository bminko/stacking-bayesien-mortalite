# Graphiques sobres en base R, directement reutilisables dans le memoire.

plot_metric_by_horizon <- function(by_horizon, metric, path,
                                   ylab = metric) {
  methods <- unique(by_horizon$method)
  horizons <- sort(unique(by_horizon$horizon))
  colors <- grDevices::hcl.colors(length(methods), "Dark 3")
  grDevices::png(path, width = 1800, height = 1100, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  y_range <- range(by_horizon[[metric]], finite = TRUE)
  graphics::plot(
    horizons, rep(NA_real_, length(horizons)),
    type = "n", ylim = y_range,
    xlab = "Horizon (annees)", ylab = ylab,
    main = paste(ylab, "selon l'horizon")
  )
  for (i in seq_along(methods)) {
    part <- by_horizon[by_horizon$method == methods[i], ]
    part <- part[order(part$horizon), ]
    graphics::lines(
      part$horizon, part[[metric]],
      type = "b", pch = 16, lwd = 2, col = colors[i]
    )
  }
  graphics::legend(
    "topright", legend = methods, col = colors,
    lty = 1, pch = 16, cex = 0.75, bg = "white"
  )
  invisible(path)
}

plot_metric_heatmap <- function(cells, method, metric, path) {
  part <- cells[cells$method == method, ]
  aggregated <- stats::aggregate(
    part[[metric]],
    by = list(age = part$age, horizon = part$horizon),
    FUN = mean
  )
  names(aggregated)[3L] <- "value"
  ages <- sort(unique(aggregated$age))
  horizons <- sort(unique(aggregated$horizon))
  matrix_values <- matrix(
    aggregated$value[
      match(
        paste(rep(ages, times = length(horizons)),
              rep(horizons, each = length(ages))),
        paste(aggregated$age, aggregated$horizon)
      )
    ],
    nrow = length(ages),
    ncol = length(horizons)
  )

  grDevices::png(path, width = 1400, height = 1200, res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::image(
    x = horizons,
    y = ages,
    z = t(matrix_values),
    col = grDevices::hcl.colors(30, "YlOrRd", rev = TRUE),
    xlab = "Horizon (annees)",
    ylab = "Age",
    main = paste(metric, "-", method)
  )
  invisible(path)
}

plot_contextual_weights <- function(weight_long, path,
                                    value = "weight_mean") {
  models <- unique(weight_long$model)
  horizons <- sort(unique(weight_long$horizon))
  ages <- sort(unique(weight_long$age))
  colors <- grDevices::hcl.colors(length(models), "Dark 3")
  selected_horizons <- unique(round(stats::quantile(
    horizons, c(0, 0.5, 1), names = FALSE
  )))

  grDevices::png(path, width = 1800, height = 600 * length(selected_horizons),
                 res = 180)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(length(selected_horizons), 1L), mar = c(4, 4, 3, 1))
  for (horizon in selected_horizons) {
    graphics::plot(
      range(ages), c(0, 1), type = "n",
      xlab = "Age", ylab = "Poids",
      main = paste("Horizon", horizon)
    )
    for (k in seq_along(models)) {
      part <- weight_long[
        weight_long$model == models[k] &
          weight_long$horizon == horizon,
      ]
      part <- part[order(part$age), ]
      graphics::lines(
        part$age, part[[value]], lwd = 2.5, col = colors[k]
      )
    }
    graphics::legend(
      "topright", legend = models, col = colors, lty = 1,
      lwd = 2.5, cex = 0.8, bg = "white"
    )
  }
  invisible(path)
}
