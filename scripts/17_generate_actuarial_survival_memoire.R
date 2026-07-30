#!/usr/bin/env Rscript

# Regénère uniquement la figure de survie utilisée dans le mémoire.
# Les trajectoires et leurs quantiles sont lus dans le CSV final déjà validé :
# aucun modèle Stan et aucune simulation actuarielle ne sont relancés.

input_file <- file.path(
  "results", "full", "chapter3_results",
  "46_actuarial_survival_all_rules.csv"
)
output_dir <- "figures_memoire"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

survival <- utils::read.csv(input_file, stringsAsFactors = FALSE)
survival <- survival[
  survival$analysis == "primary_h_le_10" &
    survival$rule == "within_lfo",
]

series <- list(
  global = list(
    method = "stacking_global",
    weight_mode = "fixed_weights",
    label = "Global",
    color = "#3B6FB6",
    lty = 1,
    pch = 16
  ),
  contextual = list(
    method = "stacking_contextual",
    weight_mode = "fixed_weights",
    label = "Contextuel non régularisé",
    color = "#E97824",
    lty = 2,
    pch = 17
  ),
  hierarchical_mean = list(
    method = "stacking_hierarchical",
    weight_mode = "posterior_mean_weights",
    label = "Hiérarchique, poids moyens",
    color = "#238A3B",
    lty = 3,
    pch = 15
  ),
  hierarchical_propagated = list(
    method = "stacking_hierarchical",
    weight_mode = "propagated_weights",
    label = "Hiérarchique, poids propagés",
    color = "#8A3D8C",
    lty = 4,
    pch = 18
  )
)

initial_ages <- sort(unique(survival$initial_age))

draw_figure <- function() {
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(2, 3),
    mar = c(4.2, 4.4, 2.8, 1.0),
    oma = c(1.4, 0.4, 0.6, 0.2),
    mgp = c(2.6, 0.8, 0),
    tcl = -0.25,
    las = 1,
    cex.axis = 0.90,
    cex.lab = 1.00,
    cex.main = 1.10,
    family = "sans"
  )

  for (age_index in seq_along(initial_ages)) {
    initial_age <- initial_ages[[age_index]]
    age_data <- survival[survival$initial_age == initial_age,]
    x_range <- range(age_data$horizon)
    y_range <- range(age_data$q025, age_data$q975)
    padding <- 0.04 * diff(y_range)
    if (!is.finite(padding) || padding == 0) {
      padding <- 0.001
    }

    graphics::plot(
      x_range,
      y_range + c(-padding, padding),
      type = "n",
      xlab = "Horizon (années)",
      ylab = "Probabilité de survie",
      main = paste("Âge initial", initial_age),
      axes = FALSE
    )
    graphics::axis(1, at = c(2, 4, 6, 8, 10))
    graphics::axis(2)
    graphics::box(bty = "l")

    for (definition in series) {
      part <- age_data[
        age_data$method == definition$method &
          age_data$weight_mode == definition$weight_mode,
      ]
      part <- part[order(part$horizon),]
      graphics::polygon(
        c(part$horizon, rev(part$horizon)),
        c(part$q025, rev(part$q975)),
        col = grDevices::adjustcolor(definition$color, alpha.f = 0.09),
        border = NA
      )
    }

    for (definition in series) {
      part <- age_data[
        age_data$method == definition$method &
          age_data$weight_mode == definition$weight_mode,
      ]
      part <- part[order(part$horizon),]
      graphics::lines(
        part$horizon,
        part$mean,
        col = definition$color,
        lty = definition$lty,
        lwd = 1.8,
        type = "b",
        pch = definition$pch,
        cex = 0.45
      )
    }

    if (age_index == 1L) {
      graphics::legend(
        "bottomleft",
        legend = vapply(series, `[[`, character(1), "label"),
        col = vapply(series, `[[`, character(1), "color"),
        lty = vapply(series, `[[`, numeric(1), "lty"),
        pch = vapply(series, `[[`, numeric(1), "pch"),
        lwd = 1.8,
        pt.cex = 0.65,
        cex = 0.68,
        bty = "n"
      )
    }
  }
}

pdf_file <- file.path(output_dir, "actuarial_primary_survival.pdf")
grDevices::pdf(pdf_file, width = 12, height = 8, useDingbats = FALSE)
draw_figure()
grDevices::dev.off()

png_file <- file.path(output_dir, "actuarial_primary_survival.png")
grDevices::png(
  png_file,
  width = 3600,
  height = 2400,
  res = 300,
  type = "cairo"
)
draw_figure()
grDevices::dev.off()

message("Figure de survie régénérée : ", pdf_file, " et ", png_file)
