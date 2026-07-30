#!/usr/bin/env Rscript

# Figure 6 du rapport technique de simulation :
# surfaces de poids vraies et estimees dans le scenario age-horizon.
#
# Ce script ne relance aucune simulation. Il lit uniquement les resultats
# deja calcules, puis produit une version de la figure avec :
#   - une meme echelle de couleur, fixee entre 0 et 1, pour tous les panneaux ;
#   - une barre de couleur commune et explicitement graduee ;
#   - des libelles francais accentues et suffisamment espaces.

options(stringsAsFactors = FALSE)

input_file <- file.path(
  "results", "full", "simulation", "main", "weights_by_repetition.csv"
)
output_dir <- "figures_memoire"
output_stem <- "figure_06_surfaces_poids_age_horizon_echelle_commune"

if (!file.exists(input_file)) {
  stop("Fichier de resultats introuvable : ", input_file)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Seules les huit premieres colonnes sont necessaires. Les colonnes
# d'intervalles sont ignorees afin de limiter la memoire utilisee.
weights <- utils::read.csv(
  input_file,
  colClasses = c(
    "character", "integer", "character", "numeric",
    "numeric", "character", "numeric", "numeric",
    "NULL", "NULL", "NULL", "NULL"
  ),
  check.names = FALSE
)

weights <- weights[weights$scenario == "age_horizon", ]
if (nrow(weights) == 0L) {
  stop("Aucun resultat trouve pour le scenario 'age_horizon'.")
}

number_of_repetitions <- length(unique(weights$repetition))
method_order <- c(
  "stacking_global",
  "stacking_contextual",
  "stacking_hierarchical"
)
model_order <- c("lc", "rh", "apc", "cbd", "m6")

expected_methods <- setdiff(method_order, unique(weights$method))
expected_models <- setdiff(model_order, unique(weights$model))
if (length(expected_methods) > 0L || length(expected_models) > 0L) {
  stop(
    "Resultats incomplets. Methodes absentes : ",
    paste(expected_methods, collapse = ", "),
    " ; modeles absents : ",
    paste(expected_models, collapse = ", ")
  )
}

# Moyenne de chaque surface sur les repetitions de l'etude principale.
surface_summary <- stats::aggregate(
  weights[c("true_weight", "estimated_weight")],
  by = weights[c("method", "age", "horizon", "model")],
  FUN = mean
)

row_methods <- c("true", method_order)
row_labels <- c("Poids vrais", "Global", "Non rég.", "Hiér.")
model_labels <- c(lc = "LC", rh = "RH", apc = "APC", cbd = "CBD", m6 = "M6")

# Palette sequentielle lisible, identique a celle du rapport original.
weight_palette <- grDevices::hcl.colors(60, "YlGnBu")
common_limits <- c(0, 1)

draw_common_scale <- function() {
  graphics::par(mar = c(0.2, 0.5, 1.3, 0.5))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

  left <- 0.20
  right <- 0.80
  bottom <- 0.38
  top <- 0.62
  edges <- seq(left, right, length.out = length(weight_palette) + 1L)

  for (index in seq_along(weight_palette)) {
    graphics::rect(
      edges[index], bottom, edges[index + 1L], top,
      col = weight_palette[index], border = NA
    )
  }
  graphics::rect(left, bottom, right, top, border = "#374151", lwd = 0.8)

  ticks <- seq(0, 1, by = 0.25)
  tick_positions <- left + (right - left) * ticks
  graphics::segments(
    tick_positions, bottom, tick_positions, bottom - 0.06,
    col = "#374151", lwd = 0.8
  )
  graphics::text(
    tick_positions, bottom - 0.11,
    labels = format(ticks, trim = TRUE, nsmall = 2),
    cex = 0.70, col = "#374151", adj = c(0.5, 1)
  )
  graphics::text(
    0.5, 0.90, "Échelle commune des poids",
    cex = 0.82, font = 2, col = "#111827"
  )
}

draw_figure <- function() {
  layout_matrix <- rbind(
    1:5,
    6:10,
    11:15,
    16:20,
    rep(21, 5)
  )
  graphics::layout(
    layout_matrix,
    heights = c(1, 1, 1, 1, 0.34)
  )
  graphics::par(
    oma = c(0.3, 2.4, 4.2, 0.5),
    family = "sans",
    fg = "#1F2937",
    col.axis = "#374151",
    col.lab = "#111827"
  )

  for (row_index in seq_along(row_methods)) {
    for (column_index in seq_along(model_order)) {
      model <- model_order[column_index]
      method <- row_methods[row_index]

      model_part <- surface_summary[surface_summary$model == model, ]
      if (method == "true") {
        # Les vrais poids sont identiques pour les trois methodes. Une seule
        # copie est donc gardee afin de ne pas compter trois fois chaque valeur.
        model_part <- model_part[
          model_part$method == "stacking_hierarchical",
        ]
        value_column <- "true_weight"
      } else {
        model_part <- model_part[model_part$method == method, ]
        value_column <- "estimated_weight"
      }

      matrix_values <- stats::xtabs(
        model_part[[value_column]] ~ model_part$age + model_part$horizon
      )
      age_values <- as.numeric(rownames(matrix_values))
      horizon_values <- as.numeric(colnames(matrix_values))

      graphics::par(
        mar = c(
          if (row_index == length(row_methods)) 2.7 else 0.65,
          if (column_index == 1L) 3.9 else 0.65,
          if (row_index == 1L) 2.1 else 0.65,
          0.65
        )
      )
      graphics::image(
        age_values,
        horizon_values,
        matrix_values,
        col = weight_palette,
        zlim = common_limits,
        xlab = "",
        ylab = "",
        axes = FALSE,
        useRaster = TRUE
      )

      if (row_index == 1L) {
        graphics::title(
          main = unname(model_labels[model]),
          cex.main = 0.92,
          font.main = 2,
          line = 0.55
        )
      }
      if (row_index == length(row_methods)) {
        graphics::axis(
          1,
          at = c(50, 70, 90),
          labels = c("50", "70", "90"),
          cex.axis = 0.69,
          tck = -0.025,
          mgp = c(1.5, 0.35, 0)
        )
        if (column_index == 3L) {
          graphics::mtext("Âge", side = 1, line = 1.55, cex = 0.78)
        }
      }
      if (column_index == 1L) {
        graphics::axis(
          2,
          at = c(1, 5, 10),
          labels = c("1", "5", "10"),
          cex.axis = 0.69,
          las = 1,
          tck = -0.025,
          mgp = c(1.5, 0.35, 0)
        )
        graphics::mtext(
          row_labels[row_index],
          side = 2,
          line = 2.55,
          cex = 0.77,
          font = 2
        )
      }
      graphics::box(col = "#374151", lwd = 0.75)
    }
  }

  draw_common_scale()

  graphics::mtext(
    "Horizon",
    side = 2,
    outer = TRUE,
    line = 1.0,
    cex = 0.84,
    font = 2
  )
  graphics::mtext(
    "Poids selon l’âge et l’horizon",
    side = 3,
    outer = TRUE,
    line = 2.35,
    cex = 1.22,
    font = 2,
    col = "#111827"
  )
  graphics::mtext(
    sprintf(
      paste0(
        "Surfaces vraies et estimations moyennes sur %d répétitions ; ",
        "tous les panneaux utilisent l’échelle 0–1."
      ),
      number_of_repetitions
    ),
    side = 3,
    outer = TRUE,
    line = 0.95,
    cex = 0.76,
    col = "#4B5563"
  )
}

pdf_file <- file.path(output_dir, paste0(output_stem, ".pdf"))
png_file <- file.path(output_dir, paste0(output_stem, ".png"))

grDevices::cairo_pdf(
  filename = pdf_file,
  width = 13.5,
  height = 9.5,
  family = "sans",
  bg = "white"
)
draw_figure()
grDevices::dev.off()

grDevices::png(
  filename = png_file,
  width = 2700,
  height = 1900,
  res = 200,
  type = "cairo",
  bg = "white"
)
draw_figure()
grDevices::dev.off()

message("Figure PDF : ", normalizePath(pdf_file, winslash = "/"))
message("Figure PNG : ", normalizePath(png_file, winslash = "/"))
