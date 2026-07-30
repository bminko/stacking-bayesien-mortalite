#!/usr/bin/env Rscript

# Figure 2 du rapport empirique :
# heatmap des log-taux bruts de mortalite observes entre 1970 et 2024.
#
# Le script ne relance aucun modele. Il lit les donnees deja preparees et
# ajoute une barre de couleur verticale couvrant exactement l'etendue
# observee. La palette HCL est sequentielle et perceptuellement uniforme :
# les couleurs claires representent les taux faibles et les couleurs foncees
# les taux eleves.

options(stringsAsFactors = FALSE)

input_file <- file.path(
  "data", "processed", "full", "mortality_data.rds"
)
output_dir <- "figures_memoire"
output_stem <- "figure_02_empirique_log_taux_mortalite_echelle_couleur"

if (!file.exists(input_file)) {
  stop("Fichier de donnees introuvable : ", input_file)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

processed <- readRDS(input_file)
if (is.null(processed$long)) {
  stop("L'objet de donnees ne contient pas la table 'long'.")
}

mortality <- processed$long
required_columns <- c("year", "age", "deaths", "exposure")
missing_columns <- setdiff(required_columns, names(mortality))
if (length(missing_columns) > 0L) {
  stop(
    "Colonnes absentes : ",
    paste(missing_columns, collapse = ", ")
  )
}
if (any(!is.finite(mortality$exposure)) || any(mortality$exposure <= 0)) {
  stop("Les expositions doivent etre finies et strictement positives.")
}

mortality$crude_rate <- mortality$deaths / mortality$exposure
if (any(!is.finite(mortality$crude_rate)) ||
    any(mortality$crude_rate <= 0)) {
  stop(
    "Les taux bruts doivent etre finis et strictement positifs ",
    "pour calculer leur logarithme."
  )
}
mortality$log_crude_rate <- log(mortality$crude_rate)

years <- sort(unique(mortality$year))
ages <- sort(unique(mortality$age))
log_rate_matrix <- stats::xtabs(
  log_crude_rate ~ age + year,
  data = mortality
)

if (nrow(log_rate_matrix) != length(ages) ||
    ncol(log_rate_matrix) != length(years)) {
  stop("La grille age-annee des taux de mortalite est incomplete.")
}

observed_limits <- range(mortality$log_crude_rate)

# Les bornes exactes sont affichees, avec des valeurs entieres intermediaires.
# Un espacement minimal evite deux etiquettes presque superposees.
integer_ticks <- seq(
  ceiling(observed_limits[1]),
  floor(observed_limits[2]),
  by = 1
)
integer_ticks <- integer_ticks[
  integer_ticks > observed_limits[1] + 0.30 &
    integer_ticks < observed_limits[2] - 0.30
]
legend_ticks <- c(
  observed_limits[1],
  integer_ticks,
  observed_limits[2]
)
legend_labels <- ifelse(
  legend_ticks %in% integer_ticks,
  as.character(as.integer(round(legend_ticks))),
  formatC(
    legend_ticks,
    format = "f",
    digits = 2,
    decimal.mark = ","
  )
)

# YlOrRd est une palette HCL a luminance ordonnee. Avec rev = TRUE,
# les petites valeurs sont jaune clair et les grandes valeurs rouge fonce.
mortality_palette <- grDevices::hcl.colors(
  160,
  palette = "YlOrRd",
  rev = TRUE
)

draw_color_bar <- function() {
  value_edges <- seq(
    observed_limits[1],
    observed_limits[2],
    length.out = length(mortality_palette) + 1L
  )
  value_padding <- diff(observed_limits) * 0.12

  graphics::par(mar = c(4.6, 0.3, 3.2, 0.3))
  graphics::plot.new()
  graphics::plot.window(
    xlim = c(0, 1.95),
    ylim = c(
      observed_limits[1] - value_padding,
      observed_limits[2] + value_padding
    ),
    xaxs = "i",
    yaxs = "i"
  )

  bar_left <- 0.08
  bar_right <- 0.42
  for (index in seq_along(mortality_palette)) {
    graphics::rect(
      bar_left,
      value_edges[index],
      bar_right,
      value_edges[index + 1L],
      col = mortality_palette[index],
      border = NA
    )
  }
  graphics::rect(
    bar_left,
    observed_limits[1],
    bar_right,
    observed_limits[2],
    border = "#374151",
    lwd = 0.8
  )

  graphics::segments(
    bar_right,
    legend_ticks,
    bar_right + 0.08,
    legend_ticks,
    col = "#374151",
    lwd = 0.8
  )
  graphics::text(
    bar_right + 0.13,
    legend_ticks,
    labels = legend_labels,
    adj = c(0, 0.5),
    cex = 0.72,
    col = "#374151"
  )

  graphics::text(
    0.88,
    observed_limits[2],
    labels = "Taux les plus élevés",
    adj = c(0, 0.5),
    cex = 0.70,
    font = 2,
    col = "#7D0025"
  )
  graphics::text(
    0.88,
    observed_limits[1],
    labels = "Taux les plus faibles",
    adj = c(0, 0.5),
    cex = 0.70,
    font = 2,
    col = "#6B5A00"
  )
  graphics::title(
    main = "Log-taux brut\nde mortalité",
    cex.main = 0.88,
    font.main = 2,
    line = 0.55,
    col.main = "#111827"
  )
}

draw_figure <- function() {
  graphics::layout(
    matrix(c(1, 2), nrow = 1L),
    widths = c(5.6, 1.75)
  )
  graphics::par(
    oma = c(0.4, 0.4, 4.4, 0.4),
    family = "sans",
    fg = "#30343B",
    col.axis = "#30343B",
    col.lab = "#30343B",
    col.main = "#111827"
  )

  graphics::par(mar = c(4.6, 4.8, 1.0, 0.5), las = 1)
  graphics::image(
    years,
    ages,
    t(log_rate_matrix),
    col = mortality_palette,
    zlim = observed_limits,
    xlab = "Année",
    ylab = "",
    axes = FALSE,
    useRaster = TRUE
  )
  graphics::axis(
    1,
    at = c(1970, 1980, 1990, 2000, 2010, 2020),
    cex.axis = 0.82,
    tck = -0.018,
    mgp = c(2.5, 0.65, 0)
  )
  graphics::axis(
    2,
    at = c(50, 60, 70, 80, 90),
    cex.axis = 0.82,
    tck = -0.018,
    mgp = c(2.5, 0.65, 0)
  )
  graphics::mtext(
    "Âge",
    side = 2,
    line = 3.0,
    cex = 0.90,
    las = 0
  )
  graphics::box(col = "#374151", lwd = 0.8)

  draw_color_bar()

  graphics::mtext(
    "Log-taux bruts de mortalité, 1970-2024",
    side = 3,
    outer = TRUE,
    line = 2.45,
    cex = 1.22,
    font = 2,
    col = "#111827"
  )
  graphics::mtext(
    paste0(
      "Population totale, âges 50 à 90 ans ; ",
      "le rouge foncé correspond aux taux les plus élevés."
    ),
    side = 3,
    outer = TRUE,
    line = 1.05,
    cex = 0.76,
    col = "#4B5563"
  )
}

pdf_file <- file.path(output_dir, paste0(output_stem, ".pdf"))
png_file <- file.path(output_dir, paste0(output_stem, ".png"))

grDevices::cairo_pdf(
  filename = pdf_file,
  width = 12.2,
  height = 7.6,
  family = "sans",
  bg = "white"
)
draw_figure()
grDevices::dev.off()

grDevices::png(
  filename = png_file,
  width = 2440,
  height = 1520,
  res = 200,
  type = "cairo",
  bg = "white"
)
draw_figure()
grDevices::dev.off()

message("Etendue observee : ", paste(round(observed_limits, 4), collapse = " a "))
message("Figure PDF : ", normalizePath(pdf_file, winslash = "/"))
message("Figure PNG : ", normalizePath(png_file, winslash = "/"))
