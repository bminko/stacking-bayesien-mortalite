#!/usr/bin/env Rscript

# Etape 8 : tableaux et figures directement utilisables au chapitre 3.
#
# Cette etape lit uniquement les ajustements et poids deja sauvegardes.
# Elle ne lance aucun ajustement Stan et ne modifie pas le meta-jeu LFO.

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
source(file.path("R", "model_fitting.R"))
source(file.path("R", "forecasting.R"))
source(file.path("R", "aggregation.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "reporting.R"))
cfg <- load_config()
# Les tirages deja archives sont relus par la voie de secours documentee
# dans R/forecasting.R : aucune recompilation Stan n'est necessaire ici.

results_dir <- file.path(cfg$paths$results, "chapter3_results")
figures_dir <- file.path(cfg$root, "output", "pdf", "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

processed <- readRDS(file.path(
  cfg$paths$processed, "mortality_data.rds"
))
meta_object <- readRDS(file.path(
  cfg$paths$processed, "lfo_meta_data.rds"
))
meta <- meta_object$data
aggregation <- readRDS(file.path(
  cfg$paths$weights, "aggregation_results.rds"
))
evaluation <- readRDS(file.path(
  cfg$paths$metrics, "test_evaluation.rds"
))
training <- subset_training_data(processed, cfg$validation_end)
method_labels <- chapter_method_labels()
colors <- chapter_model_colors()
primary_methods <- names(method_labels)
primary_cells <- evaluation$cells[
  evaluation$cells$method %in% primary_methods,
  ,
  drop = FALSE
]

write_result <- function(data, filename) {
  write_csv_atomic(data, file.path(results_dir, filename))
}

# 01 - Protocole.
protocol <- data.frame(
  element = c(
    "population", "ages", "periode_disponible",
    "entrainement_principal", "validation_lfo", "test_initial",
    "extension", "test_etendu", "origines_lfo", "horizons_lfo",
    "horizon_test_maximal", "tirages_predictifs", "chaines_stan",
    "warmup_par_chaine", "iterations_conservees_par_chaine",
    "iterations_conservees_totales", "graine", "modeles",
    "origine_test", "type_test", "convention_logs"
  ),
  value = c(
    cfg$sex,
    paste(cfg$age_min, cfg$age_max, sep = "-"),
    paste(cfg$year_start, cfg$data_end, sep = "-"),
    paste(cfg$year_start, cfg$validation_end, sep = "-"),
    "2001-2015",
    "2016-2020",
    "2021-2024",
    "2016-2024",
    paste(cfg$lfo_origins, collapse = ","),
    paste(cfg$lfo_horizons, collapse = ","),
    cfg$test_end - cfg$validation_end,
    cfg$forecast_draws,
    cfg$mcmc$chains,
    cfg$mcmc$iter_warmup,
    cfg$mcmc$iter_sampling,
    cfg$mcmc$chains * cfg$mcmc$iter_sampling,
    cfg$seed,
    paste(toupper(cfg$models), collapse = ","),
    cfg$validation_end,
    "origine fixe",
    "negative log-densite predictive moyenne, a minimiser"
  ),
  stringsAsFactors = FALSE
)
write_result(protocol, "01_protocol.csv")

# 02 - Resume des donnees.
data_long <- processed$long
data_long$crude_rate <- data_long$deaths / data_long$exposure
data_summary <- data.frame(
  population = cfg$sex,
  first_year = min(data_long$year),
  last_year = max(data_long$year),
  years = length(unique(data_long$year)),
  min_age = min(data_long$age),
  max_age = max(data_long$age),
  ages = length(unique(data_long$age)),
  cells = nrow(data_long),
  deaths_min = min(data_long$deaths),
  deaths_mean = mean(data_long$deaths),
  deaths_max = max(data_long$deaths),
  exposure_min = min(data_long$exposure),
  exposure_mean = mean(data_long$exposure),
  exposure_max = max(data_long$exposure),
  crude_rate_min = min(data_long$crude_rate),
  crude_rate_mean = mean(data_long$crude_rate),
  crude_rate_max = max(data_long$crude_rate),
  missing_deaths = sum(is.na(data_long$deaths)),
  missing_exposure = sum(is.na(data_long$exposure)),
  nonpositive_exposure = sum(data_long$exposure <= 0),
  negative_deaths = sum(data_long$deaths < 0),
  fractional_deaths = sum(
    abs(data_long$deaths - round(data_long$deaths)) > 1e-9
  ),
  duplicate_cells = sum(duplicated(data_long[c("year", "age")])),
  complete_grid = nrow(data_long) ==
    length(unique(data_long$year)) * length(unique(data_long$age)),
  stringsAsFactors = FALSE
)
write_result(data_summary, "02_data_summary.csv")

# 03 - Diagnostics des cinq ajustements utilises pour le test.
parse_clock <- function(line) {
  value <- sub("^\\[([0-9]{2}):([0-9]{2}):([0-9]{2})\\].*$",
               "\\1:\\2:\\3", line)
  parts <- as.integer(strsplit(value, ":", fixed = TRUE)[[1L]])
  sum(parts * c(3600L, 60L, 1L))
}
pipeline_log <- file.path(cfg$paths$results, "pipeline_stderr.log")
log_lines <- if (file.exists(pipeline_log)) {
  readLines(pipeline_log, warn = FALSE, encoding = "UTF-8")
} else {
  character()
}
model_starts <- setNames(rep(NA_integer_, length(cfg$models)), cfg$models)
model_ends <- setNames(rep(NA_integer_, length(cfg$models)), cfg$models)
for (index in seq_along(cfg$models)) {
  model <- cfg$models[[index]]
  start_pattern <- paste0(
    "test_training_", cfg$validation_end, "/", model,
    " - ajustement Stan, essai 1"
  )
  start_line <- grep(start_pattern, log_lines, fixed = TRUE)
  if (length(start_line)) {
    model_starts[[model]] <- parse_clock(log_lines[start_line[[1L]]])
    next_pattern <- if (index < length(cfg$models)) {
      paste0(
        "Compilation/chargement du modele ",
        toupper(cfg$models[[index + 1L]])
      )
    } else {
      "Lancement de scripts/03_lfo_validation_loop.R"
    }
    end_line <- grep(next_pattern, log_lines, fixed = TRUE)
    end_line <- end_line[end_line > start_line[[1L]]]
    if (length(end_line)) {
      model_ends[[model]] <- parse_clock(log_lines[end_line[[1L]]])
    }
  }
}

diagnostic_rows <- lapply(cfg$models, function(model) {
  fit_dir <- file.path(
    cfg$paths$fits,
    paste0("test_training_", cfg$validation_end),
    model
  )
  overview <- utils::read.csv(
    file.path(fit_dir, "diagnostics_overview.csv"),
    stringsAsFactors = FALSE
  )
  sampler <- utils::read.csv(
    file.path(fit_dir, "diagnostics_sampler.csv"),
    stringsAsFactors = FALSE
  )
  elapsed <- model_ends[[model]] - model_starts[[model]]
  if (is.finite(elapsed) && elapsed < 0) elapsed <- elapsed + 24 * 3600
  data.frame(
    model = toupper(model),
    elapsed_seconds = elapsed,
    elapsed_minutes = elapsed / 60,
    chains = cfg$mcmc$chains,
    warmup_per_chain = overview$rhat_limit * 0 + cfg$mcmc$iter_warmup,
    kept_per_chain = cfg$mcmc$iter_sampling,
    max_rhat = overview$max_rhat,
    min_ess_bulk = overview$min_ess_bulk,
    min_ess_tail = overview$min_ess_tail,
    divergences = overview$divergences,
    max_treedepth_hits = overview$max_treedepth_hits,
    min_ebfmi = min(sampler$ebfmi),
    status = if (isTRUE(overview$pass)) "valide" else "non valide",
    stringsAsFactors = FALSE
  )
})
model_diagnostics <- do.call(rbind, diagnostic_rows)
write_result(model_diagnostics, "03_model_diagnostics.csv")

# 04 - Stacking global.
global_weights <- data.frame(
  model = toupper(cfg$models),
  weight = as.numeric(aggregation$global$weights[cfg$models]),
  criterion = aggregation$global$optimization$value,
  convergence = aggregation$global$optimization$convergence,
  multistarts = length(aggregation$global$all_values),
  objective_range = diff(range(aggregation$global$all_values)),
  sum_weights = sum(aggregation$global$weights),
  stringsAsFactors = FALSE
)
write_result(global_weights, "04_global_weights.csv")

# 05 et 06 - Poids contextuels aux horizons du test.
context_grid <- expand.grid(
  age = cfg$age_min:cfg$age_max,
  horizon = 1:(cfg$test_end - cfg$validation_end)
)
nonregularized_grid <- contextual_weight_grid(
  context_grid,
  aggregation$contextual,
  cfg$models,
  aggregation$standardization
)
nonregularized_weights <- do.call(rbind, lapply(cfg$models, function(model) {
  data.frame(
    age = nonregularized_grid$age,
    horizon = nonregularized_grid$horizon,
    model = toupper(model),
    weight = nonregularized_grid[[model]],
    stringsAsFactors = FALSE
  )
}))
write_result(nonregularized_weights, "05_nonregularized_weights.csv")

hierarchical_fit <- readRDS(aggregation$hierarchical_fit_path)
hierarchical_weights <- hierarchical_weight_summary(
  context_grid,
  hierarchical_fit,
  cfg$models,
  aggregation$standardization,
  ndraws = cfg$forecast_draws,
  seed = stable_seed(cfg$seed, "chapter", "hierarchical_weights")
)
hierarchical_weights$model <- toupper(hierarchical_weights$model)
write_result(hierarchical_weights, "06_hierarchical_weights.csv")

# 07 - Coefficients des deux stackings contextuels.
contextual_coefficients <- as.data.frame(
  as.table(aggregation$contextual$coefficients),
  stringsAsFactors = FALSE
)
names(contextual_coefficients) <- c("model", "term", "mean")
contextual_coefficients$method <- "stacking_contextual"
contextual_coefficients$median <- contextual_coefficients$mean
contextual_coefficients$sd <- NA_real_
contextual_coefficients$q025 <- NA_real_
contextual_coefficients$q975 <- NA_real_
contextual_coefficients$model <- toupper(contextual_coefficients$model)
contextual_coefficients <- contextual_coefficients[c(
  "method", "model", "term", "mean", "median", "sd", "q025", "q975"
)]
hierarchical_coefficients <- hierarchical_coefficient_summary(
  hierarchical_fit,
  cfg$models,
  ndraws = cfg$forecast_draws,
  seed = stable_seed(cfg$seed, "chapter", "hierarchical_coefficients")
)
hierarchical_coefficients$model <- toupper(
  hierarchical_coefficients$model
)
weight_coefficients <- rbind(
  contextual_coefficients,
  hierarchical_coefficients
)
write_result(weight_coefficients, "07_weight_coefficients.csv")

# 08 a 11 - Metriques du test etendu.
period_summaries <- summarize_evaluation_periods(
  primary_cells,
  list(test_extended_2016_2024 = 2016:2024)
)
global_metrics <- period_summaries$overall
global_metrics$method_label <- unname(method_labels[global_metrics$method])
global_metrics <- global_metrics[c(
  "method", "method_label", "logs", "crps", "mae_deaths",
  "coverage80", "width80", "coverage95", "width95"
)]
write_result(global_metrics, "08_global_metrics.csv")

metrics_by_age <- period_summaries$by_age
metrics_by_age$method_label <- unname(
  method_labels[metrics_by_age$method]
)
write_result(metrics_by_age, "09_metrics_by_age.csv")

metrics_by_horizon <- period_summaries$by_horizon
metrics_by_horizon$method_label <- unname(
  method_labels[metrics_by_horizon$method]
)
write_result(metrics_by_horizon, "10_metrics_by_horizon.csv")

interval_metrics <- global_metrics[c(
  "method", "method_label",
  "coverage80", "width80", "coverage95", "width95"
)]
write_result(interval_metrics, "11_interval_metrics.csv")

# 12 et 13 - Quantites actuarielles deja calculees depuis 2015.
survival <- utils::read.csv(
  file.path(cfg$paths$metrics, "survival_probabilities.csv"),
  stringsAsFactors = FALSE
)
annuity <- utils::read.csv(
  file.path(cfg$paths$metrics, "annuity_values.csv"),
  stringsAsFactors = FALSE
)
survival$method_label <- unname(method_labels[survival$method])
annuity$method_label <- unname(method_labels[annuity$method])
write_result(survival, "12_survival_probabilities.csv")
write_result(annuity, "13_annuity_values.csv")

# 14 - Sensibilite temporelle et bootstrap par annee cible.
sensitivity <- utils::read.csv(
  file.path(cfg$paths$metrics, "sensitivity_metrics_overall.csv"),
  stringsAsFactors = FALSE
)
sensitivity$method_label <- unname(method_labels[sensitivity$method])
write_result(sensitivity, "14_sensitivity_results.csv")

comparisons <- list(
  contextual_minus_global = c(
    "stacking_contextual", "stacking_global"
  ),
  hierarchical_minus_global = c(
    "stacking_hierarchical", "stacking_global"
  ),
  hierarchical_minus_contextual = c(
    "stacking_hierarchical", "stacking_contextual"
  )
)
bootstrap_differences <- bootstrap_performance_differences(
  primary_cells,
  comparisons,
  repetitions = 2000L,
  seed = stable_seed(cfg$seed, "chapter", "bootstrap_year")
)
write_result(
  bootstrap_differences,
  "performance_differences_bootstrap.csv"
)

# Meta-jeu LFO en format long et resumes.
meta_long <- do.call(rbind, lapply(cfg$models, function(model) {
  data.frame(
    age = meta$age,
    origin = meta$origin,
    horizon = meta$horizon,
    target_year = meta$year,
    model = toupper(model),
    observed_deaths = meta$deaths,
    exposure = meta$exposure,
    log_pred_density = meta[[paste0("log_p_", model)]],
    multiplicity = 1 / meta$omega,
    validation_weight = meta$omega,
    age_z = meta$x_tilde,
    horizon_z = meta$h_tilde,
    stringsAsFactors = FALSE
  )
}))
write_result(meta_long, "meta_validation.csv")

lfo_model_summary <- do.call(rbind, lapply(cfg$models, function(model) {
  logp <- meta[[paste0("log_p_", model)]]
  data.frame(
    model = toupper(model),
    weighted_mean_log_density = weighted_column_mean(logp, meta$omega),
    unweighted_mean_log_density = mean(logp),
    minimum = min(logp),
    maximum = max(logp),
    stringsAsFactors = FALSE
  )
}))
write_result(lfo_model_summary, "lfo_model_summary.csv")

best_rows <- list()
position <- 1L
for (age in sort(unique(meta$age))) {
  for (horizon in sort(unique(meta$horizon))) {
    part <- meta[meta$age == age & meta$horizon == horizon, ]
    means <- vapply(cfg$models, function(model) {
      weighted_column_mean(
        part[[paste0("log_p_", model)]],
        part$omega
      )
    }, numeric(1))
    best_rows[[position]] <- data.frame(
      age = age,
      horizon = horizon,
      best_model = toupper(cfg$models[[which.max(means)]]),
      best_mean_log_density = max(means),
      stringsAsFactors = FALSE
    )
    position <- position + 1L
  }
}
best_lfo <- do.call(rbind, best_rows)
write_result(best_lfo, "lfo_best_model_age_horizon.csv")

# Comparaison directe des stackings sur le critere de validation.
log_p <- meta_logp_matrix(meta, cfg$models)
global_matrix <- constant_weight_matrix(
  aggregation$global$weights,
  nrow(meta),
  cfg$models
)
context_matrix <- as.matrix(contextual_weight_grid(
  meta[c("age", "horizon")],
  aggregation$contextual,
  cfg$models,
  aggregation$standardization
)[cfg$models])
hierarchical_matrix <- weight_matrix_from_long(
  aggregation$hierarchical_weights,
  meta[c("age", "horizon")],
  cfg$models
)
validation_criterion <- function(weights) {
  -sum(meta$omega * row_log_sum_exp(log(weights) + log_p))
}
stacking_comparison <- data.frame(
  method = c(
    "stacking_global",
    "stacking_contextual",
    "stacking_hierarchical"
  ),
  validation_criterion = c(
    validation_criterion(global_matrix),
    validation_criterion(context_matrix),
    validation_criterion(hierarchical_matrix)
  ),
  parameter_count = c(
    length(cfg$models) - 1L,
    (length(cfg$models) - 1L) * 5L,
    (length(cfg$models) - 1L) * 5L + 4L
  ),
  minimum_weight = c(
    min(global_matrix), min(context_matrix), min(hierarchical_matrix)
  ),
  maximum_weight = c(
    max(global_matrix), max(context_matrix), max(hierarchical_matrix)
  ),
  stability = c(
    sprintf(
      "%d initialisations; amplitude critere %.3g",
      length(aggregation$global$all_values),
      diff(range(aggregation$global$all_values))
    ),
    sprintf(
      "%d initialisations; amplitude critere %.3g",
      length(aggregation$contextual$all_values),
      diff(range(aggregation$contextual$all_values))
    ),
    "regularisation hierarchique; diagnostics Stan valides"
  ),
  stringsAsFactors = FALSE
)
write_result(stacking_comparison, "stacking_comparison.csv")

# Parametres principaux et ajustements historiques des modeles.
fits <- setNames(vector("list", length(cfg$models)), cfg$models)
fitted_objects <- setNames(vector("list", length(cfg$models)), cfg$models)
parameter_rows <- list()
parameter_position <- 1L
parameter_patterns <- list(
  lc = c("drift", "sigma_kappa", "phi"),
  rh = c(
    "drift", "sigma_kappa", "psi1", "psi2", "sigma_gamma", "phi"
  ),
  apc = c(
    "drift", "sigma_kappa", "psi1", "psi2", "sigma_gamma", "phi"
  ),
  cbd = c(
    "drift[1]", "drift[2]", "sigma_kappa[1]",
    "sigma_kappa[2]", "Omega[2,1]", "phi"
  ),
  m6 = c(
    "drift[1]", "drift[2]", "sigma_kappa[1]",
    "sigma_kappa[2]", "Omega[2,1]",
    "psi1", "psi2", "sigma_gamma", "phi"
  )
)
for (model in cfg$models) {
  fits[[model]] <- read_model_fit(
    model,
    training,
    cfg,
    paste0("test_training_", cfg$validation_end)
  )
  fitted_objects[[model]] <- posterior_fitted_summary(
    fits[[model]],
    model,
    training,
    ndraws = cfg$forecast_draws,
    seed = cfg$seed
  )
  draws <- fit_draw_matrix(
    fits[[model]],
    variables = parameter_patterns[[model]],
    ndraws = cfg$forecast_draws,
    seed = stable_seed(cfg$seed, "chapter", model, "parameters")
  )
  for (variable in colnames(draws)) {
    values <- draws[, variable]
    parameter_rows[[parameter_position]] <- data.frame(
      model = toupper(model),
      parameter = variable,
      mean = mean(values),
      median = stats::median(values),
      sd = stats::sd(values),
      q025 = stats::quantile(values, 0.025),
      q975 = stats::quantile(values, 0.975),
      stringsAsFactors = FALSE
    )
    parameter_position <- parameter_position + 1L
  }
}
model_parameters <- do.call(rbind, parameter_rows)
write_result(model_parameters, "model_parameter_summary.csv")

observed_training <- processed$long[
  processed$long$year <= cfg$validation_end,
  c("year", "age", "deaths", "exposure")
]
observed_training$observed_rate <-
  observed_training$deaths / observed_training$exposure
fitted_rows <- list()
period_rows <- list()
cohort_rows <- list()
for (model in cfg$models) {
  fitted <- merge(
    fitted_objects[[model]]$fitted,
    observed_training,
    by = c("year", "age"),
    sort = TRUE
  )
  fitted$model <- toupper(model)
  fitted$residual_rate <- fitted$observed_rate - fitted$mean_fitted_rate
  fitted_rows[[model]] <- fitted
  period_rows[[model]] <- fitted_objects[[model]]$period
  if (!is.null(fitted_objects[[model]]$cohort)) {
    cohort_rows[[model]] <- fitted_objects[[model]]$cohort
  }
}
fitted_all <- do.call(rbind, fitted_rows)
period_all <- do.call(rbind, period_rows)
cohort_all <- do.call(rbind, cohort_rows)
write_result(fitted_all, "individual_model_fitted_rates.csv")
write_result(period_all, "individual_model_period_factors.csv")
write_result(cohort_all, "individual_model_cohort_effects.csv")

# ------------------------------ Figures ------------------------------ #

render_chapter_figure(
  "mortality_heatmap",
  function() {
    chapter_plot_theme()
    ages <- sort(unique(data_long$age))
    years <- sort(unique(data_long$year))
    matrix_rate <- xtabs(log(crude_rate) ~ age + year, data_long)
    graphics::image(
      years, ages, t(matrix_rate),
      col = grDevices::hcl.colors(40, "YlOrRd", rev = TRUE),
      xlab = "Annee", ylab = "Age",
      main = "Log-taux bruts de mortalite, 1970-2024"
    )
    graphics::box()
  },
  figures_dir,
  width = 9.5,
  height = 6.2
)

render_chapter_figure(
  "mortality_by_age_selected_years",
  function() {
    chapter_plot_theme()
    selected <- c(1970L, 1990L, 2010L, 2024L)
    palette <- c("#4477AA", "#EE7733", "#228833", "#AA4499")
    part <- data_long[data_long$year %in% selected, ]
    graphics::plot(
      range(part$age), range(part$crude_rate),
      type = "n", log = "y",
      xlab = "Age", ylab = "Taux brut de mortalite (echelle log)",
      main = "Mortalite selon l'age pour quatre annees"
    )
    for (index in seq_along(selected)) {
      values <- part[part$year == selected[[index]], ]
      graphics::lines(
        values$age, values$crude_rate,
        col = palette[[index]], lwd = 2
      )
    }
    graphics::legend(
      "topleft", legend = selected, col = palette,
      lty = 1, lwd = 2, bty = "n", horiz = TRUE
    )
  },
  figures_dir
)

render_chapter_figure(
  "mortality_by_year_selected_ages",
  function() {
    chapter_plot_theme()
    selected <- c(50L, 60L, 70L, 80L, 90L)
    palette <- unname(colors[c("lc", "rh", "apc", "cbd", "m6")])
    part <- data_long[data_long$age %in% selected, ]
    graphics::plot(
      range(part$year), range(part$crude_rate),
      type = "n", log = "y",
      xlab = "Annee", ylab = "Taux brut de mortalite (echelle log)",
      main = "Mortalite au cours du temps pour cinq ages"
    )
    for (index in seq_along(selected)) {
      values <- part[part$age == selected[[index]], ]
      graphics::lines(
        values$year, values$crude_rate,
        col = palette[[index]], lwd = 1.8
      )
    }
    graphics::legend(
      "topright",
      legend = paste(selected, "ans"),
      col = palette, lty = 1, lwd = 2, bty = "n", ncol = 2
    )
  },
  figures_dir
)

render_chapter_figure(
  "deaths_exposures_over_time",
  function() {
    totals <- stats::aggregate(
      cbind(deaths, exposure) ~ year,
      data_long,
      sum
    )
    old <- graphics::par(mfrow = c(2, 1))
    on.exit(graphics::par(old), add = TRUE)
    chapter_plot_theme(c(4.0, 7.0, 2.2, 1))
    graphics::plot(
      totals$year, totals$deaths,
      type = "l", lwd = 2.2, col = "#4477AA",
      xlab = "Annee", ylab = "",
      main = "Nombre total de deces, ages 50-90"
    )
    chapter_plot_theme(c(4.0, 7.0, 2.2, 1))
    graphics::plot(
      totals$year, totals$exposure,
      type = "l", lwd = 2.2, col = "#EE7733",
      xlab = "Annee", ylab = "",
      main = "Exposition totale, ages 50-90"
    )
  },
  figures_dir,
  height = 7.5
)

render_chapter_figure(
  "model_diagnostics",
  function() {
    old <- graphics::par(mfrow = c(2, 2))
    on.exit(graphics::par(old), add = TRUE)
    chapter_plot_theme(c(4.2, 4.2, 2.2, 1))
    positions <- seq_along(cfg$models)
    graphics::plot(
      positions,
      model_diagnostics$max_rhat,
      pch = 16, cex = 1.2, col = "#4477AA",
      xaxt = "n", xlab = "", ylab = "R-hat",
      ylim = c(0.995, 1.052),
      main = "R-hat maximal"
    )
    graphics::axis(1, at = positions, labels = model_diagnostics$model)
    graphics::abline(h = 1.05, lty = 2, lwd = 1.5, col = "#EE7733")

    values <- list(
      min_ess_bulk = list(
        value = model_diagnostics$min_ess_bulk,
        reference = 400,
        title = "ESS bulk minimal", ylab = "ESS"
      ),
      min_ess_tail = list(
        value = model_diagnostics$min_ess_tail,
        reference = 400,
        title = "ESS tail minimal", ylab = "ESS"
      ),
      min_ebfmi = list(
        value = model_diagnostics$min_ebfmi,
        reference = 0.3,
        title = "E-BFMI minimal", ylab = "E-BFMI"
      )
    )
    for (item in values) {
      chapter_plot_theme(c(4.2, 4.2, 2.2, 1))
      graphics::barplot(
        item$value,
        names.arg = model_diagnostics$model,
        col = "#4477AA", border = "#30343B",
        ylab = item$ylab, main = item$title,
        ylim = range(c(0, item$value, item$reference)) * c(1, 1.08)
      )
      graphics::abline(
        h = item$reference,
        lty = 2, lwd = 1.5, col = "#EE7733"
      )
    }
  },
  figures_dir,
  width = 9.5,
  height = 7.5
)

render_chapter_figure(
  "observed_vs_fitted",
  function() {
    old <- graphics::par(mfrow = c(2, 3))
    on.exit(graphics::par(old), add = TRUE)
    limits <- range(
      fitted_all$observed_rate,
      fitted_all$mean_fitted_rate,
      finite = TRUE
    )
    for (model in cfg$models) {
      part <- fitted_all[fitted_all$model == toupper(model), ]
      chapter_plot_theme(c(4.2, 4.2, 2.2, 1))
      graphics::plot(
        part$observed_rate,
        part$mean_fitted_rate,
        log = "xy", pch = 16, cex = 0.35,
        col = grDevices::adjustcolor(colors[[model]], alpha.f = 0.35),
        xlim = limits, ylim = limits,
        xlab = "Taux observe", ylab = "Taux ajuste",
        main = toupper(model)
      )
      graphics::abline(0, 1, lty = 2, col = "#30343B")
    }
    graphics::plot.new()
    graphics::text(
      0.5, 0.55,
      "1970-2015\nAxe logarithmique\nLa diagonale indique l'accord parfait",
      cex = 1.1
    )
  },
  figures_dir,
  width = 10,
  height = 7.4
)

render_chapter_figure(
  "period_factors",
  function() {
    old <- graphics::par(mfrow = c(4, 2))
    on.exit(graphics::par(old), add = TRUE)
    for (model in cfg$models) {
      part <- period_all[period_all$model == model, ]
      factors <- unique(part$factor)
      for (index in seq_along(factors)) {
        values <- part[part$factor == factors[[index]], ]
        chapter_plot_theme(c(4.2, 4.2, 2.2, 1))
        graphics::plot(
          values$year, values$mean,
          type = "n", xlab = "Annee", ylab = "Facteur",
          main = paste(toupper(model), factors[[index]], sep = " - ")
        )
        graphics::polygon(
          c(values$year, rev(values$year)),
          c(values$q025, rev(values$q975)),
          col = grDevices::adjustcolor(
            colors[[model]], alpha.f = 0.16
          ),
          border = NA
        )
        graphics::lines(
          values$year, values$mean,
          col = colors[[model]], lwd = 2
        )
      }
    }
    graphics::plot.new()
  },
  figures_dir,
  width = 10,
  height = 11
)

render_chapter_figure(
  "cohort_effects",
  function() {
    old <- graphics::par(mfrow = c(1, 3))
    on.exit(graphics::par(old), add = TRUE)
    for (model in c("rh", "apc", "m6")) {
      part <- cohort_all[cohort_all$model == model, ]
      chapter_plot_theme(c(4.2, 4.2, 2.2, 1))
      graphics::plot(
        part$cohort, part$mean,
        type = "l", lwd = 2, col = colors[[model]],
        xlab = "Cohorte de naissance", ylab = "Effet de cohorte",
        main = toupper(model)
      )
      graphics::polygon(
        c(part$cohort, rev(part$cohort)),
        c(part$q025, rev(part$q975)),
        col = grDevices::adjustcolor(colors[[model]], alpha.f = 0.18),
        border = NA
      )
      graphics::lines(part$cohort, part$mean, lwd = 2, col = colors[[model]])
      graphics::abline(h = 0, lty = 3, col = "#60656D")
    }
  },
  figures_dir,
  width = 10,
  height = 4.2
)

render_chapter_figure(
  "global_weights",
  function() {
    chapter_plot_theme()
    graphics::barplot(
      global_weights$weight,
      names.arg = global_weights$model,
      col = unname(colors[c("lc", "rh", "apc", "cbd", "m6")]),
      border = "#30343B",
      ylim = c(0, 1),
      ylab = "Poids",
      main = "Poids du stacking global"
    )
    graphics::abline(h = 0, col = "#30343B")
  },
  figures_dir,
  width = 7.5,
  height = 5.2
)

draw_weight_heatmaps <- function(weight_data, value_column, title_prefix,
                                 zlim = c(0, 1), difference = FALSE) {
  old <- graphics::par(mfrow = c(2, 3))
  on.exit(graphics::par(old), add = TRUE)
  palette <- if (difference) {
    grDevices::hcl.colors(41, "Blue-Red 3", rev = TRUE)
  } else {
    grDevices::hcl.colors(40, "Blues 3", rev = TRUE)
  }
  for (model in toupper(cfg$models)) {
    part <- weight_data[weight_data$model == model, ]
    ages <- sort(unique(part$age))
    horizons <- sort(unique(part$horizon))
    matrix_value <- matrix(
      part[[value_column]][match(
        paste(
          rep(horizons, each = length(ages)),
          rep(ages, times = length(horizons))
        ),
        paste(part$horizon, part$age)
      )],
      nrow = length(horizons),
      ncol = length(ages)
    )
    chapter_plot_theme(c(4.2, 4.2, 2.2, 1))
    graphics::image(
      horizons, ages, matrix_value,
      zlim = zlim, col = palette,
      xlab = "Horizon", ylab = "Age",
      main = paste(title_prefix, model)
    )
    graphics::box()
  }
  graphics::plot.new()
}

render_chapter_figure(
  "nonregularized_weight_surfaces",
  function() {
    draw_weight_heatmaps(
      nonregularized_weights, "weight",
      "Non regul.", zlim = c(0, 1)
    )
  },
  figures_dir,
  width = 10,
  height = 7
)

render_chapter_figure(
  "hierarchical_weight_surfaces",
  function() {
    draw_weight_heatmaps(
      hierarchical_weights, "mean",
      "Hierarchique", zlim = c(0, 1)
    )
  },
  figures_dir,
  width = 10,
  height = 7
)

weight_difference <- merge(
  hierarchical_weights[c("age", "horizon", "model", "mean")],
  nonregularized_weights,
  by = c("age", "horizon", "model")
)
weight_difference$difference <- weight_difference$mean -
  weight_difference$weight
weight_limit <- max(abs(weight_difference$difference))
render_chapter_figure(
  "weights_comparison",
  function() {
    draw_weight_heatmaps(
      weight_difference, "difference",
      "Hier. - non regul.",
      zlim = c(-weight_limit, weight_limit),
      difference = TRUE
    )
  },
  figures_dir,
  width = 10,
  height = 7
)

render_chapter_figure(
  "lfo_best_model_heatmap",
  function() {
    chapter_plot_theme(c(4.2, 4.2, 2.2, 7))
    model_index <- match(best_lfo$best_model, toupper(cfg$models))
    matrix_value <- matrix(
      model_index[match(
        paste(
          rep(1:10, each = 41),
          rep(50:90, times = 10)
        ),
        paste(best_lfo$horizon, best_lfo$age)
      )],
      nrow = 10,
      ncol = 41
    )
    graphics::image(
      1:10, 50:90, matrix_value,
      breaks = seq(0.5, 5.5, by = 1),
      col = unname(colors[c("lc", "rh", "apc", "cbd", "m6")]),
      xlim = c(1, 12.4),
      xaxt = "n",
      xlab = "Horizon", ylab = "Age",
      main = "Meilleur modele LFO par age et horizon"
    )
    graphics::axis(1, at = c(1, 2, 4, 6, 8, 10))
    graphics::legend(
      x = 10.5, y = 90, legend = toupper(cfg$models),
      fill = unname(colors[c("lc", "rh", "apc", "cbd", "m6")]),
      bty = "n", cex = 0.8, xpd = NA
    )
  },
  figures_dir,
  width = 8.5,
  height = 6
)

plot_metric_lines <- function(data, x, metrics, titles, filename) {
  render_chapter_figure(
    filename,
    function() {
      old <- graphics::par(mfrow = c(length(metrics), 1L))
      on.exit(graphics::par(old), add = TRUE)
      stacking <- c(
        "stacking_global",
        "stacking_contextual",
        "stacking_hierarchical"
      )
      for (index in seq_along(metrics)) {
        metric <- metrics[[index]]
        part <- data[data$method %in% stacking, ]
        chapter_plot_theme(c(4.2, 4.5, 2.2, 1))
        graphics::plot(
          range(part[[x]]), range(part[[metric]]),
          type = "n", xlab = if (x == "age") "Age" else "Horizon",
          ylab = titles[[index]], main = titles[[index]]
        )
        for (method in stacking) {
          values <- part[part$method == method, ]
          values <- values[order(values[[x]]), ]
          graphics::lines(
            values[[x]], values[[metric]],
            type = "b", pch = c(1, 17, 16)[match(method, stacking)],
            lty = c(2, 3, 1)[match(method, stacking)],
            col = colors[[method]], lwd = 1.8, cex = 0.65
          )
        }
        if (index == 1L) {
          graphics::legend(
            "topright",
            legend = unname(method_labels[stacking]),
            col = unname(colors[stacking]),
            lty = c(2, 3, 1), pch = c(1, 17, 16),
            bty = "n", cex = 0.75
          )
        }
      }
    },
    figures_dir,
    width = 9,
    height = 9
  )
}

plot_metric_lines(
  metrics_by_age,
  "age",
  c("logs", "crps", "mae_deaths"),
  c("LogS moyen", "CRPS moyen", "MAE des deces"),
  "scores_by_age"
)
plot_metric_lines(
  metrics_by_horizon,
  "horizon",
  c("logs", "crps", "mae_deaths"),
  c("LogS moyen", "CRPS moyen", "MAE des deces"),
  "scores_by_horizon"
)

plot_coverage <- function(data, x, filename) {
  render_chapter_figure(
    filename,
    function() {
      old <- graphics::par(mfrow = c(2, 1))
      on.exit(graphics::par(old), add = TRUE)
      stacking <- c(
        "stacking_global",
        "stacking_contextual",
        "stacking_hierarchical"
      )
      for (metric in c("coverage80", "coverage95")) {
        nominal <- if (metric == "coverage80") 0.80 else 0.95
        part <- data[data$method %in% stacking, ]
        chapter_plot_theme(c(4.2, 4.5, 2.2, 1))
        graphics::plot(
          range(part[[x]]), c(0, 1),
          type = "n", xlab = if (x == "age") "Age" else "Horizon",
          ylab = "Couverture",
          main = paste(
            "Couverture predictive",
            if (metric == "coverage80") "80 %" else "95 %"
          )
        )
        graphics::abline(h = nominal, lty = 2, col = "#30343B")
        for (method in stacking) {
          values <- part[part$method == method, ]
          values <- values[order(values[[x]]), ]
          graphics::lines(
            values[[x]], values[[metric]],
            type = "b", pch = 16, cex = 0.55,
            col = colors[[method]], lwd = 1.8
          )
        }
        graphics::legend(
          "bottomleft",
          legend = unname(method_labels[stacking]),
          col = unname(colors[stacking]),
          lty = 1, pch = 16, bty = "n", cex = 0.75
        )
      }
    },
    figures_dir,
    width = 9,
    height = 7.5
  )
}
plot_coverage(metrics_by_age, "age", "coverage_by_age")
plot_coverage(metrics_by_horizon, "horizon", "coverage_by_horizon")

render_chapter_figure(
  "survival_curves",
  function() {
    old <- graphics::par(mfrow = c(2, 3))
    on.exit(graphics::par(old), add = TRUE)
    stacking <- c(
      "stacking_global",
      "stacking_contextual",
      "stacking_hierarchical"
    )
    for (initial_age in sort(unique(survival$initial_age))) {
      part <- survival[
        survival$initial_age == initial_age &
          survival$method %in% stacking,
      ]
      chapter_plot_theme(c(4.2, 4.2, 2.2, 1))
      graphics::plot(
        range(part$horizon), c(0, 1),
        type = "n", xlab = "Horizon", ylab = "Probabilite de survie",
        main = paste("Age initial", initial_age)
      )
      hierarchical_part <- part[
        part$method == "stacking_hierarchical",
      ]
      graphics::polygon(
        c(
          hierarchical_part$horizon,
          rev(hierarchical_part$horizon)
        ),
        c(
          hierarchical_part$q025,
          rev(hierarchical_part$q975)
        ),
        col = grDevices::adjustcolor(
          colors[["stacking_hierarchical"]], alpha.f = 0.15
        ),
        border = NA
      )
      for (method in stacking) {
        values <- part[part$method == method, ]
        graphics::lines(
          values$horizon, values$mean,
          col = colors[[method]],
          lty = c(2, 3, 1)[match(method, stacking)],
          lwd = 1.8
        )
      }
      if (initial_age == min(survival$initial_age)) {
        graphics::legend(
          "bottomleft",
          legend = unname(method_labels[stacking]),
          col = unname(colors[stacking]),
          lty = c(2, 3, 1), lwd = 1.8,
          bty = "n", cex = 0.65
        )
      }
    }
  },
  figures_dir,
  width = 10,
  height = 7.2
)

render_chapter_figure(
  "annuity_comparison",
  function() {
    chapter_plot_theme()
    stacking <- c("stacking_contextual", "stacking_hierarchical")
    part <- annuity[annuity$method %in% stacking, ]
    limits <- range(part$relative_difference_vs_global_pct)
    graphics::plot(
      range(part$initial_age), limits,
      type = "n", xlab = "Age initial",
      ylab = "Ecart relatif au stacking global (%)",
      main = "Valeur des rentes par rapport au stacking global"
    )
    graphics::abline(h = 0, lty = 2, col = "#30343B")
    for (method in stacking) {
      values <- part[part$method == method, ]
      graphics::lines(
        values$initial_age,
        values$relative_difference_vs_global_pct,
        type = "b", pch = if (method == stacking[[1L]]) 17 else 16,
        lty = if (method == stacking[[1L]]) 3 else 1,
        lwd = 2, col = colors[[method]]
      )
    }
    graphics::legend(
      "topright",
      legend = unname(method_labels[stacking]),
      col = unname(colors[stacking]),
      lty = c(3, 1), pch = c(17, 16), lwd = 2, bty = "n"
    )
  },
  figures_dir,
  width = 8.5,
  height = 5.5
)

render_chapter_figure(
  "performance_differences_bootstrap",
  function() {
    part <- bootstrap_differences[
      bootstrap_differences$metric == "crps",
    ]
    labels <- c(
      contextual_minus_global = "Contextuel - global",
      hierarchical_minus_global = "Hierarchique - global",
      hierarchical_minus_contextual = "Hierarchique - contextuel"
    )
    y <- rev(seq_len(nrow(part)))
    xlim <- range(c(part$ci025, part$ci975, 0))
    chapter_plot_theme(c(4.2, 11, 2.2, 1))
    graphics::plot(
      part$mean_difference, y,
      xlim = xlim, ylim = c(0.5, nrow(part) + 0.5),
      yaxt = "n", pch = 16, col = "#4477AA",
      xlab = "Difference moyenne de CRPS",
      ylab = "", main = "Bootstrap par annee cible"
    )
    graphics::segments(
      part$ci025, y, part$ci975, y,
      lwd = 2, col = "#4477AA"
    )
    graphics::axis(
      2, at = y,
      labels = unname(labels[part$comparison]),
      las = 1
    )
    graphics::abline(v = 0, lty = 2, col = "#30343B")
  },
  figures_dir,
  width = 9.5,
  height = 4.8
)

message(
  "Resultats du chapitre 3 prepares : ",
  results_dir, " ; figures : ", figures_dir, "."
)
