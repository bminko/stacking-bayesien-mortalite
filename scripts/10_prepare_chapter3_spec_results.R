#!/usr/bin/env Rscript

# Complements reproductibles demandes par le cahier des charges du chapitre 3.
#
# Le script relit exclusivement les donnees, ajustements et predictions deja
# archives. Il ne reajuste aucun modele de mortalite. Les deux ablations du
# predicteur contextuel sont des optimisations BFGS peu couteuses. Les
# sensibilites de prior hierarchique reposent sur une reponderation
# d'importance des 4 000 tirages existants, avec ESS explicitement rapporte.

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
source(file.path("R", "model_fitting.R"))
source(file.path("R", "forecasting.R"))
source(file.path("R", "aggregation.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "reporting.R"))

cfg <- load_config()
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
predictions <- readRDS(file.path(
  cfg$paths$metrics, "test_predictions.rds"
))

primary_methods <- c(
  cfg$models,
  "stacking_global",
  "stacking_contextual",
  "stacking_hierarchical"
)
model_labels <- c(lc = "LC", rh = "RH", apc = "APC", cbd = "CBD", m6 = "M6")
method_labels <- chapter_method_labels()
colors <- chapter_model_colors()

write_spec <- function(data, filename) {
  write_csv_atomic(data, file.path(results_dir, filename))
  invisible(data)
}

# ---------------------------------------------------------------------------
# 1. Protocole, source et descriptif des donnees
# ---------------------------------------------------------------------------

timeline <- data.frame(
  block = c(
    "Apprentissage des modeles finaux",
    "Origines LFO",
    "Annees cibles LFO",
    "Test final a origine fixe"
  ),
  start = c(1970L, 2000L, 2001L, 2016L),
  end = c(2015L, 2014L, 2015L, 2024L),
  role = c(
    "Ajustement final, aucune observation apres 2015",
    "15 origines; apprentissage 1970:r",
    "Validation; horizons 1 a 10; multiplicite corrigee",
    "Evaluation uniquement; origine 2015; horizons 1 a 9"
  ),
  stringsAsFactors = FALSE
)
write_spec(timeline, "15_timeline_protocol.csv")

death_info <- file.info(cfg$paths$deaths)
exposure_info <- file.info(cfg$paths$exposure)
source_description <- data.frame(
  item = c(
    "source", "country", "population", "death_definition",
    "exposure_definition", "upstream_last_modified",
    "local_download_timestamp_deaths", "local_download_timestamp_exposure",
    "methods_protocol", "integer_tolerance", "rounding_applied",
    "cells", "median_deaths", "median_exposure", "median_crude_rate"
  ),
  value = c(
    "Human Mortality Database (fichiers period 1x1 death.txt et exposure.txt)",
    "Belgique", "Total = Female + Male",
    "Nombre de deces au cours de l'annee civile et a l'age exact indique",
    "Exposition au risque en personnes-annees",
    "21 octobre 2025",
    format(death_info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
    format(exposure_info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
    "HMD Methods Protocol v6 (2017)",
    "abs(D-round(D)) <= 1e-9",
    "Non; le pipeline s'arrete si un deces est fractionnaire",
    nrow(processed$long),
    stats::median(processed$long$deaths),
    stats::median(processed$long$exposure),
    stats::median(processed$long$deaths / processed$long$exposure)
  ),
  stringsAsFactors = FALSE
)
write_spec(source_description, "16_data_source_and_definitions.csv")

render_chapter_figure(
  "protocol_timeline",
  function() {
    chapter_plot_theme(c(4.2, 9.2, 2.0, 1.2))
    y <- rev(seq_len(nrow(timeline)))
    graphics::plot(
      range(timeline$start, timeline$end), range(y) + c(-0.6, 0.6),
      type = "n", yaxt = "n", xlab = "Annee", ylab = "",
      main = "Decoupage temporel du protocole"
    )
    graphics::axis(2, at = y, labels = timeline$block, las = 1)
    palette <- c("#4477AA", "#EE7733", "#CCBB44", "#228833")
    for (i in seq_len(nrow(timeline))) {
      graphics::segments(
        timeline$start[i], y[i], timeline$end[i], y[i],
        lwd = 9, col = palette[i]
      )
      graphics::points(
        c(timeline$start[i], timeline$end[i]), rep(y[i], 2),
        pch = 16, col = palette[i]
      )
    }
    graphics::abline(v = 2015.5, lty = 2, col = "#444444")
    graphics::text(2015.5, max(y) + 0.45, "frontiere test", cex = 0.8)
  },
  figures_dir, width = 9.5, height = 4.6
)

# ---------------------------------------------------------------------------
# 2. Diagnostics des 75 ajustements LFO et des ajustements finaux
# ---------------------------------------------------------------------------

diagnostics <- collect_diagnostic_overviews(cfg)
diagnostics$origin <- suppressWarnings(as.integer(sub(
  ".*lfo_origin_([0-9]{4})/.*", "\\1", diagnostics$label
)))
diagnostics$model <- sub(".*/", "", diagnostics$label)
diagnostics$stage <- ifelse(
  grepl("lfo_origin_", diagnostics$label), "LFO",
  ifelse(grepl("test_training_", diagnostics$label), "final", "other")
)
lfo_diagnostics <- diagnostics[diagnostics$stage == "LFO", ]
lfo_diagnostics <- lfo_diagnostics[order(
  lfo_diagnostics$origin,
  match(lfo_diagnostics$model, cfg$models)
), ]
assert_true(
  nrow(lfo_diagnostics) == 75L,
  "L'audit attend exactement 75 ajustements LFO."
)
write_spec(lfo_diagnostics, "17_lfo_diagnostics_all_75.csv")

lfo_summary <- do.call(rbind, lapply(cfg$models, function(model) {
  part <- lfo_diagnostics[lfo_diagnostics$model == model, ]
  data.frame(
    model = model,
    fits = nrow(part),
    fits_passed = sum(as.logical(part$pass)),
    worst_rhat = max(part$max_rhat),
    minimum_ess_bulk = min(part$min_ess_bulk),
    minimum_ess_tail = min(part$min_ess_tail),
    divergences = sum(part$divergences),
    max_treedepth_hits = sum(part$max_treedepth),
    stringsAsFactors = FALSE
  )
}))
write_spec(lfo_summary, "18_lfo_diagnostics_summary.csv")

final_diagnostics <- diagnostics[
  diagnostics$stage %in% c("final", "other") &
    !grepl("lfo_origin_|simulation|actuarial_training", diagnostics$label),
]
write_spec(final_diagnostics, "19_final_and_aggregation_diagnostics.csv")

hier_dir <- file.path(cfg$paths$fits, "aggregation", "hierarchical")
hier_parameters <- utils::read.csv(
  file.path(hier_dir, "diagnostics_parameters.csv"),
  stringsAsFactors = FALSE
)
hier_sampler <- utils::read.csv(
  file.path(hier_dir, "diagnostics_sampler.csv"),
  stringsAsFactors = FALSE
)
write_spec(hier_parameters, "20_hierarchical_parameter_diagnostics.csv")
write_spec(hier_sampler, "21_hierarchical_sampler_diagnostics.csv")

fit_elapsed_seconds <- function(path) {
  fit <- readRDS(path)
  simulation <- attr(fit, "sim", exact = TRUE)
  per_chain <- vapply(simulation$samples, function(chain) {
    elapsed <- attr(chain, "elapsed_time", exact = TRUE)
    sum(as.numeric(elapsed))
  }, numeric(1))
  max(per_chain)
}

fit_paths <- list.files(
  cfg$paths$fits, pattern = "^fit[.]rds$",
  recursive = TRUE, full.names = TRUE
)
runtime <- data.frame(
  path = normalizePath(fit_paths, winslash = "/", mustWork = TRUE),
  elapsed_seconds_parallel_wall_proxy = vapply(
    fit_paths, fit_elapsed_seconds, numeric(1)
  ),
  stringsAsFactors = FALSE
)
runtime$stage <- ifelse(
  grepl("lfo_origin_", runtime$path), "LFO",
  ifelse(grepl("test_training_2015", runtime$path), "ajustements_finaux",
    ifelse(grepl("aggregation/hierarchical", runtime$path),
      "stacking_hierarchique",
      ifelse(grepl("simulation", runtime$path), "simulation", "autre")
    )
  )
)
runtime_summary <- do.call(rbind, lapply(
  split(runtime, runtime$stage),
  function(part) data.frame(
    stage = part$stage[1L],
    fits = nrow(part),
    elapsed_hours_parallel_wall_proxy =
      sum(part$elapsed_seconds_parallel_wall_proxy) / 3600,
    stringsAsFactors = FALSE
  )
))
runtime_summary <- rbind(
  runtime_summary,
  data.frame(
    stage = "total_archive",
    fits = nrow(runtime),
    elapsed_hours_parallel_wall_proxy =
      sum(runtime$elapsed_seconds_parallel_wall_proxy) / 3600
  )
)
write_spec(runtime, "22_fit_runtime_details.csv")
write_spec(runtime_summary, "23_fit_runtime_summary.csv")

software <- data.frame(
  component = c(
    "R executable used by pipeline", "Stan compiler",
    "RStan build note", "operating system", "chains",
    "parallel cores", "warmup per chain", "kept draws per chain",
    "total kept draws", "adapt_delta", "max_treedepth",
    "metric LC and CBD", "metric RH, APC and M6",
    "metric hierarchical stacking"
  ),
  value = c(
    "R 4.3.1 (command archived in README)",
    "stanc 2.32.2 (embedded in compiled model)",
    "RStan binary compiled under R 4.3.3; exact package patch not archived",
    "Windows x86_64 mingw32",
    cfg$mcmc$chains, cfg$mcmc$parallel_chains,
    cfg$mcmc$iter_warmup, cfg$mcmc$iter_sampling,
    cfg$mcmc$chains * cfg$mcmc$iter_sampling,
    cfg$mcmc$adapt_delta, cfg$mcmc$max_treedepth,
    cfg$mcmc$metric, "dense_e", cfg$hierarchical_mcmc$metric
  ),
  stringsAsFactors = FALSE
)
write_spec(software, "24_software_and_sampling_settings.csv")

# ---------------------------------------------------------------------------
# 3. Meta-LFO, multiplicite et optimisations
# ---------------------------------------------------------------------------

origin_year <- unique(meta[c("origin", "horizon", "year")])
multiplicity <- as.data.frame(table(origin_year$year))
names(multiplicity) <- c("year", "multiplicity")
multiplicity$year <- as.integer(as.character(multiplicity$year))
multiplicity$multiplicity <- as.integer(multiplicity$multiplicity)
meta_weight_check <- merge(
  unique(meta[c("origin", "horizon", "year", "omega")]),
  multiplicity, by = "year", all.x = TRUE
)
meta_weight_check$expected_omega <- 1 / meta_weight_check$multiplicity
meta_weight_check$absolute_error <-
  abs(meta_weight_check$omega - meta_weight_check$expected_omega)
write_spec(meta_weight_check, "25_lfo_multiplicity_and_omega.csv")

year_weight_totals <- stats::aggregate(
  omega ~ year, data = meta, FUN = sum
)
names(year_weight_totals)[2L] <- "sum_omega_all_ages"
per_age <- stats::aggregate(
  omega ~ year + age, data = meta, FUN = sum
)
per_age_summary <- stats::aggregate(
  omega ~ year, data = per_age,
  FUN = function(x) c(min = min(x), mean = mean(x), max = max(x))
)
per_age_values <- per_age_summary$omega
if (is.null(dim(per_age_values))) {
  per_age_values <- do.call(rbind, per_age_values)
}
year_weight_totals$sum_omega_per_age_min <- per_age_values[, "min"]
year_weight_totals$sum_omega_per_age_mean <- per_age_values[, "mean"]
year_weight_totals$sum_omega_per_age_max <- per_age_values[, "max"]
write_spec(year_weight_totals, "26_lfo_year_weight_totals.csv")

weighted_vs_unweighted <- do.call(rbind, lapply(cfg$models, function(model) {
  logp <- meta[[paste0("log_p_", model)]]
  by_year <- stats::aggregate(
    logp, by = list(year = meta$year),
    FUN = mean
  )
  data.frame(
    model = model,
    weighted_elpd = sum(meta$omega * logp),
    unweighted_elpd = sum(logp),
    weighted_mean_log_predictive =
      sum(meta$omega * logp) / sum(meta$omega),
    unweighted_mean_log_predictive = mean(logp),
    sd_year_mean_log_predictive = stats::sd(by_year$x),
    stringsAsFactors = FALSE
  )
}))
weighted_vs_unweighted$weighted_rank <- rank(
  -weighted_vs_unweighted$weighted_mean_log_predictive,
  ties.method = "min"
)
weighted_vs_unweighted$unweighted_rank <- rank(
  -weighted_vs_unweighted$unweighted_mean_log_predictive,
  ties.method = "min"
)
write_spec(weighted_vs_unweighted, "27_weighted_vs_unweighted_lfo.csv")

meta_unweighted <- meta
meta_unweighted$omega <- 1
global_unweighted <- fit_global_stacking(
  meta_unweighted, cfg$models, cfg$contextual_multistarts,
  stable_seed(cfg$seed, "audit", "global_unweighted")
)
context_unweighted <- fit_contextual_stacking(
  meta_unweighted, cfg$models, cfg$contextual_multistarts,
  stable_seed(cfg$seed, "audit", "context_unweighted")
)
weight_sensitivity <- rbind(
  data.frame(
    scheme = "weighted_multiplicity",
    method = "stacking_global",
    model = cfg$models,
    weight = aggregation$global$weights[cfg$models]
  ),
  data.frame(
    scheme = "unweighted",
    method = "stacking_global",
    model = cfg$models,
    weight = global_unweighted$weights[cfg$models]
  )
)
write_spec(weight_sensitivity, "28_global_weights_multiplicity_sensitivity.csv")

global_gradient <- global_stacking_gradient(
  aggregation$global$optimization$par,
  meta_logp_matrix(meta, cfg$models),
  meta$omega
)
context_gradient <- contextual_stacking_gradient(
  aggregation$contextual$optimization$par,
  meta_logp_matrix(meta, cfg$models),
  context_design_matrix(meta),
  meta$omega
)
optimization_audit <- data.frame(
  method = c("stacking_global", "stacking_contextual"),
  algorithm = "BFGS",
  starts = c(
    length(aggregation$global$all_values),
    length(aggregation$contextual$all_values)
  ),
  convergence_code = c(
    aggregation$global$optimization$convergence,
    aggregation$contextual$optimization$convergence
  ),
  iterations_function = c(
    aggregation$global$optimization$counts[["function"]],
    aggregation$contextual$optimization$counts[["function"]]
  ),
  objective = c(
    aggregation$global$optimization$value,
    aggregation$contextual$optimization$value
  ),
  objective_range_across_starts = c(
    diff(range(aggregation$global$all_values)),
    diff(range(aggregation$contextual$all_values))
  ),
  max_absolute_gradient = c(
    max(abs(global_gradient)), max(abs(context_gradient))
  ),
  stringsAsFactors = FALSE
)
write_spec(optimization_audit, "29_stacking_optimization_audit.csv")

validation_grid <- expand.grid(
  age = cfg$age_min:cfg$age_max,
  horizon = cfg$lfo_horizons
)
nonreg_grid <- contextual_weight_grid(
  validation_grid, aggregation$contextual, cfg$models,
  aggregation$standardization
)
extreme_nonreg <- do.call(rbind, lapply(cfg$models, function(model) {
  weight <- nonreg_grid[[model]]
  data.frame(
    model = model,
    proportion_below_001 = mean(weight < 0.01),
    proportion_above_099 = mean(weight > 0.99),
    minimum = min(weight),
    median = stats::median(weight),
    maximum = max(weight),
    stringsAsFactors = FALSE
  )
}))
write_spec(extreme_nonreg, "30_nonregularized_extreme_weights.csv")

logp_meta <- meta_logp_matrix(meta, cfg$models)
global_matrix <- constant_weight_matrix(
  aggregation$global$weights, nrow(meta), cfg$models
)
nonreg_meta <- contextual_weight_grid(
  meta[c("age", "horizon")], aggregation$contextual, cfg$models,
  aggregation$standardization
)
hier_meta_long <- hierarchical_weight_grid(
  meta[c("age", "horizon")],
  readRDS(aggregation$hierarchical_fit_path),
  cfg$models, aggregation$standardization,
  ndraws = cfg$forecast_draws,
  seed = stable_seed(cfg$seed, "audit", "validation_hier")
)
hier_meta <- weight_matrix_from_long(
  hier_meta_long, meta[c("age", "horizon")], cfg$models
)
validation_methods <- list(
  pseudo_bma = constant_weight_matrix(
    aggregation$pseudo_bma$weights, nrow(meta), cfg$models
  ),
  stacking_global = global_matrix,
  stacking_contextual = as.matrix(nonreg_meta[cfg$models]),
  stacking_hierarchical = hier_meta
)
validation_comparison <- do.call(rbind, lapply(
  names(validation_methods),
  function(method) {
    log_mix <- row_log_sum_exp(
      log(validation_methods[[method]]) + logp_meta
    )
    data.frame(
      method = method,
      weighted_negative_log_score = -sum(meta$omega * log_mix),
      weighted_mean_log_score = -sum(meta$omega * log_mix) /
        sum(meta$omega),
      unweighted_mean_log_score = -mean(log_mix),
      stringsAsFactors = FALSE
    )
  }
))
validation_comparison$weighted_rank <- rank(
  validation_comparison$weighted_mean_log_score,
  ties.method = "min"
)
write_spec(validation_comparison, "31_validation_method_comparison.csv")

# ---------------------------------------------------------------------------
# 4. Diagnostics d'ajustement, residus et PPC conditionnels
# ---------------------------------------------------------------------------

fitted <- utils::read.csv(
  file.path(results_dir, "individual_model_fitted_rates.csv"),
  stringsAsFactors = FALSE
)
parameters <- utils::read.csv(
  file.path(results_dir, "model_parameter_summary.csv"),
  stringsAsFactors = FALSE
)
phi <- parameters[parameters$parameter == "phi", c("model", "mean")]
names(phi)[2L] <- "phi"
fitted <- merge(fitted, phi, by = "model", all.x = TRUE)
fitted$expected_deaths <- fitted$mean_fitted_rate * fitted$exposure
fitted$residual_pearson <- (
  fitted$deaths - fitted$expected_deaths
) / sqrt(
  fitted$expected_deaths +
    fitted$expected_deaths^2 / fitted$phi
)
write_spec(fitted, "32_standardized_residuals_all_cells.csv")

residual_by_age <- stats::aggregate(
  residual_pearson ~ model + age,
  data = fitted,
  FUN = function(x) c(mean = mean(x), rms = sqrt(mean(x^2)))
)
age_values <- residual_by_age$residual_pearson
if (is.null(dim(age_values))) {
  age_values <- do.call(rbind, age_values)
}
residual_by_age$mean_residual <- age_values[, "mean"]
residual_by_age$rms_residual <- age_values[, "rms"]
residual_by_age$residual_pearson <- NULL
residual_by_year <- stats::aggregate(
  residual_pearson ~ model + year,
  data = fitted,
  FUN = function(x) c(mean = mean(x), rms = sqrt(mean(x^2)))
)
year_values <- residual_by_year$residual_pearson
if (is.null(dim(year_values))) {
  year_values <- do.call(rbind, year_values)
}
residual_by_year$mean_residual <- year_values[, "mean"]
residual_by_year$rms_residual <- year_values[, "rms"]
residual_by_year$residual_pearson <- NULL
write_spec(residual_by_age, "33_standardized_residuals_by_age.csv")
write_spec(residual_by_year, "34_standardized_residuals_by_year.csv")

ppc_rows <- list()
ppc_discrepancy <- list()
position <- 1L
for (model in unique(fitted$model)) {
  part_model <- fitted[fitted$model == model, ]
  model_phi <- unique(part_model$phi)
  for (year in sort(unique(part_model$year))) {
    part <- part_model[part_model$year == year, ]
    set.seed(stable_seed(cfg$seed, "ppc", model, year))
    replicated <- replicate(
      1000L,
      sum(stats::rnbinom(
        nrow(part), mu = part$expected_deaths, size = model_phi
      ))
    )
    ppc_rows[[position]] <- data.frame(
      model = model,
      year = year,
      observed_total_deaths = sum(part$deaths),
      predictive_mean_total_deaths = mean(replicated),
      predictive_q025_total_deaths =
        stats::quantile(replicated, 0.025, names = FALSE),
      predictive_q975_total_deaths =
        stats::quantile(replicated, 0.975, names = FALSE),
      predictive_tail_probability =
        mean(replicated >= sum(part$deaths)),
      stringsAsFactors = FALSE
    )
    position <- position + 1L
  }
  ppc_discrepancy[[model]] <- data.frame(
    model = model,
    mean_residual = mean(part_model$residual_pearson),
    rms_residual = sqrt(mean(part_model$residual_pearson^2)),
    proportion_abs_residual_gt_2 =
      mean(abs(part_model$residual_pearson) > 2),
    proportion_observed_rate_in_parameter_interval =
      mean(
        part_model$observed_rate >= part_model$q025_fitted_rate &
          part_model$observed_rate <= part_model$q975_fitted_rate
      ),
    stringsAsFactors = FALSE
  )
}
ppc_totals <- do.call(rbind, ppc_rows)
ppc_discrepancy <- do.call(rbind, ppc_discrepancy)
write_spec(ppc_totals, "35_ppc_total_deaths_by_year.csv")
write_spec(ppc_discrepancy, "36_ppc_residual_discrepancies.csv")

render_chapter_figure(
  "standardized_residual_heatmaps",
  function() {
    graphics::par(mfrow = c(2, 3), mar = c(3.2, 3.4, 2.2, 0.8))
    for (model in names(model_labels)) {
      label <- model_labels[[model]]
      part <- fitted[fitted$model == label, ]
      matrix_values <- xtabs(residual_pearson ~ age + year, data = part)
      graphics::image(
        as.integer(colnames(matrix_values)),
        as.integer(rownames(matrix_values)),
        t(matrix_values),
        col = grDevices::hcl.colors(31, "Blue-Red 3", rev = TRUE),
        zlim = c(-4, 4), xlab = "Annee", ylab = "Age",
        main = label
      )
    }
    graphics::plot.new()
    graphics::text(
      0.5, 0.6, "Residus de Pearson",
      cex = 1.1, font = 2
    )
    graphics::text(
      0.5, 0.45, "bleu: negatif  |  rouge: positif",
      cex = 0.85
    )
  },
  figures_dir, width = 10, height = 7
)

render_chapter_figure(
  "standardized_residual_profiles",
  function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.4, 1.0))
    palette <- setNames(
      colors[names(model_labels)], unname(model_labels)
    )
    graphics::plot(
      range(residual_by_age$age),
      range(residual_by_age$mean_residual),
      type = "n", xlab = "Age", ylab = "Residuel moyen",
      main = "Residus standardises selon l'age"
    )
    for (model in unique(residual_by_age$model)) {
      part <- residual_by_age[residual_by_age$model == model, ]
      graphics::lines(
        part$age, part$mean_residual,
        col = palette[[model]], lwd = 2
      )
    }
    graphics::abline(h = 0, lty = 2)
    graphics::plot(
      range(residual_by_year$year),
      range(residual_by_year$mean_residual),
      type = "n", xlab = "Annee", ylab = "Residuel moyen",
      main = "Residus standardises selon l'annee"
    )
    for (model in unique(residual_by_year$model)) {
      part <- residual_by_year[residual_by_year$model == model, ]
      graphics::lines(
        part$year, part$mean_residual,
        col = palette[[model]], lwd = 2
      )
    }
    graphics::abline(h = 0, lty = 2)
    graphics::legend(
      "topleft", legend = names(palette), col = palette,
      lwd = 2, bty = "n", cex = 0.75
    )
  },
  figures_dir, width = 10, height = 4.8
)

render_chapter_figure(
  "ppc_total_deaths",
  function() {
    graphics::par(mfrow = c(2, 3), mar = c(3.2, 3.7, 2.2, 0.8))
    for (model in unique(ppc_totals$model)) {
      part <- ppc_totals[ppc_totals$model == model, ]
      yrange <- range(
        part$observed_total_deaths,
        part$predictive_q025_total_deaths,
        part$predictive_q975_total_deaths
      )
      graphics::plot(
        part$year, part$observed_total_deaths,
        type = "l", lwd = 1.7, col = "#222222",
        xlab = "Annee", ylab = "Deces totaux",
        ylim = yrange, main = model
      )
      graphics::polygon(
        c(part$year, rev(part$year)),
        c(
          part$predictive_q025_total_deaths,
          rev(part$predictive_q975_total_deaths)
        ),
        border = NA, col = grDevices::adjustcolor("#4477AA", 0.22)
      )
      graphics::lines(
        part$year, part$predictive_mean_total_deaths,
        col = "#4477AA", lwd = 1.6
      )
      graphics::lines(
        part$year, part$observed_total_deaths,
        col = "#222222", lwd = 1.4
      )
    }
    graphics::plot.new()
    graphics::legend(
      "center",
      legend = c("Observe", "Predictif moyen", "Intervalle 95 %"),
      col = c("#222222", "#4477AA", grDevices::adjustcolor("#4477AA", 0.35)),
      lwd = c(2, 2, 8), bty = "n"
    )
  },
  figures_dir, width = 10, height = 7
)

# ---------------------------------------------------------------------------
# 5. Test final, classements, couverture et sensibilite pandemique
# ---------------------------------------------------------------------------

primary_cells <- evaluation$cells[
  evaluation$cells$method %in% primary_methods,
  ,
  drop = FALSE
]
periods <- list(
  full_2016_2024 = 2016:2024,
  initial_2016_2020 = 2016:2020,
  extension_2021_2024 = 2021:2024
)
period_summary <- summarize_evaluation_periods(primary_cells, periods)$overall
period_summary$coverage80_gap_points <-
  100 * (period_summary$coverage80 - 0.80)
period_summary$coverage95_gap_points <-
  100 * (period_summary$coverage95 - 0.95)
period_summary$rank_logs <- ave(
  period_summary$logs, period_summary$period,
  FUN = function(x) rank(x, ties.method = "min")
)
period_summary$rank_crps <- ave(
  period_summary$crps, period_summary$period,
  FUN = function(x) rank(x, ties.method = "min")
)
period_summary$rank_mae <- ave(
  period_summary$mae_deaths, period_summary$period,
  FUN = function(x) rank(x, ties.method = "min")
)
write_spec(period_summary, "37_test_metrics_8_methods_3_periods.csv")

relative_rows <- list()
position <- 1L
for (period in names(periods)) {
  part <- period_summary[period_summary$period == period, ]
  hierarchical <- part[
    part$method == "stacking_hierarchical", , drop = FALSE
  ]
  for (comparison in setdiff(part$method, "stacking_hierarchical")) {
    comparator <- part[part$method == comparison, , drop = FALSE]
    relative_rows[[position]] <- data.frame(
      period = period,
      comparator = comparison,
      logs_reduction_pct =
        100 * (comparator$logs - hierarchical$logs) / comparator$logs,
      crps_reduction_pct =
        100 * (comparator$crps - hierarchical$crps) / comparator$crps,
      mae_reduction_pct =
        100 * (comparator$mae_deaths - hierarchical$mae_deaths) /
          comparator$mae_deaths,
      coverage80_difference_points =
        100 * (hierarchical$coverage80 - comparator$coverage80),
      coverage95_difference_points =
        100 * (hierarchical$coverage95 - comparator$coverage95),
      stringsAsFactors = FALSE
    )
    position <- position + 1L
  }
}
write_spec(
  do.call(rbind, relative_rows),
  "38_hierarchical_relative_test_performance.csv"
)

coverage_by_age <- stats::aggregate(
  cbind(coverage80, coverage95, width80, width95) ~ method + age,
  data = primary_cells,
  FUN = sum
)
names(coverage_by_age)[3:6] <- c(
  "covered80_years_out_of_9", "covered95_years_out_of_9",
  "sum_width80", "sum_width95"
)
mean_width <- stats::aggregate(
  cbind(width80, width95) ~ method + age,
  data = primary_cells,
  FUN = mean
)
coverage_by_age <- merge(
  coverage_by_age, mean_width,
  by = c("method", "age"), all = TRUE
)
coverage_by_age$coverage80_gap_points <-
  100 * (coverage_by_age$covered80_years_out_of_9 / 9 - 0.80)
coverage_by_age$coverage95_gap_points <-
  100 * (coverage_by_age$covered95_years_out_of_9 / 9 - 0.95)
write_spec(coverage_by_age, "39_coverage_counts_by_age_out_of_9.csv")

year_horizon <- group_mean(
  primary_cells,
  evaluation_metric_names(),
  c("method", "year", "horizon")
)
names(year_horizon)[names(year_horizon) == "absolute_error_deaths"] <-
  "mae_deaths"
write_spec(year_horizon, "40_performance_by_target_year_and_horizon.csv")

pandemic_periods <- list(
  full_2016_2024 = 2016:2024,
  excluding_2020_2022 = c(2016:2019, 2023:2024),
  pre_pandemic_2016_2019 = 2016:2019,
  pandemic_2020 = 2020L,
  pandemic_2021_2022 = 2021:2022,
  post_pandemic_2023_2024 = 2023:2024
)
pandemic <- summarize_evaluation_periods(
  primary_cells, pandemic_periods
)$overall
pandemic$coverage80_gap_points <- 100 * (pandemic$coverage80 - 0.80)
pandemic$coverage95_gap_points <- 100 * (pandemic$coverage95 - 0.95)
for (metric in c("logs", "crps", "mae_deaths")) {
  pandemic[[paste0("rank_", metric)]] <- ave(
    pandemic[[metric]], pandemic$period,
    FUN = function(x) rank(x, ties.method = "min")
  )
}
write_spec(pandemic, "41_pandemic_sensitivity_6_blocks.csv")

jackknife_rows <- list()
position <- 1L
for (excluded_year in 2016:2024) {
  part <- primary_cells[primary_cells$year != excluded_year, ]
  summary <- summarize_evaluation_cells(part)$overall
  for (comparison in c(
    "stacking_global", "stacking_contextual"
  )) {
    hierarchical <- summary[
      summary$method == "stacking_hierarchical", ]
    comparator <- summary[summary$method == comparison, ]
    jackknife_rows[[position]] <- data.frame(
      excluded_year = excluded_year,
      comparison = paste(
        "stacking_hierarchical_minus", comparison, sep = "_"
      ),
      logs_difference = hierarchical$logs - comparator$logs,
      crps_difference = hierarchical$crps - comparator$crps,
      mae_difference =
        hierarchical$mae_deaths - comparator$mae_deaths,
      stringsAsFactors = FALSE
    )
    position <- position + 1L
  }
}
write_spec(
  do.call(rbind, jackknife_rows),
  "42_leave_one_target_year_out_sensitivity.csv"
)

render_chapter_figure(
  "performance_by_target_year_horizon",
  function() {
    chapter_plot_theme()
    selected <- c(
      "stacking_global", "stacking_contextual",
      "stacking_hierarchical"
    )
    part <- year_horizon[year_horizon$method %in% selected, ]
    graphics::plot(
      range(part$year), range(part$logs),
      type = "n", xlab = "Annee cible (horizon depuis 2015)",
      ylab = "LogS moyen",
      main = "Performances selon l'annee cible et l'horizon depuis 2015"
    )
    for (method in selected) {
      piece <- part[part$method == method, ]
      graphics::lines(
        piece$year, piece$logs,
        col = colors[[method]], lwd = 2, type = "b", pch = 16
      )
    }
    graphics::axis(
      3, at = 2016:2024, labels = 1:9,
      cex.axis = 0.75
    )
    graphics::mtext("Horizon h", side = 3, line = 2.2)
    graphics::legend(
      "topleft", legend = method_labels[selected],
      col = colors[selected], lwd = 2, pch = 16,
      bty = "n", cex = 0.8
    )
  },
  figures_dir, width = 9, height = 5.5
)

render_chapter_figure(
  "pandemic_sensitivity_scores",
  function() {
    selected <- c(
      "stacking_global", "stacking_contextual",
      "stacking_hierarchical"
    )
    part <- pandemic[pandemic$method %in% selected, ]
    matrix_values <- xtabs(logs ~ method + period, data = part)[selected, ]
    graphics::par(mar = c(8.5, 4.5, 2.2, 1.0))
    graphics::barplot(
      matrix_values, beside = TRUE,
      col = colors[selected], border = NA,
      ylab = "LogS moyen",
      main = "Sensibilite des scores aux annees pandemiques",
      las = 2, cex.names = 0.7
    )
    graphics::legend(
      "topleft", legend = method_labels[selected],
      fill = colors[selected], border = NA, bty = "n", cex = 0.8
    )
  },
  figures_dir, width = 9.5, height = 6
)

# ---------------------------------------------------------------------------
# 6. Ablations contextuelles et sensibilites de prior hierarchique
# ---------------------------------------------------------------------------

test_target <- predictions[[cfg$models[[1L]]]]$summary[c(
  "year", "age", "deaths", "exposure", "origin", "horizon"
)]

evaluate_custom_weights <- function(weights, label) {
  logp <- sapply(cfg$models, function(model) {
    predictions[[model]]$summary$log_predictive
  })
  cells <- vector("list", nrow(test_target))
  for (cell in seq_len(nrow(test_target))) {
    count_draws <- mixture_draws_for_cell(
      predictions, cfg$models, cell, weights[cell, ], "count",
      stable_seed(cfg$seed, "robustness", label, cell, "count")
    )
    rate_draws <- mixture_draws_for_cell(
      predictions, cfg$models, cell, weights[cell, ], "rate",
      stable_seed(cfg$seed, "robustness", label, cell, "rate")
    )
    interval80 <- stats::quantile(
      rate_draws, c(0.10, 0.90), names = FALSE
    )
    interval95 <- stats::quantile(
      rate_draws, c(0.025, 0.975), names = FALSE
    )
    observed_rate <- test_target$deaths[cell] /
      test_target$exposure[cell]
    cells[[cell]] <- data.frame(
      method = label,
      year = test_target$year[cell],
      age = test_target$age[cell],
      horizon = test_target$horizon[cell],
      deaths = test_target$deaths[cell],
      exposure = test_target$exposure[cell],
      observed_rate = observed_rate,
      logs = -log_sum_exp(log(weights[cell, ]) + logp[cell, ]),
      crps = empirical_crps(rate_draws, observed_rate),
      absolute_error_deaths =
        abs(test_target$deaths[cell] - mean(count_draws)),
      coverage80 = as.integer(
        observed_rate >= interval80[1L] &
          observed_rate <= interval80[2L]
      ),
      width80 = diff(interval80),
      coverage95 = as.integer(
        observed_rate >= interval95[1L] &
          observed_rate <= interval95[2L]
      ),
      width95 = diff(interval95),
      stringsAsFactors = FALSE
    )
  }
  cells <- do.call(rbind, cells)
  summary <- summarize_evaluation_cells(cells)$overall
  list(cells = cells, overall = summary)
}

ablation_no_age2 <- fit_contextual_stacking(
  meta, cfg$models, cfg$contextual_multistarts,
  stable_seed(cfg$seed, "ablation", "no_age2"),
  design_columns = c("age", "horizon", "age_horizon")
)
ablation_no_interaction <- fit_contextual_stacking(
  meta, cfg$models, cfg$contextual_multistarts,
  stable_seed(cfg$seed, "ablation", "no_interaction"),
  design_columns = c("age", "age2", "horizon")
)
weights_no_age2 <- as.matrix(contextual_weight_grid(
  test_target[c("age", "horizon")],
  ablation_no_age2, cfg$models,
  aggregation$standardization
)[cfg$models])
weights_no_interaction <- as.matrix(contextual_weight_grid(
  test_target[c("age", "horizon")],
  ablation_no_interaction, cfg$models,
  aggregation$standardization
)[cfg$models])

hier_fit <- readRDS(aggregation$hierarchical_fit_path)
hier_draws <- fit_draw_matrix(
  hier_fit, variables = c("alpha", "beta", "tau")
)
K <- length(cfg$models)
P <- 4L
alpha_draws <- extract_vector_parameter(hier_draws, "alpha", K - 1L)
beta_draws <- array(
  NA_real_, dim = c(nrow(hier_draws), K - 1L, P)
)
for (k in seq_len(K - 1L)) {
  for (p in seq_len(P)) {
    beta_draws[, k, p] <- extract_scalar_parameter(
      hier_draws, sprintf("beta[%d,%d]", k, p)
    )
  }
}
tau_draws <- extract_vector_parameter(hier_draws, "tau", P)
old_scale <- c(1, 0.5, 1, 0.5)

importance_weights_for_scale <- function(new_scale) {
  log_ratio <- rowSums(sapply(seq_len(P), function(p) {
    stats::dnorm(
      tau_draws[, p], 0, new_scale[p], log = TRUE
    ) - stats::dnorm(
      tau_draws[, p], 0, old_scale[p], log = TRUE
    )
  }))
  weight <- exp(log_ratio - max(log_ratio))
  weight / sum(weight)
}

hier_weight_matrix <- function(
    grid, draw_weights = rep(1 / nrow(hier_draws), nrow(hier_draws)),
    zero_terms = integer()) {
  standardized <- standardize_context(
    grid, aggregation$standardization
  )$data
  X <- context_design_matrix(standardized)
  result <- matrix(
    NA_real_, nrow(grid), K,
    dimnames = list(NULL, cfg$models)
  )
  beta_used <- beta_draws
  if (length(zero_terms)) beta_used[, , zero_terms] <- 0
  for (cell in seq_len(nrow(grid))) {
    scores <- alpha_draws
    for (p in seq_len(P)) {
      scores <- scores + beta_used[, , p] * X[cell, p]
    }
    weights <- scores_to_weights(cbind(scores, 0))
    result[cell, ] <- colSums(weights * draw_weights)
  }
  result
}

prior_strong <- importance_weights_for_scale(old_scale / 2)
prior_broad <- importance_weights_for_scale(old_scale * 2)
robust_weight_sets <- list(
  nonreg_no_age2 = weights_no_age2,
  nonreg_no_interaction = weights_no_interaction,
  hier_ablation_no_age2 = hier_weight_matrix(
    test_target[c("age", "horizon")], zero_terms = 2L
  ),
  hier_ablation_no_interaction = hier_weight_matrix(
    test_target[c("age", "horizon")], zero_terms = 4L
  ),
  hier_prior_stronger = hier_weight_matrix(
    test_target[c("age", "horizon")], prior_strong
  ),
  hier_prior_broader = hier_weight_matrix(
    test_target[c("age", "horizon")], prior_broad
  )
)
robust_evaluations <- lapply(
  names(robust_weight_sets),
  function(label) evaluate_custom_weights(
    robust_weight_sets[[label]], label
  )
)
names(robust_evaluations) <- names(robust_weight_sets)
robustness <- do.call(rbind, lapply(
  names(robust_evaluations),
  function(label) robust_evaluations[[label]]$overall
))
rownames(robustness) <- NULL
prior_ess <- data.frame(
  sensitivity = c("stronger_half_scale", "broader_double_scale"),
  importance_ess = c(
    1 / sum(prior_strong^2),
    1 / sum(prior_broad^2)
  ),
  total_draws = nrow(hier_draws),
  maximum_normalized_weight = c(max(prior_strong), max(prior_broad)),
  stringsAsFactors = FALSE
)
write_spec(robustness, "43_context_and_prior_robustness_metrics.csv")
write_spec(prior_ess, "44_prior_importance_reweighting_diagnostics.csv")

nonreg_existing <- utils::read.csv(
  file.path(results_dir, "05_nonregularized_weights.csv"),
  stringsAsFactors = FALSE
)
hier_existing <- utils::read.csv(
  file.path(results_dir, "06_hierarchical_weights.csv"),
  stringsAsFactors = FALSE
)
weight_difference <- merge(
  hier_existing[c("age", "horizon", "model", "mean", "q025", "q975")],
  nonreg_existing[c("age", "horizon", "model", "weight")],
  by = c("age", "horizon", "model")
)
weight_difference$difference_hier_minus_nonreg <-
  weight_difference$mean - weight_difference$weight
weight_difference$posterior_interval_width <-
  weight_difference$q975 - weight_difference$q025
write_spec(weight_difference, "45_hierarchical_vs_nonregularized_weights.csv")

render_chapter_figure(
  "hierarchical_nonregularized_weight_difference",
  function() {
    graphics::par(mfrow = c(2, 3), mar = c(3.3, 3.5, 2.1, 0.8))
    for (model in cfg$models) {
      label <- model_labels[[model]]
      part <- weight_difference[weight_difference$model == label, ]
      z <- xtabs(
        difference_hier_minus_nonreg ~ age + horizon,
        data = part
      )
      graphics::image(
        as.integer(rownames(z)), as.integer(colnames(z)), z,
        col = grDevices::hcl.colors(31, "Blue-Red 3", rev = TRUE),
        zlim = c(-1, 1),
        xlab = "Age", ylab = "Horizon",
        main = label
      )
    }
    graphics::plot.new()
    graphics::text(
      0.5, 0.55,
      "Poids hierarchique -\npoids non regularise",
      font = 2
    )
  },
  figures_dir, width = 10, height = 7
)

render_chapter_figure(
  "hierarchical_weight_uncertainty",
  function() {
    graphics::par(mfrow = c(2, 3), mar = c(3.3, 3.5, 2.1, 0.8))
    for (model in cfg$models) {
      label <- model_labels[[model]]
      part <- weight_difference[weight_difference$model == label, ]
      z <- xtabs(
        posterior_interval_width ~ age + horizon,
        data = part
      )
      graphics::image(
        as.integer(rownames(z)), as.integer(colnames(z)), z,
        col = grDevices::hcl.colors(25, "YlOrRd", rev = TRUE),
        zlim = c(0, 1),
        xlab = "Age", ylab = "Horizon",
        main = label
      )
    }
    graphics::plot.new()
    graphics::text(
      0.5, 0.55,
      "Largeur de l'intervalle\nposterieur 95 % du poids",
      font = 2
    )
  },
  figures_dir, width = 10, height = 7
)

# ---------------------------------------------------------------------------
# 7. Figures actuarielles completees
# ---------------------------------------------------------------------------

actuarial_survival <- rbind(
  utils::read.csv(
    file.path(
      cfg$paths$metrics,
      "actuarial_survival_summary_force_aggregation_fixed_weights.csv"
    ),
    stringsAsFactors = FALSE
  ),
  utils::read.csv(
    file.path(
      cfg$paths$metrics,
      paste0(
        "actuarial_survival_summary_",
        "force_aggregation_propagated_weights.csv"
      )
    ),
    stringsAsFactors = FALSE
  )
)
actuarial_annuity <- rbind(
  utils::read.csv(
    file.path(
      cfg$paths$metrics,
      "actuarial_annuity_summary_force_aggregation_fixed_weights.csv"
    ),
    stringsAsFactors = FALSE
  ),
  utils::read.csv(
    file.path(
      cfg$paths$metrics,
      paste0(
        "actuarial_annuity_summary_",
        "force_aggregation_propagated_weights.csv"
      )
    ),
    stringsAsFactors = FALSE
  )
)
actuarial_paired <- utils::read.csv(
  file.path(
    cfg$paths$metrics,
    "actuarial_annuity_paired_differences_force_aggregation.csv"
  ),
  stringsAsFactors = FALSE
)
write_spec(actuarial_survival, "46_actuarial_survival_all_rules.csv")
write_spec(actuarial_annuity, "47_actuarial_annuity_all_rules.csv")
write_spec(actuarial_paired, "48_actuarial_paired_differences.csv")
write_spec(
  utils::read.csv(
    file.path(
      cfg$paths$metrics,
      paste0(
        "actuarial_hierarchical_uncertainty_",
        "force_aggregation_propagated_weights.csv"
      )
    ),
    stringsAsFactors = FALSE
  ),
  "51_actuarial_hierarchical_weight_uncertainty.csv"
)
write_spec(
  utils::read.csv(
    file.path(
      cfg$paths$metrics, "actuarial_force_aggregation_checks.csv"
    ),
    stringsAsFactors = FALSE
  ),
  "52_actuarial_force_aggregation_checks.csv"
)

copy_actuarial_figure <- function(source_name, target_name = source_name) {
  for (extension in c("pdf", "png")) {
    source_path <- file.path(
      cfg$paths$figures, paste0(source_name, ".", extension)
    )
    target_path <- file.path(
      figures_dir, paste0(target_name, ".", extension)
    )
    assert_true(
      file.exists(source_path),
      paste("Figure actuarielle corrigee absente :", source_path)
    )
    assert_true(
      file.copy(source_path, target_path, overwrite = TRUE),
      paste("Impossible de copier la figure actuarielle :", target_path)
    )
  }
}
copy_actuarial_figure(
  "actuarial_survival_force_aggregation_fixed_and_propagated",
  "actuarial_primary_survival"
)
copy_actuarial_figure(
  "actuarial_annuity_distributions_force_aggregation"
)
copy_actuarial_figure(
  "actuarial_annuity_paired_differences_force_aggregation"
)
copy_actuarial_figure(
  "actuarial_hierarchical_weight_uncertainty_force_aggregation"
)

render_chapter_figure(
  "actuarial_annuity_interest_sensitivity",
  function() {
    chapter_plot_theme()
    selected <- c("stacking_global", "stacking_hierarchical")
    part <- actuarial_annuity[
      actuarial_annuity$rule == "within_lfo" &
        actuarial_annuity$weight_mode %in%
          c("fixed_weights", "posterior_mean_weights") &
        actuarial_annuity$method %in% selected,
    ]
    graphics::plot(
      range(part$initial_age), range(part$mean),
      type = "n", xlab = "Age initial",
      ylab = "Valeur actuarielle",
      main = "Rentes temporaires h<=10 selon le taux d'actualisation"
    )
    line_types <- c("0.01" = 1, "0.02" = 2, "0.03" = 3)
    for (method in selected) {
      for (rate in sort(unique(part$discount_rate))) {
        piece <- part[
          part$method == method &
            abs(part$discount_rate - rate) < 1e-12,
        ]
        graphics::lines(
          piece$initial_age, piece$mean,
          col = colors[[method]],
          lty = line_types[[sprintf("%.2f", rate)]],
          lwd = 2, type = "b", pch = 16
        )
      }
    }
    graphics::legend(
      "bottomleft",
      legend = c(
        "Global 1/2/3 %", "Hierarchique 1/2/3 %"
      ),
      col = colors[selected], lwd = 2, bty = "n"
    )
  },
  figures_dir, width = 9, height = 5.5
)

message(
  "Complements du cahier des charges prepares dans ", results_dir,
  " ; figures dans ", figures_dir, "."
)
