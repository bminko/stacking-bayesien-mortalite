#!/usr/bin/env Rscript

# Etape 5 : evaluation hors echantillon depuis l'origine fixe 2015.
#
# Les modeles et les poids sont ceux ajustes avant le test. L'ajout de
# nouvelles annees ne doit donc relancer ni Stan ni les aggregations.

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
source(file.path("R", "model_fitting.R"))
source(file.path("R", "forecasting.R"))
source(file.path("R", "aggregation.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "plotting.R"))
cfg <- load_config()
require_stan_backend()

processed <- readRDS(file.path(
  cfg$paths$processed, "mortality_data.rds"
))
aggregation_path <- file.path(
  cfg$paths$weights, "aggregation_results.rds"
)
assert_true(file.exists(aggregation_path),
            "Executer d'abord scripts/04_fit_aggregation_methods.R.")
aggregation <- readRDS(aggregation_path)

training <- subset_training_data(processed, cfg$validation_end)
test_horizons <- seq_len(cfg$test_end - cfg$validation_end)
target <- target_grid_from_observed(
  processed, cfg$validation_end, test_horizons, cfg$test_end
)
fit_context <- paste0("test_training_", cfg$validation_end)
predictions <- setNames(vector("list", length(cfg$models)), cfg$models)

for (model in cfg$models) {
  fit <- read_model_fit(model, training, cfg, fit_context)
  message_step("Previsions de test - ", toupper(model))
  predictions[[model]] <- forecast_model(
    fit, model, training, target, cfg,
    keep_draws = TRUE,
    seed = stable_seed(cfg$seed, cfg$profile, "test", model)
  )
}

evaluation <- evaluate_predictions(
  predictions, aggregation, cfg$models, cfg
)

primary_methods <- c(
  cfg$models,
  "stacking_global",
  "stacking_contextual",
  "stacking_hierarchical"
)
periods <- list(
  test_initial_2016_2020 = 2016:2020,
  extension_2021_2024 = 2021:2024,
  test_extended_2016_2024 = 2016:2024
)
sensitivity_periods <- list(
  test_extended_2016_2024 = 2016:2024,
  without_2020_2022 = c(2016:2019, 2023:2024)
)

# Controle de non-regression avant d'ecraser les sorties canoniques.
checks <- list()
record_check <- function(check, pass, details) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = check,
    pass = isTRUE(pass),
    details = as.character(details),
    stringsAsFactors = FALSE
  )
}

reference_dir <- file.path(
  cfg$paths$results, "reference_2016_2020"
)
reference_cells_path <- file.path(
  reference_dir, "test_metrics_cells.csv"
)
assert_true(
  file.exists(reference_cells_path),
  paste("Reference 2016-2020 absente :", reference_cells_path)
)
reference_cells <- utils::read.csv(
  reference_cells_path, stringsAsFactors = FALSE
)
extended_initial <- evaluation$cells[
  evaluation$cells$year %in% 2016:2020,
  names(reference_cells),
  drop = FALSE
]
key <- c("method", "year", "age", "horizon")
order_reference <- do.call(order, unname(reference_cells[key]))
order_extended <- do.call(order, unname(extended_initial[key]))
reference_cells <- reference_cells[order_reference, , drop = FALSE]
extended_initial <- extended_initial[order_extended, , drop = FALSE]
rownames(reference_cells) <- NULL
rownames(extended_initial) <- NULL
same_keys <- all(vapply(key, function(column) {
  identical(reference_cells[[column]], extended_initial[[column]])
}, logical(1)))
numeric_columns <- names(reference_cells)[vapply(
  reference_cells, is.numeric, logical(1)
)]
metric_differences <- vapply(numeric_columns, function(column) {
  max(
    abs(reference_cells[[column]] - extended_initial[[column]]),
    na.rm = TRUE
  )
}, numeric(1))
max_metric_difference <- max(metric_differences)
record_check(
  "resultats_2016_2020_inchanges",
  same_keys && is.finite(max_metric_difference) &&
    max_metric_difference <= 1e-12,
  sprintf("ecart numerique maximal = %.17g", max_metric_difference)
)

immutable_files <- c(
  "aggregation_results.rds",
  "weights_global.csv",
  "weights_contextual.csv",
  "contextual_coefficients.csv"
)
for (filename in immutable_files) {
  current_path <- file.path(cfg$paths$weights, filename)
  reference_path <- file.path(reference_dir, filename)
  same_hash <- identical(
    unname(tools::md5sum(current_path)),
    unname(tools::md5sum(reference_path))
  )
  record_check(
    paste0("fichier_de_poids_inchange_", filename),
    same_hash,
    paste("MD5 =", unname(tools::md5sum(current_path)))
  )
}

meta_object <- readRDS(file.path(
  cfg$paths$processed, "lfo_meta_data.rds"
))
meta <- meta_object$data
record_check(
  "validation_limitee_a_2001_2015",
  min(meta$year) == 2001L && max(meta$year) == 2015L &&
    !any(meta$year > 2015L),
  paste(min(meta$year), max(meta$year), sep = "-")
)
record_check(
  "origines_lfo_2000_2014",
  identical(sort(unique(meta$origin)), 2000:2014),
  paste(sort(unique(meta$origin)), collapse = ",")
)
record_check(
  "horizons_lfo_1_10",
  identical(sort(unique(meta$horizon)), 1:10),
  paste(sort(unique(meta$horizon)), collapse = ",")
)
record_check(
  "standardisation_inchangee",
  identical(meta_object$standardization, aggregation$standardization),
  paste(
    "age_mean", aggregation$standardization$age_mean,
    "horizon_mean", aggregation$standardization$horizon_mean
  )
)
record_check(
  "apprentissage_termine_en_2015",
  max(training$years) == 2015L,
  paste(min(training$years), max(training$years), sep = "-")
)
record_check(
  "origine_de_test_fixe_2015",
  all(target$origin == 2015L),
  paste(unique(target$origin), collapse = ",")
)
record_check(
  "correspondance_annee_horizon",
  all(target$year - 2015L == target$horizon) &&
    identical(sort(unique(target$horizon)), 1:9),
  paste(
    min(target$year), max(target$year), "-> horizons",
    min(target$horizon), max(target$horizon)
  )
)
record_check(
  "expositions_observees_utilisees",
  all(is.finite(target$exposure)) && all(target$exposure > 0),
  paste("cellules =", nrow(target))
)

draw_dimensions <- lapply(predictions, function(prediction) {
  dim(prediction$count_draws)
})
same_draw_dimensions <- length(unique(vapply(
  draw_dimensions, paste, character(1), collapse = "x"
))) == 1L
record_check(
  "dimensions_predictives_identiques",
  same_draw_dimensions &&
    identical(draw_dimensions[[1L]], c(cfg$forecast_draws, nrow(target))),
  paste(draw_dimensions[[1L]], collapse = "x")
)
finite_log_predictive <- all(vapply(
  predictions,
  function(prediction) {
    all(is.finite(prediction$summary$log_predictive))
  },
  logical(1)
))
record_check(
  "log_densites_finies",
  finite_log_predictive,
  paste("modeles =", paste(cfg$models, collapse = ","))
)

weight_sums <- unlist(lapply(
  evaluation$weight_sets[primary_methods],
  rowSums
))
weight_values <- unlist(evaluation$weight_sets[primary_methods])
record_check(
  "poids_non_negatifs",
  all(is.finite(weight_values)) && all(weight_values >= 0),
  sprintf("minimum = %.17g", min(weight_values))
)
record_check(
  "poids_somment_a_un",
  max(abs(weight_sums - 1)) <= 1e-12,
  sprintf("ecart maximal = %.17g", max(abs(weight_sums - 1)))
)

audit <- do.call(rbind, checks)
rownames(audit) <- NULL
write_csv_atomic(
  audit,
  file.path(cfg$paths$diagnostics, "extension_2024_audit_checks.csv")
)
assert_true(
  all(audit$pass),
  paste(
    "Extension 2024 rejetee :",
    paste(audit$check[!audit$pass], collapse = ", ")
  )
)

# L'audit est passe : les sorties canoniques peuvent etre remplacees.
save_rds_atomic(
  predictions,
  file.path(cfg$paths$metrics, "test_predictions.rds"),
  compress = FALSE
)
save_rds_atomic(
  evaluation,
  file.path(cfg$paths$metrics, "test_evaluation.rds"),
  compress = FALSE
)
write_csv_atomic(
  evaluation$cells,
  file.path(cfg$paths$metrics, "test_metrics_cells.csv")
)
write_csv_atomic(
  evaluation$overall,
  file.path(cfg$paths$metrics, "test_metrics_overall.csv")
)
write_csv_atomic(
  evaluation$by_horizon,
  file.path(cfg$paths$metrics, "test_metrics_by_horizon.csv")
)
write_csv_atomic(
  evaluation$by_age,
  file.path(cfg$paths$metrics, "test_metrics_by_age.csv")
)
write_csv_atomic(
  evaluation$by_age_group,
  file.path(cfg$paths$metrics, "test_metrics_by_age_group.csv")
)

primary_cells <- evaluation$cells[
  evaluation$cells$method %in% primary_methods,
  ,
  drop = FALSE
]
period_summaries <- summarize_evaluation_periods(
  primary_cells, periods
)
sensitivity_summaries <- summarize_evaluation_periods(
  primary_cells, sensitivity_periods
)
write_csv_atomic(
  primary_cells,
  file.path(
    cfg$paths$metrics, "test_metrics_cells_2016_2024.csv"
  )
)
write_csv_atomic(
  period_summaries$overall,
  file.path(
    cfg$paths$metrics, "test_metrics_overall_by_period.csv"
  )
)
write_csv_atomic(
  period_summaries$by_age,
  file.path(
    cfg$paths$metrics, "test_metrics_by_age_and_period.csv"
  )
)
write_csv_atomic(
  period_summaries$by_horizon,
  file.path(
    cfg$paths$metrics, "test_metrics_by_horizon_and_period.csv"
  )
)
write_csv_atomic(
  period_summaries$by_age_group,
  file.path(
    cfg$paths$metrics, "test_metrics_by_age_group_and_period.csv"
  )
)
write_csv_atomic(
  sensitivity_summaries$overall,
  file.path(
    cfg$paths$metrics, "sensitivity_metrics_overall.csv"
  )
)
for (period in names(periods)) {
  period_rows <- period_summaries$overall[
    period_summaries$overall$period == period,
    ,
    drop = FALSE
  ]
  write_csv_atomic(
    period_rows,
    file.path(
      cfg$paths$metrics,
      paste0("test_metrics_overall_", period, ".csv")
    )
  )
}

stacking_methods <- c(
  "stacking_global",
  "stacking_contextual",
  "stacking_hierarchical"
)
weight_target <- unique(target[c("age", "horizon")])
weight_long <- evaluation_weight_long(
  weight_target,
  evaluation$weight_sets,
  stacking_methods,
  cfg$models
)
weight_by_horizon <- summarize_weights_by_horizon(weight_long)
write_csv_atomic(
  weight_long,
  file.path(
    cfg$paths$weights, "stacking_weights_age_horizon_1_9.csv"
  )
)
write_csv_atomic(
  weight_by_horizon,
  file.path(
    cfg$paths$weights, "stacking_weights_by_horizon_1_9.csv"
  )
)

plot_metric_by_horizon(
  period_summaries$by_horizon[
    period_summaries$by_horizon$period ==
      "test_extended_2016_2024",
  ], "logs",
  file.path(cfg$paths$figures, "test_logs_by_horizon.png"),
  "Score logarithmique moyen"
)
plot_metric_by_horizon(
  period_summaries$by_horizon[
    period_summaries$by_horizon$period ==
      "test_extended_2016_2024",
  ], "crps",
  file.path(cfg$paths$figures, "test_crps_by_horizon.png"),
  "CRPS moyen sur les taux"
)
for (method in unique(evaluation$cells$method)) {
  plot_metric_heatmap(
    evaluation$cells, method, "crps",
    file.path(
      cfg$paths$figures,
      paste0("heatmap_crps_", method, ".png")
    )
  )
}
message(
  "Evaluation etendue terminee : origine 2015, ",
  min(target$year), "-", max(target$year), "."
)
