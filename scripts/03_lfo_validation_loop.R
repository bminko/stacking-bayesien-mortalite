#!/usr/bin/env Rscript

# Etape 3 : validation Leave-Future-Out a origines multiples.

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
source(file.path("R", "model_fitting.R"))
source(file.path("R", "forecasting.R"))
source(file.path("R", "aggregation.R"))
cfg <- load_config()
require_stan_backend()

processed_path <- file.path(cfg$paths$processed, "mortality_data.rds")
assert_true(file.exists(processed_path),
            "Executer d'abord scripts/01_data_preprocessing.R.")
processed <- readRDS(processed_path)

compiled_models <- lapply(cfg$models, compile_stan_model, cfg = cfg)
names(compiled_models) <- cfg$models
prediction_rows <- list()
position <- 1L

for (origin in cfg$lfo_origins) {
  target <- target_grid_from_observed(
    processed, origin, cfg$lfo_horizons, cfg$validation_end
  )
  if (is.null(target)) next
  training <- subset_training_data(processed, origin)
  context <- paste0("lfo_origin_", origin)

  for (model in cfg$models) {
    message_step("LFO origine ", origin, " - ", toupper(model))
    fit_result <- fit_mortality_model(
      model, training, compiled_models[[model]], cfg, context,
      mcmc = cfg$lfo_mcmc
    )
    prediction <- forecast_model(
      fit_result$fit, model, training, target, cfg,
      keep_draws = FALSE,
      seed = stable_seed(cfg$seed, "lfo", origin, model)
    )
    row <- prediction$summary[c(
      "origin", "horizon", "year", "age", "deaths",
      "exposure", "model", "log_predictive"
    )]
    prediction_rows[[position]] <- row
    position <- position + 1L
  }
}

long <- do.call(rbind, prediction_rows)
rownames(long) <- NULL
key_columns <- c(
  "origin", "horizon", "year", "age", "deaths", "exposure"
)
wide <- reshape(
  long,
  idvar = key_columns,
  timevar = "model",
  direction = "wide"
)
names(wide) <- sub("^log_predictive[.]", "log_p_", names(wide))
wide <- wide[order(wide$origin, wide$horizon, wide$age), ]
rownames(wide) <- NULL

expected_log_columns <- paste0("log_p_", cfg$models)
assert_true(
  all(expected_log_columns %in% names(wide)),
  "Le meta-jeu LFO ne contient pas les cinq modeles."
)
assert_true(
  all(is.finite(as.matrix(wide[expected_log_columns]))),
  "Le meta-jeu LFO contient des densites non finies."
)

origin_horizon <- unique(wide[c("origin", "horizon", "year")])
multiplicity <- table(origin_horizon$year)
wide$omega <- 1 / as.numeric(multiplicity[as.character(wide$year)])
standardized <- standardize_context(wide)
wide <- standardized$data

save_rds_atomic(
  list(data = wide, standardization = standardized$constants),
  file.path(cfg$paths$processed, "lfo_meta_data.rds")
)
write_csv_atomic(
  wide,
  file.path(cfg$paths$processed, "lfo_meta_data.csv")
)
write_csv_atomic(
  long,
  file.path(cfg$paths$processed, "lfo_predictive_densities_long.csv")
)

diagnostics <- collect_diagnostic_overviews(cfg)
write_csv_atomic(
  diagnostics,
  file.path(cfg$paths$diagnostics, "all_mcmc_diagnostics.csv")
)
message("Meta-jeu LFO termine : ", nrow(wide), " cellules.")
