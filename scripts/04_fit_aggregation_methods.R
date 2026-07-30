#!/usr/bin/env Rscript

# Etape 4 : BMA, pseudo-BMA et trois formes de stacking.

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
source(file.path("R", "model_fitting.R"))
source(file.path("R", "forecasting.R"))
source(file.path("R", "aggregation.R"))
source(file.path("R", "plotting.R"))
cfg <- load_config()
require_stan_backend()
require_packages("bridgesampling")

meta_path <- file.path(cfg$paths$processed, "lfo_meta_data.rds")
assert_true(file.exists(meta_path),
            "Executer d'abord scripts/03_lfo_validation_loop.R.")
meta_object <- readRDS(meta_path)
meta <- meta_object$data
constants <- meta_object$standardization

pseudo_bma <- fit_pseudo_bma(meta, cfg$models)
global <- fit_global_stacking(
  meta, cfg$models,
  multistarts = cfg$contextual_multistarts,
  seed = stable_seed(cfg$seed, "global")
)
contextual <- fit_contextual_stacking(
  meta, cfg$models,
  multistarts = cfg$contextual_multistarts,
  seed = stable_seed(cfg$seed, "contextual")
)
hierarchical_result <- fit_hierarchical_stacking(meta, cfg$models, cfg)

fit_context <- paste0("test_training_", cfg$validation_end)
processed <- readRDS(file.path(
  cfg$paths$processed, "mortality_data.rds"
))
training <- subset_training_data(processed, cfg$validation_end)
individual_fits <- setNames(vector("list", length(cfg$models)), cfg$models)
for (model in cfg$models) {
  individual_fits[[model]] <- read_model_fit(
    model, training, cfg, fit_context
  )
}
bma <- if (cfg$compute_bma) {
  compute_bma_weights(
    individual_fits, cfg$models, cfg, training = training
  )
} else {
  list(
    available = FALSE,
    weights = setNames(rep(NA_real_, length(cfg$models)), cfg$models),
    log_marginal = setNames(rep(NA_real_, length(cfg$models)), cfg$models),
    bridges = list(),
    error = "BMA desactive dans la configuration."
  )
}

max_horizon <- max(
  cfg$lfo_horizons,
  cfg$test_end - cfg$validation_end,
  cfg$annuity_duration
)
grid <- expand.grid(
  age = seq.int(cfg$age_min, cfg$age_max),
  horizon = seq_len(max_horizon)
)
nonregularized_grid <- contextual_weight_grid(
  grid, contextual, cfg$models, constants
)
nonregularized_long <- do.call(rbind, lapply(cfg$models, function(model) {
  data.frame(
    age = nonregularized_grid$age,
    horizon = nonregularized_grid$horizon,
    model = model,
    weight_mean = nonregularized_grid[[model]],
    weight_q05 = NA_real_,
    weight_q50 = nonregularized_grid[[model]],
    weight_q95 = NA_real_,
    method = "stacking_contextual",
    stringsAsFactors = FALSE
  )
}))
hierarchical_long <- hierarchical_weight_grid(
  grid,
  hierarchical_result$fit,
  cfg$models,
  constants,
  ndraws = cfg$forecast_draws,
  seed = stable_seed(cfg$seed, "hierarchical_weights")
)
hierarchical_long$method <- "stacking_hierarchical"
contextual_weights <- rbind(nonregularized_long, hierarchical_long)

global_rows <- rbind(
  data.frame(
    method = "pseudo_bma", model = cfg$models,
    weight = as.numeric(pseudo_bma$weights)
  ),
  data.frame(
    method = "stacking_global", model = cfg$models,
    weight = as.numeric(global$weights)
  ),
  data.frame(
    method = "bma", model = cfg$models,
    weight = as.numeric(bma$weights)
  )
)
aggregation <- list(
  standardization = constants,
  pseudo_bma = pseudo_bma,
  global = global,
  contextual = contextual,
  hierarchical_fit_path = fit_cache_path(file.path(
    cfg$paths$fits, "aggregation", "hierarchical"
  )),
  hierarchical_weights = hierarchical_long,
  bma = bma
)
save_rds_atomic(
  aggregation,
  file.path(cfg$paths$weights, "aggregation_results.rds"),
  compress = FALSE
)
write_csv_atomic(
  global_rows,
  file.path(cfg$paths$weights, "weights_global.csv")
)
write_csv_atomic(
  contextual_weights,
  file.path(cfg$paths$weights, "weights_contextual.csv")
)
write_csv_atomic(
  as.data.frame(contextual$coefficients),
  file.path(cfg$paths$weights, "contextual_coefficients.csv"),
  row.names = TRUE
)
plot_contextual_weights(
  hierarchical_long,
  file.path(cfg$paths$figures, "hierarchical_weights_by_age.png")
)
message("Methodes d'agregation estimees.")
