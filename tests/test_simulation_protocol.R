#!/usr/bin/env Rscript

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
source(file.path("R", "forecasting.R"))
source(file.path("R", "aggregation.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "simulation_study.R"))

cfg <- load_config("full")
stopifnot(
  cfg$simulation$repetitions == 30L,
  cfg$simulation$protocol_version ==
    "stan_hmd_calibrated_repetitions_v4",
  cfg$simulation$calibration_end == 2015L,
  cfg$simulation$exposure_reference_year == 2015L,
  cfg$simulation$low_information_exposure_scale == 0.10,
  cfg$simulation$pilot_repetitions == 5L,
  cfg$simulation$light_mcmc$iter_warmup == 500L,
  cfg$simulation$light_mcmc$iter_sampling == 500L,
  cfg$simulation$predictive_draws == 1000L,
  cfg$simulation$full_mcmc$iter_warmup == 1000L,
  cfg$simulation$full_mcmc$iter_sampling == 1000L,
  cfg$simulation$validation_thresholds$
    mean_weight_absolute_difference == 0.02
)

processed <- prepare_hmd_data(cfg)
reference_target <- simulation_reference_target(processed, cfg)
observed_2015 <- processed$long[
  processed$long$year == cfg$simulation$exposure_reference_year,
  c("age", "deaths", "exposure")
]
age_index <- match(reference_target$age, observed_2015$age)
observed_rate <-
  observed_2015$deaths[age_index] / observed_2015$exposure[age_index]
S <- 30L
model_factors <- setNames(
  c(0.92, 0.98, 1.00, 1.04, 1.08), cfg$models
)
test_force_draws <- setNames(lapply(cfg$models, function(model) {
  center <- observed_rate * exp(-0.01 * reference_target$horizon) *
    model_factors[[model]]
  outer(exp(seq(-0.03, 0.03, length.out = S)), center)
}), cfg$models)
test_phi_draws <- setNames(lapply(
  seq_along(cfg$models),
  function(index) rep(100 + 10 * index, S)
), cfg$models)
test_reference <- list(
  cache_key = "unit_test_hmd_reference",
  target = reference_target,
  force_draws = test_force_draws,
  phi_draws = test_phi_draws
)
stopifnot(
  nrow(reference_target) == 410L,
  identical(sort(unique(reference_target$year)), 2016:2025),
  length(unique(reference_target$exposure)) > 30L,
  all(reference_target$exposure > 0)
)

pilot_seed <- simulation_repetition_seed(
  cfg, "constant_weights", 1L, "validation", cohort = "pilot"
)
main_seed <- simulation_repetition_seed(
  cfg, "constant_weights", 1L, "validation", cohort = "main"
)
stopifnot(
  pilot_seed != main_seed,
  grepl(
    "pilot_repetitions",
    simulation_repetition_dir(
      cfg, "constant_weights", 1L, cohort = "pilot"
    ),
    fixed = TRUE
  ),
  grepl(
    "repetitions",
    simulation_repetition_dir(
      cfg, "constant_weights", 1L, cohort = "main"
    ),
    fixed = TRUE
  ),
  simulation_fit_context(
    "constant_weights", 1L, "legere", cohort = "pilot"
  ) != simulation_fit_context(
    "constant_weights", 1L, "legere", cohort = "main"
  )
)

for (scenario in cfg$simulation$scenarios) {
  generated <- simulate_study_dataset(
    scenario, 123L, cfg, reference = test_reference
  )
  expected_scale <- if (scenario == "low_information") 0.10 else 1
  stopifnot(
    nrow(generated$data) == 410L,
    all(abs(rowSums(generated$true_weights) - 1) < 1e-10),
    all(is.finite(meta_logp_matrix(generated$data, cfg$models))),
    all(generated$data$omega == 1),
    all(generated$data$deaths >= 0),
    all(generated$data$exposure ==
          reference_target$exposure * expected_scale),
    all(generated$selected_models %in% cfg$models),
    generated$reference_cache_key == test_reference$cache_key
  )
}

main_diagnostics <- file.path(
  simulation_output_root(cfg), "main", "diagnostics_by_repetition.csv"
)
if (file.exists(main_diagnostics)) {
  validation_selection <- select_simulation_validation_repetitions(cfg)
  expected_counts <- c(
    constant_weights = 2L,
    horizon_only = 2L,
    age_horizon = 2L,
    low_information = 4L
  )
  observed_counts <- table(validation_selection$scenario)
  stopifnot(
    nrow(validation_selection) == 10L,
    all(observed_counts[names(expected_counts)] == expected_counts),
    any(validation_selection$prior_retry),
    any(grepl(
      "^diagnostic_|^relance_",
      validation_selection$selection_reason
    )),
    any(grepl(
      "^couverture_rmse_",
      validation_selection$selection_reason
    ))
  )
}

message("Protocole et DGP de simulation : OK.")
