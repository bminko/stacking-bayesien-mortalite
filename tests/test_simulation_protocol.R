#!/usr/bin/env Rscript

source(file.path("R", "utils.R"))
source(file.path("R", "forecasting.R"))
source(file.path("R", "aggregation.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "simulation_study.R"))

cfg <- load_config("full")
stopifnot(
  cfg$simulation$repetitions == 30L,
  cfg$simulation$pilot_repetitions == 5L,
  cfg$simulation$light_mcmc$iter_warmup == 500L,
  cfg$simulation$light_mcmc$iter_sampling == 500L,
  cfg$simulation$predictive_draws == 1000L,
  cfg$simulation$full_mcmc$iter_warmup == 1000L,
  cfg$simulation$full_mcmc$iter_sampling == 1000L,
  cfg$simulation$validation_thresholds$
    mean_weight_absolute_difference == 0.02
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
  generated <- simulate_study_dataset(scenario, 123L, cfg)
  stopifnot(
    nrow(generated$data) == 410L,
    all(abs(rowSums(generated$true_weights) - 1) < 1e-10),
    all(is.finite(meta_logp_matrix(generated$data, cfg$models))),
    all(generated$data$omega == 1)
  )
}

main_diagnostics <- file.path(
  cfg$paths$simulation, "main", "diagnostics_by_repetition.csv"
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
