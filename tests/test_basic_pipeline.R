#!/usr/bin/env Rscript

# Contrôle déterministe et rapide de l'étude de simulation.
# Aucun ajustement Stan et aucune donnée HMD ne sont nécessaires.

source(file.path("R", "utils.R"))
source(file.path("R", "forecasting.R"))
source(file.path("R", "aggregation.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "simulation_study.R"))

cfg <- load_config("full")

stopifnot(
  cfg$simulation$repetitions == 30L,
  cfg$simulation$pilot_repetitions == 5L,
  identical(
    cfg$simulation$scenarios,
    c(
      "constant_weights",
      "horizon_only",
      "age_horizon",
      "low_information"
    )
  )
)

for (scenario in cfg$simulation$scenarios) {
  generated <- simulate_study_dataset(scenario, 123L, cfg)
  stopifnot(
    nrow(generated$data) == 410L,
    ncol(generated$true_weights) == 5L,
    all(is.finite(meta_logp_matrix(generated$data, cfg$models))),
    max(abs(rowSums(generated$true_weights) - 1)) < 1e-10
  )
}

pilot_seed <- simulation_repetition_seed(
  cfg, "constant_weights", 1L, "validation", cohort = "pilot"
)
main_seed <- simulation_repetition_seed(
  cfg, "constant_weights", 1L, "validation", cohort = "main"
)
stopifnot(pilot_seed != main_seed)

message("Contrôle de base réussi : configuration, grille, poids et graines.")
