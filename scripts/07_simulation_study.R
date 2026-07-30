#!/usr/bin/env Rscript

# Etape 7 : simulation du stacking contextuel hierarchique sous Stan.
#
# Exemples :
#   Rscript scripts/07_simulation_study.R --stage=pilot --profile=full
#   Rscript scripts/07_simulation_study.R --stage=main --profile=full
#   Rscript scripts/07_simulation_study.R --stage=main --profile=full \
#     --scenario=age_horizon
#   Rscript scripts/07_simulation_study.R --stage=main --profile=full \
#     --scenario=age_horizon --repetition=1
#   Rscript scripts/07_simulation_study.R --stage=validation --profile=full
#   Rscript scripts/07_simulation_study.R --stage=finalize --profile=full
#
# Le pilote et l'etude principale sont deux cohortes independantes. Les
# 20 ajustements pilotes ne sont jamais inclus dans les 120 ajustements
# principaux, meme si leurs numeros de repetition sont identiques.

source(file.path("R", "utils.R"))
source(file.path("R", "model_fitting.R"))
source(file.path("R", "forecasting.R"))
source(file.path("R", "aggregation.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "reporting.R"))
source(file.path("R", "simulation_study.R"))

args <- commandArgs(trailingOnly = TRUE)
profile <- parse_profile_argument(args)
stage_argument <- grep("^--stage=", args, value = TRUE)
stage <- if (length(stage_argument)) {
  sub("^--stage=", "", stage_argument[[1L]])
} else {
  "pilot"
}
allowed_stages <- c("pilot", "main", "validation", "finalize")
assert_true(
  stage %in% allowed_stages,
  paste(
    "Etape de simulation inconnue :", stage,
    "- utiliser pilot, main, validation ou finalize."
  )
)

cfg <- load_config(profile)
scenario_argument <- grep("^--scenario=", args, value = TRUE)
repetition_argument <- grep("^--repetition=", args, value = TRUE)
scenarios_to_run <- if (length(scenario_argument)) {
  sub("^--scenario=", "", scenario_argument[[1L]])
} else {
  cfg$simulation$scenarios
}
assert_true(
  all(scenarios_to_run %in% cfg$simulation$scenarios),
  paste(
    "Scenario inconnu :",
    paste(setdiff(scenarios_to_run, cfg$simulation$scenarios), collapse = ", ")
  )
)
selected_repetition <- if (length(repetition_argument)) {
  suppressWarnings(as.integer(
    sub("^--repetition=", "", repetition_argument[[1L]])
  ))
} else {
  integer(0)
}
if (length(repetition_argument)) {
  assert_true(
    length(selected_repetition) == 1L &&
      is.finite(selected_repetition) &&
      selected_repetition >= 1L,
    "La repetition doit etre un entier strictement positif."
  )
}
if (stage %in% c("validation", "finalize")) {
  assert_true(
    !length(scenario_argument) && !length(repetition_argument),
    paste(
      "--scenario et --repetition sont reserves aux etapes pilot et main.",
      "La validation et la finalisation utilisent leur selection enregistree."
    )
  )
}
stan_backend_name <- require_stan_backend()
message_step(
  "Simulation Stan independante - profil ", cfg$profile,
  ", etape ", stage, ", backend ", stan_backend_name
)

if (stage == "finalize") {
  summary <- summarize_simulation_main_results(cfg)
  message(
    "Finalisation de la simulation terminee : ",
    nrow(summary$diagnostics), " repetitions diagnostiquees."
  )
  quit(save = "no", status = 0L)
}

if (stage == "validation") {
  compiled_hierarchical <- compile_stan_model("hierarchical", cfg)
  validation <- run_simulation_validation(cfg, compiled_hierarchical)
  message(
    "Validation allegee/complete terminee : ",
    validation$summary$repetitions_valid,
    " repetitions valides sur ",
    validation$summary$repetitions_planned,
    ". Examiner validation_summary.csv avant la finalisation."
  )
  quit(save = "no", status = 0L)
}

repetitions <- if (stage == "pilot") {
  seq_len(cfg$simulation$pilot_repetitions)
} else {
  seq_len(cfg$simulation$repetitions)
}
if (length(selected_repetition)) {
  assert_true(
    selected_repetition <= max(repetitions),
    paste(
      "Repetition hors plage pour l'etape", stage, ":",
      selected_repetition, ">", max(repetitions)
    )
  )
  repetitions <- selected_repetition
}
predictive_draws <- if (
    stage == "pilot" && cfg$profile == "smoke") {
  cfg$simulation$test_predictive_draws
} else {
  cfg$simulation$predictive_draws
}

# Un ajustement occupe deja les quatre coeurs physiques avec ses quatre
# chaines. Les repetitions sont donc enchainees, ce qui evite de lancer des
# chaines et des repetitions concurrentes sur les memes coeurs.
compiled_hierarchical <- compile_stan_model("hierarchical", cfg)
for (scenario in scenarios_to_run) {
  for (repetition in repetitions) {
    message_step(
      "Debut ", stage, " : ", scenario,
      ", repetition ", repetition, "/", max(repetitions)
    )
    result <- run_simulation_repetition(
      scenario,
      repetition,
      compiled_hierarchical,
      cfg,
      predictive_draws = predictive_draws,
      cohort = stage
    )
    message_step(
      "Fin ", stage, " : ", scenario,
      ", repetition ", repetition,
      ", statut ", result$status
    )
  }
}

stage_results <- write_simulation_stage_results(
  cfg,
  repetitions,
  stage,
  scenarios = scenarios_to_run
)
if (stage == "pilot") {
  assessment_path <- file.path(
    cfg$paths$simulation, "pilot", "pilot_assessment.csv"
  )
  message(
    "Etude pilote terminee : ",
    nrow(stage_results$diagnostics),
    " ajustements sur ",
    length(scenarios_to_run) * length(repetitions),
    ". Examiner ", assessment_path,
    " avant de lancer l'etude principale."
  )
} else {
  message(
    "Etude principale terminee : ",
    nrow(stage_results$diagnostics),
    " repetitions diagnostiquees. Lancer ensuite --stage=finalize."
  )
}
