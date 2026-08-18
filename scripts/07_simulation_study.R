#!/usr/bin/env Rscript

# Etape 7 : simulation du stacking contextuel hierarchique sous Stan.
#
# Exemples :
#   Rscript scripts/07_simulation_study.R --stage=pilot --profile=full
#   Rscript scripts/07_simulation_study.R --stage=reference --profile=full
#   Rscript scripts/07_simulation_study.R --stage=main --profile=full
#   Rscript scripts/07_simulation_study.R --stage=validation --profile=full
#   Rscript scripts/07_simulation_study.R --stage=finalize --profile=full
#
# Le pilote et l'etude principale sont deux cohortes independantes. Les
# 20 ajustements pilotes ne sont jamais inclus dans les 120 ajustements
# principaux, meme si leurs numeros de repetition sont identiques.

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
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
allowed_stages <- c(
  "reference", "pilot", "main", "validation", "finalize"
)
assert_true(
  stage %in% allowed_stages,
  paste(
    "Etape de simulation inconnue :", stage,
    "- utiliser pilot, main, validation ou finalize."
  )
)

cfg <- load_config(profile)
stan_backend_name <- require_stan_backend()
message_step(
  "Simulation Stan independante - profil ", cfg$profile,
  ", etape ", stage, ", backend ", stan_backend_name
)

if (stage == "reference") {
  reference <- load_or_build_simulation_reference(cfg)
  message(
    "Calibration predictive HMD terminee : ",
    length(reference$force_draws), " modeles, ",
    nrow(reference$force_draws[[1L]]), " tirages par modele."
  )
  quit(save = "no", status = 0L)
}

if (stage == "finalize") {
  summary <- summarize_simulation_main_results(cfg)
  message(
    "Finalisation de la simulation terminee : ",
    nrow(summary$diagnostics), " repetitions diagnostiquees."
  )
  quit(save = "no", status = 0L)
}

if (stage == "validation") {
  simulation_reference <- load_or_build_simulation_reference(cfg)
  compiled_hierarchical <- compile_stan_model("hierarchical", cfg)
  validation <- run_simulation_validation(
    cfg, compiled_hierarchical, reference = simulation_reference
  )
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
scenario_argument <- grep("^--scenario=", args, value = TRUE)
selected_scenarios <- if (length(scenario_argument)) {
  sub("^--scenario=", "", scenario_argument[[1L]])
} else {
  cfg$simulation$scenarios
}
assert_true(
  all(selected_scenarios %in% cfg$simulation$scenarios),
  paste(
    "Scenario de controle inconnu :",
    paste(setdiff(selected_scenarios, cfg$simulation$scenarios), collapse = ",")
  )
)
repetition_argument <- grep("^--repetition=", args, value = TRUE)
if (length(repetition_argument)) {
  selected_repetition <- suppressWarnings(as.integer(sub(
    "^--repetition=", "", repetition_argument[[1L]]
  )))
  assert_true(
    is.finite(selected_repetition) && selected_repetition %in% repetitions,
    "Numero de repetition de controle invalide."
  )
  repetitions <- selected_repetition
}
partial_run <-
  length(selected_scenarios) != length(cfg$simulation$scenarios) ||
  length(repetitions) != if (stage == "pilot") {
    cfg$simulation$pilot_repetitions
  } else {
    cfg$simulation$repetitions
  }
predictive_draws <- if (
    stage == "pilot" && cfg$profile == "smoke") {
  cfg$simulation$test_predictive_draws
} else {
  cfg$simulation$predictive_draws
}
simulation_reference <- load_or_build_simulation_reference(cfg)

# Un ajustement occupe deja les quatre coeurs physiques avec ses quatre
# chaines. Les repetitions sont donc enchainees, ce qui evite de lancer des
# chaines et des repetitions concurrentes sur les memes coeurs.
compiled_hierarchical <- compile_stan_model("hierarchical", cfg)
for (scenario in selected_scenarios) {
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
      cohort = stage,
      reference = simulation_reference
    )
    message_step(
      "Fin ", stage, " : ", scenario,
      ", repetition ", repetition,
      ", statut ", result$status
    )
  }
}

if (partial_run) {
  message(
    "Controle partiel termine : scenarios ",
    paste(selected_scenarios, collapse = ","),
    ", repetitions ", paste(repetitions, collapse = ","), "."
  )
  quit(save = "no", status = 0L)
}

stage_results <- write_simulation_stage_results(
  cfg,
  repetitions,
  stage
)
if (stage == "pilot") {
  assessment_path <- file.path(
    simulation_output_root(cfg), "pilot", "pilot_assessment.csv"
  )
  message(
    "Etude pilote terminee : ",
    nrow(stage_results$diagnostics),
    " ajustements sur ",
    length(cfg$simulation$scenarios) * length(repetitions),
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
