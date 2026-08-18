#!/usr/bin/env Rscript

# Lance les six etapes empiriques dans des processus R separes.
# L'etude de simulation (etape 07) se lance explicitement avec --stage.

source(file.path("R", "utils.R"))
args <- commandArgs(trailingOnly = TRUE)
profile <- parse_profile_argument(args)
Sys.setenv(MEMOIRE_PROFILE = profile)

steps_argument <- grep("^--steps=", args, value = TRUE)
default_steps <- sprintf("%02d", 1:6)
steps <- if (length(steps_argument)) {
  strsplit(sub("^--steps=", "", steps_argument[[1L]]), ",")[[1L]]
} else {
  default_steps
}
steps <- sprintf("%02d", as.integer(steps))

scripts <- c(
  "01" = "scripts/01_data_preprocessing.R",
  "02" = "scripts/02_fit_individual_models.R",
  "03" = "scripts/03_lfo_validation_loop.R",
  "04" = "scripts/04_fit_aggregation_methods.R",
  "05" = "scripts/05_evaluate_test_metrics.R",
  "06" = "scripts/06_actuarial_quantities.R"
)
unknown <- setdiff(steps, names(scripts))
assert_true(!length(unknown), paste(
  "Etapes inconnues :", paste(unknown, collapse = ", ")
))

rscript <- file.path(R.home("bin"), "Rscript")
for (step in steps) {
  message_step("Lancement de ", scripts[[step]], " (profil ", profile, ")")
  status <- system2(
    rscript,
    args = scripts[[step]]
  )
  if (!identical(status, 0L)) {
    stop("Echec de l'etape ", step, " (code ", status, ").", call. = FALSE)
  }
}
message("Pipeline termine pour le profil ", profile, ".")
