#!/usr/bin/env Rscript

# Outil de validation ciblee : ajuste un seul modele sur la periode de test.
# L'ajustement utilise exactement le meme cache que l'etape 02.

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
source(file.path("R", "model_fitting.R"))

args <- commandArgs(trailingOnly = TRUE)
model_arg <- grep("^--model=", args, value = TRUE)
assert_true(length(model_arg) == 1L, "Utiliser --model=lc|rh|apc|cbd|m6")
model <- sub("^--model=", "", model_arg[[1L]])
profile <- parse_profile_argument(args)
Sys.setenv(MEMOIRE_PROFILE = profile)
cfg <- load_config(profile)
assert_true(model %in% cfg$models, paste("Modele inconnu :", model))
require_stan_backend()

processed <- readRDS(file.path(
  cfg$paths$processed, "mortality_data.rds"
))
training <- subset_training_data(processed, cfg$validation_end)
compiled <- compile_stan_model(model, cfg)
result <- fit_mortality_model(
  model,
  training,
  compiled,
  cfg,
  context = paste0("test_training_", cfg$validation_end),
  strict = FALSE
)
print(result$diagnostics$overview)
if (!result$diagnostics$pass) quit(save = "no", status = 2L)
