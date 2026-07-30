#!/usr/bin/env Rscript

# Etape 2 : ajustement des cinq modeles sur la periode precedant le test.

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
source(file.path("R", "model_fitting.R"))
cfg <- load_config()
require_stan_backend()

processed_path <- file.path(cfg$paths$processed, "mortality_data.rds")
assert_true(file.exists(processed_path),
            "Executer d'abord scripts/01_data_preprocessing.R.")
processed <- readRDS(processed_path)
training <- subset_training_data(processed, cfg$validation_end)

diagnostics <- vector("list", length(cfg$models))
names(diagnostics) <- cfg$models
for (model in cfg$models) {
  message_step("Compilation/chargement du modele ", toupper(model))
  compiled <- compile_stan_model(model, cfg)
  result <- fit_mortality_model(
    model = model,
    training = training,
    compiled_model = compiled,
    cfg = cfg,
    context = paste0("test_training_", cfg$validation_end)
  )
  diagnostics[[model]] <- result$diagnostics$overview
}

overview <- do.call(rbind, diagnostics)
rownames(overview) <- NULL
write_csv_atomic(
  overview,
  file.path(cfg$paths$diagnostics, "individual_models_test_training.csv")
)
message("Ajustements individuels termines.")
