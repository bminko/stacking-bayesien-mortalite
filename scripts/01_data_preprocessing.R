#!/usr/bin/env Rscript

# Etape 1 : lecture, controles stricts et sauvegarde des donnees nettoyees.

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
cfg <- load_config()

message_step("Pretraitement HMD : ", cfg$year_start, "-", cfg$data_end,
             ", ages ", cfg$age_min, "-", cfg$age_max,
             ", population ", cfg$sex)
processed <- prepare_hmd_data(cfg)
report <- data_quality_report(processed, cfg)

save_rds_atomic(
  processed,
  file.path(cfg$paths$processed, "mortality_data.rds")
)
write_csv_atomic(
  processed$long,
  file.path(cfg$paths$processed, "mortality_long.csv")
)
write_csv_atomic(
  report,
  file.path(cfg$paths$processed, "data_quality_report.csv")
)
message("Donnees valides : ", nrow(processed$long), " cellules.")
