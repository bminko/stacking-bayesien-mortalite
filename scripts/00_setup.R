#!/usr/bin/env Rscript

# Installe les dependances R dans le projet et verifie CmdStan.

source(file.path("R", "utils.R"))
root <- find_project_root()
library_path <- activate_project_library(root)
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(library_path, .libPaths())))

repos <- c(
  CRAN = "https://cloud.r-project.org",
  Stan = "https://stan-dev.r-universe.dev"
)
needed <- c("posterior", "bridgesampling")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installation des packages : ", paste(missing, collapse = ", "))
  utils::install.packages(missing, repos = repos, lib = library_path)
}
require_packages(needed)

if (requireNamespace("rstan", quietly = TRUE)) {
  message("RStan est deja disponible : aucune installation CmdStan requise.")
  message("Bibliotheque R du projet : ", library_path)
  quit(save = "no", status = 0L)
}

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  utils::install.packages("cmdstanr", repos = repos, lib = library_path)
}
require_packages("cmdstanr")

cmdstan_available <- tryCatch(
  nzchar(cmdstanr::cmdstan_path()),
  error = function(error) FALSE
)
if (!cmdstan_available) {
  message("Verification de la chaine de compilation CmdStan...")
  try(cmdstanr::check_cmdstan_toolchain(fix = TRUE), silent = FALSE)
  tools_dir <- file.path(root, "tools")
  dir.create(tools_dir, recursive = TRUE, showWarnings = FALSE)
  cmdstanr::install_cmdstan(
    dir = tools_dir,
    cores = max(1L, parallel::detectCores(logical = FALSE) - 1L)
  )
}

message("CmdStan : ", cmdstanr::cmdstan_path())
message("Bibliotheque R du projet : ", library_path)
message("Installation terminee.")
