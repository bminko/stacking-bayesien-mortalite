# Fonctions generales sans dependance externe.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

is_rstan_fit <- function(x) {
  # `inherits()` essaie de charger la definition S4 de RStan. La lecture
  # seule des ajustements archives doit aussi fonctionner sans ce package.
  stored_class <- attr(x, "class", exact = TRUE)
  any(as.character(stored_class) == "stanfit")
}

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    code_markers <- c(
      file.path(current, "config", "config.R"),
      file.path(current, "scripts", "00_run_pipeline.R")
    )
    legacy_marker <- file.path(
      current, "Plan_de_travail_Memoire_Minko.md"
    )
    if (all(file.exists(code_markers)) || file.exists(legacy_marker)) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Racine du projet introuvable depuis : ", start)
    }
    current <- parent
  }
}

load_config <- function(profile = Sys.getenv("MEMOIRE_PROFILE", "full")) {
  root <- find_project_root()
  activate_project_library(root)
  env <- new.env(parent = environment())
  sys.source(file.path(root, "config", "config.R"), envir = env)
  cfg <- env$minkos_config(profile = profile, root = root)
  ensure_project_dirs(cfg)
  cfg
}

activate_project_library <- function(root = find_project_root()) {
  version <- paste(R.version$major, strsplit(R.version$minor, "[.]")[[1L]][1L],
                   sep = ".")
  library_path <- file.path(root, "R", "library", version)
  if (dir.exists(library_path)) {
    .libPaths(unique(c(library_path, .libPaths())))
  }
  invisible(library_path)
}

ensure_project_dirs <- function(cfg) {
  dirs <- unique(unlist(cfg$paths[c(
    "processed", "results", "fits", "diagnostics", "weights",
    "metrics", "figures", "simulation"
  )], use.names = FALSE))
  invisible(vapply(dirs, dir.create, logical(1), recursive = TRUE,
                   showWarnings = FALSE))
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1),
                              quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Packages R manquants : ", paste(missing, collapse = ", "),
      ". Executer d'abord scripts/00_setup.R.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

stan_backend <- function() {
  requested <- tolower(Sys.getenv("MEMOIRE_STAN_BACKEND", "auto"))
  if (!requested %in% c("auto", "cmdstanr", "rstan")) {
    stop("MEMOIRE_STAN_BACKEND doit valoir auto, cmdstanr ou rstan.")
  }

  cmdstan_ready <- FALSE
  if (requested %in% c("auto", "cmdstanr") &&
      requireNamespace("cmdstanr", quietly = TRUE)) {
    cmdstan_ready <- tryCatch(
      nzchar(cmdstanr::cmdstan_path()),
      error = function(error) FALSE
    )
  }
  if (cmdstan_ready) return("cmdstanr")
  if (requested == "cmdstanr") {
    stop("CmdStanR est demande mais CmdStan n'est pas disponible.")
  }
  if (requireNamespace("rstan", quietly = TRUE)) return("rstan")
  stop(
    "Aucun moteur Stan utilisable. Installer CmdStanR ou RStan.",
    call. = FALSE
  )
}

require_stan_backend <- function() {
  backend <- stan_backend()
  require_packages("posterior")
  message("Moteur Stan : ", backend)
  invisible(backend)
}

log_sum_exp <- function(x) {
  maximum <- max(x)
  if (!is.finite(maximum)) return(maximum)
  maximum + log(sum(exp(x - maximum)))
}

log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

row_log_sum_exp <- function(x) {
  x <- as.matrix(x)
  maxima <- apply(x, 1L, max)
  result <- maxima
  finite_rows <- is.finite(maxima)
  if (any(finite_rows)) {
    shifted <- sweep(
      x[finite_rows, , drop = FALSE],
      1L,
      maxima[finite_rows],
      `-`
    )
    result[finite_rows] <- maxima[finite_rows] +
      log(rowSums(exp(shifted)))
  }
  result
}

softmax <- function(x) {
  shifted <- x - max(x)
  exp_shifted <- exp(shifted)
  exp_shifted / sum(exp_shifted)
}

stable_seed <- function(...) {
  key <- paste(..., collapse = "|")
  ints <- utf8ToInt(key)
  value <- sum((seq_along(ints) + 17L) * ints)
  as.integer((value %% (.Machine$integer.max - 1L)) + 1L)
}

write_csv_atomic <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  utils::write.csv(x, temporary, row.names = row.names, na = "")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) {
    stop("Impossible d'ecrire le fichier : ", path)
  }
  invisible(path)
}

save_rds_atomic <- function(x, path, compress = TRUE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  saveRDS(x, temporary, compress = compress)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) {
    stop("Impossible d'ecrire le fichier : ", path)
  }
  invisible(path)
}

force_recompute <- function() {
  identical(Sys.getenv("MEMOIRE_FORCE", "0"), "1")
}

parse_profile_argument <- function(args = commandArgs(trailingOnly = TRUE)) {
  match <- grep("^--profile=", args, value = TRUE)
  if (!length(match)) return(Sys.getenv("MEMOIRE_PROFILE", "full"))
  sub("^--profile=", "", match[[1L]])
}

message_step <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                  paste0(..., collapse = "")))
}

group_mean <- function(data, value_columns, by_columns) {
  stats::aggregate(
    data[value_columns],
    by = data[by_columns],
    FUN = function(x) mean(x, na.rm = TRUE)
  )
}
