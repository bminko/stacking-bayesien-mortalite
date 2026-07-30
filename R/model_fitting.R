# Compilation, ajustement et diagnostics des modeles Stan.

model_stan_file <- function(model, cfg) {
  mapping <- c(
    lc = "lc_negbin.stan",
    rh = "rh_negbin.stan",
    apc = "apc_negbin.stan",
    cbd = "cbd_negbin.stan",
    m6 = "m6_negbin.stan",
    hierarchical = "hierarchical_stacking.stan"
  )
  assert_true(model %in% names(mapping), paste("Modele inconnu :", model))
  file.path(cfg$paths$stan, unname(mapping[[model]]))
}

compile_stan_model <- function(model, cfg) {
  backend <- require_stan_backend()
  stan_file <- model_stan_file(model, cfg)
  assert_true(file.exists(stan_file), paste("Fichier Stan absent :", stan_file))
  if (backend == "cmdstanr") {
    return(cmdstanr::cmdstan_model(
      stan_file,
      compile = TRUE,
      stanc_options = list("O1")
    ))
  }
  rstan::rstan_options(auto_write = TRUE)
  options(mc.cores = cfg$mcmc$parallel_chains)
  rstan::stan_model(
    file = stan_file,
    model_name = paste0("memoire_", model),
    auto_write = TRUE,
    save_dso = TRUE
  )
}

fit_cache_path <- function(output_dir) {
  file.path(output_dir, "fit.rds")
}

fit_cache_metadata_path <- function(output_dir) {
  file.path(output_dir, "cache_metadata.rds")
}

stan_fit_cache_key <- function(model, cfg, stan_data, sampling_settings) {
  stan_file <- model_stan_file(model, cfg)
  source_md5 <- unname(tools::md5sum(stan_file))
  temporary <- tempfile("memoire_stan_data_", fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(
    list(data = stan_data, sampling_settings = sampling_settings),
    temporary,
    version = 3
  )
  data_md5 <- unname(tools::md5sum(temporary))
  paste(source_md5, data_md5, sep = ":")
}

read_cached_fit <- function(output_dir, expected_key = NULL) {
  path <- fit_cache_path(output_dir)
  if (!file.exists(path) || force_recompute()) return(NULL)
  metadata_path <- fit_cache_metadata_path(output_dir)
  if (!file.exists(metadata_path)) return(NULL)
  metadata <- tryCatch(readRDS(metadata_path), error = function(error) NULL)
  if (is.null(metadata) ||
      (!is.null(expected_key) && !identical(metadata$key, expected_key))) {
    return(NULL)
  }
  fit <- tryCatch(readRDS(path), error = function(error) NULL)
  if (is.null(fit)) return(NULL)
  if (!is_rstan_fit(fit)) {
    files <- tryCatch(fit$output_files(), error = function(error) character())
    if (!length(files) || !all(file.exists(files))) return(NULL)
  }
  list(
    fit = fit,
    sampling_settings = metadata$sampling_settings
  )
}

fit_once <- function(compiled_model, stan_data, output_dir, mcmc,
                     seed, basename) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (inherits(compiled_model, "stanmodel")) {
    return(rstan::sampling(
      object = compiled_model,
      data = stan_data,
      seed = seed,
      chains = mcmc$chains,
      cores = min(mcmc$parallel_chains, mcmc$chains),
      iter = mcmc$iter_warmup + mcmc$iter_sampling,
      warmup = mcmc$iter_warmup,
      control = list(
        adapt_delta = mcmc$adapt_delta,
        max_treedepth = mcmc$max_treedepth,
        metric = mcmc$metric
      ),
      refresh = mcmc$refresh
    ))
  }
  compiled_model$sample(
    data = stan_data,
    seed = seed,
    chains = mcmc$chains,
    parallel_chains = min(mcmc$parallel_chains, mcmc$chains),
    iter_warmup = mcmc$iter_warmup,
    iter_sampling = mcmc$iter_sampling,
    adapt_delta = mcmc$adapt_delta,
    max_treedepth = mcmc$max_treedepth,
    metric = mcmc$metric,
    refresh = mcmc$refresh,
    output_dir = output_dir,
    output_basename = basename,
    save_warmup = FALSE
  )
}

diagnose_fit <- function(fit, mcmc, label = "") {
  if (is_rstan_fit(fit)) {
    require_packages(c("rstan", "posterior"))
    raw_draws <- rstan::extract(
      fit, permuted = FALSE, inc_warmup = FALSE
    )
    draws <- posterior::as_draws_array(raw_draws)
    summary <- as.data.frame(posterior::summarise_draws(draws))
    sampler_list <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    sampler <- do.call(rbind, lapply(seq_along(sampler_list), function(chain) {
      values <- sampler_list[[chain]]
      energy <- values[, "energy__"]
      data.frame(
        chain_id = chain,
        num_divergent = sum(values[, "divergent__"]),
        num_max_treedepth = sum(
          values[, "treedepth__"] >= mcmc$max_treedepth
        ),
        ebfmi = mean(diff(energy)^2) / stats::var(energy)
      )
    }))
  } else {
    summary <- fit$summary()
    sampler <- fit$diagnostic_summary()
  }
  # lp__ est une quantite auxiliaire de Stan, pas un parametre du modele.
  # Les constantes transformees sont deja retirees par les tests is.finite.
  keep <- summary$variable != "lp__" &
    is.finite(summary$rhat) &
    is.finite(summary$ess_bulk) &
    is.finite(summary$ess_tail)
  checked <- summary[keep, ]

  divergent_column <- intersect(
    c("num_divergent", "num_divergent__"), names(sampler)
  )
  depth_column <- intersect(
    c("num_max_treedepth", "num_max_treedepth__"), names(sampler)
  )
  divergences <- if (length(divergent_column)) {
    sum(sampler[[divergent_column[[1L]]]])
  } else {
    NA_integer_
  }
  max_depth <- if (length(depth_column)) {
    sum(sampler[[depth_column[[1L]]]])
  } else {
    NA_integer_
  }

  max_rhat <- if (nrow(checked)) max(checked$rhat) else Inf
  min_ess_bulk <- if (nrow(checked)) min(checked$ess_bulk) else 0
  min_ess_tail <- if (nrow(checked)) min(checked$ess_tail) else 0
  pass <- max_rhat < mcmc$rhat_max &&
    min_ess_bulk > mcmc$ess_min &&
    min_ess_tail > mcmc$ess_min &&
    isTRUE(as.numeric(divergences) == 0) &&
    isTRUE(as.numeric(max_depth) == 0)

  overview <- data.frame(
    label = label,
    max_rhat = max_rhat,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail,
    divergences = divergences,
    max_treedepth_hits = max_depth,
    rhat_limit = mcmc$rhat_max,
    ess_limit = mcmc$ess_min,
    pass = pass,
    stringsAsFactors = FALSE
  )
  list(pass = pass, overview = overview, parameters = summary,
       sampler = sampler)
}

retry_mcmc_settings <- function(mcmc) {
  retry <- mcmc
  retry$iter_warmup <- as.integer(max(1500L, 2L * mcmc$iter_warmup))
  retry$iter_sampling <- as.integer(max(1500L, 2L * mcmc$iter_sampling))
  retry$adapt_delta <- max(0.99, mcmc$adapt_delta)
  retry$max_treedepth <- max(14L, mcmc$max_treedepth + 2L)
  retry
}

fit_with_diagnostics <- function(compiled_model, stan_data, output_dir,
                                 mcmc, seed, label,
                                 cache_key = NULL,
                                 strict = TRUE) {
  cached <- read_cached_fit(output_dir, cache_key)
  if (!is.null(cached)) {
    diagnostic_settings <- cached$sampling_settings %||% mcmc
    diagnostics <- diagnose_fit(cached$fit, diagnostic_settings, label)
    return(list(
      fit = cached$fit,
      diagnostics = diagnostics,
      cached = TRUE,
      sampling_settings = diagnostic_settings
    ))
  }

  attempts <- seq_len(mcmc$max_retries + 1L)
  current_mcmc <- mcmc
  final <- NULL
  for (attempt in attempts) {
    message_step(label, " - ajustement Stan, essai ", attempt)
    fit <- fit_once(
      compiled_model, stan_data, output_dir, current_mcmc,
      seed + attempt - 1L, paste0("draws_attempt_", attempt)
    )
    diagnostics <- diagnose_fit(fit, current_mcmc, label)
    final <- list(
      fit = fit,
      diagnostics = diagnostics,
      cached = FALSE,
      sampling_settings = current_mcmc
    )
    if (diagnostics$pass) break
    if (attempt < max(attempts)) current_mcmc <- retry_mcmc_settings(mcmc)
  }

  save_rds_atomic(final$fit, fit_cache_path(output_dir), compress = FALSE)
  save_rds_atomic(
    list(key = cache_key, sampling_settings = final$sampling_settings),
    fit_cache_metadata_path(output_dir)
  )
  write_csv_atomic(
    final$diagnostics$overview,
    file.path(output_dir, "diagnostics_overview.csv")
  )
  write_csv_atomic(
    final$diagnostics$parameters,
    file.path(output_dir, "diagnostics_parameters.csv")
  )
  write_csv_atomic(
    final$diagnostics$sampler,
    file.path(output_dir, "diagnostics_sampler.csv")
  )

  if (strict && !final$diagnostics$pass) {
    stop(
      "L'ajustement ", label,
      " ne satisfait pas les seuils de convergence apres reprise. ",
      "Consulter ", file.path(output_dir, "diagnostics_overview.csv"),
      call. = FALSE
    )
  }
  final
}

fit_mortality_model <- function(model, training, compiled_model, cfg,
                                context, strict = cfg$profile == "full",
                                mcmc = cfg$mcmc) {
  output_dir <- file.path(cfg$paths$fits, context, model)
  stan_data <- stan_data_for_model(training, model)
  if (model %in% cfg$dense_metric_models) mcmc$metric <- "dense_e"
  cache_key <- stan_fit_cache_key(model, cfg, stan_data, mcmc)
  fit_with_diagnostics(
    compiled_model = compiled_model,
    stan_data = stan_data,
    output_dir = output_dir,
    mcmc = mcmc,
    seed = stable_seed(cfg$seed, cfg$profile, context, model),
    label = paste(context, model, sep = "/"),
    cache_key = cache_key,
    strict = strict
  )
}

read_model_fit <- function(model, training, cfg, context) {
  stan_data <- stan_data_for_model(training, model)
  mcmc <- cfg$mcmc
  if (model %in% cfg$dense_metric_models) mcmc$metric <- "dense_e"
  cache_key <- stan_fit_cache_key(model, cfg, stan_data, mcmc)
  output_dir <- file.path(cfg$paths$fits, context, model)
  cached <- read_cached_fit(output_dir, cache_key)
  assert_true(!is.null(cached), paste(
    "Ajustement absent ou perime pour", model,
    "- relancer scripts/02_fit_individual_models.R."
  ))
  cached$fit
}

collect_diagnostic_overviews <- function(cfg) {
  files <- list.files(
    cfg$paths$fits,
    pattern = "diagnostics_overview[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(files)) return(data.frame())
  rows <- lapply(files, utils::read.csv, stringsAsFactors = FALSE)
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
