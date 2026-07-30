# Configuration centrale du pipeline.
#
# Le profil "full" correspond a l'analyse definitive. Le profil "smoke"
# conserve les ages 50-90 et le debut en 1970, mais reduit le cout de calcul.

minkos_config <- function(profile = Sys.getenv("MEMOIRE_PROFILE", "full"),
                          root = find_project_root()) {
  profile <- tolower(profile)
  if (!profile %in% c("full", "smoke")) {
    stop("Profil inconnu : ", profile, ". Utiliser 'full' ou 'smoke'.")
  }
  local_paths <- list()
  local_paths_file <- file.path(root, "config", "paths_local.R")
  if (file.exists(local_paths_file)) {
    local_environment <- new.env(parent = baseenv())
    sys.source(local_paths_file, envir = local_environment)
    if (exists("paths_local", envir = local_environment, inherits = FALSE)) {
      local_paths <- get(
        "paths_local", envir = local_environment, inherits = FALSE
      )
    }
  }
  deaths_path <- if (!is.null(local_paths$deaths)) {
    local_paths$deaths
  } else {
    file.path(root, "data", "raw", "death.txt")
  }
  exposure_path <- if (!is.null(local_paths$exposure)) {
    local_paths$exposure
  } else {
    file.path(root, "data", "raw", "exposure.txt")
  }

  cfg <- list(
    profile = profile,
    root = root,
    sex = "Total",
    models = c("lc", "rh", "apc", "cbd", "m6"),
    dense_metric_models = c("rh", "apc", "m6"),
    model_labels = c(
      lc = "Lee-Carter",
      rh = "Renshaw-Haberman",
      apc = "APC",
      cbd = "CBD",
      m6 = "M6"
    ),
    age_min = 50L,
    age_max = 90L,
    year_start = 1970L,
    # Les observations 2016-2024 restent strictement hors apprentissage.
    # Elles servent uniquement a l'evaluation depuis l'origine fixe 2015.
    data_end = 2024L,
    validation_start = 2001L,
    validation_end = 2015L,
    test_start = 2016L,
    test_end = 2024L,
    lfo_origins = 2000:2014,
    lfo_horizons = 1:10,
    seed = 26052026L,
    forecast_draws = 1000L,
    compute_bma = TRUE,
    bridge_repetitions = 3L,
    contextual_multistarts = 12L,
    discount_rate = 0.02,
    annuity_age = 65L,
    annuity_duration = 25L,
    mcmc = list(
      chains = 4L,
      parallel_chains = 4L,
      iter_warmup = 1000L,
      iter_sampling = 1000L,
      adapt_delta = 0.95,
      max_treedepth = 12L,
      metric = "diag_e",
      refresh = 100L,
      rhat_max = 1.05,
      ess_min = 400,
      max_retries = 1L
    ),
    hierarchical_mcmc = list(
      chains = 4L,
      parallel_chains = 4L,
      iter_warmup = 1000L,
      iter_sampling = 1000L,
      adapt_delta = 0.98,
      max_treedepth = 12L,
      metric = "dense_e",
      refresh = 100L,
      rhat_max = 1.05,
      ess_min = 400,
      max_retries = 1L
    ),
    simulation = list(
      protocol_version = "stan_independent_repetitions_v3",
      scenarios = c(
        "constant_weights",
        "horizon_only",
        "age_horizon",
        "low_information"
      ),
      repetitions = 30L,
      pilot_repetitions = 5L,
      predictive_draws = 1000L,
      test_predictive_draws = 500L,
      posterior_weight_draws = 2000L,
      multistarts = 2L,
      ebfmi_min = 0.30,
      rhat_max = 1.05,
      rhat_target = 1.01,
      ess_min = 400L,
      light_mcmc = list(
        chains = 4L,
        parallel_chains = 4L,
        iter_warmup = 500L,
        iter_sampling = 500L,
        adapt_delta = 0.98,
        max_treedepth = 12L,
        metric = "dense_e",
        refresh = 100L,
        rhat_max = 1.05,
        ess_min = 400L,
        max_retries = 0L
      ),
      full_mcmc = list(
        chains = 4L,
        parallel_chains = 4L,
        iter_warmup = 1000L,
        iter_sampling = 1000L,
        adapt_delta = 0.99,
        max_treedepth = 12L,
        metric = "dense_e",
        refresh = 100L,
        rhat_max = 1.05,
        ess_min = 400L,
        max_retries = 0L
      ),
      deep_mcmc = list(
        chains = 4L,
        parallel_chains = 4L,
        iter_warmup = 1000L,
        iter_sampling = 1000L,
        adapt_delta = 0.99,
        max_treedepth = 15L,
        metric = "dense_e",
        refresh = 100L,
        rhat_max = 1.05,
        ess_min = 400L,
        max_retries = 0L
      ),
      validation_thresholds = list(
        mean_weight_absolute_difference = 0.02,
        maximum_weight_absolute_difference = 0.05,
        mean_relative_logs_difference = 0.02,
        mean_relative_crps_difference = 0.02,
        minimum_ranking_agreement = 0.80
      )
    )
  )

  cfg$lfo_mcmc <- cfg$mcmc

  if (profile == "smoke") {
    cfg$validation_start <- 2009L
    cfg$validation_end <- 2012L
    cfg$test_start <- 2013L
    cfg$test_end <- 2015L
    cfg$data_end <- 2015L
    cfg$lfo_origins <- c(2008L, 2010L, 2011L)
    cfg$lfo_horizons <- 1:2
    cfg$forecast_draws <- 500L
    cfg$bridge_repetitions <- 1L
    cfg$contextual_multistarts <- 5L
    cfg$annuity_duration <- 15L
    cfg$mcmc$iter_warmup <- 750L
    cfg$mcmc$iter_sampling <- 500L
    cfg$mcmc$refresh <- 50L
    cfg$mcmc$ess_min <- 100
    cfg$mcmc$max_retries <- 0L
    cfg$lfo_mcmc <- cfg$mcmc
    cfg$hierarchical_mcmc <- cfg$mcmc
    cfg$hierarchical_mcmc$adapt_delta <- 0.98
    cfg$hierarchical_mcmc$metric <- "dense_e"
    # Le profil smoke verifie l'orchestration, pas la convergence finale.
    cfg$simulation$pilot_repetitions <- 1L
    cfg$simulation$repetitions <- 1L
    cfg$simulation$predictive_draws <- 100L
    cfg$simulation$test_predictive_draws <- 50L
    cfg$simulation$posterior_weight_draws <- 100L
    cfg$simulation$ess_min <- 20L
    for (setting in c("light_mcmc", "full_mcmc", "deep_mcmc")) {
      cfg$simulation[[setting]]$chains <- 2L
      cfg$simulation[[setting]]$parallel_chains <- 2L
      cfg$simulation[[setting]]$iter_warmup <- 100L
      cfg$simulation[[setting]]$iter_sampling <- 100L
      cfg$simulation[[setting]]$refresh <- 20L
      cfg$simulation[[setting]]$rhat_max <- 1.20
      cfg$simulation[[setting]]$ess_min <- 20L
    }
  }

  profile_root <- file.path(root, "results", profile)
  cfg$paths <- list(
    deaths = deaths_path,
    exposure = exposure_path,
    processed = file.path(root, "data", "processed", profile),
    results = profile_root,
    fits = file.path(profile_root, "MCMC_draws"),
    diagnostics = file.path(profile_root, "diagnostics"),
    weights = file.path(profile_root, "weights"),
    metrics = file.path(profile_root, "metrics"),
    figures = file.path(profile_root, "figures"),
    simulation = file.path(profile_root, "simulation"),
    stan = file.path(root, "stan")
  )
  cfg
}
