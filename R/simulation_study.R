# Etude de simulation du stacking contextuel hierarchique.
#
# Une repetition correspond ici a un jeu de validation independant, un jeu
# de test independant et un ajustement Stan distinct. Les resultats et les
# diagnostics sont sauvegardes repetition par repetition afin de permettre
# une reprise sans recommencer les ajustements deja valides.

simulation_output_root <- function(cfg) {
  subdir <- cfg$simulation$output_subdir %||% ""
  if (!nzchar(subdir)) return(cfg$paths$simulation)
  file.path(cfg$paths$simulation, subdir)
}

simulation_grid <- function(cfg) {
  grid <- expand.grid(
    age = cfg$age_min:cfg$age_max,
    horizon = cfg$lfo_horizons
  )
  grid[order(grid$horizon, grid$age), , drop = FALSE]
}

simulation_reference_dir <- function(cfg) {
  file.path(simulation_output_root(cfg), "reference_hmd_1970_2015")
}

simulation_reference_path <- function(cfg) {
  file.path(
    simulation_reference_dir(cfg),
    "reference_predictive_distributions.rds"
  )
}

simulation_reference_cache_key <- function(cfg) {
  context <- cfg$simulation$reference_fit_context
  fit_metadata <- file.path(
    cfg$paths$fits, context, cfg$models, "cache_metadata.rds"
  )
  input_paths <- c(
    cfg$paths$deaths,
    cfg$paths$exposure,
    file.path(cfg$root, "R", "data.R"),
    file.path(cfg$root, "R", "forecasting.R"),
    file.path(cfg$root, "R", "simulation_study.R"),
    fit_metadata
  )
  missing <- input_paths[!file.exists(input_paths)]
  assert_true(
    !length(missing),
    paste(
      "Calibration HMD impossible, fichiers absents :",
      paste(missing, collapse = ", ")
    )
  )
  paste(
    cfg$simulation$protocol_version,
    cfg$sex,
    cfg$year_start,
    cfg$simulation$calibration_end,
    cfg$age_min,
    cfg$age_max,
    paste(cfg$lfo_horizons, collapse = ","),
    cfg$simulation$exposure_reference_year,
    cfg$simulation$reference_draws,
    paste(unname(tools::md5sum(input_paths)), collapse = ":"),
    sep = "|"
  )
}

simulation_reference_target <- function(processed, cfg) {
  grid <- simulation_grid(cfg)
  reference_year <- cfg$simulation$exposure_reference_year
  exposure_profile <- processed$long[
    processed$long$year == reference_year,
    c("age", "exposure"),
    drop = FALSE
  ]
  exposure_index <- match(grid$age, exposure_profile$age)
  assert_true(
    !anyNA(exposure_index),
    paste("Profil d'exposition HMD incomplet pour", reference_year)
  )
  exposure <- exposure_profile$exposure[exposure_index]
  assert_true(
    all(is.finite(exposure) & exposure > 0),
    "Les expositions HMD de reference doivent etre positives et finies."
  )
  data.frame(
    year = as.integer(cfg$simulation$calibration_end + grid$horizon),
    age = as.integer(grid$age),
    deaths = NA_real_,
    exposure = as.numeric(exposure),
    origin = as.integer(cfg$simulation$calibration_end),
    horizon = as.integer(grid$horizon)
  )
}

read_simulation_reference_fit <- function(model, cfg, context) {
  fit_dir <- file.path(cfg$paths$fits, context, model)
  fit_path <- file.path(fit_dir, "fit.rds")
  metadata_path <- file.path(fit_dir, "cache_metadata.rds")
  diagnostics_path <- file.path(fit_dir, "diagnostics_overview.csv")
  required <- c(fit_path, metadata_path, diagnostics_path)
  assert_true(
    all(file.exists(required)),
    paste("Ajustement HMD 1970-2015 incomplet pour", model)
  )
  metadata <- readRDS(metadata_path)
  stored_stan_md5 <- strsplit(metadata$key, ":", fixed = TRUE)[[1L]][1L]
  current_stan_md5 <- unname(tools::md5sum(model_stan_file(model, cfg)))
  assert_true(
    identical(stored_stan_md5, current_stan_md5),
    paste(
      "Le code Stan a change depuis l'ajustement HMD de", model,
      ": une reestimation est necessaire."
    )
  )
  diagnostics <- utils::read.csv(
    diagnostics_path, stringsAsFactors = FALSE
  )
  assert_true(
    nrow(diagnostics) == 1L && isTRUE(diagnostics$pass[[1L]]) &&
      diagnostics$divergences[[1L]] == 0 &&
      diagnostics$max_treedepth_hits[[1L]] == 0,
    paste("Ajustement HMD non convergent pour", model)
  )
  fit <- readRDS(fit_path)
  assert_true(is_rstan_fit(fit), paste("Objet Stan invalide pour", model))
  fit
}

build_simulation_reference <- function(cfg) {
  processed_path <- file.path(
    cfg$paths$processed, "mortality_data.rds"
  )
  assert_true(
    file.exists(processed_path),
    paste("Donnees HMD preparees absentes :", processed_path)
  )
  processed <- readRDS(processed_path)
  freshly_prepared <- prepare_hmd_data(cfg)
  assert_true(
    identical(processed$long, freshly_prepared$long),
    paste(
      "Les donnees preparees ne correspondent plus aux fichiers HMD :",
      "relancer la preparation avant la calibration."
    )
  )
  calibration_end <- cfg$simulation$calibration_end
  assert_true(
    calibration_end == cfg$validation_end,
    "La calibration de simulation doit rester alignee sur la fin 2015."
  )
  training <- subset_training_data(processed, calibration_end)
  target <- simulation_reference_target(processed, cfg)
  context <- cfg$simulation$reference_fit_context
  reference_cfg <- cfg
  reference_cfg$forecast_draws <- cfg$simulation$reference_draws

  force_draws <- setNames(vector("list", length(cfg$models)), cfg$models)
  phi_draws <- setNames(vector("list", length(cfg$models)), cfg$models)
  rate_summaries <- list()
  phi_summaries <- list()
  for (model in cfg$models) {
    message_step(
      "Reference HMD 1970-", calibration_end, " - ", toupper(model)
    )
    fit <- read_simulation_reference_fit(model, cfg, context)
    model_seed <- stable_seed(
      cfg$seed, cfg$profile, "simulation_reference_hmd", model
    )
    prediction <- forecast_model(
      fit, model, training, target, reference_cfg,
      keep_draws = TRUE,
      seed = model_seed
    )
    posterior <- fit_draw_matrix(
      fit,
      ndraws = reference_cfg$forecast_draws,
      seed = stable_seed(model_seed, model, "posterior_subset")
    )
    force_draws[[model]] <- prediction$force_draws
    phi_draws[[model]] <- extract_scalar_parameter(posterior, "phi")
    assert_true(
      identical(dim(force_draws[[model]]), c(
        reference_cfg$forecast_draws, nrow(target)
      )),
      paste("Dimensions predictives inattendues pour", model)
    )
    assert_true(
      all(is.finite(force_draws[[model]]) & force_draws[[model]] > 0) &&
        all(is.finite(phi_draws[[model]]) & phi_draws[[model]] > 0),
      paste("Tirages predictifs invalides pour", model)
    )
    rate_summaries[[model]] <- data.frame(
      model = model,
      age = target$age,
      horizon = target$horizon,
      mean_force = colMeans(force_draws[[model]]),
      q025_force = apply(
        force_draws[[model]], 2L, stats::quantile,
        probs = 0.025, names = FALSE
      ),
      q50_force = apply(
        force_draws[[model]], 2L, stats::quantile,
        probs = 0.50, names = FALSE
      ),
      q975_force = apply(
        force_draws[[model]], 2L, stats::quantile,
        probs = 0.975, names = FALSE
      ),
      stringsAsFactors = FALSE
    )
    phi_summaries[[model]] <- data.frame(
      model = model,
      mean_phi = mean(phi_draws[[model]]),
      q025_phi = stats::quantile(
        phi_draws[[model]], 0.025, names = FALSE
      ),
      q50_phi = stats::quantile(
        phi_draws[[model]], 0.50, names = FALSE
      ),
      q975_phi = stats::quantile(
        phi_draws[[model]], 0.975, names = FALSE
      ),
      stringsAsFactors = FALSE
    )
  }

  metadata <- data.frame(
    protocol_version = cfg$simulation$protocol_version,
    data_source = "Human Mortality Database",
    country = processed$metadata$country,
    sex = processed$metadata$sex,
    calibration_years = paste0(cfg$year_start, "-", calibration_end),
    ages = paste0(cfg$age_min, "-", cfg$age_max),
    horizons = paste(range(cfg$lfo_horizons), collapse = "-"),
    candidate_models = paste(cfg$models, collapse = ","),
    fit_context = context,
    posterior_predictive_draws = reference_cfg$forecast_draws,
    observation_distribution = "negative_binomial_2",
    exposure_construction = paste0(
      "profil HMD observe en ",
      cfg$simulation$exposure_reference_year,
      " repete aux horizons 1-10"
    ),
    low_information_exposure_scale =
      cfg$simulation$low_information_exposure_scale,
    low_information_component_contraction =
      cfg$simulation$low_information_component_contraction,
    test_years_used_for_calibration = "aucune",
    stringsAsFactors = FALSE
  )
  reference <- list(
    protocol_version = cfg$simulation$protocol_version,
    cache_key = simulation_reference_cache_key(cfg),
    metadata = metadata,
    target = target,
    force_draws = force_draws,
    phi_draws = phi_draws
  )
  reference_dir <- simulation_reference_dir(cfg)
  dir.create(reference_dir, recursive = TRUE, showWarnings = FALSE)
  save_rds_atomic(
    reference, simulation_reference_path(cfg), compress = FALSE
  )
  write_csv_atomic(
    metadata, file.path(reference_dir, "dgp_metadata.csv")
  )
  write_csv_atomic(
    do.call(rbind, rate_summaries),
    file.path(reference_dir, "candidate_force_summaries.csv")
  )
  write_csv_atomic(
    do.call(rbind, phi_summaries),
    file.path(reference_dir, "candidate_dispersion_summaries.csv")
  )
  reference
}

load_or_build_simulation_reference <- function(cfg) {
  path <- simulation_reference_path(cfg)
  expected_key <- simulation_reference_cache_key(cfg)
  if (file.exists(path) && !force_recompute()) {
    cached <- tryCatch(readRDS(path), error = function(error) NULL)
    if (!is.null(cached) &&
        identical(cached$protocol_version, cfg$simulation$protocol_version) &&
        identical(cached$cache_key, expected_key)) {
      message_step("Reference predictive HMD - resultat en cache")
      return(cached)
    }
  }
  build_simulation_reference(cfg)
}

simulation_predictive_components <- function(reference, scenario, cfg) {
  force_draws <- reference$force_draws[cfg$models]
  phi_draws <- reference$phi_draws[cfg$models]
  if (scenario == "low_information") {
    contraction <- cfg$simulation$low_information_component_contraction
    assert_true(
      is.finite(contraction) && contraction > 0 && contraction < 1,
      "La contraction low_information doit appartenir a ]0,1[."
    )
    log_forces <- lapply(force_draws, log)
    log_reference <- Reduce(`+`, log_forces) / length(log_forces)
    force_draws <- lapply(
      log_forces,
      function(log_force) {
        exp(log_reference + contraction * (log_force - log_reference))
      }
    )
  }
  list(force_draws = force_draws, phi_draws = phi_draws)
}

simulation_log_predictive_counts <- function(
    deaths, exposure, force_draws, phi_draws) {
  force_draws <- as.matrix(force_draws)
  S <- nrow(force_draws)
  N <- ncol(force_draws)
  assert_true(length(deaths) == N, "Nombre de deces incompatible.")
  assert_true(length(exposure) == N, "Nombre d'expositions incompatible.")
  assert_true(length(phi_draws) == S, "Nombre de dispersions incompatible.")
  means <- sweep(force_draws, 2L, exposure, `*`)
  log_likelihood <- matrix(
    stats::dnbinom(
      rep(deaths, each = S),
      mu = as.vector(means),
      size = rep(phi_draws, times = N),
      log = TRUE
    ),
    nrow = S,
    ncol = N
  )
  apply(log_likelihood, 2L, log_mean_exp)
}

simulation_true_weights <- function(scenario, standardized, models) {
  N <- nrow(standardized)
  if (scenario == "constant_weights") {
    constant <- c(
      lc = 0.06, rh = 0.54, apc = 0.24, cbd = 0.10, m6 = 0.06
    )
    return(matrix(
      constant[models], N, length(models), byrow = TRUE,
      dimnames = list(NULL, models)
    ))
  }
  if (scenario == "low_information") {
    constant <- c(
      lc = 0.18, rh = 0.24, apc = 0.22, cbd = 0.18, m6 = 0.18
    )
    return(matrix(
      constant[models], N, length(models), byrow = TRUE,
      dimnames = list(NULL, models)
    ))
  }
  if (scenario == "horizon_only") {
    scores <- cbind(
      lc = -0.8 - 0.7 * standardized$h_tilde,
      rh = 0.4 + 1.1 * standardized$h_tilde,
      apc = 0.3 - 0.4 * standardized$h_tilde,
      cbd = -0.2 + 0.25 * standardized$h_tilde,
      m6 = 0
    )
    return(scores_to_weights(scores)[, models, drop = FALSE])
  }
  scores <- cbind(
    lc = -0.7 - 0.9 * standardized$x_tilde -
      0.35 * standardized$h_tilde,
    rh = 0.3 + 0.8 * standardized$h_tilde +
      0.55 * standardized$x_h,
    apc = 0.5 - 0.55 * standardized$q_age -
      0.25 * standardized$h_tilde,
    cbd = -0.2 + 0.65 * standardized$x_tilde -
      0.45 * standardized$x_h,
    m6 = 0
  )
  scores_to_weights(scores)[, models, drop = FALSE]
}

simulate_study_dataset <- function(
    scenario, seed, cfg, reference, constants = NULL,
    components = NULL) {
  grid <- simulation_grid(cfg)
  standardized <- standardize_context(grid, constants)
  data <- standardized$data
  constants <- standardized$constants
  true_weights <- simulation_true_weights(scenario, data, cfg$models)
  components <- components %||%
    simulation_predictive_components(reference, scenario, cfg)
  reference_key <- paste(reference$target$age, reference$target$horizon)
  data_key <- paste(data$age, data$horizon)
  reference_index <- match(data_key, reference_key)
  assert_true(
    !anyNA(reference_index),
    "La grille de simulation differe de la reference predictive HMD."
  )
  exposure_scale <- if (scenario == "low_information") {
    cfg$simulation$low_information_exposure_scale
  } else {
    1
  }
  data$exposure <-
    reference$target$exposure[reference_index] * exposure_scale
  data$origin <- as.integer(cfg$simulation$calibration_end)
  data$year <- as.integer(data$origin + data$horizon)
  data$omega <- 1

  set.seed(seed)
  surface_rows <- setNames(vapply(
    cfg$models,
    function(model) sample.int(
      nrow(components$force_draws[[model]]), 1L
    ),
    integer(1)
  ), cfg$models)
  generating_forces <- vapply(
    cfg$models,
    function(model) {
      components$force_draws[[model]][
        surface_rows[[model]], reference_index
      ]
    },
    numeric(nrow(data))
  )
  colnames(generating_forces) <- cfg$models
  uniforms <- stats::runif(nrow(data))
  cumulative <- t(apply(true_weights, 1L, cumsum))
  selected <- 1L + rowSums(
    matrix(
      uniforms, nrow(data), length(cfg$models) - 1L
    ) > cumulative[, -length(cfg$models), drop = FALSE]
  )
  selected_models <- cfg$models[selected]
  selected_force <- generating_forces[
    cbind(seq_len(nrow(data)), selected)
  ]
  selected_phi <- vapply(
    seq_len(nrow(data)),
    function(index) {
      model <- selected_models[[index]]
      components$phi_draws[[model]][surface_rows[[model]]]
    },
    numeric(1)
  )
  mean_deaths <- data$exposure * selected_force
  data$deaths <- stats::rnbinom(
    nrow(data), mu = mean_deaths, size = selected_phi
  )
  for (k in seq_along(cfg$models)) {
    model <- cfg$models[[k]]
    data[[paste0("log_p_", model)]] <-
      simulation_log_predictive_counts(
      data$deaths,
      data$exposure,
      components$force_draws[[model]][, reference_index, drop = FALSE],
      components$phi_draws[[model]]
    )
  }
  list(
    data = data,
    constants = constants,
    true_weights = true_weights,
    predictive_force_draws = components$force_draws,
    predictive_phi_draws = components$phi_draws,
    generating_log_rates = log(generating_forces),
    generating_phi = selected_phi,
    selected_models = selected_models,
    surface_draw_rows = surface_rows,
    exposure_scale = exposure_scale,
    reference_cache_key = reference$cache_key
  )
}

simulation_repetition_seed <- function(
    cfg, scenario, repetition, part, cohort = "main") {
  assert_true(
    cohort %in% c("pilot", "main", "validation"),
    paste("Cohorte de simulation inconnue :", cohort)
  )
  stable_seed(
    cfg$seed,
    cfg$simulation$protocol_version,
    cohort,
    scenario,
    sprintf("rep_%03d", repetition),
    part
  )
}

simulation_repetition_dir <- function(
    cfg, scenario, repetition, cohort = "main") {
  assert_true(
    cohort %in% c("pilot", "main"),
    paste("Cohorte de repetition inconnue :", cohort)
  )
  repetition_root <- if (cohort == "pilot") {
    "pilot_repetitions"
  } else {
    "repetitions"
  }
  file.path(
    simulation_output_root(cfg),
    repetition_root,
    scenario,
    sprintf("rep_%03d", repetition)
  )
}

simulation_fit_context <- function(
    scenario, repetition, configuration, cohort = "main") {
  file.path(
    "simulation_v4_hmd",
    cohort,
    scenario,
    sprintf("rep_%03d", repetition),
    configuration
  )
}

simulation_weight_draw_array <- function(
    fit, grid, models, constants, ndraws = NULL, seed = 1L,
    design_columns = c("age", "age2", "horizon", "age_horizon")) {
  standardized <- standardize_context(grid, constants)$data
  X <- context_design_matrix(standardized, design_columns)
  K <- length(models)
  P <- ncol(X)
  draws <- fit_draw_matrix(
    fit,
    variables = c("alpha", "beta"),
    ndraws = ndraws,
    seed = seed
  )
  alpha <- extract_vector_parameter(draws, "alpha", K - 1L)
  beta <- array(NA_real_, dim = c(nrow(draws), K - 1L, P))
  for (k in seq_len(K - 1L)) {
    for (p in seq_len(P)) {
      beta[, k, p] <- extract_scalar_parameter(
        draws, sprintf("beta[%d,%d]", k, p)
      )
    }
  }
  result <- array(
    NA_real_,
    dim = c(nrow(draws), nrow(grid), K),
    dimnames = list(NULL, NULL, models)
  )
  for (cell in seq_len(nrow(grid))) {
    scores <- alpha
    for (p in seq_len(P)) {
      scores <- scores + beta[, , p] * X[cell, p]
    }
    result[, cell, ] <- scores_to_weights(cbind(scores, 0))
  }
  result
}

simulation_weight_summary <- function(
    weight_draws, grid, models, true_weights) {
  means <- apply(weight_draws, c(2L, 3L), mean)
  q025 <- apply(
    weight_draws, c(2L, 3L), stats::quantile, probs = 0.025
  )
  q50 <- apply(
    weight_draws, c(2L, 3L), stats::quantile, probs = 0.50
  )
  q975 <- apply(
    weight_draws, c(2L, 3L), stats::quantile, probs = 0.975
  )
  data.frame(
    age = rep(grid$age, times = length(models)),
    horizon = rep(grid$horizon, times = length(models)),
    model = rep(models, each = nrow(grid)),
    true_weight = as.vector(true_weights),
    estimated_weight = as.vector(means),
    weight_q025 = as.vector(q025),
    weight_q50 = as.vector(q50),
    weight_q975 = as.vector(q975),
    interval_width = as.vector(q975 - q025),
    stringsAsFactors = FALSE
  )
}

simulation_representative_grid <- function(grid) {
  ages <- unique(round(stats::quantile(
    sort(unique(grid$age)), c(0, 0.25, 0.5, 0.75, 1),
    names = FALSE
  )))
  horizons <- unique(round(stats::quantile(
    sort(unique(grid$horizon)), c(0, 0.5, 1),
    names = FALSE
  )))
  expand.grid(age = ages, horizon = horizons)
}

simulation_raw_parameter_draws <- function(fit) {
  if (is_rstan_fit(fit)) {
    return(rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE))
  }
  posterior::as_draws_array(
    fit$draws(variables = c("alpha", "beta", "tau"))
  )
}

simulation_weight_chain_diagnostics <- function(
    fit, grid, models, constants,
    design_columns = c("age", "age2", "horizon", "age_horizon")) {
  require_packages("posterior")
  raw <- simulation_raw_parameter_draws(fit)
  variable_names <- dimnames(raw)[[3L]]
  K <- length(models)
  standardized <- standardize_context(grid, constants)$data
  X <- context_design_matrix(standardized, design_columns)
  P <- ncol(X)
  iterations <- dim(raw)[1L]
  chains <- dim(raw)[2L]
  sample_count <- iterations * chains

  alpha <- matrix(NA_real_, sample_count, K - 1L)
  beta <- array(NA_real_, dim = c(sample_count, K - 1L, P))
  for (k in seq_len(K - 1L)) {
    alpha[, k] <- as.vector(
      raw[, , match(sprintf("alpha[%d]", k), variable_names)]
    )
    for (p in seq_len(P)) {
      beta[, k, p] <- as.vector(
        raw[, , match(sprintf("beta[%d,%d]", k, p), variable_names)]
      )
    }
  }

  variable_count <- nrow(grid) * K
  flat <- matrix(NA_real_, sample_count, variable_count)
  output_names <- character(variable_count)
  position <- 1L
  for (cell in seq_len(nrow(grid))) {
    scores <- alpha
    for (p in seq_len(P)) {
      scores <- scores + beta[, , p] * X[cell, p]
    }
    weights <- scores_to_weights(cbind(scores, 0))
    columns <- position:(position + K - 1L)
    flat[, columns] <- weights
    output_names[columns] <- sprintf(
      "weight_age%d_h%d_%s",
      grid$age[cell],
      grid$horizon[cell],
      models
    )
    position <- position + K
  }
  weight_array <- array(
    flat,
    dim = c(iterations, chains, variable_count),
    dimnames = list(
      iteration = seq_len(iterations),
      chain = seq_len(chains),
      variable = output_names
    )
  )
  summary <- as.data.frame(posterior::summarise_draws(
    posterior::as_draws_array(weight_array)
  ))

  core <- grepl("^(alpha|beta|tau)\\[", variable_names)
  chain_stuck <- any(vapply(seq_len(chains), function(chain) {
    values <- raw[, chain, core, drop = FALSE]
    standard_deviations <- apply(values, 3L, stats::sd)
    any(!is.finite(values)) ||
      !length(standard_deviations) ||
      all(!is.finite(standard_deviations) |
            standard_deviations < 1e-10)
  }, logical(1)))
  list(summary = summary, chain_stuck = chain_stuck)
}

simulation_fit_diagnostics <- function(
    hierarchical, grid, models, constants, settings, cfg) {
  parameters <- hierarchical$diagnostics$parameters
  relevant <- grepl("^(alpha|beta|tau)\\[", parameters$variable)
  relevant_parameters <- parameters[relevant, , drop = FALSE]
  representative <- simulation_representative_grid(grid)
  weight_diagnostics <- simulation_weight_chain_diagnostics(
    hierarchical$fit,
    representative,
    models,
    constants
  )
  weights <- weight_diagnostics$summary
  sampler <- hierarchical$diagnostics$sampler
  minimum_ebfmi <- if ("ebfmi" %in% names(sampler)) {
    min(sampler$ebfmi, na.rm = TRUE)
  } else {
    NA_real_
  }
  overview <- hierarchical$diagnostics$overview[1L, , drop = FALSE]
  max_rhat <- max(
    relevant_parameters$rhat,
    weights$rhat,
    na.rm = TRUE
  )
  min_ess_bulk <- min(
    relevant_parameters$ess_bulk,
    weights$ess_bulk,
    na.rm = TRUE
  )
  min_ess_tail <- min(
    relevant_parameters$ess_tail,
    weights$ess_tail,
    na.rm = TRUE
  )
  valid <- is.finite(max_rhat) &&
    max_rhat <= cfg$simulation$rhat_max &&
    is.finite(min_ess_bulk) &&
    min_ess_bulk >= cfg$simulation$ess_min &&
    is.finite(min_ess_tail) &&
    min_ess_tail >= cfg$simulation$ess_min &&
    isTRUE(overview$divergences == 0) &&
    isTRUE(overview$max_treedepth_hits == 0) &&
    is.finite(minimum_ebfmi) &&
    minimum_ebfmi >= cfg$simulation$ebfmi_min &&
    !weight_diagnostics$chain_stuck

  list(
    overview = data.frame(
      max_rhat = max_rhat,
      min_ess_bulk = min_ess_bulk,
      min_ess_tail = min_ess_tail,
      parameter_max_rhat = max(relevant_parameters$rhat, na.rm = TRUE),
      parameter_min_ess_bulk =
        min(relevant_parameters$ess_bulk, na.rm = TRUE),
      parameter_min_ess_tail =
        min(relevant_parameters$ess_tail, na.rm = TRUE),
      weight_max_rhat = max(weights$rhat, na.rm = TRUE),
      weight_min_ess_bulk = min(weights$ess_bulk, na.rm = TRUE),
      weight_min_ess_tail = min(weights$ess_tail, na.rm = TRUE),
      divergences = overview$divergences,
      max_treedepth_hits = overview$max_treedepth_hits,
      min_ebfmi = minimum_ebfmi,
      chain_stuck = weight_diagnostics$chain_stuck,
      rhat_target_1_01 = max_rhat <= cfg$simulation$rhat_target,
      valid = valid,
      chains = settings$chains,
      iter_warmup = settings$iter_warmup,
      iter_sampling = settings$iter_sampling,
      kept_draws = settings$chains * settings$iter_sampling,
      adapt_delta = settings$adapt_delta,
      max_treedepth = settings$max_treedepth,
      metric = settings$metric,
      stringsAsFactors = FALSE
    ),
    parameters = relevant_parameters,
    weights = weights,
    sampler = sampler
  )
}

simulation_attempt <- function(
    scenario, repetition, configuration, validation,
    compiled_model, settings, cfg, cohort = "main") {
  context <- simulation_fit_context(
    scenario, repetition, configuration, cohort
  )
  attempt_dir <- file.path(
    cfg$paths$fits, context, "hierarchical"
  )
  attempt_diagnostic_path <- file.path(
    attempt_dir, "simulation_diagnostics.csv"
  )
  started <- Sys.time()
  hierarchical <- fit_hierarchical_stacking(
    validation$data,
    cfg$models,
    cfg,
    compiled_model = compiled_model,
    strict = FALSE,
    context = context,
    mcmc = settings
  )
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (isTRUE(hierarchical$cached) &&
      file.exists(attempt_diagnostic_path)) {
    previous <- tryCatch(
      utils::read.csv(
        attempt_diagnostic_path, stringsAsFactors = FALSE
      ),
      error = function(error) NULL
    )
    if (!is.null(previous) && nrow(previous) == 1L &&
        is.finite(previous$elapsed_seconds)) {
      elapsed <- previous$elapsed_seconds
    }
  }
  diagnostics <- simulation_fit_diagnostics(
    hierarchical,
    simulation_grid(cfg),
    cfg$models,
    validation$constants,
    hierarchical$sampling_settings,
    cfg
  )
  overview <- diagnostics$overview
  overview$scenario <- scenario
  overview$repetition <- repetition
  overview$seed <- stable_seed(
    cfg$seed, cfg$profile, "hierarchical", context
  )
  overview$configuration <- configuration
  overview$elapsed_seconds <- elapsed
  overview$cached <- isTRUE(hierarchical$cached)
  overview$status <- if (overview$valid) "valide" else "a_relancer"
  overview <- overview[c(
    "scenario", "repetition", "seed", "configuration",
    "elapsed_seconds", "max_rhat", "min_ess_bulk", "min_ess_tail",
    "parameter_max_rhat", "parameter_min_ess_bulk",
    "parameter_min_ess_tail", "weight_max_rhat",
    "weight_min_ess_bulk", "weight_min_ess_tail",
    "divergences", "max_treedepth_hits", "min_ebfmi",
    "chain_stuck", "rhat_target_1_01", "valid", "status",
    "chains", "iter_warmup", "iter_sampling", "kept_draws",
    "adapt_delta", "max_treedepth", "metric", "cached"
  )]
  write_csv_atomic(overview, attempt_diagnostic_path)
  list(
    hierarchical = hierarchical,
    diagnostics = diagnostics,
    overview = overview
  )
}

simulation_parameter_summary <- function(fit) {
  draws <- fit_draw_matrix(
    fit,
    variables = c("alpha", "beta", "tau")
  )
  keep <- grepl("^(alpha|beta|tau)\\[", colnames(draws))
  draws <- draws[, keep, drop = FALSE]
  rows <- lapply(seq_len(ncol(draws)), function(column) {
    values <- draws[, column]
    data.frame(
      parameter = colnames(draws)[column],
      mean = mean(values),
      sd = stats::sd(values),
      q025 = stats::quantile(values, 0.025, names = FALSE),
      q50 = stats::quantile(values, 0.50, names = FALSE),
      q975 = stats::quantile(values, 0.975, names = FALSE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

simulation_select_models <- function(probabilities, uniforms) {
  cumulative <- t(apply(probabilities, 1L, cumsum))
  1L + rowSums(
    matrix(
      uniforms,
      nrow(probabilities),
      ncol(probabilities) - 1L
    ) > cumulative[, -ncol(probabilities), drop = FALSE]
  )
}

simulation_draw_predictive_counts <- function(
    test, cell, selected, models) {
  counts <- numeric(length(selected))
  for (model_index in unique(selected)) {
    positions <- which(selected == model_index)
    model <- models[[model_index]]
    force <- test$predictive_force_draws[[model]]
    phi <- test$predictive_phi_draws[[model]]
    rows <- sample.int(nrow(force), length(positions), replace = TRUE)
    counts[positions] <- stats::rnbinom(
      length(positions),
      mu = test$data$exposure[cell] * force[rows, cell],
      size = phi[rows]
    )
  }
  counts
}

simulation_evaluate_method <- function(
    test, scenario, repetition, method, models, predictive_draws,
    seed, deterministic_weights = NULL, posterior_weights = NULL) {
  log_p <- meta_logp_matrix(test$data, models)
  cells <- nrow(test$data)
  crps <- numeric(cells)
  logs <- numeric(cells)
  set.seed(seed)
  if (!is.null(posterior_weights)) {
    available <- dim(posterior_weights)[1L]
    selected_draws <- sample.int(
      available,
      predictive_draws,
      replace = available < predictive_draws
    )
    posterior_weights <- posterior_weights[
      selected_draws, , , drop = FALSE
    ]
  }

  for (cell in seq_len(cells)) {
    if (is.null(posterior_weights)) {
      weights <- deterministic_weights[cell, ]
      logs[cell] <- -log_sum_exp(log(weights) + log_p[cell, ])
      selected <- sample.int(
        length(models),
        predictive_draws,
        replace = TRUE,
        prob = weights
      )
    } else {
      weights <- posterior_weights[, cell, , drop = FALSE]
      weights <- matrix(weights, nrow = predictive_draws)
      log_density <- apply(
        log(weights) +
          matrix(log_p[cell, ], predictive_draws, length(models),
                 byrow = TRUE),
        1L,
        log_sum_exp
      )
      logs[cell] <- -log_mean_exp(log_density)
      selected <- simulation_select_models(
        weights, stats::runif(predictive_draws)
      )
    }
    predictive <- simulation_draw_predictive_counts(
      test, cell, selected, models
    ) / test$data$exposure[cell]
    observed_rate <- test$data$deaths[cell] / test$data$exposure[cell]
    crps[cell] <- empirical_crps(predictive, observed_rate)
  }
  data.frame(
    scenario = scenario,
    repetition = repetition,
    method = method,
    mean_logs = mean(logs),
    mean_crps = mean(crps),
    predictive_draws = predictive_draws,
    stringsAsFactors = FALSE
  )
}

simulation_deterministic_weight_long <- function(
    weights, grid, models, true_weights, method) {
  data.frame(
    age = rep(grid$age, times = length(models)),
    horizon = rep(grid$horizon, times = length(models)),
    model = rep(models, each = nrow(grid)),
    method = method,
    true_weight = as.vector(true_weights),
    estimated_weight = as.vector(weights),
    weight_q025 = NA_real_,
    weight_q50 = NA_real_,
    weight_q975 = NA_real_,
    interval_width = NA_real_,
    stringsAsFactors = FALSE
  )
}

simulation_recovery_row <- function(weights_long) {
  difference <- weights_long$estimated_weight - weights_long$true_weight
  hierarchical <- all(is.finite(weights_long$weight_q025))
  data.frame(
    method = unique(weights_long$method),
    weight_bias = mean(difference),
    weight_rmse = sqrt(mean(difference^2)),
    mean_absolute_weight_error = mean(abs(difference)),
    extreme_low_proportion =
      mean(weights_long$estimated_weight < 0.01),
    extreme_high_proportion =
      mean(weights_long$estimated_weight > 0.99),
    extreme_weight_proportion =
      mean(weights_long$estimated_weight < 0.01 |
             weights_long$estimated_weight > 0.99),
    posterior_95_coverage = if (hierarchical) {
      mean(
        weights_long$true_weight >= weights_long$weight_q025 &
          weights_long$true_weight <= weights_long$weight_q975
      )
    } else {
      NA_real_
    },
    mean_interval_width = if (hierarchical) {
      mean(weights_long$interval_width)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

simulation_bind_rows <- function(items, name) {
  rows <- lapply(items, function(item) item[[name]])
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(data.frame())
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

run_simulation_repetition <- function(
    scenario, repetition, compiled_model, cfg,
    predictive_draws = cfg$simulation$predictive_draws,
    cohort = "main", reference = NULL) {
  repetition_dir <- simulation_repetition_dir(
    cfg, scenario, repetition, cohort
  )
  dir.create(repetition_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(repetition_dir, "result.rds")
  if (file.exists(result_path) && !force_recompute()) {
    cached <- tryCatch(readRDS(result_path), error = function(error) NULL)
    if (!is.null(cached) &&
        identical(
          cached$protocol_version,
          cfg$simulation$protocol_version
        ) &&
        identical(
          cached$cohort,
          cohort
        )) {
      message_step(
        "Simulation ", scenario, " repetition ", repetition,
        " - resultat en cache"
      )
      return(cached)
    }
  }

  message_step(
    "Simulation ", scenario, " repetition ", repetition,
    " - generation des donnees"
  )
  reference <- reference %||% load_or_build_simulation_reference(cfg)
  components <- simulation_predictive_components(reference, scenario, cfg)
  validation <- simulate_study_dataset(
    scenario,
    simulation_repetition_seed(
      cfg, scenario, repetition, "validation", cohort
    ),
    cfg,
    reference = reference,
    components = components
  )
  test <- simulate_study_dataset(
    scenario,
    simulation_repetition_seed(
      cfg, scenario, repetition, "test", cohort
    ),
    cfg,
    reference = reference,
    constants = validation$constants,
    components = components
  )
  grid <- simulation_grid(cfg)
  global <- fit_global_stacking(
    validation$data,
    cfg$models,
    cfg$simulation$multistarts,
    simulation_repetition_seed(
      cfg, scenario, repetition, "global", cohort
    )
  )
  contextual <- fit_contextual_stacking(
    validation$data,
    cfg$models,
    cfg$simulation$multistarts,
    simulation_repetition_seed(
      cfg, scenario, repetition, "contextual", cohort
    ),
    allow_nonconvergence = TRUE
  )
  contextual_converged <- identical(contextual$status, "converge")
  contextual_status <- data.frame(
    scenario = scenario,
    repetition = repetition,
    status = contextual$status,
    error_message = contextual$error_message,
    starts = nrow(contextual$optimization_diagnostics),
    converged_starts = sum(
      contextual$optimization_diagnostics$convergence == 0L &
        is.finite(contextual$optimization_diagnostics$value)
    ),
    best_objective = min(
      contextual$optimization_diagnostics$value,
      na.rm = TRUE
    ),
    max_abs_parameter = max(
      contextual$optimization_diagnostics$max_abs_parameter,
      na.rm = TRUE
    ),
    minimum_gradient_norm = min(
      contextual$optimization_diagnostics$gradient_norm,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
  write_csv_atomic(
    contextual_status,
    file.path(repetition_dir, "contextual_optimization_status.csv")
  )
  write_csv_atomic(
    contextual$optimization_diagnostics,
    file.path(repetition_dir, "contextual_optimization_diagnostics.csv")
  )
  if (!contextual_converged) {
    message_step(
      "Simulation ", scenario, " repetition ", repetition,
      " - stacking contextuel non regularise non convergent; ",
      "poursuite du stacking global et hierarchique"
    )
  }

  attempts <- list()
  attempts[[1L]] <- simulation_attempt(
    scenario, repetition, "legere", validation,
    compiled_model, cfg$simulation$light_mcmc, cfg, cohort
  )
  if (!attempts[[1L]]$overview$valid) {
    attempts[[2L]] <- simulation_attempt(
      scenario, repetition, "complete", validation,
      compiled_model, cfg$simulation$full_mcmc, cfg, cohort
    )
  }
  last_attempt <- attempts[[length(attempts)]]
  if (!last_attempt$overview$valid &&
      last_attempt$overview$max_treedepth_hits > 0) {
    attempts[[length(attempts) + 1L]] <- simulation_attempt(
      scenario, repetition, "complete_treedepth15", validation,
      compiled_model, cfg$simulation$deep_mcmc, cfg, cohort
    )
  }

  attempts_overview <- do.call(
    rbind, lapply(attempts, `[[`, "overview")
  )
  rownames(attempts_overview) <- NULL
  write_csv_atomic(
    attempts_overview,
    file.path(repetition_dir, "diagnostics_attempts.csv")
  )
  valid_indices <- which(vapply(
    attempts,
    function(attempt) isTRUE(attempt$overview$valid),
    logical(1)
  ))
  if (!length(valid_indices)) {
    final_overview <- attempts_overview[nrow(attempts_overview), ]
    final_overview$status <- "echec_definitif"
    write_csv_atomic(
      final_overview,
      file.path(repetition_dir, "diagnostics_final.csv")
    )
    result <- list(
      protocol_version = cfg$simulation$protocol_version,
      cohort = cohort,
      scenario = scenario,
      repetition = repetition,
      status = "echec_definitif",
      diagnostics = final_overview,
      diagnostics_attempts = attempts_overview,
      contextual_optimization = contextual_status,
      contextual_optimization_attempts =
        contextual$optimization_diagnostics
    )
    save_rds_atomic(result, result_path, compress = FALSE)
    return(result)
  }

  selected <- attempts[[valid_indices[[1L]]]]
  final_overview <- selected$overview
  final_overview$status <- "valide"
  write_csv_atomic(
    final_overview,
    file.path(repetition_dir, "diagnostics_final.csv")
  )
  # Les diagnostics sont ecrits avant les scores, conformement au protocole.
  write_csv_atomic(
    selected$diagnostics$parameters,
    file.path(repetition_dir, "diagnostics_parameters.csv")
  )
  write_csv_atomic(
    selected$diagnostics$weights,
    file.path(repetition_dir, "diagnostics_weights.csv")
  )
  write_csv_atomic(
    selected$diagnostics$sampler,
    file.path(repetition_dir, "diagnostics_sampler.csv")
  )

  weight_draws <- simulation_weight_draw_array(
    selected$hierarchical$fit,
    grid,
    cfg$models,
    validation$constants,
    ndraws = NULL,
    seed = simulation_repetition_seed(
      cfg, scenario, repetition, "weight_draws", cohort
    )
  )
  true_weights <- simulation_true_weights(
    scenario,
    standardize_context(grid, validation$constants)$data,
    cfg$models
  )
  weight_matrices <- list(
    stacking_global = constant_weight_matrix(
      global$weights, nrow(grid), cfg$models
    )
  )
  if (contextual_converged) {
    contextual_grid <- contextual_weight_grid(
      grid, contextual, cfg$models, validation$constants
    )
    weight_matrices$stacking_contextual <-
      as.matrix(contextual_grid[cfg$models])
  }
  hierarchical_long <- simulation_weight_summary(
    weight_draws, grid, cfg$models, true_weights
  )
  hierarchical_long$method <- "stacking_hierarchical"
  weights_long <- list(
    simulation_deterministic_weight_long(
      weight_matrices$stacking_global,
      grid,
      cfg$models,
      true_weights,
      "stacking_global"
    )
  )
  if (contextual_converged) {
    weights_long[[length(weights_long) + 1L]] <-
      simulation_deterministic_weight_long(
        weight_matrices$stacking_contextual,
        grid,
        cfg$models,
        true_weights,
        "stacking_contextual"
      )
  }
  weights_long[[length(weights_long) + 1L]] <- hierarchical_long
  weights_long <- do.call(rbind, weights_long)
  weights_long$scenario <- scenario
  weights_long$repetition <- repetition
  weights_long <- weights_long[c(
    "scenario", "repetition", "method", "age", "horizon", "model",
    "true_weight", "estimated_weight", "weight_q025", "weight_q50",
    "weight_q975", "interval_width"
  )]

  recovery <- do.call(rbind, lapply(
    split(weights_long, weights_long$method),
    simulation_recovery_row
  ))
  recovery$scenario <- scenario
  recovery$repetition <- repetition
  recovery <- recovery[c(
    "scenario", "repetition", "method", "weight_bias",
    "weight_rmse", "mean_absolute_weight_error",
    "extreme_low_proportion", "extreme_high_proportion",
    "extreme_weight_proportion", "posterior_95_coverage",
    "mean_interval_width"
  )]

  performance <- list()
  for (method in names(weight_matrices)) {
    performance[[method]] <- simulation_evaluate_method(
      test,
      scenario,
      repetition,
      method,
      cfg$models,
      predictive_draws,
      simulation_repetition_seed(
        cfg, scenario, repetition, paste0("predictive_", method),
        cohort
      ),
      deterministic_weights = weight_matrices[[method]]
    )
  }
  performance$stacking_hierarchical <- simulation_evaluate_method(
    test,
    scenario,
    repetition,
    "stacking_hierarchical",
    cfg$models,
    predictive_draws,
    simulation_repetition_seed(
      cfg, scenario, repetition, "predictive_hierarchical", cohort
    ),
    posterior_weights = weight_draws
  )
  performance <- do.call(rbind, performance)

  parameters <- simulation_parameter_summary(
    selected$hierarchical$fit
  )
  parameters$scenario <- scenario
  parameters$repetition <- repetition
  parameters$configuration <- selected$overview$configuration
  parameters <- parameters[c(
    "scenario", "repetition", "configuration", "parameter",
    "mean", "sd", "q025", "q50", "q975"
  )]
  tau <- parameters[grepl("^tau\\[", parameters$parameter), ]

  write_csv_atomic(
    performance, file.path(repetition_dir, "performance.csv")
  )
  write_csv_atomic(
    recovery, file.path(repetition_dir, "weight_recovery.csv")
  )
  write_csv_atomic(
    weights_long, file.path(repetition_dir, "weights.csv")
  )
  write_csv_atomic(
    parameters, file.path(repetition_dir, "parameters.csv")
  )
  write_csv_atomic(tau, file.path(repetition_dir, "tau.csv"))

  result <- list(
    protocol_version = cfg$simulation$protocol_version,
    cohort = cohort,
    scenario = scenario,
    repetition = repetition,
    status = "valide",
    configuration = selected$overview$configuration,
    diagnostics = final_overview,
    diagnostics_attempts = attempts_overview,
    contextual_optimization = contextual_status,
    contextual_optimization_attempts =
      contextual$optimization_diagnostics,
    performance = performance,
    recovery = recovery,
    weights = weights_long,
    parameters = parameters,
    tau = tau
  )
  save_rds_atomic(result, result_path, compress = FALSE)
  result
}

collect_simulation_results <- function(
    cfg, repetitions, scenarios = cfg$simulation$scenarios,
    cohort = "main") {
  results <- list()
  position <- 1L
  for (scenario in scenarios) {
    for (repetition in repetitions) {
      path <- file.path(
        simulation_repetition_dir(
          cfg, scenario, repetition, cohort
        ),
        "result.rds"
      )
      if (!file.exists(path)) next
      result <- tryCatch(readRDS(path), error = function(error) NULL)
      if (is.null(result) ||
          !identical(
            result$protocol_version,
            cfg$simulation$protocol_version
          ) ||
          !identical(
            result$cohort,
            cohort
          )) {
        next
      }
      results[[position]] <- result
      position <- position + 1L
    }
  }
  results
}

write_simulation_stage_results <- function(
    cfg, repetitions, stage,
    scenarios = cfg$simulation$scenarios) {
  results <- collect_simulation_results(
    cfg, repetitions, scenarios, cohort = stage
  )
  assert_true(length(results) > 0L, "Aucun resultat de simulation.")
  output_dir <- file.path(simulation_output_root(cfg), stage)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  diagnostics <- simulation_bind_rows(results, "diagnostics")
  attempts <- simulation_bind_rows(results, "diagnostics_attempts")
  contextual_optimization <- simulation_bind_rows(
    results, "contextual_optimization"
  )
  contextual_optimization_attempts <- simulation_bind_rows(
    results, "contextual_optimization_attempts"
  )
  performance <- simulation_bind_rows(results, "performance")
  recovery <- simulation_bind_rows(results, "recovery")
  weights <- simulation_bind_rows(results, "weights")
  parameters <- simulation_bind_rows(results, "parameters")
  tau <- simulation_bind_rows(results, "tau")

  write_csv_atomic(
    diagnostics,
    file.path(output_dir, "diagnostics_by_repetition.csv")
  )
  write_csv_atomic(
    attempts,
    file.path(output_dir, "diagnostics_all_attempts.csv")
  )
  if (nrow(contextual_optimization)) {
    write_csv_atomic(
      contextual_optimization,
      file.path(output_dir, "contextual_optimization_by_repetition.csv")
    )
    write_csv_atomic(
      contextual_optimization_attempts,
      file.path(output_dir, "contextual_optimization_all_attempts.csv")
    )
  }
  if (nrow(performance)) {
    write_csv_atomic(
      performance,
      file.path(output_dir, "performance_by_repetition.csv")
    )
    write_csv_atomic(
      recovery,
      file.path(output_dir, "weight_recovery_by_repetition.csv")
    )
    write_csv_atomic(
      weights,
      file.path(output_dir, "weights_by_repetition.csv")
    )
    write_csv_atomic(
      parameters,
      file.path(output_dir, "parameters_by_repetition.csv")
    )
    write_csv_atomic(
      tau,
      file.path(output_dir, "tau_by_repetition.csv")
    )
  }

  scenario_rows <- lapply(scenarios, function(scenario) {
    part <- diagnostics[diagnostics$scenario == scenario, , drop = FALSE]
    data.frame(
      scenario = scenario,
      repetitions_planned = length(repetitions),
      repetitions_completed = nrow(part),
      repetitions_valid = sum(part$status == "valide"),
      reestimations = sum(part$configuration != "legere"),
      final_failures = sum(part$status == "echec_definitif"),
      total_elapsed_seconds = sum(
        attempts$elapsed_seconds[attempts$scenario == scenario],
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })
  scenario_summary <- do.call(rbind, scenario_rows)
  write_csv_atomic(
    scenario_summary,
    file.path(output_dir, "diagnostics_by_scenario.csv")
  )

  if (stage == "pilot") {
    light <- attempts[
      attempts$configuration == "legere", , drop = FALSE
    ]
    pilot_assessment <- data.frame(
      adjustments_expected = length(scenarios) * length(repetitions),
      adjustments_completed = nrow(light),
      light_valid = sum(light$valid),
      light_failures = sum(!light$valid),
      divergences = sum(light$divergences),
      treedepth_hits = sum(light$max_treedepth_hits),
      stuck_chains = sum(light$chain_stuck),
      worst_rhat = max(light$max_rhat, na.rm = TRUE),
      minimum_ess_bulk = min(light$min_ess_bulk, na.rm = TRUE),
      minimum_ess_tail = min(light$min_ess_tail, na.rm = TRUE),
      minimum_ebfmi = min(light$min_ebfmi, na.rm = TRUE),
      all_light_adjustments_valid =
        nrow(light) == length(scenarios) * length(repetitions) &&
          all(light$valid),
      decision = "a_examiner_avant_etude_principale",
      stringsAsFactors = FALSE
    )
    write_csv_atomic(
      pilot_assessment,
      file.path(output_dir, "pilot_assessment.csv")
    )
  }
  invisible(list(
    results = results,
    diagnostics = diagnostics,
    attempts = attempts,
    contextual_optimization = contextual_optimization,
    contextual_optimization_attempts = contextual_optimization_attempts,
    performance = performance,
    recovery = recovery,
    weights = weights,
    parameters = parameters,
    tau = tau,
    scenario_summary = scenario_summary
  ))
}

select_simulation_validation_repetitions <- function(cfg) {
  main_dir <- file.path(simulation_output_root(cfg), "main")
  diagnostics <- utils::read.csv(
    file.path(main_dir, "diagnostics_by_repetition.csv"),
    stringsAsFactors = FALSE
  )
  attempts <- utils::read.csv(
    file.path(main_dir, "diagnostics_all_attempts.csv"),
    stringsAsFactors = FALSE
  )
  recovery <- utils::read.csv(
    file.path(main_dir, "weight_recovery_by_repetition.csv"),
    stringsAsFactors = FALSE
  )
  recovery <- recovery[
    recovery$method == "stacking_hierarchical", ,
    drop = FALSE
  ]
  initial_light <- attempts[
    attempts$configuration == "legere",
    c(
      "scenario", "repetition", "max_rhat", "min_ess_bulk",
      "min_ess_tail", "divergences", "max_treedepth_hits",
      "min_ebfmi", "chain_stuck", "valid"
    ),
    drop = FALSE
  ]
  names(initial_light)[-(1:2)] <- paste0(
    "initial_", names(initial_light)[-(1:2)]
  )
  retry_keys <- unique(attempts[
    attempts$configuration != "legere",
    c("scenario", "repetition"),
    drop = FALSE
  ])
  retry_keys$prior_retry <- TRUE
  eligible <- merge(
    diagnostics[
      diagnostics$status == "valide",
      c("scenario", "repetition", "configuration"),
      drop = FALSE
    ],
    recovery[c("scenario", "repetition", "weight_rmse")],
    by = c("scenario", "repetition")
  )
  eligible <- merge(
    eligible, initial_light,
    by = c("scenario", "repetition"), all.x = TRUE
  )
  eligible <- merge(
    eligible, retry_keys,
    by = c("scenario", "repetition"), all.x = TRUE
  )
  eligible$prior_retry[is.na(eligible$prior_retry)] <- FALSE
  names(eligible)[names(eligible) == "configuration"] <-
    "main_configuration"
  counts <- c(
    constant_weights = 2L,
    horizon_only = 2L,
    age_horizon = 2L,
    low_information = 4L
  )

  # La selection est deterministe. Elle combine un cas numeriquement
  # defavorable (en priorite une repetition deja relancee) et des cas situes
  # a plusieurs niveaux de RMSE. Le score de stress repose uniquement sur la
  # premiere tentative legere, afin de ne pas masquer une difficulte corrigee
  # par la relance complete.
  selected <- lapply(names(counts), function(scenario) {
    part <- eligible[eligible$scenario == scenario, , drop = FALSE]
    part <- part[order(part$repetition), , drop = FALSE]
    required <- counts[[scenario]]
    assert_true(
      nrow(part) >= required,
      paste(
        "Pas assez de repetitions valides pour la validation :",
        scenario
      )
    )

    n_part <- nrow(part)
    scaled_rank <- function(values, decreasing = FALSE) {
      ranked <- rank(
        if (decreasing) -values else values,
        ties.method = "average", na.last = "keep"
      )
      ifelse(is.na(ranked), 0, ranked / n_part)
    }
    part$diagnostic_stress_score <-
      scaled_rank(part$initial_max_rhat) +
      scaled_rank(part$initial_min_ess_bulk, decreasing = TRUE) +
      scaled_rank(part$initial_min_ess_tail, decreasing = TRUE) +
      scaled_rank(part$initial_min_ebfmi, decreasing = TRUE) +
      2 * as.numeric(part$initial_divergences > 0) +
      2 * as.numeric(part$initial_max_treedepth_hits > 0) +
      2 * as.numeric(part$initial_chain_stuck) +
      2 * as.numeric(!part$initial_valid)

    retry_candidates <- which(part$prior_retry)
    if (length(retry_candidates)) {
      first <- retry_candidates[
        which.max(part$diagnostic_stress_score[retry_candidates])
      ]
      first_reason <- "relance_precedente_et_diagnostic_defavorable"
    } else {
      first <- which.max(part$diagnostic_stress_score)
      first_reason <- "diagnostic_numerique_defavorable"
    }
    chosen <- first
    reasons <- first_reason
    targets_used <- NA_real_

    remaining <- required - 1L
    if (remaining > 0L) {
      if (remaining == 1L) {
        median_rmse <- stats::median(part$weight_rmse)
        probabilities <- if (
          part$weight_rmse[first] <= median_rmse
        ) 0.80 else 0.20
      } else {
        probabilities <- seq(0.20, 0.80, length.out = remaining)
      }
      targets <- stats::quantile(
        part$weight_rmse, probabilities, names = FALSE
      )
      for (target_index in seq_along(targets)) {
        available <- setdiff(seq_len(nrow(part)), chosen)
        next_index <- available[
          which.min(abs(part$weight_rmse[available] - targets[target_index]))
        ]
        chosen <- c(chosen, next_index)
        reasons <- c(
          reasons,
          sprintf(
            "couverture_rmse_q%02d",
            round(100 * probabilities[target_index])
          )
        )
        targets_used <- c(targets_used, probabilities[target_index])
      }
    }

    result <- part[chosen, , drop = FALSE]
    result$selection_reason <- reasons
    result$selection_target <- targets_used
    result <- result[
      ,
      c(
        "scenario", "repetition", "main_configuration", "prior_retry",
        "weight_rmse", "selection_reason", "selection_target",
        "diagnostic_stress_score", "initial_max_rhat",
        "initial_min_ess_bulk", "initial_min_ess_tail",
        "initial_divergences", "initial_max_treedepth_hits",
        "initial_min_ebfmi", "initial_chain_stuck", "initial_valid"
      ),
      drop = FALSE
    ]
    result
  })
  result <- do.call(rbind, selected)
  rownames(result) <- NULL
  result
}

run_simulation_validation_repetition <- function(
    scenario, repetition, compiled_model, cfg, reference = NULL) {
  validation_dir <- file.path(
    simulation_output_root(cfg),
    "validation",
    scenario,
    sprintf("rep_%03d", repetition)
  )
  dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
  result_path <- file.path(validation_dir, "result.rds")
  validation_version <- paste0(
    cfg$simulation$protocol_version, "_full_validation_v2"
  )
  if (file.exists(result_path) && !force_recompute()) {
    cached <- tryCatch(readRDS(result_path), error = function(error) NULL)
    if (!is.null(cached) &&
        identical(cached$protocol_version, validation_version)) {
      return(cached)
    }
  }

  original_path <- file.path(
    simulation_repetition_dir(
      cfg, scenario, repetition, cohort = "main"
    ),
    "result.rds"
  )
  assert_true(
    file.exists(original_path),
    paste("Resultat leger absent :", scenario, repetition)
  )
  original <- readRDS(original_path)
  assert_true(
    original$status == "valide" &&
      original$configuration %in%
        c("legere", "complete", "complete_treedepth15"),
    paste(
      "La validation complete requiert une repetition principale valide :",
      scenario, repetition
    )
  )
  reference <- reference %||% load_or_build_simulation_reference(cfg)
  components <- simulation_predictive_components(reference, scenario, cfg)
  validation <- simulate_study_dataset(
    scenario,
    simulation_repetition_seed(
      cfg, scenario, repetition, "validation", cohort = "main"
    ),
    cfg,
    reference = reference,
    components = components
  )
  test <- simulate_study_dataset(
    scenario,
    simulation_repetition_seed(
      cfg, scenario, repetition, "test", cohort = "main"
    ),
    cfg,
    reference = reference,
    constants = validation$constants,
    components = components
  )
  attempts <- list(simulation_attempt(
    scenario,
    repetition,
    "validation_complete",
    validation,
    compiled_model,
    cfg$simulation$full_mcmc,
    cfg,
    cohort = "validation"
  ))
  if (!attempts[[1L]]$overview$valid &&
      attempts[[1L]]$overview$max_treedepth_hits > 0) {
    attempts[[2L]] <- simulation_attempt(
      scenario,
      repetition,
      "validation_complete_treedepth15",
      validation,
      compiled_model,
      cfg$simulation$deep_mcmc,
      cfg,
      cohort = "validation"
    )
  }
  valid <- which(vapply(
    attempts,
    function(attempt) isTRUE(attempt$overview$valid),
    logical(1)
  ))
  diagnostics <- do.call(rbind, lapply(attempts, `[[`, "overview"))
  write_csv_atomic(
    diagnostics,
    file.path(validation_dir, "diagnostics.csv")
  )
  if (!length(valid)) {
    result <- list(
      protocol_version = validation_version,
      scenario = scenario,
      repetition = repetition,
      status = "echec_definitif",
      diagnostics = diagnostics
    )
    save_rds_atomic(result, result_path, compress = FALSE)
    return(result)
  }
  selected <- attempts[[valid[[1L]]]]
  grid <- simulation_grid(cfg)
  weight_draws <- simulation_weight_draw_array(
    selected$hierarchical$fit,
    grid,
    cfg$models,
    validation$constants,
    ndraws = NULL,
    seed = simulation_repetition_seed(
      cfg, scenario, repetition, "validation_full_weights",
      cohort = "validation"
    )
  )
  true_weights <- simulation_true_weights(
    scenario,
    standardize_context(grid, validation$constants)$data,
    cfg$models
  )
  full_weights <- simulation_weight_summary(
    weight_draws, grid, cfg$models, true_weights
  )
  full_weights$method <- "stacking_hierarchical"
  full_recovery <- simulation_recovery_row(full_weights)
  full_performance <- simulation_evaluate_method(
    test,
    scenario,
    repetition,
    "stacking_hierarchical",
    cfg$models,
    cfg$simulation$predictive_draws,
    simulation_repetition_seed(
      cfg, scenario, repetition, "validation_full_predictive",
      cohort = "validation"
    ),
    posterior_weights = weight_draws
  )
  full_parameters <- simulation_parameter_summary(
    selected$hierarchical$fit
  )

  baseline_weights <- original$weights[
    original$weights$method == "stacking_hierarchical", ,
    drop = FALSE
  ]
  key_baseline <- paste(
    baseline_weights$age,
    baseline_weights$horizon,
    baseline_weights$model
  )
  key_full <- paste(
    full_weights$age,
    full_weights$horizon,
    full_weights$model
  )
  full_weights <- full_weights[
    match(key_baseline, key_full), ,
    drop = FALSE
  ]
  baseline_recovery <- original$recovery[
    original$recovery$method == "stacking_hierarchical", ,
    drop = FALSE
  ]
  baseline_performance <- original$performance[
    original$performance$method == "stacking_hierarchical", ,
    drop = FALSE
  ]
  baseline_parameters <- original$parameters
  parameter_key <- baseline_parameters$parameter
  full_parameters <- full_parameters[
    match(parameter_key, full_parameters$parameter), ,
    drop = FALSE
  ]
  contextual_parameters <- grepl(
    "^(alpha|beta)\\[", parameter_key
  )
  tau_parameters <- grepl("^tau\\[", parameter_key)

  original_ranks_logs <- original$performance$method[
    order(original$performance$mean_logs)
  ]
  original_ranks_crps <- original$performance$method[
    order(original$performance$mean_crps)
  ]
  full_method_performance <- original$performance
  hierarchical_row <- full_method_performance$method ==
    "stacking_hierarchical"
  full_method_performance$mean_logs[hierarchical_row] <-
    full_performance$mean_logs
  full_method_performance$mean_crps[hierarchical_row] <-
    full_performance$mean_crps

  comparison <- data.frame(
    scenario = scenario,
    repetition = repetition,
    comparison_type = if (
      original$configuration == "legere"
    ) "legere_vs_complete" else "relance_complete_vs_complete",
    baseline_configuration = original$configuration,
    configuration_full = selected$overview$configuration,
    mean_weight_rmse_between = sqrt(mean(
      (baseline_weights$estimated_weight -
         full_weights$estimated_weight)^2
    )),
    mean_weight_absolute_difference = mean(abs(
      baseline_weights$estimated_weight -
        full_weights$estimated_weight
    )),
    q025_absolute_difference = mean(abs(
      baseline_weights$weight_q025 - full_weights$weight_q025
    )),
    q975_absolute_difference = mean(abs(
      baseline_weights$weight_q975 - full_weights$weight_q975
    )),
    contextual_coefficient_mean_absolute_difference = mean(abs(
      baseline_parameters$mean[contextual_parameters] -
        full_parameters$mean[contextual_parameters]
    )),
    tau_mean_absolute_difference = mean(abs(
      baseline_parameters$mean[tau_parameters] -
        full_parameters$mean[tau_parameters]
    )),
    weight_rmse_baseline = baseline_recovery$weight_rmse,
    weight_rmse_full = full_recovery$weight_rmse,
    logs_baseline = baseline_performance$mean_logs,
    logs_full = full_performance$mean_logs,
    crps_baseline = baseline_performance$mean_crps,
    crps_full = full_performance$mean_crps,
    relative_logs_difference = abs(
      full_performance$mean_logs - baseline_performance$mean_logs
    ) / max(abs(baseline_performance$mean_logs), .Machine$double.eps),
    relative_crps_difference = abs(
      full_performance$mean_crps - baseline_performance$mean_crps
    ) / max(abs(baseline_performance$mean_crps), .Machine$double.eps),
    extreme_frequency_baseline =
      baseline_recovery$extreme_weight_proportion,
    extreme_frequency_full =
      full_recovery$extreme_weight_proportion,
    coverage_baseline = baseline_recovery$posterior_95_coverage,
    coverage_full = full_recovery$posterior_95_coverage,
    interval_width_baseline = baseline_recovery$mean_interval_width,
    interval_width_full = full_recovery$mean_interval_width,
    logs_ranking_identical = identical(
      original_ranks_logs,
      full_method_performance$method[
        order(full_method_performance$mean_logs)
      ]
    ),
    crps_ranking_identical = identical(
      original_ranks_crps,
      full_method_performance$method[
        order(full_method_performance$mean_crps)
      ]
    ),
    stringsAsFactors = FALSE
  )
  write_csv_atomic(
    comparison, file.path(validation_dir, "comparison.csv")
  )
  write_csv_atomic(
    full_weights, file.path(validation_dir, "weights_full.csv")
  )
  write_csv_atomic(
    full_parameters, file.path(validation_dir, "parameters_full.csv")
  )
  result <- list(
    protocol_version = validation_version,
    scenario = scenario,
    repetition = repetition,
    status = "valide",
    diagnostics = selected$overview,
    comparison = comparison
  )
  save_rds_atomic(result, result_path, compress = FALSE)
  result
}

run_simulation_validation <- function(
    cfg, compiled_model, reference = NULL) {
  reference <- reference %||% load_or_build_simulation_reference(cfg)
  selection <- select_simulation_validation_repetitions(cfg)
  validation_root <- file.path(simulation_output_root(cfg), "validation")
  dir.create(validation_root, recursive = TRUE, showWarnings = FALSE)
  write_csv_atomic(
    selection, file.path(validation_root, "selection.csv")
  )
  results <- lapply(seq_len(nrow(selection)), function(index) {
    run_simulation_validation_repetition(
      selection$scenario[index],
      selection$repetition[index],
      compiled_model,
      cfg,
      reference = reference
    )
  })
  diagnostics <- simulation_bind_rows(results, "diagnostics")
  comparisons <- simulation_bind_rows(results, "comparison")
  write_csv_atomic(
    diagnostics,
    file.path(validation_root, "diagnostics.csv")
  )
  write_csv_atomic(
    comparisons,
    file.path(validation_root, "comparisons.csv")
  )
  light_comparisons <- comparisons[
    comparisons$comparison_type == "legere_vs_complete", ,
    drop = FALSE
  ]
  retry_comparisons <- comparisons[
    comparisons$comparison_type == "relance_complete_vs_complete", ,
    drop = FALSE
  ]
  thresholds <- cfg$simulation$validation_thresholds
  repetitions_valid <- sum(vapply(
    results, function(result) result$status == "valide", logical(1)
  ))
  logs_ranking_agreement <- mean(
    light_comparisons$logs_ranking_identical
  )
  crps_ranking_agreement <- mean(
    light_comparisons$crps_ranking_identical
  )
  diagnostics_stable <-
    repetitions_valid == nrow(selection) &&
    all(diagnostics$valid) &&
    all(diagnostics$divergences == 0) &&
    all(diagnostics$max_treedepth_hits == 0) &&
    !any(diagnostics$chain_stuck) &&
    max(diagnostics$max_rhat) <= cfg$simulation$rhat_max &&
    min(diagnostics$min_ess_bulk) >= cfg$simulation$ess_min &&
    min(diagnostics$min_ess_tail) >= cfg$simulation$ess_min &&
    min(diagnostics$min_ebfmi) >= cfg$simulation$ebfmi_min
  weights_stable <-
    mean(light_comparisons$mean_weight_absolute_difference) <=
      thresholds$mean_weight_absolute_difference &&
    max(light_comparisons$mean_weight_absolute_difference) <=
      thresholds$maximum_weight_absolute_difference
  scores_stable <-
    mean(light_comparisons$relative_logs_difference) <=
      thresholds$mean_relative_logs_difference &&
    mean(light_comparisons$relative_crps_difference) <=
      thresholds$mean_relative_crps_difference
  rankings_stable <-
    logs_ranking_agreement >= thresholds$minimum_ranking_agreement &&
    crps_ranking_agreement >= thresholds$minimum_ranking_agreement
  validation_passed <-
    diagnostics_stable && weights_stable &&
    scores_stable && rankings_stable
  summary <- data.frame(
    repetitions_planned = nrow(selection),
    repetitions_valid = repetitions_valid,
    light_to_full_comparisons = nrow(light_comparisons),
    retry_full_to_full_controls = nrow(retry_comparisons),
    mean_weight_absolute_difference =
      mean(light_comparisons$mean_weight_absolute_difference),
    maximum_weight_absolute_difference =
      max(light_comparisons$mean_weight_absolute_difference),
    mean_relative_logs_difference =
      mean(light_comparisons$relative_logs_difference),
    mean_relative_crps_difference =
      mean(light_comparisons$relative_crps_difference),
    logs_ranking_agreement = logs_ranking_agreement,
    crps_ranking_agreement = crps_ranking_agreement,
    maximum_full_rhat = max(diagnostics$max_rhat),
    minimum_full_ess_bulk = min(diagnostics$min_ess_bulk),
    minimum_full_ess_tail = min(diagnostics$min_ess_tail),
    full_divergences = sum(diagnostics$divergences),
    full_treedepth_hits = sum(diagnostics$max_treedepth_hits),
    minimum_full_ebfmi = min(diagnostics$min_ebfmi),
    diagnostics_stable = diagnostics_stable,
    weights_stable = weights_stable,
    scores_stable = scores_stable,
    rankings_stable = rankings_stable,
    threshold_mean_weight_absolute_difference =
      thresholds$mean_weight_absolute_difference,
    threshold_maximum_weight_absolute_difference =
      thresholds$maximum_weight_absolute_difference,
    threshold_mean_relative_logs_difference =
      thresholds$mean_relative_logs_difference,
    threshold_mean_relative_crps_difference =
      thresholds$mean_relative_crps_difference,
    threshold_minimum_ranking_agreement =
      thresholds$minimum_ranking_agreement,
    decision = if (
      validation_passed
    ) "configuration_legere_stable" else "examen_requis",
    stringsAsFactors = FALSE
  )
  write_csv_atomic(
    summary, file.path(validation_root, "validation_summary.csv")
  )
  invisible(list(
    selection = selection,
    diagnostics = diagnostics,
    comparisons = comparisons,
    summary = summary
  ))
}

summarize_simulation_main_results <- function(cfg) {
  main_dir <- file.path(simulation_output_root(cfg), "main")
  performance <- utils::read.csv(
    file.path(main_dir, "performance_by_repetition.csv"),
    stringsAsFactors = FALSE
  )
  recovery <- utils::read.csv(
    file.path(main_dir, "weight_recovery_by_repetition.csv"),
    stringsAsFactors = FALSE
  )
  diagnostics <- utils::read.csv(
    file.path(main_dir, "diagnostics_by_repetition.csv"),
    stringsAsFactors = FALSE
  )
  weights <- utils::read.csv(
    file.path(main_dir, "weights_by_repetition.csv"),
    stringsAsFactors = FALSE
  )

  performance_summary <- do.call(rbind, lapply(
    split(performance, list(performance$scenario, performance$method)),
    function(part) {
      data.frame(
        scenario = part$scenario[1L],
        method = part$method[1L],
        repetitions = nrow(part),
        mean_logs = mean(part$mean_logs),
        sd_logs = stats::sd(part$mean_logs),
        mean_crps = mean(part$mean_crps),
        sd_crps = stats::sd(part$mean_crps),
        stringsAsFactors = FALSE
      )
    }
  ))
  performance_summary$rank_logs <- ave(
    performance_summary$mean_logs,
    performance_summary$scenario,
    FUN = function(x) rank(x, ties.method = "min")
  )
  performance_summary$rank_crps <- ave(
    performance_summary$mean_crps,
    performance_summary$scenario,
    FUN = function(x) rank(x, ties.method = "min")
  )

  recovery_summary <- do.call(rbind, lapply(
    split(recovery, list(recovery$scenario, recovery$method)),
    function(part) {
      data.frame(
        scenario = part$scenario[1L],
        method = part$method[1L],
        repetitions = nrow(part),
        weight_bias = mean(part$weight_bias),
        weight_rmse = mean(part$weight_rmse),
        mean_absolute_weight_error =
          mean(part$mean_absolute_weight_error),
        weight_rmse_mean = mean(part$weight_rmse),
        weight_rmse_median = stats::median(part$weight_rmse),
        weight_rmse_q025 = stats::quantile(
          part$weight_rmse, 0.025, names = FALSE
        ),
        weight_rmse_q975 = stats::quantile(
          part$weight_rmse, 0.975, names = FALSE
        ),
        weight_bias_mean = mean(part$weight_bias),
        extreme_weight_frequency =
          mean(part$extreme_weight_proportion),
        posterior_95_coverage =
          mean(part$posterior_95_coverage, na.rm = TRUE),
        mean_interval_width =
          mean(part$mean_interval_width, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  ))
  recovery_summary$posterior_95_coverage[
    is.nan(recovery_summary$posterior_95_coverage)
  ] <- NA_real_
  recovery_summary$mean_interval_width[
    is.nan(recovery_summary$mean_interval_width)
  ] <- NA_real_
  surface_summary <- stats::aggregate(
    weights[c("true_weight", "estimated_weight")],
    by = weights[c(
      "scenario", "method", "age", "horizon", "model"
    )],
    FUN = mean
  )

  write_csv_atomic(
    performance_summary,
    file.path(simulation_output_root(cfg), "simulation_summary.csv")
  )
  write_csv_atomic(
    recovery_summary,
    file.path(
      simulation_output_root(cfg), "simulation_weight_recovery.csv"
    )
  )
  write_csv_atomic(
    diagnostics,
    file.path(simulation_output_root(cfg), "simulation_diagnostics.csv")
  )
  write_csv_atomic(
    surface_summary[
      surface_summary$method == "stacking_hierarchical", ,
      drop = FALSE
    ],
    file.path(
      simulation_output_root(cfg),
      "simulation_hierarchical_weights.csv"
    )
  )
  save_rds_atomic(
    list(
      performance = performance_summary,
      recovery = recovery_summary,
      diagnostics = diagnostics,
      surfaces = surface_summary
    ),
    file.path(simulation_output_root(cfg), "simulation_results.rds"),
    compress = FALSE
  )

  methods <- c(
    "stacking_global",
    "stacking_contextual",
    "stacking_hierarchical"
  )
  scenarios <- cfg$simulation$scenarios
  figures_dir <- if (cfg$profile == "full") {
    file.path(cfg$root, "output", "pdf", "figures")
  } else {
    file.path(cfg$paths$figures, "simulation_smoke")
  }
  render_chapter_figure(
    "simulation_weight_recovery",
    function() {
      graphics::par(mfrow = c(1, 2), mar = c(8, 4.5, 2, 1))
      rmse <- xtabs(
        weight_rmse ~ method + scenario,
        data = recovery_summary
      )[methods, scenarios]
      graphics::barplot(
        rmse,
        beside = TRUE,
        col = chapter_model_colors()[methods],
        border = NA,
        ylab = "RMSE moyenne des poids",
        las = 2,
        cex.names = 0.72,
        main = "Recuperation des poids"
      )
      logs <- xtabs(
        mean_logs ~ method + scenario,
        data = performance_summary
      )[methods, scenarios]
      graphics::barplot(
        logs,
        beside = TRUE,
        col = chapter_model_colors()[methods],
        border = NA,
        ylab = "LogS moyen",
        las = 2,
        cex.names = 0.72,
        main = "Performance predictive"
      )
      graphics::legend(
        "topright",
        legend = chapter_method_labels()[methods],
        fill = chapter_model_colors()[methods],
        border = NA,
        bty = "n",
        cex = 0.70
      )
    },
    figures_dir,
    width = 11,
    height = 6
  )
  render_chapter_figure(
    "simulation_coverage_and_extremes",
    function() {
      graphics::par(mfrow = c(1, 2), mar = c(8, 4.5, 2, 1))
      hierarchical <- recovery_summary[
        recovery_summary$method == "stacking_hierarchical", ,
        drop = FALSE
      ]
      hierarchical <- hierarchical[
        match(scenarios, hierarchical$scenario), ,
        drop = FALSE
      ]
      graphics::barplot(
        hierarchical$posterior_95_coverage,
        names.arg = hierarchical$scenario,
        col = "#228833",
        border = NA,
        ylim = c(0, 1),
        las = 2,
        cex.names = 0.72,
        ylab = "Couverture moyenne",
        main = "Intervalles credibles a 95 %"
      )
      graphics::abline(h = 0.95, lty = 2)
      extremes <- xtabs(
        extreme_weight_frequency ~ method + scenario,
        data = recovery_summary
      )[methods, scenarios]
      graphics::barplot(
        extremes,
        beside = TRUE,
        col = chapter_model_colors()[methods],
        border = NA,
        las = 2,
        cex.names = 0.72,
        ylab = "Frequence des poids extremes",
        main = "Poids inferieurs a 0,01 ou superieurs a 0,99"
      )
    },
    figures_dir,
    width = 11,
    height = 6
  )
  render_chapter_figure(
    "simulation_weight_error_distributions",
    function() {
      errors <- weights$estimated_weight - weights$true_weight
      labels <- paste(weights$scenario, weights$method, sep = "\n")
      graphics::par(mar = c(11, 4.5, 2, 1))
      graphics::boxplot(
        split(errors, labels),
        las = 2,
        outline = FALSE,
        col = "#BBCCEE",
        border = "#4477AA",
        ylab = "Poids estime - poids vrai",
        main = "Distribution des erreurs de poids"
      )
      graphics::abline(h = 0, lty = 2)
    },
    figures_dir,
    width = 13,
    height = 7
  )

  hierarchical_surface <- surface_summary[
    surface_summary$method == "stacking_hierarchical", ,
    drop = FALSE
  ]
  for (scenario in scenarios) {
    part <- hierarchical_surface[
      hierarchical_surface$scenario == scenario, ,
      drop = FALSE
    ]
    render_chapter_figure(
      paste0("simulation_weight_surfaces_", scenario),
      function() {
        graphics::par(
          mfrow = c(2, length(cfg$models)),
          mar = c(2.5, 2.5, 2.5, 1)
        )
        for (value in c("true_weight", "estimated_weight")) {
          for (model in cfg$models) {
            model_part <- part[part$model == model, ]
            matrix_values <- xtabs(
              model_part[[value]] ~
                model_part$age + model_part$horizon
            )
            graphics::image(
              x = as.numeric(rownames(matrix_values)),
              y = as.numeric(colnames(matrix_values)),
              z = matrix_values,
              col = grDevices::hcl.colors(30, "YlGnBu"),
              zlim = c(0, 1),
              xlab = "Age",
              ylab = "Horizon",
              main = paste(
                if (value == "true_weight") "Vrai" else "Estime",
                toupper(model)
              )
            )
          }
        }
      },
      figures_dir,
      width = 14,
      height = 6
    )
  }
  invisible(list(
    performance = performance_summary,
    recovery = recovery_summary,
    diagnostics = diagnostics,
    surfaces = surface_summary
  ))
}
