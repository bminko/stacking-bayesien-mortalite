# Pseudo-BMA, stacking global, stacking contextuel et BMA classique.

standardize_context <- function(data, constants = NULL) {
  if (is.null(constants)) {
    unique_ages <- sort(unique(data$age))
    unique_horizons <- sort(unique(data$horizon))
    constants <- list(
      age_mean = mean(unique_ages),
      age_sd = stats::sd(unique_ages),
      age_square_mean = mean(
        ((unique_ages - mean(unique_ages)) / stats::sd(unique_ages))^2
      ),
      horizon_mean = mean(unique_horizons),
      horizon_sd = stats::sd(unique_horizons)
    )
  }
  assert_true(constants$age_sd > 0, "Ecart-type d'age nul.")
  assert_true(constants$horizon_sd > 0, "Ecart-type d'horizon nul.")

  data$x_tilde <- (data$age - constants$age_mean) / constants$age_sd
  data$q_age <- data$x_tilde^2 - constants$age_square_mean
  data$h_tilde <- (data$horizon - constants$horizon_mean) /
    constants$horizon_sd
  data$x_h <- data$x_tilde * data$h_tilde
  list(data = data, constants = constants)
}

context_design_matrix <- function(
    data,
    columns = c("age", "age2", "horizon", "age_horizon")) {
  design <- cbind(
    age = data$x_tilde,
    age2 = data$q_age,
    horizon = data$h_tilde,
    age_horizon = data$x_h
  )
  unknown <- setdiff(columns, colnames(design))
  assert_true(
    !length(unknown),
    paste("Colonnes contextuelles inconnues :", paste(unknown, collapse = ", "))
  )
  design[, columns, drop = FALSE]
}

meta_logp_matrix <- function(meta, models) {
  columns <- paste0("log_p_", models)
  missing <- setdiff(columns, names(meta))
  assert_true(!length(missing), paste(
    "Colonnes predictives absentes :", paste(missing, collapse = ", ")
  ))
  matrix <- as.matrix(meta[columns])
  assert_true(all(is.finite(matrix)),
              "Le meta-jeu contient des log-densites non finies.")
  colnames(matrix) <- models
  matrix
}

scores_to_weights <- function(scores) {
  maxima <- apply(scores, 1L, max)
  exponentials <- exp(scores - maxima)
  exponentials / rowSums(exponentials)
}

theta_to_global_weights <- function(theta, K) {
  softmax(c(theta, 0))
}

select_best_optimization <- function(fits, label) {
  converged <- vapply(
    fits,
    function(fit) identical(fit$convergence, 0L) && is.finite(fit$value),
    logical(1)
  )
  assert_true(
    any(converged),
    paste("Aucune optimisation convergente pour", label)
  )
  candidates <- which(converged)
  best_index <- candidates[which.min(vapply(
    fits[candidates], `[[`, numeric(1), "value"
  ))]
  fits[[best_index]]
}

validate_weight_vector <- function(weights, label) {
  assert_true(
    all(is.finite(weights)) && all(weights >= 0),
    paste("Poids invalides pour", label)
  )
  assert_true(
    abs(sum(weights) - 1) < 1e-8,
    paste("Les poids ne somment pas a 1 pour", label)
  )
  invisible(weights)
}

stacking_responsibilities <- function(log_p, weights) {
  log_joint <- if (is.matrix(weights)) {
    log_p + log(weights)
  } else {
    sweep(log_p, 2L, log(weights), `+`)
  }
  exp(sweep(log_joint, 1L, row_log_sum_exp(log_joint), `-`))
}

global_stacking_value <- function(theta, log_p, omega) {
  weights <- theta_to_global_weights(theta, ncol(log_p))
  log_mixture <- row_log_sum_exp(
    sweep(log_p, 2L, log(weights), `+`)
  )
  -sum(omega * log_mixture)
}

global_stacking_gradient <- function(theta, log_p, omega) {
  K <- ncol(log_p)
  weights <- theta_to_global_weights(theta, K)
  responsibilities <- stacking_responsibilities(log_p, weights)
  difference <- sweep(-responsibilities, 2L, weights, `+`)
  colSums(sweep(difference, 1L, omega, `*`))[seq_len(K - 1L)]
}

fit_pseudo_bma <- function(meta, models) {
  log_p <- meta_logp_matrix(meta, models)
  omega <- meta$omega %||% rep(1, nrow(meta))
  elpd <- colSums(log_p * omega)
  weights <- softmax(elpd)
  names(weights) <- models
  validate_weight_vector(weights, "pseudo-BMA")
  list(weights = weights, elpd = elpd)
}

fit_global_stacking <- function(meta, models, multistarts = 8L,
                                seed = 1L) {
  log_p <- meta_logp_matrix(meta, models)
  omega <- meta$omega %||% rep(1, nrow(meta))
  K <- length(models)
  objective <- function(theta) {
    global_stacking_value(theta, log_p, omega)
  }
  gradient <- function(theta) {
    global_stacking_gradient(theta, log_p, omega)
  }

  set.seed(seed)
  starts <- c(
    list(rep(0, K - 1L)),
    replicate(multistarts - 1L, stats::rnorm(K - 1L), simplify = FALSE)
  )
  fits <- lapply(starts, function(start) {
    stats::optim(
      start, objective, gr = gradient, method = "BFGS",
      control = list(maxit = 3000L, reltol = 1e-11)
    )
  })
  values <- vapply(fits, `[[`, numeric(1), "value")
  best <- select_best_optimization(fits, "stacking global")
  weights <- theta_to_global_weights(best$par, K)
  names(weights) <- models
  validate_weight_vector(weights, "stacking global")
  list(weights = weights, optimization = best, all_values = values)
}

unpack_context_theta <- function(theta, K, P) {
  matrix(theta, nrow = K - 1L, ncol = P + 1L, byrow = TRUE)
}

context_weights_from_coefficients <- function(coefficients, X) {
  K_minus_one <- nrow(coefficients)
  scores_reference_free <- cbind(1, X) %*% t(coefficients)
  scores <- cbind(scores_reference_free, reference = 0)
  weights <- scores_to_weights(scores)
  colnames(weights) <- c(seq_len(K_minus_one), K_minus_one + 1L)
  weights
}

contextual_stacking_value <- function(theta, log_p, X, omega) {
  K <- ncol(log_p)
  P <- ncol(X)
  coefficients <- unpack_context_theta(theta, K, P)
  weights <- context_weights_from_coefficients(coefficients, X)
  log_mixture <- row_log_sum_exp(log(weights) + log_p)
  -sum(omega * log_mixture)
}

contextual_stacking_gradient <- function(theta, log_p, X, omega) {
  K <- ncol(log_p)
  P <- ncol(X)
  coefficients <- unpack_context_theta(theta, K, P)
  weights <- context_weights_from_coefficients(coefficients, X)
  responsibilities <- stacking_responsibilities(log_p, weights)
  difference <- sweep(weights - responsibilities, 1L, omega, `*`)
  gradient <- t(difference[, seq_len(K - 1L), drop = FALSE]) %*%
    cbind(intercept = 1, X)
  as.vector(t(gradient))
}

fit_contextual_stacking <- function(
    meta, models, multistarts = 12L, seed = 1L,
    design_columns = c("age", "age2", "horizon", "age_horizon"),
    allow_nonconvergence = FALSE) {
  log_p <- meta_logp_matrix(meta, models)
  X <- context_design_matrix(meta, design_columns)
  omega <- meta$omega %||% rep(1, nrow(meta))
  K <- length(models)
  P <- ncol(X)

  objective <- function(theta) {
    contextual_stacking_value(theta, log_p, X, omega)
  }
  gradient <- function(theta) {
    contextual_stacking_gradient(theta, log_p, X, omega)
  }

  set.seed(seed)
  dimension <- (K - 1L) * (P + 1L)
  starts <- c(
    list(rep(0, dimension)),
    replicate(
      multistarts - 1L,
      stats::rnorm(dimension, sd = 0.25),
      simplify = FALSE
    )
  )
  fits <- lapply(starts, function(start) {
    stats::optim(
      start, objective, gr = gradient, method = "BFGS",
      control = list(maxit = 5000L, reltol = 1e-10)
    )
  })
  values <- vapply(fits, `[[`, numeric(1), "value")
  optimization_diagnostics <- do.call(rbind, lapply(
    seq_along(fits),
    function(index) {
      fit <- fits[[index]]
      gradient_norm <- tryCatch(
        sqrt(sum(gradient(fit$par)^2)),
        error = function(error) NA_real_
      )
      data.frame(
        start = index,
        convergence = fit$convergence,
        value = fit$value,
        function_evaluations = unname(fit$counts[["function"]]),
        gradient_evaluations = unname(fit$counts[["gradient"]]),
        message = fit$message %||% "",
        max_abs_parameter = max(abs(fit$par)),
        gradient_norm = gradient_norm,
        stringsAsFactors = FALSE
      )
    }
  ))
  converged <- vapply(
    fits,
    function(fit) identical(fit$convergence, 0L) && is.finite(fit$value),
    logical(1)
  )
  if (!any(converged) && isTRUE(allow_nonconvergence)) {
    return(list(
      status = "echec_non_convergence",
      error_message =
        "Aucune optimisation convergente pour stacking contextuel",
      coefficients = NULL,
      design_columns = design_columns,
      optimization = NULL,
      all_values = values,
      optimization_diagnostics = optimization_diagnostics
    ))
  }
  best <- select_best_optimization(fits, "stacking contextuel")
  coefficients <- unpack_context_theta(best$par, K, P)
  rownames(coefficients) <- models[seq_len(K - 1L)]
  colnames(coefficients) <- c("intercept", design_columns)
  list(
    status = "converge",
    error_message = "",
    coefficients = coefficients,
    design_columns = design_columns,
    optimization = best,
    all_values = values,
    optimization_diagnostics = optimization_diagnostics
  )
}

contextual_weight_grid <- function(grid, contextual_fit, models,
                                   constants) {
  standardized <- standardize_context(grid, constants)$data
  design_columns <- contextual_fit$design_columns %||%
    c("age", "age2", "horizon", "age_horizon")
  X <- context_design_matrix(standardized, design_columns)
  weights <- context_weights_from_coefficients(
    contextual_fit$coefficients, X
  )
  colnames(weights) <- models
  cbind(standardized, as.data.frame(weights, check.names = FALSE))
}

fit_hierarchical_stacking <- function(
    meta, models, cfg, compiled_model = NULL,
    strict = cfg$profile == "full", context = "aggregation",
    design_columns = c("age", "age2", "horizon", "age_horizon"),
    tau_scale = NULL, mcmc = cfg$hierarchical_mcmc) {
  if (is.null(compiled_model)) {
    compiled_model <- compile_stan_model("hierarchical", cfg)
  }
  X <- context_design_matrix(meta, design_columns)
  if (is.null(tau_scale)) {
    default_tau <- c(age = 1, age2 = 0.5, horizon = 1, age_horizon = 0.5)
    tau_scale <- unname(default_tau[design_columns])
  }
  assert_true(
    length(tau_scale) == ncol(X) && all(is.finite(tau_scale)) &&
      all(tau_scale > 0),
    "Les echelles a priori de tau doivent etre positives et de taille P."
  )
  stan_data <- list(
    N = nrow(meta),
    K = length(models),
    P = ncol(X),
    log_p = unname(meta_logp_matrix(meta, models)),
    X = unname(X),
    obs_weight = as.numeric(meta$omega),
    tau_scale = as.numeric(tau_scale)
  )
  output_dir <- file.path(cfg$paths$fits, context, "hierarchical")
  cache_key <- stan_fit_cache_key(
    "hierarchical", cfg, stan_data, mcmc
  )
  fit_with_diagnostics(
    compiled_model = compiled_model,
    stan_data = stan_data,
    output_dir = output_dir,
    mcmc = mcmc,
    seed = stable_seed(cfg$seed, cfg$profile, "hierarchical", context),
    label = paste(context, "hierarchical", sep = "/"),
    cache_key = cache_key,
    strict = strict
  )
}

hierarchical_weight_grid <- function(
    grid, hierarchical_fit, models, constants, ndraws = NULL,
    seed = 1L,
    design_columns = c("age", "age2", "horizon", "age_horizon")) {
  standardized <- standardize_context(grid, constants)$data
  X <- context_design_matrix(standardized, design_columns)
  K <- length(models)
  P <- ncol(X)
  draws <- fit_draw_matrix(
    hierarchical_fit,
    variables = c("alpha", "beta", "tau"),
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

  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    scores <- alpha
    for (p in seq_len(P)) {
      scores <- scores + beta[, , p] * X[i, p]
    }
    weights <- scores_to_weights(cbind(scores, 0))
    rows[[i]] <- data.frame(
      age = grid$age[i],
      horizon = grid$horizon[i],
      model = rep(models, each = 1L),
      weight_mean = colMeans(weights),
      weight_q025 = apply(
        weights, 2L, stats::quantile, probs = 0.025
      ),
      weight_q05 = apply(weights, 2L, stats::quantile, probs = 0.05),
      weight_q50 = apply(weights, 2L, stats::quantile, probs = 0.50),
      weight_q95 = apply(weights, 2L, stats::quantile, probs = 0.95),
      weight_q975 = apply(
        weights, 2L, stats::quantile, probs = 0.975
      ),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

weight_matrix_from_long <- function(weight_long, grid, models,
                                    value = "weight_mean") {
  key_grid <- paste(grid$age, grid$horizon, sep = "|")
  result <- matrix(NA_real_, nrow(grid), length(models),
                   dimnames = list(NULL, models))
  for (k in seq_along(models)) {
    part <- weight_long[weight_long$model == models[k], ]
    key_part <- paste(part$age, part$horizon, sep = "|")
    result[, k] <- part[[value]][match(key_grid, key_part)]
  }
  assert_true(!anyNA(result), "Poids contextuels manquants sur la grille.")
  result
}

prepare_rstan_fit_for_bridge <- function(fit, model, training) {
  if (!is_rstan_fit(fit)) return(fit)

  # Sous Windows, le pointeur C++ d'un stanfit n'est plus valide apres
  # saveRDS/readRDS. Une instance vide, construite avec les memes donnees,
  # restaure ce pointeur sans refaire l'echantillonnage.
  model_shell <- suppressMessages(rstan::sampling(
    fit@stanmodel,
    data = stan_data_for_model(training, model),
    chains = 0L,
    iter = 1L,
    warmup = 0L,
    refresh = 0L
  ))
  fit@.MISC <- model_shell@.MISC
  assert_true(
    rstan:::is_sfinstance_valid(fit),
    paste("Instance RStan invalide pour le bridge sampling de", model)
  )
  fit
}

compute_bma_weights <- function(fits, models, cfg, training = NULL,
                                strict = cfg$profile == "full") {
  require_packages("bridgesampling")
  log_marginal <- setNames(rep(NA_real_, length(models)), models)
  bridges <- vector("list", length(models))
  names(bridges) <- models
  error_message <- NULL

  for (model in models) {
    message_step("Bridge sampling : ", model)
    bridge_fit <- fits[[model]]
    if (is_rstan_fit(bridge_fit)) {
      assert_true(
        !is.null(training),
        "Les donnees d'apprentissage sont requises pour le BMA avec RStan."
      )
      bridge_fit <- prepare_rstan_fit_for_bridge(
        bridge_fit, model, training
      )
    }
    bridge <- tryCatch(
      bridgesampling::bridge_sampler(
        bridge_fit,
        repetitions = cfg$bridge_repetitions,
        method = "warp3",
        cores = if (.Platform$OS.type == "windows") {
          1L
        } else {
          min(cfg$mcmc$parallel_chains, cfg$mcmc$chains)
        },
        silent = FALSE
      ),
      error = function(error) error
    )
    if (inherits(bridge, "error")) {
      if (strict) {
        stop("Bridge sampling echoue pour ", model, " : ",
             conditionMessage(bridge), call. = FALSE)
      }
      error_message <- paste(
        "Bridge sampling indisponible pour", model, ":",
        conditionMessage(bridge)
      )
      warning(error_message, call. = FALSE)
      return(list(
        available = FALSE,
        weights = setNames(rep(NA_real_, length(models)), models),
        log_marginal = log_marginal,
        bridges = bridges,
        error = error_message
      ))
    }
    bridges[[model]] <- bridge
    log_marginal[model] <- as.numeric(bridgesampling::logml(bridge))
  }

  weights <- softmax(log_marginal)
  names(weights) <- models
  validate_weight_vector(weights, "BMA")
  list(
    available = TRUE,
    weights = weights,
    log_marginal = log_marginal,
    bridges = bridges,
    error = NULL
  )
}
