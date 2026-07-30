# Projection des facteurs et construction des distributions predictives.

stanfit_draw_matrix_base <- function(fit, variables = NULL) {
  # Extraction de secours ne faisant appel qu'a la structure serialisee du
  # stanfit. Elle permet de relire les tirages deja produits meme si RStan
  # n'est pas charge dans la session qui ne fait que construire le rapport.
  simulation <- attr(fit, "sim", exact = TRUE)
  assert_true(
    is.list(simulation) && length(simulation$samples) == simulation$chains,
    "Structure serialisee du stanfit invalide."
  )
  available <- simulation$fnames_oi
  if (!is.null(variables)) {
    keep <- vapply(available, function(column) {
      any(column == variables | startsWith(column, paste0(variables, "[")))
    }, logical(1))
    available <- available[keep]
  }
  assert_true(length(available) > 0L, "Aucun tirage Stan selectionne.")

  chains <- lapply(seq_len(simulation$chains), function(chain) {
    samples <- simulation$samples[[chain]]
    start <- simulation$warmup2[[chain]] + 1L
    stop <- simulation$n_save[[chain]]
    values <- vapply(
      available,
      function(parameter) as.numeric(samples[[parameter]][start:stop]),
      numeric(stop - start + 1L)
    )
    if (is.null(dim(values))) {
      values <- matrix(values, ncol = 1L)
    }
    colnames(values) <- available
    values
  })
  do.call(rbind, chains)
}

fit_draw_matrix <- function(fit, variables = NULL, ndraws = NULL,
                            seed = 1L) {
  if (is_rstan_fit(fit)) {
    draws <- if (requireNamespace("rstan", quietly = TRUE)) {
      if (is.null(variables)) {
        as.matrix(fit)
      } else {
        as.matrix(fit, pars = variables)
      }
    } else {
      stanfit_draw_matrix_base(fit, variables)
    }
  } else {
    draws <- fit$draws(variables = variables, format = "matrix")
    draws <- as.matrix(draws)
  }
  if (!is.null(ndraws) && nrow(draws) > ndraws) {
    set.seed(seed)
    draws <- draws[sample.int(nrow(draws), ndraws), , drop = FALSE]
  }
  draws
}

extract_vector_parameter <- function(draws, name, size) {
  columns <- sprintf("%s[%d]", name, seq_len(size))
  missing <- setdiff(columns, colnames(draws))
  assert_true(!length(missing), paste(
    "Parametres absents :", paste(missing, collapse = ", ")
  ))
  draws[, columns, drop = FALSE]
}

extract_scalar_parameter <- function(draws, name) {
  assert_true(name %in% colnames(draws), paste("Parametre absent :", name))
  as.numeric(draws[, name])
}

project_period_factors <- function(draws, model, training, max_horizon,
                                   seed) {
  set.seed(seed)
  S <- nrow(draws)
  T <- length(training$years)

  if (model %in% c("lc", "rh", "apc")) {
    previous <- extract_scalar_parameter(draws, sprintf("kappa[%d]", T))
    drift <- extract_scalar_parameter(draws, "drift")
    sigma <- extract_scalar_parameter(draws, "sigma_kappa")
    projected <- matrix(NA_real_, S, max_horizon)
    for (h in seq_len(max_horizon)) {
      previous <- previous + drift + sigma * stats::rnorm(S)
      projected[, h] <- previous
    }
    return(list(kappa1 = projected))
  }

  previous1 <- extract_scalar_parameter(
    draws, sprintf("kappa[1,%d]", T)
  )
  previous2 <- extract_scalar_parameter(
    draws, sprintf("kappa[2,%d]", T)
  )
  drift <- extract_vector_parameter(draws, "drift", 2L)
  sigma <- extract_vector_parameter(draws, "sigma_kappa", 2L)
  rho <- extract_scalar_parameter(draws, "Omega[2,1]")
  rho <- pmax(-0.999999, pmin(0.999999, rho))

  projected1 <- matrix(NA_real_, S, max_horizon)
  projected2 <- matrix(NA_real_, S, max_horizon)
  for (h in seq_len(max_horizon)) {
    z1 <- stats::rnorm(S)
    z2 <- stats::rnorm(S)
    innovation1 <- sigma[, 1L] * z1
    innovation2 <- sigma[, 2L] * (
      rho * z1 + sqrt(1 - rho^2) * z2
    )
    previous1 <- previous1 + drift[, 1L] + innovation1
    previous2 <- previous2 + drift[, 2L] + innovation2
    projected1[, h] <- previous1
    projected2[, h] <- previous2
  }
  list(kappa1 = projected1, kappa2 = projected2)
}

project_cohort_effects <- function(draws, training, target, seed) {
  set.seed(seed)
  S <- nrow(draws)
  cohort_train <- outer(
    training$ages, training$years,
    function(age, year) year - age
  )
  cohort_min <- min(cohort_train)
  cohort_max <- max(cohort_train)
  C <- cohort_max - cohort_min + 1L
  gamma <- extract_vector_parameter(draws, "gamma", C)

  target_cohort <- as.integer(target$year - target$age)
  assert_true(
    min(target_cohort) >= cohort_min,
    "Une cohorte cible precede la plage d'apprentissage."
  )

  extra <- max(0L, max(target_cohort) - cohort_max)
  if (extra > 0L) {
    psi1 <- extract_scalar_parameter(draws, "psi1")
    psi2 <- extract_scalar_parameter(draws, "psi2")
    sigma <- extract_scalar_parameter(draws, "sigma_gamma")
    extended <- matrix(NA_real_, S, C + extra)
    extended[, seq_len(C)] <- gamma
    for (index in (C + 1L):(C + extra)) {
      extended[, index] <- psi1 * extended[, index - 1L] +
        psi2 * extended[, index - 2L] +
        sigma * stats::rnorm(S)
    }
    gamma <- extended
  }
  indices <- target_cohort - cohort_min + 1L
  gamma[, indices, drop = FALSE]
}

forecast_model <- function(fit, model, training, target, cfg,
                           keep_draws = TRUE, seed = cfg$seed) {
  assert_true(nrow(target) > 0L, "Grille de prevision vide.")
  assert_true(all(target$horizon >= 1L), "Les horizons doivent etre positifs.")
  max_horizon <- max(target$horizon)
  draws <- fit_draw_matrix(
    fit,
    ndraws = cfg$forecast_draws,
    seed = stable_seed(seed, model, "posterior_subset")
  )
  S <- nrow(draws)
  N <- nrow(target)
  age_index <- match(target$age, training$ages)
  assert_true(!anyNA(age_index), "Age cible absent de l'ajustement.")

  period <- project_period_factors(
    draws, model, training, max_horizon,
    seed = stable_seed(seed, model, "period")
  )
  cohort <- if (model %in% c("rh", "apc", "m6")) {
    project_cohort_effects(
      draws, training, target,
      seed = stable_seed(seed, model, "cohort")
    )
  } else {
    NULL
  }

  eta <- matrix(NA_real_, S, N)
  if (model %in% c("lc", "rh")) {
    alpha <- extract_vector_parameter(draws, "alpha", length(training$ages))
    beta <- extract_vector_parameter(draws, "beta", length(training$ages))
    for (i in seq_len(N)) {
      eta[, i] <- alpha[, age_index[i]] +
        beta[, age_index[i]] * period$kappa1[, target$horizon[i]]
      if (model == "rh") eta[, i] <- eta[, i] + cohort[, i]
    }
  } else if (model == "apc") {
    alpha <- extract_vector_parameter(draws, "alpha", length(training$ages))
    for (i in seq_len(N)) {
      eta[, i] <- alpha[, age_index[i]] +
        period$kappa1[, target$horizon[i]] + cohort[, i]
    }
  } else {
    centered_age <- target$age - mean(training$ages)
    for (i in seq_len(N)) {
      eta[, i] <- period$kappa1[, target$horizon[i]] +
        centered_age[i] * period$kappa2[, target$horizon[i]]
      if (model == "m6") eta[, i] <- eta[, i] + cohort[, i]
    }
  }

  force <- exp(eta)
  phi <- extract_scalar_parameter(draws, "phi")
  expected_deaths <- sweep(force, 2L, target$exposure, `*`)
  log_predictive <- rep(NA_real_, N)
  if ("deaths" %in% names(target)) {
    for (i in seq_len(N)) {
      if (is.finite(target$deaths[i])) {
        log_likelihood <- stats::dnbinom(
          target$deaths[i],
          mu = expected_deaths[, i],
          size = phi,
          log = TRUE
        )
        log_predictive[i] <- log_mean_exp(log_likelihood)
      }
    }
  }

  summary <- target
  summary$model <- model
  summary$log_predictive <- log_predictive
  summary$mean_force <- colMeans(force)
  summary$mean_deaths <- colMeans(expected_deaths)

  result <- list(summary = summary)
  if (keep_draws) {
    set.seed(stable_seed(seed, model, "observation_noise"))
    count_draws <- matrix(
      stats::rnbinom(
        S * N,
        size = rep(phi, times = N),
        mu = as.vector(expected_deaths)
      ),
      nrow = S,
      ncol = N
    )
    rate_draws <- sweep(count_draws, 2L, target$exposure, `/`)
    result$count_draws <- count_draws
    result$rate_draws <- rate_draws
    result$force_draws <- force
  }
  result
}

make_future_force_target <- function(ages, origin, horizons) {
  assert_true(length(ages) == length(horizons),
              "ages et horizons doivent avoir la meme longueur.")
  data.frame(
    year = as.integer(origin + horizons),
    age = as.integer(ages),
    deaths = NA_real_,
    exposure = 1,
    origin = as.integer(origin),
    horizon = as.integer(horizons)
  )
}
