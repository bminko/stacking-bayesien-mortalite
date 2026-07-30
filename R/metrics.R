# Regles de score, intervalles predictifs et melanges de tirages.

empirical_crps <- function(draws, observation) {
  draws <- sort(as.numeric(draws))
  S <- length(draws)
  assert_true(S >= 2L, "Le CRPS requiert au moins deux tirages.")
  first_term <- mean(abs(draws - observation))
  coefficients <- 2 * seq_len(S) - S - 1
  second_term <- sum(coefficients * draws) / (S * (S - 1))
  first_term - second_term
}

mixture_draws_for_cell <- function(predictions, models, cell, weights,
                                   type = c("count", "rate"), seed = 1L) {
  type <- match.arg(type)
  field <- if (type == "count") "count_draws" else "rate_draws"
  S <- nrow(predictions[[models[[1L]]]][[field]])
  set.seed(seed)
  selected_model <- sample.int(
    length(models), S, replace = TRUE, prob = weights
  )
  result <- numeric(S)
  for (k in seq_along(models)) {
    rows <- which(selected_model == k)
    if (length(rows)) {
      result[rows] <- predictions[[models[k]]][[field]][rows, cell]
    }
  }
  result
}

constant_weight_matrix <- function(weights, N, models) {
  weights <- weights[models]
  matrix(
    rep(weights, each = N),
    nrow = N,
    ncol = length(models),
    dimnames = list(NULL, models)
  )
}

build_evaluation_weight_sets <- function(target, aggregation, models) {
  N <- nrow(target)
  result <- list()
  for (model in models) {
    weight <- setNames(rep(0, length(models)), models)
    weight[model] <- 1
    result[[model]] <- constant_weight_matrix(weight, N, models)
  }

  if (isTRUE(aggregation$bma$available)) {
    result$bma <- constant_weight_matrix(
      aggregation$bma$weights, N, models
    )
  }
  result$pseudo_bma <- constant_weight_matrix(
    aggregation$pseudo_bma$weights, N, models
  )
  result$stacking_global <- constant_weight_matrix(
    aggregation$global$weights, N, models
  )

  nonregularized_grid <- contextual_weight_grid(
    target[c("age", "horizon")],
    aggregation$contextual,
    models,
    aggregation$standardization
  )
  result$stacking_contextual <- as.matrix(nonregularized_grid[models])

  result$stacking_hierarchical <- weight_matrix_from_long(
    aggregation$hierarchical_weights,
    target[c("age", "horizon")],
    models
  )
  result
}

evaluate_predictions <- function(predictions, aggregation, models, cfg) {
  reference <- predictions[[models[[1L]]]]$summary
  target <- reference[c(
    "year", "age", "deaths", "exposure", "origin", "horizon"
  )]
  N <- nrow(target)
  for (model in models[-1L]) {
    candidate <- predictions[[model]]$summary
    assert_true(
      identical(target$year, candidate$year) &&
        identical(target$age, candidate$age),
      paste("Grille predictive non alignee pour", model)
    )
  }

  log_p <- sapply(models, function(model) {
    predictions[[model]]$summary$log_predictive
  })
  colnames(log_p) <- models
  weight_sets <- build_evaluation_weight_sets(target, aggregation, models)

  output <- vector("list", length(weight_sets) * N)
  position <- 1L
  for (method in names(weight_sets)) {
    weights <- weight_sets[[method]]
    for (i in seq_len(N)) {
      log_mixture <- log_sum_exp(log(weights[i, ]) + log_p[i, ])
      count_draws <- mixture_draws_for_cell(
        predictions, models, i, weights[i, ], "count",
        seed = stable_seed(cfg$seed, cfg$profile, method, i, "count")
      )
      rate_draws <- mixture_draws_for_cell(
        predictions, models, i, weights[i, ], "rate",
        seed = stable_seed(cfg$seed, cfg$profile, method, i, "rate")
      )
      interval80 <- stats::quantile(
        rate_draws, c(0.10, 0.90), names = FALSE
      )
      interval95 <- stats::quantile(
        rate_draws, c(0.025, 0.975), names = FALSE
      )
      observed_rate <- target$deaths[i] / target$exposure[i]

      output[[position]] <- data.frame(
        method = method,
        year = target$year[i],
        age = target$age[i],
        horizon = target$horizon[i],
        deaths = target$deaths[i],
        exposure = target$exposure[i],
        observed_rate = observed_rate,
        logs = -log_mixture,
        crps = empirical_crps(rate_draws, observed_rate),
        absolute_error_deaths = abs(target$deaths[i] - mean(count_draws)),
        coverage80 = as.integer(
          observed_rate >= interval80[1L] &&
            observed_rate <= interval80[2L]
        ),
        width80 = interval80[2L] - interval80[1L],
        coverage95 = as.integer(
          observed_rate >= interval95[1L] &&
            observed_rate <= interval95[2L]
        ),
        width95 = interval95[2L] - interval95[1L],
        stringsAsFactors = FALSE
      )
      position <- position + 1L
    }
  }
  cells <- do.call(rbind, output)
  rownames(cells) <- NULL

  summaries <- summarize_evaluation_cells(cells)

  list(
    cells = cells,
    overall = summaries$overall,
    by_horizon = summaries$by_horizon,
    by_age = summaries$by_age,
    by_age_group = summaries$by_age_group,
    weight_sets = weight_sets
  )
}

evaluation_metric_names <- function() {
  c(
    "logs", "crps", "absolute_error_deaths",
    "coverage80", "width80", "coverage95", "width95"
  )
}

add_age_group <- function(cells) {
  cells$age_group <- cut(
    cells$age,
    breaks = c(49, 59, 69, 79, 90),
    labels = c("50-59", "60-69", "70-79", "80-90"),
    include.lowest = TRUE,
    right = TRUE
  )
  cells
}

summarize_evaluation_cells <- function(cells) {
  assert_true(nrow(cells) > 0L, "Aucune cellule a resumer.")
  metrics <- evaluation_metric_names()
  missing <- setdiff(c("method", "age", "horizon", metrics), names(cells))
  assert_true(
    !length(missing),
    paste("Colonnes de metriques absentes :", paste(missing, collapse = ", "))
  )
  cells <- add_age_group(cells)
  overall <- group_mean(cells, metrics, "method")
  by_horizon <- group_mean(cells, metrics, c("method", "horizon"))
  by_age <- group_mean(cells, metrics, c("method", "age"))
  by_age_group <- group_mean(cells, metrics, c("method", "age_group"))
  for (object in c("overall", "by_horizon", "by_age", "by_age_group")) {
    value <- get(object)
    names(value)[names(value) == "absolute_error_deaths"] <- "mae_deaths"
    assign(object, value)
  }
  list(
    overall = overall,
    by_horizon = by_horizon,
    by_age = by_age,
    by_age_group = by_age_group
  )
}

summarize_evaluation_periods <- function(cells, periods) {
  rows <- list(
    overall = list(),
    by_horizon = list(),
    by_age = list(),
    by_age_group = list()
  )
  for (label in names(periods)) {
    part <- cells[cells$year %in% periods[[label]], , drop = FALSE]
    assert_true(
      nrow(part) > 0L,
      paste("Aucune observation pour la periode", label)
    )
    summary <- summarize_evaluation_cells(part)
    for (level in names(rows)) {
      summary[[level]]$period <- label
      rows[[level]][[label]] <- summary[[level]]
    }
  }
  lapply(rows, function(parts) {
    result <- do.call(rbind, parts)
    rownames(result) <- NULL
    result
  })
}

evaluation_weight_long <- function(target, weight_sets, methods, models) {
  rows <- vector("list", length(methods))
  for (index in seq_along(methods)) {
    method <- methods[[index]]
    weights <- weight_sets[[method]]
    assert_true(!is.null(weights), paste("Poids absents pour", method))
    assert_true(
      nrow(weights) == nrow(target) && ncol(weights) == length(models),
      paste("Dimensions de poids invalides pour", method)
    )
    rows[[index]] <- data.frame(
      method = rep(method, each = nrow(target) * length(models)),
      age = rep(target$age, times = length(models)),
      horizon = rep(target$horizon, times = length(models)),
      model = rep(models, each = nrow(target)),
      weight = as.vector(weights),
      stringsAsFactors = FALSE
    )
  }
  result <- unique(do.call(rbind, rows))
  result <- result[order(
    result$method, result$horizon, result$age, result$model
  ), ]
  rownames(result) <- NULL
  result
}

summarize_weights_by_horizon <- function(weight_long) {
  mean_weight <- stats::aggregate(
    weight ~ method + horizon + model,
    data = weight_long,
    FUN = mean
  )
  min_weight <- stats::aggregate(
    weight ~ method + horizon + model,
    data = weight_long,
    FUN = min
  )
  max_weight <- stats::aggregate(
    weight ~ method + horizon + model,
    data = weight_long,
    FUN = max
  )
  names(mean_weight)[4L] <- "weight_mean_over_ages"
  names(min_weight)[4L] <- "weight_min_over_ages"
  names(max_weight)[4L] <- "weight_max_over_ages"
  result <- Reduce(
    function(left, right) {
      merge(left, right, by = c("method", "horizon", "model"))
    },
    list(mean_weight, min_weight, max_weight)
  )
  result[order(result$method, result$horizon, result$model), ]
}

summarize_distribution <- function(x, prefix = "") {
  quantiles <- stats::quantile(x, c(0.025, 0.05, 0.5, 0.95, 0.975))
  result <- c(
    mean = mean(x),
    sd = stats::sd(x),
    q025 = quantiles[[1L]],
    q05 = quantiles[[2L]],
    median = quantiles[[3L]],
    q95 = quantiles[[4L]],
    q975 = quantiles[[5L]]
  )
  names(result) <- paste0(prefix, names(result))
  result
}
