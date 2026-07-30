#!/usr/bin/env Rscript

# Stabilite Monte Carlo des predictions finales.
#
# Comparaison demandee : 1 000 contre 4 000 tirages, trois graines, sans
# aucun reajustement Stan. Une sensibilite distincte propage les tirages des
# poids hierarchiques dans la construction du melange predictif.

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
source(file.path("R", "model_fitting.R"))
source(file.path("R", "forecasting.R"))
source(file.path("R", "aggregation.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "reporting.R"))
cfg <- load_config()

results_dir <- file.path(cfg$paths$results, "chapter3_results")
figures_dir <- file.path(cfg$root, "output", "pdf", "figures")
processed <- readRDS(file.path(
  cfg$paths$processed, "mortality_data.rds"
))
aggregation <- readRDS(file.path(
  cfg$paths$weights, "aggregation_results.rds"
))
training <- subset_training_data(processed, cfg$validation_end)
target <- target_grid_from_observed(
  processed, cfg$validation_end,
  1:(cfg$test_end - cfg$validation_end),
  cfg$test_end
)
methods <- c(
  "stacking_global", "stacking_contextual",
  "stacking_hierarchical"
)
settings <- expand.grid(
  forecast_draws = c(1000L, 4000L),
  seed_offset = 0:2
)
settings$seed <- cfg$seed + settings$seed_offset

hierarchical_weight_draw_array <- function(target, S, seed) {
  standardized <- standardize_context(
    target, aggregation$standardization
  )$data
  X <- context_design_matrix(standardized)
  fit <- readRDS(aggregation$hierarchical_fit_path)
  draws <- fit_draw_matrix(
    fit, variables = c("alpha", "beta"),
    ndraws = S, seed = seed
  )
  K <- length(cfg$models)
  alpha <- extract_vector_parameter(draws, "alpha", K - 1L)
  beta <- array(NA_real_, dim = c(S, K - 1L, ncol(X)))
  for (k in seq_len(K - 1L)) {
    for (p in seq_len(ncol(X))) {
      beta[, k, p] <- extract_scalar_parameter(
        draws, sprintf("beta[%d,%d]", k, p)
      )
    }
  }
  result <- array(
    NA_real_, c(S, nrow(target), K),
    dimnames = list(NULL, NULL, cfg$models)
  )
  for (cell in seq_len(nrow(target))) {
    score <- alpha
    for (p in seq_len(ncol(X))) {
      score <- score + beta[, , p] * X[cell, p]
    }
    result[, cell, ] <- scores_to_weights(cbind(score, 0))
  }
  result
}

draw_mixture <- function(predictions, cell, weights, field, seed) {
  S <- nrow(predictions[[cfg$models[[1L]]]][[field]])
  set.seed(seed)
  uniform <- stats::runif(S)
  if (is.matrix(weights)) {
    cumulative <- t(apply(weights, 1L, cumsum))
    selected <- 1L + rowSums(
      matrix(uniform, S, length(cfg$models) - 1L) >
        cumulative[, -length(cfg$models), drop = FALSE]
    )
  } else {
    selected <- 1L + rowSums(outer(
      uniform, cumsum(weights)[-length(cfg$models)], `>`
    ))
  }
  result <- numeric(S)
  for (model_index in seq_along(cfg$models)) {
    rows <- which(selected == model_index)
    if (length(rows)) {
      result[rows] <- predictions[[cfg$models[[model_index]]]][[field]][
        rows, cell
      ]
    }
  }
  result
}

evaluate_weight_mode <- function(
    predictions, weights, method, weight_mode, seed) {
  logp <- sapply(cfg$models, function(model) {
    predictions[[model]]$summary$log_predictive
  })
  cells <- vector("list", nrow(target))
  for (cell in seq_len(nrow(target))) {
    cell_weights <- if (length(dim(weights)) == 3L) {
      weights[, cell, , drop = FALSE][, 1L, ]
    } else {
      weights[cell, ]
    }
    mean_weights <- if (is.matrix(cell_weights)) {
      colMeans(cell_weights)
    } else {
      cell_weights
    }
    count <- draw_mixture(
      predictions, cell, cell_weights, "count_draws",
      stable_seed(seed, method, weight_mode, cell, "count")
    )
    rate <- draw_mixture(
      predictions, cell, cell_weights, "rate_draws",
      stable_seed(seed, method, weight_mode, cell, "rate")
    )
    interval80 <- stats::quantile(rate, c(0.10, 0.90), names = FALSE)
    interval95 <- stats::quantile(rate, c(0.025, 0.975), names = FALSE)
    observed_rate <- target$deaths[cell] / target$exposure[cell]
    cells[[cell]] <- data.frame(
      method = method,
      weight_mode = weight_mode,
      year = target$year[cell],
      age = target$age[cell],
      horizon = target$horizon[cell],
      deaths = target$deaths[cell],
      exposure = target$exposure[cell],
      observed_rate = observed_rate,
      logs = -log_sum_exp(log(mean_weights) + logp[cell, ]),
      crps = empirical_crps(rate, observed_rate),
      absolute_error_deaths = abs(target$deaths[cell] - mean(count)),
      coverage80 = as.integer(
        observed_rate >= interval80[1L] &
          observed_rate <= interval80[2L]
      ),
      width80 = diff(interval80),
      coverage95 = as.integer(
        observed_rate >= interval95[1L] &
          observed_rate <= interval95[2L]
      ),
      width95 = diff(interval95),
      stringsAsFactors = FALSE
    )
  }
  cells <- do.call(rbind, cells)
  overall <- summarize_evaluation_cells(cells)$overall
  overall$weight_mode <- weight_mode
  list(cells = cells, overall = overall)
}

summary_rows <- list()
propagation_rows <- list()
position <- 1L
propagation_position <- 1L
for (setting_index in seq_len(nrow(settings))) {
  draw_count <- settings$forecast_draws[setting_index]
  seed <- settings$seed[setting_index]
  message_step(
    "Stabilite predictive : ", draw_count,
    " tirages, graine ", seed
  )
  local_cfg <- cfg
  local_cfg$forecast_draws <- draw_count
  predictions <- setNames(
    vector("list", length(cfg$models)), cfg$models
  )
  for (model in cfg$models) {
    fit <- read_model_fit(
      model, training, cfg,
      paste0("test_training_", cfg$validation_end)
    )
    predictions[[model]] <- forecast_model(
      fit, model, training, target, local_cfg,
      keep_draws = TRUE,
      seed = stable_seed(seed, "mc_stability", model)
    )
    rm(fit)
    gc(verbose = FALSE)
  }
  mean_weight_sets <- build_evaluation_weight_sets(
    target, aggregation, cfg$models
  )[methods]
  for (method in methods) {
    evaluation_result <- evaluate_weight_mode(
      predictions, mean_weight_sets[[method]], method,
      "posterior_mean", seed
    )
    row <- evaluation_result$overall
    row$forecast_draws <- draw_count
    row$seed <- seed
    summary_rows[[position]] <- row
    position <- position + 1L
  }

  # La propagation integrale des poids est comparee aux poids moyens pour
  # chacune des six configurations; elle ne change aucun log-predictif modele.
  propagated <- hierarchical_weight_draw_array(
    target, draw_count,
    stable_seed(seed, "mc_stability", "weight_draws")
  )
  propagated_result <- evaluate_weight_mode(
    predictions, propagated, "stacking_hierarchical",
    "posterior_draw", seed
  )
  propagated_row <- propagated_result$overall
  propagated_row$forecast_draws <- draw_count
  propagated_row$seed <- seed
  propagation_rows[[propagation_position]] <- propagated_row
  propagation_position <- propagation_position + 1L
  rm(predictions, propagated)
  gc(verbose = FALSE)
}

stability <- do.call(rbind, summary_rows)
propagation <- do.call(rbind, propagation_rows)
all_results <- rbind(stability, propagation)
rownames(all_results) <- NULL
write_csv_atomic(
  all_results,
  file.path(results_dir, "49_prediction_mc_stability.csv")
)

metric_names <- c(
  "logs", "crps", "mae_deaths",
  "coverage80", "width80", "coverage95", "width95"
)
stability_summary <- do.call(rbind, lapply(
  split(
    all_results,
    interaction(
      all_results$method, all_results$weight_mode,
      all_results$forecast_draws, drop = TRUE
    )
  ),
  function(part) {
    row <- part[1L, c("method", "weight_mode", "forecast_draws")]
    for (metric in metric_names) {
      row[[paste0(metric, "_mean_over_seeds")]] <- mean(part[[metric]])
      row[[paste0(metric, "_sd_over_seeds")]] <- stats::sd(part[[metric]])
      row[[paste0(metric, "_range_over_seeds")]] <-
        diff(range(part[[metric]]))
    }
    row
  }
))
rownames(stability_summary) <- NULL
write_csv_atomic(
  stability_summary,
  file.path(results_dir, "50_prediction_mc_stability_summary.csv")
)

render_chapter_figure(
  "prediction_mc_stability",
  function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.2, 1.0))
    colors <- chapter_model_colors()
    methods <- c(
      "stacking_global", "stacking_contextual",
      "stacking_hierarchical"
    )
    part <- stability
    graphics::plot(
      range(part$forecast_draws), range(part$logs),
      type = "n", xlab = "Nombre de tirages predictifs",
      ylab = "LogS moyen",
      main = "Stabilite du LogS"
    )
    for (method in methods) {
      piece <- part[part$method == method, ]
      graphics::points(
        jitter(piece$forecast_draws, factor = 0.4),
        piece$logs, col = colors[[method]], pch = 16
      )
      means <- stats::aggregate(
        logs ~ forecast_draws, data = piece, FUN = mean
      )
      graphics::lines(
        means$forecast_draws, means$logs,
        col = colors[[method]], lwd = 2
      )
    }
    graphics::plot(
      range(part$forecast_draws), range(part$crps),
      type = "n", xlab = "Nombre de tirages predictifs",
      ylab = "CRPS moyen",
      main = "Stabilite du CRPS"
    )
    for (method in methods) {
      piece <- part[part$method == method, ]
      graphics::points(
        jitter(piece$forecast_draws, factor = 0.4),
        piece$crps, col = colors[[method]], pch = 16
      )
      means <- stats::aggregate(
        crps ~ forecast_draws, data = piece, FUN = mean
      )
      graphics::lines(
        means$forecast_draws, means$crps,
        col = colors[[method]], lwd = 2
      )
    }
    graphics::legend(
      "topright",
      legend = chapter_method_labels()[methods],
      col = colors[methods], pch = 16, lwd = 2,
      bty = "n", cex = 0.75
    )
  },
  figures_dir, width = 10, height = 5
)

render_chapter_figure(
  "hierarchical_weight_propagation",
  function() {
    graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.2, 1.0))
    part <- all_results[
      all_results$method == "stacking_hierarchical", ]
    group <- interaction(
      part$forecast_draws, part$weight_mode, drop = TRUE
    )
    labels <- levels(group)
    values_logs <- split(part$logs, group)
    values_crps <- split(part$crps, group)
    graphics::boxplot(
      values_logs, names = labels, las = 2,
      col = rep(c("#66C2A5", "#FC8D62"), 2),
      ylab = "LogS moyen",
      main = "Poids moyens ou propages"
    )
    graphics::boxplot(
      values_crps, names = labels, las = 2,
      col = rep(c("#66C2A5", "#FC8D62"), 2),
      ylab = "CRPS moyen",
      main = "Sensibilite predictive"
    )
  },
  figures_dir, width = 10, height = 5
)

message("Analyse de stabilite Monte Carlo terminee.")
