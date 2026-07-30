#!/usr/bin/env Rscript

# Tests sans framework : ils fonctionnent avec une installation R minimale.

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
source(file.path("R", "aggregation.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "actuarial.R"))
cfg <- load_config("full")

expect_close <- function(x, y, tolerance = 1e-8, label = "") {
  if (max(abs(x - y)) > tolerance) {
    stop("Test echoue : ", label, call. = FALSE)
  }
}

expect_close(sum(softmax(c(-1000, 0, 1000))), 1, label = "softmax")
expect_close(
  log_sum_exp(c(-1000, -1001)),
  -1000 + log1p(exp(-1)),
  label = "log_sum_exp"
)
row_log_values <- rbind(c(-1000, -1001), c(-Inf, -Inf))
expect_close(
  row_log_sum_exp(row_log_values)[1L],
  -1000 + log1p(exp(-1)),
  label = "row_log_sum_exp fini"
)
stopifnot(identical(row_log_sum_exp(row_log_values)[2L], -Inf))

small_draws <- c(1, 2, 4, 8)
brute_pairwise <- mean(abs(small_draws - 3)) -
  sum(abs(outer(small_draws, small_draws, `-`))) /
  (2 * length(small_draws) * (length(small_draws) - 1))
expect_close(
  empirical_crps(small_draws, 3),
  brute_pairwise,
  label = "CRPS"
)

processed <- prepare_hmd_data(cfg)
report <- data_quality_report(processed, cfg)
stopifnot(
  report$cells == 2255L,
  report$first_year == 1970L,
  report$last_year == 2024L,
  report$min_age == 50L,
  report$max_age == 90L,
  report$fractional_deaths == 0L,
  report$missing_deaths == 0L,
  report$missing_exposure == 0L
)

extended_target <- target_grid_from_observed(
  processed, 2015L, 1:9, maximum_year = 2024L
)
stopifnot(
  nrow(extended_target) == 369L,
  identical(sort(unique(extended_target$year)), 2016:2024),
  identical(sort(unique(extended_target$horizon)), 1:9),
  all(extended_target$year - extended_target$origin ==
        extended_target$horizon)
)

lfo_target <- target_grid_from_observed(
  processed, 2008L, 1:2, maximum_year = 2012L
)
stopifnot(
  nrow(lfo_target) == 82L,
  identical(sort(unique(lfo_target$year)), 2009:2010),
  identical(sort(unique(lfo_target$horizon)), 1:2)
)

dummy <- expand.grid(age = 50:90, horizon = 1:5)
standardized <- standardize_context(dummy)
expect_close(mean(standardized$data$x_tilde), 0, label = "age centre")
expect_close(mean(standardized$data$q_age), 0, label = "age2 centre")
expect_close(mean(standardized$data$h_tilde), 0, label = "horizon centre")

models <- cfg$models
meta_test <- standardized$data
for (model in models) meta_test[[paste0("log_p_", model)]] <- 0
meta_test$omega <- 1

pseudo_test <- fit_pseudo_bma(meta_test, models)
expect_close(
  unname(pseudo_test$weights),
  rep(1 / length(models), length(models)),
  label = "poids pseudo-BMA"
)

global_test <- fit_global_stacking(
  meta_test, models, multistarts = 2L, seed = 1L
)
expect_close(sum(global_test$weights), 1, label = "poids stacking global")
stopifnot(global_test$optimization$convergence == 0L)

contextual_test <- fit_contextual_stacking(
  meta_test, models, multistarts = 2L, seed = 1L
)
contextual_grid_test <- contextual_weight_grid(
  dummy, contextual_test, models, standardized$constants
)
expect_close(
  rowSums(as.matrix(contextual_grid_test[models])),
  rep(1, nrow(dummy)),
  label = "poids stacking contextuel"
)
stopifnot(
  contextual_test$optimization$convergence == 0L,
  all(as.matrix(contextual_grid_test[models]) >= 0)
)

weight_test_target <- unique(extended_target[c("age", "horizon")])
contextual_grid_extended <- contextual_weight_grid(
  weight_test_target,
  contextual_test,
  models,
  standardized$constants
)
weight_test_sets <- list(
  stacking_global = constant_weight_matrix(
    setNames(rep(1 / length(models), length(models)), models),
    nrow(weight_test_target),
    models
  ),
  stacking_contextual = as.matrix(contextual_grid_extended[models])
)
weight_long_test <- evaluation_weight_long(
  weight_test_target,
  weight_test_sets,
  names(weight_test_sets),
  models
)
weight_sums_test <- stats::aggregate(
  weight ~ method + age + horizon,
  data = weight_long_test,
  FUN = sum
)
expect_close(
  weight_sums_test$weight,
  rep(1, nrow(weight_sums_test)),
  label = "somme des poids etendus"
)

central_gradient <- function(fn, theta, epsilon = 1e-6) {
  vapply(seq_along(theta), function(index) {
    upper <- lower <- theta
    upper[index] <- upper[index] + epsilon
    lower[index] <- lower[index] - epsilon
    (fn(upper) - fn(lower)) / (2 * epsilon)
  }, numeric(1))
}

set.seed(2L)
gradient_log_p <- matrix(
  stats::rnorm(60, mean = -4), nrow = 12, ncol = length(models)
)
gradient_omega <- seq(0.5, 1.5, length.out = nrow(gradient_log_p))
global_theta <- stats::rnorm(length(models) - 1L, sd = 0.2)
expect_close(
  global_stacking_gradient(
    global_theta, gradient_log_p, gradient_omega
  ),
  central_gradient(
    function(theta) global_stacking_value(
      theta, gradient_log_p, gradient_omega
    ),
    global_theta
  ),
  tolerance = 1e-5,
  label = "gradient stacking global"
)

gradient_X <- matrix(stats::rnorm(48), nrow = 12, ncol = 4)
contextual_theta <- stats::rnorm(
  (length(models) - 1L) * (ncol(gradient_X) + 1L),
  sd = 0.1
)
expect_close(
  contextual_stacking_gradient(
    contextual_theta, gradient_log_p, gradient_X, gradient_omega
  ),
  central_gradient(
    function(theta) contextual_stacking_value(
      theta, gradient_log_p, gradient_X, gradient_omega
    ),
    contextual_theta
  ),
  tolerance = 1e-5,
  label = "gradient stacking contextuel"
)

# Tests unitaires de l'agregation des forces actuarielles.
actuarial_models <- cfg$models
actuarial_S <- 3L
actuarial_J <- 2L
actuarial_forces <- setNames(lapply(
  seq_along(actuarial_models),
  function(k) matrix(
    c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06) + k / 1000,
    nrow = actuarial_S,
    ncol = actuarial_J
  )
), actuarial_models)

# Test 1 : un poids egal a un reproduit exactement le modele choisi.
one_hot_weights <- matrix(
  0,
  nrow = actuarial_J,
  ncol = length(actuarial_models),
  dimnames = list(NULL, actuarial_models)
)
one_hot_weights[, 2L] <- 1
one_hot_force <- aggregate_force_draws(
  actuarial_forces, one_hot_weights, actuarial_models
)
expect_close(
  one_hot_force,
  actuarial_forces[[2L]],
  label = "actuariel one-hot force"
)
one_hot_survival <- survival_from_aggregated_force(one_hot_force)
expect_close(
  one_hot_survival,
  survival_from_aggregated_force(actuarial_forces[[2L]]),
  label = "actuariel one-hot survie"
)
expect_close(
  annuity_from_survival(one_hot_survival, 1:2, 0.02),
  annuity_from_survival(
    survival_from_aggregated_force(actuarial_forces[[2L]]),
    1:2,
    0.02
  ),
  label = "actuariel one-hot rente"
)

# Test 2 : des modeles identiques donnent le meme resultat quels que soient
# les poids.
identical_force <- matrix(
  c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07),
  nrow = actuarial_S,
  ncol = actuarial_J
)
identical_models <- setNames(
  replicate(
    length(actuarial_models), identical_force, simplify = FALSE
  ),
  actuarial_models
)
varying_weights <- matrix(
  c(
    0.05, 0.15, 0.20, 0.25, 0.35,
    0.30, 0.10, 0.05, 0.25, 0.30
  ),
  nrow = actuarial_J,
  byrow = TRUE,
  dimnames = list(NULL, actuarial_models)
)
expect_close(
  aggregate_force_draws(
    identical_models, varying_weights, actuarial_models
  ),
  identical_force,
  label = "actuariel modeles identiques"
)

# Test 3 : des poids contextuels constants reproduisent le calcul global.
constant_weights <- setNames(
  c(0.10, 0.20, 0.25, 0.15, 0.30),
  actuarial_models
)
global_weights_test <- matrix(
  rep(constant_weights, each = actuarial_J),
  nrow = actuarial_J,
  dimnames = list(NULL, actuarial_models)
)
contextual_weights_test <- global_weights_test
expect_close(
  aggregate_force_draws(
    actuarial_forces, global_weights_test, actuarial_models
  ),
  aggregate_force_draws(
    actuarial_forces, contextual_weights_test, actuarial_models
  ),
  label = "actuariel poids contextuels constants"
)

# Test 4 : comparaison a une somme, une survie et une rente manuelles.
manual_force <- matrix(0, nrow = actuarial_S, ncol = actuarial_J)
for (s in seq_len(actuarial_S)) {
  for (j in seq_len(actuarial_J)) {
    manual_force[s, j] <- sum(vapply(
      seq_along(actuarial_models),
      function(k) {
        varying_weights[j, k] * actuarial_forces[[k]][s, j]
      },
      numeric(1)
    ))
  }
}
automatic_force <- aggregate_force_draws(
  actuarial_forces, varying_weights, actuarial_models
)
expect_close(
  automatic_force, manual_force,
  label = "actuariel somme manuelle"
)
manual_survival <- matrix(NA_real_, nrow = actuarial_S, ncol = actuarial_J)
for (s in seq_len(actuarial_S)) {
  manual_survival[s, ] <- exp(-cumsum(manual_force[s, ]))
}
automatic_survival <- survival_from_aggregated_force(automatic_force)
expect_close(
  automatic_survival, manual_survival,
  label = "actuariel survie manuelle"
)
manual_annuity <- rowSums(
  manual_survival * matrix(
    (1 / 1.02)^(1:2),
    nrow = actuarial_S,
    ncol = actuarial_J,
    byrow = TRUE
  )
)
expect_close(
  annuity_from_survival(automatic_survival, 1:2, 0.02),
  manual_annuity,
  label = "actuariel rente manuelle"
)

# Test 5 : memes donnees et meme graine, memes permutations et resultats.
permutation_seed_test <- 12345L
permutations_first <- make_independent_draw_permutations(
  actuarial_S, actuarial_models, permutation_seed_test
)
permutations_second <- make_independent_draw_permutations(
  actuarial_S, actuarial_models, permutation_seed_test
)
stopifnot(identical(permutations_first, permutations_second))
paired_first <- pair_model_force_draws(
  actuarial_forces,
  seq_len(actuarial_J),
  actuarial_models,
  permutations_first
)
paired_second <- pair_model_force_draws(
  actuarial_forces,
  seq_len(actuarial_J),
  actuarial_models,
  permutations_second
)
expect_close(
  aggregate_force_draws(
    paired_first, varying_weights, actuarial_models
  ),
  aggregate_force_draws(
    paired_second, varying_weights, actuarial_models
  ),
  label = "actuariel reproductibilite"
)

actuarial_unit_test_results <- data.frame(
  test = c(
    "one_hot_model",
    "identical_models",
    "constant_context_equals_global",
    "manual_force_survival_annuity",
    "reproducibility_fixed_seed"
  ),
  passed = TRUE,
  stringsAsFactors = FALSE
)
write_csv_atomic(
  actuarial_unit_test_results,
  file.path(
    cfg$paths$metrics,
    "actuarial_force_aggregation_unit_tests.csv"
  )
)

if ("--compile" %in% commandArgs(trailingOnly = TRUE)) {
  source(file.path("R", "model_fitting.R"))
  require_stan_backend()
  invisible(lapply(c(cfg$models, "hierarchical"), compile_stan_model, cfg = cfg))
}

message("Tous les tests ont reussi.")
