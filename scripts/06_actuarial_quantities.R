#!/usr/bin/env Rscript

# Etape 6 corrigee : survie et rentes temporaires depuis l'origine 2015.
#
# Pour les quantites actuarielles, les forces de mortalite des modeles sont
# d'abord combinees par une moyenne ponderee. La survie et la rente sont
# ensuite calculees sur cette trajectoire de force agregee :
#
#   mu_agr(s, j) = somme_k w_k(j) * mu_k(s, j)
#   p_x(s, t)    = exp(-somme_{j=1}^t mu_agr(s, j))
#   a_x:n(s)     = somme_{t=1}^n v^t * p_x(s, t)
#
# Les calculs predictifs marginaux (LogS, CRPS, couverture) ne sont pas
# modifies. Les anciens fichiers actuariels sont conserves et servent
# uniquement a la comparaison diagnostique.

source(file.path("R", "utils.R"))
source(file.path("R", "data.R"))
source(file.path("R", "model_fitting.R"))
source(file.path("R", "forecasting.R"))
source(file.path("R", "aggregation.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "actuarial.R"))
source(file.path("R", "reporting.R"))
cfg <- load_config()

metrics_dir <- cfg$paths$metrics
figures_dir <- cfg$paths$figures

# Les anciens resultats, lorsqu'ils existent, servent uniquement au diagnostic
# comparatif. Ils ne sont pas requis pour une reproduction depuis zero.
old_survival_path <- file.path(
  metrics_dir, "actuarial_survival_all_rules.csv"
)
old_annuity_path <- file.path(
  metrics_dir, "actuarial_annuity_all_rules.csv"
)
has_old_actuarial_diagnostic <-
  file.exists(old_survival_path) && file.exists(old_annuity_path)
old_survival <- if (has_old_actuarial_diagnostic) {
  utils::read.csv(
    old_survival_path, stringsAsFactors = FALSE, check.names = FALSE
  )
} else {
  NULL
}
old_annuity <- if (has_old_actuarial_diagnostic) {
  utils::read.csv(
    old_annuity_path, stringsAsFactors = FALSE, check.names = FALSE
  )
} else {
  NULL
}
if (!has_old_actuarial_diagnostic) {
  message(
    "Diagnostic ancien/nouveau omis : aucun resultat actuariel categoriel ",
    "anterieur n'est disponible."
  )
}

processed <- readRDS(file.path(
  cfg$paths$processed, "mortality_data.rds"
))
aggregation <- readRDS(file.path(
  cfg$paths$weights, "aggregation_results.rds"
))
training <- subset_training_data(processed, cfg$validation_end)

initial_ages <- c(50L, 60L, 65L, 70L, 75L, 80L)
initial_ages <- initial_ages[
  initial_ages >= cfg$age_min & initial_ages < cfg$age_max
]
target_rows <- lapply(initial_ages, function(initial_age) {
  duration <- min(cfg$annuity_duration, cfg$age_max - initial_age)
  horizons <- seq_len(duration)
  result <- make_future_force_target(
    initial_age + horizons - 1L,
    cfg$validation_end,
    horizons
  )
  result$initial_age <- initial_age
  result
})
target <- do.call(rbind, target_rows)
rownames(target) <- NULL

assert_true(
  all(target$age == target$initial_age + target$horizon - 1L),
  "L'age atteint ne correspond pas a x+j-1."
)

# Les ajustements existants sont relus. Seules les projections deterministes
# a graine fixee sont reconstruites ; aucun MCMC ou LFO n'est relance.
predictions <- setNames(vector("list", length(cfg$models)), cfg$models)
context <- paste0("test_training_", cfg$validation_end)
for (model in cfg$models) {
  fit <- read_model_fit(model, training, cfg, context)
  predictions[[model]] <- forecast_model(
    fit, model, training, target, cfg,
    keep_draws = TRUE,
    seed = stable_seed(
      cfg$seed, cfg$profile, "actuarial_origin_2015", model
    )
  )
}

model_force_draws <- setNames(
  lapply(predictions, function(prediction) prediction$force_draws),
  cfg$models
)
draw_dimensions <- lapply(model_force_draws, dim)
assert_true(
  all(vapply(
    draw_dimensions,
    function(candidate) identical(candidate, draw_dimensions[[1L]]),
    logical(1)
  )),
  "Les dimensions des tirages de force different entre les modeles."
)
S <- draw_dimensions[[1L]][1L]
K <- length(cfg$models)
discount_rates <- c(0.01, 0.02, 0.03)

# Une permutation independante est appliquee a chaque modele. Les memes
# permutations sont ensuite reutilisees pour toutes les methodes afin de
# conserver l'appariement des comparaisons de rente.
permutation_seed <- stable_seed(
  cfg$seed, cfg$profile, "actuarial_independent_model_pairing"
)
permutations <- make_independent_draw_permutations(
  S, cfg$models, permutation_seed
)

context_target_for_rule <- function(target, rule) {
  result <- target
  if (rule == "freeze_h10") {
    result$horizon <- pmin(result$horizon, max(cfg$lfo_horizons))
  }
  result
}

mean_weight_sets_for_rule <- function(target, rule) {
  context_target <- context_target_for_rule(target, rule)
  result <- build_evaluation_weight_sets(
    context_target, aggregation, cfg$models
  )[c(
    "stacking_global",
    "stacking_contextual",
    "stacking_hierarchical"
  )]
  if (rule == "global_after_h10") {
    beyond <- target$horizon > max(cfg$lfo_horizons)
    global <- constant_weight_matrix(
      aggregation$global$weights, sum(beyond), cfg$models
    )
    result$stacking_contextual[beyond, ] <- global
    result$stacking_hierarchical[beyond, ] <- global
  }
  result
}

hierarchical_weight_draws_for_rule <- function(target, rule) {
  context_target <- context_target_for_rule(target, rule)
  standardized <- standardize_context(
    context_target, aggregation$standardization
  )$data
  X <- context_design_matrix(standardized)
  fit <- readRDS(aggregation$hierarchical_fit_path)
  weight_seed <- stable_seed(
    cfg$seed, cfg$profile, "actuarial_hierarchical_weight_surface"
  )
  draws <- fit_draw_matrix(
    fit,
    variables = c("alpha", "beta"),
    ndraws = S,
    seed = weight_seed
  )
  assert_true(
    nrow(draws) == S,
    "Le nombre de tirages hierarchiques ne correspond pas aux forces."
  )
  alpha <- extract_vector_parameter(draws, "alpha", K - 1L)
  beta <- array(NA_real_, dim = c(S, K - 1L, ncol(X)))
  for (k in seq_len(K - 1L)) {
    for (p in seq_len(ncol(X))) {
      beta[, k, p] <- extract_scalar_parameter(
        draws, sprintf("beta[%d,%d]", k, p)
      )
    }
  }

  # Pour une ligne s, les memes coefficients alpha[s,] et beta[s,,] sont
  # utilises sur toutes les cellules : une surface complete par trajectoire.
  result <- array(
    NA_real_,
    dim = c(S, nrow(target), K),
    dimnames = list(NULL, NULL, cfg$models)
  )
  for (cell in seq_len(nrow(target))) {
    scores <- alpha
    for (p in seq_len(ncol(X))) {
      scores <- scores + beta[, , p] * X[cell, p]
    }
    result[, cell, ] <- scores_to_weights(cbind(scores, 0))
  }
  if (rule == "global_after_h10") {
    beyond <- which(target$horizon > max(cfg$lfo_horizons))
    for (cell in beyond) {
      result[, cell, ] <- matrix(
        aggregation$global$weights[cfg$models],
        nrow = S,
        ncol = K,
        byrow = TRUE
      )
    }
  }
  validate_actuarial_weights(result, cfg$models, S)
  result
}

analysis_label <- function(rule) {
  if (rule == "within_lfo") "primary_h_le_10" else "exploratory"
}

variant_specs <- list(
  global = list(
    method = "stacking_global",
    weight_mode = "fixed_weights",
    output_group = "fixed"
  ),
  contextual = list(
    method = "stacking_contextual",
    weight_mode = "fixed_weights",
    output_group = "fixed"
  ),
  hierarchical_mean = list(
    method = "stacking_hierarchical",
    weight_mode = "posterior_mean_weights",
    output_group = "fixed"
  ),
  hierarchical_propagated = list(
    method = "stacking_hierarchical",
    weight_mode = "propagated_weights",
    output_group = "propagated"
  )
)

rules <- c("within_lfo", "extrapolate", "freeze_h10", "global_after_h10")
operational_rules <- c(
  within_lfo = "freeze_h10",
  extrapolate = "extrapolate",
  freeze_h10 = "freeze_h10",
  global_after_h10 = "global_after_h10"
)
fixed_weight_cache <- list()
propagated_weight_cache <- list()
for (operational_rule in unique(unname(operational_rules))) {
  fixed_weight_cache[[operational_rule]] <- mean_weight_sets_for_rule(
    target, operational_rule
  )
  propagated_weight_cache[[operational_rule]] <-
    hierarchical_weight_draws_for_rule(target, operational_rule)
  invisible(lapply(
    fixed_weight_cache[[operational_rule]],
    validate_actuarial_weights,
    models = cfg$models,
    S = S
  ))
}

survival_rows <- list()
annuity_rows <- list()
survival_draw_rows <- list(fixed = list(), propagated = list())
annuity_draw_rows <- list(fixed = list(), propagated = list())
distribution_store <- list()
fixed_weight_rows <- list()
propagated_weight_summary_rows <- list()
survival_position <- 1L
annuity_position <- 1L
survival_draw_position <- c(fixed = 1L, propagated = 1L)
annuity_draw_position <- c(fixed = 1L, propagated = 1L)
weight_position <- 1L
propagated_weight_position <- 1L

diagnostic_state <- list(
  minimum_force = Inf,
  maximum_force = -Inf,
  minimum_survival = Inf,
  maximum_survival = -Inf,
  maximum_survival_increase = -Inf,
  minimum_annuity = Inf,
  maximum_annuity_to_certain_ratio = -Inf
)

for (rule in rules) {
  operational_rule <- operational_rules[[rule]]
  fixed_weights <- fixed_weight_cache[[operational_rule]]
  propagated_weights <- propagated_weight_cache[[operational_rule]]

  for (method in names(fixed_weights)) {
    weights <- fixed_weights[[method]]
    N <- nrow(target)
    fixed_weight_rows[[weight_position]] <- data.frame(
      rule = rule,
      method = method,
      weight_mode = if (method == "stacking_hierarchical") {
        "posterior_mean_weights"
      } else {
        "fixed_weights"
      },
      age = target$age[rep(seq_len(N), times = K)],
      horizon = target$horizon[rep(seq_len(N), times = K)],
      model = rep(cfg$models, each = N),
      weight = as.vector(weights),
      stringsAsFactors = FALSE
    )
    weight_position <- weight_position + 1L
  }

  N <- nrow(target)
  for (k in seq_along(cfg$models)) {
    values <- propagated_weights[, , k, drop = FALSE][, , 1L]
    propagated_weight_summary_rows[[propagated_weight_position]] <-
      data.frame(
        rule = rule,
        method = "stacking_hierarchical",
        weight_mode = "propagated_weights",
        age = target$age,
        horizon = target$horizon,
        model = cfg$models[[k]],
        mean = colMeans(values),
        sd = apply(values, 2L, stats::sd),
        q025 = apply(values, 2L, stats::quantile, probs = 0.025),
        q975 = apply(values, 2L, stats::quantile, probs = 0.975),
        stringsAsFactors = FALSE
      )
    propagated_weight_position <- propagated_weight_position + 1L
  }

  for (initial_age in initial_ages) {
    all_indices <- which(target$initial_age == initial_age)
    indices <- if (rule == "within_lfo") {
      all_indices[
        target$horizon[all_indices] <= max(cfg$lfo_horizons)
      ]
    } else {
      all_indices
    }
    horizons <- target$horizon[indices]
    paired_forces <- pair_model_force_draws(
      model_force_draws, indices, cfg$models, permutations
    )
    store_key <- paste(rule, initial_age, sep = "|")
    distribution_store[[store_key]] <- list()

    for (variant in names(variant_specs)) {
      specification <- variant_specs[[variant]]
      weights <- if (variant == "hierarchical_propagated") {
        propagated_weights[, indices, , drop = FALSE]
      } else {
        fixed_weights[[specification$method]][
          indices, , drop = FALSE
        ]
      }
      force <- aggregate_force_draws(
        paired_forces, weights, cfg$models
      )
      survival <- survival_from_aggregated_force(force)
      diagnostic_state$minimum_force <- min(
        diagnostic_state$minimum_force, force
      )
      diagnostic_state$maximum_force <- max(
        diagnostic_state$maximum_force, force
      )
      diagnostic_state$minimum_survival <- min(
        diagnostic_state$minimum_survival, survival
      )
      diagnostic_state$maximum_survival <- max(
        diagnostic_state$maximum_survival, survival
      )
      if (ncol(survival) > 1L) {
        diagnostic_state$maximum_survival_increase <- max(
          diagnostic_state$maximum_survival_increase,
          survival[, -1L, drop = FALSE] -
            survival[, -ncol(survival), drop = FALSE]
        )
      }

      distribution_store[[store_key]][[variant]] <- list(
        survival = survival,
        annuity = list()
      )
      group <- specification$output_group
      draw_position <- survival_draw_position[[group]]
      survival_draw_rows[[group]][[draw_position]] <- data.frame(
        analysis = analysis_label(rule),
        rule = rule,
        method = specification$method,
        weight_mode = specification$weight_mode,
        initial_age = initial_age,
        draw = rep(seq_len(S), times = length(horizons)),
        horizon = rep(horizons, each = S),
        attained_age = rep(initial_age + horizons, each = S),
        survival = as.vector(survival),
        stringsAsFactors = FALSE
      )
      survival_draw_position[[group]] <- draw_position + 1L

      for (column in seq_along(horizons)) {
        summary <- summarize_actuarial_draws(survival[, column])
        survival_rows[[survival_position]] <- data.frame(
          analysis = analysis_label(rule),
          rule = rule,
          method = specification$method,
          weight_mode = specification$weight_mode,
          initial_age = initial_age,
          horizon = horizons[[column]],
          attained_age = initial_age + horizons[[column]],
          as.list(summary),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        survival_position <- survival_position + 1L
      }

      for (rate in discount_rates) {
        annuity <- annuity_from_survival(survival, horizons, rate)
        certain_annuity <- sum((1 / (1 + rate))^horizons)
        diagnostic_state$minimum_annuity <- min(
          diagnostic_state$minimum_annuity, annuity
        )
        diagnostic_state$maximum_annuity_to_certain_ratio <- max(
          diagnostic_state$maximum_annuity_to_certain_ratio,
          annuity / certain_annuity
        )
        rate_key <- sprintf("%.2f", rate)
        distribution_store[[store_key]][[variant]]$annuity[[
          rate_key
        ]] <- annuity

        draw_position <- annuity_draw_position[[group]]
        annuity_draw_rows[[group]][[draw_position]] <- data.frame(
          analysis = analysis_label(rule),
          rule = rule,
          method = specification$method,
          weight_mode = specification$weight_mode,
          initial_age = initial_age,
          duration = length(horizons),
          discount_rate = rate,
          draw = seq_len(S),
          annuity = annuity,
          stringsAsFactors = FALSE
        )
        annuity_draw_position[[group]] <- draw_position + 1L

        summary <- summarize_actuarial_draws(annuity)
        annuity_rows[[annuity_position]] <- data.frame(
          analysis = analysis_label(rule),
          rule = rule,
          method = specification$method,
          weight_mode = specification$weight_mode,
          initial_age = initial_age,
          duration = length(horizons),
          discount_rate = rate,
          as.list(summary),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        annuity_position <- annuity_position + 1L
      }
    }
  }
}

survival_summary <- do.call(rbind, survival_rows)
annuity_summary <- do.call(rbind, annuity_rows)
survival_draws_fixed <- do.call(rbind, survival_draw_rows$fixed)
survival_draws_propagated <- do.call(
  rbind, survival_draw_rows$propagated
)
annuity_draws_fixed <- do.call(rbind, annuity_draw_rows$fixed)
annuity_draws_propagated <- do.call(
  rbind, annuity_draw_rows$propagated
)
fixed_weights_long <- do.call(rbind, fixed_weight_rows)
propagated_weight_summary <- do.call(
  rbind, propagated_weight_summary_rows
)

fixed_mode <- survival_summary$weight_mode != "propagated_weights"
survival_summary_fixed <- survival_summary[fixed_mode, , drop = FALSE]
survival_summary_propagated <- survival_summary[
  !fixed_mode, , drop = FALSE
]
fixed_mode <- annuity_summary$weight_mode != "propagated_weights"
annuity_summary_fixed <- annuity_summary[fixed_mode, , drop = FALSE]
annuity_summary_propagated <- annuity_summary[
  !fixed_mode, , drop = FALSE
]

# Differences appariées : le meme indice de tirage et les memes permutations
# sont utilises pour les deux termes de chaque difference.
comparisons <- list(
  hierarchical_mean_minus_global =
    c("hierarchical_mean", "global"),
  hierarchical_propagated_minus_global =
    c("hierarchical_propagated", "global"),
  hierarchical_mean_minus_contextual =
    c("hierarchical_mean", "contextual"),
  hierarchical_propagated_minus_contextual =
    c("hierarchical_propagated", "contextual"),
  hierarchical_propagated_minus_hierarchical_mean =
    c("hierarchical_propagated", "hierarchical_mean")
)
paired_rows <- list()
paired_draw_rows <- list()
paired_position <- 1L
paired_draw_position <- 1L
for (rule in rules) {
  for (initial_age in initial_ages) {
    store <- distribution_store[[paste(rule, initial_age, sep = "|")]]
    for (comparison in names(comparisons)) {
      pair <- comparisons[[comparison]]
      for (rate in discount_rates) {
        rate_key <- sprintf("%.2f", rate)
        difference <- store[[pair[[1L]]]]$annuity[[rate_key]] -
          store[[pair[[2L]]]]$annuity[[rate_key]]
        summary <- summarize_paired_difference(difference)
        duration <- ncol(store[[pair[[1L]]]]$survival)
        paired_rows[[paired_position]] <- data.frame(
          analysis = analysis_label(rule),
          rule = rule,
          comparison = comparison,
          method = pair[[1L]],
          comparator = pair[[2L]],
          initial_age = initial_age,
          duration = duration,
          discount_rate = rate,
          as.list(summary),
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        paired_position <- paired_position + 1L
        paired_draw_rows[[paired_draw_position]] <- data.frame(
          analysis = analysis_label(rule),
          rule = rule,
          comparison = comparison,
          method = pair[[1L]],
          comparator = pair[[2L]],
          initial_age = initial_age,
          duration = duration,
          discount_rate = rate,
          draw = seq_len(S),
          difference = difference,
          stringsAsFactors = FALSE
        )
        paired_draw_position <- paired_draw_position + 1L
      }
    }
  }
}
paired_differences <- do.call(rbind, paired_rows)
paired_difference_draws <- do.call(rbind, paired_draw_rows)

# Diagnostic ancien/nouveau. Un identifiant commun rapproche les deux modes
# hierarchiques tout en distinguant poids moyens et poids propages.
variant_from_columns <- function(method, weight_mode) {
  result <- rep(NA_character_, length(method))
  result[method == "stacking_global"] <- "global"
  result[method == "stacking_contextual"] <- "contextual"
  result[
    method == "stacking_hierarchical" &
      weight_mode %in% c("posterior_mean", "posterior_mean_weights")
  ] <- "hierarchical_mean"
  result[
    method == "stacking_hierarchical" &
      weight_mode %in% c("posterior_draw", "propagated_weights")
  ] <- "hierarchical_propagated"
  result
}

if (has_old_actuarial_diagnostic) {
old_annuity$variant <- variant_from_columns(
  old_annuity$method, old_annuity$weight_mode
)
new_annuity <- annuity_summary
new_annuity$variant <- variant_from_columns(
  new_annuity$method, new_annuity$weight_mode
)
annuity_keys <- c(
  "analysis", "rule", "variant", "initial_age", "duration",
  "discount_rate"
)
old_annuity_diagnostic <- old_annuity[
  !is.na(old_annuity$variant),
  c(annuity_keys, "mean", "sd", "q025", "q975"),
  drop = FALSE
]
new_annuity_diagnostic <- new_annuity[
  !is.na(new_annuity$variant),
  c(annuity_keys, "mean", "sd", "q025", "q975"),
  drop = FALSE
]
annuity_old_new <- merge(
  old_annuity_diagnostic,
  new_annuity_diagnostic,
  by = annuity_keys,
  suffixes = c("_old_categorical", "_new_force"),
  all = TRUE,
  sort = FALSE
)
annuity_old_new$mean_difference <-
  annuity_old_new$mean_new_force -
  annuity_old_new$mean_old_categorical
annuity_old_new$mean_relative_difference_pct <- 100 *
  annuity_old_new$mean_difference /
  annuity_old_new$mean_old_categorical
annuity_old_new$width95_old_categorical <-
  annuity_old_new$q975_old_categorical -
  annuity_old_new$q025_old_categorical
annuity_old_new$width95_new_force <-
  annuity_old_new$q975_new_force -
  annuity_old_new$q025_new_force
annuity_old_new$width95_difference <-
  annuity_old_new$width95_new_force -
  annuity_old_new$width95_old_categorical

old_survival$variant <- variant_from_columns(
  old_survival$method, old_survival$weight_mode
)
new_survival <- survival_summary
new_survival$variant <- variant_from_columns(
  new_survival$method, new_survival$weight_mode
)
survival_keys <- c(
  "analysis", "rule", "variant", "initial_age", "horizon",
  "attained_age"
)
old_survival_diagnostic <- old_survival[
  !is.na(old_survival$variant),
  c(survival_keys, "mean", "sd", "q025", "q975"),
  drop = FALSE
]
new_survival_diagnostic <- new_survival[
  !is.na(new_survival$variant),
  c(survival_keys, "mean", "sd", "q025", "q975"),
  drop = FALSE
]
survival_old_new <- merge(
  old_survival_diagnostic,
  new_survival_diagnostic,
  by = survival_keys,
  suffixes = c("_old_categorical", "_new_force"),
  all = TRUE,
  sort = FALSE
)
survival_old_new$mean_difference <-
  survival_old_new$mean_new_force -
  survival_old_new$mean_old_categorical
survival_old_new$mean_relative_difference_pct <- 100 *
  survival_old_new$mean_difference /
  survival_old_new$mean_old_categorical
survival_old_new$width95_old_categorical <-
  survival_old_new$q975_old_categorical -
  survival_old_new$q025_old_categorical
survival_old_new$width95_new_force <-
  survival_old_new$q975_new_force -
  survival_old_new$q025_new_force
survival_old_new$width95_difference <-
  survival_old_new$width95_new_force -
  survival_old_new$width95_old_categorical
} else {
  annuity_old_new <- NULL
  survival_old_new <- NULL
}

# Contribution de l'incertitude des poids hierarchiques.
hierarchical_mean <- annuity_summary[
  annuity_summary$method == "stacking_hierarchical" &
    annuity_summary$weight_mode == "posterior_mean_weights",
  ,
  drop = FALSE
]
hierarchical_propagated <- annuity_summary[
  annuity_summary$method == "stacking_hierarchical" &
    annuity_summary$weight_mode == "propagated_weights",
  ,
  drop = FALSE
]
hierarchy_keys <- c(
  "analysis", "rule", "initial_age", "duration", "discount_rate"
)
hierarchy_uncertainty <- merge(
  hierarchical_mean[
    c(hierarchy_keys, "mean", "sd", "q10", "q90", "q025", "q975",
      "width80", "width95")
  ],
  hierarchical_propagated[
    c(hierarchy_keys, "mean", "sd", "q10", "q90", "q025", "q975",
      "width80", "width95")
  ],
  by = hierarchy_keys,
  suffixes = c("_mean_weights", "_propagated_weights"),
  sort = FALSE
)
hierarchy_uncertainty$mean_difference <-
  hierarchy_uncertainty$mean_propagated_weights -
  hierarchy_uncertainty$mean_mean_weights
hierarchy_uncertainty$sd_difference <-
  hierarchy_uncertainty$sd_propagated_weights -
  hierarchy_uncertainty$sd_mean_weights
hierarchy_uncertainty$width80_difference <-
  hierarchy_uncertainty$width80_propagated_weights -
  hierarchy_uncertainty$width80_mean_weights
hierarchy_uncertainty$width95_difference <-
  hierarchy_uncertainty$width95_propagated_weights -
  hierarchy_uncertainty$width95_mean_weights

# Controles automatiques.
check_rows <- list()
check_position <- 1L
add_check <- function(check, passed, value = NA_real_, tolerance = NA_real_,
                      detail = "") {
  check_rows[[check_position]] <<- data.frame(
    check = check,
    passed = isTRUE(passed),
    value = value,
    tolerance = tolerance,
    detail = detail,
    stringsAsFactors = FALSE
  )
  check_position <<- check_position + 1L
}

all_fixed_weights <- fixed_weights_long$weight
fixed_sums <- stats::aggregate(
  weight ~ rule + method + weight_mode + age + horizon,
  data = fixed_weights_long,
  FUN = sum
)
maximum_fixed_sum_error <- max(abs(fixed_sums$weight - 1))
maximum_propagated_sum_error <- max(vapply(
  propagated_weight_cache,
  function(values) max(abs(apply(values, c(1L, 2L), sum) - 1)),
  numeric(1)
))
add_check(
  "fixed_weights_sum_to_one",
  maximum_fixed_sum_error < 1e-8,
  maximum_fixed_sum_error,
  1e-8
)
add_check(
  "propagated_weights_sum_to_one",
  maximum_propagated_sum_error < 1e-8,
  maximum_propagated_sum_error,
  1e-8
)
add_check(
  "fixed_weights_in_zero_one",
  all(is.finite(all_fixed_weights)) &&
    all(all_fixed_weights >= -1e-8 & all_fixed_weights <= 1 + 1e-8),
  min(all_fixed_weights),
  1e-8,
  paste("maximum =", signif(max(all_fixed_weights), 8))
)
add_check(
  "model_order_matches_weights",
  identical(unique(fixed_weights_long$model), cfg$models),
  detail = paste(cfg$models, collapse = ",")
)
add_check(
  "attained_age_equals_x_plus_j_minus_one",
  all(target$age == target$initial_age + target$horizon - 1L)
)
add_check(
  "force_dimensions_identical",
  all(vapply(
    draw_dimensions,
    function(candidate) identical(candidate, draw_dimensions[[1L]]),
    logical(1)
  )),
  detail = paste(draw_dimensions[[1L]], collapse = " x ")
)
add_check(
  "forces_finite_and_nonnegative",
  is.finite(diagnostic_state$minimum_force) &&
    diagnostic_state$minimum_force >= -1e-10,
  diagnostic_state$minimum_force,
  1e-10,
  paste("maximum =", signif(diagnostic_state$maximum_force, 8))
)
add_check(
  "survival_in_zero_one",
  diagnostic_state$minimum_survival >= -1e-10 &&
    diagnostic_state$maximum_survival <= 1 + 1e-10,
  diagnostic_state$minimum_survival,
  1e-10,
  paste("maximum =", signif(diagnostic_state$maximum_survival, 8))
)
add_check(
  "survival_nonincreasing",
  diagnostic_state$maximum_survival_increase <= 1e-10,
  diagnostic_state$maximum_survival_increase,
  1e-10
)
add_check(
  "annuities_positive",
  diagnostic_state$minimum_annuity > 0,
  diagnostic_state$minimum_annuity
)
add_check(
  "annuities_below_certain_annuity",
  diagnostic_state$maximum_annuity_to_certain_ratio <= 1 + 1e-10,
  diagnostic_state$maximum_annuity_to_certain_ratio,
  1e-10
)
add_check(
  "first_payment_at_t_one",
  all(vapply(
    split(target$horizon, target$initial_age),
    function(horizons) min(horizons) == 1L,
    logical(1)
  )),
  detail = "discount factor v^1 times one-year survival"
)
permutations_repeated <- make_independent_draw_permutations(
  S, cfg$models, permutation_seed
)
add_check(
  "independent_pairing_reproducible",
  identical(permutations, permutations_repeated),
  detail = paste("seed", permutation_seed)
)
add_check(
  "hierarchical_surface_uses_one_coefficient_draw_per_trajectory",
  all(vapply(
    propagated_weight_cache,
    function(values) dim(values)[1L] == S &&
      dim(values)[2L] == nrow(target),
    logical(1)
  )),
  detail = "one posterior coefficient row s reused across every cell"
)

script_lines <- readLines(
  file.path("scripts", "06_actuarial_quantities.R"),
  warn = FALSE
)
forbidden_tokens <- vapply(
  list(
    c("sample", ".int("),
    c("run", "if("),
    c("select_models", "_with_uniforms"),
    c("mixture", "_force")
  ),
  paste0,
  collapse = "",
  FUN.VALUE = character(1)
)
residual_calls <- sum(vapply(
  forbidden_tokens,
  function(token) any(grepl(token, script_lines, fixed = TRUE)),
  logical(1)
))
add_check(
  "no_categorical_selection_in_actuarial_script",
  residual_calls == 0L,
  residual_calls,
  0
)

checks <- do.call(rbind, check_rows)

# Metadonnees, nombre exact de tirages et permutations conservees.
run_timestamp <- format(
  Sys.time(), "%Y-%m-%dT%H:%M:%S%z", tz = "Europe/Brussels"
)
metadata <- data.frame(
  run_timestamp = rep(run_timestamp, 2L),
  profile = cfg$profile,
  sex = cfg$sex,
  training_end = cfg$validation_end,
  projection_origin = cfg$validation_end,
  force_construction = "weighted_sum_of_model_forces",
  weight_mode = c("fixed_or_posterior_mean", "propagated_weights"),
  weight_uncertainty_propagated = c(FALSE, TRUE),
  model_pairing = "independent_saved_permutation_per_model",
  seed_pipeline = cfg$seed,
  seed_permutations = permutation_seed,
  seed_hierarchical_weight_draws = stable_seed(
    cfg$seed, cfg$profile, "actuarial_hierarchical_weight_surface"
  ),
  draws = S,
  models = paste(cfg$models, collapse = ","),
  initial_ages = paste(initial_ages, collapse = ","),
  primary_horizons = paste(
    seq_len(max(cfg$lfo_horizons)), collapse = ","
  ),
  exploratory_max_horizon = max(target$horizon),
  durations = paste(
    vapply(
      initial_ages,
      function(age) min(cfg$annuity_duration, cfg$age_max - age),
      integer(1)
    ),
    collapse = ","
  ),
  discount_rates = paste(discount_rates, collapse = ","),
  weights_source = normalizePath(
    file.path(cfg$paths$weights, "aggregation_results.rds"),
    winslash = "/",
    mustWork = TRUE
  ),
  excluded_cells = 0L,
  stringsAsFactors = FALSE
)
draw_counts <- do.call(rbind, lapply(
  names(variant_specs),
  function(variant) {
    specification <- variant_specs[[variant]]
    data.frame(
      variant = variant,
      method = specification$method,
      weight_mode = specification$weight_mode,
      survival_draws_per_age_rule_horizon = S,
      annuity_draws_per_age_rule_rate = S,
      stringsAsFactors = FALSE
    )
  }
))
permutations_long <- do.call(rbind, lapply(
  cfg$models,
  function(model) data.frame(
    model = model,
    aggregated_draw = seq_len(S),
    original_model_draw = permutations[[model]],
    seed = permutation_seed,
    stringsAsFactors = FALSE
  )
))
warnings_and_exclusions <- data.frame(
  level = character(),
  category = character(),
  message = character(),
  affected_cells = integer(),
  stringsAsFactors = FALSE
)

# Ecriture sous des noms nouveaux. Les anciens fichiers ne sont pas remplaces.
write_csv_atomic(
  survival_summary_fixed,
  file.path(
    metrics_dir,
    "actuarial_survival_summary_force_aggregation_fixed_weights.csv"
  )
)
write_csv_atomic(
  survival_summary_propagated,
  file.path(
    metrics_dir,
    "actuarial_survival_summary_force_aggregation_propagated_weights.csv"
  )
)
write_csv_gz_atomic(
  survival_draws_fixed,
  file.path(
    metrics_dir,
    "actuarial_survival_draws_force_aggregation_fixed_weights.csv.gz"
  )
)
write_csv_gz_atomic(
  survival_draws_propagated,
  file.path(
    metrics_dir,
    paste0(
      "actuarial_survival_draws_",
      "force_aggregation_propagated_weights.csv.gz"
    )
  )
)
write_csv_atomic(
  annuity_summary_fixed,
  file.path(
    metrics_dir,
    "actuarial_annuity_summary_force_aggregation_fixed_weights.csv"
  )
)
write_csv_atomic(
  annuity_summary_propagated,
  file.path(
    metrics_dir,
    "actuarial_annuity_summary_force_aggregation_propagated_weights.csv"
  )
)
write_csv_gz_atomic(
  annuity_draws_fixed,
  file.path(
    metrics_dir,
    "actuarial_annuity_draws_force_aggregation_fixed_weights.csv.gz"
  )
)
write_csv_gz_atomic(
  annuity_draws_propagated,
  file.path(
    metrics_dir,
    paste0(
      "actuarial_annuity_draws_",
      "force_aggregation_propagated_weights.csv.gz"
    )
  )
)
write_csv_atomic(
  paired_differences,
  file.path(
    metrics_dir,
    "actuarial_annuity_paired_differences_force_aggregation.csv"
  )
)
write_csv_gz_atomic(
  paired_difference_draws,
  file.path(
    metrics_dir,
    "actuarial_annuity_paired_difference_draws_force_aggregation.csv.gz"
  )
)
write_csv_atomic(
  hierarchy_uncertainty,
  file.path(
    metrics_dir,
    paste0(
      "actuarial_hierarchical_uncertainty_",
      "force_aggregation_propagated_weights.csv"
    )
  )
)
if (has_old_actuarial_diagnostic) {
  write_csv_atomic(
    annuity_old_new,
    file.path(
      metrics_dir,
      "actuarial_annuity_comparison_categorical_old_diagnostic.csv"
    )
  )
  write_csv_atomic(
    survival_old_new,
    file.path(
      metrics_dir,
      "actuarial_survival_comparison_categorical_old_diagnostic.csv"
    )
  )
}
write_csv_atomic(
  fixed_weights_long,
  file.path(
    metrics_dir,
    "actuarial_weights_force_aggregation_fixed_weights.csv"
  )
)
write_csv_atomic(
  propagated_weight_summary,
  file.path(
    metrics_dir,
    paste0(
      "actuarial_weight_summary_",
      "force_aggregation_propagated_weights.csv"
    )
  )
)
write_csv_atomic(
  checks,
  file.path(
    metrics_dir, "actuarial_force_aggregation_checks.csv"
  )
)
write_csv_atomic(
  metadata,
  file.path(
    metrics_dir, "actuarial_force_aggregation_metadata.csv"
  )
)
write_csv_atomic(
  draw_counts,
  file.path(
    metrics_dir, "actuarial_force_aggregation_draw_counts.csv"
  )
)
write_csv_atomic(
  permutations_long,
  file.path(
    metrics_dir, "actuarial_model_draw_permutations.csv"
  )
)
write_csv_atomic(
  warnings_and_exclusions,
  file.path(
    metrics_dir,
    "actuarial_force_aggregation_warnings_and_exclusions.csv"
  )
)
save_rds_atomic(
  list(
    distributions = distribution_store,
    permutations = permutations,
    metadata = metadata
  ),
  file.path(
    metrics_dir, "actuarial_distributions_force_aggregation.rds"
  ),
  compress = TRUE
)
save_rds_atomic(
  propagated_weight_cache,
  file.path(
    metrics_dir,
    paste0(
      "actuarial_hierarchical_weight_draws_",
      "force_aggregation_propagated_weights.rds"
    )
  ),
  compress = TRUE
)

# Figures principales : survie, distributions de rente, differences,
# incertitude des poids et diagnostic ancien/nouveau.
variant_labels <- c(
  global = "Global",
  contextual = "Contextuel",
  hierarchical_mean = "Hier. moyen",
  hierarchical_propagated = "Hier. propage"
)
variant_colors <- c(
  global = "#4477AA",
  contextual = "#EE7733",
  hierarchical_mean = "#228833",
  hierarchical_propagated = "#AA3377"
)
summary_variant <- variant_from_columns(
  survival_summary$method, survival_summary$weight_mode
)

render_chapter_figure(
  "actuarial_survival_force_aggregation_fixed_and_propagated",
  function() {
    graphics::layout(matrix(seq_along(initial_ages), nrow = 2L, byrow = TRUE))
    for (age_index in seq_along(initial_ages)) {
      initial_age <- initial_ages[[age_index]]
      part <- survival_summary[
        survival_summary$rule == "within_lfo" &
          survival_summary$initial_age == initial_age,
        ,
        drop = FALSE
      ]
      part$variant <- summary_variant[
        survival_summary$rule == "within_lfo" &
          survival_summary$initial_age == initial_age
      ]
      chapter_plot_theme(c(4.0, 4.0, 2.3, 0.8))
      graphics::plot(
        NA,
        xlim = range(part$horizon),
        ylim = range(part$q025, part$q975),
        xlab = "Horizon (annees)",
        ylab = "Probabilite de survie",
        main = paste("Age initial", initial_age)
      )
      for (variant in names(variant_labels)) {
        curve <- part[part$variant == variant, , drop = FALSE]
        curve <- curve[order(curve$horizon), , drop = FALSE]
        graphics::polygon(
          c(curve$horizon, rev(curve$horizon)),
          c(curve$q025, rev(curve$q975)),
          border = NA,
          col = grDevices::adjustcolor(
            variant_colors[[variant]], alpha.f = 0.10
          )
        )
        graphics::lines(
          curve$horizon, curve$mean,
          col = variant_colors[[variant]], lwd = 2
        )
      }
      if (age_index == 1L) {
        graphics::legend(
          "bottomleft",
          legend = unname(variant_labels),
          col = unname(variant_colors),
          lwd = 2,
          bty = "n",
          cex = 0.70
        )
      }
    }
  },
  figures_dir,
  width = 10,
  height = 7
)

render_chapter_figure(
  "actuarial_annuity_distributions_force_aggregation",
  function() {
    graphics::layout(matrix(seq_along(initial_ages), nrow = 2L, byrow = TRUE))
    for (age_index in seq_along(initial_ages)) {
      initial_age <- initial_ages[[age_index]]
      store <- distribution_store[[
        paste("within_lfo", initial_age, sep = "|")
      ]]
      draws <- lapply(
        names(variant_labels),
        function(variant) store[[variant]]$annuity[["0.02"]]
      )
      names(draws) <- unname(variant_labels)
      chapter_plot_theme(c(7.0, 4.0, 2.3, 0.8))
      graphics::boxplot(
        draws,
        col = grDevices::adjustcolor(
          unname(variant_colors), alpha.f = 0.55
        ),
        border = unname(variant_colors),
        outline = FALSE,
        las = 2,
        ylab = "Valeur de rente",
        main = paste("Age initial", initial_age),
        cex.axis = 0.72
      )
    }
  },
  figures_dir,
  width = 10,
  height = 8
)

comparison_labels <- c(
  hierarchical_mean_minus_global = "Hier. moyen - global",
  hierarchical_propagated_minus_global = "Hier. propage - global",
  hierarchical_mean_minus_contextual = "Hier. moyen - contextuel",
  hierarchical_propagated_minus_contextual =
    "Hier. propage - contextuel",
  hierarchical_propagated_minus_hierarchical_mean =
    "Hier. propage - hier. moyen"
)
render_chapter_figure(
  "actuarial_annuity_paired_differences_force_aggregation",
  function() {
    graphics::layout(matrix(1:6, nrow = 2L, byrow = TRUE))
    part_all <- paired_differences[
      paired_differences$rule == "within_lfo" &
        abs(paired_differences$discount_rate - 0.02) < 1e-12,
      ,
      drop = FALSE
    ]
    for (comparison in names(comparison_labels)) {
      part <- part_all[
        part_all$comparison == comparison, , drop = FALSE
      ]
      part <- part[order(part$initial_age), , drop = FALSE]
      chapter_plot_theme(c(4.0, 4.4, 2.5, 0.8))
      limits <- range(part$q025, part$q975, 0)
      graphics::plot(
        part$initial_age, part$mean,
        ylim = limits,
        type = "n",
        xlab = "Age initial",
        ylab = "Difference de rente",
        main = comparison_labels[[comparison]]
      )
      graphics::abline(h = 0, lty = 2, col = "#777777")
      graphics::segments(
        part$initial_age, part$q025,
        part$initial_age, part$q975,
        col = "#4477AA", lwd = 2
      )
      graphics::points(
        part$initial_age, part$mean,
        pch = 19, col = "#AA3377"
      )
    }
    graphics::plot.new()
  },
  figures_dir,
  width = 10,
  height = 7
)

render_chapter_figure(
  "actuarial_hierarchical_weight_uncertainty_force_aggregation",
  function() {
    part <- hierarchy_uncertainty[
      hierarchy_uncertainty$rule == "within_lfo" &
        abs(hierarchy_uncertainty$discount_rate - 0.02) < 1e-12,
      ,
      drop = FALSE
    ]
    part <- part[order(part$initial_age), , drop = FALSE]
    graphics::layout(matrix(1:2, nrow = 1L))
    chapter_plot_theme(c(5.2, 7.0, 2.8, 0.8))
    graphics::par(
      mgp = c(3.5, 1.2, 0),
      cex.lab = 0.85,
      cex.axis = 0.85,
      cex.main = 0.90
    )
    graphics::plot(
      part$initial_age, part$sd_difference,
      type = "b", pch = 19, col = "#AA3377",
      xlab = "Age initial", ylab = "Delta ecart-type",
      main = "Effet sur l'ecart-type"
    )
    graphics::abline(h = 0, lty = 2, col = "#777777")
    chapter_plot_theme(c(5.2, 7.0, 2.8, 0.8))
    graphics::par(
      mgp = c(3.5, 1.2, 0),
      cex.lab = 0.85,
      cex.axis = 0.85,
      cex.main = 0.90
    )
    graphics::plot(
      part$initial_age, part$width95_difference,
      type = "b", pch = 19, col = "#4477AA",
      xlab = "Age initial", ylab = "Delta largeur 95 %",
      main = "Effet sur l'intervalle a 95 %"
    )
    graphics::abline(h = 0, lty = 2, col = "#777777")
  },
  figures_dir,
  width = 10,
  height = 4.8
)

if (has_old_actuarial_diagnostic) {
render_chapter_figure(
  "actuarial_annuity_categorical_old_diagnostic",
  function() {
    selected_variants <- names(variant_labels)
    graphics::layout(matrix(seq_along(selected_variants), nrow = 2L))
    for (variant in selected_variants) {
      part <- annuity_old_new[
        annuity_old_new$rule == "extrapolate" &
          annuity_old_new$variant == variant &
          abs(annuity_old_new$discount_rate - 0.02) < 1e-12,
        ,
        drop = FALSE
      ]
      part <- part[order(part$initial_age), , drop = FALSE]
      chapter_plot_theme()
      limits <- range(
        part$mean_old_categorical, part$mean_new_force,
        finite = TRUE
      )
      graphics::plot(
        part$initial_age, part$mean_old_categorical,
        type = "b", pch = 1, lty = 2, col = "#777777",
        ylim = limits,
        xlab = "Age initial", ylab = "Rente moyenne",
        main = variant_labels[[variant]]
      )
      graphics::lines(
        part$initial_age, part$mean_new_force,
        type = "b", pch = 19, col = variant_colors[[variant]]
      )
      if (variant == selected_variants[[1L]]) {
        graphics::legend(
          "bottomleft",
          legend = c("Ancien categoriel", "Nouvelle force agregee"),
          col = c("#777777", variant_colors[[variant]]),
          lty = c(2, 1), pch = c(1, 19), bty = "n", cex = 0.8
        )
      }
    }
  },
  figures_dir,
  width = 9,
  height = 7
)
}

assert_true(
  all(checks$passed),
  paste(
    "Un controle actuariel a echoue :",
    paste(checks$check[!checks$passed], collapse = ", ")
  )
)

message(
  "Quantites actuarielles corrigees terminees : ", S,
  " tirages par methode, aucune cellule exclue, ",
  sum(checks$passed), "/", nrow(checks), " controles reussis."
)
