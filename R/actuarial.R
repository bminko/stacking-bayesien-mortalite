# Fonctions pures pour les trajectoires actuarielles.
#
# Les distributions predictives marginales restent des melanges de modeles
# dans R/metrics.R. Ici, une trajectoire actuarielle est construite en
# agregeant d'abord les forces de mortalite, puis en calculant la survie et
# la rente. Aucune selection categorielle de modele n'intervient.

validate_actuarial_weights <- function(weights, models, S = NULL,
                                       tolerance = 1e-8) {
  dimensions <- dim(weights)
  assert_true(
    length(dimensions) %in% c(2L, 3L),
    "Les poids actuariels doivent etre une matrice ou un tableau 3D."
  )
  model_dimension <- if (length(dimensions) == 2L) {
    dimensions[2L]
  } else {
    dimensions[3L]
  }
  assert_true(
    model_dimension == length(models),
    "Le nombre de colonnes de poids ne correspond pas aux modeles."
  )
  model_names <- if (length(dimensions) == 2L) {
    colnames(weights)
  } else {
    dimnames(weights)[[3L]]
  }
  assert_true(
    identical(model_names, models),
    "L'ordre des modeles et l'ordre des poids ne correspondent pas."
  )
  if (length(dimensions) == 3L && !is.null(S)) {
    assert_true(
      dimensions[1L] == S,
      "Le nombre de tirages de poids ne correspond pas aux forces."
    )
  }
  assert_true(
    all(is.finite(weights)),
    "Les poids actuariels contiennent une valeur non finie."
  )
  assert_true(
    all(weights >= -tolerance & weights <= 1 + tolerance),
    "Les poids actuariels doivent etre compris entre zero et un."
  )
  sums <- if (length(dimensions) == 2L) {
    rowSums(weights)
  } else {
    apply(weights, c(1L, 2L), sum)
  }
  assert_true(
    max(abs(sums - 1)) < tolerance,
    "La somme des poids actuariels n'est pas egale a un."
  )
  invisible(list(
    minimum = min(weights),
    maximum = max(weights),
    maximum_sum_error = max(abs(sums - 1))
  ))
}

make_independent_draw_permutations <- function(S, models, seed) {
  assert_true(S >= 1L, "Le nombre de tirages doit etre positif.")
  set.seed(seed)
  result <- setNames(
    lapply(models, function(model) sample.int(S, S, replace = FALSE)),
    models
  )
  assert_true(
    all(vapply(
      result,
      function(indices) identical(sort(indices), seq_len(S)),
      logical(1)
    )),
    "Une permutation de tirages est invalide."
  )
  result
}

pair_model_force_draws <- function(force_draws, indices, models,
                                   permutations) {
  assert_true(
    identical(names(force_draws), models),
    "L'ordre des matrices de forces ne correspond pas aux modeles."
  )
  assert_true(
    identical(names(permutations), models),
    "L'ordre des permutations ne correspond pas aux modeles."
  )
  dimensions <- lapply(force_draws, dim)
  assert_true(
    all(vapply(dimensions, length, integer(1)) == 2L),
    "Chaque jeu de forces doit etre une matrice."
  )
  reference <- dimensions[[1L]]
  assert_true(
    all(vapply(
      dimensions,
      function(candidate) identical(candidate, reference),
      logical(1)
    )),
    "Les dimensions des forces different entre les modeles."
  )
  assert_true(
    all(indices >= 1L & indices <= reference[2L]),
    "Un indice de cellule actuarielle est hors de la grille predictive."
  )
  assert_true(
    all(vapply(permutations, length, integer(1)) == reference[1L]),
    "Une permutation ne contient pas le bon nombre de tirages."
  )
  result <- setNames(lapply(models, function(model) {
    values <- force_draws[[model]][
      permutations[[model]], indices, drop = FALSE
    ]
    assert_true(
      all(is.finite(values) & values >= 0),
      paste("Forces invalides pour le modele", model)
    )
    values
  }), models)
  result
}

aggregate_force_draws <- function(paired_forces, weights, models,
                                  tolerance = 1e-8) {
  assert_true(
    identical(names(paired_forces), models),
    "L'ordre des forces appariees ne correspond pas aux modeles."
  )
  reference <- dim(paired_forces[[1L]])
  assert_true(
    all(vapply(
      paired_forces,
      function(values) identical(dim(values), reference),
      logical(1)
    )),
    "Les forces appariees n'ont pas les memes dimensions."
  )
  S <- reference[1L]
  J <- reference[2L]
  validate_actuarial_weights(weights, models, S, tolerance)
  assert_true(
    if (length(dim(weights)) == 2L) {
      nrow(weights) == J
    } else {
      dim(weights)[2L] == J
    },
    "Le nombre de cellules des poids ne correspond pas aux forces."
  )

  result <- matrix(0, nrow = S, ncol = J)
  for (k in seq_along(models)) {
    cell_weights <- if (length(dim(weights)) == 2L) {
      matrix(
        weights[, k],
        nrow = S,
        ncol = J,
        byrow = TRUE
      )
    } else {
      weights[, , k, drop = FALSE][, , 1L]
    }
    assert_true(
      identical(dim(cell_weights), reference),
      "Recyclage implicite detecte dans l'agregation des forces."
    )
    result <- result + paired_forces[[models[[k]]]] * cell_weights
  }
  assert_true(
    all(is.finite(result) & result >= -tolerance),
    "La force agregee contient une valeur invalide."
  )
  result[result < 0 & result >= -tolerance] <- 0
  result
}

survival_from_aggregated_force <- function(force, tolerance = 1e-10) {
  force <- as.matrix(force)
  assert_true(
    all(is.finite(force) & force >= -tolerance),
    "La survie requiert des forces finies et non negatives."
  )
  survival <- exp(-t(apply(force, 1L, cumsum)))
  assert_true(
    all(is.finite(survival)),
    "La survie contient une valeur non finie."
  )
  assert_true(
    all(survival >= -tolerance & survival <= 1 + tolerance),
    "Une probabilite de survie est hors de [0,1]."
  )
  if (ncol(survival) > 1L) {
    assert_true(
      max(survival[, -1L, drop = FALSE] -
            survival[, -ncol(survival), drop = FALSE]) <= tolerance,
      "Une trajectoire de survie n'est pas non croissante."
    )
  }
  survival
}

annuity_from_survival <- function(survival, horizons, interest_rate,
                                  tolerance = 1e-10) {
  survival <- as.matrix(survival)
  horizons <- as.integer(horizons)
  assert_true(
    ncol(survival) == length(horizons),
    "La duree de survie et les horizons de rente different."
  )
  assert_true(
    identical(horizons, seq_len(length(horizons))),
    "Les paiements de rente doivent commencer a t=1 sans lacune."
  )
  assert_true(
    is.finite(interest_rate) && interest_rate > -1,
    "Le taux d'actualisation est invalide."
  )
  discount_factors <- (1 / (1 + interest_rate))^horizons
  result <- as.numeric(survival %*% discount_factors)
  certain <- sum(discount_factors)
  assert_true(
    all(is.finite(result) & result > 0),
    "Une valeur de rente est non positive ou non finie."
  )
  assert_true(
    all(result <= certain + tolerance),
    "Une rente depasse la rente certaine de meme duree."
  )
  result
}

summarize_actuarial_draws <- function(values) {
  values <- as.numeric(values)
  assert_true(
    length(values) >= 2L && all(is.finite(values)),
    "La synthese requiert au moins deux tirages finis."
  )
  quantiles <- stats::quantile(
    values,
    probs = c(0.025, 0.05, 0.10, 0.90, 0.95, 0.975),
    names = FALSE
  )
  c(
    mean = mean(values),
    median = stats::median(values),
    sd = stats::sd(values),
    q025 = quantiles[1L],
    q05 = quantiles[2L],
    q10 = quantiles[3L],
    q90 = quantiles[4L],
    q95 = quantiles[5L],
    q975 = quantiles[6L],
    width80 = quantiles[4L] - quantiles[3L],
    width95 = quantiles[6L] - quantiles[1L]
  )
}

summarize_paired_difference <- function(values) {
  summary <- summarize_actuarial_draws(values)
  c(
    summary,
    probability_positive = mean(values > 0),
    probability_negative = mean(values < 0),
    probability_zero = mean(values == 0)
  )
}

write_csv_gz_atomic <- function(x, path, row.names = FALSE) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp")
  connection <- gzfile(temporary, open = "wt")
  on.exit({
    try(close(connection), silent = TRUE)
    if (file.exists(temporary)) file.remove(temporary)
  }, add = TRUE)
  utils::write.csv(
    x, connection, row.names = row.names, na = "", fileEncoding = "UTF-8"
  )
  close(connection)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(temporary, path)) {
    stop("Impossible d'ecrire le fichier compresse : ", path)
  }
  invisible(path)
}
