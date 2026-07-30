# Outils de restitution pour le chapitre 3.
#
# Ce fichier ne contient aucun ajustement statistique. Il transforme les
# objets deja valides du pipeline en tableaux et figures reproductibles.

chapter_method_labels <- function() {
  c(
    lc = "LC",
    rh = "RH",
    apc = "APC",
    cbd = "CBD",
    m6 = "M6",
    stacking_global = "Stacking global",
    stacking_contextual = "Contextuel non regularise",
    stacking_hierarchical = "Contextuel hierarchique"
  )
}

chapter_model_colors <- function() {
  c(
    lc = "#4477AA",
    rh = "#CC6677",
    apc = "#DDCC77",
    cbd = "#117733",
    m6 = "#AA4499",
    stacking_global = "#4477AA",
    stacking_contextual = "#EE7733",
    stacking_hierarchical = "#228833"
  )
}

render_chapter_figure <- function(name, draw, output_dir,
                                  width = 9, height = 6) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(output_dir, paste0(name, ".pdf"))
  png_path <- file.path(output_dir, paste0(name, ".png"))

  grDevices::pdf(
    pdf_path,
    width = width,
    height = height,
    useDingbats = FALSE,
    paper = "special"
  )
  draw()
  grDevices::dev.off()

  grDevices::png(
    png_path,
    width = round(width * 180),
    height = round(height * 180),
    res = 180
  )
  draw()
  grDevices::dev.off()
  invisible(c(pdf = pdf_path, png = png_path))
}

chapter_plot_theme <- function(mar = c(4.2, 4.5, 2.4, 1.2)) {
  graphics::par(
    mar = mar,
    las = 1,
    bty = "l",
    col.axis = "#30343B",
    col.lab = "#30343B",
    col.main = "#20242A",
    fg = "#30343B",
    family = "sans"
  )
}

weighted_column_mean <- function(values, weights) {
  sum(values * weights) / sum(weights)
}

posterior_fitted_summary <- function(fit, model, training,
                                     ndraws = 1000L, seed = 1L) {
  draws <- fit_draw_matrix(
    fit,
    ndraws = ndraws,
    seed = stable_seed(seed, model, "chapter_fitted")
  )
  S <- nrow(draws)
  A <- length(training$ages)
  T <- length(training$years)
  grid <- expand.grid(
    age = training$ages,
    year = training$years
  )
  eta <- matrix(NA_real_, nrow = S, ncol = nrow(grid))

  period_rows <- list()
  if (model %in% c("lc", "rh", "apc")) {
    kappa <- extract_vector_parameter(draws, "kappa", T)
    period_rows[[1L]] <- data.frame(
      model = model,
      factor = "kappa",
      year = training$years,
      mean = colMeans(kappa),
      q025 = apply(kappa, 2L, stats::quantile, probs = 0.025),
      q975 = apply(kappa, 2L, stats::quantile, probs = 0.975),
      stringsAsFactors = FALSE
    )
  } else {
    kappa1 <- do.call(cbind, lapply(seq_len(T), function(t) {
      extract_scalar_parameter(draws, sprintf("kappa[1,%d]", t))
    }))
    kappa2 <- do.call(cbind, lapply(seq_len(T), function(t) {
      extract_scalar_parameter(draws, sprintf("kappa[2,%d]", t))
    }))
    period_rows[[1L]] <- data.frame(
      model = model,
      factor = "kappa1",
      year = training$years,
      mean = colMeans(kappa1),
      q025 = apply(kappa1, 2L, stats::quantile, probs = 0.025),
      q975 = apply(kappa1, 2L, stats::quantile, probs = 0.975),
      stringsAsFactors = FALSE
    )
    period_rows[[2L]] <- data.frame(
      model = model,
      factor = "kappa2",
      year = training$years,
      mean = colMeans(kappa2),
      q025 = apply(kappa2, 2L, stats::quantile, probs = 0.025),
      q975 = apply(kappa2, 2L, stats::quantile, probs = 0.975),
      stringsAsFactors = FALSE
    )
  }

  cohort <- NULL
  cohort_values <- NULL
  if (model %in% c("rh", "apc", "m6")) {
    cohort_grid <- outer(
      training$ages,
      training$years,
      function(age, year) year - age
    )
    cohort_values <- seq.int(min(cohort_grid), max(cohort_grid))
    gamma <- extract_vector_parameter(
      draws, "gamma", length(cohort_values)
    )
    cohort <- data.frame(
      model = model,
      cohort = cohort_values,
      mean = colMeans(gamma),
      q025 = apply(gamma, 2L, stats::quantile, probs = 0.025),
      q975 = apply(gamma, 2L, stats::quantile, probs = 0.975),
      stringsAsFactors = FALSE
    )
  }

  if (model %in% c("lc", "rh")) {
    alpha <- extract_vector_parameter(draws, "alpha", A)
    beta <- extract_vector_parameter(draws, "beta", A)
    kappa <- extract_vector_parameter(draws, "kappa", T)
    for (column in seq_len(nrow(grid))) {
      age_index <- match(grid$age[column], training$ages)
      year_index <- match(grid$year[column], training$years)
      eta[, column] <- alpha[, age_index] +
        beta[, age_index] * kappa[, year_index]
      if (model == "rh") {
        cohort_index <- grid$year[column] - grid$age[column] -
          min(cohort_values) + 1L
        eta[, column] <- eta[, column] + gamma[, cohort_index]
      }
    }
  } else if (model == "apc") {
    alpha <- extract_vector_parameter(draws, "alpha", A)
    kappa <- extract_vector_parameter(draws, "kappa", T)
    gamma <- extract_vector_parameter(draws, "gamma", length(cohort_values))
    for (column in seq_len(nrow(grid))) {
      age_index <- match(grid$age[column], training$ages)
      year_index <- match(grid$year[column], training$years)
      cohort_index <- grid$year[column] - grid$age[column] -
        min(cohort_values) + 1L
      eta[, column] <- alpha[, age_index] +
        kappa[, year_index] + gamma[, cohort_index]
    }
  } else {
    kappa1 <- do.call(cbind, lapply(seq_len(T), function(t) {
      extract_scalar_parameter(draws, sprintf("kappa[1,%d]", t))
    }))
    kappa2 <- do.call(cbind, lapply(seq_len(T), function(t) {
      extract_scalar_parameter(draws, sprintf("kappa[2,%d]", t))
    }))
    gamma <- if (model == "m6") {
      extract_vector_parameter(draws, "gamma", length(cohort_values))
    } else {
      NULL
    }
    for (column in seq_len(nrow(grid))) {
      year_index <- match(grid$year[column], training$years)
      eta[, column] <- kappa1[, year_index] +
        (grid$age[column] - mean(training$ages)) *
        kappa2[, year_index]
      if (model == "m6") {
        cohort_index <- grid$year[column] - grid$age[column] -
          min(cohort_values) + 1L
        eta[, column] <- eta[, column] + gamma[, cohort_index]
      }
    }
  }

  fitted_force <- exp(eta)
  grid$mean_fitted_rate <- colMeans(fitted_force)
  grid$q025_fitted_rate <- apply(
    fitted_force, 2L, stats::quantile, probs = 0.025
  )
  grid$q975_fitted_rate <- apply(
    fitted_force, 2L, stats::quantile, probs = 0.975
  )
  list(
    fitted = grid,
    period = do.call(rbind, period_rows),
    cohort = cohort
  )
}

hierarchical_weight_summary <- function(grid, hierarchical_fit, models,
                                        constants, ndraws = 1000L,
                                        seed = 1L) {
  standardized <- standardize_context(grid, constants)$data
  X <- context_design_matrix(standardized)
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
      model = models,
      mean = colMeans(weights),
      median = apply(weights, 2L, stats::median),
      sd = apply(weights, 2L, stats::sd),
      q025 = apply(weights, 2L, stats::quantile, probs = 0.025),
      q975 = apply(weights, 2L, stats::quantile, probs = 0.975),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

hierarchical_coefficient_summary <- function(fit, models,
                                             ndraws = 1000L,
                                             seed = 1L) {
  draws <- fit_draw_matrix(
    fit,
    variables = c("alpha", "beta", "tau"),
    ndraws = ndraws,
    seed = seed
  )
  keep <- grep("^(alpha|beta|tau)\\[", colnames(draws), value = TRUE)
  rows <- lapply(keep, function(variable) {
    values <- draws[, variable]
    data.frame(
      method = "stacking_hierarchical",
      model = if (grepl("^(alpha|beta)\\[", variable)) {
        index <- as.integer(sub("^[^[]+\\[([0-9]+).*$", "\\1", variable))
        models[index]
      } else {
        NA_character_
      },
      term = variable,
      mean = mean(values),
      median = stats::median(values),
      sd = stats::sd(values),
      q025 = stats::quantile(values, 0.025),
      q975 = stats::quantile(values, 0.975),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

bootstrap_performance_differences <- function(cells, comparisons,
                                              repetitions = 2000L,
                                              seed = 1L) {
  metrics <- c("logs", "crps", "absolute_error_deaths")
  years <- sort(unique(cells$year))
  set.seed(seed)
  output <- list()
  position <- 1L

  for (comparison in names(comparisons)) {
    methods <- comparisons[[comparison]]
    first <- cells[cells$method == methods[[1L]], ]
    second <- cells[cells$method == methods[[2L]], ]
    key <- c("year", "age", "horizon")
    first <- first[do.call(order, unname(first[key])), ]
    second <- second[do.call(order, unname(second[key])), ]
    rownames(first) <- NULL
    rownames(second) <- NULL
    assert_true(
      all(vapply(key, function(column) {
        identical(first[[column]], second[[column]])
      }, logical(1))),
      paste("Cellules non alignees pour", comparison)
    )

    for (metric in metrics) {
      difference <- first[[metric]] - second[[metric]]
      bootstrap <- replicate(repetitions, {
        sampled_years <- sample(years, length(years), replace = TRUE)
        sampled_difference <- unlist(lapply(sampled_years, function(year) {
          difference[first$year == year]
        }))
        mean(sampled_difference)
      })
      output[[position]] <- data.frame(
        comparison = comparison,
        first_method = methods[[1L]],
        second_method = methods[[2L]],
        metric = metric,
        mean_difference = mean(difference),
        ci025 = stats::quantile(bootstrap, 0.025),
        ci975 = stats::quantile(bootstrap, 0.975),
        proportion_first_better = mean(difference < 0),
        bootstrap_unit = "target_year",
        bootstrap_repetitions = repetitions,
        stringsAsFactors = FALSE
      )
      position <- position + 1L
    }
  }
  do.call(rbind, output)
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\", "__LATEX_BACKSLASH__", x, fixed = TRUE)
  replacements <- c(
    "&" = "\\&",
    "%" = "\\%",
    "$" = "\\$",
    "#" = "\\#",
    "_" = "\\_",
    "{" = "\\{",
    "}" = "\\}",
    "~" = "\\textasciitilde{}",
    "^" = "\\textasciicircum{}"
  )
  for (pattern in names(replacements)) {
    x <- gsub(
      pattern,
      replacements[[pattern]],
      x,
      fixed = TRUE
    )
  }
  x <- gsub(
    "__LATEX_BACKSLASH__",
    "\\textbackslash{}",
    x,
    fixed = TRUE
  )
  x
}

write_latex_longtable <- function(data, path, headers = names(data),
                                  align = NULL, caption = NULL,
                                  label = NULL) {
  if (is.null(align)) {
    align <- paste0("l", paste(rep("r", ncol(data) - 1L), collapse = ""))
  }
  rows <- apply(data, 1L, function(row) {
    paste(latex_escape(row), collapse = " & ")
  })
  lines <- c(
    sprintf("\\begin{longtable}{%s}", align),
    if (!is.null(caption)) {
      paste0("\\caption{", latex_escape(caption), "}")
    },
    if (!is.null(label)) paste0("\\label{", label, "}\\\\"),
    paste(latex_escape(headers), collapse = " & "),
    "\\\\ \\toprule",
    "\\endfirsthead",
    paste(latex_escape(headers), collapse = " & "),
    "\\\\ \\toprule",
    "\\endhead",
    paste0(rows, " \\\\"),
    "\\bottomrule",
    "\\end{longtable}"
  )
  writeLines(lines[!is.na(lines)], path, useBytes = TRUE)
  invisible(path)
}
