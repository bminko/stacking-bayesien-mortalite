#!/usr/bin/env Rscript

# Genere les tableaux, figures et macros du rapport technique de simulation.
#
# Cette etape est strictement documentaire :
# - elle lit les 120 repetitions Stan deja validees ;
# - elle lit les dix controles complets deja valides ;
# - elle ne compile aucun modele et ne lance aucun echantillonnage MCMC.

options(stringsAsFactors = FALSE, scipen = 999)

source(file.path("R", "utils.R"))
source(file.path("R", "forecasting.R"))
source(file.path("R", "aggregation.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "simulation_study.R"))

cfg <- load_config("full")
root <- cfg$root
main_dir <- file.path(cfg$paths$simulation, "main")
validation_dir <- file.path(cfg$paths$simulation, "validation")
report_dir <- file.path(root, "output", "pdf", "simulation_technical_report")
tables_dir <- file.path(report_dir, "tables")
figures_dir <- file.path(report_dir, "figures")
derived_dir <- file.path(report_dir, "derived")
invisible(vapply(
  c(report_dir, tables_dir, figures_dir, derived_dir),
  dir.create, logical(1), recursive = TRUE, showWarnings = FALSE
))

read_csv <- function(path) {
  assert_true(file.exists(path), paste("Fichier absent :", path))
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

diagnostics <- read_csv(
  file.path(main_dir, "diagnostics_by_repetition.csv")
)
attempts <- read_csv(
  file.path(main_dir, "diagnostics_all_attempts.csv")
)
diagnostics_scenario <- read_csv(
  file.path(main_dir, "diagnostics_by_scenario.csv")
)
performance <- read_csv(
  file.path(main_dir, "performance_by_repetition.csv")
)
recovery <- read_csv(
  file.path(main_dir, "weight_recovery_by_repetition.csv")
)
parameters <- read_csv(
  file.path(main_dir, "parameters_by_repetition.csv")
)
tau <- read_csv(
  file.path(main_dir, "tau_by_repetition.csv")
)
weights <- read_csv(
  file.path(main_dir, "weights_by_repetition.csv")
)
validation_selection <- read_csv(
  file.path(validation_dir, "selection.csv")
)
validation_diagnostics <- read_csv(
  file.path(validation_dir, "diagnostics.csv")
)
validation_comparisons <- read_csv(
  file.path(validation_dir, "comparisons.csv")
)
validation_summary <- read_csv(
  file.path(validation_dir, "validation_summary.csv")
)

scenario_order <- c(
  "constant_weights", "horizon_only", "age_horizon", "low_information"
)
method_order <- c(
  "stacking_global", "stacking_contextual", "stacking_hierarchical"
)
model_order <- cfg$models
scenario_labels <- c(
  constant_weights = "Poids constants",
  horizon_only = "Poids selon l'horizon",
  age_horizon = "Poids selon l'age et l'horizon",
  low_information = "Faible information"
)
method_labels <- c(
  stacking_global = "Global",
  stacking_contextual = "Contextuel non regularise",
  stacking_hierarchical = "Contextuel hierarchique"
)
model_labels <- c(
  lc = "LC", rh = "RH", apc = "APC", cbd = "CBD", m6 = "M6"
)
feature_labels <- c(
  "Age", "Age quadratique", "Horizon", "Age x horizon"
)

# ---------------------------------------------------------------------------
# Audit de provenance et d'exhaustivite
# ---------------------------------------------------------------------------

result_files <- list.files(
  file.path(cfg$paths$simulation, "repetitions"),
  pattern = "result[.]rds$", recursive = TRUE, full.names = TRUE
)
assert_true(length(result_files) == 120L, "Le nombre de RDS principaux n'est pas 120.")
result_objects <- lapply(result_files, readRDS)
protocols <- vapply(
  result_objects, function(x) x$protocol_version, character(1)
)
cohorts <- vapply(result_objects, function(x) x$cohort, character(1))
statuses <- vapply(result_objects, function(x) x$status, character(1))
configurations <- vapply(
  result_objects, function(x) x$configuration, character(1)
)
object_methods <- sort(unique(unlist(lapply(
  result_objects, function(x) x$performance$method
))))

assert_true(
  all(protocols == cfg$simulation$protocol_version),
  "Un RDS principal n'appartient pas au protocole Stan v3."
)
assert_true(all(cohorts == "main"), "Une repetition n'appartient pas a la cohorte main.")
assert_true(all(statuses == "valide"), "Une repetition principale n'est pas valide.")
assert_true(
  identical(object_methods, sort(method_order)),
  "Les methodes des RDS principaux ne correspondent pas aux trois stackings."
)
assert_true(
  !any(vapply(
    result_objects,
    function(x) any(grepl("laplace", names(x), ignore.case = TRUE)),
    logical(1)
  )),
  "Une structure nommee Laplace apparait dans un RDS principal."
)
assert_true(nrow(diagnostics) == 120L, "Diagnostics principaux incomplets.")
assert_true(nrow(performance) == 360L, "Performances principales incompletes.")
assert_true(nrow(recovery) == 360L, "Recuperation des poids incomplete.")
assert_true(nrow(parameters) == 2880L, "Parametres hierarchiques incomplets.")
assert_true(nrow(tau) == 480L, "Parametres tau incomplets.")
assert_true(
  nrow(validation_selection) == 10L &&
    nrow(validation_diagnostics) == 10L &&
    validation_summary$repetitions_valid[[1L]] == 10L,
  "Controles complets incomplets."
)

audit_manifest <- data.frame(
  controle = c(
    "RDS principaux", "Protocole", "Cohorte", "Statut",
    "Configurations finales", "Methodes", "Reference Laplace",
    "Controles complets"
  ),
  resultat = c(
    "120/120 presents",
    unique(protocols),
    unique(cohorts),
    "120/120 valides",
    paste(names(table(configurations)), table(configurations), collapse = "; "),
    paste(object_methods, collapse = "; "),
    "Aucune dans les RDS utilises",
    "10/10 valides"
  ),
  source = c(
    "results/full/simulation/repetitions",
    "champ protocol_version des RDS",
    "champ cohort des RDS",
    "champ status des RDS",
    "champ configuration des RDS",
    "performance des RDS",
    "noms des structures des RDS",
    "results/full/simulation/validation"
  )
)
utils::write.csv(
  audit_manifest,
  file.path(derived_dir, "audit_manifest.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

# ---------------------------------------------------------------------------
# Fonctions de formatage et d'ecriture LaTeX
# ---------------------------------------------------------------------------

latex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\\\", "<<BS>>", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x <- gsub("<<BS>>", "\\\\textbackslash{}", x, fixed = TRUE)
  x
}

fmt_num <- function(x, digits = 3L) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    is.na(x), "",
    formatC(
      x, format = "f", digits = digits,
      big.mark = " ", decimal.mark = ","
    )
  )
}

fmt_sci <- function(x, digits = 3L) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    is.na(x), "",
    ifelse(
      abs(x) > 0 & abs(x) < 0.001,
      gsub("e", "\\\\,10^{", paste0(
        sub("e.*", "", formatC(x, format = "e", digits = digits)),
        "e",
        as.integer(sub(".*e", "", formatC(x, format = "e", digits = digits))),
        "}"
      ), fixed = TRUE),
      fmt_num(x, digits)
    )
  )
}

fmt_int <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(
    is.na(x), "",
    formatC(round(x), format = "d", big.mark = " ")
  )
}

fmt_pct <- function(x, digits = 1L) {
  ifelse(is.na(x), "", paste0(fmt_num(100 * x, digits), "%"))
}

fmt_pct_macro <- function(x, digits = 1L) {
  gsub("%", "\\%", fmt_pct(x, digits), fixed = TRUE)
}

label_scenario <- function(x) {
  unname(ifelse(x %in% names(scenario_labels), scenario_labels[x], x))
}

label_method <- function(x) {
  unname(ifelse(x %in% names(method_labels), method_labels[x], x))
}

label_configuration <- function(x) {
  labels <- c(
    legere = "Legere",
    complete = "Complete",
    complete_treedepth15 = "Complete, profondeur 15",
    validation_complete = "Controle complet",
    validation_complete_treedepth15 = "Controle complet, profondeur 15"
  )
  unname(ifelse(x %in% names(labels), labels[x], x))
}

write_longtable <- function(
    data, filename, caption, label, align = NULL,
    font_size = "\\scriptsize") {
  if (is.null(align)) {
    align <- paste(rep("l", ncol(data)), collapse = "")
  }
  headers <- latex_escape(names(data))
  body <- apply(data, 1L, function(row) {
    paste0(paste(latex_escape(row), collapse = " & "), " \\\\")
  })
  header_line <- paste0(paste(headers, collapse = " & "), " \\\\")
  lines <- c(
    font_size,
    paste0("\\begin{longtable}{", align, "}"),
    paste0("\\caption{", caption, "}\\label{", label, "}\\\\"),
    "\\toprule",
    header_line,
    "\\midrule",
    "\\endfirsthead",
    paste0(
      "\\multicolumn{", ncol(data),
      "}{c}{\\tablename\\ \\thetable\\ -- suite}\\\\"
    ),
    "\\toprule",
    header_line,
    "\\midrule",
    "\\endhead",
    paste0(
      "\\midrule\\multicolumn{", ncol(data),
      "}{r}{Suite page suivante}\\\\"
    ),
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot",
    body,
    "\\end{longtable}",
    "\\normalsize"
  )
  writeLines(
    lines, file.path(tables_dir, filename),
    useBytes = TRUE
  )
}

summary_stats <- function(x) {
  c(
    mean = mean(x),
    median = stats::median(x),
    sd = stats::sd(x),
    q025 = stats::quantile(x, 0.025, names = FALSE),
    q25 = stats::quantile(x, 0.25, names = FALSE),
    q75 = stats::quantile(x, 0.75, names = FALSE),
    q975 = stats::quantile(x, 0.975, names = FALSE)
  )
}

bind_group_summary <- function(data, groups, metrics) {
  split_key <- interaction(data[groups], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(data, split_key), function(part) {
    base <- part[1L, groups, drop = FALSE]
    values <- unlist(lapply(metrics, function(metric) {
      stats <- summary_stats(part[[metric]])
      names(stats) <- paste(metric, names(stats), sep = "_")
      stats
    }))
    cbind(base, as.data.frame(as.list(values), check.names = FALSE))
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

# ---------------------------------------------------------------------------
# Agregats principaux
# ---------------------------------------------------------------------------

performance_summary <- bind_group_summary(
  performance, c("scenario", "method"), c("mean_logs", "mean_crps")
)
recovery_summary <- bind_group_summary(
  recovery, c("scenario", "method"),
  c("weight_rmse", "mean_absolute_weight_error")
)
recovery_extra <- do.call(rbind, lapply(
  split(
    recovery,
    interaction(recovery$scenario, recovery$method, drop = TRUE)
  ),
  function(part) {
    data.frame(
      scenario = part$scenario[1L],
      method = part$method[1L],
      weight_bias_mean = mean(part$weight_bias),
      extreme_low_mean = mean(part$extreme_low_proportion),
      extreme_high_mean = mean(part$extreme_high_proportion),
      extreme_total_mean = mean(part$extreme_weight_proportion),
      coverage_mean = if (
        all(is.na(part$posterior_95_coverage))
      ) NA_real_ else mean(part$posterior_95_coverage, na.rm = TRUE),
      interval_width_mean = if (
        all(is.na(part$mean_interval_width))
      ) NA_real_ else mean(part$mean_interval_width, na.rm = TRUE)
    )
  }
))
recovery_summary <- merge(
  recovery_summary, recovery_extra,
  by = c("scenario", "method"), sort = FALSE
)

best_rows <- do.call(rbind, lapply(
  split(
    performance,
    interaction(performance$scenario, performance$repetition, drop = TRUE)
  ),
  function(part) {
    data.frame(
      scenario = part$scenario[1L],
      repetition = part$repetition[1L],
      best_logs = part$method[which.min(part$mean_logs)],
      best_crps = part$method[which.min(part$mean_crps)]
    )
  }
))
best_frequency <- do.call(rbind, lapply(scenario_order, function(scenario) {
  part <- best_rows[best_rows$scenario == scenario, ]
  do.call(rbind, lapply(method_order, function(method) {
    data.frame(
      scenario = scenario,
      method = method,
      best_logs_frequency = mean(part$best_logs == method),
      best_crps_frequency = mean(part$best_crps == method)
    )
  }))
}))

tau$feature_index <- as.integer(sub("tau\\[([0-9]+)\\]", "\\1", tau$parameter))
tau$feature <- feature_labels[tau$feature_index]
tau_summary <- bind_group_summary(
  tau, c("scenario", "feature"), c("mean")
)

beta <- parameters[grepl("^beta\\[", parameters$parameter), ]
beta$model_index <- as.integer(sub(
  "beta\\[([0-9]+),([0-9]+)\\]", "\\1", beta$parameter
))
beta$feature_index <- as.integer(sub(
  "beta\\[([0-9]+),([0-9]+)\\]", "\\2", beta$parameter
))
beta$model <- model_order[beta$model_index]
beta$feature <- feature_labels[beta$feature_index]
coefficient_summary <- do.call(rbind, lapply(
  split(
    beta,
    interaction(beta$scenario, beta$model, beta$feature, drop = TRUE)
  ),
  function(part) {
    stats <- summary_stats(part$mean)
    data.frame(
      scenario = part$scenario[1L],
      model = part$model[1L],
      feature = part$feature[1L],
      mean = stats["mean"],
      median = stats["median"],
      sd = stats["sd"],
      q025 = stats["q025"],
      q975 = stats["q975"]
    )
  }
))

initial_attempts <- attempts[attempts$configuration == "legere", ]
retry_keys <- unique(attempts[
  attempts$configuration != "legere",
  c("scenario", "repetition")
])
retry_initial <- merge(
  initial_attempts, retry_keys,
  by = c("scenario", "repetition")
)
retry_final <- merge(
  diagnostics, retry_keys,
  by = c("scenario", "repetition")
)
diagnostic_reason <- function(row) {
  reasons <- character()
  if (row[["max_rhat"]] > cfg$simulation$rhat_max) {
    reasons <- c(reasons, "R-hat")
  }
  if (row[["min_ess_bulk"]] < cfg$simulation$ess_min) {
    reasons <- c(reasons, "ESS bulk")
  }
  if (row[["min_ess_tail"]] < cfg$simulation$ess_min) {
    reasons <- c(reasons, "ESS tail")
  }
  if (row[["divergences"]] > 0) reasons <- c(reasons, "divergence")
  if (row[["max_treedepth_hits"]] > 0) reasons <- c(reasons, "treedepth")
  if (row[["min_ebfmi"]] < cfg$simulation$ebfmi_min) {
    reasons <- c(reasons, "E-BFMI")
  }
  if (isTRUE(row[["chain_stuck"]])) reasons <- c(reasons, "chaine bloquee")
  if (!length(reasons)) "Diagnostic global" else paste(reasons, collapse = " + ")
}
retry_initial$reason <- apply(retry_initial, 1L, diagnostic_reason)
retries <- merge(
  retry_initial[
    c(
      "scenario", "repetition", "reason", "max_rhat",
      "min_ess_bulk", "min_ess_tail", "min_ebfmi"
    )
  ],
  retry_final[
    c(
      "scenario", "repetition", "configuration", "max_rhat",
      "min_ess_bulk", "min_ess_tail", "min_ebfmi"
    )
  ],
  by = c("scenario", "repetition"),
  suffixes = c("_initial", "_final")
)

worst_diagnostics <- do.call(rbind, lapply(
  c("global", scenario_order),
  function(scenario) {
    part <- if (scenario == "global") {
      diagnostics
    } else {
      diagnostics[diagnostics$scenario == scenario, ]
    }
    data.frame(
      scenario = scenario,
      max_rhat = max(part$max_rhat),
      min_ess_bulk = min(part$min_ess_bulk),
      min_ess_tail = min(part$min_ess_tail),
      divergences = sum(part$divergences),
      treedepth = sum(part$max_treedepth_hits),
      min_ebfmi = min(part$min_ebfmi),
      elapsed_hours = sum(part$elapsed_seconds) / 3600
    )
  }
))

# ---------------------------------------------------------------------------
# Controles complete contre configuration principale acceptee
# ---------------------------------------------------------------------------

validation_weight_rows <- list()
validation_parameter_rows <- list()
validation_metrics <- list()

for (index in seq_len(nrow(validation_selection))) {
  selected <- validation_selection[index, ]
  scenario <- selected$scenario
  repetition <- selected$repetition
  original_path <- file.path(
    cfg$paths$simulation, "repetitions", scenario,
    sprintf("rep_%03d", repetition), "result.rds"
  )
  original <- readRDS(original_path)
  baseline <- original$weights[
    original$weights$method == "stacking_hierarchical", ,
    drop = FALSE
  ]
  full_path <- file.path(
    validation_dir, scenario, sprintf("rep_%03d", repetition),
    "weights_full.csv"
  )
  full <- read_csv(full_path)
  key <- c("age", "horizon", "model")
  joined <- merge(
    baseline, full, by = key,
    suffixes = c("_baseline", "_full"), sort = FALSE
  )
  comparison_row <- validation_comparisons[
    validation_comparisons$scenario == scenario &
      validation_comparisons$repetition == repetition, ,
    drop = FALSE
  ]
  joined$scenario <- scenario
  joined$repetition <- repetition
  joined$comparison_type <- comparison_row$comparison_type
  joined$difference <- joined$estimated_weight_full -
    joined$estimated_weight_baseline
  validation_weight_rows[[index]] <- joined

  baseline_parameters <- original$parameters
  full_parameters <- read_csv(file.path(
    validation_dir, scenario, sprintf("rep_%03d", repetition),
    "parameters_full.csv"
  ))
  parameter_join <- merge(
    baseline_parameters, full_parameters,
    by = "parameter", suffixes = c("_baseline", "_full")
  )
  parameter_join$scenario <- scenario
  parameter_join$repetition <- repetition
  parameter_join$comparison_type <- comparison_row$comparison_type
  validation_parameter_rows[[index]] <- parameter_join

  baseline_recovery <- original$recovery[
    original$recovery$method == "stacking_hierarchical", ,
    drop = FALSE
  ]
  full_difference <- joined$estimated_weight_full -
    joined$true_weight_full
  validation_metrics[[index]] <- data.frame(
    scenario = scenario,
    repetition = repetition,
    comparison_type = comparison_row$comparison_type,
    baseline_configuration = original$configuration,
    surface_correlation = stats::cor(
      joined$estimated_weight_baseline,
      joined$estimated_weight_full
    ),
    surface_rmse = sqrt(mean(joined$difference^2)),
    surface_mae = mean(abs(joined$difference)),
    surface_max_abs = max(abs(joined$difference)),
    bias_baseline = baseline_recovery$weight_bias,
    bias_full = mean(full_difference),
    rmse_baseline = baseline_recovery$weight_rmse,
    rmse_full = sqrt(mean(full_difference^2)),
    extreme_baseline = baseline_recovery$extreme_weight_proportion,
    extreme_full = mean(
      joined$estimated_weight_full < 0.01 |
        joined$estimated_weight_full > 0.99
    ),
    coverage_baseline = baseline_recovery$posterior_95_coverage,
    coverage_full = mean(
      joined$true_weight_full >= joined$weight_q025_full &
        joined$true_weight_full <= joined$weight_q975_full
    ),
    width_baseline = baseline_recovery$mean_interval_width,
    width_full = mean(
      joined$weight_q975_full - joined$weight_q025_full
    ),
    q025_mae = mean(abs(
      joined$weight_q025_full - joined$weight_q025_baseline
    )),
    q975_mae = mean(abs(
      joined$weight_q975_full - joined$weight_q975_baseline
    ))
  )
}

validation_weights <- do.call(rbind, validation_weight_rows)
validation_parameters <- do.call(rbind, validation_parameter_rows)
validation_metrics <- do.call(rbind, validation_metrics)
validation_metrics <- merge(
  validation_metrics,
  validation_comparisons,
  by = c("scenario", "repetition", "comparison_type"),
  sort = FALSE
)

validation_selection_report <- merge(
  validation_selection,
  diagnostics[
    c("scenario", "repetition", "seed")
  ],
  by = c("scenario", "repetition")
)
names(validation_selection_report)[
  names(validation_selection_report) == "seed"
] <- "main_seed"
validation_selection_report <- merge(
  validation_selection_report,
  validation_diagnostics[
    c("scenario", "repetition", "seed")
  ],
  by = c("scenario", "repetition")
)
names(validation_selection_report)[
  names(validation_selection_report) == "seed"
] <- "full_seed"

utils::write.csv(
  validation_metrics,
  file.path(derived_dir, "validation_metrics_extended.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

# ---------------------------------------------------------------------------
# Tableaux LaTeX
# ---------------------------------------------------------------------------

dgp_table <- data.frame(
  Scenario = unname(scenario_labels[scenario_order]),
  Objectif = c(
    "Verifier la recuperation de poids invariants",
    "Verifier une surface variant uniquement avec l'horizon",
    "Verifier une surface non separable en age et horizon",
    "Tester la regularisation lorsque les modeles sont peu discernables"
  ),
  `Poids vrais` = c(
    "LC 0,06; RH 0,54; APC 0,24; CBD 0,10; M6 0,06",
    "Softmax de scores lineaires en horizon",
    "Softmax de scores en age, age2, horizon et interaction",
    "LC 0,18; RH 0,24; APC 0,22; CBD 0,18; M6 0,18"
  ),
  Exposition = c("50 000", "50 000", "50 000", "5 000"),
  Dispersion = c("150", "150", "150", "55"),
  `Echelle du signal` = c("1", "1", "1", "0,12"),
  check.names = FALSE
)
write_longtable(
  dgp_table, "dgp_scenarios.tex",
  "Definition exacte des quatre scenarios simules.",
  "tab:dgp-scenarios", "p{2.7cm}p{4.3cm}p{5.5cm}rrr", "\\footnotesize"
)

method_table <- data.frame(
  Methode = unname(method_labels[method_order]),
  `Poids` = c(
    "Un vecteur constant de cinq poids",
    "Softmax d'un predicteur age-horizon sans regularisation",
    "Softmax d'un predicteur age-horizon avec regularisation hierarchique"
  ),
  Estimation = c(
    "BFGS, 2 demarrages dans la simulation",
    "BFGS, 2 demarrages dans la simulation",
    "MCMC Stan, 4 chaines"
  ),
  Incertitude = c(
    "Pas d'intervalle posterieur",
    "Pas d'intervalle posterieur",
    "Tirages posterieurs des poids, coefficients et tau"
  ),
  check.names = FALSE
)
write_longtable(
  method_table, "methods.tex",
  "Methodes de stacking effectivement comparees.",
  "tab:methods", "p{3.3cm}p{5.2cm}p{4.0cm}p{4.3cm}", "\\footnotesize"
)

diag_scenario_table <- diagnostics_scenario
diag_scenario_table$scenario <- label_scenario(diag_scenario_table$scenario)
names(diag_scenario_table) <- c(
  "Scenario", "Prevues", "Terminees", "Valides",
  "Relances", "Echecs", "Temps (h)"
)
diag_scenario_table$`Temps (h)` <- fmt_num(
  diag_scenario_table$`Temps (h)` / 3600, 2
)
write_longtable(
  diag_scenario_table, "diagnostics_by_scenario.tex",
  "Diagnostics synthetiques par scenario pour les 120 repetitions.",
  "tab:diagnostics-scenario", "lrrrrrr"
)

worst_table <- worst_diagnostics
worst_table$scenario <- ifelse(
  worst_table$scenario == "global",
  "Ensemble", label_scenario(worst_table$scenario)
)
worst_table <- data.frame(
  Scenario = worst_table$scenario,
  `R-hat max.` = fmt_num(worst_table$max_rhat, 4),
  `ESS bulk min.` = fmt_int(worst_table$min_ess_bulk),
  `ESS tail min.` = fmt_int(worst_table$min_ess_tail),
  Divergences = fmt_int(worst_table$divergences),
  Treedepth = fmt_int(worst_table$treedepth),
  `E-BFMI min.` = fmt_num(worst_table$min_ebfmi, 3),
  `Temps final (h)` = fmt_num(worst_table$elapsed_hours, 2),
  check.names = FALSE
)
write_longtable(
  worst_table, "worst_diagnostics.tex",
  "Pires diagnostics finaux et temps des configurations retenues.",
  "tab:worst-diagnostics", "lrrrrrrr"
)

retry_table <- data.frame(
  Scenario = label_scenario(retries$scenario),
  Rep. = fmt_int(retries$repetition),
  Raison = retries$reason,
  `R-hat init.` = fmt_num(retries$max_rhat_initial, 4),
  `ESS bulk init.` = fmt_int(retries$min_ess_bulk_initial),
  `ESS tail init.` = fmt_int(retries$min_ess_tail_initial),
  `E-BFMI init.` = fmt_num(retries$min_ebfmi_initial, 3),
  `Config. finale` = label_configuration(retries$configuration),
  `R-hat final` = fmt_num(retries$max_rhat_final, 4),
  `ESS bulk final` = fmt_int(retries$min_ess_bulk_final),
  `ESS tail final` = fmt_int(retries$min_ess_tail_final),
  `E-BFMI final` = fmt_num(retries$min_ebfmi_final, 3),
  check.names = FALSE
)
write_longtable(
  retry_table, "retries.tex",
  "Dix repetitions relancees : cause initiale et diagnostic final.",
  "tab:retries", "llp{2.4cm}rrrrlrrrr", "\\tiny"
)

performance_logs_table <- data.frame(
  Scenario = label_scenario(performance_summary$scenario),
  Methode = label_method(performance_summary$method),
  Moyenne = fmt_num(performance_summary$mean_logs_mean, 4),
  Mediane = fmt_num(performance_summary$mean_logs_median, 4),
  `Ecart-type` = fmt_num(performance_summary$mean_logs_sd, 4),
  `Q2,5` = fmt_num(performance_summary$mean_logs_q025, 4),
  Q25 = fmt_num(performance_summary$mean_logs_q25, 4),
  Q75 = fmt_num(performance_summary$mean_logs_q75, 4),
  `Q97,5` = fmt_num(performance_summary$mean_logs_q975, 4),
  check.names = FALSE
)
write_longtable(
  performance_logs_table, "performance_logs_summary.tex",
  "Distribution du LogS moyen sur 30 repetitions par scenario.",
  "tab:logs-summary", "llrrrrrrr"
)

performance_crps_table <- data.frame(
  Scenario = label_scenario(performance_summary$scenario),
  Methode = label_method(performance_summary$method),
  Moyenne = fmt_num(performance_summary$mean_crps_mean, 7),
  Mediane = fmt_num(performance_summary$mean_crps_median, 7),
  `Ecart-type` = fmt_num(performance_summary$mean_crps_sd, 7),
  `Q2,5` = fmt_num(performance_summary$mean_crps_q025, 7),
  Q25 = fmt_num(performance_summary$mean_crps_q25, 7),
  Q75 = fmt_num(performance_summary$mean_crps_q75, 7),
  `Q97,5` = fmt_num(performance_summary$mean_crps_q975, 7),
  check.names = FALSE
)
write_longtable(
  performance_crps_table, "performance_crps_summary.tex",
  "Distribution du CRPS moyen sur 30 repetitions par scenario.",
  "tab:crps-summary", "llrrrrrrr"
)

best_frequency_table <- data.frame(
  Scenario = label_scenario(best_frequency$scenario),
  Methode = label_method(best_frequency$method),
  `Meilleur LogS` = fmt_pct(best_frequency$best_logs_frequency, 1),
  `Meilleur CRPS` = fmt_pct(best_frequency$best_crps_frequency, 1),
  check.names = FALSE
)
write_longtable(
  best_frequency_table, "best_frequency.tex",
  "Frequence a laquelle chaque methode obtient le meilleur score.",
  "tab:best-frequency", "llrr"
)

recovery_table <- data.frame(
  Scenario = label_scenario(recovery_summary$scenario),
  Methode = label_method(recovery_summary$method),
  `RMSE moy.` = fmt_num(recovery_summary$weight_rmse_mean, 4),
  `RMSE med.` = fmt_num(recovery_summary$weight_rmse_median, 4),
  `RMSE sd` = fmt_num(recovery_summary$weight_rmse_sd, 4),
  `RMSE Q2,5` = fmt_num(recovery_summary$weight_rmse_q025, 4),
  `RMSE Q97,5` = fmt_num(recovery_summary$weight_rmse_q975, 4),
  `Biais moy.` = fmt_num(recovery_summary$weight_bias_mean, 5),
  `<0,01` = fmt_pct(recovery_summary$extreme_low_mean, 1),
  `>0,99` = fmt_pct(recovery_summary$extreme_high_mean, 1),
  Couverture = fmt_pct(recovery_summary$coverage_mean, 1),
  Largeur = fmt_num(recovery_summary$interval_width_mean, 3),
  check.names = FALSE
)
write_longtable(
  recovery_table, "recovery_summary.tex",
  paste(
    "Recuperation des poids. La couverture et la largeur ne sont definies",
    "que pour le modele hierarchique."
  ),
  "tab:recovery-summary", "llrrrrrrrrrr", "\\tiny"
)

tau_table <- data.frame(
  Scenario = label_scenario(tau_summary$scenario),
  Effet = tau_summary$feature,
  Moyenne = fmt_num(tau_summary$mean_mean, 3),
  Mediane = fmt_num(tau_summary$mean_median, 3),
  `Ecart-type` = fmt_num(tau_summary$mean_sd, 3),
  `Q2,5` = fmt_num(tau_summary$mean_q025, 3),
  `Q97,5` = fmt_num(tau_summary$mean_q975, 3),
  check.names = FALSE
)
write_longtable(
  tau_table, "tau_summary.tex",
  "Distribution entre repetitions des moyennes posterieures de tau.",
  "tab:tau-summary", "llrrrrr"
)

coefficient_table <- data.frame(
  Scenario = label_scenario(coefficient_summary$scenario),
  Modele = unname(model_labels[coefficient_summary$model]),
  Effet = coefficient_summary$feature,
  Moyenne = fmt_num(coefficient_summary$mean, 3),
  Mediane = fmt_num(coefficient_summary$median, 3),
  `Ecart-type` = fmt_num(coefficient_summary$sd, 3),
  `Q2,5` = fmt_num(coefficient_summary$q025, 3),
  `Q97,5` = fmt_num(coefficient_summary$q975, 3),
  check.names = FALSE
)
write_longtable(
  coefficient_table, "coefficient_summary.tex",
  paste(
    "Resume entre repetitions des coefficients contextuels beta.",
    "M6 est le modele de reference."
  ),
  "tab:coefficient-summary", "lllrrrrr", "\\tiny"
)

selection_table <- data.frame(
  Scenario = label_scenario(validation_selection_report$scenario),
  Rep. = fmt_int(validation_selection_report$repetition),
  `Graine principale` = fmt_int(validation_selection_report$main_seed),
  `Graine controle` = fmt_int(validation_selection_report$full_seed),
  `Config. principale` = label_configuration(
    validation_selection_report$main_configuration
  ),
  `Relance anterieure` = ifelse(
    validation_selection_report$prior_retry, "Oui", "Non"
  ),
  `RMSE poids` = fmt_num(validation_selection_report$weight_rmse, 4),
  `Raison de selection` = gsub(
    "_", " ", validation_selection_report$selection_reason
  ),
  check.names = FALSE
)
write_longtable(
  selection_table, "validation_selection.tex",
  "Dix repetitions retenues pour les controles de robustesse numerique.",
  "tab:validation-selection", "llrrllrp{4.0cm}", "\\tiny"
)

validation_surface_table <- data.frame(
  Scenario = label_scenario(validation_metrics$scenario),
  Rep. = fmt_int(validation_metrics$repetition),
  Comparaison = ifelse(
    validation_metrics$comparison_type == "legere_vs_complete",
    "Legere vs complete", "Complete vs complete"
  ),
  Correlation = fmt_num(validation_metrics$surface_correlation, 5),
  `RMSE surfaces` = fmt_num(validation_metrics$surface_rmse, 5),
  `Ecart abs. moyen` = fmt_num(validation_metrics$surface_mae, 5),
  `Ecart abs. max.` = fmt_num(validation_metrics$surface_max_abs, 5),
  `RMSE base` = fmt_num(validation_metrics$rmse_baseline, 5),
  `RMSE complete` = fmt_num(validation_metrics$rmse_full, 5),
  check.names = FALSE
)
write_longtable(
  validation_surface_table, "validation_surfaces.tex",
  "Stabilite des surfaces de poids entre configuration principale et controle.",
  "tab:validation-surfaces", "llp{2.5cm}rrrrrr", "\\tiny"
)

validation_recovery_table <- data.frame(
  Scenario = label_scenario(validation_metrics$scenario),
  Rep. = fmt_int(validation_metrics$repetition),
  `Biais base` = fmt_num(validation_metrics$bias_baseline, 6),
  `Biais complet` = fmt_num(validation_metrics$bias_full, 6),
  `Extremes base` = fmt_pct(validation_metrics$extreme_baseline, 2),
  `Extremes complet` = fmt_pct(validation_metrics$extreme_full, 2),
  `Couv. base` = fmt_pct(validation_metrics$coverage_baseline.x, 2),
  `Couv. complete` = fmt_pct(validation_metrics$coverage_full.x, 2),
  `Largeur base` = fmt_num(validation_metrics$width_baseline, 4),
  `Largeur complete` = fmt_num(validation_metrics$width_full, 4),
  check.names = FALSE
)
write_longtable(
  validation_recovery_table, "validation_recovery.tex",
  paste(
    "Stabilite du biais, des poids extremes, de la couverture",
    "et de la largeur des intervalles credibles."
  ),
  "tab:validation-recovery", "llrrrrrrrr", "\\tiny"
)

validation_parameter_table <- data.frame(
  Scenario = label_scenario(validation_metrics$scenario),
  Rep. = fmt_int(validation_metrics$repetition),
  `MAE borne 2,5` = fmt_num(validation_metrics$q025_mae, 5),
  `MAE borne 97,5` = fmt_num(validation_metrics$q975_mae, 5),
  `MAE beta` = fmt_num(
    validation_metrics$contextual_coefficient_mean_absolute_difference, 5
  ),
  `MAE tau` = fmt_num(
    validation_metrics$tau_mean_absolute_difference, 5
  ),
  check.names = FALSE
)
write_longtable(
  validation_parameter_table, "validation_parameter_stability.tex",
  paste(
    "Stabilite des bornes des intervalles credibles,",
    "des coefficients contextuels et des parametres tau."
  ),
  "tab:validation-parameters", "llrrrr", "\\tiny"
)

validation_predictive_table <- data.frame(
  Scenario = label_scenario(validation_metrics$scenario),
  Rep. = fmt_int(validation_metrics$repetition),
  `LogS base` = fmt_num(validation_metrics$logs_baseline, 5),
  `LogS complet` = fmt_num(validation_metrics$logs_full, 5),
  `Delta LogS rel.` = fmt_pct(
    validation_metrics$relative_logs_difference, 3
  ),
  `CRPS base` = fmt_num(validation_metrics$crps_baseline, 7),
  `CRPS complet` = fmt_num(validation_metrics$crps_full, 7),
  `Delta CRPS rel.` = fmt_pct(
    validation_metrics$relative_crps_difference, 3
  ),
  `Rang LogS` = ifelse(
    validation_metrics$logs_ranking_identical, "Stable", "Modifie"
  ),
  `Rang CRPS` = ifelse(
    validation_metrics$crps_ranking_identical, "Stable", "Modifie"
  ),
  check.names = FALSE
)
write_longtable(
  validation_predictive_table, "validation_predictive.tex",
  "Stabilite des scores predictifs et des classements.",
  "tab:validation-predictive", "llrrrrrrll", "\\tiny"
)

validation_diag_table <- data.frame(
  Scenario = label_scenario(validation_diagnostics$scenario),
  Rep. = fmt_int(validation_diagnostics$repetition),
  `R-hat max.` = fmt_num(validation_diagnostics$max_rhat, 4),
  `ESS bulk min.` = fmt_int(validation_diagnostics$min_ess_bulk),
  `ESS tail min.` = fmt_int(validation_diagnostics$min_ess_tail),
  Divergences = fmt_int(validation_diagnostics$divergences),
  Treedepth = fmt_int(validation_diagnostics$max_treedepth_hits),
  `E-BFMI min.` = fmt_num(validation_diagnostics$min_ebfmi, 3),
  Statut = ifelse(validation_diagnostics$valid, "Valide", "Invalide"),
  check.names = FALSE
)
write_longtable(
  validation_diag_table, "validation_diagnostics.tex",
  "Diagnostics MCMC des dix controles complets.",
  "tab:validation-diagnostics", "llrrrrrrl"
)

audit_table <- audit_manifest
names(audit_table) <- c("Controle", "Resultat", "Source auditee")
write_longtable(
  audit_table, "audit_manifest.tex",
  "Audit de provenance des resultats utilises dans le rapport.",
  "tab:audit", "p{3.5cm}p{5.0cm}p{7.0cm}", "\\footnotesize"
)

# Annexes detaillees.
diag_detail <- data.frame(
  Scenario = label_scenario(diagnostics$scenario),
  Rep. = fmt_int(diagnostics$repetition),
  Graine = fmt_int(diagnostics$seed),
  Configuration = label_configuration(diagnostics$configuration),
  `Temps (min)` = fmt_num(diagnostics$elapsed_seconds / 60, 1),
  `R-hat` = fmt_num(diagnostics$max_rhat, 4),
  `ESS bulk` = fmt_int(diagnostics$min_ess_bulk),
  `ESS tail` = fmt_int(diagnostics$min_ess_tail),
  Div. = fmt_int(diagnostics$divergences),
  Prof. = fmt_int(diagnostics$max_treedepth_hits),
  `E-BFMI` = fmt_num(diagnostics$min_ebfmi, 3),
  Statut = ifelse(diagnostics$valid, "Valide", "Invalide"),
  check.names = FALSE
)
write_longtable(
  diag_detail, "diagnostics_120.tex",
  "Diagnostics individuels des 120 repetitions principales.",
  "tab:diagnostics-120", "llrlrrrrrrrl", "\\tiny"
)

performance_detail <- merge(
  performance,
  best_rows,
  by = c("scenario", "repetition")
)
performance_detail <- data.frame(
  Scenario = label_scenario(performance_detail$scenario),
  Rep. = fmt_int(performance_detail$repetition),
  Methode = label_method(performance_detail$method),
  LogS = fmt_num(performance_detail$mean_logs, 5),
  CRPS = fmt_num(performance_detail$mean_crps, 7),
  `Meilleur LogS` = ifelse(
    performance_detail$method == performance_detail$best_logs, "Oui", ""
  ),
  `Meilleur CRPS` = ifelse(
    performance_detail$method == performance_detail$best_crps, "Oui", ""
  ),
  check.names = FALSE
)
write_longtable(
  performance_detail, "performance_360.tex",
  "Scores individuels des trois methodes dans les 120 repetitions.",
  "tab:performance-360", "lllrlll", "\\tiny"
)

recovery_detail <- data.frame(
  Scenario = label_scenario(recovery$scenario),
  Rep. = fmt_int(recovery$repetition),
  Methode = label_method(recovery$method),
  RMSE = fmt_num(recovery$weight_rmse, 5),
  MAE = fmt_num(recovery$mean_absolute_weight_error, 5),
  Biais = fmt_num(recovery$weight_bias, 6),
  `<0,01` = fmt_pct(recovery$extreme_low_proportion, 1),
  `>0,99` = fmt_pct(recovery$extreme_high_proportion, 1),
  Couverture = fmt_pct(recovery$posterior_95_coverage, 1),
  Largeur = fmt_num(recovery$mean_interval_width, 3),
  check.names = FALSE
)
write_longtable(
  recovery_detail, "recovery_360.tex",
  "Recuperation individuelle des poids dans les 120 repetitions.",
  "tab:recovery-360", "lllrrrrrrr", "\\tiny"
)

tau_detail <- data.frame(
  Scenario = label_scenario(tau$scenario),
  Rep. = fmt_int(tau$repetition),
  Configuration = label_configuration(tau$configuration),
  Effet = tau$feature,
  Moyenne = fmt_num(tau$mean, 4),
  `Q2,5` = fmt_num(tau$q025, 4),
  Mediane = fmt_num(tau$q50, 4),
  `Q97,5` = fmt_num(tau$q975, 4),
  check.names = FALSE
)
write_longtable(
  tau_detail, "tau_480.tex",
  "Parametres tau individuels des 120 repetitions.",
  "tab:tau-480", "lllrrrrr", "\\tiny"
)

# ---------------------------------------------------------------------------
# Macros numeriques
# ---------------------------------------------------------------------------

main_attempt_hours <- sum(attempts$elapsed_seconds) / 3600
main_retries <- sum(attempts$configuration != "legere")
light_validation <- validation_metrics[
  validation_metrics$comparison_type == "legere_vs_complete", ,
  drop = FALSE
]

macro_lines <- c(
  sprintf("\\newcommand{\\MainFits}{%d}", nrow(diagnostics)),
  sprintf("\\newcommand{\\MainRetries}{%d}", main_retries),
  sprintf("\\newcommand{\\MainRhat}{%s}", fmt_num(max(diagnostics$max_rhat), 4)),
  sprintf("\\newcommand{\\MainEssBulk}{%s}", fmt_int(min(diagnostics$min_ess_bulk))),
  sprintf("\\newcommand{\\MainEssTail}{%s}", fmt_int(min(diagnostics$min_ess_tail))),
  sprintf("\\newcommand{\\MainEbfmi}{%s}", fmt_num(min(diagnostics$min_ebfmi), 3)),
  sprintf("\\newcommand{\\MainHours}{%s}", fmt_num(main_attempt_hours, 2)),
  sprintf("\\newcommand{\\ValidationRhat}{%s}", fmt_num(max(validation_diagnostics$max_rhat), 4)),
  sprintf("\\newcommand{\\ValidationEssBulk}{%s}", fmt_int(min(validation_diagnostics$min_ess_bulk))),
  sprintf("\\newcommand{\\ValidationEssTail}{%s}", fmt_int(min(validation_diagnostics$min_ess_tail))),
  sprintf("\\newcommand{\\ValidationEbfmi}{%s}", fmt_num(min(validation_diagnostics$min_ebfmi), 3)),
  sprintf(
    "\\newcommand{\\ValidationWeightMae}{%s}",
    fmt_num(mean(light_validation$surface_mae), 5)
  ),
  sprintf(
    "\\newcommand{\\ValidationWeightMaxMean}{%s}",
    fmt_num(max(light_validation$surface_mae), 5)
  ),
  sprintf(
    "\\newcommand{\\ValidationPointMax}{%s}",
    fmt_num(max(light_validation$surface_max_abs), 5)
  ),
  sprintf(
    "\\newcommand{\\ValidationCorrelation}{%s}",
    fmt_num(min(light_validation$surface_correlation), 4)
  ),
  sprintf(
    "\\newcommand{\\ValidationLogsRelative}{%s}",
    fmt_pct_macro(mean(light_validation$relative_logs_difference), 3)
  ),
  sprintf(
    "\\newcommand{\\ValidationCrpsRelative}{%s}",
    fmt_pct_macro(mean(light_validation$relative_crps_difference), 3)
  )
)
writeLines(
  macro_lines,
  file.path(report_dir, "report_macros.tex"),
  useBytes = TRUE
)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

method_colors <- c(
  stacking_global = "#4477AA",
  stacking_contextual = "#EE7733",
  stacking_hierarchical = "#AA3377"
)
scenario_colors <- c(
  constant_weights = "#4477AA",
  horizon_only = "#CCBB44",
  age_horizon = "#EE7733",
  low_information = "#AA3377"
)
comparison_colors <- c(
  legere_vs_complete = "#4477AA",
  relance_complete_vs_complete = "#EE7733"
)

open_figure <- function(filename, width = 11, height = 7) {
  grDevices::cairo_pdf(
    file.path(figures_dir, filename),
    width = width, height = height, family = "sans"
  )
}

close_figure <- function() {
  grDevices::dev.off()
}

plot_group_box <- function(values, groups, labels, colors, ylab, main) {
  graphics::boxplot(
    split(values, factor(groups, levels = names(labels))),
    names = unname(labels),
    col = unname(colors[names(labels)]),
    border = "#333333", outline = FALSE,
    ylab = ylab, main = main, las = 2, cex.axis = 0.75
  )
  graphics::grid(nx = NA, ny = NULL, col = "#E5E7EB")
}

# Carte des graphiques, conservee pour l'audit de conception.
chart_map <- data.frame(
  figure = c(
    "diagnostics_overview.pdf", "rmse_distributions.pdf",
    "predictive_performance.pdf", "recovery_quality.pdf",
    "best_frequency.pdf", "tau_distributions.pdf",
    "coefficient_magnitudes.pdf", "low_information_detail.pdf",
    paste0("surfaces_", scenario_order, ".pdf"),
    paste0("surface_errors_", scenario_order, ".pdf"),
    "validation_weights_scatter.pdf", "validation_weight_heatmaps.pdf",
    "validation_metrics.pdf", "validation_tau.pdf"
  ),
  famille = c(
    "diagnostic", "distribution", "distribution", "comparaison",
    "comparaison", "distribution", "comparaison", "distribution",
    rep("matrice", 4), rep("matrice", 4),
    "relation", "matrice", "comparaison", "relation"
  ),
  question = c(
    "Les 120 ajustements convergent-ils ?",
    "Quelle methode recupere le mieux les poids ?",
    "Quelle methode predit le mieux ?",
    "Comment regularisation, couverture et extremes se comparent-ils ?",
    "Quelle methode gagne le plus souvent ?",
    "Quels effets contextuels sont le plus regularises ?",
    "Quelle est l'amplitude des coefficients contextuels ?",
    "Que se passe-t-il lorsque l'information est faible ?",
    rep("Les surfaces vraies sont-elles recuperees ?", 4),
    rep("Ou les erreurs de surface se concentrent-elles ?", 4),
    "Les poids principaux et complets concordent-ils ?",
    "Ou se trouvent les ecarts des dix controles ?",
    "Les ecarts respectent-ils les seuils de stabilite ?",
    "Les parametres tau sont-ils stables ?"
  )
)
utils::write.csv(
  chart_map, file.path(derived_dir, "chart_map.csv"),
  row.names = FALSE, fileEncoding = "UTF-8"
)

# Diagnostics.
open_figure("diagnostics_overview.pdf", 11, 7.5)
graphics::par(mfrow = c(2, 2), mar = c(6, 4.2, 2.5, 1))
graphics::barplot(
  diagnostics_scenario$reestimations,
  names.arg = label_scenario(diagnostics_scenario$scenario),
  col = scenario_colors[diagnostics_scenario$scenario],
  border = "#333333", las = 2, ylab = "Nombre de relances",
  main = "Relances ciblees par scenario"
)
graphics::grid(nx = NA, ny = NULL, col = "#E5E7EB")
plot_group_box(
  diagnostics$max_rhat, diagnostics$scenario,
  scenario_labels, scenario_colors,
  "R-hat maximal", "Distribution des R-hat maximaux"
)
graphics::abline(h = c(1.01, 1.05), lty = c(2, 3), col = c("#555555", "#AA3377"))
ess_long <- rbind(
  data.frame(type = "ESS bulk", value = diagnostics$min_ess_bulk),
  data.frame(type = "ESS tail", value = diagnostics$min_ess_tail)
)
graphics::boxplot(
  value ~ type, data = ess_long,
  col = c("#4477AA", "#EE7733"), border = "#333333",
  outline = FALSE, ylab = "ESS minimal",
  main = "ESS minimaux sur les parametres et les poids"
)
graphics::abline(h = 400, lty = 2, col = "#AA3377")
plot_group_box(
  diagnostics$min_ebfmi, diagnostics$scenario,
  scenario_labels, scenario_colors,
  "E-BFMI minimal", "Distribution des E-BFMI minimaux"
)
graphics::abline(h = 0.30, lty = 2, col = "#AA3377")
close_figure()

# RMSE des poids.
open_figure("rmse_distributions.pdf", 11, 7.5)
graphics::par(mfrow = c(2, 2), mar = c(4.5, 4.2, 2.5, 1))
for (scenario in scenario_order) {
  part <- recovery[recovery$scenario == scenario, ]
  plot_group_box(
    part$weight_rmse, part$method,
    method_labels, method_colors,
    "RMSE des poids", scenario_labels[[scenario]]
  )
}
close_figure()

# Scores predictifs.
open_figure("predictive_performance.pdf", 12, 7)
graphics::par(mfrow = c(2, 4), mar = c(4.8, 4.0, 2.2, 0.8))
for (metric in c("mean_logs", "mean_crps")) {
  for (scenario in scenario_order) {
    part <- performance[performance$scenario == scenario, ]
    plot_group_box(
      part[[metric]], part$method,
      method_labels, method_colors,
      if (metric == "mean_logs") "LogS moyen" else "CRPS moyen",
      scenario_labels[[scenario]]
    )
  }
}
close_figure()

# Extremes, couverture et largeur.
open_figure("recovery_quality.pdf", 11, 4.2)
graphics::par(mfrow = c(1, 3), mar = c(6.5, 4.2, 2.5, 1))
means_extreme <- stats::aggregate(
  extreme_weight_proportion ~ scenario + method,
  recovery, mean
)
for (metric in c("extreme_weight_proportion")) {
  matrix_values <- xtabs(
    means_extreme[[metric]] ~ means_extreme$method + means_extreme$scenario
  )[method_order, scenario_order]
  graphics::barplot(
    matrix_values, beside = TRUE,
    col = method_colors[method_order], border = "#333333",
    names.arg = label_scenario(scenario_order),
    las = 2, ylab = "Frequence moyenne",
    main = "Poids inferieurs a 0,01 ou superieurs a 0,99"
  )
}
hierarchical_recovery <- recovery[
  recovery$method == "stacking_hierarchical", ]
plot_group_box(
  hierarchical_recovery$posterior_95_coverage,
  hierarchical_recovery$scenario,
  scenario_labels, scenario_colors,
  "Couverture", "Couverture des intervalles a 95 %"
)
graphics::abline(h = 0.95, lty = 2, col = "#555555")
plot_group_box(
  hierarchical_recovery$mean_interval_width,
  hierarchical_recovery$scenario,
  scenario_labels, scenario_colors,
  "Largeur moyenne", "Largeur des intervalles credibles"
)
close_figure()

# Frequence du meilleur score.
open_figure("best_frequency.pdf", 11, 4.5)
graphics::par(mfrow = c(1, 2), mar = c(6.5, 4.2, 2.5, 1))
for (metric in c("best_logs_frequency", "best_crps_frequency")) {
  matrix_values <- xtabs(
    best_frequency[[metric]] ~
      best_frequency$method + best_frequency$scenario
  )[method_order, scenario_order]
  graphics::barplot(
    matrix_values, beside = TRUE,
    col = method_colors[method_order], border = "#333333",
    names.arg = label_scenario(scenario_order),
    las = 2, ylim = c(0, 1),
    ylab = "Frequence",
    main = if (
      metric == "best_logs_frequency"
    ) "Meilleur LogS" else "Meilleur CRPS"
  )
  graphics::legend(
    "topright", legend = method_labels[method_order],
    fill = method_colors[method_order], bty = "n", cex = 0.75
  )
}
close_figure()

# Tau.
open_figure("tau_distributions.pdf", 11, 7.5)
graphics::par(mfrow = c(2, 2), mar = c(6, 4.2, 2.5, 1))
for (scenario in scenario_order) {
  part <- tau[tau$scenario == scenario, ]
  graphics::boxplot(
    split(part$mean, factor(part$feature, levels = feature_labels)),
    names = feature_labels, col = scenario_colors[[scenario]],
    border = "#333333", outline = FALSE, las = 2,
    ylab = "Moyenne posterieure de tau",
    main = scenario_labels[[scenario]], cex.axis = 0.8
  )
  graphics::grid(nx = NA, ny = NULL, col = "#E5E7EB")
}
close_figure()

# Amplitude des coefficients beta.
coefficient_magnitude <- stats::aggregate(
  abs(mean) ~ scenario + feature, beta, mean
)
open_figure("coefficient_magnitudes.pdf", 11, 7.5)
graphics::par(mfrow = c(2, 2), mar = c(6, 4.2, 2.5, 1))
for (scenario in scenario_order) {
  part <- coefficient_magnitude[
    coefficient_magnitude$scenario == scenario, ]
  part <- part[match(feature_labels, part$feature), ]
  graphics::barplot(
    part[[3L]], names.arg = part$feature,
    col = scenario_colors[[scenario]], border = "#333333",
    las = 2, ylab = "Moyenne de |beta|",
    main = scenario_labels[[scenario]]
  )
  graphics::grid(nx = NA, ny = NULL, col = "#E5E7EB")
}
close_figure()

# Focus low_information.
low_recovery <- recovery[recovery$scenario == "low_information", ]
low_performance <- performance[performance$scenario == "low_information", ]
open_figure("low_information_detail.pdf", 11, 7.5)
graphics::par(mfrow = c(2, 2), mar = c(5.2, 4.2, 2.5, 1))
plot_group_box(
  low_recovery$weight_rmse, low_recovery$method,
  method_labels, method_colors,
  "RMSE des poids", "Recuperation des poids"
)
plot_group_box(
  low_recovery$extreme_weight_proportion, low_recovery$method,
  method_labels, method_colors,
  "Frequence", "Poids extremes"
)
plot_group_box(
  low_performance$mean_logs, low_performance$method,
  method_labels, method_colors,
  "LogS moyen", "Performance LogS"
)
plot_group_box(
  low_performance$mean_crps, low_performance$method,
  method_labels, method_colors,
  "CRPS moyen", "Performance CRPS"
)
close_figure()

# Surfaces moyennes vraies, estimees et erreurs.
surface_summary <- stats::aggregate(
  weights[c("true_weight", "estimated_weight")],
  by = weights[c("scenario", "method", "age", "horizon", "model")],
  FUN = mean
)
weight_palette <- grDevices::hcl.colors(40, "YlGnBu")
error_palette <- grDevices::colorRampPalette(
  c("#2166AC", "#F7F7F7", "#B2182B")
)(41)

for (scenario in scenario_order) {
  part <- surface_summary[surface_summary$scenario == scenario, ]
  open_figure(paste0("surfaces_", scenario, ".pdf"), 11.7, 8.2)
  graphics::par(
    mfrow = c(4, 5), mar = c(1.8, 1.8, 2.2, 0.8),
    oma = c(2.3, 2.8, 2.4, 1)
  )
  row_methods <- c("true", method_order)
  row_labels <- c("Poids vrais", unname(method_labels[method_order]))
  for (row_index in seq_along(row_methods)) {
    for (model in model_order) {
      model_part <- part[part$model == model, ]
      if (row_methods[row_index] == "true") {
        model_part <- model_part[
          model_part$method == "stacking_hierarchical", ]
        value <- "true_weight"
      } else {
        model_part <- model_part[
          model_part$method == row_methods[row_index], ]
        value <- "estimated_weight"
      }
      matrix_values <- xtabs(
        model_part[[value]] ~ model_part$age + model_part$horizon
      )
      graphics::image(
        as.numeric(rownames(matrix_values)),
        as.numeric(colnames(matrix_values)),
        matrix_values,
        col = weight_palette, zlim = c(0, 1),
        xlab = "", ylab = "", axes = FALSE
      )
      if (row_index == 1L) {
        graphics::title(main = model_labels[[model]], cex.main = 0.9)
      }
      if (row_index == length(row_methods)) {
        graphics::axis(1, at = c(50, 70, 90), cex.axis = 0.7)
      }
      if (model == model_order[1L]) {
        graphics::axis(2, at = c(1, 5, 10), cex.axis = 0.7)
        graphics::mtext(
          row_labels[row_index], side = 2, line = 2.2, cex = 0.75
        )
      }
      graphics::box(col = "#333333")
    }
  }
  graphics::mtext("Age", side = 1, outer = TRUE, line = 0.8)
  graphics::mtext("Horizon", side = 2, outer = TRUE, line = 1.1)
  graphics::mtext(
    scenario_labels[[scenario]], side = 3, outer = TRUE,
    line = 0.6, font = 2
  )
  close_figure()

  error_part <- part
  error_part$error <- error_part$estimated_weight - error_part$true_weight
  maximum <- max(abs(error_part$error))
  open_figure(paste0("surface_errors_", scenario, ".pdf"), 11.7, 6.8)
  graphics::par(
    mfrow = c(3, 5), mar = c(1.8, 1.8, 2.2, 0.8),
    oma = c(2.3, 2.8, 2.4, 1)
  )
  for (row_index in seq_along(method_order)) {
    for (model in model_order) {
      model_part <- error_part[
        error_part$method == method_order[row_index] &
          error_part$model == model, ]
      matrix_values <- xtabs(
        model_part$error ~ model_part$age + model_part$horizon
      )
      graphics::image(
        as.numeric(rownames(matrix_values)),
        as.numeric(colnames(matrix_values)),
        matrix_values,
        col = error_palette, zlim = c(-maximum, maximum),
        xlab = "", ylab = "", axes = FALSE
      )
      if (row_index == 1L) {
        graphics::title(main = model_labels[[model]], cex.main = 0.9)
      }
      if (row_index == length(method_order)) {
        graphics::axis(1, at = c(50, 70, 90), cex.axis = 0.7)
      }
      if (model == model_order[1L]) {
        graphics::axis(2, at = c(1, 5, 10), cex.axis = 0.7)
        graphics::mtext(
          method_labels[[method_order[row_index]]],
          side = 2, line = 2.2, cex = 0.75
        )
      }
      graphics::box(col = "#333333")
    }
  }
  graphics::mtext("Age", side = 1, outer = TRUE, line = 0.8)
  graphics::mtext("Horizon", side = 2, outer = TRUE, line = 1.1)
  graphics::mtext(
    paste("Erreur estime - vrai :", scenario_labels[[scenario]]),
    side = 3, outer = TRUE, line = 0.6, font = 2
  )
  close_figure()
}

# Concordance des poids lors des controles.
open_figure("validation_weights_scatter.pdf", 11, 5.2)
graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
for (comparison_type in names(comparison_colors)) {
  part <- validation_weights[
    validation_weights$comparison_type == comparison_type, ]
  graphics::plot(
    part$estimated_weight_baseline,
    part$estimated_weight_full,
    pch = 16, cex = 0.25,
    col = grDevices::adjustcolor(
      comparison_colors[[comparison_type]], alpha.f = 0.18
    ),
    xlim = c(0, 1), ylim = c(0, 1),
    xlab = "Poids de la configuration principale",
    ylab = "Poids du controle complet",
    main = if (
      comparison_type == "legere_vs_complete"
    ) "Legere contre complete (7 cas)" else
      "Complete contre complete (3 relances)"
  )
  graphics::abline(0, 1, lty = 2, col = "#333333")
  graphics::grid(col = "#E5E7EB")
}
close_figure()

# Heatmaps des ecarts RMS entre modeles pour chaque cellule.
cell_differences <- stats::aggregate(
  difference ~ scenario + repetition + comparison_type + age + horizon,
  validation_weights,
  function(x) sqrt(mean(x^2))
)
max_cell_difference <- max(cell_differences$difference)
open_figure("validation_weight_heatmaps.pdf", 12, 5.4)
graphics::par(
  mfrow = c(2, 5), mar = c(2.0, 2.0, 2.4, 0.8),
  oma = c(2.2, 2.5, 1.2, 1)
)
for (index in seq_len(nrow(validation_selection))) {
  selected <- validation_selection[index, ]
  part <- cell_differences[
    cell_differences$scenario == selected$scenario &
      cell_differences$repetition == selected$repetition, ]
  matrix_values <- xtabs(part$difference ~ part$age + part$horizon)
  graphics::image(
    as.numeric(rownames(matrix_values)),
    as.numeric(colnames(matrix_values)),
    matrix_values,
    col = grDevices::hcl.colors(40, "YlOrRd"),
    zlim = c(0, max_cell_difference),
    xlab = "", ylab = "", axes = FALSE
  )
  graphics::axis(1, at = c(50, 70, 90), cex.axis = 0.65)
  graphics::axis(2, at = c(1, 5, 10), cex.axis = 0.65)
  graphics::title(
    main = paste0(
      selected$scenario, " - rep. ", selected$repetition
    ),
    cex.main = 0.68
  )
  graphics::box(col = "#333333")
}
graphics::mtext("Age", side = 1, outer = TRUE, line = 0.8)
graphics::mtext("Horizon", side = 2, outer = TRUE, line = 0.8)
close_figure()

# Metriques de stabilite.
validation_metrics$label <- paste0(
  validation_metrics$scenario, "-", validation_metrics$repetition
)
open_figure("validation_metrics.pdf", 12, 7)
graphics::par(mfrow = c(2, 2), mar = c(7.2, 4.2, 2.5, 1))
metric_specs <- list(
  list(
    field = "surface_mae", title = "Ecart absolu moyen des poids",
    ylab = "Ecart absolu", threshold = 0.02
  ),
  list(
    field = "relative_logs_difference", title = "Difference relative de LogS",
    ylab = "Difference relative", threshold = 0.02
  ),
  list(
    field = "relative_crps_difference", title = "Difference relative de CRPS",
    ylab = "Difference relative", threshold = 0.02
  )
)
for (spec in metric_specs) {
  colors <- comparison_colors[validation_metrics$comparison_type]
  graphics::barplot(
    validation_metrics[[spec$field]],
    names.arg = validation_metrics$label,
    col = colors, border = "#333333", las = 2,
    ylab = spec$ylab, main = spec$title
  )
  graphics::abline(h = spec$threshold, lty = 2, col = "#AA3377")
  graphics::grid(nx = NA, ny = NULL, col = "#E5E7EB")
}
ranking_values <- c(
  mean(validation_metrics$logs_ranking_identical),
  mean(validation_metrics$crps_ranking_identical)
)
graphics::barplot(
  ranking_values, names.arg = c("LogS", "CRPS"),
  col = c("#4477AA", "#EE7733"), border = "#333333",
  ylim = c(0, 1), ylab = "Proportion",
  main = "Classements integralement stables"
)
graphics::abline(h = 0.80, lty = 2, col = "#AA3377")
close_figure()

# Tau principal contre tau complet.
validation_tau <- validation_parameters[
  grepl("^tau\\[", validation_parameters$parameter), ]
open_figure("validation_tau.pdf", 11, 5.2)
graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
tau_range <- range(
  validation_tau$mean_baseline,
  validation_tau$mean_full
)
for (comparison_type in names(comparison_colors)) {
  part <- validation_tau[
    validation_tau$comparison_type == comparison_type, ]
  graphics::plot(
    part$mean_baseline, part$mean_full,
    pch = 21, bg = comparison_colors[[comparison_type]],
    col = "#333333", xlim = tau_range, ylim = tau_range,
    xlab = "Moyenne posterieure principale",
    ylab = "Moyenne posterieure complete",
    main = if (
      comparison_type == "legere_vs_complete"
    ) "Tau : legere contre complete" else
      "Tau : complete contre complete"
  )
  graphics::abline(0, 1, lty = 2, col = "#333333")
  graphics::grid(col = "#E5E7EB")
}
close_figure()

message(
  "Rapport technique : tables, macros et figures generees dans ",
  report_dir
)
