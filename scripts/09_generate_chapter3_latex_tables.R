#!/usr/bin/env Rscript

# Génère les tableaux LaTeX du rapport du chapitre 3 à partir des CSV finaux.
# Cette étape est purement documentaire : elle ne réestime aucun modèle et
# ne recalcule aucun poids de stacking.

options(stringsAsFactors = FALSE, scipen = 999)

root_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
results_dir <- file.path(root_dir, "results", "full", "chapter3_results")
metrics_dir <- file.path(root_dir, "results", "full", "metrics")
weights_dir <- file.path(root_dir, "results", "full", "weights")
diagnostics_dir <- file.path(root_dir, "results", "full", "diagnostics")
tables_dir <- file.path(root_dir, "output", "pdf", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(path) {
  read.csv(path, check.names = FALSE, fileEncoding = "UTF-8")
}

latex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

fmt_num <- function(x, digits = 3L) {
  x <- suppressWarnings(as.numeric(x))
  out <- ifelse(
    is.na(x),
    "",
    ifelse(
      abs(x) > 0 & abs(x) < 10^(-digits),
      formatC(x, format = "e", digits = 2),
      formatC(x, format = "f", digits = digits, big.mark = " ", decimal.mark = ",")
    )
  )
  out
}

fmt_int <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", formatC(round(x), format = "d", big.mark = " "))
}

fmt_pct <- function(x, digits = 1L) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", paste0(fmt_num(100 * x, digits), "%"))
}

method_label <- function(x) {
  labels <- c(
    lc = "LC",
    rh = "RH",
    apc = "APC",
    cbd = "CBD",
    m6 = "M6",
    stacking_global = "Stacking global",
    stacking_contextual = "Contextuel non régularisé",
    stacking_hierarchical = "Contextuel hiérarchique"
  )
  unname(ifelse(x %in% names(labels), labels[x], x))
}

period_label <- function(x) {
  labels <- c(
    test_initial_2016_2020 = "2016--2020",
    extension_2021_2024 = "2021--2024",
    test_extended_2016_2024 = "2016--2024",
    without_2020_2022 = "2016--2019 et 2023--2024"
  )
  unname(ifelse(x %in% names(labels), labels[x], x))
}

write_longtable <- function(df, filename, caption, label, align = NULL,
                            font_size = "\\scriptsize") {
  if (is.null(align)) {
    align <- paste0("l", paste(rep("r", ncol(df) - 1L), collapse = ""))
  }
  headers <- latex_escape(names(df))
  body <- apply(df, 1L, function(row) {
    paste0(paste(latex_escape(row), collapse = " & "), " \\\\")
  })
  head_line <- paste0(paste(headers, collapse = " & "), " \\\\")
  lines <- c(
    font_size,
    paste0("\\begin{longtable}{", align, "}"),
    paste0("\\caption{", caption, "}\\label{", label, "}\\\\"),
    "\\toprule",
    head_line,
    "\\midrule",
    "\\endfirsthead",
    paste0("\\multicolumn{", ncol(df), "}{c}{\\tablename\\ \\thetable\\ -- suite}\\\\"),
    "\\toprule",
    head_line,
    "\\midrule",
    "\\endhead",
    paste0("\\midrule\\multicolumn{", ncol(df), "}{r}{Suite page suivante}\\\\"),
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot",
    body,
    "\\end{longtable}",
    "\\normalsize"
  )
  writeLines(lines, file.path(tables_dir, filename), useBytes = TRUE)
}

# 1. Protocole
protocol <- read_csv(file.path(results_dir, "01_protocol.csv"))
protocol_labels <- c(
  population = "Population",
  ages = "Âges",
  periode_disponible = "Période disponible",
  entrainement_principal = "Entraînement principal",
  validation_lfo = "Années cibles de validation LFO",
  test_initial = "Test initial",
  extension = "Extension temporelle",
  test_etendu = "Test étendu",
  origines_lfo = "Origines LFO",
  horizons_lfo = "Horizons LFO",
  horizon_test_maximal = "Horizon maximal du test",
  tirages_predictifs = "Tirages prédictifs",
  chaines_stan = "Chaînes Stan",
  warmup_par_chaine = "Warm-up par chaîne",
  iterations_conservees_par_chaine = "Itérations conservées par chaîne",
  iterations_conservees_totales = "Tirages postérieurs totaux",
  graine = "Graine aléatoire",
  modeles = "Modèles",
  origine_test = "Origine du test",
  type_test = "Type de test",
  convention_logs = "Convention LogS"
)
protocol$element <- unname(protocol_labels[protocol$element])
protocol$value <- gsub(",", ", ", protocol$value, fixed = TRUE)
names(protocol) <- c("Élément", "Valeur")
write_longtable(
  protocol, "protocol.tex", "Protocole empirique définitif.",
  "tab:protocol", align = "p{0.34\\textwidth}p{0.60\\textwidth}", font_size = "\\small"
)

# 2. Résumé des données
data_summary <- read_csv(file.path(results_dir, "02_data_summary.csv"))[1, ]
data_table <- data.frame(
  Indicateur = c(
    "Population", "Période", "Nombre d'années", "Âges", "Nombre d'âges",
    "Cellules âge--année", "Décès (min. / moyenne / max.)",
    "Exposition (min. / moyenne / max.)",
    "Taux brut (min. / moyenne / max.)", "Valeurs manquantes",
    "Expositions non positives", "Décès négatifs", "Décès fractionnaires",
    "Doublons", "Grille complète"
  ),
  Valeur = c(
    data_summary$population,
    paste0(data_summary$first_year, "--", data_summary$last_year),
    fmt_int(data_summary$years),
    paste0(data_summary$min_age, "--", data_summary$max_age),
    fmt_int(data_summary$ages),
    fmt_int(data_summary$cells),
    paste(fmt_int(data_summary$deaths_min), fmt_num(data_summary$deaths_mean, 2),
          fmt_int(data_summary$deaths_max), sep = " / "),
    paste(fmt_num(data_summary$exposure_min, 2), fmt_num(data_summary$exposure_mean, 2),
          fmt_num(data_summary$exposure_max, 2), sep = " / "),
    paste(fmt_num(data_summary$crude_rate_min, 6), fmt_num(data_summary$crude_rate_mean, 6),
          fmt_num(data_summary$crude_rate_max, 6), sep = " / "),
    fmt_int(as.numeric(data_summary$missing_deaths) + as.numeric(data_summary$missing_exposure)),
    fmt_int(data_summary$nonpositive_exposure),
    fmt_int(data_summary$negative_deaths),
    fmt_int(data_summary$fractional_deaths),
    fmt_int(data_summary$duplicate_cells),
    ifelse(data_summary$complete_grid == "TRUE", "Oui", "Non")
  )
)
write_longtable(
  data_table, "data_summary.tex", "Résumé et contrôles de qualité des données.",
  "tab:data-summary", align = "p{0.47\\textwidth}p{0.47\\textwidth}", font_size = "\\small"
)

# 3. Diagnostics MCMC
diag <- read_csv(file.path(results_dir, "03_model_diagnostics.csv"))
diag_table <- data.frame(
  Modèle = diag$model,
  `Temps (min)` = fmt_num(diag$elapsed_minutes, 2),
  `Rhat max.` = fmt_num(diag$max_rhat, 3),
  `ESS bulk min.` = fmt_int(diag$min_ess_bulk),
  `ESS tail min.` = fmt_int(diag$min_ess_tail),
  Divergences = fmt_int(diag$divergences),
  Treedepth = fmt_int(diag$max_treedepth_hits),
  `E-BFMI min.` = fmt_num(diag$min_ebfmi, 3),
  Statut = diag$status,
  check.names = FALSE
)
write_longtable(
  diag_table, "model_diagnostics.tex", "Diagnostics des ajustements finaux sur 1970--2015.",
  "tab:model-diagnostics", align = "lrrrrrrrl", font_size = "\\tiny"
)

# 4. Paramètres principaux
params <- read_csv(file.path(results_dir, "model_parameter_summary.csv"))
params_table <- data.frame(
  Modèle = params$model,
  Paramètre = params$parameter,
  Moyenne = fmt_num(params$mean, 5),
  Médiane = fmt_num(params$median, 5),
  Écart.type = fmt_num(params$sd, 5),
  `Q2,5 %` = fmt_num(params$q025, 5),
  `Q97,5 %` = fmt_num(params$q975, 5),
  check.names = FALSE
)
write_longtable(
  params_table, "model_parameters.tex", "Résumé postérieur des principaux paramètres.",
  "tab:model-parameters", align = "llrrrrr"
)

# 5. Validation LFO
meta <- read_csv(file.path(results_dir, "meta_validation.csv"))
lfo_summary <- data.frame(
  Information = c(
    "Origines distinctes", "Origines", "Années cibles", "Horizons",
    "Âges", "Cellules de validation hors modèle", "Log-densités (5 modèles)",
    "Valeurs manquantes", "Valeurs non finies"
  ),
  Résultat = c(
    fmt_int(length(unique(meta$origin))),
    paste(range(meta$origin), collapse = "--"),
    paste(range(meta$target_year), collapse = "--"),
    paste(sort(unique(meta$horizon)), collapse = ", "),
    fmt_int(length(unique(meta$age))),
    fmt_int(nrow(meta) / length(unique(meta$model))),
    fmt_int(nrow(meta)),
    fmt_int(sum(is.na(meta$log_pred_density))),
    fmt_int(sum(!is.finite(meta$log_pred_density)))
  )
)
write_longtable(
  lfo_summary, "lfo_summary.tex", "Structure et intégrité du méta-jeu LFO.",
  "tab:lfo-summary", align = "p{0.58\\textwidth}p{0.36\\textwidth}", font_size = "\\small"
)

lfo_models <- read_csv(file.path(results_dir, "lfo_model_summary.csv"))
lfo_model_table <- data.frame(
  Modèle = lfo_models$model,
  `Log-densité pondérée` = fmt_num(lfo_models$weighted_mean_log_density, 4),
  `Log-densité non pondérée` = fmt_num(lfo_models$unweighted_mean_log_density, 4),
  Minimum = fmt_num(lfo_models$minimum, 3),
  Maximum = fmt_num(lfo_models$maximum, 3),
  check.names = FALSE
)
write_longtable(
  lfo_model_table, "lfo_model_summary.tex", "Log-densités prédictives dans la validation LFO (plus élevée = meilleure).",
  "tab:lfo-model-summary", align = "lrrrr"
)

# 6. Poids et comparaison des méthodes
global_weights <- read_csv(file.path(results_dir, "04_global_weights.csv"))
global_table <- data.frame(
  Modèle = global_weights$model,
  Poids = fmt_num(global_weights$weight, 6)
)
write_longtable(
  global_table, "global_weights.tex", "Poids du stacking global.",
  "tab:global-weights", align = "lr", font_size = "\\small"
)

stack_comp <- read_csv(file.path(results_dir, "stacking_comparison.csv"))
stack_names <- c(
  stacking_global = "Global",
  stacking_contextual = "Contextuel non régularisé",
  stacking_hierarchical = "Contextuel hiérarchique"
)
stack_table <- data.frame(
  Méthode = unname(stack_names[stack_comp$method]),
  Critère = fmt_num(stack_comp$validation_criterion, 3),
  Paramètres = fmt_int(stack_comp$parameter_count),
  `Poids min.` = fmt_num(stack_comp$minimum_weight, 6),
  `Poids max.` = fmt_num(stack_comp$maximum_weight, 6),
  Stabilité = stack_comp$stability,
  check.names = FALSE
)
write_longtable(
  stack_table, "stacking_comparison.tex", "Comparaison des trois formes de stacking sur la validation.",
  "tab:stacking-comparison", align = "p{0.22\\textwidth}rrrrp{0.34\\textwidth}"
)

weight_h <- read_csv(file.path(weights_dir, "stacking_weights_by_horizon_1_9.csv"))
weight_h_table <- data.frame(
  Méthode = method_label(weight_h$method),
  Horizon = fmt_int(weight_h$horizon),
  Modèle = toupper(weight_h$model),
  Moyenne = fmt_num(weight_h$weight_mean_over_ages, 4),
  Minimum = fmt_num(weight_h$weight_min_over_ages, 4),
  Maximum = fmt_num(weight_h$weight_max_over_ages, 4)
)
write_longtable(
  weight_h_table, "weights_by_horizon.tex",
  "Poids aux horizons 1 à 9, résumés sur les 41 âges.",
  "tab:weights-horizon", align = "llrrrr"
)

coef <- read_csv(file.path(results_dir, "07_weight_coefficients.csv"))
coef_table <- data.frame(
  Méthode = method_label(coef$method),
  Modèle = coef$model,
  Terme = coef$term,
  Moyenne = fmt_num(coef$mean, 4),
  Médiane = fmt_num(coef$median, 4),
  Écart.type = fmt_num(coef$sd, 4),
  `Q2,5 %` = fmt_num(coef$q025, 4),
  `Q97,5 %` = fmt_num(coef$q975, 4),
  check.names = FALSE
)
write_longtable(
  coef_table, "weight_coefficients.tex", "Coefficients des stackings contextuels.",
  "tab:weight-coefficients", align = "lllrrrrr"
)

# 7. Performances globales et par période
overall <- read_csv(file.path(results_dir, "08_global_metrics.csv"))
overall_table <- data.frame(
  Méthode = method_label(overall$method),
  LogS = fmt_num(overall$logs, 3),
  CRPS = fmt_num(overall$crps, 6),
  `MAE décès` = fmt_num(overall$mae_deaths, 1),
  `Couv. 80 %` = fmt_pct(overall$coverage80, 1),
  `Larg. 80 %` = fmt_num(overall$width80, 5),
  `Couv. 95 %` = fmt_pct(overall$coverage95, 1),
  `Larg. 95 %` = fmt_num(overall$width95, 5),
  check.names = FALSE
)
write_longtable(
  overall_table, "overall_metrics.tex", "Performances sur le test étendu 2016--2024.",
  "tab:overall-metrics", align = "lrrrrrrr"
)

period_metrics <- read_csv(file.path(metrics_dir, "test_metrics_overall_by_period.csv"))
period_order <- c(
  test_initial_2016_2020 = 1L,
  extension_2021_2024 = 2L,
  test_extended_2016_2024 = 3L
)
method_order <- c(
  lc = 1L, rh = 2L, apc = 3L, cbd = 4L, m6 = 5L,
  stacking_global = 6L, stacking_contextual = 7L, stacking_hierarchical = 8L
)
period_metrics <- period_metrics[order(period_order[period_metrics$period],
                                       method_order[period_metrics$method]), ]
period_table <- data.frame(
  Période = period_label(period_metrics$period),
  Méthode = method_label(period_metrics$method),
  LogS = fmt_num(period_metrics$logs, 3),
  CRPS = fmt_num(period_metrics$crps, 6),
  `MAE décès` = fmt_num(period_metrics$mae_deaths, 1),
  `Couv. 80 %` = fmt_pct(period_metrics$coverage80, 1),
  `Couv. 95 %` = fmt_pct(period_metrics$coverage95, 1),
  check.names = FALSE
)
write_longtable(
  period_table, "metrics_by_period.tex",
  "Performances des huit méthodes dans les trois blocs temporels.",
  "tab:metrics-period", align = "llrrrrr"
)

# 8. Détails par âge et horizon
by_age <- read_csv(file.path(results_dir, "09_metrics_by_age.csv"))
by_age <- by_age[order(method_order[by_age$method], by_age$age), ]
by_age_table <- data.frame(
  Méthode = method_label(by_age$method),
  Âge = fmt_int(by_age$age),
  LogS = fmt_num(by_age$logs, 3),
  CRPS = fmt_num(by_age$crps, 6),
  `MAE décès` = fmt_num(by_age$mae_deaths, 1),
  `Couv. 80 %` = fmt_pct(by_age$coverage80, 1),
  `Couv. 95 %` = fmt_pct(by_age$coverage95, 1),
  check.names = FALSE
)
write_longtable(
  by_age_table, "metrics_by_age.tex",
  "Performances détaillées par âge sur 2016--2024.",
  "tab:metrics-age", align = "lrrrrrr"
)

by_horizon <- read_csv(file.path(results_dir, "10_metrics_by_horizon.csv"))
by_horizon <- by_horizon[order(method_order[by_horizon$method], by_horizon$horizon), ]
by_horizon_table <- data.frame(
  Méthode = method_label(by_horizon$method),
  Horizon = fmt_int(by_horizon$horizon),
  LogS = fmt_num(by_horizon$logs, 3),
  CRPS = fmt_num(by_horizon$crps, 6),
  `MAE décès` = fmt_num(by_horizon$mae_deaths, 1),
  `Couv. 80 %` = fmt_pct(by_horizon$coverage80, 1),
  `Couv. 95 %` = fmt_pct(by_horizon$coverage95, 1),
  check.names = FALSE
)
write_longtable(
  by_horizon_table, "metrics_by_horizon.tex",
  "Performances détaillées par horizon sur 2016--2024.",
  "tab:metrics-horizon", align = "lrrrrrr"
)

# 9. Incertitude, intervalles et sensibilité
boot <- read_csv(file.path(results_dir, "performance_differences_bootstrap.csv"))
metric_names <- c(
  logs = "LogS",
  crps = "CRPS",
  absolute_error_deaths = "MAE décès"
)
comparison_names <- c(
  contextual_minus_global = "Contextuel -- global",
  hierarchical_minus_global = "Hiérarchique -- global",
  hierarchical_minus_contextual = "Hiérarchique -- contextuel"
)
boot_table <- data.frame(
  Comparaison = unname(comparison_names[boot$comparison]),
  Métrique = unname(metric_names[boot$metric]),
  Différence = fmt_num(boot$mean_difference, 6),
  `IC 2,5 %` = fmt_num(boot$ci025, 6),
  `IC 97,5 %` = fmt_num(boot$ci975, 6),
  `Cellules gagnées` = fmt_pct(boot$proportion_first_better, 1),
  Répétitions = fmt_int(boot$bootstrap_repetitions),
  check.names = FALSE
)
write_longtable(
  boot_table, "bootstrap_differences.tex",
  "Différences de performance et bootstrap par année cible. Une différence négative favorise la première méthode.",
  "tab:bootstrap", align = "llrrrrr"
)

interval <- read_csv(file.path(results_dir, "11_interval_metrics.csv"))
interval_table <- data.frame(
  Méthode = method_label(interval$method),
  `Couv. 80 %` = fmt_pct(interval$coverage80, 1),
  `Larg. 80 %` = fmt_num(interval$width80, 5),
  `Couv. 95 %` = fmt_pct(interval$coverage95, 1),
  `Larg. 95 %` = fmt_num(interval$width95, 5),
  check.names = FALSE
)
write_longtable(
  interval_table, "interval_metrics.tex",
  "Couverture et largeur des intervalles prédictifs sur 2016--2024.",
  "tab:intervals", align = "lrrrr"
)

sensitivity <- read_csv(file.path(results_dir, "14_sensitivity_results.csv"))
sensitivity <- sensitivity[order(
  match(sensitivity$period, c("test_extended_2016_2024", "without_2020_2022")),
  method_order[sensitivity$method]
), ]
sensitivity_table <- data.frame(
  Période = period_label(sensitivity$period),
  Méthode = method_label(sensitivity$method),
  LogS = fmt_num(sensitivity$logs, 3),
  CRPS = fmt_num(sensitivity$crps, 6),
  `MAE décès` = fmt_num(sensitivity$mae_deaths, 1),
  `Couv. 80 %` = fmt_pct(sensitivity$coverage80, 1),
  `Couv. 95 %` = fmt_pct(sensitivity$coverage95, 1),
  check.names = FALSE
)
write_longtable(
  sensitivity_table, "sensitivity.tex",
  "Sensibilité à l'exclusion des années 2020--2022.",
  "tab:sensitivity", align = "llrrrrr"
)

# 10. Quantités actuarielles
survival <- read_csv(file.path(results_dir, "12_survival_probabilities.csv"))
survival <- survival[
  survival$method %in% c("stacking_global", "stacking_contextual", "stacking_hierarchical") &
    survival$horizon %in% c(1, 5, 10, 15, 20, 25),
]
survival <- survival[order(survival$initial_age, survival$horizon,
                           method_order[survival$method]), ]
survival_table <- data.frame(
  Âge = fmt_int(survival$initial_age),
  Horizon = fmt_int(survival$horizon),
  Méthode = method_label(survival$method),
  Moyenne = fmt_num(survival$mean, 5),
  Médiane = fmt_num(survival$median, 5),
  `Q2,5 %` = fmt_num(survival$q025, 5),
  `Q97,5 %` = fmt_num(survival$q975, 5),
  check.names = FALSE
)
write_longtable(
  survival_table, "survival_selected.tex",
  "Probabilités de survie sélectionnées pour les trois formes de stacking.",
  "tab:survival-selected", align = "rrlrrrr"
)

annuity <- read_csv(file.path(results_dir, "13_annuity_values.csv"))
annuity <- annuity[
  annuity$method %in% c("stacking_global", "stacking_contextual", "stacking_hierarchical"),
]
annuity <- annuity[order(annuity$initial_age, method_order[annuity$method]), ]
annuity_duration <- if ("horizon" %in% names(annuity)) {
  annuity$horizon
} else {
  annuity$duration
}
annuity_table <- data.frame(
  Âge = fmt_int(annuity$initial_age),
  Durée = fmt_int(annuity_duration),
  Méthode = method_label(annuity$method),
  Moyenne = fmt_num(annuity$mean, 4),
  Médiane = fmt_num(annuity$median, 4),
  Écart.type = fmt_num(annuity$sd, 4),
  `Q2,5 %` = fmt_num(annuity$q025, 4),
  `Q97,5 %` = fmt_num(annuity$q975, 4),
  `Écart au global` = fmt_num(annuity$difference_vs_global, 4),
  `Écart relatif` = paste0(fmt_num(annuity$relative_difference_vs_global_pct, 2), "%"),
  check.names = FALSE
)
write_longtable(
  annuity_table, "annuity_values.tex",
  "Valeurs de rentes temporaires immédiates à terme échu, taux annuel de 2\\%.",
  "tab:annuity", align = "rrlrrrrrrr"
)

# 11. Audit de non-régression
audit <- read_csv(file.path(diagnostics_dir, "extension_2024_audit_checks.csv"))
audit_table <- data.frame(
  Contrôle = gsub("_", " ", audit$check),
  Statut = ifelse(audit$pass == "TRUE", "Réussi", "Échec"),
  Détails = gsub(",", ", ", audit$details, fixed = TRUE)
)
write_longtable(
  audit_table, "extension_audit.tex",
  "Contrôles de non-régression et d'absence de fuite d'information.",
  "tab:extension-audit",
  align = "p{0.33\\textwidth}p{0.10\\textwidth}p{0.49\\textwidth}",
  font_size = "\\scriptsize"
)

message("Tableaux LaTeX écrits dans : ", tables_dir)
