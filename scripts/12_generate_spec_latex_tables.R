#!/usr/bin/env Rscript

# Tableaux LaTeX supplémentaires exigés par le cahier des charges.
# Ce script est documentaire : il ne calcule ni modèle, ni poids, ni score.

options(stringsAsFactors = FALSE, scipen = 999)
root_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
results_dir <- file.path(root_dir, "results", "full", "chapter3_results")
simulation_dir <- file.path(root_dir, "results", "full", "simulation")
weights_dir <- file.path(root_dir, "results", "full", "weights")
tables_dir <- file.path(root_dir, "output", "pdf", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

read_csv <- function(filename, directory = results_dir) {
  utils::read.csv(
    file.path(directory, filename),
    check.names = FALSE, fileEncoding = "UTF-8",
    stringsAsFactors = FALSE
  )
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
  ifelse(
    is.na(x), "",
    ifelse(
      abs(x) > 0 & abs(x) < 10^(-digits),
      formatC(x, format = "e", digits = 2),
      formatC(
        x, format = "f", digits = digits,
        big.mark = " ", decimal.mark = ","
      )
    )
  )
}

fmt_int <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", formatC(round(x), format = "d", big.mark = " "))
}

fmt_pct <- function(x, digits = 1L) {
  x <- suppressWarnings(as.numeric(x))
  # write_latex_table() échappe ensuite le signe %, il ne faut donc pas
  # injecter ici une commande LaTeX qui serait échappée une seconde fois.
  ifelse(is.na(x), "", paste0(fmt_num(100 * x, digits), "%"))
}

method_label <- function(x) {
  labels <- c(
    lc = "LC", rh = "RH", apc = "APC", cbd = "CBD", m6 = "M6",
    pseudo_bma = "Pseudo-BMA",
    stacking_global = "Stacking global",
    stacking_contextual = "Contextuel non régularisé",
    stacking_hierarchical = "Contextuel hiérarchique"
  )
  unname(ifelse(x %in% names(labels), labels[x], gsub("_", " ", x)))
}

actuarial_method_label <- function(method, weight_mode) {
  result <- method_label(method)
  result[
    method == "stacking_hierarchical" &
      weight_mode == "posterior_mean_weights"
  ] <- "Hiérarchique - poids moyens"
  result[
    method == "stacking_hierarchical" &
      weight_mode == "propagated_weights"
  ] <- "Hiérarchique - poids propagés"
  result
}

actuarial_comparison_label <- function(comparison) {
  labels <- c(
    hierarchical_mean_minus_global = "Hiér. moyen - global",
    hierarchical_propagated_minus_global = "Hiér. propagé - global",
    hierarchical_mean_minus_contextual = "Hiér. moyen - contextuel",
    hierarchical_propagated_minus_contextual =
      "Hiér. propagé - contextuel",
    hierarchical_propagated_minus_hierarchical_mean =
      "Hiér. propagé - hiér. moyen"
  )
  unname(ifelse(
    comparison %in% names(labels),
    labels[comparison],
    gsub("_", " ", comparison)
  ))
}

write_longtable <- function(
    df, filename, caption, label, align = NULL,
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
    "\\toprule", head_line, "\\midrule", "\\endfirsthead",
    paste0(
      "\\multicolumn{", ncol(df),
      "}{c}{\\tablename\\ \\thetable\\ -- suite}\\\\"
    ),
    "\\toprule", head_line, "\\midrule", "\\endhead",
    paste0(
      "\\midrule\\multicolumn{", ncol(df),
      "}{r}{Suite page suivante}\\\\"
    ),
    "\\endfoot", "\\bottomrule", "\\endlastfoot",
    body, "\\end{longtable}", "\\normalsize"
  )
  writeLines(
    lines, file.path(tables_dir, filename),
    useBytes = TRUE
  )
}

# Données et protocole.
source_data <- read_csv("16_data_source_and_definitions.csv")
names(source_data) <- c("Élément", "Valeur")
# Les identifiants techniques très longs deviennent lisibles et peuvent
# revenir à la ligne normalement dans la première colonne.
source_data$Élément <- gsub("_", " ", source_data$Élément, fixed = TRUE)
write_longtable(
  source_data, "source_data_definitions.tex",
  "Source, définitions et contrôles numériques des données.",
  "tab:source-definitions",
  align = "p{0.33\\textwidth}p{0.61\\textwidth}",
  font_size = "\\small"
)

software <- read_csv("24_software_and_sampling_settings.csv")
names(software) <- c("Composant", "Valeur")
write_longtable(
  software, "software_sampling.tex",
  "Logiciels et réglages d'échantillonnage.",
  "tab:software-sampling",
  align = "p{0.42\\textwidth}p{0.52\\textwidth}",
  font_size = "\\small"
)

runtime <- read_csv("23_fit_runtime_summary.csv")
write_longtable(
  data.frame(
    Étape = runtime$stage,
    Ajustements = fmt_int(runtime$fits),
    `Heures (proxy parallèle)` =
      fmt_num(runtime$elapsed_hours_parallel_wall_proxy, 2),
    check.names = FALSE
  ),
  "runtime_summary.tex",
  "Temps de calcul reconstitués depuis les objets Stan.",
  "tab:runtime", align = "lrr", font_size = "\\small"
)

# Diagnostics MCMC.
lfo_summary <- read_csv("18_lfo_diagnostics_summary.csv")
write_longtable(
  data.frame(
    Modèle = toupper(lfo_summary$model),
    Ajustements = fmt_int(lfo_summary$fits),
    Validés = fmt_int(lfo_summary$fits_passed),
    `Rhat max.` = fmt_num(lfo_summary$worst_rhat, 4),
    `ESS bulk min.` = fmt_int(lfo_summary$minimum_ess_bulk),
    `ESS tail min.` = fmt_int(lfo_summary$minimum_ess_tail),
    Divergences = fmt_int(lfo_summary$divergences),
    Treedepth = fmt_int(lfo_summary$max_treedepth_hits),
    check.names = FALSE
  ),
  "lfo_diagnostics_summary.tex",
  "Diagnostics des 75 ajustements LFO, résumés par modèle.",
  "tab:lfo-diag-summary", align = "lrrrrrrr", font_size = "\\tiny"
)

lfo_all <- read_csv("17_lfo_diagnostics_all_75.csv")
write_longtable(
  data.frame(
    Origine = fmt_int(lfo_all$origin),
    Modèle = toupper(lfo_all$model),
    `Rhat max.` = fmt_num(lfo_all$max_rhat, 4),
    `ESS bulk min.` = fmt_int(lfo_all$min_ess_bulk),
    `ESS tail min.` = fmt_int(lfo_all$min_ess_tail),
    Divergences = fmt_int(lfo_all$divergences),
    Treedepth = fmt_int(lfo_all$max_treedepth),
    Statut = ifelse(lfo_all$pass == "TRUE", "Valide", "Échec"),
    check.names = FALSE
  ),
  "lfo_diagnostics_all_75.tex",
  "Diagnostics détaillés des 75 ajustements LFO.",
  "tab:lfo-diag-all", align = "rlrrrrrl", font_size = "\\tiny"
)

hier_sampler <- read_csv("21_hierarchical_sampler_diagnostics.csv")
write_longtable(
  data.frame(
    Chaîne = fmt_int(hier_sampler$chain_id),
    Divergences = fmt_int(hier_sampler$num_divergent),
    Treedepth = fmt_int(hier_sampler$num_max_treedepth),
    `E-BFMI` = fmt_num(hier_sampler$ebfmi, 4),
    check.names = FALSE
  ),
  "hierarchical_sampler.tex",
  "Diagnostics par chaîne du stacking hiérarchique.",
  "tab:hier-sampler", align = "rrrr", font_size = "\\small"
)

hier_params <- read_csv("20_hierarchical_parameter_diagnostics.csv")
write_longtable(
  data.frame(
    Paramètre = hier_params$variable,
    Moyenne = fmt_num(hier_params$mean, 5),
    Écart.type = fmt_num(hier_params$sd, 5),
    Rhat = fmt_num(hier_params$rhat, 4),
    `ESS bulk` = fmt_int(hier_params$ess_bulk),
    `ESS tail` = fmt_int(hier_params$ess_tail),
    check.names = FALSE
  ),
  "hierarchical_parameter_diagnostics.tex",
  "Paramètres et diagnostics du stacking hiérarchique.",
  "tab:hier-parameters", align = "lrrrrr", font_size = "\\tiny"
)

# Méta-LFO et optimisation.
year_weights <- read_csv("26_lfo_year_weight_totals.csv")
write_longtable(
  data.frame(
    Année = fmt_int(year_weights$year),
    `Somme tous âges` = fmt_num(year_weights$sum_omega_all_ages, 3),
    `Somme par âge min.` =
      fmt_num(year_weights$sum_omega_per_age_min, 3),
    `Somme par âge moy.` =
      fmt_num(year_weights$sum_omega_per_age_mean, 3),
    `Somme par âge max.` =
      fmt_num(year_weights$sum_omega_per_age_max, 3),
    check.names = FALSE
  ),
  "lfo_year_weight_totals.tex",
  "Somme des poids de multiplicité par année cible.",
  "tab:lfo-year-weights", align = "rrrrr", font_size = "\\small"
)

weighted <- read_csv("27_weighted_vs_unweighted_lfo.csv")
write_longtable(
  data.frame(
    Modèle = toupper(weighted$model),
    `Log-pred. pond.` =
      fmt_num(weighted$weighted_mean_log_predictive, 4),
    `Rang pond.` = fmt_int(weighted$weighted_rank),
    `Log-pred. non pond.` =
      fmt_num(weighted$unweighted_mean_log_predictive, 4),
    `Rang non pond.` = fmt_int(weighted$unweighted_rank),
    `Éc.-type annuel` =
      fmt_num(weighted$sd_year_mean_log_predictive, 4),
    check.names = FALSE
  ),
  "weighted_unweighted_lfo.tex",
  "Validation LFO pondérée et non pondérée.",
  "tab:weighted-unweighted", align = "lrrrrr", font_size = "\\scriptsize"
)

optimization <- read_csv("29_stacking_optimization_audit.csv")
write_longtable(
  data.frame(
    Méthode = method_label(optimization$method),
    Algorithme = optimization$algorithm,
    Départs = fmt_int(optimization$starts),
    Convergence = fmt_int(optimization$convergence_code),
    Évaluations = fmt_int(optimization$iterations_function),
    Objectif = fmt_num(optimization$objective, 4),
    `Étendue départs` =
      fmt_num(optimization$objective_range_across_starts, 6),
    `Gradient max.` =
      fmt_num(optimization$max_absolute_gradient, 6),
    check.names = FALSE
  ),
  "stacking_optimization_audit.tex",
  "Audit des optimisations de stacking.",
  "tab:optimization-audit", align = "llrrrrrr", font_size = "\\tiny"
)

extreme <- read_csv("30_nonregularized_extreme_weights.csv")
write_longtable(
  data.frame(
    Modèle = toupper(extreme$model),
    `Poids <1 %` = fmt_pct(extreme$proportion_below_001, 1),
    `Poids >99 %` = fmt_pct(extreme$proportion_above_099, 1),
    Minimum = fmt_num(extreme$minimum, 5),
    Médiane = fmt_num(extreme$median, 5),
    Maximum = fmt_num(extreme$maximum, 5),
    check.names = FALSE
  ),
  "nonregularized_extreme_weights.tex",
  "Fréquence des poids extrêmes du stacking non régularisé.",
  "tab:extreme-weights", align = "lrrrrr", font_size = "\\small"
)

validation <- read_csv("31_validation_method_comparison.csv")
write_longtable(
  data.frame(
    Méthode = method_label(validation$method),
    `Nég. LogS pondéré` =
      fmt_num(validation$weighted_negative_log_score, 3),
    `LogS moyen pondéré` =
      fmt_num(validation$weighted_mean_log_score, 4),
    `LogS moyen non pondéré` =
      fmt_num(validation$unweighted_mean_log_score, 4),
    Rang = fmt_int(validation$weighted_rank),
    check.names = FALSE
  ),
  "validation_method_comparison.tex",
  "Comparaison des agrégateurs sur le méta-jeu LFO.",
  "tab:validation-comparison", align = "lrrrr", font_size = "\\scriptsize"
)

global_methods <- read_csv("weights_global.csv", weights_dir)
write_longtable(
  data.frame(
    Méthode = method_label(global_methods$method),
    Modèle = toupper(global_methods$model),
    Poids = fmt_num(global_methods$weight, 8),
    check.names = FALSE
  ),
  "global_bma_pseudobma_weights.tex",
  "Poids BMA, pseudo-BMA et stacking global.",
  "tab:global-method-weights", align = "llr", font_size = "\\small"
)

# Ajustement et performances finales.
ppc <- read_csv("36_ppc_residual_discrepancies.csv")
write_longtable(
  data.frame(
    Modèle = ppc$model,
    `Résidu moyen` = fmt_num(ppc$mean_residual, 3),
    `RMS résiduel` = fmt_num(ppc$rms_residual, 3),
    `|résidu|>2` = fmt_pct(ppc$proportion_abs_residual_gt_2, 1),
    `Taux dans IC param. 95 %` =
      fmt_pct(ppc$proportion_observed_rate_in_parameter_interval, 1),
    check.names = FALSE
  ),
  "ppc_residual_summary.tex",
  "Diagnostics résiduels et contrôles prédictifs conditionnels.",
  "tab:ppc-residual", align = "lrrrr", font_size = "\\small"
)

test_periods <- read_csv("37_test_metrics_8_methods_3_periods.csv")
write_longtable(
  data.frame(
    Période = test_periods$period,
    Méthode = method_label(test_periods$method),
    LogS = fmt_num(test_periods$logs, 3),
    CRPS = fmt_num(test_periods$crps, 6),
    MAE = fmt_num(test_periods$mae_deaths, 1),
    `Couv. 80 %` = fmt_pct(test_periods$coverage80, 1),
    `Écart 80 (pts)` =
      fmt_num(test_periods$coverage80_gap_points, 1),
    `Couv. 95 %` = fmt_pct(test_periods$coverage95, 1),
    `Écart 95 (pts)` =
      fmt_num(test_periods$coverage95_gap_points, 1),
    check.names = FALSE
  ),
  "test_metrics_8_methods_3_periods.tex",
  "Huit méthodes, trois périodes et écarts aux couvertures nominales.",
  "tab:test-three-periods", align = "llrrrrrrr", font_size = "\\tiny"
)

relative <- read_csv("38_hierarchical_relative_test_performance.csv")
write_longtable(
  data.frame(
    Période = relative$period,
    Comparateur = method_label(relative$comparator),
    `Gain LogS` = paste0(fmt_num(relative$logs_reduction_pct, 1), "%"),
    `Gain CRPS` = paste0(fmt_num(relative$crps_reduction_pct, 1), "%"),
    `Gain MAE` = paste0(fmt_num(relative$mae_reduction_pct, 1), "%"),
    `Écart 80` =
      fmt_num(relative$coverage80_difference_points, 1),
    `Écart 95` =
      fmt_num(relative$coverage95_difference_points, 1),
    check.names = FALSE
  ),
  "hierarchical_relative_performance.tex",
  "Performance relative du stacking hiérarchique.",
  "tab:relative-performance",
  align = "p{0.18\\textwidth}p{0.22\\textwidth}rrrrr",
  font_size = "\\tiny"
)

coverage <- read_csv("39_coverage_counts_by_age_out_of_9.csv")
coverage <- coverage[coverage$method == "stacking_hierarchical", ]
write_longtable(
  data.frame(
    Âge = fmt_int(coverage$age),
    `Couverts 80 / 9` =
      fmt_int(coverage$covered80_years_out_of_9),
    `Écart 80 (pts)` =
      fmt_num(coverage$coverage80_gap_points, 1),
    `Largeur 80 moy.` = fmt_num(coverage$width80, 5),
    `Couverts 95 / 9` =
      fmt_int(coverage$covered95_years_out_of_9),
    `Écart 95 (pts)` =
      fmt_num(coverage$coverage95_gap_points, 1),
    `Largeur 95 moy.` = fmt_num(coverage$width95, 5),
    check.names = FALSE
  ),
  "coverage_counts_hierarchical.tex",
  "Nombre d'années couvertes sur neuf, par âge, stacking hiérarchique.",
  "tab:coverage-counts", align = "rrrrrrr", font_size = "\\scriptsize"
)

pandemic <- read_csv("41_pandemic_sensitivity_6_blocks.csv")
write_longtable(
  data.frame(
    Bloc = pandemic$period,
    Méthode = method_label(pandemic$method),
    LogS = fmt_num(pandemic$logs, 3),
    Rang.LogS = fmt_int(pandemic$rank_logs),
    CRPS = fmt_num(pandemic$crps, 6),
    Rang.CRPS = fmt_int(pandemic$rank_crps),
    MAE = fmt_num(pandemic$mae_deaths, 1),
    Rang.MAE = fmt_int(pandemic$rank_mae_deaths),
    `Couv. 80 %` = fmt_pct(pandemic$coverage80, 1),
    `Couv. 95 %` = fmt_pct(pandemic$coverage95, 1),
    check.names = FALSE
  ),
  "pandemic_sensitivity_6_blocks.tex",
  "Sensibilité pandémique en six blocs et classements.",
  "tab:pandemic-six-blocks", align = "llrrrrrrrr", font_size = "\\tiny"
)

jackknife <- read_csv("42_leave_one_target_year_out_sensitivity.csv")
write_longtable(
  data.frame(
    `Année exclue` = fmt_int(jackknife$excluded_year),
    Comparaison = gsub("_", " ", jackknife$comparison),
    `Diff. LogS` = fmt_num(jackknife$logs_difference, 5),
    `Diff. CRPS` = fmt_num(jackknife$crps_difference, 7),
    `Diff. MAE` = fmt_num(jackknife$mae_difference, 2),
    check.names = FALSE
  ),
  "jackknife_target_year.tex",
  "Sensibilité leave-one-target-year-out.",
  "tab:jackknife", align = "rlrrr", font_size = "\\scriptsize"
)

# Robustesse.
robustness <- read_csv("43_context_and_prior_robustness_metrics.csv")
write_longtable(
  data.frame(
    Variante = gsub("_", " ", robustness$method),
    LogS = fmt_num(robustness$logs, 3),
    CRPS = fmt_num(robustness$crps, 6),
    MAE = fmt_num(robustness$mae_deaths, 1),
    `Couv. 80 %` = fmt_pct(robustness$coverage80, 1),
    `Couv. 95 %` = fmt_pct(robustness$coverage95, 1),
    check.names = FALSE
  ),
  "context_prior_robustness.tex",
  "Ablations du prédicteur et sensibilités de prior.",
  "tab:context-prior-robustness", align = "lrrrrr", font_size = "\\small"
)

prior <- read_csv("44_prior_importance_reweighting_diagnostics.csv")
write_longtable(
  data.frame(
    Variante = prior$sensitivity,
    `ESS d'importance` = fmt_int(prior$importance_ess),
    Tirages = fmt_int(prior$total_draws),
    `Poids normalisé max.` =
      fmt_num(prior$maximum_normalized_weight, 6),
    check.names = FALSE
  ),
  "prior_reweighting_diagnostics.tex",
  "Diagnostics de repondération d'importance des priors.",
  "tab:prior-reweighting", align = "lrrr", font_size = "\\small"
)

# Actuariat.
actuarial <- read_csv("47_actuarial_annuity_all_rules.csv")
actuarial <- actuarial[
  actuarial$rule == "within_lfo" &
    actuarial$method %in%
      c("stacking_global", "stacking_contextual", "stacking_hierarchical") &
    abs(actuarial$discount_rate - 0.02) < 1e-12,
]
write_longtable(
  data.frame(
    Âge = fmt_int(actuarial$initial_age),
    Durée = fmt_int(actuarial$duration),
    Méthode = actuarial_method_label(
      actuarial$method, actuarial$weight_mode
    ),
    Moyenne = fmt_num(actuarial$mean, 5),
    `Écart-type` = fmt_num(actuarial$sd, 5),
    `Q2,5 %` = fmt_num(actuarial$q025, 5),
    `Q97,5 %` = fmt_num(actuarial$q975, 5),
    `Largeur 95 %` = fmt_num(actuarial$width95, 5),
    check.names = FALSE
  ),
  "actuarial_primary_annuity.tex",
  paste0(
    "Rentes principales par agrégation des forces, ",
    "limitées à h<=10, taux 2\\%."
  ),
  "tab:actuarial-primary", align = "rrlrrrrr",
  font_size = "\\scriptsize"
)

paired <- read_csv("48_actuarial_paired_differences.csv")
paired <- paired[
  paired$rule == "within_lfo" &
    abs(paired$discount_rate - 0.02) < 1e-12,
]
write_longtable(
  data.frame(
    Âge = fmt_int(paired$initial_age),
    Comparaison = actuarial_comparison_label(paired$comparison),
    `Diff. moyenne` = fmt_num(paired$mean, 5),
    `Q10 %` = fmt_num(paired$q10, 5),
    `Q90 %` = fmt_num(paired$q90, 5),
    `Q2,5 %` = fmt_num(paired$q025, 5),
    `Q97,5 %` = fmt_num(paired$q975, 5),
    `Pr(diff>0)` = fmt_pct(paired$probability_positive, 1),
    `Pr(diff<0)` = fmt_pct(paired$probability_negative, 1),
    check.names = FALSE
  ),
  "actuarial_paired_differences.tex",
  "Différences actuarielles appariées après agrégation des forces, h<=10.",
  "tab:actuarial-paired", align = "rlrrrrrrr",
  font_size = "\\tiny"
)

actuarial_all <- read_csv("47_actuarial_annuity_all_rules.csv")
actuarial_all <- actuarial_all[
  actuarial_all$method %in%
    c("stacking_global", "stacking_contextual", "stacking_hierarchical"),
]
write_longtable(
  data.frame(
    Analyse = actuarial_all$analysis,
    Règle = actuarial_all$rule,
    Méthode = actuarial_method_label(
      actuarial_all$method, actuarial_all$weight_mode
    ),
    Âge = fmt_int(actuarial_all$initial_age),
    Durée = fmt_int(actuarial_all$duration),
    Taux = fmt_pct(actuarial_all$discount_rate, 0),
    Moyenne = fmt_num(actuarial_all$mean, 5),
    `Q2,5 %` = fmt_num(actuarial_all$q025, 5),
    `Q97,5 %` = fmt_num(actuarial_all$q975, 5),
    check.names = FALSE
  ),
  "actuarial_all_rules_annuity.tex",
  paste0(
    "Rentes par agrégation des forces pour les trois règles ",
    "d'extrapolation et les taux 1\\%, 2\\%, 3\\%."
  ),
  "tab:actuarial-all-rules", align = "lllrrlrrr",
  font_size = "\\tiny"
)

survival_all <- read_csv("46_actuarial_survival_all_rules.csv")
survival_all <- survival_all[
  survival_all$method %in%
    c("stacking_global", "stacking_contextual", "stacking_hierarchical") &
    survival_all$horizon %in% c(1, 5, 10, 15, 20, 25),
]
write_longtable(
  data.frame(
    Analyse = survival_all$analysis,
    Règle = survival_all$rule,
    Méthode = actuarial_method_label(
      survival_all$method, survival_all$weight_mode
    ),
    Âge = fmt_int(survival_all$initial_age),
    Horizon = fmt_int(survival_all$horizon),
    Moyenne = fmt_num(survival_all$mean, 6),
    `Q2,5 %` = fmt_num(survival_all$q025, 6),
    `Q97,5 %` = fmt_num(survival_all$q975, 6),
    check.names = FALSE
  ),
  "actuarial_survival_selected_rules.tex",
  paste0(
    "Probabilités de survie par agrégation des forces, ",
    "sélectionnées selon la règle de poids."
  ),
  "tab:actuarial-survival-rules", align = "lllrrrrr",
  font_size = "\\tiny"
)

hierarchy_uncertainty <- read_csv(
  "51_actuarial_hierarchical_weight_uncertainty.csv"
)
hierarchy_uncertainty <- hierarchy_uncertainty[
  hierarchy_uncertainty$rule == "within_lfo" &
    abs(hierarchy_uncertainty$discount_rate - 0.02) < 1e-12,
]
write_longtable(
  data.frame(
    Âge = fmt_int(hierarchy_uncertainty$initial_age),
    `Delta moy.` =
      fmt_num(hierarchy_uncertainty$mean_difference, 5),
    `SD moy.` =
      fmt_num(hierarchy_uncertainty$sd_mean_weights, 5),
    `SD prop.` =
      fmt_num(hierarchy_uncertainty$sd_propagated_weights, 5),
    `Delta SD` =
      fmt_num(hierarchy_uncertainty$sd_difference, 5),
    `L95 moy.` =
      fmt_num(hierarchy_uncertainty$width95_mean_weights, 5),
    `L95 prop.` =
      fmt_num(hierarchy_uncertainty$width95_propagated_weights, 5),
    `Delta L95` =
      fmt_num(hierarchy_uncertainty$width95_difference, 5),
    check.names = FALSE
  ),
  "actuarial_hierarchical_uncertainty.tex",
  paste0(
    "Contribution de l'incertitude des poids hiérarchiques ",
    "aux rentes principales, taux 2\\%."
  ),
  "tab:actuarial-hierarchical-uncertainty",
  align = "rrrrrrrr",
  font_size = "\\scriptsize"
)

actuarial_checks <- read_csv("52_actuarial_force_aggregation_checks.csv")
write_longtable(
  data.frame(
    Contrôle = gsub("_", " ", actuarial_checks$check),
    Réussi = ifelse(actuarial_checks$passed, "oui", "non"),
    Valeur = fmt_num(actuarial_checks$value, 8),
    Tolérance = fmt_num(actuarial_checks$tolerance, 8),
    Détail = actuarial_checks$detail,
    check.names = FALSE
  ),
  "actuarial_validation_checks.tex",
  "Contrôles automatiques du calcul actuariel par agrégation des forces.",
  "tab:actuarial-validation-checks",
  align = "llrrl",
  font_size = "\\tiny"
)

# Annexe D : 4 100 lignes, une par méthode, âge, horizon et modèle.
nonreg <- read_csv("05_nonregularized_weights.csv")
nonreg <- nonreg[nonreg$horizon <= 10, ]
nonreg$method <- "Contextuel non régularisé"
nonreg$mean <- nonreg$weight
nonreg$q025 <- NA_real_
nonreg$q975 <- NA_real_
hier <- read_csv("06_hierarchical_weights.csv")
hier <- hier[hier$horizon <= 10, ]
hier$method <- "Contextuel hiérarchique"
weights <- rbind(
  nonreg[c("method", "age", "horizon", "model", "mean", "q025", "q975")],
  hier[c("method", "age", "horizon", "model", "mean", "q025", "q975")]
)
weights <- weights[order(
  weights$method, weights$horizon, weights$age, weights$model
), ]
write_longtable(
  data.frame(
    Méthode = weights$method,
    Âge = fmt_int(weights$age),
    Horizon = fmt_int(weights$horizon),
    Modèle = toupper(weights$model),
    Poids = fmt_num(weights$mean, 6),
    `Q2,5 %` = fmt_num(weights$q025, 6),
    `Q97,5 %` = fmt_num(weights$q975, 6),
    check.names = FALSE
  ),
  "weights_age_horizon_model_all.tex",
  "Poids détaillés par âge, horizon et modèle (horizons validés 1 à 10).",
  "tab:weights-all", align = "lrrlrrr", font_size = "\\tiny"
)

# Stabilité Monte Carlo et simulations, produits par les étapes longues.
mc <- read_csv("50_prediction_mc_stability_summary.csv")
write_longtable(
  data.frame(
    Méthode = method_label(mc$method),
    Poids = mc$weight_mode,
    Tirages = fmt_int(mc$forecast_draws),
    `LogS moy.` = fmt_num(mc$logs_mean_over_seeds, 4),
    `Éc.-type LogS` = fmt_num(mc$logs_sd_over_seeds, 6),
    `CRPS moy.` = fmt_num(mc$crps_mean_over_seeds, 7),
    `Éc.-type CRPS` = fmt_num(mc$crps_sd_over_seeds, 8),
    check.names = FALSE
  ),
  "prediction_mc_stability.tex",
  "Stabilité Monte Carlo sur trois graines.",
  "tab:mc-stability", align = "llrrrrr", font_size = "\\scriptsize"
)

simulation <- read_csv("simulation_summary.csv", simulation_dir)
write_longtable(
  data.frame(
    Scénario = simulation$scenario,
    Méthode = method_label(simulation$method),
    LogS = fmt_num(simulation$mean_logs, 4),
    Rang.LogS = fmt_int(simulation$rank_logs),
    CRPS = fmt_num(simulation$mean_crps, 7),
    Rang.CRPS = fmt_int(simulation$rank_crps),
    check.names = FALSE
  ),
  "simulation_summary_full.tex",
  "Performances prédictives de l'étude de simulation.",
  "tab:simulation-summary", align = "llrrrr", font_size = "\\small"
)

recovery <- read_csv("simulation_weight_recovery.csv", simulation_dir)
write_longtable(
  data.frame(
    Scénario = recovery$scenario,
    Méthode = method_label(recovery$method),
    Biais = fmt_num(recovery$weight_bias, 5),
    RMSE = fmt_num(recovery$weight_rmse, 5),
    `Erreur abs.` =
      fmt_num(recovery$mean_absolute_weight_error, 5),
    `Poids extrêmes` =
      fmt_pct(recovery$extreme_weight_proportion, 1),
    `Couv. post. 95 %` =
      fmt_pct(recovery$posterior_95_coverage, 1),
    check.names = FALSE
  ),
  "simulation_weight_recovery.tex",
  "Récupération des poids vrais dans les quatre scénarios.",
  "tab:simulation-recovery", align = "llrrrrr", font_size = "\\scriptsize"
)

message("Tableaux complémentaires écrits dans : ", tables_dir)
