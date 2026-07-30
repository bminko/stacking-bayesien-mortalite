#!/usr/bin/env Rscript

# Reprise légère de la fin de l'étape 10.
#
# Utile lorsqu'un périphérique graphique s'interrompt après que les tableaux
# 15--45 ont déjà été produits. Aucun calcul statistique n'est recommencé.

source(file.path("R", "utils.R"))
source(file.path("R", "reporting.R"))
cfg <- load_config()
results_dir <- file.path(cfg$paths$results, "chapter3_results")
figures_dir <- file.path(cfg$root, "output", "pdf", "figures")
model_labels <- c(lc = "LC", rh = "RH", apc = "APC", cbd = "CBD", m6 = "M6")
colors <- chapter_model_colors()
method_labels <- chapter_method_labels()

software <- data.frame(
  component = c(
    "R executable used by pipeline", "Stan compiler",
    "RStan build note", "operating system", "chains",
    "parallel cores", "warmup per chain", "kept draws per chain",
    "total kept draws", "adapt_delta", "max_treedepth",
    "metric LC and CBD", "metric RH, APC and M6",
    "metric hierarchical stacking"
  ),
  value = c(
    "R 4.3.1 (command archived in README)",
    "stanc 2.32.2 (embedded in compiled model)",
    "RStan binary compiled under R 4.3.3; exact package patch not archived",
    "Windows x86_64 mingw32",
    cfg$mcmc$chains, cfg$mcmc$parallel_chains,
    cfg$mcmc$iter_warmup, cfg$mcmc$iter_sampling,
    cfg$mcmc$chains * cfg$mcmc$iter_sampling,
    cfg$mcmc$adapt_delta, cfg$mcmc$max_treedepth,
    cfg$mcmc$metric, "dense_e", cfg$hierarchical_mcmc$metric
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(
  software,
  file.path(results_dir, "24_software_and_sampling_settings.csv")
)

weight_difference <- utils::read.csv(
  file.path(results_dir, "45_hierarchical_vs_nonregularized_weights.csv"),
  stringsAsFactors = FALSE
)

render_chapter_figure(
  "hierarchical_nonregularized_weight_difference",
  function() {
    graphics::par(mfrow = c(2, 3), mar = c(3.3, 3.5, 2.1, 0.8))
    for (model in cfg$models) {
      label <- model_labels[[model]]
      part <- weight_difference[weight_difference$model == label, ]
      z <- xtabs(
        difference_hier_minus_nonreg ~ age + horizon,
        data = part
      )
      graphics::image(
        as.integer(rownames(z)), as.integer(colnames(z)), z,
        col = grDevices::hcl.colors(31, "Blue-Red 3", rev = TRUE),
        zlim = c(-1, 1), xlab = "Age", ylab = "Horizon",
        main = label
      )
    }
    graphics::plot.new()
    graphics::text(
      0.5, 0.55,
      "Poids hiérarchique -\npoids non régularisé",
      font = 2
    )
  },
  figures_dir, width = 10, height = 7
)

render_chapter_figure(
  "hierarchical_weight_uncertainty",
  function() {
    graphics::par(mfrow = c(2, 3), mar = c(3.3, 3.5, 2.1, 0.8))
    for (model in cfg$models) {
      label <- model_labels[[model]]
      part <- weight_difference[weight_difference$model == label, ]
      z <- xtabs(
        posterior_interval_width ~ age + horizon,
        data = part
      )
      graphics::image(
        as.integer(rownames(z)), as.integer(colnames(z)), z,
        col = grDevices::hcl.colors(25, "YlOrRd", rev = TRUE),
        zlim = c(0, 1), xlab = "Age", ylab = "Horizon",
        main = label
      )
    }
    graphics::plot.new()
    graphics::text(
      0.5, 0.55,
      "Largeur de l'intervalle\npostérieur 95 % du poids",
      font = 2
    )
  },
  figures_dir, width = 10, height = 7
)

actuarial_survival <- rbind(
  utils::read.csv(
    file.path(
      cfg$paths$metrics,
      "actuarial_survival_summary_force_aggregation_fixed_weights.csv"
    ),
    stringsAsFactors = FALSE
  ),
  utils::read.csv(
    file.path(
      cfg$paths$metrics,
      paste0(
        "actuarial_survival_summary_",
        "force_aggregation_propagated_weights.csv"
      )
    ),
    stringsAsFactors = FALSE
  )
)
actuarial_annuity <- rbind(
  utils::read.csv(
    file.path(
      cfg$paths$metrics,
      "actuarial_annuity_summary_force_aggregation_fixed_weights.csv"
    ),
    stringsAsFactors = FALSE
  ),
  utils::read.csv(
    file.path(
      cfg$paths$metrics,
      paste0(
        "actuarial_annuity_summary_",
        "force_aggregation_propagated_weights.csv"
      )
    ),
    stringsAsFactors = FALSE
  )
)
actuarial_paired <- utils::read.csv(
  file.path(
    cfg$paths$metrics,
    "actuarial_annuity_paired_differences_force_aggregation.csv"
  ),
  stringsAsFactors = FALSE
)
write_csv_atomic(
  actuarial_survival,
  file.path(results_dir, "46_actuarial_survival_all_rules.csv")
)
write_csv_atomic(
  actuarial_annuity,
  file.path(results_dir, "47_actuarial_annuity_all_rules.csv")
)
write_csv_atomic(
  actuarial_paired,
  file.path(results_dir, "48_actuarial_paired_differences.csv")
)

write_csv_atomic(
  utils::read.csv(
    file.path(
      cfg$paths$metrics,
      paste0(
        "actuarial_hierarchical_uncertainty_",
        "force_aggregation_propagated_weights.csv"
      )
    ),
    stringsAsFactors = FALSE
  ),
  file.path(
    results_dir, "51_actuarial_hierarchical_weight_uncertainty.csv"
  )
)
write_csv_atomic(
  utils::read.csv(
    file.path(
      cfg$paths$metrics, "actuarial_force_aggregation_checks.csv"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(results_dir, "52_actuarial_force_aggregation_checks.csv")
)

copy_actuarial_figure <- function(source_name, target_name = source_name) {
  for (extension in c("pdf", "png")) {
    source_path <- file.path(
      cfg$paths$figures, paste0(source_name, ".", extension)
    )
    target_path <- file.path(
      figures_dir, paste0(target_name, ".", extension)
    )
    assert_true(
      file.exists(source_path),
      paste("Figure actuarielle corrigée absente :", source_path)
    )
    assert_true(
      file.copy(source_path, target_path, overwrite = TRUE),
      paste("Impossible de copier la figure actuarielle :", target_path)
    )
  }
}
copy_actuarial_figure(
  "actuarial_survival_force_aggregation_fixed_and_propagated",
  "actuarial_primary_survival"
)
copy_actuarial_figure(
  "actuarial_annuity_distributions_force_aggregation"
)
copy_actuarial_figure(
  "actuarial_annuity_paired_differences_force_aggregation"
)
copy_actuarial_figure(
  "actuarial_hierarchical_weight_uncertainty_force_aggregation"
)

render_chapter_figure(
  "actuarial_annuity_interest_sensitivity",
  function() {
    chapter_plot_theme()
    selected <- c("stacking_global", "stacking_hierarchical")
    part <- actuarial_annuity[
      actuarial_annuity$rule == "within_lfo" &
        actuarial_annuity$weight_mode %in%
          c("fixed_weights", "posterior_mean_weights") &
        actuarial_annuity$method %in% selected,
    ]
    graphics::plot(
      range(part$initial_age), range(part$mean),
      type = "n", xlab = "Âge initial",
      ylab = "Valeur actuarielle",
      main = "Rentes temporaires h<=10 selon le taux d'actualisation"
    )
    line_types <- c("0.01" = 1, "0.02" = 2, "0.03" = 3)
    for (method in selected) {
      for (rate in sort(unique(part$discount_rate))) {
        piece <- part[
          part$method == method &
            abs(part$discount_rate - rate) < 1e-12,
        ]
        graphics::lines(
          piece$initial_age, piece$mean,
          col = colors[[method]],
          lty = line_types[[sprintf("%.2f", rate)]],
          lwd = 2, type = "b", pch = 16
        )
      }
    }
    graphics::legend(
      "bottomleft",
      legend = c("Global 1/2/3 %", "Hiérarchique 1/2/3 %"),
      col = colors[selected], lwd = 2, bty = "n"
    )
  },
  figures_dir, width = 9, height = 5.5
)

message("Finalisation légère de l'étape 10 terminée.")
