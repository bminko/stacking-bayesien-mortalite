# Lecture, controle et mise en forme des tables HMD.

read_hmd_table <- function(path) {
  assert_true(file.exists(path), paste("Fichier absent :", path))
  data <- utils::read.table(
    path,
    skip = 2L,
    header = TRUE,
    stringsAsFactors = FALSE,
    na.strings = c(".", "NA", "")
  )
  expected <- c("Year", "Age", "Female", "Male", "Total")
  assert_true(
    identical(names(data), expected),
    paste("Colonnes HMD inattendues dans", basename(path))
  )
  data$Age_num <- suppressWarnings(as.integer(data$Age))
  for (column in c("Female", "Male", "Total")) {
    data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
  }
  data
}

prepare_hmd_data <- function(cfg) {
  deaths <- read_hmd_table(cfg$paths$deaths)
  exposures <- read_hmd_table(cfg$paths$exposure)

  keep_deaths <- !is.na(deaths$Age_num) &
    deaths$Age_num >= cfg$age_min & deaths$Age_num <= cfg$age_max &
    deaths$Year >= cfg$year_start & deaths$Year <= cfg$data_end
  keep_exposure <- !is.na(exposures$Age_num) &
    exposures$Age_num >= cfg$age_min & exposures$Age_num <= cfg$age_max &
    exposures$Year >= cfg$year_start & exposures$Year <= cfg$data_end

  deaths <- deaths[keep_deaths, c("Year", "Age_num", cfg$sex)]
  exposures <- exposures[keep_exposure, c("Year", "Age_num", cfg$sex)]
  names(deaths) <- c("year", "age", "deaths")
  names(exposures) <- c("year", "age", "exposure")

  duplicated_deaths <- sum(duplicated(deaths[c("year", "age")]))
  duplicated_exposure <- sum(duplicated(exposures[c("year", "age")]))
  if (duplicated_deaths || duplicated_exposure) {
    stop(
      "PROBLEME DE DONNEES : cellules age-annee dupliquees. ",
      "Deces=", duplicated_deaths, ", expositions=", duplicated_exposure,
      ". Aucun traitement automatique n'a ete applique.",
      call. = FALSE
    )
  }

  data <- merge(
    deaths, exposures,
    by = c("year", "age"),
    all = TRUE,
    sort = TRUE
  )
  expected_rows <- (cfg$age_max - cfg$age_min + 1L) *
    (cfg$data_end - cfg$year_start + 1L)

  issues <- character()
  if (nrow(data) != expected_rows) {
    issues <- c(issues, sprintf(
      "grille incomplete : %d cellules au lieu de %d",
      nrow(data), expected_rows
    ))
  }
  if (anyNA(data$deaths)) issues <- c(issues, "deces manquants")
  if (anyNA(data$exposure)) issues <- c(issues, "expositions manquantes")
  if (any(data$exposure <= 0, na.rm = TRUE)) {
    issues <- c(issues, "expositions nulles ou negatives")
  }
  if (any(data$deaths < 0, na.rm = TRUE)) {
    issues <- c(issues, "deces negatifs")
  }

  fractional <- abs(data$deaths - round(data$deaths)) > 1e-9
  if (any(fractional, na.rm = TRUE)) {
    examples <- head(data[fractional, c("year", "age", "deaths")], 8L)
    example_text <- paste(
      apply(examples, 1L, paste, collapse = "/"),
      collapse = ", "
    )
    issues <- c(
      issues,
      paste0(
        sum(fractional, na.rm = TRUE),
        " deces fractionnaires (annee/age/deces : ", example_text, ")"
      )
    )
  }

  if (length(issues)) {
    stop(
      "PROBLEME DE DONNEES : ", paste(issues, collapse = " ; "),
      ". Le pipeline est arrete sans arrondi ni imputation. ",
      "Merci d'indiquer le traitement souhaite.",
      call. = FALSE
    )
  }

  data$deaths <- as.integer(round(data$deaths))
  data <- data[order(data$year, data$age), ]
  rownames(data) <- NULL

  ages <- seq.int(cfg$age_min, cfg$age_max)
  years <- seq.int(cfg$year_start, cfg$data_end)
  death_matrix <- t(vapply(
    ages,
    function(age) data$deaths[data$age == age],
    integer(length(years))
  ))
  exposure_matrix <- t(vapply(
    ages,
    function(age) data$exposure[data$age == age],
    numeric(length(years))
  ))
  dimnames(death_matrix) <- list(age = ages, year = years)
  dimnames(exposure_matrix) <- list(age = ages, year = years)

  list(
    long = data,
    ages = ages,
    years = years,
    D = death_matrix,
    E = exposure_matrix,
    metadata = list(
      country = "Belgium",
      sex = cfg$sex,
      age_min = cfg$age_min,
      age_max = cfg$age_max,
      year_start = cfg$year_start,
      year_end = cfg$data_end,
      source_deaths = basename(cfg$paths$deaths),
      source_exposure = basename(cfg$paths$exposure)
    )
  )
}

data_quality_report <- function(processed, cfg) {
  data <- processed$long
  data.frame(
    profile = cfg$profile,
    sex = cfg$sex,
    first_year = min(data$year),
    last_year = max(data$year),
    min_age = min(data$age),
    max_age = max(data$age),
    cells = nrow(data),
    missing_deaths = sum(is.na(data$deaths)),
    missing_exposure = sum(is.na(data$exposure)),
    duplicate_cells = sum(duplicated(data[c("year", "age")])),
    nonpositive_exposure = sum(data$exposure <= 0),
    fractional_deaths = sum(abs(data$deaths - round(data$deaths)) > 1e-9),
    min_deaths = min(data$deaths),
    max_deaths = max(data$deaths),
    min_exposure = min(data$exposure),
    max_exposure = max(data$exposure),
    stringsAsFactors = FALSE
  )
}

subset_training_data <- function(processed, end_year) {
  keep <- processed$years <= end_year
  assert_true(any(keep), paste("Aucune annee avant", end_year))
  list(
    ages = processed$ages,
    years = processed$years[keep],
    D = processed$D[, keep, drop = FALSE],
    E = processed$E[, keep, drop = FALSE],
    end_year = as.integer(end_year)
  )
}

stan_data_for_model <- function(training, model) {
  ages <- training$ages
  years <- training$years
  A <- length(ages)
  T <- length(years)
  common <- list(
    A = A,
    T = T,
    D = unname(training$D),
    log_E = log(unname(training$E))
  )

  if (model %in% c("cbd", "m6")) {
    common$age_centered <- as.numeric(ages - mean(ages))
    common$age_scale <- as.numeric(stats::sd(ages))
    assert_true(
      is.finite(common$age_scale) && common$age_scale > 0,
      "L'echelle des ages doit etre strictement positive."
    )
  }
  if (model %in% c("rh", "apc", "m6")) {
    cohort_values <- outer(ages, years, function(age, year) year - age)
    cohort_min <- min(cohort_values)
    cohort_max <- max(cohort_values)
    common$C <- as.integer(cohort_max - cohort_min + 1L)
    common$cohort_idx <- unname(cohort_values - cohort_min + 1L)
  }
  common
}

target_grid_from_observed <- function(processed, origin, horizons,
                                      maximum_year) {
  horizons <- horizons[(origin + horizons) <= maximum_year]
  if (!length(horizons)) return(NULL)
  target <- processed$long[
    processed$long$year %in% (origin + horizons),
    c("year", "age", "deaths", "exposure")
  ]
  target$origin <- as.integer(origin)
  target$horizon <- as.integer(target$year - origin)
  target <- target[order(target$horizon, target$age), ]
  rownames(target) <- NULL
  target
}
