#!/usr/bin/env Rscript

# NDHS PAF analysis
# Reproduces the published NDHS anxiety/depression PAF tables and writes
# outputs to outputs/ndhs_paf/.
#
# Before running:
# - Keep this script in the project root.
# - Keep the input file at Data/MergedHHwomen.dta. (fine name can be anything, make sure to replace the file name in data read section)
# - Run the script from the project root so relative paths resolve correctly.
# - Treat Data/ as input-only and outputs/ndhs_paf/ as generated output.
# - Archive old outputs if you want a clean rerun.
#
# Predictor rules:
# - full_risks contains all predictors included in the fully adjusted models.
# - paf_variables contains the candidate modifiable predictors for the PAF step.
# - A variable contributes to PAF only if it is listed in paf_variables and its
#   fully adjusted model term is statistically significant for that outcome.
# - Non-significant predictors remain in the regression model but do not enter
#   the unadjusted, adjusted, or combined PAF calculations.
#
# PAF method:
# - Unadjusted PAF is calculated as Pe * (OR - 1) / OR.
# - Pe is taken from the published Table 2 prevalence definition: weighted
#   exposure distribution within the outcome group from the full survey sample.
# - OR is taken from the fully adjusted multilevel model fitted on the
#   complete-case analytic sample.
# - Tetrachoric correlations and PCA are used to estimate communality among the
#   selected risk factors, and each factor is weighted as 1 - communality.
# - Combined PAF is calculated as 1 - prod(1 - weighted PAF).
# - Adjusted PAF for each factor is its proportional share of the combined PAF.

options(stringsAsFactors = FALSE, scipen = 999)

required_packages <- c(
  "broom.mixed",
  "data.table",
  "datawizard",
  "FactoMineR",
  "forcats",
  "glmmTMB",
  "psych",
  "survey"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Install missing packages before running this script: ",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(normalizePath(getwd()))
  }
  normalizePath(dirname(sub("^--file=", "", file_arg[1])))
}

project_dir <- get_script_dir()
data_path <- file.path(project_dir, "Data", "MergedHHwomen.dta")
output_dir <- file.path(project_dir, "outputs", "ndhs_paf")

if (!file.exists(data_path)) {
  stop(
    paste0(
      "Data file not found at ", data_path,
      ". Run the script from the project folder or update the path."
    ),
    call. = FALSE
  )
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

full_risks <- c(
  "age",
  "province",
  "residence",
  "wealth.quintile",
  "education.lvl",
  "employment",
  "smokesYN",
  "alcohol",
  "menstrual.sep",
  "wife.beating",
  "sexual.nego",
  "decisions",
  "func.diff.pres",
  "hhfs",
  "dv_phy",
  "dv_sex",
  "dv_prtnr_emot"
)

paf_variables <- c("func.diff.pres", "hhfs", "dv_phy", "dv_sex", "dv_prtnr_emot")

# Models use the complete-case analytic sample. Table 2 prevalence inputs use
# the full weighted survey sample within each outcome group.
use_complete_case_sample <- TRUE

prepare_analysis_data <- function(path) {
  raw_data <- data.table(datawizard::data_read(path, convert_factors = TRUE))
  raw_data <- forcats::as_factor(raw_data, only_labelled = TRUE)

  analysis_data <- raw_data[, .(
    cluster = hv001,
    hh.num = hv002,
    resp.line.num = hv003,
    weight = v005 / 1000000,
    age = v012,
    psu = v021,
    strata = v022,
    domain = v023,
    province = v024,
    residence = v025,
    wealth.quintile = v190,
    education.lvl = v106,
    currently.married = currently_married,
    employment = employed,
    smokesYN = smokesYN,
    alcohol = alcohol,
    menstrual.sep = menstrual_sep,
    wife.beating = wife_beating,
    sexual.nego = sexual_nego,
    decisions = decisions,
    func.diff.pres = func_diff_pres,
    dv_phy = dv_phy,
    dv_sex = dv_sex,
    dv_prtnr_emot = dv_prtnr_emot,
    anxiety.symptoms = anxiety,
    depression.symptoms = depression,
    hhfs = hfs_mod
  )]

  analysis_data[, province := factor(
    province,
    levels = c(
      "koshi",
      "madhesh province",
      "bagmati province",
      "gandaki province",
      "lumbini province",
      "karnali province",
      "sudurpashchim province"
    ),
    labels = c(
      "Koshi",
      "Madhesh",
      "Bagmati",
      "Gandaki",
      "Lumbini",
      "Karnali",
      "Sudurpaschim"
    )
  )]

  analysis_data[, province := stats::relevel(province, ref = "Bagmati")]
  analysis_data[, sexual.nego := stats::relevel(sexual.nego, ref = "yes")]

  factor_vars <- c(
    "province",
    "residence",
    "wealth.quintile",
    "education.lvl",
    "currently.married",
    "employment",
    "smokesYN",
    "alcohol",
    "menstrual.sep",
    "wife.beating",
    "sexual.nego",
    "decisions",
    "func.diff.pres",
    "dv_phy",
    "dv_sex",
    "dv_prtnr_emot",
    "anxiety.symptoms",
    "depression.symptoms",
    "hhfs"
  )

  for (var_name in factor_vars) {
    analysis_data[, (var_name) := droplevels(get(var_name))]
  }

  analysis_data[, anxiety_bin := as.integer(anxiety.symptoms == "anxiety")]
  analysis_data[, depression_bin := as.integer(depression.symptoms == "depressed")]

  analysis_data[]
}

build_term_lookup <- function(data, predictors) {
  model_frame <- stats::model.frame(
    stats::reformulate(predictors),
    data = data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )
  model_matrix <- model.matrix(
    stats::reformulate(predictors),
    data = model_frame
  )
  assign_index <- attr(model_matrix, "assign")
  term_labels <- attr(stats::terms(stats::reformulate(predictors)), "term.labels")
  keep_cols <- colnames(model_matrix) != "(Intercept)"

  lookup <- data.table(
    term = colnames(model_matrix)[keep_cols],
    variable = term_labels[assign_index[keep_cols]]
  )
  lookup[, level := substring(term, nchar(variable) + 1L)]
  list(
    lookup = lookup,
    matrix = model_matrix[, keep_cols, drop = FALSE]
  )
}

fit_full_model <- function(outcome, predictors, data) {
  formula_text <- paste(
    outcome,
    "~",
    paste(c(predictors, "(1 | cluster)"), collapse = " + ")
  )

  suppressWarnings(
    glmmTMB::glmmTMB(
      formula = as.formula(formula_text),
      data = data,
      weights = data[["weight"]],
      family = binomial(link = "logit")
    )
  )
}

tidy_or_results <- function(model, outcome_name, lookup) {
  result <- as.data.table(
    broom.mixed::tidy(
      model,
      effects = "fixed",
      conf.int = TRUE,
      exponentiate = TRUE
    )
  )

  result <- result[term != "(Intercept)"]
  result <- merge(result, lookup, by = "term", all.x = TRUE)
  result[, outcome := outcome_name]
  result[, .(
    outcome,
    term,
    variable,
    level,
    estimate,
    conf.low,
    conf.high,
    p.value
  )]
}

# Published Table 2 prevalence is the weighted exposure distribution within the
# outcome group, not the population exposure prevalence.
estimate_exposure_prevalence_in_cases <- function(design, terms, outcome_var, outcome_name) {
  prevalence_rows <- lapply(terms, function(term_name) {
    case_index <- design$variables[[outcome_var]] == 1
    case_index[is.na(case_index)] <- FALSE

    if (!any(case_index)) {
      return(data.table(
        outcome = outcome_name,
        term = term_name,
        prevalence_n = 0L,
        prevalence_pct = NA_real_,
        prevalence_lci_pct = NA_real_,
        prevalence_uci_pct = NA_real_
      ))
    }

    estimate <- survey::svymean(
      as.formula(paste0("~", term_name)),
      design[case_index, ],
      na.rm = TRUE
    )
    ci <- suppressMessages(stats::confint(estimate))

    data.table(
      outcome = outcome_name,
      term = term_name,
      prevalence_n = sum(case_index),
      prevalence_pct = as.numeric(stats::coef(estimate)[1] * 100),
      prevalence_lci_pct = as.numeric(ci[1, 1] * 100),
      prevalence_uci_pct = as.numeric(ci[1, 2] * 100)
    )
  })

  rbindlist(prevalence_rows, use.names = TRUE, fill = TRUE)
}

calculate_pca_weights <- function(data, terms) {
  dummy_data <- as.data.frame(data[, ..terms])
  tetrachoric_corr <- psych::tetrachoric(dummy_data)
  pca <- FactoMineR::PCA(tetrachoric_corr$rho, graph = FALSE)

  eigenvalues <- pca$eig[, "eigenvalue"]
  n_components <- sum(eigenvalues > 1)
  if (n_components < 1) {
    n_components <- 1
  }

  loadings <- pca$svd$V[, seq_len(n_components), drop = FALSE]
  weights <- data.table(term = rownames(tetrachoric_corr$rho))
  weights[, communality := rowSums(loadings^2)]
  weights[, pca_weight := 1 - communality]
  weights[]
}

# Only significant modifiable factors from the full model enter the PAF step.
calculate_paf <- function(result_data, pca_weights, outcome_name, sample_n) {
  paf_data <- merge(
    result_data[outcome == outcome_name & variable %in% paf_variables & p.value < 0.05],
    pca_weights,
    by = "term",
    all.x = TRUE
  )

  if (nrow(paf_data) == 0) {
    empty_summary <- data.table(
      outcome = outcome_name,
      included_terms = 0L,
      combined_paf_pct = NA_real_,
      combined_paf_lci_pct = NA_real_,
      combined_paf_uci_pct = NA_real_
    )
    return(list(detail = paf_data, summary = empty_summary))
  }

  paf_data[, paf_unadj_pct := prevalence_pct * (estimate - 1) / estimate]
  paf_data[, paf_unadj_lci_pct := prevalence_lci_pct * (conf.low - 1) / conf.low]
  paf_data[, paf_unadj_uci_pct := prevalence_uci_pct * (conf.high - 1) / conf.high]

  paf_data[, paf_unadj_prop := paf_unadj_pct / 100]
  paf_data[, paf_weighted_prop := paf_unadj_prop * pca_weight]

  combined_paf_prop <- 1 - prod(1 - paf_data$paf_weighted_prop)
  combined_paf_se <- sqrt(abs(combined_paf_prop * (1 - combined_paf_prop)) / sample_n)

  paf_data[, paf_adjusted_prop := (paf_unadj_prop / sum(paf_unadj_prop)) * combined_paf_prop]
  paf_data[, paf_adjusted_pct := paf_adjusted_prop * 100]

  paf_data[, paf_adjusted_lci_pct := (
    paf_adjusted_prop - 1.96 * sqrt(abs(paf_adjusted_prop * (1 - paf_adjusted_prop)) / sample_n)
  ) * 100]

  paf_data[, paf_adjusted_uci_pct := (
    paf_adjusted_prop + 1.96 * sqrt(abs(paf_adjusted_prop * (1 - paf_adjusted_prop)) / sample_n)
  ) * 100]

  summary <- data.table(
    outcome = outcome_name,
    included_terms = nrow(paf_data),
    combined_paf_pct = combined_paf_prop * 100,
    combined_paf_lci_pct = (combined_paf_prop - 1.96 * combined_paf_se) * 100,
    combined_paf_uci_pct = (combined_paf_prop + 1.96 * combined_paf_se) * 100
  )

  list(detail = paf_data, summary = summary)
}

full_analysis_data <- prepare_analysis_data(data_path)

term_info <- build_term_lookup(full_analysis_data, full_risks)
term_lookup <- term_info$lookup
model_terms <- term_lookup$term
term_columns_to_add <- setdiff(model_terms, names(full_analysis_data))
for (term_name in term_columns_to_add) {
  full_analysis_data[, (term_name) := as.numeric(term_info$matrix[, term_name])]
}

full_survey_design <- survey::svydesign(
  data = full_analysis_data,
  ids = ~1,
  weights = ~weight
)

analysis_data <- copy(full_analysis_data)

if (use_complete_case_sample) {
  complete_case_vars <- c(full_risks, "anxiety.symptoms", "depression.symptoms")
  analysis_data <- analysis_data[complete.cases(analysis_data[, ..complete_case_vars])]
}

analytic_sample_n <- nrow(analysis_data)
paf_ci_n <- nrow(full_analysis_data)

anxiety_model <- fit_full_model("anxiety.symptoms", full_risks, analysis_data)
depression_model <- fit_full_model("depression.symptoms", full_risks, analysis_data)

anxiety_or <- tidy_or_results(anxiety_model, "Anxiety", term_lookup)
depression_or <- tidy_or_results(depression_model, "Depression", term_lookup)
or_results <- rbindlist(list(anxiety_or, depression_or), use.names = TRUE)

# Age has no binary exposure term, so it is excluded from the prevalence table.
prevalence_lookup <- term_lookup[variable != "age"]

anxiety_prev <- estimate_exposure_prevalence_in_cases(
  design = full_survey_design,
  terms = prevalence_lookup$term,
  outcome_var = "anxiety_bin",
  outcome_name = "Anxiety"
)

depression_prev <- estimate_exposure_prevalence_in_cases(
  design = full_survey_design,
  terms = prevalence_lookup$term,
  outcome_var = "depression_bin",
  outcome_name = "Depression"
)

prevalence_results <- rbindlist(list(anxiety_prev, depression_prev), use.names = TRUE)
prevalence_results <- merge(prevalence_results, prevalence_lookup, by = "term", all.x = TRUE)
setcolorder(
  prevalence_results,
  c("outcome", "term", "variable", "level", "prevalence_n", "prevalence_pct", "prevalence_lci_pct", "prevalence_uci_pct")
)

combined_results <- merge(
  or_results,
  prevalence_results,
  by = c("outcome", "term", "variable", "level"),
  all.x = TRUE
)

paf_terms <- term_lookup[variable %in% paf_variables, term]
pca_weights <- calculate_pca_weights(analysis_data, paf_terms)

anxiety_paf <- calculate_paf(combined_results, pca_weights, "Anxiety", paf_ci_n)
depression_paf <- calculate_paf(combined_results, pca_weights, "Depression", paf_ci_n)

paf_detail <- rbindlist(
  list(anxiety_paf$detail, depression_paf$detail),
  use.names = TRUE,
  fill = TRUE
)

paf_summary <- rbindlist(
  list(anxiety_paf$summary, depression_paf$summary),
  use.names = TRUE,
  fill = TRUE
)

risk_factor_paf <- paf_detail[, .(
  outcome,
  risk_factor = variable,
  level,
  prevalence_n,
  risk_factor_prevalence_pct = prevalence_pct,
  risk_factor_prevalence_lci_pct = prevalence_lci_pct,
  risk_factor_prevalence_uci_pct = prevalence_uci_pct,
  odds_ratio = estimate,
  or_lci = conf.low,
  or_uci = conf.high,
  p_value = p.value,
  pca_weight,
  paf_unadj_pct,
  paf_unadj_lci_pct,
  paf_unadj_uci_pct,
  paf_adjusted_pct,
  paf_adjusted_lci_pct,
  paf_adjusted_uci_pct
)]

risk_factor_paf_formatted <- risk_factor_paf[, .(
  outcome,
  risk_factor,
  level,
  prevalence_n,
  prevalence = sprintf("%.1f (%.1f, %.1f)", risk_factor_prevalence_pct, risk_factor_prevalence_lci_pct, risk_factor_prevalence_uci_pct),
  odds_ratio = sprintf("%.2f (%.2f, %.2f)", odds_ratio, or_lci, or_uci),
  p_value = signif(p_value, 3),
  pca_weight = round(pca_weight, 3),
  paf_unadjusted = sprintf("%.1f (%.1f, %.1f)", paf_unadj_pct, paf_unadj_lci_pct, paf_unadj_uci_pct),
  paf_adjusted = sprintf("%.1f (%.1f, %.1f)", paf_adjusted_pct, paf_adjusted_lci_pct, paf_adjusted_uci_pct)
)]

data.table::fwrite(
  or_results,
  file.path(output_dir, "01_model_odds_ratios.csv")
)

data.table::fwrite(
  prevalence_results,
  file.path(output_dir, "02_outcome_prevalence_by_term.csv")
)

data.table::fwrite(
  combined_results,
  file.path(output_dir, "03_model_and_prevalence_combined.csv")
)

data.table::fwrite(
  pca_weights,
  file.path(output_dir, "04_pca_weights.csv")
)

data.table::fwrite(
  paf_detail[, .(
    outcome,
    term,
    variable,
    level,
    estimate,
    conf.low,
    conf.high,
    p.value,
    prevalence_pct,
    prevalence_lci_pct,
    prevalence_uci_pct,
    pca_weight,
    paf_unadj_pct,
    paf_unadj_lci_pct,
    paf_unadj_uci_pct,
    paf_adjusted_pct,
    paf_adjusted_lci_pct,
    paf_adjusted_uci_pct
  )],
  file.path(output_dir, "05_paf_detail.csv")
)

data.table::fwrite(
  paf_summary,
  file.path(output_dir, "06_paf_summary.csv")
)

data.table::fwrite(
  risk_factor_paf,
  file.path(output_dir, "07_risk_factor_paf.csv")
)

data.table::fwrite(
  risk_factor_paf_formatted,
  file.path(output_dir, "08_risk_factor_paf_formatted.csv")
)

saveRDS(
  list(
    full_analysis_data = full_analysis_data,
    analysis_data = analysis_data,
    term_lookup = term_lookup,
    odds_ratios = or_results,
    prevalence = prevalence_results,
    combined_results = combined_results,
    pca_weights = pca_weights,
    paf_detail = paf_detail,
    paf_summary = paf_summary,
    risk_factor_paf = risk_factor_paf,
    risk_factor_paf_formatted = risk_factor_paf_formatted
  ),
  file.path(output_dir, "ndhs_paf_results.rds")
)

capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))

writeLines(
  c(
    "NDHS PAF analysis outputs",
    paste("Run time:", format(Sys.time(), tz = "", usetz = TRUE)),
    paste("Project directory:", project_dir),
    paste("Data file:", data_path),
    paste("Full survey sample size:", nrow(full_analysis_data)),
    paste("Complete-case analytic sample size:", analytic_sample_n),
    "",
    "Files written:",
    "01_model_odds_ratios.csv",
    "02_outcome_prevalence_by_term.csv",
    "03_model_and_prevalence_combined.csv",
    "04_pca_weights.csv",
    "05_paf_detail.csv",
    "06_paf_summary.csv",
    "07_risk_factor_paf.csv",
    "08_risk_factor_paf_formatted.csv",
    "ndhs_paf_results.rds",
    "session_info.txt"
  ),
  file.path(output_dir, "README.txt")
)

message("Analysis completed successfully.")
message("Outputs written to: ", output_dir)
print(paf_summary)
