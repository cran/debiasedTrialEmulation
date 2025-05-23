#' @title Target Trial Emulation (TTE) Pipeline
#' @description Implements a Target Trial Emulation pipeline using propensity score methods, including matching, weighting, and stratification.
#' @import dplyr
#' @import janitor
#' @import cobalt
#' @import MatchIt
#' @import geex
#' @import glmnet
#' @import survival
#' @import ggplot2
#' @import ParallelLogger
#' @param data A dataset containing treatment assignment, covariates, and outcomes.
#' @param xvars A character vector of covariate names used for propensity score estimation.
#' @param yvars A character vector of primary outcome variable names.
#' @param ncovars Optional. A character vector of negative control outcome variable names.
#' @param ps_type The propensity score method: "Matching", "Stratification", or "Weighting".
#' @param outcome_measure The outcome measure to estimate: "RR" (Risk Ratio), "OR" (Odds Ratio), or "HR" (Hazard Ratio).
#' @return An object of class "TTE" containing the propensity score analysis results
#'
#' @examples
#' library("dplyr")
#' data(demo_data)
#' xvars <- c("eth_cat", "age_cat", "sex", "cohort_entry_month", "obese", "pmca_index", "n_ed", "n_inpatient",
#'            "n_tests", "imm_date_diff_grp", "medical_1", "medical_2", "medical_3", "medical_4", "medical_5")
#' yvars1 <- colnames(demo_data %>% select(starts_with("visits_")))
#' yvars2 <- colnames(demo_data %>% select(starts_with("event_")))
#' # without negative controls
#' TTE_pipeline(demo_data, xvars=xvars, yvars=yvars1, ps_type="Matching", outcome_measure="RR")
#' # with negative controls
#' ncovars1 <- colnames(demo_data %>% select(starts_with("nco_visits_")))
#' TTE_pipeline(demo_data, xvars=xvars, yvars=yvars1, ncovars=ncovars1, ps_type="Matching", outcome_measure="RR")
#' @export

TTE_pipeline <- function(data,
                        xvars,
                        yvars,
                        ncovars = NULL,
                        ps_type,
                        outcome_measure) {
  
  # Input validation
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  if (!is.character(xvars)) {
    stop("xvars must be a character vector")
  }
  if (!is.character(yvars)) {
    stop("yvars must be a character vector")
  }
  if (!is.null(ncovars) && !is.character(ncovars)) {
    stop("ncovars must be a character vector or NULL")
  }
  if (!ps_type %in% c("Matching", "Stratification", "Weighting")) {
    stop("ps_type must be one of: 'Matching', 'Stratification', 'Weighting'")
  }
  if (!outcome_measure %in% c("RR", "OR", "HR")) {
    stop("outcome_measure must be one of: 'RR', 'OR', 'HR'")
  }
  
  # Clean column names
  colnames(data) <- janitor::make_clean_names(colnames(data))
  
  # Create formula
  form <- as.formula(paste("treatment ~ ", paste(xvars, collapse = "+")))
  
  # Run appropriate analysis based on ps_type
  if (ps_type == "Matching") {
    res <- estEffect_matching(
      form = form,
      data = data,
      yvars = yvars,
      ncovars = ncovars,
      distance = "glm",
      outcome_measure = outcome_measure
    )
    p.SMD <- plot_SMD_matching(res$m.out)
    p.equipoise <- plot_Equipoise_matching(data, res$m.out)
  } else if (ps_type == "Stratification") {
    res <- estEffect_stratification(
      form = form,
      data = data,
      yvars = yvars,
      ncovars = ncovars,
      distance = "glm",
      outcome_measure = outcome_measure
    )
    p.SMD <- plot_SMD_stratification(res$stratifiedPop, xvars)
    p.equipoise <- plot_Equipoise_stratification(res$stratifiedPop)
  } else if (ps_type == "Weighting") {
    res <- estEffect_weighting(
      form = form,
      data = data,
      yvars = yvars,
      ncovars = ncovars,
      distance = "glm",
      outcome_measure = outcome_measure
    )
    p.SMD <- plot_SMD_weighting(res$data.trimmed)
    p.equipoise <- plot_Equipoise_weighting(res$data.trimmed)
  }
  
  # Create the result object
  result <- list(
    res = res,
    p.SMD = p.SMD,
    p.equipoise = p.equipoise,
    df_rslt = res$df_rslt,
    df_rslt_nco = res$df_rslt_nco
  )
  
  class(result) <- "TTE"
  
  return(result)
}
