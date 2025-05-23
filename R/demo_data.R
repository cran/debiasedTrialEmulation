#' @title Example Dataset for Debiased Trial Emulation
#'
#' @description
#' `demo_data` is a simulated dataset used to demonstrate the functionality of the `debiasedTrialEmulation` package.
#' It includes patient demographic information, treatment assignment, covariates, clinical outcomes,
#' and negative control outcomes for evaluating treatment effects using propensity score methods.
#'
#' The dataset contains 50,000 observations and 93 variables, including:
#'
#' - **Demographic variables**: Ethnicity, age, sex, and cohort entry month.
#' - **Treatment assignment**: Binary treatment indicator.
#' - **Covariates**: Baseline health conditions and healthcare utilization variables.
#' - **Primary outcomes**: Binary and time-to-event outcomes related to cardiovascular health.
#' - **Negative control outcomes (NCOs)**: Outcomes used for bias calibration.
#'
#' @format A data frame with 50,000 rows and 93 variables
#' @name demo_data
#' @docType data
#' @keywords datasets
#' @usage data(demo_data)
NULL
