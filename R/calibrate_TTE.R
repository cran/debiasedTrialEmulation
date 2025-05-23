#' @title Negative Control Calibration for Target Trial Emulation
#' @description Performs calibration with negative control outcomes to further reduce confounding bias in
#' Target Trial Emulation results.
#' @import EmpiricalCalibration
#' @import ggplot2
#' @param tte_obj Optional. An object of class "TTE" from the TTE_pipeline function. If provided, the function
#'               will use the negative control results from this object. Either this or both custom_results
#'               and custom_nco_results must be provided.
#' @param custom_results Optional. A data frame containing the primary outcome results if no TTE object is available.
#'                      Must contain columns 'names_outcome', 'logEst', and 'seLogEst'.
#' @param custom_nco_results Optional. A data frame containing the negative control outcome results if no TTE object is available.
#'                          Must contain columns 'names_outcome', 'logEst', and 'seLogEst'.
#' @return An object of class "dTTE" containing both the TTE results and calibration results
#'
#' @examples
#' library("dplyr")
#' data(demo_data)
#' # First run TTE pipeline
#' xvars <- c("eth_cat", "age_cat", "sex", "cohort_entry_month", "obese", "pmca_index", "n_ed", "n_inpatient",
#'            "n_tests", "imm_date_diff_grp", "medical_1", "medical_2", "medical_3", "medical_4", "medical_5")
#' yvars1 <- colnames(demo_data %>% select(starts_with("visits_")))
#' ncovars1 <- colnames(demo_data %>% select(starts_with("nco_visits_")))
#' tte_result <- TTE_pipeline(demo_data, xvars=xvars, yvars=yvars1, ncovars=ncovars1,
#'                          ps_type="Matching", outcome_measure="RR")
#'
#' # Then calibrate results
#' dtte_result <- calibrate_TTE(tte_obj = tte_result)
#'
#' # Alternatively, provide custom results
#' df_results <- data.frame(
#'   names_outcome = c("outcome1", "outcome2"),
#'   logEst = c(0.1, -0.2),
#'   seLogEst = c(0.05, 0.08)
#' )
#' df_nco_results <- data.frame(
#'   names_outcome = c("nco1", "nco2", "nco3"),
#'   logEst = c(0.02, -0.03, 0.01),
#'   seLogEst = c(0.03, 0.04, 0.02)
#' )
#' dtte_result <- calibrate_TTE(custom_results = df_results, custom_nco_results = df_nco_results)
#' @export

calibrate_TTE <- function(tte_obj = NULL, custom_results = NULL, custom_nco_results = NULL) {
  # Validate input
  if (is.null(tte_obj) && (is.null(custom_results) || is.null(custom_nco_results))) {
    stop("Either a TTE object or both custom_results and custom_nco_results must be provided")
  }

  # Get results from TTE object or custom inputs
  if (!is.null(tte_obj)) {
    if (!inherits(tte_obj, "TTE")) {
      stop("tte_obj must be of class 'TTE'")
    }
    df_rslt <- tte_obj$df_rslt
    df_rslt_nco <- tte_obj$df_rslt_nco

    # Extract other components from TTE object
    res <- tte_obj$res
  } else {
    # Use custom results
    required_cols <- c("names_outcome", "logEst", "seLogEst")

    if (!all(required_cols %in% colnames(custom_results)) ||
        !all(required_cols %in% colnames(custom_nco_results))) {
      stop("custom_results and custom_nco_results must contain columns: names_outcome, logEst, seLogEst")
    }

    df_rslt <- custom_results
    df_rslt_nco <- custom_nco_results

    # Create placeholder components if using custom results
    res <- list(df_rslt = df_rslt, df_rslt_nco = df_rslt_nco)
  }

  # Ensure numeric values
  df_rslt$logEst <- suppressWarnings(as.numeric(df_rslt$logEst))
  df_rslt$seLogEst <- suppressWarnings(as.numeric(df_rslt$seLogEst))
  df_rslt_nco$logEst <- suppressWarnings(as.numeric(df_rslt_nco$logEst))
  df_rslt_nco$seLogEst <- suppressWarnings(as.numeric(df_rslt_nco$seLogEst))

  # Remove rows with NA values
  df_rslt <- df_rslt[!is.na(df_rslt$logEst) & !is.na(df_rslt$seLogEst), ]
  df_rslt_nco <- df_rslt_nco[!is.na(df_rslt_nco$logEst) & !is.na(df_rslt_nco$seLogEst), ]

  # Check if we have enough data for calibration
  if (nrow(df_rslt) == 0) {
    stop("No valid primary outcome results available for calibration")
  }
  if (nrow(df_rslt_nco) == 0) {
    stop("No valid negative control results available for calibration")
  }

  # Perform negative control calibration
  tryCatch({
    fitnull <- fitNull(df_rslt_nco$logEst, df_rslt_nco$seLogEst)

    plot.nc <- plotCalibrationEffect(logRrNegatives = df_rslt_nco$logEst,
                                     seLogRrNegatives = df_rslt_nco$seLogEst,
                                     showExpectedAbsoluteSystematicError = TRUE,
                                     null = fitnull)

    model <- fitSystematicErrorModel(df_rslt_nco$logEst, df_rslt_nco$seLogEst, rep(0, length(df_rslt_nco$logEst)))

    rslt.cal <- calibrateConfidenceInterval(logRr = df_rslt$logEst,
                                            seLogRr = df_rslt$seLogEst,
                                            model, ciWidth = 0.95)

    # Extract outcome names from df_rslt
    out_names <- df_rslt$names_outcome

    rslt.cal <- cbind(out_names, rslt.cal)
    colnames(rslt.cal) <- c("out.names", "logEst", "ll", "ul", "selogEst")
    rslt.cal <- data.frame(rslt.cal)

    # Create the final dTTE object
    result <- list(
      res = res,
      df_rslt = df_rslt,
      df_rslt_nco = df_rslt_nco,
      plot.nc = plot.nc,
      rslt.cal = rslt.cal
    )

    class(result) <- "dTTE"

    return(result)
  }, error = function(e) {
    stop(paste("Error during calibration:", e$message))
  })
}
