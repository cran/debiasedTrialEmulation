#' @title Print Method for TTE Objects
#' @description Prints a concise summary of the TTE pipeline output.
#' @param x An object of class "TTE".
#' @param ... Additional arguments (currently ignored).
#' @export
print.TTE <- function(x, ...) {
  cat("Target Trial Emulation (TTE) Object\n")
  cat("-------------------------------------------------\n")
  cat("Primary results (first few rows):\n")
  print(utils::head(x$df_rslt))

  if (!is.null(x$df_rslt_nco) && nrow(x$df_rslt_nco) > 0) {
    cat("\nNegative control results (first few rows):\n")
    print(utils::head(x$df_rslt_nco))
  } else {
    cat("\nNo negative control results available.\n")
  }

  cat("\nFor more details, try summary() or plot().\n")
  invisible(x)
}

#' @title Summary Method for TTE Objects
#' @description Provides a detailed summary of the TTE pipeline output.
#' @param x An object of class "TTE".
#' @param ... Additional arguments (currently ignored).
#' @export
summary.TTE <- function(x, ...) {
  cat("Summary of Target Trial Emulation (TTE) Results\n")
  cat("-------------------------------------------------\n")

  cat("\nPrimary Result Data Frame:\n")
  print(summary(x$df_rslt))

  if (!is.null(x$df_rslt_nco) && nrow(x$df_rslt_nco) > 0) {
    cat("\nNegative Control Result Data Frame:\n")
    print(summary(x$df_rslt_nco))
  } else {
    cat("\nNo negative control results available.\n")
  }

  invisible(x)
}

#' @title Plot Method for TTE Objects
#' @description Plots diagnostic graphs for the TTE pipeline output.
#' @param x An object of class "TTE".
#' @param which A character vector specifying which plot(s) to display.
#'              Options are "SMD" (Standardized Mean Differences) and
#'              "Equipoise" (Equipoise plot). Default shows all plots.
#' @param ... Additional arguments passed to plotting functions.
#' @export
plot.TTE <- function(x, which = c("SMD", "Equipoise"), ...) {
  # Ensure valid plot choices:
  which <- match.arg(which, choices = c("SMD", "Equipoise"), several.ok = TRUE)

  if ("SMD" %in% which) {
    message("Displaying Standardized Mean Difference plot...")
    print(x$p.SMD)
  }
  if ("Equipoise" %in% which) {
    message("Displaying Equipoise plot...")
    print(x$p.equipoise)
  }

  invisible(x)
}

#' @title Print Method for dTTE Objects
#' @description Prints a concise summary of the dTTE pipeline output.
#' @param x An object of class "dTTE".
#' @param ... Additional arguments (currently ignored).
#' @export
print.dTTE <- function(x, ...) {
  cat("Debiased Target Trial Emulation (dTTE) Object\n")
  cat("-------------------------------------------------\n")
  cat("Calibrated results (first few rows):\n")
  print(utils::head(x$rslt.cal))
  cat("\nFor more details, try summary() or plot().\n")
  invisible(x)
}

#' @title Summary Method for dTTE Objects
#' @description Provides a detailed summary of the dTTE pipeline output.
#' @param x An object of class "dTTE".
#' @param ... Additional arguments (currently ignored).
#' @export
summary.dTTE <- function(x, ...) {
  cat("Summary of dTTE Pipeline Results\n")
  cat("-------------------------------------------------\n")

  cat("\nPrimary Result Data Frame:\n")
  print(summary(x$df_rslt))

  cat("\nNegative Control Result Data Frame:\n")
  print(summary(x$df_rslt_nco))

  cat("\nCalibrated Results (rslt.cal) - first few rows:\n")
  print(utils::head(x$rslt.cal))

  invisible(x)
}


#' @title Plot Method for dTTE Objects (Calibration Only)
#' @description Plots only the calibration graph for the dTTE pipeline output.
#' @param x An object of class "dTTE".
#' @param ... Additional arguments passed to plotting functions.
#' @export
plot.dTTE <- function(x, ...) {
  message("Displaying Calibration plot...")
  print(x$plot.nc)

  invisible(x)
}
