#' @title Compute Standardized Mean Differences (SMD)
#' @description Computes the standardized mean differences for covariates before and after adjustment.
#' @param data A data frame containing covariates.
#' @param treat A binary variable indicating treatment assignment.
#' @param weights Optional weight vector.
#' @param std Logical; whether to standardize.
#' @return A data frame with covariate names and their standardized mean differences.
#' @export
GetSMD <- function(data, treat, weights = NULL, std = TRUE){
  table <- col_w_smd(data, treat, weights, std)
  smd = t(data.frame(as.list(table)))
  smd = data.frame(row.names(smd), smd)
  rownames(smd) <- NULL
  colnames(smd) = c("Xvars", "SMD")
  #smd$SMD = abs(smd$SMD)

  return(smd)
}

#' @title Compute Weights for Stratification
#' @description Computes inverse probability weights for stratification-based adjustment.
#' @param data A data frame containing strata IDs and treatment assignments.
#' @return A data frame with row IDs and computed weights.
#' @export
Compute_weight <- function(data){
  stratumSize <- data %>% group_by(stratumId, treatment) %>% count() %>% ungroup()

  w <- stratumSize %>% mutate(weight = 1/n) %>%
    inner_join(data, by = c("stratumId", "treatment"), multiple = "all") %>% select(rowId, treatment, weight)

  wSum <- w %>% group_by(treatment) %>% summarize(wSum = sum(weight, na.rm = TRUE)) %>% ungroup()

  w_final <- w %>% inner_join(wSum, by = "treatment") %>%
    mutate(weight = weight / wSum) %>% select(rowId, treatment, weight)

  return(w_final[,c(1,3)])
}

#' @title Compute Propensity Score Weights
#' @description Computes inverse probability treatment weights (IPTW) for ATE or ATT estimation.
#' @param population A data frame containing treatment assignments and propensity scores.
#' @param estimator Type of estimator, either "ate" (average treatment effect) or "att" (average treatment effect on the treated).
#' @return A vector of computed weights.
#' @export
computeWeights <- function(population, estimator = "ate") {
  if (estimator == "ate") {
    # 'Stabilized' ATE:
    return(ifelse(population$treatment == 1,
                  mean(population$treatment == 1) / population$propensityScore,
                  mean(population$treatment == 0) / (1 - population$propensityScore)))
  } else {
    # 'Stabilized' ATT:
    return(ifelse(population$treatment == 1,
                  mean(population$treatment == 1),
                  mean(population$treatment == 0) * population$propensityScore / (1 - population$propensityScore)))
  }
}

#' @title Trim Propensity Scores
#' @description Trims propensity scores by removing extreme values at both ends of the distribution.
#' @param propensityScore A numeric vector of propensity scores.
#' @param trimFraction Fraction of extreme values to trim (default 5%).
#' @return A vector of indices indicating which scores to keep.
#' @export
trimByPsQuantile <- function(propensityScore, trimFraction = 0.05) {
  # 0.05 cutoff on both sides
  cutoffUpper <- quantile(propensityScore, 1 - trimFraction)
  cutoffLower <- quantile(propensityScore, trimFraction)
  result <- which(propensityScore >= cutoffLower &  propensityScore <= cutoffUpper)
  return(result)
}

#' @title Compute Preference Score
#' @description Computes preference scores based on the propensity scores.
#' @param data A data frame containing propensity scores.
#' @param unfilteredData Optional dataset to compute proportions from.
#' @return A data frame with added preference scores.
#' @export
computePreferenceScore <- function(data, unfilteredData = NULL) {
  if (is.null(unfilteredData)) {
    proportion <- sum(data$treatment)/nrow(data)
  } else {
    proportion <- sum(unfilteredData$treatment)/nrow(unfilteredData)
  }
  propensityScore <- data$propensityScore
  propensityScore[propensityScore > 0.9999999] <- 0.9999999
  x <- exp(log(propensityScore/(1 - propensityScore)) - log(proportion/(1 - proportion)))
  data$preferenceScore <- x/(x + 1)
  return(data)
}

#' @title Stratify Population by Propensity Score
#' @description Assigns individuals to strata based on their propensity scores.
#' @param population A data frame containing row IDs, treatment assignments, and propensity scores.
#' @param numberOfStrata Number of strata to create.
#' @param stratificationColumns Additional columns to use for stratification.
#' @param baseSelection Defines which group is used to determine strata cutoffs ("all", "target", or "comparator").
#' @return A data frame with stratum assignments.
#' @export
stratifyByPs <- function(population, numberOfStrata = 5, stratificationColumns = c(), baseSelection = "all") {
  if (!("rowId" %in% colnames(population)))
    stop("Missing column rowId in population")
  if (!("treatment" %in% colnames(population)))
    stop("Missing column treatment in population")
  if (!("propensityScore" %in% colnames(population)))
    stop("Missing column propensityScore in population")
  ParallelLogger::logTrace("Stratifying by propensity score")
  if (nrow(population) == 0) {
    return(population)
  }
  baseSelection <- tolower(baseSelection)
  if (baseSelection == "all") {
    basePop <- population$propensityScore
  } else if (baseSelection == "target") {
    basePop <- population$propensityScore[population$treatment == 1]
  } else if (baseSelection == "target") {
    basePop <- population$propensityScore[population$treatment == 0]
  } else {
    stop(paste0("Unknown base selection: '", baseSelection, "'. Please choose 'all', 'target', or 'comparator'"))
  }
  if (length(basePop) == 0) {
    psStrata <- c()
  } else {
    psStrata <- unique(quantile(basePop, (1:(numberOfStrata - 1))/numberOfStrata))
  }
  attr(population, "strata") <- psStrata
  breaks <- unique(c(0, psStrata, 1))
  breaks[1] <- -1 # So 0 is included in the left-most stratum
  if (length(breaks) - 1 < numberOfStrata) {
    warning("Specified ", numberOfStrata, " strata, but only ", length(breaks) - 1, " could be created")
  }
  if (length(stratificationColumns) == 0) {
    if (length(breaks) - 1 == 1) {
      population$stratumId <- rep(1, nrow(population))
    } else {
      population$stratumId <- as.integer(as.character(cut(population$propensityScore,
                                                          breaks = breaks,
                                                          labels = 1:(length(breaks) - 1))))
    }
    return(population)
  } else {
    f <- function(subset, psStrata, numberOfStrata) {
      if (length(breaks) - 1 == 1) {
        subset$stratumId <- rep(1, nrow(subset))
      } else {
        subset$stratumId <- as.integer(as.character(cut(subset$propensityScore,
                                                        breaks = breaks,
                                                        labels = 1:(length(breaks) - 1))))
      }
      return(subset)
    }

    results <- population %>%
      group_by(across(all_of(stratificationColumns))) %>%
      group_split() %>%
      set_names(group_keys(.)) %>%
      map(~ f(.x, psStrata = psStrata, numberOfStrata = numberOfStrata))
    maxStratumId <- 0
    for (i in 1:length(results)) {
      if (nrow(results[[i]]) > 0) {
        if (maxStratumId != 0)
          results[[i]]$stratumId <- results[[i]]$stratumId + maxStratumId + 1
        maxStratumId <- max(results[[i]]$stratumId)
      }
    }
    result <- do.call(rbind, results)
    return(result)
  }
}
