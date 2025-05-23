#' @title Estimate Odds Ratio (OR) after Propensity Score Matching
#' @description Computes the odds ratio for a binary outcome after applying propensity score matching.
#' @param data A dataset containing treatment assignment and outcome variables.
#' @param names_outcome A character vector of outcome variable names.
#' @return A data frame with the estimated log odds ratio, standard error, and p-value.
#' @rdname outcome_models
#' @export
get_OR_matching <- function(data, names_outcome){
  data_lm=data.frame(
    outcome = data[,names_outcome],
    treatment = data$treatment
  )
  s = summary(glm(outcome~treatment,data=data_lm,family = binomial(link = "logit")))
  rslt = s$coefficients[2,c(1,2,4)]

  return(data.frame(names_outcome = names_outcome, logEst = rslt[1], seLogEst = rslt[2], p = rslt[3]))
}

#' @title Estimate Risk Ratio (RR) after Propensity Score Matching
#' @description Computes the risk ratio for a binary outcome after applying propensity score matching.
#' @rdname outcome_models
#' @export
get_RR_matching <- function(data, names_outcome){
  data_lm=data.frame(
    outcome = data[,names_outcome],
    treatment = data$treatment
  )
  
  # Check if outcome is all zeros or all ones
  if (all(data_lm$outcome == 0) || all(data_lm$outcome == 1)) {
    warning(paste("Outcome", names_outcome, "is constant (all", ifelse(all(data_lm$outcome == 0), "zeros", "ones"), "). Cannot estimate risk ratio."))
    return(data.frame(names_outcome = names_outcome, logEst = NA, seLogEst = NA, p = NA))
  }
  
  # Try to fit the model
  tryCatch({
    s = summary(glm(outcome~treatment, data=data_lm, family = poisson(link = "log")))
    rslt = s$coefficients[2,c(1,2,4)]
    
    # Check if results are valid
    if (any(is.na(rslt))) {
      warning(paste("Model for outcome", names_outcome, "did not converge. Returning NA values."))
      return(data.frame(names_outcome = names_outcome, logEst = NA, seLogEst = NA, p = NA))
    }
    
    return(data.frame(names_outcome = names_outcome, logEst = rslt[1], seLogEst = rslt[2], p = rslt[3]))
  }, error = function(e) {
    warning(paste("Error fitting model for outcome", names_outcome, ":", e$message))
    return(data.frame(names_outcome = names_outcome, logEst = NA, seLogEst = NA, p = NA))
  })
}

#' @title Estimate Hazard Ratio (HR) after Propensity Score Matching
#' @description Computes the hazard ratio for a time-to-event outcome after propensity score matching.
#' @rdname outcome_models
#' @export
get_HR_matching <- function(data, names_outcome){
  names_time <- gsub("event_", "time_", names_outcome)
  data_lm=data.frame(
    outcome = data[,names_outcome],
    treatment = data$treatment,
    time = data[,names_time]
  )

  cox_fit <- coxph(Surv(time, outcome) ~ treatment,
                   data = data_lm)

  logEst <- coef(cox_fit)
  seLogEst <- sqrt(diag(vcov(cox_fit)))
  p <- 2 * (1 - pnorm(abs(logEst/seLogEst)))

  return(data.frame(names_outcome = names_outcome, logEst = logEst, seLogEst = seLogEst, p = p))
}

#' @title Estimate Treatment Effects after Propensity Score Stratification
#' @description Computes OR, RR, and HR for a binary or time-to-event outcome after stratification.
#' @rdname outcome_models
#' @export
get_OR_stratification <- function(data, names_outcome){
  data_lm=data.frame(
    outcome = data[,names_outcome],
    treatment = data$treatment,
    stratumID = data$stratumId
  )
  s = summary(glm(outcome~treatment+strata(stratumID),data=data_lm,family = binomial(link = "logit")))
  rslt = s$coefficients[2,c(1,2,4)]

  return(data.frame(names_outcome = names_outcome, logEst = rslt[1], seLogEst = rslt[2], p = rslt[3]))
}

#' @title Estimate Risk Ratio (RR) after Propensity Score Stratification
#' @description Computes the risk ratio for a binary outcome after applying propensity score stratification.
#' @rdname outcome_models
#' @export
get_RR_stratification <- function(data, names_outcome){
  data_lm=data.frame(
    outcome = data[,names_outcome],
    treatment = data$treatment,
    stratumID = data$stratumId
  )
  s = summary(glm(outcome~treatment+strata(stratumID),data=data_lm,family = poisson(link = "log")))
  rslt=s$coefficients[2,c(1,2,4)]

  return(data.frame(names_outcome = names_outcome, logEst = rslt[1], seLogEst = rslt[2], p = rslt[3]))
}

#' @title Estimate Hazard Ratio (HR) after Propensity Score Stratification
#' @description Computes the hazard ratio for a time-to-event outcome after propensity score stratification.
#' @rdname outcome_models
#' @export
get_HR_stratification <- function(data, names_outcome){
  names_time <- gsub("event_", "time_", names_outcome)
  data_lm=data.frame(
    outcome = data[,names_outcome],
    treatment = data$treatment,
    time = data[,names_time],
    stratumID = data$stratumId
  )

  cox_fit <- coxph(Surv(time, outcome) ~ treatment+strata(stratumID),
                   data = data_lm)

  logEst <- coef(cox_fit)
  seLogEst <- sqrt(diag(vcov(cox_fit)))
  p <- 2 * (1 - pnorm(abs(logEst/seLogEst)))

  return(data.frame(names_outcome = names_outcome, logEst = logEst, seLogEst = seLogEst, p = p))
}

#' @title Estimate Odds Ratio (OR) after Propensity Score Weighting
#' @description Computes the odds ratio for a binary outcome after applying propensity score weighting.
#' @param IPTW A numeric vector of inverse probability of treatment weights.
#' @rdname outcome_models
#' @export
get_OR_weighting <- function(data, names_outcome, IPTW){
  data_lm=data.frame(
    outcome = data[,names_outcome],
    treatment = data$treatment,
    IPTW = IPTW
  )

  s = summary(glm(outcome~treatment,data=data_lm,family = binomial(link = "logit"), weights = IPTW))
  rslt = s$coefficients[2,c(1,2,4)]

  return(data.frame(names_outcome = names_outcome, logEst = rslt[1], seLogEst = rslt[2], p = rslt[3]))
}

#' @title Estimate Risk Ratio (RR) after Propensity Score Weighting
#' @description Computes the risk ratio for a binary outcome after applying propensity score weighting.
#' @rdname outcome_models
#' @export
get_RR_weighting <- function(data, names_outcome, IPTW){
  data_lm=data.frame(
    outcome = data[,names_outcome],
    treatment = data$treatment,
    IPTW = IPTW
  )
  s = summary(glm(outcome~treatment,data=data_lm,family = poisson(link = "log"), weights = IPTW))
  rslt=s$coefficients[2,c(1,2,4)]

  return(data.frame(names_outcome = names_outcome, logEst = rslt[1], seLogEst = rslt[2], p = rslt[3]))
}

#' @title Estimate Hazard Ratio (HR) after Propensity Score Weighting
#' @description Computes the hazard ratio for a time-to-event outcome after applying propensity score weighting.
#' @rdname outcome_models
#' @export
get_HR_weighting <- function(data, names_outcome, IPTW){
  names_time <- gsub("event_", "time_", names_outcome)
  data_lm=data.frame(
    outcome = data[,names_outcome],
    treatment = data$treatment,
    time = data[,names_time],
    IPTW = IPTW
  )

  cox_fit <- coxph(Surv(time, outcome) ~ treatment,
                   data = data_lm,
                   weights = IPTW)

  logEst <- coef(cox_fit)
  seLogEst <- sqrt(diag(vcov(cox_fit)))
  p <- 2 * (1 - pnorm(abs(logEst/seLogEst)))

  return(data.frame(names_outcome = names_outcome, logEst = logEst, seLogEst = seLogEst, p = p))
}
