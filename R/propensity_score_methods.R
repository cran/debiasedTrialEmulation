#' @title Estimate Treatment Effects using Propensity Score Matching
#' @description Computes effect estimates using propensity score matching to reduce confounding.
#' @param form A formula specifying the treatment assignment model.
#' @param data A dataset containing covariates and treatment assignment.
#' @param yvars A character vector of outcome variable names.
#' @param ncovars A character vector of negative control outcome variable names.
#' @param distance The method for estimating propensity scores ("glm").
#' @param outcome_measure The outcome measure to estimate: "RR" (Risk Ratio), "OR" (Odds Ratio), or "HR" (Hazard Ratio).
#' @return List of components
#' @rdname propensity_score_methods
#' @export
estEffect_matching <- function(form, data, yvars, ncovars, distance, outcome_measure){
  ### PS matching
  if (distance == "lasso"){
    m.out = matchit(form,
                    data=data,
                    caliper =0.2,
                    method="nearest",
                    distance="lasso",
                    ratio=5,
                    verbose=FALSE,
                    estimand="ATT",
                    distance.options = list(s="lambda.min"))
  }else{
    m.out = matchit(form,
                    data=data,
                    caliper =0.2,
                    method="nearest",
                    distance = "glm",
                    ratio=5,
                    verbose=FALSE,
                    estimand="ATT",
                    distance.options = list(s="lambda.min"))
  }
  m.match = match.data(m.out)
  m.match = as.data.frame(m.match)

  ### Outcome model
  df_rslt <- data.frame()
  df_rslt_nco <- data.frame()
  if (outcome_measure == "RR"){
    for (yvar in yvars){
      rslt <- get_RR_matching(m.match, yvar)
      df_rslt <- rbind(df_rslt, rslt)
    }
    for (ncovar in ncovars){
      rslt_nco <- get_RR_matching(m.match, ncovar)
      df_rslt_nco <- rbind(df_rslt_nco, rslt_nco)
    }
  }else if (outcome_measure == "OR"){
    for (yvar in yvars){
      rslt <- get_OR_matching(m.match, yvar)
      df_rslt <- rbind(df_rslt, rslt)
    }
    for (ncovar in ncovars){
      rslt_nco <- get_OR_matching(m.match, ncovar)
      df_rslt_nco <- rbind(df_rslt_nco, rslt_nco)
    }
  }else{
    for (yvar in yvars){
      rslt <- get_HR_matching(m.match, yvar)
      df_rslt <- rbind(df_rslt, rslt)
    }
    for (ncovar in ncovars){
      rslt_nco <- get_HR_matching(m.match, ncovar)
      df_rslt_nco <- rbind(df_rslt_nco, rslt_nco)
    }
  }

  df_rslt$ll <- df_rslt$logEst - 1.96*df_rslt$seLogEst
  df_rslt$ul <- df_rslt$logEst + 1.96*df_rslt$seLogEst

  df_rslt$names_outcome <- factor(df_rslt$names_outcome, levels = df_rslt$names_outcome)

  return(list(m.out = m.out, df_rslt = df_rslt, df_rslt_nco = df_rslt_nco))
}

#' @title Estimate Treatment Effects using Propensity Score Stratification
#' @description Computes effect estimates using propensity score stratification to adjust for confounding.
#' @return List of components
#' @rdname propensity_score_methods
#' @export
estEffect_stratification <- function(form, data, yvars, ncovars, distance, outcome_measure){
  ### PS stratification
  if (distance == "lasso"){
    Xmat <- grab_design_matrix(data = data, rhs_formula = form)
    Y <- data$treatment
    nfolds = 5
    foldid = sample(rep(seq(nfolds), length.out = length(Y)))
    Fit_ps_cv <- cv.glmnet(Xmat, Y, alpha = 1, family = "binomial", nfolds = nfolds, foldid = foldid) # cross validation to select lambda
    Fit_ps <- glmnet(Xmat, Y, alpha = 1, family = "binomial", lambda = Fit_ps_cv$lambda.min)
    propensityScore <- predict(Fit_ps, Xmat, type = "response")
  }else{
    Fit_ps <- glm(form, family = "binomial", data = data)
    propensityScore <- predict(Fit_ps, type = "response")
  }

  data$propensityScore <- propensityScore

  nstrata = 5
  rowId = c(1:length(data$treatment))
  data_Id <- cbind(rowId, data)
  stratifiedPop <- stratifyByPs(data_Id, numberOfStrata = nstrata)
  rm(data_Id)

  ### Outcome model
  df_rslt <- data.frame()
  df_rslt_nco <- data.frame()
  if (outcome_measure == "RR"){
    for (yvar in yvars){
      rslt <- get_RR_stratification(stratifiedPop, yvar)
      df_rslt <- rbind(df_rslt, rslt)
    }
    for (ncovar in ncovars){
      rslt_nco <- get_RR_stratification(stratifiedPop, ncovar)
      df_rslt_nco <- rbind(df_rslt_nco, rslt_nco)
    }
  }else if (outcome_measure == "OR"){
    for (yvar in yvars){
      rslt <- get_OR_stratification(stratifiedPop, yvar)
      df_rslt <- rbind(df_rslt, rslt)
    }
    for (ncovar in ncovars){
      rslt_nco <- get_OR_stratification(stratifiedPop, ncovar)
      df_rslt_nco <- rbind(df_rslt_nco, rslt_nco)
    }
  }else{
    for (yvar in yvars){
      rslt <- get_HR_stratification(stratifiedPop, yvar)
      df_rslt <- rbind(df_rslt, rslt)
    }
    for (ncovar in ncovars){
      rslt_nco <- get_HR_stratification(stratifiedPop, ncovar)
      df_rslt_nco <- rbind(df_rslt_nco, rslt_nco)
    }
  }

  df_rslt$ll <- df_rslt$logEst - 1.96*df_rslt$seLogEst
  df_rslt$ul <- df_rslt$logEst + 1.96*df_rslt$seLogEst

  df_rslt$names_outcome <- factor(df_rslt$names_outcome, levels = df_rslt$names_outcome)

  return(list(stratifiedPop = stratifiedPop, df_rslt = df_rslt, df_rslt_nco = df_rslt_nco))
}

#' @title Estimate Treatment Effects using Propensity Score Weighting
#' @description Computes effect estimates using propensity score weighting to balance covariates between treatment groups.
#' @return List of components
#' @rdname propensity_score_methods
#' @export
estEffect_weighting <- function(form, data, yvars, ncovars, distance, outcome_measure){
  ### PS weighting
  if (distance == "lasso"){
    Xmat <- grab_design_matrix(data = data, rhs_formula = form)
    Y <- data$treatment
    nfolds = 10
    foldid = sample(rep(seq(nfolds), length.out = length(Y)))
    Fit_ps_cv <- cv.glmnet(Xmat, Y, alpha = 1, family = "binomial", nfolds = nfolds, foldid = foldid) # cross validation to select lambda
    Fit_ps <- glmnet(Xmat, Y, alpha = 1, family = "binomial", lambda = Fit_ps_cv$lambda.min)
    propensityScore <- predict(Fit_ps, Xmat, type = "response")
  }else{
    Fit_ps <- glm(form, family = "binomial", data = data)
    propensityScore <- predict(Fit_ps, type = "response")
  }

  data$propensityScore <- propensityScore
  trimFraction <- 0.05
  propensityScore.trimmed <- trimByPsQuantile(propensityScore = data$propensityScore, trimFraction = trimFraction)

  data.trimmed <- data[propensityScore.trimmed != "trimmed", ]

  ### Outcome model
  df_rslt <- data.frame()
  df_rslt_nco <- data.frame()
  IPTW = computeWeights(data.trimmed)
  if (outcome_measure == "RR"){
    for (yvar in yvars){
      rslt <- get_RR_weighting(data.trimmed, yvar, IPTW)
      df_rslt <- rbind(df_rslt, rslt)
    }
    for (ncovar in ncovars){
      rslt_nco <- get_RR_weighting(data.trimmed, ncovar, IPTW)
      df_rslt_nco <- rbind(df_rslt_nco, rslt_nco)
    }
  }else if (outcome_measure == "OR"){
    for (yvar in yvars){
      rslt <- get_OR_weighting(data.trimmed, yvar, IPTW)
      df_rslt <- rbind(df_rslt, rslt)
    }
    for (ncovar in ncovars){
      rslt_nco <- get_OR_weighting(data.trimmed, ncovar, IPTW)
      df_rslt_nco <- rbind(df_rslt_nco, rslt_nco)
    }
  }else{
    for (yvar in yvars){
      rslt <- get_HR_weighting(data.trimmed, yvar, IPTW)
      df_rslt <- rbind(df_rslt, rslt)
    }
    for (ncovar in ncovars){
      rslt_nco <- get_HR_weighting(data.trimmed, ncovar, IPTW)
      df_rslt_nco <- rbind(df_rslt_nco, rslt_nco)
    }
  }

  df_rslt$ll <- df_rslt$logEst - 1.96*df_rslt$seLogEst
  df_rslt$ul <- df_rslt$logEst + 1.96*df_rslt$seLogEst
  df_rslt$names_outcome <- factor(df_rslt$names_outcome, levels = df_rslt$names_outcome)

  return(list(data.trimmed = data.trimmed, df_rslt = df_rslt, df_rslt_nco = df_rslt_nco))
}
