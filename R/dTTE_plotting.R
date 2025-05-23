#' @title Plot Standardized Mean Differences (SMD) for Matching
#' @description Generates a plot of standardized mean differences before and after propensity score matching.
#' @param m.out The output from `MatchIt`, containing matched data.
#' @return A `ggplot2` object showing balance improvement after matching.
#' @rdname dTTE_plotting
#' @export
plot_SMD_matching <- function(m.out){
  s3=summary(m.out)

  smd_before <- s3[["sum.all"]][, "Std. Mean Diff."]
  smd_after <- s3[["sum.matched"]][, "Std. Mean Diff."]

  df_smd <- data.frame(smd_before,smd_after)
  df_smd = data.frame(row.names(df_smd), df_smd)
  df_smd = df_smd[-1,] # delete first row "distance"
  rownames(df_smd) <- NULL

  colnames(df_smd) = c("covariateName", "beforeMatchingStdDiff", "afterMatchingStdDiff")

  # Largest 20 imbalance
  largest.imbalance.plot <- plotCovariateBalanceOfTopVariables(df_smd, n = 20,
                                                               beforeLabel = "Before matching",
                                                               afterLabel = "After matching")
  return(largest.imbalance.plot)
}

#' @title Plot Equipoise for Matching
#' @description Generates a plot showing the distribution of preference scores to assess equipoise after matching.
#' @param data The dataset containing treatment and propensity scores.
#' @param m.out The output from `MatchIt`, containing matched data.
#' @rdname dTTE_plotting
#' @export
plot_Equipoise_matching <- function(data, m.out){
  data$propensityScore <- m.out$distance
  data <- computePreferenceScore(data)

  equipoise.plot <- plotPs(data, scale = "preference",
                           showCountsLabel = TRUE, showAucLabel = TRUE, showEquiposeLabel = TRUE)

  return(equipoise.plot)
}

#' @title Plot Standardized Mean Differences (SMD) for Stratification
#' @description Generates an SMD plot to compare balance before and after stratification.
#' @param stratifiedPop The dataset containing stratified propensity scores.
#' @param xvars The covariate names to assess balance.
#' @rdname dTTE_plotting
#' @export
plot_SMD_stratification <- function(stratifiedPop, xvars){
  smd_before <- GetSMD(stratifiedPop[, xvars], stratifiedPop$treatment)

  stratifiedPop <- merge(stratifiedPop, Compute_weight(stratifiedPop), by = "rowId", all.x = TRUE)
  smd_after <- GetSMD(stratifiedPop[, xvars], stratifiedPop$treatment, stratifiedPop$weight)

  balance_table = cbind(smd_before, smd_after[,2])
  colnames(balance_table) = c("covariateName", "beforeMatchingStdDiff", "afterMatchingStdDiff")

  largest.imbalance.plot <- plotCovariateBalanceOfTopVariables(balance_table, n = 20,
                                                               beforeLabel = "Before stratification",
                                                               afterLabel = "After stratification")

  return(largest.imbalance.plot)
}

#' @title Plot Equipoise for Stratification
#' @description Generates a plot showing the distribution of preference scores to assess equipoise after stratification.
#' @param data The dataset containing treatment and propensity scores.
#' @rdname dTTE_plotting
#' @export
plot_Equipoise_stratification <- function(data){
  data <- computePreferenceScore(data)
  equipoise.plot <- plotPs(data, scale = "preference",
                           showCountsLabel = TRUE, showAucLabel = TRUE, showEquiposeLabel = TRUE)

  return(equipoise.plot)
}

#' @title Plot Standardized Mean Differences (SMD) for Weighting
#' @description Generates an SMD plot to compare balance before and after propensity score weighting.
#' @param data The dataset containing treatment and propensity scores.
#' @rdname dTTE_plotting
#' @export
plot_SMD_weighting <- function(data){
  smd_before <- GetSMD(data, data$treatment)

  IPTW = computeWeights(data)
  smd_after <- GetSMD(data, data$treatment, IPTW)

  balance_table = cbind(smd_before, smd_after[,2])
  colnames(balance_table) = c("covariateName", "beforeMatchingStdDiff", "afterMatchingStdDiff")

  largest.imbalance.plot <- plotCovariateBalanceOfTopVariables(balance_table, n = 20,
                                                               beforeLabel = "Before weighting",
                                                               afterLabel = "After weighting")

  return(largest.imbalance.plot)
}

#' @title Plot Equipoise for Weighting
#' @description Generates a plot showing the distribution of preference scores to assess equipoise after weighting.
#' @param data The dataset containing treatment and propensity scores.
#' @rdname dTTE_plotting
#' @export
plot_Equipoise_weighting <- function(data){
  data <- computePreferenceScore(data)

  equipoise.plot <- plotPs(data, scale = "preference",
                           showCountsLabel = TRUE, showAucLabel = TRUE, showEquiposeLabel = TRUE)

  return(equipoise.plot)
}

#' @title Plot Propensity Score Distributions
#' @description Creates a plot of propensity score distributions using density or histogram visualization.
#' @param data A dataset containing treatment and propensity scores.
#' @param scale The scale to use: "preference" or "propensity".
#' @param type The type of plot: "density", "histogramCount", or "histogramProportion".
#' @param binWidth The bin width for histograms (default = 0.05).
#' @param targetLabel Label for the treated group.
#' @param comparatorLabel Label for the control group.
#' @param showCountsLabel Logical; whether to show sample counts.
#' @param showAucLabel Logical; whether to show AUC.
#' @param showEquiposeLabel Logical; whether to indicate equipoise range.
#' @param equipoiseBounds A numeric vector of two values defining the equipoise range (default = c(0.3, 0.7)).
#' @param unitOfAnalysis Unit label for counts (e.g., "subjects").
#' @param title Optional title for the plot.
#' @param fileName Optional file name to save the plot.
#' @rdname dTTE_plotting
#' @export
plotPs <- function(data,
                   unfilteredData = NULL,
                   scale = "preference",
                   type = "density",
                   binWidth = 0.05,
                   targetLabel = "Target",
                   comparatorLabel = "Comparator",
                   showCountsLabel = FALSE,
                   showAucLabel = FALSE,
                   showEquiposeLabel = FALSE,
                   equipoiseBounds = c(0.3, 0.7),
                   unitOfAnalysis = "subjects",
                   title = NULL,
                   fileName = NULL) {
  if (!("treatment" %in% colnames(data)))
    stop("Missing column treatment in data")
  if (!("propensityScore" %in% colnames(data)))
    stop("Missing column propensityScore in data")
  if (!is.null(unfilteredData)) {
    if (!("treatment" %in% colnames(unfilteredData)))
      stop("Missing column treatment in unfilteredData")
    if (!("propensityScore" %in% colnames(unfilteredData)))
      stop("Missing column propensityScore in unfilteredData")
  }
  if (type != "density" && type != "histogram" && type != "histogramCount" && type != "histogramProportion")
    stop(paste("Unknown type '", type, "', please choose either 'density', 'histogram', 'histogramCount', or 'histogramProportion'"),
         sep = "")
  if (type == "histogram")
    type <- "histogramCount"
  if (scale != "propensity" && scale != "preference")
    stop(paste("Unknown scale '", scale, "', please choose either 'propensity' or 'preference'"),
         sep = "")
  targetLabel <- as.character(targetLabel)
  comparatorLabel <- as.character(comparatorLabel)

  if (scale == "preference") {
    #data <- computePreferenceScore(data, unfilteredData)
    data$score <- data$preferenceScore
    label <- "Preference score"
  } else {
    data$score <- data$propensityScore
    label <- "Propensity score"
  }
  if (showAucLabel || showCountsLabel || showEquiposeLabel) {
    yMultiplier <- 1.25
  } else {
    yMultiplier <- 1
  }
  if (type == "density") {
    d1 <- density(data$score[data$treatment == 1], from = 0, to = 1, n = 200)
    d0 <- density(data$score[data$treatment == 0], from = 0, to = 1, n = 200)
    d <- data.frame(x = c(d1$x, d0$x), y = c(d1$y, d0$y), treatment = c(rep(as.character(targetLabel), length(d1$x)),
                                                                        rep(as.character(comparatorLabel), length(d0$x))))
    d$treatment <- factor(d$treatment, levels = c(targetLabel, comparatorLabel))
    plot <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_density(stat = "identity", ggplot2::aes(color = .data$treatment, group = .data$treatment, fill = .data$treatment)) +
      ggplot2::scale_fill_manual(values = c(rgb(0.8, 0, 0, alpha = 0.5),
                                            rgb(0, 0, 0.8, alpha = 0.5))) +
      ggplot2::scale_color_manual(values = c(rgb(0.8, 0, 0, alpha = 0.5),
                                             rgb(0, 0, 0.8, alpha = 0.5))) +
      ggplot2::scale_x_continuous(label, limits = c(0, 1)) +
      ggplot2::scale_y_continuous("Density", limits = c(0, max(d$y)*yMultiplier)) +
      ggplot2::theme(legend.title = ggplot2::element_blank(),
                     legend.position = "top",
                     legend.text = ggplot2::element_text(margin = ggplot2::margin(0, 0.5, 0, 0.1, "cm")))
    if (!is.null(attr(data, "strata"))) {
      strata <- data.frame(propensityScore = attr(data, "strata"))
      if (scale == "preference") {
        if (is.null(unfilteredData)) {
          strata <- computePreferenceScore(strata, data)
        } else {
          strata <- computePreferenceScore(strata, unfilteredData)
        }
        strata$score <- strata$preferenceScore
      } else {
        strata$score <- strata$propensityScore
      }
      plot <- plot +
        ggplot2::geom_vline(xintercept = strata$score, color = rgb(0, 0, 0, alpha = 0.5))
    }

  } else {
    x <- seq(from = 0, to = 1, by = binWidth)
    d1 <- data.frame(xmin = cut(data$score[data$treatment == 1], x, labels = x[1:(length(x) - 1)]), y = 1)
    d1 <- aggregate(y ~   xmin, d1, sum)
    d1$xmin <- as.numeric(as.character(d1$xmin))
    d0 <- data.frame(xmin = cut(data$score[data$treatment == 0], x, labels = x[1:(length(x) - 1)]), y = 1)
    d0 <- aggregate(y ~   xmin, d0, sum)
    d0$xmin <- as.numeric(as.character(d0$xmin))
    d <- data.frame(xmin = c(d1$xmin, d0$xmin), y = c(d1$y, d0$y), treatment = c(rep(as.character(targetLabel), nrow(d1)),
                                                                                 rep(as.character(comparatorLabel), nrow(d0))))
    d$xmax <- d$xmin + binWidth
    d$treatment <- factor(d$treatment, levels = c(targetLabel, comparatorLabel))
    yAxisScale <- "Number"
    if (type == "histogramProportion") {
      d$y <- d$y / sum(d$y)
      yAxisScale <- "Proportion"
    }
    plot <- ggplot2::ggplot(d, ggplot2::aes(x = .data$xmin)) +
      ggplot2::geom_rect(ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax, ymin = 0, ymax = .data$y, color = .data$treatment, group = .data$treatment, fill = .data$treatment)) +
      ggplot2::scale_fill_manual(values = c(rgb(0.8, 0, 0, alpha = 0.5),
                                            rgb(0, 0, 0.8, alpha = 0.5))) +
      ggplot2::scale_color_manual(values = c(rgb(0.8, 0, 0, alpha = 0.5),
                                             rgb(0, 0, 0.8, alpha = 0.5))) +
      ggplot2::scale_x_continuous(label, limits = c(0, 1)) +
      ggplot2::scale_y_continuous(paste(yAxisScale, "of", unitOfAnalysis), limits = c(0, max(d$y)*1.25)) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "top")
  }
  if (showAucLabel || showCountsLabel || showEquiposeLabel) {
    labelsLeft <- c()
    labelsRight <- c()
    if (showCountsLabel) {
      labelsLeft <- c(labelsLeft, sprintf("%s: %s %s", targetLabel, format(sum(data$treatment == 1), big.mark = ",", scientific = FALSE), unitOfAnalysis))
      labelsLeft <- c(labelsLeft, sprintf("%s: %s %s", comparatorLabel, format(sum(data$treatment == 0), big.mark = ",", scientific = FALSE), unitOfAnalysis))
    }

    #if (showAucLabel) {
    #  auc <- CohortMethod::computePsAuc(data, confidenceIntervals = FALSE)
    #  labelsRight <- c(labelsRight, sprintf("AUC:\t\t%0.2f", auc))
    #}
    if (showEquiposeLabel) {
      if (is.null(data$preferenceScore)) {
        data <- computePreferenceScore(data, unfilteredData)
      }
      equipoise <- mean(data$preferenceScore >= equipoiseBounds[1] & data$preferenceScore <= equipoiseBounds[2])
      labelsRight <- c(labelsRight, sprintf("%2.1f%% is in equipoise", equipoise*100))
    }
    # maxY <- ggplot2::ggplot_build(plot)$layout$panel_ranges[[1]]$y.range[2]
    if (length(labelsLeft) > 0) {
      dummy <- data.frame(text = paste(labelsLeft, collapse = "\n"))
      plot <- plot + ggplot2::geom_label(x = 0, y = max(d$y) * 1.24, hjust = "left", vjust = "top", alpha = 0.8, ggplot2::aes(label = text), data = dummy, size = 3.5)
    }
    if (length(labelsRight) > 0) {
      dummy <- data.frame(text = paste(labelsRight, collapse = "\n"))
      plot <- plot + ggplot2::geom_label(x = 1, y =  max(d$y) * 1.24, hjust = "right", vjust = "top", alpha = 0.8, ggplot2::aes(label = text), data = dummy, size = 3.5)
    }
  }
  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 5, height = 3.5, dpi = 400)
  return(plot)
}


#' @title Plot Covariate Balance of Top Variables
#' @description Creates a plot showing the top covariates with the largest standardized mean differences before and after matching.
#' @param balance A data frame containing standardized mean differences.
#' @param n Number of top covariates to display.
#' @param beforeLabel Label for pre-matching imbalance.
#' @param afterLabel Label for post-matching balance.
#' @param title Optional title for the plot.
#' @param fileName Optional file name to save the plot.
#' @rdname dTTE_plotting
#' @export
plotCovariateBalanceOfTopVariables <- function(balance,
                                               n = 20,
                                               maxNameWidth = 100,
                                               title = NULL,
                                               fileName = NULL,
                                               beforeLabel = "before matching",
                                               afterLabel = "after matching") {
  n <- min(n, nrow(balance))
  beforeLabel <- as.character(beforeLabel)
  afterLabel <- as.character(afterLabel)
  topBefore <- balance[order(-abs(balance$beforeMatchingStdDiff)), ]
  topBefore <- topBefore[1:n, ]
  topBefore$facet <- paste("Top", n, beforeLabel)
  topAfter <- balance[order(-abs(balance$afterMatchingStdDiff)), ]
  topAfter <- topAfter[1:n, ]
  topAfter$facet <- paste("Top", n, afterLabel)
  filtered <- rbind(topBefore, topAfter)

  data <- dplyr::tibble(covariate = rep(filtered$covariateName, 2),
                        difference = c(filtered$beforeMatchingStdDiff, filtered$afterMatchingStdDiff),
                        group = rep(c(beforeLabel, afterLabel), each = nrow(filtered)),
                        facet = rep(filtered$facet, 2),
                        rowId = rep(nrow(filtered):1, 2))
  #filtered$covariateName <- .truncRight(as.character(filtered$covariateName), maxNameWidth)
  data$facet <- factor(data$facet, levels = c(paste("Top", n, beforeLabel), paste("Top", n, afterLabel)))
  data$group <- factor(data$group, levels = c(beforeLabel, afterLabel))
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$difference,
                                             y = .data$rowId,
                                             color = .data$group,
                                             group = .data$group,
                                             fill = .data$group,
                                             shape = .data$group)) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_vline(xintercept = 0.1, lty = 2) +
    ggplot2::geom_vline(xintercept = -0.1, lty = 2) +
    ggplot2::scale_fill_manual(values = c(rgb(0.8, 0, 0, alpha = 0.5),
                                          rgb(0, 0, 0.8, alpha = 0.5))) +
    ggplot2::scale_color_manual(values = c(rgb(0.8, 0, 0, alpha = 0.5),
                                           rgb(0, 0, 0.8, alpha = 0.5))) +
    ggplot2::scale_x_continuous("Standardized difference of mean") +
    ggplot2::scale_y_continuous(breaks = nrow(filtered):1, labels = filtered$covariateName) +
    ggplot2::facet_grid(facet ~ ., scales = "free", space = "free") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7),
                   axis.title.y = ggplot2::element_blank(),
                   legend.position = "top",
                   legend.direction = "vertical",
                   legend.title = ggplot2::element_blank())
  if (!is.null(title)) {
    plot <- plot + ggplot2::ggtitle(title)
  }
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 10, height = max(2 + n * 0.2, 5), dpi = 400)
  return(plot)
}
