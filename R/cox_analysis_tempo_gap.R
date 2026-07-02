#' Cox Proportional Hazards Analysis on TempoGap
#'
#' Performs Cox regression using temporal gap as a predictor of recovery or outcome.
#'
#' @param tempo_gap Numeric vector of predicted - actual post-injury time.
#' @param time_to_event Numeric vector of time until event or censoring.
#' @param event Binary vector (1 = event occurred, 0 = censored).
#' @param covariates Optional data.frame of additional covariates.
#' @param plot Logical. Plot Kaplan-Meier survival stratified by TempoGap? Default TRUE.
#' @param cut_method How to group TempoGap values: "tertiles", "median", or "none". Default is "none".
#' @param return_model Logical. Return coxph object? Default FALSE.
#'
#' @return List with model summary and plot (if requested).
#' @export
cox_analysis_tempo_gap <- function(tempo_gap,
                                   time_to_event,
                                   event,
                                   covariates = NULL,
                                   plot = TRUE,
                                   cut_method = "none",
                                   return_model = FALSE){
  require(survival)
  require(survminer)
  require(tidyverse)
  surv_obj <- Surv(time_to_event, event)
  if (cut_method == "tertiles"){
    gap_group <- cut(tempo_gap, breaks = quantile(tempo_gap, probs = c(0, 1/3, 2/3, 1), na.rm = T),
                     include.lowest = TRUE, labels = c("Low", "Medium", "High"))
  } else if (cut_method == "median") {
    median_gap <- median(tempo_gap, na.rm = TRUE)
    gap_group <- ifelse(tempo_gap > median_gap, "High", "Low")
  } else {
    gap_group <- tempo_gap
  }
  df <- data.frame(surv_obj, tempo_gap = tempo_gap, gap_group = gap_group, time_to_event = time_to_event)
  if (!is.null(covariates)) {
    df <- cbind(df, covariates)
  }
  if (cut_method %in% c("tertiles", "median")) {
    cox_formula <- as.formula("surv_obj ~ gap_group")
  } else {
    cox_formula <- as.formula("surv_obj ~ tempo_gap")
  }
  if (!is.null(covariates)) {
    covars <- paste(colnames(covariates), collapse = " + ")
    cox_formula <- update(cox_formula, paste("~ . +", covars))
  }
  model <- coxph(cox_formula, data = df)
  summary_model <- summary(model)
  if (plot && cut_method != "none"){
    fit <- survfit(Surv(time_to_event, event) ~ gap_group, data = df)
    km_plot <- ggsurvplot(fit,
                          data = df,
                          pval = TRUE,
                          title = "Kaplan-Meier Curve by TempoGap Group",
                          xlab = "Time ro Recovery",
                          risk.table = TRUE)
  } else {
    km_plot <- NULL
  }
  result <- list(
    summary = summary_model,
    model = if(return_model) model else NULL,
    plot = km_plot
  )
  return(result)
}
