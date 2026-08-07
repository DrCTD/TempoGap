#' AUC and ROC Analysis for TempoGap Predictions
#'
#' Evaluates the classification performance of TempoGap using ROC curve and AUC.
#'
#' @param tempo_gap Numeric vector of TempoGap values.
#' @param recovery_class Vector of binary class labels (e.g., 0/1 or "complicated"/"uncomplicated").
#' @param positive_class The label to be considered as the "positive" class in ROC analysis.
#' @param plot Logical. Whether to return the ROC curve as a ggplot. Default: TRUE.
#' @param return_curve Logical. Whether to return the `roc` object from pROC. Default: FALSE.
#' @param title Optional plot title.
#'
#' @return A list with AUC, plot (if requested), and ROC object (if requested).
#' @export
auc_analysis_tempo_gap <- function(tempo_gap,
                                   recovery_class,
                                   positive_class = NULL,
                                   plot = TRUE,
                                   return_curve = FALSE,
                                   title = NULL){
  require(pROC)
  require(ggplot2)
  if (!is.factor(recovery_class)) {
    recovery_class <- as.factor(recovery_class)
  }
  if (length(levels(recovery_class)) != 2) {
    stop("Only binary classification supported for AUC. Provide a binary recovery_class.")
  }
  if (is.null(positive_class)) {
    positive_class <- levels(recovery_class)[2]
  }
  roc_obj <- roc(response = recovery_class,
                 predictor = tempo_gap,
                 levels = rev(levels(recovery_class)),
                 direction = "<",
                 quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  if (plot) {
    roc_df <- data.frame(
      specificity = rev(roc_obj$specificities),
      sensitivity = rev(roc_obj$sensitivities)
    )
    roc_plot <- ggplot2::ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
      ggplot2::geom_line(color = "#2C7BB6", size = 1.2) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
      ggplot2::labs(
        title = title %||% paste("ROC Curve (AUC =", round(auc_val, 3), ")"),
        x = "1 - Specificity",
        y = "Sensitivity"
      ) +
      ggplot2::theme_minimal(base_size = 14)
  } else {
    roc_plot <- NULL
  }
  return(list(
    auc = auc_val,
    plot = roc_plot,
    roc_obj = if (return_curve) roc_obj else NULL
  ))
}







