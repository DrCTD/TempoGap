#' Evaluate Prediction Performance of TempoGap Models
#'
#' This function evaluates model performance using standard regression metrics,
#' including RMSE, MAE, and R-squared. If a recovery class is provided, it will
#' also compute classification metrics such as AUC and accuracy.
#'
#' @param predicted Numeric vector of predicted disease process time.
#' @param actual Numeric vector of actual disease process time.
#' @param recovery_class Optional. Factor or character vector of recovery class labels.
#' @param predicted_class Optional. Factor or character vector of predicted recovery class (for classification metrics).
#' @param plot Logical. Whether to show a scatterplot of predicted vs actual. Default is TRUE.
#'
#' @return A list of performance metrics. Optionally, a plot if \code{plot = TRUE}.
#' @importFrom pROC auc roc
#' @import ggplot2
#' @export
#'
#' @examples
#' predicted <- runif(20, 0, 72)
#' actual <- runif(20, 0, 72)
#' evaluate_prediction_performance(predicted, actual)
evaluate_prediction_performance <- function(predicted,
                                            actual,
                                            recovery_class = NULL,
                                            predicted_class = NULL,
                                            plot = TRUE){
  if (length(predicted) != length(actual)) {
    stop("Length of predicted and actual must be equal.")
  }
  rmse <- sqrt(mean((predicted - actual)^2))
  mae <- mean(abs(predicted - actual))
  r_squared <- 1 - sum((predicted - actual)^2) / sum((actual-mean(actual))^2)
  metrics <- list(
    RMSE = rmse,
    MAE = mae,
    R2 = r_squared
  )
  
  # Optional: classification metrics if recovery class is given
  if (!is.null(recovery_class) && !is.null(predicted_class)){
    if (length(recovery_class) != length(predicted_class)) {
      stop("Length of recovery_class and predicted_class must match.")
    }
    recovery_class <- factor(recovery_class)
    predicted_class <- factor(predicted_class, levels = levels(recovery_class))
    confusion <- table(Predicted = predicted_class, Actual = recovery_class)
    accuracy <- sum(diag(confusion)) / sum(confusion)
    metrics$ConfusionMatrix <- confusion
    metrics$Accuracy <- accuracy
    if (length(levels(recovery_class)) == 2){
      roc_obj <- pROC::roc(recovery_class, as.numeric(predicted_class))
      metrics$AUC <- pROC::auc(roc_obj)
    }
  }
  if (plot) {
    df <- data.frame(Predicted = predicted, Actual = actual)
    p <- ggplot(df, aes(x = Actual, y = Predicted)) +
      geom_point(alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "steelblue") +
      labs(
        title = "Predicted vs Actual Post-Injury Time",
        x = "Actual Time (hours)",
        y = "Predicted Time (hours)"
      ) +
      theme_minimal()
    print(p)
    metrics$Plot <- p
  }
  return(metrics)
}

