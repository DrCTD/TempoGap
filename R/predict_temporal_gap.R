#' Predict the disease process time and Compute Temporal Gap
#' This function uses a trained LASSO model to predict disease process time for gene expression data, and
#' optionally calculates the Temporal Gap between predicted and actual post-injury hours.
#'
#' @param model A trained \code{cv.glmnet} object return from \code{trian_temporal_model()}
#' @param expression_matrix A numeric matrix (features x samples) of gene/protein/lipid... expression values for prediction.
#' @param actual_time Optional. A numeric vector of actual post-injury hours (for computing TempoGap).
#' @param nboot Bootstrapping numbers. Default is 500.
#' @return A data frame with predicted time, actual time (if provided), and TempoGap (if applicable).
#' @examples
#' # set.seed(123)
#' # expr_train <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' # time_train <- runif(10, 0, 72)
#' # model <- train_temporal_model(expr_train, time_train)$model
#'
#' # expr_test <- matrix(rnorm(500), nrow = 100, ncol = 5)
#' # time_test <- runif(5, 0, 72)
#' # results <- predict_temporal_model(model, expr_test, time_test)
#' @export
predict_temporal_gap <- function(model, expression_matrix, actual_time = NULL, n_boot = 500){
  if (!inherits(model[[1]], "cv.glmnet")){
    stop("model must be a trained cv.glmnet object.")
  }
  if(!is.matrix(expression_matrix)){
    stop("expression_matrix must be a numeric matrix (features x samples)")
  }
  x <- t(expression_matrix)
  n_boot <- n_boot
  pred_matrix <- matrix(0, nrow = nrow(x), ncol = n_boot)
  for (i in 1:n_boot) {
    cv_model <- model[[i]]
    pred_i <- as.numeric(predict(model[[i]], newx = x, s = "lambda.min"))
    pred_matrix[,i] <- pred_i
    print(paste0(i, " times bootstrap sampling"))
  }
  pred_time_mean <- apply(pred_matrix, 1, mean)
  result <- data.frame(PredictedTime = pred_time_mean)
  if (!is.null(actual_time)) {
    if (length(actual_time) != length(pred_time_mean)) {
      stop("Length of actual_time must match number of samples in expression_matrix.")
    }
  }
  result$ActualTime <- actual_time
  loess_fit <- loess(PredictedTime ~ ActualTime, data = result)
  loess_pred <- predict(loess_fit, newdata = result)
  result$LoessAge <- loess_pred
  result$TempoGap <- pred_time_mean - loess_pred
  result$ScaledTempoGap <- scale(result$TempoGap)
  pred_matrix <- as.data.frame(pred_matrix)
  colnames(pred_matrix) <- paste0("Time", 1:n_boot)
  rownames(pred_matrix) <- colnames(expression_matrix)
  return(list(
    predicted_temporal_gap = result,
    predicted_matrix = pred_matrix
  ))
}












