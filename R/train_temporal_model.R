#' Train Temporal Model
#' Purpose: Trains the LASSO mdoel to predict recovery class of disease based on expression data
#' @param expression_matrix A numeric matrix (features x samples) of gene/protein/lipid... expression values.
#' @param time_vextor A numeric vextor of disease process timing data for each sample. Note: It must have same units.
#' @param cv_folds Integer, number of folds for cross-validation. Default is 5.
#' @param alpha Elastic net mixing parameter. Default is 1 (LASSO). If alpha is a vector of continuous numbers, hyperparameter optimization of alpha will be performed based on the vector.
#' @param seed Optional, a random seed for reproducibility.
#' @param nboot Bootstrapping numbers. Default is 500.
#' @import glmnet
#' @export
#'
#' @examples
#' # set.seed(42)
#' # data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' # time <- runif(10, 0, 72)
#' # model <- train_temporal_model(data, time)

train_temporal_model <- function(expression_matrix, time_vector, cv_folds = 5, alpha = 1, seed = NULL, n_boot = 500){
  if (!is.null(seed)) set.seed(seed)
  # Check inputs
  if (!is.matrix(expression_matrix)) {
    stop("expression_matrix must be a numeric matrix (features x samples)")
  }
  if (length(time_vector) != ncol(expression_matrix)){
    stop("Length of time_vector must match number of samples (columns) in expression_matrix.")
  }
  x <- t(expression_matrix)
  y <- time_vector
  coef_matrix <- matrix(0, nrow = ncol(x) + 1, ncol = n_boot) # +1 for intercept
  cv_model_list <- list()
  lambda_min_list <- list()
  m <- 1
  for (i in 1:n_boot) {
    sample_idx <- sample(1:nrow(x), replace = T)
    x_boot <- x[sample_idx, ]
    y_boot <- y[sample_idx]
    # Define grid of alpha values
    alpha_grid <- alpha
    cv_errors <- numeric(length = length(alpha_grid))
    lambda_min_values <- numeric(length = length(alpha_grid))
    cv_models <- list()
    # Loop over alpha values
    for (n in seq_along(alpha_grid)) {
      alpha_val <- alpha_grid[n]
      # 5-fold CV for each alpha
      cv_fit <- cv.glmnet(x, y, alpha = alpha_val, nfolds = cv_folds)
      # Save results
      cv_errors[n] <- min(cv_fit$cvm) # lowest mean CV error
      lambda_min_values[n] <- cv_fit$lambda.min
      cv_models[[n]] <- cv_fit
    }
    # Combine results
    tune_results <- data.frame(
      alpha = alpha_grid,
      lambda_min = lambda_min_values,
      cvm_min = cv_errors
    )
    best_idx <- which.min(tune_results$cvm_min)
    best_alpha <- tune_results$alpha[best_idx]
    best_lambda <- tune_results$lambda_min[best_idx]
    print(tune_results)
    cv_models <- cv_models[[best_idx]]
    cat("\nBest alpha:", best_alpha, "\nBest lambda:", best_lambda, "\n")
    coef_i <- as.numeric(coef(cv_fit, s = "lambda.min"))
    coef_matrix[,i] <- coef_i
    cv_model_list[[m]] <- cv_models
    lambda_min_list[[m]] <- cv_models$lambda.min
    m <- m + 1
  }
  coef_matrix <- as.data.frame(coef_matrix)
  colnames(coef_matrix) <- paste0("Times", 1:n_boot)
  rownames(coef_matrix) <- c("Intercept", rownames(expression_matrix))
  return(list(
    model = cv_model_list,
    lambda_min = lambda_min_list,
    coef_matrix = coef_matrix
  ))
}
















