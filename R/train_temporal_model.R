#' Train Temporal Model (patient-level resampling)
#' Purpose: Trains the LASSO model to predict recovery-process timing from
#' expression data, using PATIENT-level bootstrap resampling and PATIENT-
#' grouped cross-validation folds, so that samples from the same patient
#' never end up split across the bootstrap train set / CV train-fold /
#' CV validation-fold in a way that would leak information.
#'
#' @param expression_matrix A numeric matrix (features x samples) of gene/protein/lipid... expression values.
#' @param time_vector A numeric vector of disease process timing data for each sample. Note: it must have the same units.
#' @param patient_id A vector (length = ncol(expression_matrix)) giving the patient/subject
#'   identifier for each sample/column. Samples sharing an ID are treated as
#'   belonging to the same patient and are always kept together (same
#'   bootstrap draw, same CV fold).
#' @param cv_folds Integer, number of folds for cross-validation. Default is 5.
#' @param alpha Elastic net mixing parameter. Default is 1 (LASSO). If alpha is a vector of continuous numbers, hyperparameter optimization of alpha will be performed based on the vector.
#' @param seed Optional, a random seed for reproducibility.
#' @param n_boot Bootstrapping numbers. Default is 500.
#' @import glmnet
#' @export
#'
#' @examples
#' # set.seed(42)
#' # data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
#' # time <- runif(10, 0, 72)
#' # patient <- rep(1:5, each = 2)  # e.g. 2 samples per patient
#' # model <- train_temporal_model(data, time, patient_id = patient)
train_temporal_model <- function(expression_matrix, time_vector, patient_id,
                                 cv_folds = 5, alpha = 1, seed = NULL, n_boot = 500) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # ---- Check inputs --------------------------------------------------------
  if (!is.matrix(expression_matrix)) {
    stop("expression_matrix must be a numeric matrix (features x samples)")
  }
  if (length(time_vector) != ncol(expression_matrix)) {
    stop("Length of time_vector must match number of samples (columns) in expression_matrix.")
  }
  if (length(patient_id) != ncol(expression_matrix)) {
    stop("Length of patient_id must match number of samples (columns) in expression_matrix.")
  }
  
  x <- t(expression_matrix)
  y <- time_vector
  unique_patients <- unique(patient_id)
  
  # ---- Helper: assign CV folds by patient, not by row -----------------------
  # All samples belonging to the same patient (even duplicated rows coming
  # from the same patient being drawn more than once in a bootstrap) are
  # always assigned to the SAME fold, so cv.glmnet never trains and
  # validates on the same patient at once.
  make_patient_foldid <- function(patient_ids_for_rows, k) {
    ids <- unique(patient_ids_for_rows)
    fold_lookup <- sample(rep(seq_len(k), length.out = length(ids)))
    names(fold_lookup) <- ids
    unname(fold_lookup[as.character(patient_ids_for_rows)])
  }
  
  coef_matrix <- matrix(0, nrow = ncol(x) + 1, ncol = n_boot) # +1 for intercept
  cv_model_list <- list()
  lambda_min_list <- list()
  
  for (i in 1:n_boot) {
    
    # ---- Patient-level bootstrap draw ---------------------------------------
    # Resample PATIENTS with replacement, then pull in every sample belonging
    # to each drawn patient. A patient drawn twice contributes its full set
    # of samples twice.
    boot_patients <- sample(unique_patients, replace = TRUE)
    sample_idx <- unlist(lapply(boot_patients, function(p) which(patient_id == p)))
    
    x_boot <- x[sample_idx, , drop = FALSE]
    y_boot <- y[sample_idx]
    patient_id_boot <- patient_id[sample_idx]
    
    # Patient-grouped fold assignment for this bootstrap draw
    foldid_boot <- make_patient_foldid(patient_id_boot, cv_folds)
    
    # ---- Hyperparameter (alpha) search --------------------------------------
    alpha_grid <- alpha
    cv_errors <- numeric(length = length(alpha_grid))
    lambda_min_values <- numeric(length = length(alpha_grid))
    cv_fits_by_alpha <- list()
    
    for (n in seq_along(alpha_grid)) {
      alpha_val <- alpha_grid[n]
      # patient-grouped CV instead of a plain random nfolds split
      cv_fit <- cv.glmnet(x_boot, y_boot, alpha = alpha_val, foldid = foldid_boot)
      cv_errors[n] <- min(cv_fit$cvm)
      lambda_min_values[n] <- cv_fit$lambda.min
      cv_fits_by_alpha[[n]] <- cv_fit
    }
    
    tune_results <- data.frame(
      alpha = alpha_grid,
      lambda_min = lambda_min_values,
      cvm_min = cv_errors
    )
    best_idx <- which.min(tune_results$cvm_min)
    best_alpha <- tune_results$alpha[best_idx]
    best_lambda <- tune_results$lambda_min[best_idx]
    best_cv_fit <- cv_fits_by_alpha[[best_idx]]
    
    print(tune_results)
    cat("\nBest alpha:", best_alpha, "\nBest lambda:", best_lambda, "\n")
    
    # BUGFIX vs. original: coefficients must come from the *best* alpha's
    # model (best_cv_fit), not from whichever alpha happened to run last
    # in the loop above (the original code mistakenly used `cv_fit` here).
    coef_i <- as.numeric(coef(best_cv_fit, s = "lambda.min"))
    coef_matrix[, i] <- coef_i
    
    cv_model_list[[i]] <- best_cv_fit
    lambda_min_list[[i]] <- best_cv_fit$lambda.min
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
