#' Plot Bootstrap LASSO Coefficient Stability
#'
#' This function visulizas the frequency and average effect size of gene selected in an ensemble of bootstrapped LASSO models.
#'
#' @param coef_matrix The coefficient matrix which produced by \code{train_temporal_model} function.
#' @param dot_col The point colors of the dot plot. Default is lightblue.
#' @param top_n Integer, number of top features to display. Default is 20.
#' @param return_data Logical, whether to return the summary table. Default is FALSE.
#' @param n_boot Bootstrapping numbers. Default is 500.
#'
#' @return A \code{ggplot2} dot plot showing selection frequency and mean coefficients.
#'         If \code{return = TRUE}, returns a summary data frame instead.
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' # boot_models <- bootstrap_temporal_models(...)
#' # plot_bootstrap_coefficients(boot_models$models)

plot_bootstrap_coefficients <- function(coef_matrix,
                                        dot_col = "lightblue",
                                        top_n = 20,
                                        return_data = FALSE,
                                        n_boot = 500){
  coef_matrix <- coef_matrix[-1, , drop = FALSE]
  freq <- rowSums(coef_matrix != 0)
  freq <- freq/n_boot
  mean_coef <- rowMeans(coef_matrix)
  df <- data.frame(Gene = rownames(coef_matrix),
                   SelectionFrequency = freq,
                   MeanCoefficient = mean_coef)
  df <- df[order(-df$SelectionFrequency), ]
  df_top <- head(df, top_n)
  p <- ggplot(df_top, aes(x = SelectionFrequency, y = MeanCoefficient)) +
    geom_point(aes(size = abs(MeanCoefficient)), fill = dot_col, color = "black", shape = 21) +
    geom_text(aes(label = Gene)) +
    labs(
      title = "Top Predictors by Bootstrap Selection Frequency",
      y = "Mean Coefficient",
      x = paste0("Selection Percentage (out of ", n_boot, " bootstraps %)")
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          axis.text = element_text(size = 15, color = "black"),
          axis.title = element_text(size = 15))
  if (return_data) {
    return(list(plot = p, data = df))
  } else {
    return(p)
  }
}








