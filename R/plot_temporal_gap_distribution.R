#' Plot Temporal Gap Distribution Across Disease Process Times
#'
#' This function visulizes the distribution of Temporal Gap (predicted - actual time)
#' across patient actual times using loess function plots.
#'
#' @param gap_data A data frame containing at least three columns: \code{ActualTime}, \code{ScaledTempoGap}, and \code{PredictedTime}
#' @param xlab Label for the x-axis. Default is "Actual Disease Process Time".
#' @param ylab Label for the y-axis. Default is "Predicted Disease Process Time".
#' @param title Plot title. Default is NULL.
#'
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @export
#'
#' @examples
#' gap_data <- data.frame(
#'   ActualTime = rnorm(30)
#'   ScaledTempoGap = rnorm(30)
#'   PredictedTime = rnorm(30)
#' )
#' plot_temporal_gap_distribution <- function(gap_data)
plot_temporal_gap_distribution <- function(gap_data,
                                           xlab = "Actual Disease Process Time",
                                           ylab = "Predicted Disease Process Time",
                                           title = NULL){
  if(!requireNamespace("ggplot2", quietly = T)){
    stop("ggplot2 package is required for this function")
  }
  if(!("ScaledTempoGap" %in% colnames(gap_data))){
    stop("gap_data must contain a 'ScaledTempoGap' column.")
  }
  if(!("ActualTime" %in% colnames(gap_data))){
    stop("gap_data must contain a 'ActualTime' column.")
  }
  if(!("PredictedTime" %in% colnames(gap_data))){
    stop("gap_data must contain a 'PredictedTime' column.")
  }
  p <- ggplot(gap_data, aes(x = ActualTime, y = PredictedTime)) +
    geom_point(aes(color = ScaledTempoGap)) +
    geom_smooth(method = "loess", color = "black", se = F) +
    xlab("factual hours from injury") +
    ylab("predicteed hours from injury") +
    scale_color_gradient2(low = "darkblue", high = "darkred", mid = "white") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15))

  return(p)
}

















