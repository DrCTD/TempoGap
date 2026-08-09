#' Plot Bootstrap LASSO Coefficient Selection Frequency
#'
#' Visualizes bootstrap LASSO regression results as a scatter plot of selection
#' frequency versus mean coefficient. Point size reflects the absolute mean
#' coefficient, point color reflects the coefficient direction (positive/negative),
#' and the top N genes by selection frequency are labeled.
#'
#' @param coef_matrix A matrix (or data.frame) with genes as rows (including an
#'   "Intercept" row, which will be automatically removed) and bootstrap
#'   iterations as columns. Typically extracted from
#'   \code{model_list[[cell_type]]$coef_matrix}.
#' @param bootstrap_times Integer. Number of bootstrap iterations used to
#'   calculate selection frequency. Default is 500.
#' @param top_n Integer. Number of top genes (by selection frequency) to label
#'   on the plot. Default is 20.
#' @param positive_color Character. Color for positive coefficients.
#'   Default is "#C44E52".
#' @param negative_color Character. Color for negative coefficients.
#'   Default is "#4C72B0".
#' @param point_size_range Numeric vector of length 2. Range of point sizes
#'   mapped to \code{AbsMeanCoefficient}. Default is c(1.5, 6).
#' @param base_size Numeric. Base font size for the theme. Default is 11.
#' @param base_family Character. Font family. Default is "Arial".
#' @param label_size Numeric. Font size for gene labels via
#'   \code{ggrepel::geom_text_repel}. Default is 3.5.
#' @param seed Integer or NULL. Random seed passed to \code{geom_text_repel}
#'   for reproducible label placement. Default is NULL (no seed set).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \dontrun{
#' p <- plot_bootstrap_coef(
#'   coef_matrix = gg2_predicted_model$`T-cell_cr_genes`$coef_matrix,
#'   bootstrap_times = 500,
#'   top_n = 20
#' )
#' p
#' }
#'
#' @export
plot_bootstrap_coef <- function(coef_matrix,
                                bootstrap_times = 500,
                                top_n = 20,
                                positive_color = "#C44E52",
                                negative_color = "#4C72B0",
                                point_size_range = c(1.5, 6),
                                base_size = 11,
                                base_family = "Arial",
                                label_size = 3.5,
                                seed = NULL) {
  
  if (!is.matrix(coef_matrix) && !is.data.frame(coef_matrix)) {
    stop("`coef_matrix` must be a matrix or data.frame.")
  }
  
  # 1. Compute selection frequency and mean coefficient --------------------
  selection_percentage <- rowSums(coef_matrix != 0) / bootstrap_times * 100
  
  coef_summary <- data.frame(
    Gene = rownames(coef_matrix),
    MeanCoefficient = rowSums(coef_matrix),
    SelectionPercentage = selection_percentage,
    stringsAsFactors = FALSE
  )
  
  plot_df <- coef_summary |>
    dplyr::mutate(
      AbsMeanCoefficient = abs(.data$MeanCoefficient),
      Direction = ifelse(.data$MeanCoefficient >= 0, "Positive", "Negative")
    )
  
  plot_df <- plot_df[plot_df$Gene != "Intercept", ]
  
  # 2. Select top N genes for labeling -------------------------------------
  label_df <- plot_df |>
    dplyr::arrange(dplyr::desc(.data$SelectionPercentage)) |>
    dplyr::slice_head(n = top_n)
  
  # 3. Color mapping --------------------------------------------------------
  direction_colors <- c(
    "Positive" = positive_color,
    "Negative" = negative_color
  )
  
  # 4. Build the plot --------------------------------------------------------
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$SelectionPercentage, y = .data$MeanCoefficient)
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linewidth = 0.4,
      color = "grey60",
      linetype = "dashed"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$AbsMeanCoefficient, fill = .data$Direction),
      shape = 21,
      color = "black",
      stroke = 0.35,
      alpha = 0.75
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = .data$Gene),
      size = label_size,
      seed = seed
    ) +
    ggplot2::scale_x_continuous(name = "Bootstrap selection frequency") +
    ggplot2::scale_y_continuous(name = "Mean LASSO coefficient") +
    ggplot2::scale_size_continuous(
      name = "|Mean coefficient|",
      range = point_size_range,
      breaks = scales::pretty_breaks(n = 4)
    ) +
    ggplot2::scale_fill_manual(values = direction_colors) +
    ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 12, color = "black"),
      axis.text = ggplot2::element_text(size = 10, color = "black"),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.4)
    )
  
  p
}
