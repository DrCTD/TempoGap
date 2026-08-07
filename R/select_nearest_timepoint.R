#' Select Time Point Closest to Target for Each Patient
#'
#' For each patient, select the sample closest to the target time (e.g., 72 hours).
#'
#' @param df A data frame containing patient IDs, time values, and sample metadata.
#' @param patient_col Name of the column with patient IDs (string).
#' @param time_col Name of the column with time values (in hours).
#' @param target_time Numeric value indicating the target time to select (default: 72).
#' @param with_ties Logical. If TRUE, return all equally-close samples; if FALSE, return one.
#'
#' @return A data frame with the nearest sample(s) to the target time for each patient.
#' @export
select_nearest_timepoint <- function(df,
                                     patient_col = "patient_id",
                                     time_col = "time",
                                     target_time = 72,
                                     with_ties = FALSE) {
  # Load necessary library
  require(dplyr)
  
  df %>%
    group_by(.data[[patient_col]]) %>%
    mutate(time_diff = abs(.data[[time_col]] - target_time)) %>%
    slice_min(order_by = time_diff, n = 1, with_ties = with_ties) %>%
    ungroup()
}