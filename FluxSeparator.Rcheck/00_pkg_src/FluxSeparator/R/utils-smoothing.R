# Internal helper: smooth concentration data using repeated running means
#
# @param data Data frame with concentration and grouping columns.
# @param concentration_values Name of the concentration column.
# @param smooth_window Window size for each running mean pass.
# @param min_obs Minimum observations per PumpCycle to process.
# @return Data frame with `concentration` and `smooth` columns added.
# @noRd
smooth_concentration <- function(data, concentration_values,
                                 smooth_window = 10, min_obs = 100) {
  data %>%
    drop_na(concentration_values) %>%
    group_by(.data$PumpCycle, .data$sensor) %>%
    rename(concentration_raw_input = any_of(concentration_values)) %>%
    add_tally(name = "obs_in_PumpCycle") %>%
    filter(.data$obs_in_PumpCycle > min_obs) %>%
    mutate(
      concentration_smooth = runMean(.data$concentration_raw_input, smooth_window),
      concentration_smooth = runMean(.data$concentration_smooth, smooth_window),
      concentration_smooth = runMean(.data$concentration_smooth, smooth_window),
      concentration_smooth = runMean(.data$concentration_smooth, smooth_window),
      concentration_smooth = runMean(.data$concentration_smooth, smooth_window),
      smooth = TRUE
    ) %>%
    rename(concentration = concentration_smooth)
}

# Internal helper: prepare concentration data without smoothing
#
# @param data Data frame.
# @param concentration_values Name of the concentration column.
# @return Data frame with `concentration` and `smooth` columns.
# @noRd
rename_concentration <- function(data, concentration_values) {
  data %>%
    rename(concentration = contains(concentration_values)) %>%
    mutate(smooth = FALSE)
}
