# Internal helper: expand indices around detected events
#
# @param x Center index.
# @param index_span Number of observations before and after to include.
# @return Integer vector of indices (positive only).
# @noRd
get_ids_before_after <- function(x, index_span) {
  v <- (x - index_span):(x + index_span)
  v[v > 0]
}

# Internal helper: detect ebullitive events via running variance
#
# Computes a running variance on concentration data. Observations where the
# running variance exceeds `runvar_cutoff` (plus a surrounding window of
# `index_span`) are flagged as ebullitive.
#
# @param data Data frame with `concentration`, `PumpCycle`, `datetime` columns.
# @param runvar_cutoff Variance threshold for bubble detection.
# @param index_span Number of observations to include before/after each event.
# @param min_obs Minimum observations per PumpCycle to process.
# @param runvar_window Window size for the running variance.
# @param time_gap_seconds Time gap (seconds) to separate event groups.
# @return A list with:
#   - `running_var`: the full dataset with running variance and row numbers
#   - `bubble_ids`: integer vector of row indices flagged as ebullitive
# @noRd
detect_bubbles <- function(data, runvar_cutoff = 0.5, index_span = 30,
                           min_obs = 100, runvar_window = 5,
                           time_gap_seconds = 30) {
  running_var <- data %>%
    add_count(.data$PumpCycle) %>%
    filter(.data$n > min_obs) %>%
    drop_na(.data$concentration) %>%
    group_by(.data$PumpCycle) %>%
    mutate(run_var5 = runVar(.data$concentration, n = runvar_window)) %>%
    ungroup() %>%
    mutate(row = row_number())

  bubble_ids <- running_var %>%
    filter(.data$run_var5 > runvar_cutoff) %>%
    mutate(
      time_diff = .data$datetime - lag(.data$datetime),
      time_diff = as.numeric(.data$time_diff)
    ) %>%
    drop_na(.data$time_diff) %>%
    mutate(group_id = 1 + cumsum(.data$time_diff > time_gap_seconds)) %>%
    group_by(.data$group_id) %>%
    pull(.data$row) %>%
    map(~get_ids_before_after(., index_span)) %>%
    unlist() %>%
    unique()

  list(running_var = running_var, bubble_ids = bubble_ids)
}
