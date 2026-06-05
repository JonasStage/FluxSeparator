#' FluxSeparator: Separation of Diffusive and Ebullitive Fluxes
#'
#' Separates diffusive and ebullitive (bubble) fluxes from continuous
#' concentration measurements using a running variance approach.
#'
#' @author Jonas Stage Sø \email{Jonassoe@biology.sdu.dk}
#'
#' @name FluxSeparator-package
#' @aliases FluxSeparator
#'
#' @importFrom dplyr filter mutate group_by select rename summarise summarize reframe
#'   arrange pull ungroup left_join full_join inner_join join_by bind_rows
#'   tibble any_of contains between first last add_count add_tally lag n
#'   tally row_number if_else across
#' @importFrom lubridate ymd_hms as_datetime
#' @importFrom TTR runVar runMean
#' @importFrom HMR HMR
#' @importFrom broom tidy
#' @importFrom readr read_csv write_csv parse_number col_double col_character
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth geom_vline geom_hline
#'   scale_color_gradient scale_color_manual scale_x_datetime scale_y_continuous
#'   labs guides guide_colourbar theme element_text margin element_rect element_blank
#' @importFrom graphics par
#' @importFrom stringr str_detect
#' @importFrom purrr map map_dbl map_df
#' @importFrom stats lm
#' @importFrom tidyr drop_na unnest nest separate
#' @importFrom ggpubr ggarrange
#' @importFrom rlang .data
NULL

# Suppress R CMD check NOTEs for non-standard evaluation variables
utils::globalVariables(c(
  ".", "concentration_smooth", "model_coef",
  "concentration_raw_input", "obs_in_PumpCycle",
  "row", "concentration", "time", "smooth",
  "bubble_sum", "n_bubbles_count", "group_id",
  "start_value", "run_var5", "station",
  # FluxSeparatorApp column/variable references
  "CH4smV", "K33_CO2", "PumpCycle", "RsR0", "SampleNumber", "V0",
  "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9", "X10", "X11", "X12", "X13",
  "a", "abs_H", "airt", "b", "ch4_smv", "concentration_per_hour", "concentration_per_time",
  "datetime", "datetime_end", "datetime_start", "diff_flux", "ebul_flux",
  "first_concentration", "g", "k", "ppm_H20", "plot_number", "rh", "s",
  "sec", "second_concentration", "sensor", "slope_concentration_hr",
  "temp", "tempC", "water"
))
