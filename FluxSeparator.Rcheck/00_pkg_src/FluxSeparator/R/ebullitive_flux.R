#' @title Ebullitive flux
#'
#' @description
#'  Separates the diffusive and ebullitive fluxes, to calculate the ebullitive
#'  flux, as the change in concentration from ebullitive events. This is done by
#'  computing a running variance; if the running variance exceeds a customizable
#'  cutoff value it is considered an ebullitive event.
#'
#'  Additional factors can be set to determine what is considered ebullitive
#'  events. Output data is converted to concentration change per hour.
#'
#' @param data Your data frame. Must contain columns \code{datetime},
#'   \code{PumpCycle}, \code{station}, and \code{tempC}.
#' @param concentration_values Name of your variable representing the
#'   concentration.
#' @param top_selection Can be set to \code{"last"} or \code{"max"} to either
#'   use the last or maximum concentration value in each ebullitive event.
#' @param runvar_cutoff Cutoff of the running variance, which is used to
#'   determine if an increase in concentration is an ebullitive event. Lower
#'   values increase the number of ebullitive events registered.
#' @param show_plots Show diagnostic plots. A logical; defaults to
#'   \code{TRUE}.
#' @param IndexSpan Number of observations which are included before and after
#'   an ebullitive event, to ensure the entire event is captured.
#' @param concentration_diffusion_cutoff Minimum concentration change that is
#'   considered an ebullitive event.
#' @param number_of_pumpcycles_in_plot Number of cycles which are plotted.
#'   Used only if \code{show_plots = TRUE}.
#' @param smooth_data Computes a running mean on the concentration data five
#'   times, to smooth data if data is low bit resolution. See
#'   Sø et al. (2023) for more information.
#' @param min_obs_per_cycle Minimum number of observations per pump cycle
#'   required for processing (default 100).
#' @param time_gap_seconds Time gap in seconds used to separate bubble event
#'   groups (default 30).
#' @param runvar_window Window size for the running variance computation
#'   (default 5).
#' @param smooth_window Window size for each running-mean smoothing pass
#'   (default 10). Only used when \code{smooth_data = TRUE}.
#'
#' @returns A data frame containing the following:
#' \itemize{
#'   \item{station \cr The station column from the input data}
#'   \item{PumpCycle \cr Cycle number}
#'   \item{datetime_start \cr Start time of the cycle}
#'   \item{datetime_end \cr End time of the cycle}
#'   \item{sum_bubbles_concentration \cr The sum of the differences in
#'     concentration caused by bubbles}
#'   \item{n_bubbles \cr Number of bubbles detected. Bear in mind that this
#'     function has difficulties detecting the number of bubbles if they are
#'     close to each other}
#'   \item{pumpcycle_duration_hr \cr Length of the cycle duration in hours}
#'   \item{temp \cr Average temperature within the chamber}
#'   \item{bubbles_per_time \cr Amount of bubbles divided by the duration of
#'     the cycle in hours}
#'   \item{concentration_per_time \cr Ebullitive flux, as the sum of
#'     concentration change divided by the duration in hours}
#' }
#' When \code{show_plots = TRUE}, a \code{"plots"} attribute is attached to
#' the result containing a list of \pkg{ggplot2} objects.
#'
#' @references Sø et al. (2024). Self-Made Equipment for Automatic Methane
#'   Diffusion and Ebullition Measurements From Aquatic Environments.
#'   \doi{10.1029/2024JG008035}.
#' @references Sø et al. (2023). Methane and carbon dioxide fluxes at high
#'   spatiotemporal resolution from a small temperate lake.
#'   \doi{10.1016/j.scitotenv.2023.162895}.
#'
#' @seealso \code{\link{diffusive_flux}}, \code{\link{ppm_to_umol}}
#'
#' @author Jonas Stage Sø \email{Jonassoe@@biology.sdu.dk}
#'
#' @examples
#' library(FluxSeparator)
#'
#' data(DIY_sensor_data)
#'
#' DIY_sensor_data %>%
#'   ebullitive_flux()
#'
#' @export
ebullitive_flux <- function(data,
                            concentration_values = "pred_CH4",
                            top_selection = "last",
                            runvar_cutoff = 0.5,
                            show_plots = TRUE,
                            IndexSpan = 30,
                            concentration_diffusion_cutoff = 1,
                            number_of_pumpcycles_in_plot = 24,
                            smooth_data = FALSE,
                            min_obs_per_cycle = 100,
                            time_gap_seconds = 30,
                            runvar_window = 5,
                            smooth_window = 10) {

  # --- input validation -------------------------------------------------------
  required_cols <- c("datetime", "PumpCycle", "station", "tempC")
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!concentration_values %in% names(data)) {
    stop("Column '", concentration_values, "' not found in data.", call. = FALSE)
  }
  match.arg(top_selection, c("last", "max"))

  # --- smoothing / renaming ---------------------------------------------------
  if (smooth_data) {
    data <- smooth_concentration(data, concentration_values,
                                 smooth_window = smooth_window,
                                 min_obs = min_obs_per_cycle)
  } else {
    data <- rename_concentration(data, concentration_values)
  }

  # --- cycle start / end times ------------------------------------------------
  data %>%
    group_by(.data$station, .data$PumpCycle) %>%
    reframe(datetime_start = min(.data$datetime, na.rm = TRUE),
            datetime_end   = max(.data$datetime, na.rm = TRUE)) -> times

  # --- bubble detection -------------------------------------------------------
  detection <- detect_bubbles(data,
                              runvar_cutoff    = runvar_cutoff,
                              index_span       = IndexSpan,
                              min_obs          = min_obs_per_cycle,
                              runvar_window    = runvar_window,
                              time_gap_seconds = time_gap_seconds)

  running_var   <- detection$running_var
  ids_to_remain <- detection$bubble_ids

  # --- baseline (no-bubble) summary per cycle ---------------------------------
  running_var %>%
    group_by(.data$station, .data$PumpCycle) %>%
    mutate(PumpCycle_Timediff = max(.data$datetime) - min(.data$datetime),
           PumpCycle_Timediff = as.numeric(.data$PumpCycle_Timediff,
                                           units = "hours")) %>%
    summarise(bubble_sum         = 0,
              n_bubbles_count    = 0,
              PumpCycle_Timediff = mean(.data$PumpCycle_Timediff),
              PumpCycle_Timediff_hr = as.numeric(.data$PumpCycle_Timediff),
              temp = mean(.data$tempC, na.rm = TRUE)) %>%
    rename(time = .data$PumpCycle_Timediff_hr) %>%
    select(-.data$PumpCycle_Timediff) -> no_bubbles

  # --- classify bubble events -------------------------------------------------
  running_var %>%
    filter(row %in% ids_to_remain) %>%
    mutate(time_diff = .data$datetime - lag(.data$datetime),
           time_diff = as.numeric(.data$time_diff)) %>%
    drop_na(.data$time_diff) %>%
    group_by(.data$station, .data$PumpCycle) %>%
    mutate(group_id = 1 + cumsum(.data$time_diff > time_gap_seconds)) %>%
    group_by(.data$station, .data$group_id, .data$PumpCycle) %>%
    mutate(first = first(.data$concentration),
           last  = last(.data$concentration),
           first = if_else(is.na(.data$first),
                           first(.data$concentration), .data$first),
           last  = if_else(is.na(.data$last),
                           last(.data$concentration), .data$last)) %>%
    filter(.data$first < .data$last) %>%
    bind_rows(filter(running_var, !row %in% ids_to_remain)) -> bubbles_check2

  # --- summarise each bubble event --------------------------------------------
  if (top_selection == "max") {
    bubbles_check2 %>%
      arrange(.data$row) %>%
      mutate(PumpCycle_Timediff = as.numeric(
        max(.data$datetime) - min(.data$datetime), units = "hours")) %>%
      summarize(
        time_diff          = max(.data$datetime) - min(.data$datetime),
        min_datetime       = .data$datetime[which.min(.data$concentration)],
        max_datetime       = .data$datetime[which.max(.data$concentration)],
        datetime           = mean(.data$datetime),
        min_concentration  = min(.data$concentration, na.rm = TRUE),
        top_concentration  = max(.data$concentration, na.rm = TRUE),
        concentration_diff = .data$top_concentration - .data$min_concentration,
        PumpCycle_Timediff = mean(.data$PumpCycle_Timediff),
        temp = mean(.data$tempC, na.rm = TRUE)
      ) %>%
      ungroup() -> bubbles_detected
  } else if (top_selection == "last") {
    bubbles_check2 %>%
      arrange(.data$row) %>%
      mutate(PumpCycle_Timediff = as.numeric(
        max(.data$datetime) - min(.data$datetime), units = "hours")) %>%
      summarize(
        time_diff          = max(.data$datetime) - min(.data$datetime),
        min_datetime       = .data$datetime[which.min(.data$concentration)],
        max_datetime       = .data$datetime[which.max(.data$concentration)],
        datetime           = mean(.data$datetime),
        min_concentration  = min(.data$concentration, na.rm = TRUE),
        top_concentration  = last(.data$concentration),
        concentration_diff = .data$top_concentration - .data$min_concentration,
        PumpCycle_Timediff = mean(.data$PumpCycle_Timediff),
        temp = mean(.data$tempC, na.rm = TRUE)
      ) %>%
      ungroup() -> bubbles_detected
  } else {
    stop("'top_selection' must be either 'max' or 'last'", call. = FALSE)
  }

  # --- filter + count bubbles per pump cycle ----------------------------------
  bubbles_detected %>%
    filter(.data$concentration_diff > concentration_diffusion_cutoff &
             .data$min_datetime < .data$max_datetime) %>%
    drop_na(.data$group_id) %>%
    add_count(.data$station, .data$PumpCycle) -> n_bubbles_per_pump

  # --- combine with no-bubble baseline and compute final summary --------------
  n_bubbles_per_pump %>%
    rename(bubble_sum      = .data$concentration_diff,
           time            = .data$PumpCycle_Timediff,
           n_bubbles_count = .data$n) %>%
    mutate(index = IndexSpan) %>%
    bind_rows(no_bubbles) %>%
    group_by(.data$station, .data$PumpCycle) %>%
    summarise(sum_bubbles_concentration = sum(.data$bubble_sum),
              n_bubbles            = max(.data$n_bubbles_count),
              pumpcycle_duration_hr = max(.data$time),
              temp                 = mean(.data$temp, na.rm = TRUE),
              bubbles_per_time     = .data$n_bubbles / .data$pumpcycle_duration_hr,
              concentration_per_time = .data$sum_bubbles_concentration /
                .data$pumpcycle_duration_hr) %>%
    full_join(times, by = c("station", "PumpCycle")) %>%
    select(.data$station, .data$PumpCycle, .data$datetime_start,
           .data$datetime_end,
           .data$sum_bubbles_concentration:.data$concentration_per_time) ->
    bubbles_found

  # --- diagnostic plots -------------------------------------------------------
  plotting_data <- running_var %>%
    mutate(plot_number = floor(.data$PumpCycle / number_of_pumpcycles_in_plot)) %>%
    full_join(n_bubbles_per_pump,
              by = c("station", "PumpCycle"),
              multiple = "all",
              suffix = c("", "_bubbles"))

  if (show_plots) {
    plot_list <- list()
    if (interactive()) {
      oldpar <- par(ask = TRUE)
      on.exit(par(oldpar), add = TRUE)
    }
    for (i in unique(plotting_data$station)) {
      wp_select <- i
      bubbles_found %>%
        filter(.data$station == wp_select) %>%
        ungroup() %>%
        summarize(n_bubles = sum(.data$n_bubbles)) -> bubles_count
      for (j in unique(filter(plotting_data,
                              .data$station == wp_select)$plot_number)) {
        plot_number_select <- j
        plotting_data %>%
          filter(.data$station == wp_select,
                 .data$plot_number == plot_number_select) -> plot1_dat
        plot1_dat %>%
          ggplot(aes(.data$datetime, .data$concentration,
                     group = .data$PumpCycle)) +
          geom_point() +
          geom_point(
            data = filter(drop_na(bubbles_check2, .data$group_id),
                          .data$station == wp_select),
            aes(.data$datetime, .data$concentration), col = "blue") +
          geom_vline(
            data = filter(plot1_dat,
                          .data$station == wp_select,
                          .data$concentration_diff > concentration_diffusion_cutoff),
            aes(xintercept = .data$datetime_bubbles), col = "red") +
          scale_x_datetime(limits = c(min(plot1_dat$datetime),
                                      max(plot1_dat$datetime))) +
          scale_y_continuous(limits = c(
            min(plot1_dat$concentration, na.rm = TRUE),
            max(plot1_dat$concentration, na.rm = TRUE))) +
          labs(y = bquote("Concentration"),
               x = "Datetime",
               title = paste0("This is station: ", i)) -> graf1
        ggplot() +
          geom_point(data = filter(plot1_dat, .data$run_var5 > 0.1),
                     aes(.data$datetime, .data$run_var5,
                         col = "run_var5 > 0.1")) +
          geom_point(data = filter(plot1_dat, .data$run_var5 > 0.2),
                     aes(.data$datetime, .data$run_var5,
                         col = "run_var5 > 0.2")) +
          geom_point(data = filter(plot1_dat, .data$run_var5 > 0.5),
                     aes(.data$datetime, .data$run_var5,
                         col = "run_var5 > 0.5")) +
          geom_point(data = filter(plot1_dat, .data$run_var5 > 1),
                     aes(.data$datetime, .data$run_var5,
                         col = "run_var5 > 1")) +
          scale_x_datetime(limits = c(min(plot1_dat$datetime),
                                      max(plot1_dat$datetime))) +
          scale_y_continuous(limits = c(
            0, max(plot1_dat$run_var5, na.rm = TRUE))) +
          scale_color_manual(
            limits = c("run_var5 > 0.1", "run_var5 > 0.2",
                        "run_var5 > 0.5", "run_var5 > 1"),
            labels = c("Variance > 0.1", "Variance > 0.2",
                        "Variance > 0.5", "Variance > 1"),
            values = c("red", "blue", "green", "black")) +
          labs(y = "Running variance", x = "Datetime", col = "") +
          geom_hline(yintercept = runvar_cutoff) +
          theme(legend.position = c(.9, .75)) -> graf2
        p <- ggarrange(graf1, graf2, ncol = 1)
        plot_list[[length(plot_list) + 1]] <- p
        if (interactive()) print(p)
      }
    }
    attr(bubbles_found, "plots") <- plot_list
  }

  return(bubbles_found)
}
