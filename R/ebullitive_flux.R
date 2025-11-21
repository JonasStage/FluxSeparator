#' @title Ebullitive flux
#'
#' @description
#'  Separates the diffusive and ebullitive fluxes, to calculate the ebullitive flux, as the change in concentration from ebullitive events. This is done by computing a running variance, if the running variance exceed a customizable cutoff value it is considered an ebullitive event. <br>Additional factors can be set to determine what is considered ebullitive events.
#'  <br><br> Ouput data is converted to concentration change per hour.
#'
#' @param data Your data frame.
#' @param concentration_values Name of your variable representing the concentration.
#' @param station Variable to distinguish between stations.
#' @param top_selection Can be set to "last" or "max" to either use the last or maximum concentration value in each ebullitive event.
#' @param runvar_cutoff Cutoff of the running variance, which is used to determine if an increase in concentration is an ebullitive event. Lower values increases number of ebullitive events registered.
#' @param show_plots Show plots which can assist in the determination of good fits for the model. A boolean variable which should be TRUE or FALSE.
#' @param IndexSpan Number of observations which are included before and after an ebullitive event, to ensure the entire event is determined.
#' @param concentration_diffusion_cutoff A variable used to set a minimum concentration change that is considered an ebullitive event.
#' @param number_of_pumpcycles_in_plot Number of cycles which are plotted. Used only if show_plots = TRUE.
#' @param smooth_data Computes a running mean on the concentration data five times, to smoothen data if data is low bit resolution. See Sø et al., (2023) for more information
#'
#' @returns  A data frame containing the following:
#' \itemize{
#'    \item{station - The station provided to the input 'station'}
#'    \item{PumpCycle - Cycle number}
#'    \item{datetime_start - Start time of the cycle}
#'    \item{datetime_end - End time of the cycle}
#'    \item{sum_bubbles_concentration - The sum of the differences in concentration caused by bubbbles}
#'    \item{n_bubbles - Number of bubbles detected. Bear in mind that this function has dificulties detecting the number of bubbles if they are close to each other}
#'    \item{pumpcycle_duration_hr - Length of the cycle duration in hours}
#'    \item{temp - Average temperature within the chamber}
#'    \item{bubbles_per_time - Amount of bubbles divided by the duration of the cycle in hours}
#'    \item{concentration_per_hour - Ebullitive flux, as the sum of concentration change divided by the duration in hours}
#'    }
#'
#' @references Sø et al. (2024). Self-Made Equipment for Automatic Methane Diffusion and Ebullition Measurements From Aquatic Environments. DOI: <a href="https://doi.org/10.1029/2024JG008035">https://doi.org/10.1029/2024JG008035</a>.
#' @references Sø et al. (2023). Methane and carbon dioxide fluxes at high spatiotemporal resolution from a small temperate lake. DOI: <a href="http://dx.doi.org/10.1016/j.scitotenv.2023.162895">http://dx.doi.org/10.1016/j.scitotenv.2023.162895</a>.
#' @author Jonas Stage Sø \email{Jonassoe@biology.sdu.dk}
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


ebullitive_flux <- function(data,concentration_values = "pred_CH4",station, top_selection = "last",runvar_cutoff = .5,
                            show_plots = TRUE, IndexSpan = 30, concentration_diffusion_cutoff = 1, number_of_pumpcycles_in_plot = 24,
                            smooth_data = FALSE) {
  par(ask=T)
  GetIDsBeforeAfter = function(x,IndexSpan) {
    v = (x-IndexSpan) : (x+IndexSpan)
    v[v > 0]
  }

  if(smooth_data) {
    data %>%
      drop_na(concentration_values) %>%
      group_by(.data$PumpCycle,.data$sensor) %>%
      rename(concentration_raw_input = any_of(concentration_values)) %>%
      add_tally(name = "obs_in_PumpCycle") %>%
      filter(.data$obs_in_PumpCycle > 100) %>%
      mutate(concentration_smooth = runMean(.data$concentration_raw_input, 10),
             concentration_smooth = runMean(.data$concentration_smooth, 10),
             concentration_smooth = runMean(.data$concentration_smooth, 10),
             concentration_smooth = runMean(.data$concentration_smooth, 10),
             concentration_smooth = runMean(.data$concentration_smooth, 10),
             smooth = T) %>%
      rename(concentration = concentration_smooth) -> data
  } else {
    data %>%
      rename(concentration = contains(concentration_values)) %>%
      mutate(smooth = F) -> data
  }

  data %>%
    group_by(.data$station,.data$PumpCycle) %>%
    reframe(datetime_start = min(.data$datetime, na.rm=T),
            datetime_end = max(.data$datetime, na.rm=T)) -> times

  data %>%
    add_count(.data$PumpCycle) %>%
    filter(.data$n > 100) %>%
    drop_na(.data$concentration) %>%
    group_by(.data$PumpCycle) %>%
    mutate(run_var5 = runVar(.data$concentration, n = 5)) %>%
    ungroup() %>%
    mutate(row = row_number()) -> running_var

  running_var %>%
    filter(.data$run_var5 > runvar_cutoff) %>%
    mutate(time_diff = .data$datetime -lag(.data$datetime),
           time_diff = as.numeric(.data$time_diff)) %>%
    drop_na(.data$time_diff) %>%
    mutate(gruppering =  1 + cumsum(.data$time_diff>30)) %>%
    group_by(.data$gruppering) %>%
    pull(.data$row) %>%
    map(~GetIDsBeforeAfter(., IndexSpan)) %>%
    unlist() %>%
    unique() -> ids_to_remain

  running_var %>%
    rename(concentration = contains(concentration_values)) %>%
    group_by(.data$station, .data$PumpCycle) %>%
    mutate(PumpCycle_Timediff = max(.data$datetime)-min(.data$datetime),
           PumpCycle_Timediff =as.numeric(.data$PumpCycle_Timediff, units = "hours")) %>%
    summarise(sum_bobler = 0,
              n_bobler = 0,
              PumpCycle_Timediff = mean(.data$PumpCycle_Timediff),
              PumpCycle_Timediff_hr = as.numeric(.data$PumpCycle_Timediff),
              temp = mean(.data$tempC,na.rm=T)) %>%
    rename(time = .data$PumpCycle_Timediff_hr) %>%
    select(-.data$PumpCycle_Timediff) ->no_bobler

  running_var %>%
    filter(row %in% ids_to_remain) %>%
    mutate(time_diff = .data$datetime -lag(.data$datetime),
           time_diff = as.numeric(.data$time_diff)) %>%
    drop_na(.data$time_diff) %>%
    group_by(.data$station, .data$PumpCycle) %>%
    mutate(gruppering =  1 + cumsum(.data$time_diff>30)) %>%
    group_by(.data$station,.data$gruppering,.data$PumpCycle) %>%
    mutate(first = first(.data$concentration),
           last = last(.data$concentration),
           first = if_else(is.na(.data$first), first(.data$concentration), .data$first),
           last = if_else(is.na(.data$last), last(.data$concentration), .data$last)) %>%
    filter(.data$first < .data$last) %>%
    bind_rows(filter(running_var, !row %in% ids_to_remain)) -> bubbles_check2


  if(top_selection == "max") {
    bubbles_check2 %>%
      arrange(.data$row) %>%
      mutate(PumpCycle_Timediff = as.numeric(max(.data$datetime)-min(.data$datetime), units = "hours")) %>%
      summarize(time_diff = max(.data$datetime)-min(.data$datetime),
                min_datetime = .data$datetime[which.min(.data$concentration)],
                max_datetime = .data$datetime[which.max(.data$concentration)],
                datetime = mean(.data$datetime),
                min_concentration = min(.data$concentration, na.rm=T),
                top_concentration = max(.data$concentration, na.rm=T),
                concentration_diff = .data$top_concentration-.data$min_concentration,
                PumpCycle_Timediff = mean(.data$PumpCycle_Timediff),
                temp = mean(.data$tempC, na.rm=T)) %>%
      ungroup() -> bubbles_detected
  } else if (top_selection == "last") {
    bubbles_check2 %>%
      arrange(.data$row) %>%
      mutate(PumpCycle_Timediff = as.numeric(max(.data$datetime)-min(.data$datetime), units = "hours")) %>%
      summarize(time_diff = max(.data$datetime)-min(.data$datetime),
                min_datetime = .data$datetime[which.min(.data$concentration)],
                max_datetime = .data$datetime[which.max(.data$concentration)],
                datetime = mean(.data$datetime),
                min_concentration = min(.data$concentration, na.rm=T),
                top_concentration = last(.data$concentration),
                concentration_diff = .data$top_concentration-.data$min_concentration,
                PumpCycle_Timediff = mean(.data$PumpCycle_Timediff),
                temp = mean(.data$tempC, na.rm=T)) %>%
      ungroup() -> bubbles_detected
  } else {
    print("top_selection can only be max or last")
    }

  bubbles_detected %>%
    filter(.data$concentration_diff > concentration_diffusion_cutoff & .data$min_datetime < .data$max_datetime) %>%
    drop_na(.data$gruppering) %>%
    add_count(.data$station,.data$PumpCycle) -> n_bubbles_per_pump

  n_bubbles_per_pump %>%
    rename(sum_bobler = .data$concentration_diff,
           time = .data$PumpCycle_Timediff,
           n_bobler = .data$n) %>%
    mutate(index = IndexSpan) %>%
    bind_rows(no_bobler) %>%
    group_by(.data$station,.data$PumpCycle) %>%
    summarise(sum_bubbles_concentration = sum(.data$sum_bobler),
              n_bubbles = max(.data$n_bobler),
              pumpcycle_duration_hr = max(.data$time),
              temp = mean(.data$temp, na.rm=T),
              bubbles_per_time = .data$n_bubbles/.data$pumpcycle_duration_hr,
              concentration_per_time = .data$sum_bubbles_concentration/.data$pumpcycle_duration_hr)  %>%
    full_join(times, by = c("station","PumpCycle")) %>%
    select(.data$station,.data$PumpCycle, .data$datetime_start, .data$datetime_end,
           .data$sum_bubbles_concentration:.data$concentration_per_time)-> bubbles_found

  plotting_data <- running_var %>%
    mutate(plot_number = floor(.data$PumpCycle/number_of_pumpcycles_in_plot)) %>%
    full_join(n_bubbles_per_pump, by = c("station","PumpCycle"), multiple = "all",
              suffix = c("","_bubbles"))

  if(show_plots) {for(i in unique(plotting_data$station)) {
    wp_select = i
    bubbles_found %>%
      filter(.data$station == wp_select) %>%
      ungroup %>%
      summarize(n_bubles = sum(.data$n_bubbles)) -> bubles_count
    #cat("waypoint ",i," has a total of ",bubles_count$n_bubles, "bubbles\n")
    for(j in unique(filter(plotting_data, .data$station == wp_select)$plot_number)){
      par(ask=T)
      plot_number_select = j;plotting_data %>%
        filter(.data$station == wp_select , .data$plot_number == plot_number_select) ->plot1_dat
      plot1_dat %>%
        ggplot(aes(.data$datetime,.data$concentration, group = .data$PumpCycle)) +
        geom_point() +
        geom_point(data = filter(drop_na(bubbles_check2, .data$gruppering), .data$station == wp_select), aes(.data$datetime, .data$concentration), col = "blue") +
        geom_vline(data = filter(plot1_dat, .data$station == wp_select, .data$concentration_diff > concentration_diffusion_cutoff),
                   aes(xintercept= .data$datetime_bubbles), col = "red")+
        scale_x_datetime(limits=c(min(plot1_dat$datetime), max(plot1_dat$datetime))) +
        scale_y_continuous(limits=c(min(plot1_dat$concentration,na.rm=T), max(plot1_dat$concentration,na.rm=T))) +
        labs(y = bquote("Concentration"),
             x = "Datetime",
             title = paste0("This is station: ",i)) -> graf1
      ggplot() +
        geom_point(data = filter(plot1_dat, .data$run_var5 > 0.1), aes(.data$datetime, .data$run_var5, col = "run_var5 > 0.1")) +
        geom_point(data = filter(plot1_dat, .data$run_var5 > 0.2), aes(.data$datetime, .data$run_var5, col = "run_var5 > 0.2")) +
        geom_point(data = filter(plot1_dat, .data$run_var5 > 0.5), aes(.data$datetime, .data$run_var5, col = "run_var5 > 0.5")) +
        geom_point(data = filter(plot1_dat, .data$run_var5 > 1), aes(.data$datetime, .data$run_var5, col = "run_var5 > 1")) +
        scale_x_datetime(limits=c(min(plot1_dat$datetime), max(plot1_dat$datetime))) +
        scale_y_continuous(limits=c(0,max(plot1_dat$run_var5,na.rm=T))) +
        scale_color_manual(limits = c("run_var5 > 0.1","run_var5 > 0.2","run_var5 > 0.5","run_var5 > 1"),
                           labels = c("Variance > 0.1","Variance > 0.2","Variance > 0.5","Variance > 1"),
                           values = c("red","blue","green","black")) +
        labs(y = "Running variance", x = "Datetime", col = "") +
        geom_hline(yintercept = runvar_cutoff) + theme(legend.position=c(.9,.75))->graf2
      ggarrange(graf1,graf2, ncol = 1) ->p
      print(p)
    }}} else {}
  par(ask=F)
  return(bubbles_found)}
