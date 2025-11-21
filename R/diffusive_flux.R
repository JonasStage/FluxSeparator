#' @title Diffusive flux
#'
#' @description
#' Separates the diffusive and ebullitive fluxes, to calculate the diffusive flux, as a linear function of concentrations over time. This is done by firstly finding all events that are considered ebullitive (for more info see ebullitive_flux).<br>Several factors can be set to determine what is considered ebullitive events, remove observations before doing diffusive flux, number of observations used in the diffusive flux calculation, cutoffs if concentrations start to high and number of observations needed in the linear model.
#' <br><br> Ouput data is converted to concentration change per hour.
#'
#' @param data Your data frame.
#' @param concentration_values Name of your variable representing the concentration.
#' @param station Variable to distinguish between stations.
#' @param runvar_cutoff Cutoff of the running variance, which is used to determine if an increase in concentration is an ebullitive event. Lower values increases number of ebullitive events registered.
#' @param remove_observations_prior Remove n number of observations before calculating the diffusive flux by a linear slope.
#' @param number_of_observations_used Number of observations used to calculate the diffusive flux by a linear slope.
#' @param show_plots Show plots which can assist in the determination of good fits for the model. A boolean variable which should be TRUE or FALSE.
#' @param IndexSpan Number of observations which are included before and after an ebullitive event, to ensure the entire event is determined.
#' @param cutoff_start_value Variable indicating what the maximum starting concentration can be.
#' @param number_of_observations_required Number of observations required in each cycle for the function to compute a linear model on the data
#' @param number_of_pumpcycles_in_plot Number of cycles which are plotted. Used only if show_plots = TRUE.
#' @param smooth_data Computes a running mean on the concentration data five times, to smoothen data if data is low bit resolution. See Sø et al., (2023) for more information
#' @param look_for_bubbles Can be used for the function to not consider ebullitive events. Can be useful when calculating diffusive CO2 flux
#' @param Hutchinson_Mosier_correction Can be used to correct flux measurements based on the Hutchinson-Mosier correction (1981). However, fluxes are only calculated of three points. A boolean variable which should be TRUE or FALSE
#' @param volume Volume of the chamber used for calculating fluxes (L). This is only needed if the calculation of the Hutchinson_Mosier_correction = TRUE
#' @param area Surface area of the chamber used for calculating fluxes m<sup>2</sup>. This is only needed if the calculation of the Hutchinson_Mosier_correction = TRUE
#'
#' @returns  A data frame containing the following:
#' \itemize{
#'    \item station - The station provided to the input 'station'
#'    \item PumpCycle - Cycle number
#'    \item datetime_start - Start time of the cycle
#'    \item datetime_end - End time of the cycle
#'    \item slope_concentration_hr - The diffusiv flux per hour (ppm h^<sup>-1</sup>)
#'    \item slope_standard_error - The standard error of the flux
#'    \item n_obs_included_in_lm - Number of observations used to calculate the diffusive flux
#'    \item r2 - Variance explained by the linear model
#'    \item temp - Average temperature within the chamber (C)
#'    \item hmr_slope - Flux calculated using the HMR package (ppm h^<sup>-1</sup>). This is only returned if the method selected by the HMR package is either \emph{Hutchinson-Mosier correction} or \emph{No flux}, in case of \emph{LR} used the function will return NA. Only given if the Hutchinson-Mosier correction is calculated
#'    \item hmr_se - Standard error of the \emph{hmr_slope}. This is only returned if the method selected by the \emph{hmr_slope} is calculated.
#'    \item hmr_pvalue - The p-value for the null hypothesis of zero flux. This is only returned if the method selected by the \emph{hmr_slope} is calculated.
#'    \item hmr_lower95 - The lower end-point of the 95%-confidence interval for the flux. This is only returned if the method selected by the \emph{hmr_slope} is calculated.
#'    \item hmr_upper95 - The upper end-point of the 95%-confidence interval for the flux. This is only returned if the method selected by the \emph{hmr_slope} is calculated.
#'    \item method - The method used for by the HMR package to calculate flux. This is only returned if the method selected by the \emph{hmr_slope} is calculated.
#'    }
#'
#' @references Sø et al. (2024). Self-Made Equipment for Automatic Methane Diffusion and Ebullition Measurements From Aquatic Environments. DOI: <a href="https://doi.org/10.1029/2024JG008035">https://doi.org/10.1029/2024JG008035</a>.
#' @references Sø et al. (2023). Methane and carbon dioxide fluxes at high spatiotemporal resolution from a small temperate lake. DOI: <a href="http://dx.doi.org/10.1016/j.scitotenv.2023.162895">http://dx.doi.org/10.1016/j.scitotenv.2023.162895</a>.
#' @references Hutchinson, G.L. and Mosier, A.R. (1981). Improved soil cover method for field measurement of nitrous oxide fluxes. Soil Science Society of America Journal, 45, pp. 311-316
#' @references Pullens, J.W.M., Abalos, D., Petersen, S.O. and Pedersen, A.R. (2023). Identifying criteria for greenhouse gas flux estimation with automatic and manual chambers: A case study for N2O. Euro-pean Journal of Soil Science, 74, e13340. <a href="https://doi.org/10.1111/ejss.13340">https://doi.org/10.1111/ejss.13340</a>
#'
#' @author Jonas Stage Sø \email{Jonassoe@biology.sdu.dk}
#'
#' @importFrom rlang .data
#'
#' @examples
#' library(FluxSeparator)
#'
#' data(DIY_sensor_data)
#'
#' DIY_sensor_data %>%
#'   diffusive_flux(cutoff_start_value = 450)
#'   # 450 would be good for CO2, while 5 could be good for CH4
#'
#' @export


diffusive_flux <- function(data, concentration_values = "pred_CH4", station, runvar_cutoff = 0.5, remove_observations_prior = 200,
                           number_of_observations_used = 400, show_plots = TRUE,IndexSpan = 30, cutoff_start_value,
                           number_of_observations_required = 50,number_of_pumpcycles_in_plot = 50, smooth_data = FALSE,
                           look_for_bubbles = TRUE, Hutchinson_Mosier_correction = FALSE, volume, area) {

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

  if(look_for_bubbles) {
    data %>%
      add_count(.data$PumpCycle) %>%
      filter(.data$n > 100) %>%
      drop_na(.data$concentration) %>%
      mutate(run_var5 = runVar(.data$concentration, n = 5)) %>%
      ungroup() %>%
      mutate(row = row_number()) -> running_var_diff

    running_var_diff %>%
      filter(.data$run_var5 > runvar_cutoff) %>%
      mutate(time_diff = .data$datetime -lag(.data$datetime),
             time_diff = as.numeric(.data$time_diff)) %>%
      drop_na(.data$time_diff) %>%
      mutate(gruppering =  1 + cumsum(.data$time_diff>30)) %>%
      group_by(.data$gruppering) %>%
      pull(.data$row) %>%
      map(~GetIDsBeforeAfter(., IndexSpan)) %>%
      unlist() %>%
      unique() -> ids_to_remain_diff

    running_var_diff %>%
      filter(!row %in% ids_to_remain_diff) %>%
      mutate(time_diff = .data$datetime -lag(.data$datetime),
             time_diff = as.numeric(.data$time_diff)) %>%
      drop_na(.data$time_diff) %>%
      mutate(gruppering =  1 + cumsum(.data$time_diff>30)) %>%
      arrange(.data$row) %>%
      group_by(.data$PumpCycle,.data$station) %>%
      mutate(first = first(.data$concentration),
             first = if_else(is.na(.data$first),.data$concentration,.data$first)) %>%
      drop_na(.data$time_diff) %>%
      mutate(min_grp = min(.data$gruppering)) %>%
      filter(.data$gruppering == .data$min_grp) %>%
      drop_na(.data$concentration) -> bubbles_removed_dataset

    bubbles_removed_dataset -> diffusive_dataset

    } else {
    data %>%
        group_by(.data$PumpCycle,.data$station) -> diffusive_dataset
      }

    diffusive_dataset %>%
      filter(between(row_number(),remove_observations_prior,(remove_observations_prior+number_of_observations_used))) %>%
      mutate(time = .data$datetime-min(.data$datetime)) -> dif_check

    dif_check %>%
      nest() %>%
      mutate(model = map(.data$data, ~lm(concentration ~ time, data = .)),
             model_summary = map(.data$model, summary),
             model = map(.data$model, tidy)) %>%
      unnest(.data$model) %>%
      filter(.data$term == "time") %>%
      mutate(start_value = map_dbl(data, ~min(.$concentration, na.rm=T)),
             n     = map_df(data, tally)$n,
             r2    = map_dbl(.data$model_summary, "r.squared"),
             temp  = map_dbl(data, ~mean(.$tempC)),
             datetime_start = map_dbl(data, ~min(.$datetime, na.rm=T)),
             datetime_start = as_datetime(.data$datetime_start),
             datetime_end = map_dbl(data, ~max(.$datetime, na.rm=T)),
             datetime_end = as_datetime(.data$datetime_end)) %>%
      rename(n_obs_included_in_lm = .data$n) %>%
      filter(.data$start_value < cutoff_start_value,
             .data$n_obs_included_in_lm > number_of_observations_required) %>%
      mutate(slope_concentration_hr = .data$estimate*3600,
             slope_standard_error = .data$std.error*3600,
             station = .data$station) %>%
      select(.data$PumpCycle, .data$station, .data$slope_concentration_hr,.data$slope_standard_error,p_value = .data$p.value, .data$temp, .data$n_obs_included_in_lm,.data$r2, .data$datetime_start,.data$datetime_end) -> diffusive_results

  if(Hutchinson_Mosier_correction) {
    diffusive_dataset %>%
      mutate(time = as.numeric(.data$datetime-min(.data$datetime))) %>%
      filter(between(row_number(),remove_observations_prior,(remove_observations_prior+number_of_observations_used))) %>%
      ungroup() %>%
      cbind(.data$volume,.data$area) %>%
      mutate(concentration = .data$concentration+abs(min(.data$concentration))+1,
             id = paste0(.data$station,"___",.data$PumpCycle)) %>%
      select(.data$id, .data$volume, .data$area, .data$time,.data$concentration) %>%
      arrange(.data$time) %>%
      write_csv(., file = paste0(tempdir(),"/FluxSeparator.csv"))

    wd <- getwd()
    setwd(tempdir())
    HMR("FluxSeparator.csv", sep = ",",IfNoValidHMR = "No flux",FollowHMR = T) -> hmr
    setwd(wd)

    hmr %>%
      tibble() %>%
      separate(.data$Series, c("station","PumpCycle", sep = "___")) %>%
      select(-"___") %>%
      cbind(.data$volume,.data$area) %>%
      mutate(across(.data$f0:.data$f0.up95, ~(parse_number(.x)*.data$area/.data$volume)*3600),
             PumpCycle = parse_number(.data$PumpCycle)) %>%
      select(.data$station:.data$Method) %>%
      rename(hmr_slope = .data$f0, hmr_se = .data$f0.se, hmr_pvalue = .data$f0.p, hmr_lower95 = .data$f0.lo95, hmr_upper95 = .data$f0.up95, method = .data$Method) %>%
      filter(!.data$method == "LR")-> hmr_results

    diffusive_results %>%
      left_join(hmr_results,
                      by = join_by(.data$station, .data$PumpCycle)) %>%
      select(.data$station,.data$datetime_start,.data$datetime_end,.data$slope_concentration_hr,.data$slope_standard_error,.data$n_obs_included_in_lm, .data$r2,.data$temp, .data$hmr_slope:.data$method)-> diff_flux
  } else {
        diffusive_results %>%
          select(.data$station,.data$datetime_start,.data$datetime_end,.data$slope_concentration_hr,.data$slope_standard_error,.data$n_obs_included_in_lm, .data$r2,.data$temp)-> diff_flux
      }

    plotting_data <- dif_check %>%
      left_join(diffusive_results, by = c("PumpCycle","station")) %>%
      mutate(plot_number = floor(.data$PumpCycle/number_of_pumpcycles_in_plot))
    data_indelt <- data %>%
      mutate(plot_number = floor(.data$PumpCycle/number_of_pumpcycles_in_plot))

  if(show_plots){
    for(i in unique(plotting_data$station)) {
      station_select = i
      for(j in unique(filter(plotting_data, station == station_select)$plot_number)){
        plot_number_select = j;
        plotting_data %>%
          filter(.data$station == station_select , .data$plot_number == plot_number_select) ->plot_lm
        data_indelt %>%
          filter(.data$station == station_select , .data$plot_number == plot_number_select) ->plot_raw
        par(ask=T)

        ggplot() +
          geom_point(data = plot_raw, aes(.data$datetime, .data$concentration, group = .data$PumpCycle)) +
          geom_smooth(data = plot_lm, aes(.data$datetime, .data$concentration, group = .data$PumpCycle, col = .data$r2),
                      method = "lm", se =F, linewidth = 2, formula = y ~ x,orientation = "y") +
          scale_color_gradient(limits = c(0,1), low = "red",high = "green") +
          labs(y = bquote("Concentration"),
               x = "Datetime",
               col = bquote("R"^2*"       "),
               title = paste0("This is station: ",i)) +
          guides(col = guide_colourbar(barwidth = 20)) +
          scale_x_datetime()+
          theme(legend.position = "bottom") ->p
        print(p)
      }}} else {}
    par(ask=F)

  return(diff_flux) }
