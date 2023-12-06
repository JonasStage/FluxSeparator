ebullitive_flux <- function(data,concentration_values = "pred_CH4",station, top_selection = "max",IndexSpan = 30,runvar_cutoff = .5,
                            show_plots = T, concentration_diffusion_cutoff = 1, number_of_pumpcycles_in_plot = 24,
                            smooth_data = F) {
  par(ask=T)
  GetIDsBeforeAfter = function(x,IndexSpan) {
    v = (x-IndexSpan) : (x+IndexSpan)
    v[v > 0]
  }

  if(smooth_data) {
    data %>%
      drop_na(concentration_values) %>%
      group_by(PumpCycle,sensor) %>%
      rename(concentration_raw = any_of(concentration_values)) %>%
      add_tally(name = "obs_in_PumpCycle") %>%
      filter(obs_in_PumpCycle > 100) %>%
      mutate(concentration_smooth = runMean(concentration_raw, 10),
             concentration_smooth = runMean(concentration_smooth, 10),
             concentration_smooth = runMean(concentration_smooth, 10),
             concentration_smooth = runMean(concentration_smooth, 10),
             concentration_smooth = runMean(concentration_smooth, 10)) -> data
  } else {
    data %>%
      rename(concentration_raw = contains(concentration_values)) -> data
  }

  data %>%
    colnames() %>%
    str_detect("concentration_smooth") %>%
    sum() -> smooth_present

  if(smooth_present == 1) {
    data %>% rename(concentration = concentration_smooth) -> data
  } else {
    data %>% rename(concentration = concentration_raw) -> data
  }

  data %>%
    group_by(station,PumpCycle) %>%
    reframe(datetime_start = min(datetime, na.rm=T),
            datetime_end = max(datetime, na.rm=T)) -> times

  data %>%
    add_count(PumpCycle) %>%
    filter(n > 100) %>%
    drop_na(concentration) %>%
    group_by(PumpCycle) %>%
    mutate(run_var5 = runVar(concentration, n = 5)) %>%
    ungroup() %>%
    mutate(row = row_number()) -> running_var

  running_var %>%
    filter(run_var5 > runvar_cutoff) %>%
    mutate(time_diff = datetime -lag(datetime),
           time_diff = as.numeric(time_diff)) %>%
    drop_na(time_diff) %>%
    mutate(gruppering =  1 + cumsum(time_diff>6)) %>%
    group_by(gruppering) %>%
    pull(row) %>%
    map(~GetIDsBeforeAfter(., IndexSpan)) %>%
    unlist() %>%
    unique() -> ids_to_remain

  running_var %>%
    rename(concentration = contains(concentration_values)) %>%
    group_by(station, PumpCycle) %>%
    mutate(PumpCycle_Timediff = max(datetime)-min(datetime),
           PumpCycle_Timediff =as.numeric(PumpCycle_Timediff, units = "hours")) %>%
    summarise(sum_bobler = 0,
              n_bobler = 0,
              PumpCycle_Timediff = mean(PumpCycle_Timediff),
              PumpCycle_Timediff_hr = as.numeric(PumpCycle_Timediff),
              temp = mean(tempC,na.rm=T)) %>%
    rename(time = PumpCycle_Timediff_hr) %>%
    select(-PumpCycle_Timediff) ->no_bobler

  running_var %>%
    filter(row %in% ids_to_remain) %>%
    mutate(time_diff = datetime -lag(datetime),
           time_diff = as.numeric(time_diff)) %>%
    drop_na(time_diff) %>%
    group_by(station, PumpCycle) %>%
    mutate(gruppering =  1 + cumsum(time_diff>6)) %>%
    group_by(station,gruppering,PumpCycle) %>%
    mutate(first = first(concentration),
           last = last(concentration),
           first = if_else(is.na(first), first(concentration), first),
           last = if_else(is.na(last), last(concentration), last)) %>%
    filter(first < last) %>%
    bind_rows(filter(running_var, !row %in% ids_to_remain)) -> bubbles_check2


  if(top_selection == "max") {
    bubbles_check2 %>%
      arrange(row) %>%
      mutate(PumpCycle_Timediff = as.numeric(max(datetime)-min(datetime), units = "hours")) %>%
      summarize(time_diff = max(datetime)-min(datetime),
                min_datetime = datetime[which.min(concentration)],
                max_datetime = datetime[which.max(concentration)],
                datetime = mean(datetime),
                min_concentration = min(concentration, na.rm=T),
                top_concentration = max(concentration, na.rm=T),
                concentration_diff = top_concentration-min_concentration,
                PumpCycle_Timediff = mean(PumpCycle_Timediff),
                temp = mean(tempC, na.rm=T)) %>%
      ungroup() -> bubbles_detected
  } else if (top_selection == "last") {
    bubbles_check2 %>%
      arrange(row) %>%
      mutate(PumpCycle_Timediff = as.numeric(max(datetime)-min(datetime), units = "hours")) %>%
      summarize(time_diff = max(datetime)-min(datetime),
                min_datetime = datetime[which.min(concentration)],
                max_datetime = datetime[which.max(concentration)],
                datetime = mean(datetime),
                min_concentration = min(concentration, na.rm=T),
                top_concentration = last(concentration),
                concentration_diff = top_concentration-min_concentration,
                PumpCycle_Timediff = mean(PumpCycle_Timediff),
                temp = mean(tempC, na.rm=T)) %>%
      ungroup() -> bubbles_detected
  } else {
    print("top_selection can only be max or last")
    }

  bubbles_detected %>%
    filter(concentration_diff > concentration_diffusion_cutoff & min_datetime < max_datetime) %>%
    drop_na(gruppering) %>%
    add_count(station,PumpCycle) -> n_bubbles_per_pump

  n_bubbles_per_pump %>%
    rename(sum_bobler = concentration_diff,
           time = PumpCycle_Timediff,
           n_bobler = n) %>%
    mutate(index = IndexSpan) %>%
    bind_rows(no_bobler) %>%
    group_by(station,PumpCycle) %>%
    summarise(sum_bubbles_concentration = sum(sum_bobler),
              n_bubbles = max(n_bobler),
              pumpcycle_duration_hr = max(time),
              temp = mean(temp, na.rm=T),
              bubbles_per_time = n_bubbles/pumpcycle_duration_hr,
              concentration_per_time = sum_bubbles_concentration/pumpcycle_duration_hr)  %>%
    full_join(times, by = c("station","PumpCycle")) %>%
    select(station,PumpCycle, datetime_start, datetime_end,
           sum_bubbles_concentration:concentration_per_time)-> bubbles_found

  plotting_data <- running_var %>%
    mutate(plot_number = floor(PumpCycle/number_of_pumpcycles_in_plot)) %>%
    full_join(n_bubbles_per_pump, by = c("station","PumpCycle"), multiple = "all",
              suffix = c("","_bubbles"))

  if(show_plots) {for(i in unique(plotting_data$station)) {
    wp_select = i
    bubbles_found %>%
      filter(station == wp_select) %>%
      ungroup %>%
      summarize(n_bubles = sum(n_bubbles)) -> bubles_count
    #cat("waypoint ",i," has a total of ",bubles_count$n_bubles, "bubbles\n")
    for(j in unique(filter(plotting_data, station == wp_select)$plot_number)){
      par(ask=T)
      plot_number_select = j;plotting_data %>%
        filter(station == wp_select , plot_number == plot_number_select) ->plot1_dat
      plot1_dat %>%
        ggplot(aes(datetime,concentration, group = PumpCycle)) +
        geom_point() +
        geom_point(data = filter(drop_na(bubbles_check2, gruppering), station == wp_select), aes(datetime, concentration), col = "blue") +
        geom_vline(data = filter(plot1_dat, station == wp_select, concentration_diff > concentration_diffusion_cutoff),
                   aes(xintercept= datetime_bubbles), col = "red")+
        scale_x_datetime(limits=c(min(plot1_dat$datetime), max(plot1_dat$datetime))) +
        scale_y_continuous(limits=c(min(plot1_dat$concentration,na.rm=T), max(plot1_dat$concentration,na.rm=T))) +
        labs(y = bquote("CH"[4]*" concentration (ppm)"),
             x = "Datetime",
             title = paste0("This is station: ",i)) -> graf1
      ggplot() +
        geom_point(data = filter(plot1_dat, run_var5 > 0.1), aes(datetime, run_var5, col = "run_var5 > 0.1")) +
        geom_point(data = filter(plot1_dat, run_var5 > 0.2), aes(datetime, run_var5, col = "run_var5 > 0.2")) +
        geom_point(data = filter(plot1_dat, run_var5 > 0.5), aes(datetime, run_var5, col = "run_var5 > 0.5")) +
        geom_point(data = filter(plot1_dat, run_var5 > 1), aes(datetime, run_var5, col = "run_var5 > 1")) +
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
