read_CH4_files <- function(data, files, pump_present = T, join_model_coef = T,
                           model_coef_data = model_coef) {

  lookup <- c(RH = "RH%")

  if(join_model_coef) {
    data <- data %>%
      inner_join(model_coef, by ="sensor") %>%
      rename(files = path)
  }
  else {
    data <- data %>%
      rename(files = path)
  }
  data %>%
    mutate(data = lapply(files,read_csv, show_col_types = T, col_types = list(
      col_double(),col_double(),col_character(),col_double(),
      col_double(),col_double(),col_double(),col_double(),
      col_double(),col_double(),col_double(),col_double(),col_double()))) %>%
    unnest(data) %>%
    rename(any_of(lookup)) %>%
    mutate(datetime = ymd_hms(datetime),
           abs_H = (6.112*exp((17.67*tempC)/(tempC+243.5))*RH*18.02)/((273.15+tempC)*100*0.08314),
           V0 = abs_H*g+S,
           RsR0 = ((5000/CH4smV)-1)/((5000/V0)-1),
           pred_CH4 = a*(RsR0^b)+c*abs_H*(a*RsR0^b) + K)  %>%
    select(files, datetime, RH:tempC,K33_RH:ncol(.), pred_CH4, sensor, abs_H, contains("volumen"), contains("station")) -> done_data
  return(done_data) }
