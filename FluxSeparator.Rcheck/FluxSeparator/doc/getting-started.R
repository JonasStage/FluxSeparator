## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----load-data----------------------------------------------------------------
library(FluxSeparator)

data(DIY_sensor_data)
str(DIY_sensor_data)

## ----ebullitive---------------------------------------------------------------
ebul_flux <- ebullitive_flux(DIY_sensor_data, show_plots = FALSE)
ebul_flux

## ----diffusive----------------------------------------------------------------
diff_flux <- diffusive_flux(DIY_sensor_data,
                            cutoff_start_value = 5,
                            show_plots = FALSE)
diff_flux

## ----conversion---------------------------------------------------------------
ppm_to_umol(
  pressure = 101325,     # Pa (standard atmospheric)
  concentration = 10,    # ppm/hr
  volume = 0.01,         # m³
  temperature_C = 20,    # °C
  area = 0.05            # m²
)

## ----plots, eval = FALSE------------------------------------------------------
# result <- ebullitive_flux(DIY_sensor_data, show_plots = TRUE)
# 
# # Access individual plots
# plots <- attr(result, "plots")

## ----citation-----------------------------------------------------------------
citation("FluxSeparator")

