# FluxSeparator - Separation of diffusive and ebullitive fluxes


This R-package is dedicated to separation and calculation of diffusive and ebullitive (bubble) fluxes. 
Determination of ebullitive events is characterized by sudden increase in concentration over a short period of time. Thus, to determine ebullitive events a running variance approach is used. If the running variance is above a user-set threshold value, the data is considered an ebullitive event. Ebullitive events can then be calculated as the difference from the lowest concentration value to either the highest or the last value. Furthermore, correct identification of ebullitive events can be visually inspected when using the function. The functions also allows for several parameters to be changed, for full list of parameters see the R documentation for each function (called by ?_function_). The function for diffusive flux also looks for ebullitive events, to avoid diffusive flux being calculated post an ebullitive event. This function can likewise be turned off.   

For citing this package follow the instructions on https://zenodo.org/doi/10.5281/zenodo.8297153

Sø. J.S. (2023). JonasStage/FluxSeparator: v1.0.0 (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.8297154

For further information on package and usage see [Sø et al 2024](https://doi.org/10.1029/2024JG008035)

```
remotes::install_github('JonasStage/FluxSeparator')
```

## Functions
As of now the package consists of four functions and one dataset. A quick walkthrough of the functions is listed here, for more information and examples on usage see the R documentation of each function:

### Read DIY sensor files
```
library(FluxSeparator)

# read in model coef
model_coef <- read_csv("model_coef.csv")

# path to DIY sensors files
path_to_files <- list.files(pattern = ".csv")

# create data frame for path, sensor, and station.
data_path <- tibble(path = path_to_files,
                    sensor = c(1,2,3,4),
                    station = c(1,2,4,3))


# join with model_coef and calculate CH4 in ppm.
data <- read_CH4_files(data_path,
                       path)
```
The _read_CH4_files_ function allows for easy import and transformation of data files measured using the the low-cost methane and carbon dioxide sensors propsed by Sø et al., (2023). The function imports data and uses the individual calibration values from the dataframe _model_coef_ to transform the sensor voltage, temperature and humidity into methane concentration (ppm).

### Load data
```
data(DIY_sensor_data)
```
A dataset is included which holds two runs (PumpCycle), with the first only experiencing diffusive methane flux, whereas the second experiences ebullitive flux. This is data from an older version of the sensor, which has a low analog-to-digital converter, and thus the concentration changes occurs in steps. If using a newer version of the sensor this is not the case.

### Ebullitive flux
```
# Load dataset
DIY_sensor_data <-  data(DIY_sensor_data)

ebul_flux <- DIY_sensor_data %>%
               ebullitive_flux()
```
The _ebullitive_flux_ function allows for identification of ebullitive events in multiple data seperated by _station_ and _PumpCycle_. The function determines ebullitive events based on a running variance approach, in which a user-set threshold values is set (_runvar_cutoff_). The resulting dataframe is the ebullitive flux based on the ebullitive events determined. The function show plots of all data and areas which are considered ebullitive events, if _show_plots_ is set to _TRUE_.

### Diffusive flux
```
# Load dataset
DIY_sensor_data <-  data(DIY_sensor_data)

diff_flux <- DIY_sensor_data %>%
               diffusive_flux()
```
The _diffusive_flux_ function allows for calculation of diffusive flux before any ebullitive events or disregarding ebullitive events (_look_for_bubbles = FALSE_). The diffusive flux is calculated by a linear model and the corresponding slope, R<sup>2</sup> and number of observations is the result. The _runvar_cutoff_ for each sensor should be set to the same value as the corresponding _ebullitive_flux_.

### Conversion from ppm h<sup>-1</sup> to µmol m<sup>-2</sup> h<sup>-1</sup>
```
# Using the ebullitive flux function
DIY_sensor_data %>%
  ebullitive_flux() %>%
  mutate(umol_m2_h1_flux = ppm_to_umol(pressure,
                                 concentration_per_time,
                                 volume,
                                 temp,
                                 area))

# Using the diffusive flux function
DIY_sensor_data %>%
  diffusive_flux() %>%
  mutate(umol_m2_h1_flux = ppm_to_umol(pressure,
                                 concentration_per_time,
                                 volume,
                                 temp,
                                 area))
```
The _ppm_to_umol_ function allows for conversion of ppmV h<sup>-1</sup> to µmol m<sup>-2</sup> h<sup>-1</sup>. The function is called within the _mutate_ function and atmospheric pressure, volume and area should be given to the function. 


## Citations

Sø, J. S., Sand-Jensen, K., Martinsen, K. T., Polauke, E., Kjær, J. E., Reitzel, K., & Kragh, T. (2023). Methane and carbon dioxide fluxes at high spatiotemporal resolution from a small temperate lake. Science of the Total Environment, 878, 162895. doi:https://doi.org/10.1016/j.scitotenv.2023.162895

