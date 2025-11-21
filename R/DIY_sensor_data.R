#' @title Data from Lake Lyng
#'
#' @description
#' A subset of data from the paper Methane and carbon dioxide fluxes at high spatiotemporal resolution from a small temperate lake.
#' One measurement containing only diffusive flux and one containing ebullitive events
#'
#' @name DIY_sensor_data
#'
#' @docType data
#'
#' @author Jonas Stage Sø \email{Jonassoe@biology.sdu.dk}
#'
#' @param datetime Date and time of measurement
#' @param RH\% Relative humidity
#' @param tempC Temperature in degree celcius
#' @param CH4smV Methane sensor voltage
#' @param K33_RH Relative humidity measured by the CO2 sensor
#' @param K33_Temp Temperature in degree celcius measured by the CO2 sensor
#' @param K33_CO2 CO2 concentration in ppm
#' @param SampleNumber Sample number in this pump cycle
#' @param PumpCycle Pump cycle which counts upwards after the chamber has been flushed
#' @param pred_CH4 Predicted methane concentration, this is automaticly calculated using the read_CH4_files function, and calibration values
#' @param station Station name. Needed for diffusive and ebullitive fluxes calculations
#' @param sensor Sensor name
#'
#' @return A data frame with 2,201 rows and 10 variables.
#'
#' @references Sø et al. (2023). Methane and carbon dioxide fluxes at high spatiotemporal resolution from a small temperate lake. DOI: <a href="http://dx.doi.org/10.1016/j.scitotenv.2023.162895">http://dx.doi.org/10.1016/j.scitotenv.2023.162895</a>.
#'
#' @source {Load a subset of the dataset.}
#' @examples data(DIY_sensor_data)
data(DIY_sensor_data, envir = environment())

