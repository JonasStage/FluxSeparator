#' Data from Lake Lyng
#'
#' A subset of data from the paper Methane and carbon dioxide fluxes at high
#' spatiotemporal resolution from a small temperate lake. One measurement
#' containing only diffusive flux and one containing ebullitive events.
#'
#' @format A data frame with 2,201 rows and 12 variables:
#' \describe{
#'   \item{datetime}{Date and time of measurement}
#'   \item{RH}{Relative humidity (percent)}
#'   \item{tempC}{Temperature in degrees Celsius}
#'   \item{CH4smV}{Methane sensor voltage}
#'   \item{K33_RH}{Relative humidity measured by the CO2 sensor}
#'   \item{K33_Temp}{Temperature in degrees Celsius measured by the CO2 sensor}
#'   \item{K33_CO2}{CO2 concentration in ppm}
#'   \item{SampleNumber}{Sample number in this pump cycle}
#'   \item{PumpCycle}{Pump cycle which counts upwards after the chamber has been flushed}
#'   \item{pred_CH4}{Predicted methane concentration, calculated using the read_CH4_files function and calibration values}
#'   \item{station}{Station name, needed for diffusive and ebullitive flux calculations}
#'   \item{sensor}{Sensor identification number}
#' }
#'
#' @source Sø et al. (2023). Methane and carbon dioxide fluxes at high
#'   spatiotemporal resolution from a small temperate lake.
#'   \doi{10.1016/j.scitotenv.2023.162895}
#'
#' @examples
#' data(DIY_sensor_data)
#' head(DIY_sensor_data)
"DIY_sensor_data"

