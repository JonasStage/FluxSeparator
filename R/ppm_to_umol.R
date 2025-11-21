#' @title ppm_to_µmol
#'
#' @description
#'  Convertion of \eqn{ppmV} to \eqn{µmol} \eqn{m^{-2}} \eqn{h^{-1}} using the ideal gas law.
#'
#' @param pressure Air pressure during measurement (Pa).
#' @param concentration Concentration of the gas in ppm (µmol/mol).
#' @param volume Volume of the chamber used for measuring in \eqn{m^{3}}.
#' @param temperature_C Temperature in degrees celcius in the chamber.
#' @param area Surface area of the chamber used in \eqn{m^{2}}
#'
#' @author Jonas Stage Sø \email{Jonassoe@biology.sdu.dk}
#'
#' @examples
#' \dontrun{
#' library(FluxSeparator)
#'
#'# See ?ebullitive_flux or ?diffusive_flux for help with the respective functions
#'
#'# Using the ebullitive flux function
#' DIY_sensor_data %>%
#'   ebullitive_flux() %>%
#'   mutate(umol_m2_h1_flux = ppm_to_umol(pressure
#'                                        concentration
#'                                        volume
#'                                        temperature_C
#'                                        area))
#'
#'# Using the diffusive flux function
#' DIY_sensor_data %>%
#'  diffusive_flux(cutoff_start_value = 5) %>%
#'# 450 would be good for CO2, while 5 could be good for CH4
#'  mutate(umol_m2_h1_flux = ppm_to_umol(pressure
#'                                       concentration
#'                                       volume
#'                                       temperature_C
#'                                       area))
#' }
#'@export
ppm_to_umol <- function(pressure, concentration, volume, temperature_C, area) {
  ((pressure*concentration*volume)/(8.314*(temperature_C+ 273.15)))/area
}
