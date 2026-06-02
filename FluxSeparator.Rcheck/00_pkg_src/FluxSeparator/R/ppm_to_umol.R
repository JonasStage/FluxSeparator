#' @title ppm_to_µmol
#'
#' @description
#'  Conversion of \eqn{ppmV} to \eqn{\mu mol} \eqn{m^{-2}} \eqn{h^{-1}} using the ideal gas law.
#'
#' @param pressure Air pressure during measurement (Pa).
#' @param concentration Concentration of the gas in ppm (µmol/mol).
#' @param volume Volume of the chamber used for measuring in \eqn{m^{3}}.
#' @param temperature_C Temperature in degrees Celsius in the chamber.
#' @param area Surface area of the chamber used in \eqn{m^{2}}.
#'
#' @returns A numeric vector of flux values in \eqn{\mu mol} \eqn{m^{-2}} \eqn{h^{-1}}.
#'
#' @author Jonas Stage Sø \email{Jonassoe@biology.sdu.dk}
#'
#' @examples
#' # Convert a single value
#' ppm_to_umol(pressure = 101325, concentration = 10,
#'             volume = 0.01, temperature_C = 20, area = 0.05)
#'
#' @export
ppm_to_umol <- function(pressure, concentration, volume, temperature_C, area) {
  ((pressure*concentration*volume)/(8.314*(temperature_C+ 273.15)))/area
}
