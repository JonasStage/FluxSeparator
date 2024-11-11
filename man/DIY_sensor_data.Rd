\name{DIY_sensor_data}
\alias{DIY_sensor_data}
\docType{data}
\title{
Load DIY_sensor_data dataset
}
\description{
A subset of data from the paper Methane and carbon dioxide fluxes at high spatiotemporal resolution from a small temperate lake.
One measurement containing only diffusive flux and one containing ebullitive events
}
\usage{DIY_sensor_data}
\format{
  A data frame with 2400 observations on the following 10 variables.
  \describe{
    \item{\code{datetime}}{Date and time of measurement}
    \item{\samp{RH\%}}{Relative humidity}
    \item{\code{tempC}}{Temperature in degree celcius}
    \item{\code{CH4smV}}{Methane sensor voltage}
    \item{\code{K33_RH}}{Relative humidity measured by the CO2 sensor}
    \item{\code{K33_Temp}}{Temperature in degree celcius measured by the CO2 sensor}
    \item{\code{K33_CO2}}{CO2 concentration in ppm}
    \item{\code{SampleNumber}}{Sample number in this pump cycle}
    \item{\code{PumpCycle}}{Pump cycle which counts upwards after the chamber has been flushed}
    \item{\code{pred_CH4}}{Predicted methane concentration in ppm, this is automaticly calculated using the read_CH4_files function, and
calibration values}
    \item{\code{station}}{Station measured at. This is used to seperate stations when calling diffusive flux or ebullitive flux functions}
    \item{\code{sensor}}{Station measured with. This is used to seperate sensors when calling diffusive flux or ebullitive flux functions}
  }
}
\references{
Methane and carbon dioxide fluxes at high spatiotemporal resolution from a small temperate lake. DOI: http://dx.doi.org/10.1016/j.scitotenv.2023.162895.

Self-Made Equipment for Automatic Methane Diffusion and Ebullition Measurements From Aquatic Environments. DOI: https://doi.org/10.1029/2024JG008035.}
\examples{
data(DIY_sensor_data)
}
\keyword{datasets}
