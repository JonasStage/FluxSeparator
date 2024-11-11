\name{diffusive_flux}
\alias{diffusive_flux}
\title{
diffusive_flux
}
\description{
Separates the diffusive and ebullitive fluxes, to calculate the diffusive flux, as a linear function of concentrations over time. This is done by firstly finding all events that are considered ebullitive (for more info see ebullitive_flux). Several factors can be set to determine what is considered ebullitive events, remove observations before doing diffusive flux, number of observations used in the diffusive flux calculation, cutoffs if concentrations start to high and number of observations needed in the linear model.
}
\usage{
diffusive_flux(data, concentration_values = "pred_CH4",
station, runvar_cutoff = 0.5,
remove_observations_prior = 200,
number_of_observations_used = 400,
show_plots = TRUE, IndexSpan = 30,
cutoff_start_value,
number_of_observations_required = 50,
number_of_pumpcycles_in_plot = 50,
smooth_data = FALSE,
look_for_bubbles = TRUE,
Hutchinson_Mosier_correction = FALSE,
volume,
area)}
\arguments{
  \item{data}{Your data frame.}
  \item{concentration_values}{Name of your variable representing the concentration.}
  \item{station}{Variable to distinguish between stations.}
  \item{runvar_cutoff}{Cutoff of the running variance, which is used to determine if an increase in concentration is an ebullitive event. Lower values increases number of ebullitive events registered.}
  \item{remove_observations_prior}{Remove n number of observations before calculating the diffusive flux by a linear slope.}
  \item{number_of_observations_used}{Number of observations used to calculate the diffusive flux by a linear slope.}
  \item{show_plots}{Show plots which can assist in the determination of good fits for the model. A boolean variable which should be TRUE or FALSE.}
  \item{IndexSpan}{Number of observations which are included before and after an ebullitive event, to ensure the entire event is determined.}
  \item{cutoff_start_value}{Variable indicating what the maximum starting concentration can be.}
  \item{number_of_observations_required}{Number of observations required in each cycle for the function to compute a linear model on the data}
  \item{number_of_pumpcycles_in_plot}{Number of cycles which are plotted. Used only if show_plots = TRUE.}
  \item{smooth_data}{Computes a running mean on the concentration data five times, to smoothen data if data is low bit resolution. See Sø et al., (2023) for more information}
  \item{look_for_bubbles}{Can be used for the function to not consider ebullitive events. Can be useful when calculating diffusive CO2 flux}
  \item{Hutchinson_Mosier_correction}{Can be used to correct flux measurements based on the Hutchinson-Mosier correction (1981). However, fluxes are only calculated of three points. A boolean variable which should be TRUE or FALSE}
  \item{volume}{Volume of the chamber used for calculating fluxes (L). This is only needed if the calculation of the Hutchinson_Mosier_correction = TRUE}
  \item{area}{Surface area of the chamber used for calculating fluxes \eqn{m^{2}}. This is only needed if the calculation of the Hutchinson_Mosier_correction = TRUE}

}
\details{
This code is similar to the code used in the scientific peer-reviewed article Sø et al., (2023).
}
\value{
A data frame containing the following
  \item{station}{The station provided to the input 'station'}
  \item{PumpCycle}{Cycle number}
  \item{datetime_start}{Start time of the cycle}
  \item{datetime_end}{End time of the cycle}
  \item{slope_concentration_hr}{The diffusiv flux per hour (ppm \eqn{h^{-1}})}
  \item{n_obs_included_in_lm}{Number of observations used to calculate the diffusive flux}
  \item{r2}{Variance explained by the linear model}
  \item{temp}{Average temperature within the chamber (C)}
  \item{hmr_slope}{Flux calculated using the HMR package (ppm \eqn{h^{-1}}). This is only returned if the method selected by the HMR package is either \emph{Hutchinson-Mosier correction} or \emph{No flux}, in case of \emph{LR} used the function will return NA. Only given if the Hutchinson-Mosier correction is calculated}
  \item{hmr_se}{Standard error of the \emph{hmr_slope}. This is only returned if the method selected by the \emph{hmr_slope} is calculated.}
  \item{hmr_pvalue}{The p-value for the null hypothesis of zero flux. This is only returned if the method selected by the \emph{hmr_slope} is calculated.}
  \item{hmr_lower95}{The lower end-point of the 95\%-confidence interval for the flux. This is only returned if the method selected by the \emph{hmr_slope} is calculated.}
  \item{hmr_upper95}{The upper end-point of the 95\%-confidence interval for the flux. This is only returned if the method selected by the \emph{hmr_slope} is calculated.}
  \item{method}{The method used for by the HMR package to calculate flux. This is only returned if the method selected by the \emph{hmr_slope} is calculated.}
  }
}
\references{
Methane and carbon dioxide fluxes at high spatiotemporal resolution from a small temperate lake. DOI: http://dx.doi.org/10.1016/j.scitotenv.2023.162895.

Self-Made Equipment for Automatic Methane Diffusion and Ebullition Measurements From Aquatic Environments. DOI: https://doi.org/10.1029/2024JG008035

Hutchinson, G.L. and Mosier, A.R. (1981). Improved soil cover method for field measurement of nitrous oxide fluxes. Soil Science Society of America Journal, 45, pp. 311-316

Pullens, J.W.M., Abalos, D., Petersen, S.O. and Pedersen, A.R. (2023). Identifying criteria for greenhouse gas flux estimation with automatic and manual chambers: A case study for N2O. Euro-pean Journal of Soil Science, 74, e13340. https://doi.org/10.1111/ejss.13340
}
\author{
Jonas Stage Sø.
}
\note{
For information on the Hutchinson-Mosier correction see the HMR-package
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(FluxSeparator)

# Load dataset
DIY_sensor_data <-  data(DIY_sensor_data)

DIY_sensor_data %>%
  diffusive_flux()
}

