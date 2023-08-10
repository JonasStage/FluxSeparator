\name{ebullitive_flux}
\alias{ebullitive_flux}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ebullitive_flux
}
\description{
Separates the diffusive and ebullitive fluxes, to calculate the ebullitive flux, as the change in concentration from ebullitive events. This is done by computing a running variance, if the running variance exceed a customizable cutoff value it is considered an ebullitive event. Additional factors can be set to determine what is considered ebullitive events.

}
\usage{
ebullitive_flux(data,
concentration_values = "pred_CH4",
station, top_selection = "max",
runvar_cutoff = .5,
show_plots = TRUE, IndexSpan = 30,
CH4_diffusion_cutoff = 1,
number_of_pumpcycles_in_plot = 24,
smooth_data = FALSE)}

\arguments{
  \item{data}{Your data frame.}
  \item{concentration_values}{Name of your variable representing the concentration.}
  \item{station}{Variable to distinguish between stations.}
  \item{top_selection}{Can be set to "max" or "last" to either use the maximum concentration value or the last in each ebullitive event.}
  \item{runvar_cutoff}{Cutoff of the running variance, which is used to determine if an increase in concentration is an ebullitive event. Lower values increases number of ebullitive events registered.}
  \item{show_plots}{Show plots which can assist in the determination of good fits for the model. A boolean variable which should be TRUE or FALSE.}
  \item{IndexSpan}{Number of observations which are included before and after an ebullitive event, to ensure the entire event is determined.}
  \item{CH4_diffusion_cutoff}{A variable used to set a minimum concentration change that is considered an ebullitive event.}
  \item{number_of_pumpcycles_in_plot}{Number of cycles which are plotted. Used only if show_plots = TRUE.}
  \item{smooth_data}{Computes a running mean on the concentration data five times, to smoothen data if data is low bit resolution. See Sø et al., (2023) for more information}

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
  \item{sum_bubbles_concentration}{The sum of the differences in concentration caused by bubbbles}
  \item{n_bubbles}{Number of bubbles detected. Bear in mind that this function has dificulties detecting the number of bubbles if they are close to each other}
  \item{pumpcycle_duration_hr}{Length of the cycle duration in hours}
  \item{temp}{Average temperature within the chamber}
  \item{bubbles_per_time}{Amount of bubbles divided by the duration of the cycle in hours}
  \item{concentration_per_time}{Ebullitive flux, as the sum of concentration change divided by the duration in hours}
}
\references{
Methane and carbon dioxide fluxes at high spatiotemporal resolution from a small temperate lake. DOI: http://dx.doi.org/10.1016/j.scitotenv.2023.162895.
}
\author{
Jonas Stage Sø
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(FluxSeparator)

# Load dummy dataset
DIY_sensor_data <-  data(DIY_sensor_data)

DIY_sensor_data %>%
  ebullitive_flux()
}

