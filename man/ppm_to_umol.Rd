\name{ppm_to_umol}
\alias{ppm_to_umol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ppm_to_umol}
\description{
Convertion of \eqn{ppmV} to \eqn{µmol} \eqn{m^{-2}} \eqn{h^{-1}} using the ideal gas law.}
\usage{
ppm_to_umol(pressure, concentration, volume, temperature_C, area)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pressure}{Air pressure during measurement (Pa).}
  \item{concentration}{Concentration of the gas in ppm (µmol/mol).}
  \item{volume}{Volume of the chamber used for measuring in \eqn{m^{3}}.}
  \item{temperature_C}{Temperature in degrees celcius in the chamber.}
  \item{area}{Surface area of the chamber used in \eqn{m^{2}}}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Methane and carbon dioxide fluxes at high spatiotemporal resolution from a small temperate lake. DOI: http://dx.doi.org/10.1016/j.scitotenv.2023.162895.
}
\author{
Jonas Stage Sø.
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

# Conversion of ppm h-1 to µmol m-2 h-1
# See ?ebullitive_flux or ?diffusive_flux for help with the respective functions

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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
