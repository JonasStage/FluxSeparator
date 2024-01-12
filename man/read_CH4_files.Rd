\name{read_CH4_files}
\alias{read_CH4_files}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
read_CH4_files
}
\description{
A function to ease the import of data from DIY sensors, which reads a csv file, calculates the absolute humidity, V0, RsR0, and the concentration following Bastviken et al., (2020).
}
\usage{
read_CH4_files(data, files, pump_present = TRUE, join_model_coef = TRUE,
model_coef_data = model_coef)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{A data frame containing the the path, sensor identification, and model coefficients for this specific sensor. Model coefficients can also be read in as a seperate data frame and defined in the model_coef_data variable.}
  \item{files}{A vector supplying the path to the file being read}
  \item{pump_present}{A boolean value to express whether or not a pump is added present on the sensor. Standards to TRUE.}
}
\details{
This code is similar to the code used in the scientific peer-reviewed article Sø et al., (2023).
}
\value{
A data frame output including all the original values, with the exception of the model coefficients. pred_CH4 is computed based on the calibration model.
\item{pred_CH4}{Computed from the calibration model. The CH4 concentration calculated from the sensor resistance and expressed in ppm.}
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

# read in model coef
model_coef <- read_csv("model_coef.csv")

# path to DIY sensors files
path_to_files <- list.files(pattern = ".csv")

# create data frame for path, sensor and station.
data_path <- tibble(path = path_to_files,
                    sensor = c(1,2,3,4),
                    station = c(1,2,4,3))


# join with model_coef and calculate CH4 in ppm.
read_CH4_files(data_path,
               path)

#### Example using join_model_coef = FALSE ####

# join with model_coef.
joined_data_path <- left_join(data_path, model_coef, by = join_by(sensor))

# calculate CH4 in ppm.
read_CH4_files(joined_data_path,
               path,
               join_model_coef = FALSE)

}
