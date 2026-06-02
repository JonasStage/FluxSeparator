#' @title read_CH4_files
#'
#' @description
#' A function to ease the import of data from DIY sensors, which reads a CSV
#' file, calculates the absolute humidity, V0, RsR0, and the methane
#' concentration.
#'
#' @param data A data frame containing the path, sensor identification, and
#'   model coefficients for this specific sensor. Model coefficients can also
#'   be read in as a separate data frame and defined in the
#'   \code{model_coef_data} variable.
#' @param files A vector supplying the path to the file being read.
#' @param join_model_coef Boolean variable. Join data with dataframe
#'   \code{model_coef_data} to convert sensor voltage signal to methane
#'   concentration.
#' @param model_coef_data Data frame consisting of the calibration values used
#'   to convert sensor voltage signal to methane concentration. Required when
#'   \code{join_model_coef = TRUE}.
#'
#' @returns A data frame output including all the original values, with the
#'   exception of the model coefficients.
#' \describe{
#'   \item{pred_CH4}{Computed from the calibration model. The CH4 concentration
#'     calculated from the sensor resistance, expressed in ppm.}
#' }
#'
#' @author Jonas Stage Sø \email{Jonassoe@biology.sdu.dk}
#'
#' @references Sø et al. (2023). Methane and carbon dioxide fluxes at high
#'   spatiotemporal resolution from a small temperate lake.
#'   \doi{10.1016/j.scitotenv.2023.162895}
#' @references Sø et al. (2024). Self-Made Equipment for Automatic Methane
#'   Diffusion and Ebullition Measurements From Aquatic Environments.
#'   \doi{10.1029/2024JG008035}
#'
#' @seealso \code{\link{ebullitive_flux}}, \code{\link{diffusive_flux}},
#'   \code{\link{ppm_to_umol}}
#'
#' @examples
#' \dontrun{
#' library(FluxSeparator)
#'
#' # read in model coef
#' model_coef <- read_csv("model_coef.csv")
#'
#' # path to DIY sensors files
#' path_to_files <- list.files(pattern = ".csv")
#'
#' # create data frame for path, sensor and station.
#' data_path <- tibble(path = path_to_files,
#'                     sensor = c(1, 2, 3, 4),
#'                     station = c(1, 2, 4, 3))
#'
#'
#' # join with model_coef and calculate CH4 in ppm.
#' read_CH4_files(data_path, path,
#'                model_coef_data = model_coef)
#'
#' #### Example using join_model_coef = FALSE ####
#'
#' # join with model_coef.
#' joined_data_path <- left_join(data_path, model_coef, by = join_by(sensor))
#'
#' # calculate CH4 in ppm.
#' read_CH4_files(joined_data_path,
#'                path,
#'                join_model_coef = FALSE)
#' }
#' @export
read_CH4_files <- function(data, files, join_model_coef = TRUE,
                           model_coef_data = NULL) {

  if (join_model_coef && is.null(model_coef_data)) {
    stop("'model_coef_data' must be provided when join_model_coef = TRUE",
         call. = FALSE)
  }

  lookup <- c(RH = "RH%")

  if (join_model_coef) {
    data <- data %>%
      inner_join(model_coef_data, by = "sensor") %>%
      rename(files = .data$path)
  } else {
    data <- data %>%
      rename(files = .data$path)
  }
  data %>%
    mutate(data = lapply(.data$files, read_csv, show_col_types = TRUE, col_types = list(
      col_double(), col_double(), col_character(), col_double(),
      col_double(), col_double(), col_double(), col_double(),
      col_double(), col_double(), col_double(), col_double(), col_double()))) %>%
    unnest(.data$data) %>%
    rename(any_of(lookup)) %>%
    mutate(datetime = ymd_hms(.data$datetime),
           abs_H = (6.112 * exp((17.67 * .data$tempC) / (.data$tempC + 243.5)) * .data$RH * 18.02) / ((273.15 + .data$tempC) * 100 * 0.08314),
           V0 = .data$abs_H * .data$g + .data$S,
           RsR0 = ((5000 / .data$CH4smV) - 1) / ((5000 / .data$V0) - 1),
           pred_CH4 = .data$a * (.data$RsR0^.data$b) + .data$c * .data$abs_H * (.data$a * .data$RsR0^.data$b) + .data$K) %>%
    select(.data$files, .data$datetime, .data$RH:.data$tempC, .data$K33_RH:.data$SampleNumber, contains("PumpCycle"), .data$pred_CH4, .data$sensor, .data$abs_H, contains("volumen"), contains("station")) -> done_data
  return(done_data)
}
