pkgname <- "FluxSeparator"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('FluxSeparator')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("DIY_sensor_data")
### * DIY_sensor_data

flush(stderr()); flush(stdout())

### Name: DIY_sensor_data
### Title: Data from Lake Lyng
### Aliases: DIY_sensor_data
### Keywords: datasets

### ** Examples

data(DIY_sensor_data)
head(DIY_sensor_data)



cleanEx()
nameEx("FluxSeparatorApp")
### * FluxSeparatorApp

flush(stderr()); flush(stdout())

### Name: FluxSeparatorApp
### Title: FluxSeparator Shiny App
### Aliases: FluxSeparatorApp

### ** Examples




cleanEx()
nameEx("diffusive_flux")
### * diffusive_flux

flush(stderr()); flush(stdout())

### Name: diffusive_flux
### Title: Diffusive flux
### Aliases: diffusive_flux

### ** Examples





cleanEx()
nameEx("ebullitive_flux")
### * ebullitive_flux

flush(stderr()); flush(stdout())

### Name: ebullitive_flux
### Title: Ebullitive flux
### Aliases: ebullitive_flux

### ** Examples

library(FluxSeparator)

data(DIY_sensor_data)

DIY_sensor_data %>%
  ebullitive_flux()




cleanEx()
nameEx("ppm_to_umol")
### * ppm_to_umol

flush(stderr()); flush(stdout())

### Name: ppm_to_umol
### Title: ppm_to_µmol
### Aliases: ppm_to_umol

### ** Examples

# Convert a single value
ppm_to_umol(pressure = 101325, concentration = 10,
            volume = 0.01, temperature_C = 20, area = 0.05)




cleanEx()
nameEx("read_CH4_files")
### * read_CH4_files

flush(stderr()); flush(stdout())

### Name: read_CH4_files
### Title: read_CH4_files
### Aliases: read_CH4_files

### ** Examples

## Not run: 
##D library(FluxSeparator)
##D 
##D # read in model coef
##D model_coef <- read_csv("model_coef.csv")
##D 
##D # path to DIY sensors files
##D path_to_files <- list.files(pattern = ".csv")
##D 
##D # create data frame for path, sensor and station.
##D data_path <- tibble(path = path_to_files,
##D                     sensor = c(1, 2, 3, 4),
##D                     station = c(1, 2, 4, 3))
##D 
##D 
##D # join with model_coef and calculate CH4 in ppm.
##D read_CH4_files(data_path, path,
##D                model_coef_data = model_coef)
##D 
##D #### Example using join_model_coef = FALSE ####
##D 
##D # join with model_coef.
##D joined_data_path <- left_join(data_path, model_coef, by = join_by(sensor))
##D 
##D # calculate CH4 in ppm.
##D read_CH4_files(joined_data_path,
##D                path,
##D                join_model_coef = FALSE)
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
