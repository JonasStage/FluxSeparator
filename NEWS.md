# FluxSeparator 2.0.1
* Minor update to citations

# FluxSeparator 2.0.0
* Included a shiny app `FluxSeparatorApp()`
  * The app helps visualize the process of determining diffusive and ebullitive fluxes
  * Allows the user to mark areas of diffusive fluxes by dragging the mouse. 
  * Can be used with DIY sensors or common datatypes like .csv and .data. 

# FluxSeparator 1.0.1

* Moved all dependencies from `Depends` to `Imports` for CRAN compliance.
* Replaced broad `@import` directives with specific `@importFrom` statements.
* Fixed interactive plotting to properly save/restore graphical parameters.
* Removed `setwd()` usage in `diffusive_flux()`.
* Improved documentation and examples across all exported functions.
* Added unit tests with `testthat`.

# FluxSeparator 1.0.0

* Initial release.
* Core functions: `ebullitive_flux()`, `diffusive_flux()`, `ppm_to_umol()`,
  `read_CH4_files()`.
* Included `DIY_sensor_data` example dataset from Lake Lyng.
