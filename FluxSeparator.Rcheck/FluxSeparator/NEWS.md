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
