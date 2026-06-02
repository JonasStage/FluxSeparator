test_that("diffusive_flux runs on example data", {
  data("DIY_sensor_data", package = "FluxSeparator")
  result <- diffusive_flux(DIY_sensor_data,
                           cutoff_start_value = 5,
                           show_plots = FALSE)
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expected_cols <- c("station", "datetime_start", "datetime_end",
                     "slope_concentration_hr", "slope_standard_error",
                     "n_obs_included_in_lm", "r2", "temp")
  for (col in expected_cols) {
    expect_true(col %in% names(result),
                info = paste("Missing column:", col))
  }
})

test_that("diffusive_flux works with look_for_bubbles = FALSE", {
  data("DIY_sensor_data", package = "FluxSeparator")
  result <- diffusive_flux(DIY_sensor_data,
                           cutoff_start_value = 5,
                           look_for_bubbles = FALSE,
                           show_plots = FALSE)
  expect_s3_class(result, "data.frame")
})

test_that("diffusive_flux validates required columns", {
  bad_data <- data.frame(x = 1:10, y = 1:10)
  expect_error(
    diffusive_flux(bad_data, show_plots = FALSE),
    "Missing required columns"
  )
})

test_that("diffusive_flux validates concentration column exists", {
  data("DIY_sensor_data", package = "FluxSeparator")
  expect_error(
    diffusive_flux(DIY_sensor_data,
                   concentration_values = "nonexistent_col",
                   show_plots = FALSE),
    "not found"
  )
})

test_that("diffusive_flux uses cutoff_start_value = Inf as default", {
  data("DIY_sensor_data", package = "FluxSeparator")
  # With Inf default, should not filter out any cycles by start value
  result <- diffusive_flux(DIY_sensor_data, show_plots = FALSE)
  expect_s3_class(result, "data.frame")
})
