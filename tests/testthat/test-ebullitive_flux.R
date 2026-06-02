test_that("ebullitive_flux runs on example data", {
  data("DIY_sensor_data", package = "FluxSeparator")
  result <- ebullitive_flux(DIY_sensor_data, show_plots = FALSE)
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expected_cols <- c("station", "PumpCycle", "datetime_start", "datetime_end",
                     "sum_bubbles_concentration", "n_bubbles",
                     "pumpcycle_duration_hr", "temp",
                     "bubbles_per_time", "concentration_per_time")
  for (col in expected_cols) {
    expect_true(col %in% names(result),
                info = paste("Missing column:", col))
  }
})

test_that("ebullitive_flux top_selection works for both options", {
  data("DIY_sensor_data", package = "FluxSeparator")
  result_last <- ebullitive_flux(DIY_sensor_data,
                                  top_selection = "last",
                                  show_plots = FALSE)
  result_max <- ebullitive_flux(DIY_sensor_data,
                                 top_selection = "max",
                                 show_plots = FALSE)
  expect_s3_class(result_last, "data.frame")
  expect_s3_class(result_max, "data.frame")
})

test_that("ebullitive_flux rejects invalid top_selection", {
  data("DIY_sensor_data", package = "FluxSeparator")
  expect_error(
    ebullitive_flux(DIY_sensor_data, top_selection = "invalid",
                    show_plots = FALSE)
  )
})

test_that("ebullitive_flux validates required columns", {
  bad_data <- data.frame(x = 1:10, y = 1:10)
  expect_error(
    ebullitive_flux(bad_data, show_plots = FALSE),
    "Missing required columns"
  )
})

test_that("ebullitive_flux validates concentration column exists", {
  data("DIY_sensor_data", package = "FluxSeparator")
  expect_error(
    ebullitive_flux(DIY_sensor_data,
                    concentration_values = "nonexistent_col",
                    show_plots = FALSE),
    "not found"
  )
})
