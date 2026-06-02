test_that("ppm_to_umol computes correct conversion", {
  # At standard conditions: 101325 Pa, 20°C
  result <- ppm_to_umol(
    pressure = 101325,
    concentration = 10,
    volume = 0.01,
    temperature_C = 20,
    area = 0.05
  )
  expect_type(result, "double")
  expect_length(result, 1)
  # Manual calculation: (101325 * 10 * 0.01) / (8.314 * 293.15) / 0.05
  expected <- (101325 * 10 * 0.01) / (8.314 * 293.15) / 0.05
  expect_equal(result, expected, tolerance = 1e-6)
})

test_that("ppm_to_umol is vectorized", {
  result <- ppm_to_umol(
    pressure = c(101325, 101325),
    concentration = c(10, 20),
    volume = 0.01,
    temperature_C = 20,
    area = 0.05
  )
  expect_length(result, 2)
  expect_equal(result[2], result[1] * 2, tolerance = 1e-6)
})

test_that("DIY_sensor_data loads correctly", {
  data("DIY_sensor_data", package = "FluxSeparator")
  expect_s3_class(DIY_sensor_data, "data.frame")
  expect_true(nrow(DIY_sensor_data) > 0)
  expect_true("pred_CH4" %in% names(DIY_sensor_data))
  expect_true("PumpCycle" %in% names(DIY_sensor_data))
  expect_true("station" %in% names(DIY_sensor_data))
  expect_true("datetime" %in% names(DIY_sensor_data))
})
