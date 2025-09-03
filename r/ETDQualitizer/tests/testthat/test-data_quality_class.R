test_that("constructor with degrees stores data correctly", {
  duration <- 1000
  freq <- 100
  timestamps <- seq(0, duration - 1 / freq, by = 1 / freq)
  n_samples <- length(timestamps)
  azi <- rnorm(n_samples)
  ele <- rnorm(n_samples)

  dq <- DataQuality$new(azi, ele, timestamps, "degrees")
  expect_equal(dq$azi, azi)
  expect_equal(dq$ele, ele)
  expect_equal(dq$timestamps, timestamps)
})

test_that("constructor with pixels converts correctly", {
  duration <- 1000
  freq <- 100
  timestamps <- seq(0, duration - 1 / freq, by = 1 / freq)
  n_samples <- length(timestamps)
  azi <- rnorm(n_samples)
  ele <- rnorm(n_samples)
  screen <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)

  x_pix <- 960 + azi * 10
  y_pix <- 540 + ele * 10
  dq <- DataQuality$new(x_pix, y_pix, timestamps, "pixels", screen)
  deg <- screen$pix_to_deg(x_pix, y_pix)
  expect_equal(dq$azi, deg$azi, tolerance = 1e-8)
  expect_equal(dq$ele, deg$ele, tolerance = 1e-8)
})

test_that("accuracy returns near-zero offset for symmetric data", {
  dq <- DataQuality$new(c(0, 1, -1), c(0, 1, -1), c(0, 1, 2), "degrees")
  result <- dq$accuracy(0, 0)
  expect_equal(result$offset, 0, tolerance = 0.1)
  expect_equal(result$offset_azi, 0, tolerance = 0.1)
  expect_equal(result$offset_ele, 0, tolerance = 0.1)
})

test_that("precision RMS with large dataset", {
  dq <- DataQuality$new(rnorm(100000), rnorm(100000), seq(0, 1000, length.out = 100000), "degrees")
  result <- dq$precision_RMS_S2S()
  expect_equal(result$rms, 2, tolerance = 0.1)
  expect_equal(result$rms_azi, sqrt(2), tolerance = 0.1)
  expect_equal(result$rms_ele, sqrt(2), tolerance = 0.1)
})

test_that("precision RMS with one axis varying", {
  dq <- DataQuality$new(c(1, 2, 4), c(1, 1, 1), c(0, 1, 2), "degrees")
  result <- dq$precision_RMS_S2S()
  expected_rms <- sqrt(mean(c(1, 2)^2))
  expect_gte(result$rms, expected_rms)
  expect_gte(result$rms_azi, expected_rms)
  expect_gte(result$rms_ele, 0)
})

test_that("precision RMS with both axes varying", {
  dq <- DataQuality$new(c(1, 2, 4), c(1, 2, 4), c(0, 1, 2), "degrees")
  result <- dq$precision_RMS_S2S()
  expected_rms_xy <- sqrt(mean(c(1, 2)^2))
  expected_rms <- sqrt(2) * expected_rms_xy
  expect_gte(result$rms, expected_rms)
  expect_gte(result$rms_azi, expected_rms_xy)
  expect_gte(result$rms_ele, expected_rms_xy)
})

test_that("precision STD returns expected values", {
  dq <- DataQuality$new(rnorm(100000, sd = 1), rnorm(100000, sd = 1), seq(0, 1000, length.out = 100000), "degrees")
  result <- dq$precision_STD()
  expect_equal(result$std, sqrt(2), tolerance = 0.1)
  expect_equal(result$std_azi, 1, tolerance = 0.1)
  expect_equal(result$std_ele, 1, tolerance = 0.1)
})

test_that("precision BCEA returns valid values", {
  dq <- DataQuality$new(rnorm(10000), rnorm(10000), seq(0, 100, length.out = 10000), "degrees")
  result <- dq$precision_BCEA()
  expect_gt(result$area, 0)
  expect_equal(result$aspect_ratio, 1, tolerance = 0.1)
  expect_equal(result$area, 2 * pi * result$ax1 * result$ax2, tolerance = 1e-3)
})

test_that("precision STD using moving window returns reasonable value", {
  dq <- DataQuality$new(rnorm(10000), rnorm(10000), seq(0, 100, length.out = 10000), "degrees")
  s <- dq$precision_using_moving_window(50, "STD")
  expect_lte(s, sqrt(2))
})

test_that("data loss is zero for complete data", {
  dq <- DataQuality$new(rnorm(1000), rnorm(1000), seq(0, 10, length.out = 1000), "degrees")
  expect_equal(dq$data_loss(), 0)
})

test_that("data loss from expected is zero for complete data", {
  freq <- 100
  timestamps <- seq(0, 10, by = 1 / freq)
  dq <- DataQuality$new(rnorm(length(timestamps)), rnorm(length(timestamps)), timestamps, "degrees")
  expect_equal(dq$data_loss_from_expected(freq), 0)
})

test_that("effective frequency matches expected", {
  freq <- 100
  timestamps <- seq(0, 10 - 1 / freq, by = 1 / freq)
  dq <- DataQuality$new(rnorm(length(timestamps)), rnorm(length(timestamps)), timestamps, "degrees")
  expect_equal(dq$effective_frequency(), freq)
})

test_that("get_duration returns correct value", {
  duration <- 1000
  freq <- 100
  timestamps <- seq(0, duration - 1 / freq, by = 1 / freq)
  dq <- DataQuality$new(rnorm(length(timestamps)), rnorm(length(timestamps)), timestamps, "degrees")
  expect_equal(dq$get_duration(), duration, tolerance = 1e-6)
})