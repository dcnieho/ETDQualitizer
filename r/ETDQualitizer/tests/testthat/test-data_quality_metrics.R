test_that("accuracy returns zero offset for symmetric gaze around target", {
  x <- c(0, 1, -1)
  y <- c(0, 1, -1)
  result <- accuracy(x, y, 0, 0)
  expect_equal(result$offset, 0, tolerance = 1e-8)
  expect_equal(result$offset_azi, 0, tolerance = 1e-8)
  expect_equal(result$offset_ele, 0, tolerance = 1e-8)
})

test_that("accuracy works with custom central tendency function", {
  x <- c(0, 1, -1)
  y <- c(0, 1, -1)
  result <- accuracy(x, y, 0, 0, central_tendency_fun = median)
  expect_equal(result$offset, 0, tolerance = 1e-8)
  expect_equal(result$offset_azi, 0, tolerance = 1e-8)
  expect_equal(result$offset_ele, 0, tolerance = 1e-8)
})

test_that("std returns correct values", {
  x <- c(1, 2, 3)
  y <- c(4, 5, 6)
  result <- std(x, y)
  expect_equal(result$std_azi, sqrt(mean((x - mean(x))^2)))
  expect_equal(result$std_ele, sqrt(mean((y - mean(y))^2)))
  expect_equal(result$std, sqrt(result$std_a^2 + result$std_e^2))
})

test_that("std returns correct values (with NA)", {
  x <- c(1, 2, NA, 3)
  y <- c(4, 5, NA, 6)
  result <- std(x, y)
  expect_equal(result$std_azi, sqrt(mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE)))
  expect_equal(result$std_ele, sqrt(mean((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)))
  expect_equal(result$std, sqrt(result$std_a^2 + result$std_e^2))
})

test_that("bcea returns valid area and aspect ratio", {
  set.seed(42)
  x <- rnorm(100000)
  y <- rnorm(100000)
  result <- bcea(x, y)
  expect_gt(result$area, 0)
  expect_equal(result$aspect_ratio, 1, tolerance = 0.1)
  expect_equal(result$area, 2 * pi * result$ax1 * result$ax2, tolerance = 1e-3)
})

test_that("rms_s2s returns correct RMS values", {
  x <- c(1, 2, 3)
  y <- c(4, 5, 6)
  result <- rms_s2s(x, y)
  expect_gte(result$rms, 0)
  expect_equal(result$rms_azi, sqrt(mean(diff(x)^2)))
  expect_equal(result$rms_ele, sqrt(mean(diff(y)^2)))
})

test_that("data_loss_from_invalid computes correct percentage", {
  x <- c(1, NA, 3)
  y <- c(4, 5, NA)
  loss <- data_loss_from_invalid(x, y)
  expect_equal(loss, 2 / 3 * 100)
})

test_that("data_loss_from_expected computes correct loss", {
  x <- c(1, NA, 3)
  y <- c(4, 5, NA)
  loss <- data_loss_from_expected(x, y, duration = 1, frequency = 3)
  expect_equal(loss, (1 - 1 / 3) * 100)
})

test_that("effective_frequency computes correct value", {
  x <- c(1, NA, 3)
  y <- c(4, 5, NA)
  freq <- effective_frequency(x, y, duration = 1)
  expect_equal(freq, 1)
})
