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

test_that("accuracy with median uses the Frechet median on the sphere", {
  x <- c(7.42, 73.96, 53.70, -84.53)
  y <- c(28.59, -37.31, 18.37, -25.95)
  grid_azi <- seq(-85, 85, length.out = 121)
  grid_ele <- seq(-40, 40, length.out = 101)
  gaze_vectors <- Fick_to_vector(x, y)
  sample_vectors <- cbind(gaze_vectors$x, gaze_vectors$y, gaze_vectors$z)

  objective <- function(target_azi, target_ele) {
    target_vector <- unlist(Fick_to_vector(target_azi, target_ele), use.names = FALSE)
    dots <- pmin(1, pmax(-1, as.vector(sample_vectors %*% target_vector)))
    sum(acos(dots))
  }

  best_value <- Inf
  best_target <- c(NA_real_, NA_real_)
  for (target_azi in grid_azi) {
    for (target_ele in grid_ele) {
      value <- objective(target_azi, target_ele)
      if (value < best_value) {
        best_value <- value
        best_target <- c(target_azi, target_ele)
      }
    }
  }

  result <- accuracy(x, y, best_target[1], best_target[2], central_tendency_fun = median)
  expect_lt(result$offset, 1)
  expect_lt(abs(result$offset_azi), 1)
  expect_lt(abs(result$offset_ele), 1)

  legacy_vector <- c(median(sample_vectors[, 1]), median(sample_vectors[, 2]), median(sample_vectors[, 3]))
  legacy_vector <- legacy_vector / sqrt(sum(legacy_vector^2))
  legacy_value <- sum(acos(pmin(1, pmax(-1, as.vector(sample_vectors %*% legacy_vector)))))
  expect_lt(best_value, legacy_value)
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
