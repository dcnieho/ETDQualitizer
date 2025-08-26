# Unit tests based on Python output -- not ideal and should be replaced.

# library(readr)
#
# mypath <- file.path(dirname(dirname(dirname(dirname(getwd())))), "test", "data", "Tobii_Spectrum_1200Hz.tsv")
#
# dat = read_tsv(file=mypath, n_max=2000)
# dat <- dat[dat$target_id==7,]
# #print(paste("Number of rows for this target: ", nrow(dat)))
#
#
# test_that("dummy test works", {
#   expect_equal(2 * 2, 4)
# })

test_that("pix_to_mm and mm_to_pix are consistent", {
  config <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
  mm <- config$pix_to_mm(960, 540)
  pix <- config$mm_to_pix(mm$x, mm$y)
  expect_equal(pix$x, 960, tolerance = 1e-10)
  expect_equal(pix$y, 540, tolerance = 1e-10)
})
test_that("mm_to_deg and deg_to_mm are consistent", {
  config <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
  deg <- config$mm_to_deg(250, 0)
  mm <- config$deg_to_mm(deg$azi, deg$ele)
  expect_equal(mm$x, 250, tolerance = 1e-10)
  expect_equal(mm$y, 0, tolerance = 1e-10)
})
test_that("screen_extents returns positive values", {
  config <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
  extents <- config$screen_extents()
  expect_gt(extents$x_deg, 0)
  expect_gt(extents$y_deg, 0)
})

test_that("Fick to vector and back is consistent", {
  azi <- 45
  ele <- 30
  vec <- Fick_to_vector(azi, ele)
  angles <- vector_to_Fick(vec$x, vec$y, vec$z)
  expect_equal(angles$azi, azi, tolerance = 1e-10)
  expect_equal(angles$ele, ele, tolerance = 1e-10)
})

test_that("accuracy returns valid offsets", {
  x <- c(0, 1, -1)
  y <- c(0, 1, -1)
  result <- accuracy(x, y, 0, 0)
  expect_gte(result$offset, 0)
})
test_that("std_ returns correct values", {
  x <- c(1, 2, 3)
  y <- c(4, 5, 6)
  result <- std_(x, y)
  expect_equal(result$std_, sqrt(var(x, na.rm = TRUE) + var(y, na.rm = TRUE)))
})
test_that("bcea returns positive area", {
  x <- rnorm(100)
  y <- rnorm(100)
  result <- bcea(x, y)
  expect_gt(result$area, 0)
})
test_that("rms_s2s returns valid RMS", {
  x <- c(1, 2, 3, 4)
  y <- c(4, 5, 6, 7)
  result <- rms_s2s(x, y)
  expect_gte(result$rms, 0)
})
test_that("data_loss returns correct percentage", {
  x <- c(1, NA, 3)
  y <- c(4, 5, NA)
  loss <- data_loss(x, y)
  expect_equal(loss, 2/3 * 100)
})
test_that("effective_frequency returns correct value", {
  x <- c(1, NA, 3)
  y <- c(4, 5, NA)
  freq <- effective_frequency(x, y, duration = 1)
  expect_equal(freq, 1)
})

test_that("DataQuality methods work", {
  timestamps <- seq(0, 0.99, by = 0.01)
  azi <- rnorm(100)
  ele <- rnorm(100)
  dq <- DataQuality$new(azi, ele, timestamps, unit = "degrees")
  expect_equal(length(dq$azi), 100)
  expect_equal(length(dq$ele), 100)
  acc <- dq$accuracy(0, 0)
  expect_gte(acc$offset, 0)
  rms <- dq$precision_RMS_S2S()
  expect_gte(rms$rms, 0)
  std <- dq$precision_STD()
  expect_equal(std$std__, sqrt(var(azi) + var(ele)))
  bcea <- dq$precision_BCEA()
  expect_gt(bcea$area, 0)
  loss <- dq$data_loss()
  expect_equal(loss, 0)
  freq <- dq$effective_frequency()
  expect_equal(freq, 100)
  duration <- dq$get_duration()
  expect_equal(duration, 0.99 + 0.01)
})
