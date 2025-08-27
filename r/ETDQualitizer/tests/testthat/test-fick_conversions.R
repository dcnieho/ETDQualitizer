test_that("Fick_to_vector with default radius", {
  vec <- Fick_to_vector(0, 0)
  expect_equal(vec$x, 0, tolerance = 1e-8)
  expect_equal(vec$y, 0, tolerance = 1e-8)
  expect_equal(vec$z, 1, tolerance = 1e-8)
})

test_that("Fick_to_vector with custom radius", {
  vec <- Fick_to_vector(90, 0, 2)
  expect_equal(vec$x, 2, tolerance = 1e-8)
  expect_equal(vec$y, 0, tolerance = 1e-8)
  expect_equal(vec$z, 0, tolerance = 1e-8)
})

test_that("vector_to_Fick basic conversion", {
  fick <- vector_to_Fick(0, 0, 1)
  expect_equal(fick$azi, 0, tolerance = 1e-8)
  expect_equal(fick$ele, 0, tolerance = 1e-8)
})

test_that("Fick round-trip conversion", {
  azi <- 45
  ele <- 30
  r <- 1.5
  vec <- Fick_to_vector(azi, ele, r)
  fick2 <- vector_to_Fick(vec$x, vec$y, vec$z)
  expect_equal(fick2$azi, azi, tolerance = 1e-5)
  expect_equal(fick2$ele, ele, tolerance = 1e-5)
})

test_that("vector_to_Fick with negative components", {
  fick <- vector_to_Fick(-1, -1, -1)
  expected_azi <- atan2(-1, -1) * 180 / pi
  expected_ele <- atan2(-1, sqrt((-1)^2 + (-1)^2)) * 180 / pi
  expect_equal(fick$azi, expected_azi, tolerance = 1e-8)
  expect_equal(fick$ele, expected_ele, tolerance = 1e-8)
})
