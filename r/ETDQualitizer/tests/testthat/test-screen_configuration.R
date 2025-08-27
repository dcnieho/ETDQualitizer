test_that("constructor sets fields correctly", {
  sc <- ScreenConfiguration$new(400, 250, 1600, 900, 500)
  expect_equal(sc$screen_size_x_mm, 400)
  expect_equal(sc$screen_size_y_mm, 250)
  expect_equal(sc$screen_res_x_pix, 1600)
  expect_equal(sc$screen_res_y_pix, 900)
  expect_equal(sc$viewing_distance_mm, 500)
})

test_that("pix_to_mm conversion is correct", {
  config <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
  mm <- config$pix_to_mm(960, 540)
  expect_equal(mm$x, 250, tolerance = 1e-8)
  expect_equal(mm$y, 150, tolerance = 1e-8)
})

test_that("mm_to_pix conversion is correct", {
  config <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
  pix <- config$mm_to_pix(250, 150)
  expect_equal(pix$x, 960, tolerance = 1e-8)
  expect_equal(pix$y, 540, tolerance = 1e-8)
})

test_that("mm_to_deg conversion is correct", {
  config <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
  deg <- config$mm_to_deg(250, 0)
  expect_equal(deg$azi, atan2(250, 600) * 180 / pi, tolerance = 1e-8)
  expect_equal(deg$ele, 0, tolerance = 1e-8)
})

test_that("deg_to_mm conversion is correct", {
  config <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
  azi <- atan2(250, 600) * 180 / pi
  ele <- 0
  mm <- config$deg_to_mm(azi, ele)
  expect_equal(mm$x, 250, tolerance = 1e-8)
  expect_equal(mm$y, 0, tolerance = 1e-8)
})

test_that("pix_to_deg returns correct azimuth and elevation", {
  config <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
  deg <- config$pix_to_deg(960, 540)
  expected_azi <- atan2(250, 600) * 180 / pi
  expected_ele <- atan2(150, sqrt(600^2 + 250^2)) * 180 / pi
  expect_equal(deg$azi, expected_azi, tolerance = 1e-10)
  expect_equal(deg$ele, expected_ele, tolerance = 1e-10)
})

test_that("deg_to_pix conversion is correct", {
  config <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
  azi <- atan2(250, 600) * 180 / pi
  ele <- atan2(150, sqrt(600^2 + 250^2)) * 180 / pi
  pix <- config$deg_to_pix(azi, ele)
  expect_equal(pix$x, 960, tolerance = 1e-1)
  expect_equal(pix$y, 540, tolerance = 1e-1)
})

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
  expect_equal(mm$y,   0, tolerance = 1e-10)
})

test_that("screen_extents returns correct values", {
  config <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
  extents <- config$screen_extents()
  expected_x_deg <- 2 * atan2(250, 600) * 180 / pi
  expected_y_deg <- 2 * atan2(150, 600) * 180 / pi
  expect_equal(extents$width, expected_x_deg, tolerance = 1e-8)
  expect_equal(extents$height, expected_y_deg, tolerance = 1e-8)
})
