library(R6)
library(stats)


#' Get ETDQualitizer Version
#'
#' Returns the current version string of the ETDQualitizer tool.
#'
#' @return A character string representing the version number.
#'
#' @examples
#' ETDQ_version()
#'
#' @export
ETDQ_version <- function() {
  return("0.9.0")
}


#' Convert Fick Angles to 3D Vector
#'
#' Converts azimuth and elevation angles (in degrees) to a 3D unit vector.
#'
#' @param azi Azimuth angle in degrees.
#' @param ele Elevation angle in degrees.
#' @param rho Radius (default is 1.0).
#'
#' @return A list with components \code{x}, \code{y}, and \code{z}.
#' @examples
#' Fick_to_vector(30, 10)
#' @export
Fick_to_vector <- function(azi, ele, rho=1.0) {
  # Convert degrees to radians
  azi_rad <- pi/180 * azi
  ele_rad <- pi/180 * ele
  r_cos_ele <- rho * cos(ele_rad)
  x <- r_cos_ele * sin(azi_rad)
  y <-       rho * sin(ele_rad)
  z <- r_cos_ele * cos(azi_rad)
  list(x = x, y = y, z = z)
}

#' Convert 3D Vector to Fick Angles
#'
#' Converts a 3D vector to azimuth and elevation angles (in degrees).
#'
#' @param x X component of the vector.
#' @param y Y component of the vector.
#' @param z Z component of the vector.
#'
#' @return A list with components \code{azi} and \code{ele}.
#' @examples
#' vector_to_Fick(0.5, 0.2, 0.8)
#' @export
vector_to_Fick <- function(x, y, z) {
  azi <- 180/pi * atan2(x, z)
  ele <- 180/pi * atan2(y, sqrt(x^2 + z^2))
  list(azi = azi, ele = ele)
}


#' Compute Gaze Accuracy
#'
#' Calculates the angular offset between gaze and target directions.
#'
#' @param azi Gaze azimuth in degrees.
#' @param ele Gaze elevation in degrees.
#' @param target_azi Target azimuth in degrees.
#' @param target_ele Target elevation in degrees.
#' @param central_tendency_fun Function to compute central tendency (default: \code{mean}).
#'
#' @return A list with \code{offset}, \code{offset_azi}, and \code{offset_ele}, the total, horizontal and vertical offset of gaze from the target (in degrees).
#' @examples
#' accuracy(c(1, 2), c(1, 2), 0, 0)
#' @export
accuracy <- function(azi, ele, target_azi, target_ele, central_tendency_fun = mean) {
  # Get unit vectors from gaze directions
  g <- Fick_to_vector(azi, ele)
  t <- Fick_to_vector(target_azi, target_ele)
  # calculate angular offset for each sample using dot product
  dot_products <- g$x*t$x + g$y*t$y + g$z*t$z
  dot_products <- pmin(pmax(dot_products, -1), 1)  # Clamp to [-1,1]
  offsets      <- acos(dot_products)
  # calculate on-screen orientation so we can decompose offset into x and y
  direction    <- atan2(g$y/g$z-t$y/t$z, g$x/g$z-t$x/t$z)
  offsets_2D   <- offsets*cbind(cos(direction), sin(direction)) * 180/pi
  # calculate mean horizontal and vertical offsets
  offset_azi   <- central_tendency_fun(offsets_2D[, 1], na.rm = TRUE)
  offset_ele   <- central_tendency_fun(offsets_2D[, 2], na.rm = TRUE)
  # calculate offset of centroid
  offset       <- sqrt(offset_azi^2 + offset_ele^2)
  list(offset = offset, offset_azi = offset_azi, offset_ele = offset_ele)
}

#' RMS of Sample-to-Sample Differences
#'
#' Computes root mean square of differences between successive gaze samples.
#'
#' @param azi Azimuth values in degrees.
#' @param ele Elevation values in degrees.
#' @param central_tendency_fun Function to compute central tendency (default: \code{mean}).
#'
#' @return A list with \code{rms}, \code{rms_azi}, and \code{rms_ele}, the total RMS of sample-to-sample distances and that of the azimuthal and elevation components (all in degrees).
#' @examples
#' rms_s2s(c(1, 2, 3), c(1, 2, 3))
#' @export
rms_s2s <- function(azi, ele, central_tendency_fun = mean) {
  a_diff  <- diff(azi)^2
  e_diff  <- diff(ele)^2
  rms_azi <- sqrt(central_tendency_fun(a_diff, na.rm = TRUE))
  rms_ele <- sqrt(central_tendency_fun(e_diff, na.rm = TRUE))
  rms     <- sqrt(central_tendency_fun(a_diff + e_diff, na.rm = TRUE))
  list(rms = rms, rms_azi = rms_azi, rms_ele = rms_ele)
}

#' Standard Deviation of Gaze Samples
#'
#' Computes standard deviation of azimuth and elevation.
#'
#' @param azi Azimuth values in degrees.
#' @param ele Elevation values in degrees.
#'
#' @return A list with \code{std}, \code{std_azi}, and \code{std_ele}, the total STD of sample-to-sample distances and that of the azimuthal and elevation components (all in degrees).
#' @examples
#' std(c(1, 2, 3), c(1, 2, 3))
#' @export
std <- function(azi, ele) {
  std_azi <- pop_sd(azi)
  std_ele <- pop_sd(ele)
  std     <- sqrt(std_azi^2 + std_ele^2)
  list(std = std, std_azi = std_azi, std_ele = std_ele)
}
#' Population Standard Deviation
#'
#' Computes population standard deviation of a numeric vector.
#'
#' @param x Numeric vector.
#'
#' @return Population standard deviation.
#' @examples
#' pop_sd(c(1, 2, 3, 4))
#' @noRd
pop_sd <- function(x) {
  x <- x[!is.na(x)]
  sqrt(mean((x - mean(x))^2))
}

#' Bivariate Contour Ellipse Area (BCEA)
#'
#' Computes BCEA and ellipse parameters for gaze precision.
#'
#' @param azi Azimuth values in degrees.
#' @param ele Elevation values in degrees.
#' @param P Cumulative probability (default: 0.68).
#'
#' @return A list with the BCEA (\code{area}) and additional info about the BCEA ellipse: \code{orientation}, \code{ax1}, \code{ax2}, and \code{aspect_ratio}.
#' @examples
#' bcea(rnorm(100), rnorm(100))
#' @export
#' @importFrom stats cor cov sd
bcea <- function(azi, ele, P = 0.68) {
  # turn cumulative probability into scale factor
  k <- log(1 / (1 - P))
  # remove NaNs
  valid <- !(is.na(azi) | is.na(ele))
  azi_valid <- azi[valid]
  ele_valid <- ele[valid]
  # compute sample standard deviations
  std_a <- sd(azi_valid)
  std_e <- sd(ele_valid)
  # compute correlation
  rho <- cor(azi_valid, ele_valid)
  # compute BCEA
  area <- 2 * k * pi * std_a * std_e * sqrt(1 - rho^2)

  # compute major and minor axis radii, and orientation, of the BCEA ellipse
  cov_matrix <- cov(cbind(azi_valid, ele_valid))
  eig <- eigen(cov_matrix)
  i <- which.max(eig$values)
  orientation <- atan2(eig$vectors[2, i], eig$vectors[1, i]) * 180 / pi
  ax1 <- sqrt(k * eig$values[ i])
  ax2 <- sqrt(k * eig$values[-i])
  aspect_ratio <- max(ax1, ax2) / min(ax1, ax2)
  # sanity check: this (formula for area of ellipse) should
  # closely match directly computed area from above
  # 2*pi*ax1*ax2
  list(
    area = area,
    orientation = orientation,
    ax1 = ax1,
    ax2 = ax2,
    aspect_ratio = aspect_ratio
  )
}

#' Compute Data Loss from number of invalid samples.
#'
#' Calculates percentage of missing gaze samples.
#'
#' @param a Horizontal gaze values (e.g. azimuth or horizontal coordinate in pixels or mm).
#' @param b Vertical gaze values (e.g. azimuth or horizontal coordinate in pixels or mm).
#'
#' @return Percentage of missing samples.
#' @examples
#' data_loss_from_invalid(c(1, NA, 3), c(1, 2, NA))
#' @export
data_loss_from_invalid <- function(a, b) {
  missing <- is.na(a) | is.na(b)
  sum(missing)/length(missing)*100
}

#' Compute Data Loss from Expected Sample Count
#'
#' Calculates data loss based on expected number of samples.
#'
#' @param a Horizontal gaze values (e.g. azimuth or horizontal coordinate in pixels or mm).
#' @param b Vertical gaze values (e.g. azimuth or horizontal coordinate in pixels or mm).
#' @param duration Duration in seconds.
#' @param frequency Sampling frequency in Hz.
#'
#' @return Percentage of data loss.
#' @examples
#' data_loss_from_expected(c(1, NA, 3), c(1, 2, NA), duration = 1, frequency = 3)
#' @export
data_loss_from_expected <- function(a, b, duration, frequency) {
  N_valid <- sum(!is.na(a) & !is.na(b))
  (1 - N_valid/(duration* frequency))*100
}

#' Compute Effective Sampling Frequency
#'
#' Calculates effective frequency based on valid samples.
#'
#' @param a Horizontal gaze values (e.g. azimuth or horizontal coordinate in pixels or mm).
#' @param b Vertical gaze values (e.g. azimuth or horizontal coordinate in pixels or mm).
#' @param duration Duration in seconds.
#'
#' @return Effective frequency in Hz.
#' @examples
#' effective_frequency(c(1, NA, 3), c(1, 2, NA), duration = 1)
#' @export
effective_frequency <- function(a, b, duration) {
  N_valid <- sum(!is.na(a) & !is.na(b))
  N_valid/duration
}


#' Precision Using Moving Window
#'
#' Computes gaze precision using a moving window and selected metric.
#'
#' @param azi Azimuth values.
#' @param ele Elevation values.
#' @param window_length Window size in samples.
#' @param metric Precision metric: \code{"RMS-S2S"}, \code{"STD"}, or \code{"BCEA"}.
#' @param aggregation_fun Function to aggregate precision values across the windows (default: \code{median}).
#' @param ... Additional arguments passed to metric function.
#'
#' @return Aggregated precision value.
#' @examples
#' precision_using_moving_window(rnorm(100), rnorm(100), 10, "STD")
#' @export
#' @importFrom stats median
precision_using_moving_window <- function(azi, ele, window_length, metric, aggregation_fun = median, ...) {
  # Select the appropriate precision metric function
  fun <- switch(metric,
    "RMS-S2S" = rms_s2s,
    "STD"     = std,
    "BCEA"    = bcea,
    stop(sprintf('metric "%s" is not understood', metric))
  )

  # Get number of samples in data
  ns <- length(azi)

  if (window_length < ns) {
    # If number of samples in data exceeds window size
    values <- rep(NA_real_, ns - window_length + 1)  # pre-allocate
    for (p in seq_len(ns - window_length + 1)) {
      result <- fun(azi[p:(p + window_length - 1)], ele[p:(p + window_length - 1)], ...)
      values[p] <- result[[1]]  # extract first element (e.g., rms, std, or area)
    }
    precision <- aggregation_fun(values, na.rm = TRUE)
  } else {
    # If too few samples in data
    precision <- NA_real_
  }

  precision
}


#' @title R6 Screen Configuration Class
#'
#' @description Provides methods for converting between pixel, millimeter, and degree units.
#'
#' @examples
#' sc <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
#' sc$pix_to_deg(960, 540)
#' @export
#' @importFrom R6 R6Class
ScreenConfiguration <- R6Class("ScreenConfiguration",
  public = list(
    #' @field screen_size_x_mm Screen width in mm.
    #' @field screen_size_y_mm Screen height in mm.
    #' @field screen_res_x_pix Horizontal screen resolution in pixels.
    #' @field screen_res_y_pix Vertical screen resolution in pixels.
    #' @field viewing_distance_mm Viewing distance in mm.
    screen_size_x_mm = NULL,
    screen_size_y_mm = NULL,
    screen_res_x_pix = NULL,
    screen_res_y_pix = NULL,
    viewing_distance_mm = NULL,

    #' @description Creates a new ScreenConfiguration object with screen and viewing distance parameters.
    #' @param screen_size_x_mm Screen width in millimeters.
    #' @param screen_size_y_mm Screen height in millimeters.
    #' @param screen_res_x_pix Horizontal screen resolution in pixels.
    #' @param screen_res_y_pix Vertical screen resolution in pixels.
    #' @param viewing_distance_mm Viewing distance in millimeters.
    #' @return A new ScreenConfiguration object.
    #' @examples
    #' sc <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
    initialize = function(screen_size_x_mm, screen_size_y_mm, screen_res_x_pix, screen_res_y_pix, viewing_distance_mm) {
      self$screen_size_x_mm <- screen_size_x_mm
      self$screen_size_y_mm <- screen_size_y_mm
      self$screen_res_x_pix <- screen_res_x_pix
      self$screen_res_y_pix <- screen_res_y_pix
      self$viewing_distance_mm <- viewing_distance_mm
    },

    #' @description Converts pixel coordinates to millimeter coordinates on the screen.
    #' @param x Horizontal pixel coordinate.
    #' @param y Vertical pixel coordinate.
    #' @return A list with x and y in millimeters.
    #' @examples
    #' sc$pix_to_mm(960, 540)
    pix_to_mm = function(x, y) {
      x_mm <- x / self$screen_res_x_pix * self$screen_size_x_mm
      y_mm <- y / self$screen_res_y_pix * self$screen_size_y_mm
      list(x = x_mm, y = y_mm)
    },

    #' @description Converts pixel coordinates to an angular gaze direction in degrees.
    #' @param x Horizontal pixel coordinate.
    #' @param y Vertical pixel coordinate.
    #' @return A list with azimuth (\code{"azi"}) and elevation (\code{"ele"}) in degrees.
    #' @examples
    #' sc$pix_to_deg(960, 540)
    pix_to_deg = function(x, y) {
      mm <- self$pix_to_mm(x, y)
      self$mm_to_deg(mm$x, mm$y)
    },

    #' @description Converts millimeter coordinates to an angular gaze direction in degrees.
    #' @param x Horizontal position in millimeters.
    #' @param y Vertical position in millimeters.
    #' @return A list with azimuth (\code{"azi"}) and elevation (\code{"ele"}) in degrees.
    #' @examples
    #' sc$mm_to_deg(100, 50)
    mm_to_deg = function(x, y) {
      azi <- atan2(x, self$viewing_distance_mm) * 180 / pi
      ele <- atan2(y, sqrt(self$viewing_distance_mm^2 + x^2)) * 180 / pi
      list(azi = azi, ele = ele)
    },

    #' @description Converts millimeter coordinates on the screen to pixel coordinates.
    #' @param x Horizontal position in millimeters.
    #' @param y Vertical position in millimeters.
    #' @return A list with x and y in pixels.
    #' @examples
    #' sc$mm_to_pix(100, 50)
    mm_to_pix = function(x, y) {
      x_pix <- x / self$screen_size_x_mm * self$screen_res_x_pix
      y_pix <- y / self$screen_size_y_mm * self$screen_res_y_pix
      list(x = x_pix, y = y_pix)
    },

    #' @description Converts an angular gaze direction in degrees to millimeter coordinates on the screen.
    #' @param azi Azimuth in degrees (Fick angles).
    #' @param ele Elevation in degrees (Fick angles).
    #' @return A list with x and y in millimeters.
    #' @examples
    #' sc$deg_to_mm(2, 1)
    deg_to_mm = function(azi, ele) {
      x_mm <- self$viewing_distance_mm * tan(azi * pi / 180)
      y_mm <- self$viewing_distance_mm * tan(ele * pi / 180) / cos(azi * pi / 180)
      list(x = x_mm, y = y_mm)
    },

    #' @description Converts an angular gaze direction in degrees to pixel coordinates.
    #' @param azi Azimuth in degrees (Fick angles).
    #' @param ele Elevation in degrees (Fick angles).
    #' @return A list with x and y in pixels.
    #' @examples
    #' sc$deg_to_pix(2, 1)
    deg_to_pix = function(azi, ele) {
      mm <- self$deg_to_mm(azi, ele)
      self$mm_to_pix(mm$x, mm$y)
    },

    #' @description Computes the horizontal and vertical extents of the screen (in degrees).
    #' @return A list with width and height in degrees.
    #' @examples
    #' sc$screen_extents()
    screen_extents = function() {
      w <- self$mm_to_deg(self$screen_size_x_mm / 2, 0.)
      h <- self$mm_to_deg(0., self$screen_size_y_mm / 2)
      list(width = w$azi * 2, height = h$ele * 2)
    }
  )
)


#' @title R6 class for calculating Data Quality from a gaze data segment
#'
#' @description Provides methods for assessing the quality of gaze data, including accuracy, precision, data loss, and effective sampling frequency.
#'
#' @examples
#' sc <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
#' gaze_x <- c(0, 1, -1)
#' gaze_y <- c(0, 1, -1)
#' timestamps <- c(0, 1, 2)
#' dq <- DataQuality$new(gaze_x, gaze_y, timestamps, unit = "pixels", screen = sc)
#' dq$accuracy(0, 0)
#' dq$precision_RMS_S2S()
#' dq$data_loss_from_invalid()
#'
#' @export
#' @importFrom R6 R6Class
DataQuality <- R6Class("DataQuality",
  public = list(
    #' @field timestamps Vector of timestamps in seconds. Samples with missing data should not be removed, or the RMS calculation would be incorrect.
    #' @field azi Vector of azimuth angles in degrees (Fick angles). Missing data should be coded as NA, not using some special value such as (0,0) or (-xres,-yres).
    #' @field ele Vector of elevation angles in degrees (Fick angles). Missing data should be coded as NA, not using some special value such as (0,0) or (-xres,-yres).
    timestamps = NULL,
    azi = NULL,
    ele = NULL,

    #' @description Creates a new DataQuality object from gaze data and timestamps.
    #' @param gaze_x Horizontal gaze positions (pixels or degrees).
    #' @param gaze_y Vertical gaze positions (pixels or degrees).
    #' @param timestamps Vector of timestamps in seconds.
    #' @param unit Unit of gaze data: either \code{"pixels"} or \code{"degrees"}.
    #' @param screen Optional \code{ScreenConfiguration} object, required if unit is \code{"pixels"}.
    #' @return A new DataQuality object.
    #' @examples
    #' dq <- DataQuality$new(gaze_x, gaze_y, timestamps, unit = "pixels", screen = sc)
    initialize = function(gaze_x, gaze_y, timestamps, unit, screen = NULL) {
      self$timestamps <- as.numeric(timestamps)

      gaze_x <- as.numeric(gaze_x)
      gaze_y <- as.numeric(gaze_y)

      if (unit == "pixels") {
        if (is.null(screen)) {
          stop('If unit is "pixels", a screen configuration must be supplied')
        }
        deg <- screen$pix_to_deg(gaze_x, gaze_y)
        gaze_x <- deg$azi
        gaze_y <- deg$ele
      } else if (unit != "degrees") {
        stop('unit should be "pixels" or "degrees"')
      }

      self$azi <- gaze_x
      self$ele <- gaze_y
    },

    #' @description Calculates the accuracy of gaze data relative to a known target location.
    #' @param target_azi Target azimuth in degrees.
    #' @param target_ele Target elevation in degrees.
    #' @param central_tendency_fun Function to compute central tendency (e.g., \code{mean}, \code{median}).
    #' @return Accuracy in degrees.
    #' @examples
    #' dq$accuracy(0, 0)
    accuracy = function(target_azi, target_ele, central_tendency_fun = mean) {
      accuracy(self$azi, self$ele, target_azi, target_ele, central_tendency_fun)
    },

    #' @description Calculates precision as root mean square of sample-to-sample distances
    #' @param central_tendency_fun Function to compute central tendency (e.g., \code{mean}, \code{median}).
    #' @return Precision in degrees.
    #' @examples
    #' dq$precision_RMS_S2S()
    precision_RMS_S2S = function(central_tendency_fun = mean) {
      rms_s2s(self$azi, self$ele, central_tendency_fun)
    },

    #' @description Calculates precision as standard deviation of gaze positions.
    #' @return Standard deviation in degrees.
    #' @examples
    #' dq$precision_STD()
    precision_STD = function() {
      std(self$azi, self$ele)
    },

    #' @description Calculates the Bivariate Contour Ellipse Area (BCEA) and ellipse parameters for gaze precision.
    #' @param P Proportion of data to include in the ellipse (default is 0.68).
    #' @return BCEA in degrees-squared.
    #' @examples
    #' dq$precision_BCEA()
    precision_BCEA = function(P = 0.68) {
      bcea(self$azi, self$ele, P)
    },

    #' @description Calculates the proportion of missing data (coded as NA).
    #' @return Proportion of missing samples.
    #' @examples
    #' dq$data_loss_from_invalid()
    data_loss_from_invalid = function() {
      data_loss_from_invalid(self$azi, self$ele)
    },

    #' @description Estimates data loss based on expected number of samples given the duration and sampling frequency.
    #' @param frequency Expected sampling frequency in Hz.
    #' @return Proportion of missing samples.
    #' @examples
    #' dq$data_loss_from_expected(500)
    data_loss_from_expected = function(frequency) {
      data_loss_from_expected(self$azi, self$ele, self$get_duration(), frequency)
    },

    #' @description Calculates the effective sampling frequency based on timestamps.
    #' @return Effective frequency in Hz.
    #' @examples
    #' dq$effective_frequency()
    effective_frequency = function() {
      effective_frequency(self$azi, self$ele, self$get_duration())
    },

    #' @description Computes the total duration of the gaze recording, including the last sample.
    #' @return Duration in seconds.
    #' @examples
    #' dq$get_duration()
    get_duration = function() {
      # to get duration right, we need to include duration of last sample
      isi <- median(diff(self$timestamps), na.rm = TRUE)
      self$timestamps[length(self$timestamps)] - self$timestamps[1] + isi
    },

    #' @description Calculates precision using a moving window approach.
    #' @param window_length Length of the moving window in number of samples.
    #' @param metric Precision metric to use (\code{"RMS-S2S"}, \code{"STD"}, or \code{"BCEA"}).
    #' @param aggregation_fun Function to aggregate windowed precision values (e.g., \code{median}).
    #' @param ... Additional arguments passed to the precision metric function.
    #' @return Precision value.
    #' @examples
    #' dq$precision_using_moving_window(0.2, "RMS-S2S")
    precision_using_moving_window = function(window_length, metric, aggregation_fun = median, ...) {
      precision_using_moving_window(self$azi, self$ele, window_length, metric, aggregation_fun, ...)
    }
  )
)


#' Compute Data Quality Metrics from Validation Data
#'
#' This function computes a set of data quality metrics for gaze data collected during the PsychoPy validation procedure
#' that is provided in the ETDQualitizer repository on github
#' (https://github.com/dcnieho/ETDQualitizer/tree/master/python/ETDQualitizer/stim).
#' It evaluates accuracy, precision, and optionally data loss and effective sampling frequency, per eye and per target.
#'
#' @param gaze A `data.frame` containing gaze data. Must include columns `target_id`, `tar_x`, `tar_y`, `timestamp`,
#'   and eye-specific columns such as `left_x`, `left_y`, `right_x`, `right_y`. Timestamps should be provided in milliseconds.
#' @param unit A character string specifying the unit of measurement for gaze and target coordinates in the gaze data.frame.
#'   Must be either `"pixels"` or `"degrees"`.
#' @param screen An optional `ScreenConfiguration` object or numeric scalar used to convert pixel coordinates to degrees.
#'   Required if `unit == "pixels"`.
#' @param advanced Logical. If `TRUE`, all available metrics are returned. If `FALSE`, only a simplified subset is included (default is FALSE).
#' @param include_data_loss Logical. If `TRUE`, includes data loss and effective frequency metrics in the output (default is FALSE).
#'
#' @return A `data.frame` with one row per eye-target combination, containing computed metrics:
#'   - `eye`, `target_id`: identifiers
#'   - `accuracy`, `accuracy_x`, `accuracy_y`: accuracy metrics (`accuracy_x`, `accuracy_y` only if `advanced` is `TRUE`)
#'   - `rms_s2s`, `rms_s2s_x`, `rms_s2s_y`: precision (RMS sample-to-sample) (`rms_s2s_x`, `rms_s2s_y` only if `advanced` is `TRUE`)
#'   - `std`, `std_x`, `std_y`: precision (standard deviation) (`std_x`, `std_y` only if `advanced` is `TRUE`)
#'   - `bcea`, `bcea_orientation`, `bcea_ax1`, `bcea_ax2`, `bcea_aspect_ratio`: precision (BCEA metrics) (`bcea_orientation`, `bcea_ax1`, `bcea_ax2`, `bcea_aspect_ratio` only if `advanced` is `TRUE`)
#'   - `data_loss`, `effective_frequency`: optional metrics if `include_data_loss = TRUE`
#'
#' @details
#' This function uses the following methods in the `DataQuality` class to compute the returned results:
#' `accuracy()`, `precision_RMS_S2S()`, `precision_STD()`, `precision_BCEA()`, `data_loss_from_invalid()`, and `effective_frequency()`.
#'
#' @examples
#' \dontrun{
#' # NB: this example requires a gaze data table to run. See the complete example at
#' # https://github.com/dcnieho/ETDQualitizer/blob/master/example/R.R for how to prepare
#' # the input data for this function
#' dq <- compute_data_quality_from_validation(gaze_data, unit = "pixels", screen = my_screen_config)
#' }
#'
#' @export
compute_data_quality_from_validation <- function(gaze, unit, screen = NULL, advanced = FALSE, include_data_loss = FALSE) {
  stopifnot(is.data.frame(gaze))
  stopifnot(unit %in% c("pixels", "degrees"))

  # Get all targets
  targets <- unique(gaze$target_id)
  targets <- targets[targets != -1]
  target_locations <- matrix(NA, nrow = length(targets), ncol = 2)

  for (t in seq_along(targets)) {
    qTarget <- gaze$target_id == targets[t]
    iTarget <- which(qTarget)[1]
    target_locations[t, ] <- c(gaze$tar_x[iTarget], gaze$tar_y[iTarget])
  }

  # Convert target locations to degrees if needed
  if (unit == "pixels") {
    if (is.null(screen)) {
      stop('If unit is "pixels", a screen configuration must be supplied')
    }
    deg_coords <- screen$pix_to_deg(target_locations[, 1], target_locations[, 2])
    target_locations[, 1] <- deg_coords$azi
    target_locations[, 2] <- deg_coords$ele
  } else if (unit != "degrees") {
    stop('unit should be "pixels" or "degrees"')
  }

  # Determine which eyes are present
  eyes <- c("left", "right")
  have_eye <- sapply(eyes, function(e) paste0(e, "_x") %in% names(gaze))
  eyes <- eyes[have_eye]

  # Prepare output table
  vars <- c("eye", "target_id", "accuracy", "accuracy_x", "accuracy_y",
            "rms_s2s", "rms_s2s_x", "rms_s2s_y",
            "std", "std_x", "std_y",
            "bcea", "bcea_orientation", "bcea_ax1", "bcea_ax2", "bcea_aspect_ratio")
  if (include_data_loss) {
    vars <- c(vars, "data_loss", "effective_frequency")
  }

  dq <- data.frame(matrix(NA, nrow = length(targets) * length(eyes), ncol = length(vars)))
  names(dq) <- vars
  dq$eye <- ""
  dq$target_id <- NA

  # Compute metrics per eye and target
  for (e in seq_along(eyes)) {
    eye <- eyes[e]
    x_col <- paste0(eye, "_x")
    y_col <- paste0(eye, "_y")

    for (t in seq_along(targets)) {
      oi <- (e - 1) * length(targets) + t
      t_id <- targets[t]
      is_target <- gaze$target_id == t_id

      dq_calc <- DataQuality$new(
        gaze_x = gaze[[x_col]][is_target],
        gaze_y = gaze[[y_col]][is_target],
        timestamps = gaze$timestamp[is_target] / 1000,
        unit = unit,
        screen = screen
      )

      dq$eye[oi] <- eye
      dq$target_id[oi] <- t_id
      dq[oi, c("accuracy", "accuracy_x", "accuracy_y")] <- dq_calc$accuracy(target_locations[t, 1], target_locations[t, 2])
      dq[oi, c("rms_s2s", "rms_s2s_x", "rms_s2s_y")] <- dq_calc$precision_RMS_S2S()
      dq[oi, c("std", "std_x", "std_y")] <- dq_calc$precision_STD()
      dq[oi, c("bcea", "bcea_orientation", "bcea_ax1", "bcea_ax2", "bcea_aspect_ratio")] <- dq_calc$precision_BCEA()

      if (include_data_loss) {
        dq$data_loss[oi] <- dq_calc$data_loss_from_invalid()
        dq$effective_frequency[oi] <- dq_calc$effective_frequency()
      }
    }
  }

  # Drop advanced metrics if not requested
  if (!advanced) {
    keep_vars <- c("eye", "target_id", "accuracy", "rms_s2s", "std", "bcea")
    if (include_data_loss) {
      keep_vars <- c(keep_vars, "data_loss", "effective_frequency")
    }
    dq <- dq[, keep_vars, drop = FALSE]
  }

  dq
}


#' Summarize and Report Data Quality Metrics
#'
#' This function summarizes data quality metrics from a validation procedure by computing averages per participant and generating descriptive statistics across participants.
#' It also returns a formatted textual summary suitable for reporting.
#'
#' @param dq_table A `data.frame` containing data quality metrics. Must include columns `file`, `eye`, `target_id`, and relevant numeric metrics such as `accuracy`, `rms_s2s`, and `std`.
#'  This would generally be created by concatenating the output of the compute_data_quality_from_validation() for multiple files.
#'
#' @return A named list with two elements:
#' \describe{
#'   \item{txt}{A character string summarizing key metrics (accuracy, RMS-S2S precision, STD precision).}
#'   \item{measures}{A list containing:
#'     \itemize{
#'       \item{\code{all}: A data frame with per-participant averages (grouped by `file`).}
#'       \item{\code{mean}, \code{std}, \code{min}, \code{max}: Named numeric vectors with summary statistics across participants.}
#'     }
#'   }
#' }
#'
#' @details
#' The summary text excludes BCEA and data loss metrics. BCEA is considered a niche metric and data loss is best reported across the full dataset rather than just the validation subset.
#'
#' @examples
#' \dontrun{
#' # NB: this example requires a gaze data table to run. See the complete example at
#' # https://github.com/dcnieho/ETDQualitizer/blob/master/example/R.R for how to prepare
#' # the input data for this function
#' result <- report_data_quality_table(dq_table)
#' cat(result$txt)
#' head(result$measures$all)
#' }
#'
#' @export
#' @importFrom stats aggregate
report_data_quality_table <- function(dq_table) {
  stopifnot(is.data.frame(dq_table))

  measures <- list()

  # Average over targets and eyes, grouped by file
  grouped <- aggregate(. ~ file, data = dq_table[, !(names(dq_table) %in% c("eye", "target_id"))], FUN = mean, na.rm = TRUE)
  measures$all <- grouped

  # Summary statistics
  numeric_cols <- names(grouped)[-1]  # exclude 'file'
  measures$mean <- sapply(grouped[numeric_cols], mean, na.rm = TRUE)
  measures$std  <- sapply(grouped[numeric_cols], sd, na.rm = TRUE)
  measures$min  <- sapply(grouped[numeric_cols], min, na.rm = TRUE)
  measures$max  <- sapply(grouped[numeric_cols], max, na.rm = TRUE)

  # Text summary (excluding BCEA and data loss)
  n_target <- length(unique(dq_table$target_id))
  n_subj   <- nrow(measures$all)
  version  <- ETDQ_version()  # Assumes this function exists

  txt <- sprintf(
    "For %d participants, the average inaccuracy in the data determined from a %d-point validation procedure using ETDQualitizer v%s (Niehorster et al., in prep) was %.2f\u00b0 (SD=%.2f\u00b0, range=%.2f\u00b0--%.2f\u00b0). Average RMS-S2S precision was %.3f\u00b0 (SD=%.3f\u00b0, range=%.3f\u00b0--%.3f\u00b0) and STD precision %.3f\u00b0 (SD=%.3f\u00b0, range=%.3f\u00b0--%.3f\u00b0).",
    n_subj, n_target, version,
    measures$mean["accuracy"], measures$std["accuracy"], measures$min["accuracy"], measures$max["accuracy"],
    measures$mean["rms_s2s"], measures$std["rms_s2s"], measures$min["rms_s2s"], measures$max["rms_s2s"],
    measures$mean["std"], measures$std["std"], measures$min["std"], measures$max["std"]
  )

  list(txt = txt, measures = measures)
}