library(R6)
library(dplyr)
library(tidyr)
library(stats)


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
#' @param x Gaze azimuth in degrees.
#' @param y Gaze elevation in degrees.
#' @param target_x_deg Target azimuth in degrees.
#' @param target_y_deg Target elevation in degrees.
#' @param central_tendency_fun Function to compute central tendency (default: \code{mean}).
#'
#' @return A list with \code{offset}, \code{offset_x}, and \code{offset_y}, the total, horizontal and vertical offset of gaze from the target (in degrees).
#' @examples
#' accuracy(c(1, 2), c(1, 2), 0, 0)
#' @export
accuracy <- function(x, y, target_x_deg, target_y_deg, central_tendency_fun = mean) {
  # Get unit vectors from gaze directions
  g <- Fick_to_vector(x, y)
  t <- Fick_to_vector(target_x_deg, target_y_deg)
  # calculate angular offset for each sample using dot product
  dot_products <- g$x*t$x + g$y*t$y + g$z*t$z
  dot_products <- pmin(pmax(dot_products, -1), 1)  # Clamp to [-1,1]
  offsets      <- acos(dot_products)
  # calculate on-screen orientation so we can decompose offset into x and y
  direction    <- atan2(g$y/g$z-t$y/t$z, g$x/g$z-t$x/t$z)
  offsets_2D   <- offsets*cbind(cos(direction), sin(direction)) * 180/pi
  # calculate mean horizontal and vertical offsets
  offset_x     <- central_tendency_fun(offsets_2D[, 1], na.rm = TRUE)
  offset_y     <- central_tendency_fun(offsets_2D[, 2], na.rm = TRUE)
  # calculate offset of centroid
  offset       <- sqrt(offset_x^2 + offset_y^2)
  list(offset = offset, offset_x = offset_x, offset_y = offset_y)
}

#' RMS of Sample-to-Sample Differences
#'
#' Computes root mean square of differences between successive gaze samples.
#'
#' @param azi Azimuth values in degrees.
#' @param ele Elevation values in degrees.
#' @param central_tendency_fun Function to compute central tendency (default: \code{mean}).
#'
#' @return A list with \code{rms}, \code{rms_a}, and \code{rms_e}, the total RMS of sample-to-sample distances and that of the azimuthal and elevation components (all in degrees).
#' @examples
#' rms_s2s(c(1, 2, 3), c(1, 2, 3))
#' @export
rms_s2s <- function(azi, ele, central_tendency_fun = mean) {
  a_diff <- diff(azi)^2
  e_diff <- diff(ele)^2
  rms_a  <- sqrt(central_tendency_fun(a_diff, na.rm = TRUE))
  rms_e  <- sqrt(central_tendency_fun(e_diff, na.rm = TRUE))
  rms    <- sqrt(central_tendency_fun(a_diff + e_diff, na.rm = TRUE))
  list(rms = rms, rms_a = rms_a, rms_e = rms_e)
}

#' Standard Deviation of Gaze Samples
#'
#' Computes standard deviation of azimuth and elevation.
#'
#' @param azi Azimuth values in degrees.
#' @param ele Elevation values in degrees.
#'
#' @return A list with \code{std}, \code{std_a}, and \code{std_e}, the total STD of sample-to-sample distances and that of the azimuthal and elevation components (all in degrees).
#' @examples
#' std(c(1, 2, 3), c(1, 2, 3))
#' @export
std <- function(azi, ele) {
  std_a <- pop_sd(azi)
  std_e <- pop_sd(ele)
  std   <- sqrt(std_a^2 + std_e^2)
  list(std = std, std_a = std_a, std_e = std_e)
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
  sqrt(mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))
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
  cov_matrix <- cov(cbind(azi, ele))
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

#' Compute Data Loss Percentage
#'
#' Calculates percentage of missing gaze samples.
#'
#' @param x Azimuth values.
#' @param y Elevation values.
#'
#' @return Percentage of missing samples.
#' @examples
#' data_loss(c(1, NA, 3), c(1, 2, NA))
#' @export
data_loss <- function(x, y) {
  missing <- is.na(x) | is.na(y)
  sum(missing)/length(missing)*100
}

#' Compute Data Loss from Expected Sample Count
#'
#' Calculates data loss based on expected number of samples.
#'
#' @param x Azimuth values.
#' @param y Elevation values.
#' @param duration Duration in seconds.
#' @param frequency Sampling frequency in Hz.
#'
#' @return Percentage of data loss.
#' @examples
#' data_loss_from_expected(c(1, NA, 3), c(1, 2, NA), duration = 1, frequency = 3)
#' @export
data_loss_from_expected <- function(x, y, duration, frequency) {
  N_valid <- sum(!is.na(x) & !is.na(y))
  (1 - N_valid/(duration* frequency))*100
}

#' Compute Effective Sampling Frequency
#'
#' Calculates effective frequency based on valid samples.
#'
#' @param x Azimuth values.
#' @param y Elevation values.
#' @param duration Duration in seconds.
#'
#' @return Effective frequency in Hz.
#' @examples
#' effective_frequency(c(1, NA, 3), c(1, 2, NA), duration = 1)
#' @export
effective_frequency <- function(x, y, duration) {
  N_valid <- sum(!is.na(x) & !is.na(y))
  N_valid/duration
}


#' Precision Using Moving Window
#'
#' Computes gaze precision using a moving window and selected metric.
#'
#' @param x Azimuth values.
#' @param y Elevation values.
#' @param window_length Window size in samples.
#' @param metric Precision metric: \code{"RMS_S2S"}, \code{"STD"}, or \code{"BCEA"}.
#' @param aggregation_fun Function to aggregate precision values across the windows (default: \code{median}).
#' @param ... Additional arguments passed to metric function.
#'
#' @return Aggregated precision value.
#' @examples
#' precision_using_moving_window(rnorm(100), rnorm(100), 10, "STD")
#' @export
precision_using_moving_window <- function(x, y, window_length, metric, aggregation_fun = median, ...) {
  # Select the appropriate precision metric function
  fun <- switch(metric,
    "RMS_S2S" = rms_s2s,
    "STD"     = std,
    "BCEA"    = bcea,
    stop(sprintf('metric "%s" is not understood', metric))
  )

  # Get number of samples in data
  ns <- length(x)

  if (window_length < ns) {
    # If number of samples in data exceeds window size
    values <- rep(NA_real_, ns - window_length + 1)  # pre-allocate
    for (p in seq_len(ns - window_length + 1)) {
      result <- fun(x[p:(p + window_length - 1)], y[p:(p + window_length - 1)], ...)
      values[p] <- result[[1]]  # extract first element (e.g., rms, std, or area)
    }
    precision <- aggregation_fun(values, na.rm = TRUE)
  } else {
    # If too few samples in data
    precision <- NA_real_
  }

  precision
}


#' Screen Configuration Class
#'
#' Provides methods for converting between pixel, millimeter, and degree units.
#'
#' @examples
#' sc <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
#' sc$pix_to_deg(960, 540)
#' @export
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

    #' Initialize Screen Configuration
    #'
    #' Creates a new ScreenConfiguration object with screen and viewing distance parameters.
    #'
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

    #' Convert Pixels to Millimeters
    #'
    #' Converts pixel coordinates to millimeter coordinates on the screen.
    #'
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

    #' Convert Pixels to Degrees
    #'
    #' Converts pixel coordinates to an angular gaze direction in degrees.
    #'
    #' @param x Horizontal pixel coordinate.
    #' @param y Vertical pixel coordinate.
    #' @return A list with azimuth (\code{"azi"}) and elevation (\code{"ele"}) in degrees.
    #' @examples
    #' sc$pix_to_deg(960, 540)
    pix_to_deg = function(x, y) {
      mm <- self$pix_to_mm(x, y)
      self$mm_to_deg(mm$x, mm$y)
    },

    #' Convert Millimeters to Degrees
    #'
    #' Converts millimeter coordinates to an angular gaze direction in degrees.
    #'
    #' @param x Horizontal position in millimeters.
    #' @param y Vertical position in millimeters.
    #' @return A list with azimuth (\code{"azi"}) and elevation (\code{"ele"}) in degrees.
    #' @examples
    #' sc$mm_to_deg(100, 50)
    mm_to_deg = function(x, y) {
      azi <- atan2(x, self$viewing_distance_mm)
      ele <- atan2(y, sqrt(self$viewing_distance_mm^2 + x^2))
      list(azi = azi * 180 / pi, ele = ele * 180 / pi)
    },

    #' Convert Millimeters to Pixels
    #'
    #' Converts millimeter coordinates on the screen to pixel coordinates.
    #'
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

    #' Convert Degrees to Millimeters
    #'
    #' Converts an angular gaze direction in degrees to millimeter coordinates on the screen.
    #'
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

    #' Convert Degrees to Pixels
    #'
    #' Converts an angular gaze direction in degrees to pixel coordinates.
    #'
    #' @param x Azimuth in degrees (Fick angles).
    #' @param y Elevation in degrees (Fick angles).
    #' @return A list with x and y in pixels.
    #' @examples
    #' sc$deg_to_pix(2, 1)
    deg_to_pix = function(x, y) {
      mm <- self$deg_to_mm(x, y)
      self$mm_to_pix(mm$x, mm$y)
    },

    #' Get Screen Extents in Degrees
    #'
    #' Computes the horizontal and vertical extents of the screen (in degrees).
    #'
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

#' Data Quality Class
#'
#' Provides methods for assessing the quality of gaze data, including accuracy, precision, data loss, and effective sampling frequency.
#'
#' @examples
#' sc <- ScreenConfiguration$new(500, 300, 1920, 1080, 600)
#' dq <- DataQuality$new(gaze_x, gaze_y, timestamps, unit = "pixels", screen = sc)
#' dq$accuracy(0, 0)
#' dq$precision_RMS_S2S()
#' dq$data_loss()
#' @export
DataQuality <- R6Class("DataQuality",
  public = list(
    #' @field timestamps Vector of timestamps in seconds. Samples with missing data should not be removed, or the RMS calculation would be incorrect.
    #' @field azi Vector of azimuth angles in degrees (Fick coordinates). Missing data should be coded as NA, not using some special value such as (0,0) or (-xres,-yres).
    #' @field ele Vector of elevation angles in degrees (Fick coordinates). Missing data should be coded as NA, not using some special value such as (0,0) or (-xres,-yres).
    timestamps = NULL,
    azi = NULL,
    ele = NULL,

    #' Initialize Data Quality
    #'
    #' Creates a new DataQuality object from gaze data and timestamps.
    #'
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

    #' Compute Accuracy
    #'
    #' Calculates the accuracy of gaze data relative to a known target location.
    #'
    #' @param target_x_deg Target azimuth in degrees.
    #' @param target_y_deg Target elevation in degrees.
    #' @param central_tendency_fun Function to compute central tendency (e.g., \code{mean}, \code{median}).
    #' @return Accuracy in degrees.
    #' @examples
    #' dq$accuracy(0, 0)
    accuracy = function(target_x_deg, target_y_deg, central_tendency_fun = mean) {
      accuracy(self$azi, self$ele, target_x_deg, target_y_deg, central_tendency_fun)
    },

    #' Compute Precision (RMS S2S)
    #'
    #' Calculates precision using root mean square of sample-to-sample differences.
    #'
    #' @param central_tendency_fun Function to compute central tendency (e.g., \code{mean}, \code{median}).
    #' @return Precision in degrees.
    #' @examples
    #' dq$precision_RMS_S2S()
    precision_RMS_S2S = function(central_tendency_fun = mean) {
      rms_s2s(self$azi, self$ele, central_tendency_fun)
    },

    #' Compute Precision (Standard Deviation)
    #'
    #' Calculates precision using standard deviation of gaze positions.
    #'
    #' @return Standard deviation in degrees.
    #' @examples
    #' dq$precision_STD()
    precision_STD = function() {
      std(self$azi, self$ele)
    },

    #' Compute Precision (BCEA)
    #'
    #' Calculates the Bivariate Contour Ellipse Area (BCEA) and ellipse parameters for gaze precision.
    #'
    #' @param P Proportion of data to include in the ellipse (default is 0.68).
    #' @return BCEA in degrees-squared.
    #' @examples
    #' dq$precision_BCEA()
    precision_BCEA = function(P = 0.68) {
      bcea(self$azi, self$ele, P)
    },

    #' Compute Data Loss
    #'
    #' Calculates the proportion of missing data (coded as NA).
    #'
    #' @return Proportion of missing samples.
    #' @examples
    #' dq$data_loss()
    data_loss = function() {
      data_loss(self$azi, self$ele)
    },

    #' Compute Data Loss from Expected Number of Samples
    #'
    #' Estimates data loss based on expected number of samples given the duration and sampling frequency.
    #'
    #' @param frequency Expected sampling frequency in Hz.
    #' @return Proportion of missing samples.
    #' @examples
    #' dq$data_loss_from_expected(500)
    data_loss_from_expected = function(frequency) {
      data_loss_from_expected(self$azi, self$ele, self$get_duration(), frequency)
    },

    #' Compute Effective Sampling Frequency
    #'
    #' Calculates the effective sampling frequency based on timestamps.
    #'
    #' @return Effective frequency in Hz.
    #' @examples
    #' dq$effective_frequency()
    effective_frequency = function() {
      effective_frequency(self$azi, self$ele, self$get_duration())
    },

    #' Get Duration of Gaze Data
    #'
    #' Computes the total duration of the gaze recording, including the last sample.
    #'
    #' @return Duration in seconds.
    #' @examples
    #' dq$get_duration()
    get_duration = function() {
      # to get duration right, we need to include duration of last sample
      isi <- median(diff(self$timestamps), na.rm = TRUE)
      self$timestamps[length(self$timestamps)] - self$timestamps[1] + isi
    },

    #' Compute Precision Using Moving Window
    #'
    #' Calculates precision using a moving window approach.
    #'
    #' @param window_length Length of the moving window in number of samples.
    #' @param metric Precision metric to use (e.g., \code{"RMS"}, \code{"STD"}).
    #' @param aggregation_fun Function to aggregate windowed precision values (e.g., \code{median}).
    #' @param ... Additional arguments passed to the precision metric function.
    #' @return Precision value.
    #' @examples
    #' dq$precision_using_moving_window(0.2, "RMS")
    precision_using_moving_window = function(window_length, metric, aggregation_fun = median, ...) {
      precision_using_moving_window(self$azi, self$ele, window_length, metric, aggregation_fun, ...)
    }
  )
)

compute_data_quality_from_validation <- function(gaze, unit, screen = NULL, advanced = FALSE, include_data_loss = FALSE) {
  targets <- sort(unique(gaze$target_id[gaze$target_id != -1]))
  target_locations <- sapply(targets, function(t) {
    idx <- which(gaze$target_id == t)[1]
    c(gaze$tar_x[idx], gaze$tar_y[idx])
  })
  target_locations <- t(target_locations)

  if (unit == "pixels") {
    if (is.null(screen)) stop("If unit is 'pixels', a screen configuration must be supplied")
    deg_coords <- screen$pix_to_deg(target_locations[,1], target_locations[,2])
    target_locations[,1] <- deg_coords[1]
    target_locations[,2] <- deg_coords[2]
  } else if (unit != "degrees") {
    stop("unit should be 'pixels' or 'degrees'")
  }

  rows <- list()
  for (e in c("left", "right")) {
    if (!(paste0(e, "_x") %in% colnames(gaze))) next
    for (i in seq_along(targets)) {
      t_id <- targets[i]
      is_target <- gaze$target_id == t_id
      dq <- DataQuality$new(
        gaze_x = gaze[[paste0(e, "_x")]][is_target],
        gaze_y = gaze[[paste0(e, "_y")]][is_target],
        timestamps = gaze$timestamp[is_target] / 1000,
        unit = unit,
        screen = screen
      )
      row <- list(eye = e, target_id = t_id)
      row <- c(row, setNames(as.list(dq$accuracy(target_locations[i,1], target_locations[i,2])), c("offset", "offset_x", "offset_y")))
      row <- c(row, setNames(as.list(dq$precision_RMS_S2S()), c("rms_s2s", "rms_s2s_x", "rms_s2s_y")))
      row <- c(row, setNames(as.list(dq$precision_STD()), c("std", "std_x", "std_y")))
      row <- c(row, setNames(as.list(dq$precision_BCEA()), c("bcea", "bcea_orientation", "bcea_ax1", "bcea_ax2", "bcea_aspect_ratio")))
      if (include_data_loss) {
        row$data_loss <- dq$data_loss()
        row$effective_frequency <- dq$effective_frequency()
      }
      rows[[length(rows) + 1]] <- row
    }
  }

  dq_df <- bind_rows(rows)
  if (!advanced) {
    keep_cols <- c("eye", "target_id", "offset", "rms_s2s", "std", "bcea", "data_loss", "effective_frequency")
    dq_df <- dq_df[, intersect(names(dq_df), keep_cols)]
  }
  return(dq_df)
}
