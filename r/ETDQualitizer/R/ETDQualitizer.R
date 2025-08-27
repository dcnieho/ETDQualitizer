library(R6)
library(dplyr)
library(tidyr)
library(stats)

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

vector_to_Fick <- function(x, y, z) {
  azi <- 180/pi * atan2(x, z)
  ele <- 180/pi * atan2(y, sqrt(x^2 + z^2))
  list(azi = azi, ele = ele)
}


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

rms_s2s <- function(azi, ele, central_tendency_fun = mean) {
  a_diff <- diff(azi)^2
  e_diff <- diff(ele)^2
  rms_a  <- sqrt(central_tendency_fun(a_diff, na.rm = TRUE))
  rms_e  <- sqrt(central_tendency_fun(e_diff, na.rm = TRUE))
  rms    <- sqrt(central_tendency_fun(a_diff + e_diff, na.rm = TRUE))
  list(rms = rms, rms_a = rms_a, rms_e = rms_e)
}

std <- function(azi, ele) {
  std_a <- pop_sd(azi)
  std_e <- pop_sd(ele)
  std   <- sqrt(std_a^2 + std_e^2)
  list(std = std, std_a = std_a, std_e = std_e)
}
pop_sd <- function(x) {
  x <- x[!is.na(x)]
  sqrt(mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))
}

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

data_loss <- function(x, y) {
  missing <- is.na(x) | is.na(y)
  sum(missing)/length(missing)*100
}

data_loss_from_expected <- function(x, y, duration, frequency) {
  N_valid <- sum(!is.na(x) & !is.na(y))
  (1 - N_valid/(duration* frequency))*100
}

effective_frequency <- function(x, y, duration) {
  N_valid <- sum(!is.na(x) & !is.na(y))
  N_valid/duration
}


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


ScreenConfiguration <- R6Class("ScreenConfiguration",
  public = list(
    screen_size_x_mm = NULL,
    screen_size_y_mm = NULL,
    screen_res_x_pix = NULL,
    screen_res_y_pix = NULL,
    viewing_distance_mm = NULL,

    initialize = function(screen_size_x_mm, screen_size_y_mm, screen_res_x_pix, screen_res_y_pix, viewing_distance_mm) {
      self$screen_size_x_mm <- screen_size_x_mm
      self$screen_size_y_mm <- screen_size_y_mm
      self$screen_res_x_pix <- screen_res_x_pix
      self$screen_res_y_pix <- screen_res_y_pix
      self$viewing_distance_mm <- viewing_distance_mm
    },

    pix_to_mm = function(x, y) {
      x_mm <- x / self$screen_res_x_pix * self$screen_size_x_mm
      y_mm <- y / self$screen_res_y_pix * self$screen_size_y_mm
      list(x = x_mm, y = y_mm)
    },

    pix_to_deg = function(x, y) {
      mm <- self$pix_to_mm(x, y)
      self$mm_to_deg(mm$x, mm$y)
    },

    mm_to_deg = function(x, y) {
      azi <- atan2(x, self$viewing_distance_mm)
      ele <- atan2(y, sqrt(self$viewing_distance_mm^2 + x^2))
      list(azi = azi * 180 / pi, ele = ele * 180 / pi)
    },

    mm_to_pix = function(x, y) {
      x_pix <- x / self$screen_size_x_mm * self$screen_res_x_pix
      y_pix <- y / self$screen_size_y_mm * self$screen_res_y_pix
      list(x = x_pix, y = y_pix)
    },

    deg_to_mm = function(x, y) {
      cart <- Fick_to_vector(x, y)
      x_mm <- cart$x / cart$z * self$viewing_distance_mm
      y_mm <- cart$y / cart$z * self$viewing_distance_mm
      list(x = x_mm, y = y_mm)
    },

    deg_to_pix = function(x, y) {
      mm <- self$deg_to_mm(x, y)
      self$mm_to_pix(mm$x, mm$y)
    },

    screen_extents = function() {
      w <- self$mm_to_deg(self$screen_size_x_mm / 2, 0.)
      h <- self$mm_to_deg(0., self$screen_size_y_mm / 2)
      list(width = w$azi * 2, height = h$ele * 2)
    }
  )
)

DataQuality <- R6::R6Class("DataQuality",
  public = list(
    timestamps = NULL,
    azi = NULL,
    ele = NULL,

    # N.B: for this module it is assumed that any missing data are not coded
    # with some special value such as (0,0) or (-xres,-yres) but as NA. Missing
    # data should also not be removed, or the RMS calculation would be incorrect.
    #
    # timestamps should be in seconds.
    #
    # all angular positions are expected to be expressed in Fick angles.
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

    accuracy = function(target_x_deg, target_y_deg, central_tendency_fun = mean) {
      accuracy(self$azi, self$ele, target_x_deg, target_y_deg, central_tendency_fun)
    },

    precision_RMS_S2S = function(central_tendency_fun = mean) {
      rms_s2s(self$azi, self$ele, central_tendency_fun)
    },

    precision_STD = function() {
      std(self$azi, self$ele)
    },

    precision_BCEA = function(P = 0.68) {
      bcea(self$azi, self$ele, P)
    },

    data_loss = function() {
      data_loss(self$azi, self$ele)
    },

    data_loss_from_expected = function(frequency) {
      data_loss_from_expected(self$azi, self$ele, self$get_duration(), frequency)
    },

    effective_frequency = function() {
      effective_frequency(self$azi, self$ele, self$get_duration())
    },

    get_duration = function() {
      # to get duration right, we need to include duration of last sample
      isi <- median(diff(self$timestamps), na.rm = TRUE)
      self$timestamps[length(self$timestamps)] - self$timestamps[1] + isi
    },

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
