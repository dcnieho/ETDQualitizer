library(R6)
library(dplyr)
library(tidyr)
library(stats)

Fick_to_cartesian <- function(azi, ele, r = 1) {
  azi <- azi * pi / 180
  ele <- ele * pi / 180
  r_cos_ele <- r * cos(ele)
  x <- r_cos_ele * sin(azi)
  y <- r * sin(ele)
  z <- r_cos_ele * cos(azi)
  return(list(x = x, y = y, z = z))
}


accuracy <- function(x, y, target_x_deg, target_y_deg, central_tendency_fun = mean) {
  g <- Fick_to_cartesian(x, y)
  t <- Fick_to_cartesian(target_x_deg, target_y_deg)
  offsets <- acos(g$x * t$x + g$y * t$y + g$z * t$z)
  direction <- atan2(g$y / g$z - t$y / t$z, g$x / g$z - t$x / t$z)
  offsets_2D <- offsets * cbind(cos(direction), sin(direction)) * 180 / pi
  offset_x <- central_tendency_fun(offsets_2D[, 1], na.rm = TRUE)
  offset_y <- central_tendency_fun(offsets_2D[, 2], na.rm = TRUE)
  return(c(sqrt(offset_x^2 + offset_y^2), offset_x, offset_y))
}

rms_s2s <- function(x, y, central_tendency_fun = mean) {
  x_diff <- diff(x)^2
  y_diff <- diff(y)^2
  return(c(
    sqrt(central_tendency_fun(x_diff + y_diff, na.rm = TRUE)),
    sqrt(central_tendency_fun(x_diff, na.rm = TRUE)),
    sqrt(central_tendency_fun(y_diff, na.rm = TRUE))
  ))
}

std <- function(x, y) {
  std_x <- sd(x, na.rm = TRUE)
  std_y <- sd(y, na.rm = TRUE)
  return(c(sqrt(std_x^2 + std_y^2), std_x, std_y))
}

bcea <- function(x, y, P = 0.68) {
  k <- log(1 / (1 - P))
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  std_x <- sd(x)
  std_y <- sd(y)
  rho <- cor(x, y)
  area <- 2 * k * pi * std_x * std_y * sqrt(1 - rho^2)
  cov_matrix <- cov(cbind(x, y))
  eig <- eigen(cov_matrix)
  i <- which.max(eig$values)
  orientation <- atan2(eig$vectors[2, i], eig$vectors[1, i]) * 180 / pi
  ax1 <- sqrt(k * eig$values[i])
  ax2 <- sqrt(k * eig$values[-i])
  aspect_ratio <- max(ax1, ax2) / min(ax1, ax2)
  return(c(area, orientation, ax1, ax2, aspect_ratio))
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
                                   return(c(x_mm, y_mm))
                                 },

                                 pix_to_deg = function(x, y) {
                                   mm <- self$pix_to_mm(x, y)
                                   return(self$mm_to_deg(mm[1], mm[2]))
                                 },

                                 mm_to_deg = function(x, y) {
                                   azi <- atan2(x, self$viewing_distance_mm)
                                   ele <- atan2(y, sqrt(self$viewing_distance_mm^2 + x^2))
                                   return(c(azi * 180 / pi, ele * 180 / pi))
                                 },

                                 mm_to_pix = function(x, y) {
                                   x_pix <- x / self$screen_size_x_mm * self$screen_res_x_pix
                                   y_pix <- y / self$screen_size_y_mm * self$screen_res_y_pix
                                   return(c(x_pix, y_pix))
                                 },

                                 deg_to_mm = function(x, y) {
                                   cart <- Fick_to_cartesian(x, y)
                                   x_mm <- cart$x / cart$z * self$viewing_distance_mm
                                   y_mm <- cart$y / cart$z * self$viewing_distance_mm
                                   return(c(x_mm, y_mm))
                                 },

                                 deg_to_pix = function(x, y) {
                                   mm <- self$deg_to_mm(x, y)
                                   return(self$mm_to_pix(mm[1], mm[2]))
                                 },

                                 screen_extents = function() {
                                   deg <- self$mm_to_deg(self$screen_size_x_mm / 2, self$screen_size_y_mm / 2)
                                   return(c(deg[1] * 2, deg[2] * 2))
                                 }
                               )
)

DataQuality <- R6Class("DataQuality",
                       public = list(
                         x = NULL,
                         y = NULL,
                         timestamps = NULL,

                         initialize = function(gaze_x, gaze_y, timestamps, unit, screen = NULL) {
                           self$timestamps <- timestamps
                           if (unit == "pixels") {
                             if (is.null(screen)) stop("If unit is 'pixels', a screen configuration must be supplied")
                             coords <- screen$pix_to_deg(gaze_x, gaze_y)
                             gaze_x <- coords[1]
                             gaze_y <- coords[2]
                           } else if (unit != "degrees") {
                             stop("unit should be 'pixels' or 'degrees'")
                           }
                           self$x <- gaze_x
                           self$y <- gaze_y
                         },

                         accuracy = function(target_x_deg, target_y_deg, central_tendency_fun = mean) {
                           accuracy(self$x, self$y, target_x_deg, target_y_deg, central_tendency_fun)
                         },

                         precision_RMS_S2S = function(central_tendency_fun = mean) {
                           rms_s2s(self$x, self$y, central_tendency_fun)
                         },

                         precision_STD = function() {
                           std(self$x, self$y)
                         },

                         precision_BCEA = function(P = 0.68) {
                           bcea(self$x, self$y, P)
                         },

                         data_loss = function() {
                           data_loss(self$x, self$y)
                         },

                         data_loss_from_expected = function(frequency) {
                           data_loss_from_expected(self$x, self$y, self$get_duration(), frequency)
                         },

                         effective_frequency = function() {
                           effective_frequency(self$x, self$y, self$get_duration())
                         },

                         get_duration = function() {
                           isi <- median(diff(self$timestamps))
                           return(tail(self$timestamps, 1) - self$timestamps[1] + isi)
                         },

                         precision_using_moving_window = function(window_length, metric, aggregation_fun = median, ...) {
                           precision_using_moving_window(self$x, self$y, window_length, metric, aggregation_fun, ...)
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
