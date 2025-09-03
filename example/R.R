library(readr)
library(dplyr)
library(stringr)
library(ETDQualitizer)

get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    normalizePath(dirname(sub("--file=", "", file_arg)))
  } else if (!is.null(sys.frame(1)$ofile)) {
    normalizePath(dirname(sys.frame(1)$ofile))
  } else {
    getwd()
  }
}

# Set up screen configuration
screen <- ScreenConfiguration$new(
  screen_size_x_mm = 528.0,
  screen_size_y_mm = 296.9997253417969,
  screen_res_x_pix = 1920,
  screen_res_y_pix = 1080,
  viewing_distance_mm = 650
)

# Locate data folder
script_path <- get_script_path()  # or use getwd() if running interactively
if (basename(script_path)=="ETDQualitizer") {
    script_path <- file.path(script_path, 'example')
}
data_path <- file.path(script_path, "data")
tsv_files <- list.files(data_path, pattern = "\\.tsv$", full.names = TRUE)

all_dq <- list()

for (file in tsv_files) {
  cat("----------\n", basename(file), "\n")
  gaze <- read_tsv(file, show_col_types = FALSE)

  # automatically compute data quality measures per target
  dq <- compute_data_quality_from_validation(gaze, "pixels", screen, advanced = FALSE, include_data_loss = TRUE)
  dq$file <- basename(file)
  all_dq[[length(all_dq) + 1]] <- dq

  # manually perform some further data quality computations on data from the whole validation instead of per target
  eyes <- c("left", "right")
  for (eye in eyes) {
    x_col <- paste0(eye, "_x")
    y_col <- paste0(eye, "_y")
    if (!(x_col %in% names(gaze))) next

    dq_calc <- DataQuality$new(
      gaze_x = gaze[[x_col]],
      gaze_y = gaze[[y_col]],
      timestamps = gaze$timestamp / 1000,
      unit = "pixels",
      screen = screen
    )

    # determine sampling frequency from the filename (assumes our test files which end in '<xxx>Hz')
    parts <- str_split(basename(file), "_", simplify = TRUE)
    fs <- as.numeric(str_remove(parts[ncol(parts)], "Hz\\.tsv"))
    window_len <- round(0.2 * fs)

    # and RMS-S2S calculated in two ways over the whole datafile
    cat(sprintf("RMS-S2S using median (%s eye): %.4f deg\n", eye,
                dq_calc$precision_RMS_S2S(median)[1]))
    cat(sprintf("RMS-S2S using moving window (%s eye): %.4f deg\n", eye,
                dq_calc$precision_using_moving_window(window_len, "RMS-S2S")))

    # data loss and effective frequency
    cat(sprintf("Data loss from invalid samples (%s eye): %.1f%%\n", eye, dq_calc$data_loss_from_invalid()))
    cat(sprintf("Data loss from expected #samples (%s eye): %.1f%%\n", eye, dq_calc$data_loss_from_expected(fs)))
    cat(sprintf("Effective frequency (%s eye): %.1f Hz\n", eye, dq_calc$effective_frequency()))
  }
}

# Combine all results
all_dq_df <- bind_rows(all_dq)

# Generate report
report <- report_data_quality_table(all_dq_df)
cat("--------\nAll:\n")
cat(report$txt, "\n")
print(report$measures$all)
