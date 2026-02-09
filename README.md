[![PyPI Latest Release](https://img.shields.io/pypi/v/ETDQualitizer.svg)](https://pypi.org/project/ETDQualitizer/)
[![Python Downloads](https://static.pepy.tech/badge/ETDQualitizer)](https://pepy.tech/project/ETDQualitizer)
[![File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/181899-etdqualitizer)
[![CRAN](https://www.r-pkg.org/badges/version/ETDQualitizer?color=green)](https://cran.r-project.org/package=ETDQualitizer)
[![](http://cranlogs.r-pkg.org/badges/grand-total/ETDQualitizer?color=green)](https://cran.r-project.org/package=ETDQualitizer)

# ETDQualitizer v0.10.0
ETDQualitizer is a toolbox for automated eye tracking data quality determination for screen-based eye trackers. This repository consists of two parts:
1. The Validator module: A [Python module for collecting validation data](#collecting-validation-data) in a format that is suitable for automated determination and reporting of data quality for screen-based eye trackers. This module can both be run as a standalone Python script using [PsychoPy](https://www.psychopy.org/) or an integratable module that can be added to an exisiting PsychoPy experiment using a single line of code.
2. ETDQualitizer: a [tool for determining and reporting eye-tracking data quality](#determining-data-quality) using recordings made by the Validator script, or any data formatted according to the same format. ETDQualitizer is available for MATLAB, Python and R, and can furthermore be run [as a webpage](https://dcnieho.github.io/ETDQualitizer/).

Note that for determining data quality for recordings made with wearable eye trackers, the [glassesValidator](https://github.com/dcnieho/glassesValidator) tool is available (which is also integrated in [gazeMapper](https://github.com/dcnieho/gazeMapper)).

Please cite:
Niehorster, D.C., Nyström, M., Hessels, R.S., Benjamins, J.S., Andersson, R. & Hooge, I.T.C. (in prep). The fundamentals of eye tracking part 7: Determining data quality

For questions, bug reports or to check for updates, please visit
www.github.com/dcnieho/ETDQualitizer.

# Content
Below, the following topics are discussed:
1. [The Validator module, its configuration, and its output data format](#collecting-validation-data).
2. [The ETDQualitizer tool for determining and reporting eye-tracking data quality](#determining-data-quality).
3. [A complete walkthrough using the Validator script and ETDQualitizer](#walkthrough).

# Collecting validation data
The procedure for collecting validation data is available from [the `/python/ETDQualitizer/validator` subfolder](/python/ETDQualitizer/validator). In this manual, we will refer to this script as the validator.

## How to acquire
The validator is not available as an installable package. The user is recommended to download [the `/python/ETDQualitizer/validator` subfolder](/python/ETDQualitizer/validator) from github and add it to their experiment.

## Documentation
Several files are contained in [the `/python/ETDQualitizer/validator` subfolder](/python/ETDQualitizer/validator) that together comprise the validator. The validator can be run in two ways, either as a standalone script, or integrated into an existing PsychoPy script. By default, support is included for EyeLink, SMI and Tobii eye trackers. To run as a standalone script, simply run the `run_for_eyelink.py`, `run_for_smi.py`, or `run_for_tobii.py` scripts. To integrate the validator into an existing script, simply call the `run_validation()` function inside the `validator.py` module in your existing PsychoPy script at appropriate time(s).

## Configuration files
The procedure (e.g., position of the validation targets, and duration for which they are shown) is configured using the following two configuration files.

### `targetPositions.csv`
The `targetPositions.csv` file contains the positions (and associated ID) for the validation targets. It contains the following columns:
|name|description|
| --- | --- |
| `ID` | Unique identifier of the validation target (integer). |
| `x` | Horizontal position of the center of the validation target. |
| `y` | Vertical position of the center of the validation target. |
| `color` | Optional column allowing to set the color of the target place holder for this target. |

Validation target coordinates are in PsychoPy's `'pix'` coordinate system. That is, they are in pixels with respect to the center of the screen, with the positive y direction being upward.

### `setup.json`
The `setup.json` file contains information about the screen and configuration for the validation procedure. The below table describes the options. It is a nested json file; nesting levels are indicated with `.` in the table.
|name|description|
| --- | --- |
| `screen.which_monitor` | Which monitor the validation screen will be shown on (only used for standalone operation) |
| `screen.width` | The width of the monitor (its display area to be exact) in cm |
| `screen.viewing_distance` | The distance from which the monitor is viewed |
| `screen.refresh_rate` | The refresh rate of the monitor |
| `screen.resolution` | The pixel resolution of the monitor |
| `screen.background_color` | The background color to use for the validation display |
|||
| `instruction_text.height` | The size of lines of the instruction text |
| `instruction_text.color` | The color of the instruction text |
|||
| `validation.n_repetitions` | The number of times each validation target should be shown during a single validation session |
| `validation.targets.file` | The name of the file containing the validation targets, probably [`targetPositions.csv`](#targetpositionscsv) |
| `validation.targets.duration` | The duration of gaze data to capture for each target |
| `validation.targets.randomize_order` | Boolean indicating whether the target presentation order should be shuffled or as in the file |
| `validation.targets.placeholder.color` | Color of the placeholder to be shown at all fixation target locations (`null`, Python's `None` to disable) |
| `validation.targets.placeholder.diameter` | Diameter of the placeholder circle |
| `validation.targets.cue.color` | Color of the cue to be shown as the fixation marker is traveling to a next fixation target location (`null`, Python's `None` to disable) |
| `validation.targets.cue.diameter` | Diameter of the cue circle |
| `validation.targets.move.duration` | Duration it would take for the fixation target to move from one end of the screen to the other, horizontally. This duration is scaled by the actual distance between the current and the next fixation target location |
| `validation.targets.move.diameter` | Diameter of the fixation target as it is moving |
| `validation.targets.move.min_duration` | Minimum duration of move from the current to next fixation target location |
| `validation.targets.move.move_with_acceleration` | Boolean indicating whether the target should accelerate and decelerate as it moves between the current and next fixation target locations, or move at a constant speed |
| `validation.targets.shrink.duration` | Duration of interval during which the target shrinks after it has arrived at the locaiton of the upcoming fixation target |
| `validation.targets.look.diameter_min` | Diameter at the start of the shrinking interval |
| `validation.targets.look.diameter_max` | Diameter at the end of the shrinking interval|
| `validation.targets.look.outer_color` | Color of the large and small circles in the fixation target |
| `validation.targets.look.inner_color` | Color of the cross in the fixation target |

All sizes are in pixels, and all durations in seconds. Ensure that the screen information is correctly filled out or the timing or positioning of validation targets may be incorrect.

## Output data format
The data files output by the validator contain raw gaze data (indvidual samples) with the following columns:
|name|description|
| --- | --- |
| `timestamp` | Timestamp of the sample, in ms |
| `left_x` | Horizontal position of the gaze point of the left eye, in pixels |
| `left_y` | Vertical position of the gaze point of the left eye, in pixels |
| `right_x` | Horizontal position of the gaze point of the right eye, in pixels |
| `right_y` | Vertical position of the gaze point of the right eye, in pixels |
| `target_id` | ID of the target at which the participant looked. -1 for samples that should not be included in the data quality calculations, an integer ID (see description of [`targetPositions.csv`](#targetpositionscsv)) for intervals to be included |
| `tar_x` | Horizontal screen position of the center of the target, in pixels |
| `tar_y` | Vertical screen position of the center of the target, in pixels |

Note that positions in this file _do not_ use PsychoPy's `'pix'` coordinate system where the positive y direction is upward. Like in most other conventions, the positive y direction is downward. Positions of both gaze and targets are expressed with respect to the center of the screen. Note that it is possible a file contains data from only the left or only the right eye, in which case the other columns should be missing.

## Adding your own eye tracker
The validator provides out of the box support for EyeLink, SMI and Tobii eye trackers. To add support for other eye trackers, users should add a class wrapping communication with their eye tracker in [`validator.py`](/python/ETDQualitizer/validator/validator.py). This class should derive from the `EyeTrackerBase` class in that file, and implement the methods `start_recording`, `send_message`, `stop_recording`, and `save_data`. See the `EyeLinkTracker`, `SMIEyeTracker`, and `TobiiEyeTracker` classes for examples.

# Determining Data Quality
Example scripts for using the MATLAB, Python and R implementations, as well as example data are available from [the `/example` subfolder](/example).

## How to acquire
While the [MATLAB](/matlab), [Python](/python/ETDQualitizer) and [R](/r/ETDQualitizer) implementations can be directly downloaded from their respective subfolders in this repository, the recommended installation method is to use the respective platform's package manager:
### MATLAB
Search for and install the ETDQualitizer package from the MATLAB Add-Ons manager.
### Python
Execute `pip install ETDQualitizer`
### R
Execute `install.packages("ETDQualitizer")`

## Documentation
The documentation is integrated with the code and accessible using the `help` function in MATLAB, Python and R. As such, it will not be repeated here. Here we provide a high-level overview of the ETDQualitizer package's design.
- Basic functions: the package provides the following functions used for determining data quality of a segment of eye tracking data: `accuracy`, `rms_s2s`, `std`, `bcea`, `data_loss_from_invalid`, `data_loss_from_expected`, `effective_frequency`. Furthermore, the function `precision_using_moving_window` is provided for determining RMS-S2S, STD or BCEA using a moving window approach. The algorithm implemented in each of these functions is described in Niehorster et al. (in prep).
- Coordinate transformation class and functions: the package provides the `ScreenConfiguration` class, which can be used to convert data represented as pixels on a monitor to physical distances in the world (e.g., mm) and to angular gaze directions (using Fick angles), and vice-versa. Further available are the functions `Fick_to_vector` and `vector_to_Fick` that turn angular gaze directions into gaze vectors, and vice versa.
- Data quality computation class: the package provides the `DataQuality` class, which is a wrapper around the basic functions above. It is constructed using a segment of eye tracking data (optionally converted to gaze angles if they are not already expressed as such). The various data quality measures can then be directly computed on the segment by using the class methods.
- Wrapper functions: finally, the package provides two higher-level wrapper functions, `compute_data_quality_from_validation` and `report_data_quality_table`. `compute_data_quality_from_validation` can be used to take a segment of eye tracking formatted according to [the format produced by the Validator module](#output-data-format) and calculate the various data quality metrics for each presented target in this segment. The `report_data_quality_table` then takes one or multiple tables output by `compute_data_quality_from_validation` and provides an overall summary of data quality across targets and (if multiple tables are provided) recordings/participants, along with a textual report that can be directly pasted into a paper.

# Walkthrough
This walkthrough will help you become familiar with how to use the Validator script and the ETDQualitizer package to determine the data quality of an eye tracking recording. It uses the [standalone Validator script for a Tobii] as an example, along with the [ETDQualitizer tool packaged as a webpage](https://dcnieho.github.io/ETDQualitizer/). Users who prefer to use ETDQualitizer from their own MATLAB, Python or R script are referred to the [respective example scripts](/example). The walkthrough is organized in two parts, 1) acquiring gaze data that is suitable for determining data quality and 2) determining eye tracking data quality from this gaze data using ETDQualitizer. If one only wishes to learn how to determine data quality using ETDQualitizer, one can skip to this section and use the provided [example data](/example/data).

## Collecting gaze data for data quality assessment using the Validator script

## Determining data quality for screen-based eye tracker recordings using ETDQualitizer

# Citation
If you use this tool or any of the code in this repository, please cite:<br>
Niehorster, D.C., Nyström, M., Hessels, R.S., Benjamins, J.S., Andersson, R. & Hooge, I.T.C. (in prep). The fundamentals of eye tracking part 7: Determining data quality

## BibTeX
```latex
@article{niehorsterFundamentals7,
    Author = {Niehorster, Diederick C. and
              Nystr{\"o}m, Marcus and
              Hessels, Roy S. and
              Benjamins, Jeroen S. and
              Andersson, Richard and
              Hooge, Ignace T. C.},
    Journal = {},
    Number = {},
    Title = {The fundamentals of eye tracking part 7: Determining data quality},
    Year = {in prep},
    doi = {}
}
