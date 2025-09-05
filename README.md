[![PyPI Latest Release](https://img.shields.io/pypi/v/ETDQualitizer.svg)](https://pypi.org/project/ETDQualitizer/)
[![Python Downloads](https://static.pepy.tech/badge/ETDQualitizer)](https://pepy.tech/project/ETDQualitizer)
[![File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/181328-etdqualitizer)
[![CRAN](https://www.r-pkg.org/badges/version/ETDQualitizer?color=green)](https://cran.r-project.org/package=ETDQualitizer)
[![](http://cranlogs.r-pkg.org/badges/grand-total/ETDQualitizer?color=green)](https://cran.r-project.org/package=ETDQualitizer)

# ETDQualitizer v0.9.0
ETDQualitizer is a toolbox for automated eye tracking data quality determination for screen-based eye trackers. This repository offers two tools:
1. A [tool for collecting validation data](#collecting-validation-data) in a format that is suitable for automated determination and reporting of data quality for screen-based eye trackers. This tool is available in the form of a standalone Python script using [PsychoPy](https://www.psychopy.org/) or an integratable module that can be added to an exisiting PsychoPy experiment using a single line of code.
2. A [tool for determining and reporting eye-tracking data quality](#determining-data-quality) using recordings made by the first tool, or any data formatted according to the same format. These tools are available for MATLAB, Python and R. This tool is furthermore available [as a webpage](https://dcnieho.github.io/ETDQualitizer/).

Note that for wearable eye trackers, the [glassesValidator](https://github.com/dcnieho/glassesValidator) and [gazeMapper](https://github.com/dcnieho/gazeMapper) tools are available.

Please cite:
Niehorster, D.C., Nyström, M., Hessels, R.S., Benjamins, J.S., Andersson, R. & Hooge, I.T.C. (in prep). The fundamentals of eye tracking part 7: Data quality

For questions, bug reports or to check for updates, please visit
www.github.com/dcnieho/ETDQualitizer.

# Collecting validation data
The procedure for collecting validation data is available from [the `/python/ETDQualitizer/validator` subfolder](/python/ETDQualitizer/validator). In this manual, we will refer to this script as the validator.

## How to acquire
The validator is not available as an installable package. The user is recommended to download [the `/python/ETDQualitizer/validator` subfolder](/python/ETDQualitizer/validator) from github and add it to their experiment.

## Documentation
Several files are contained in [the `/python/ETDQualitizer/validator` subfolder](/python/ETDQualitizer/validator) that together comprise the validator. First, the validator can be run in two ways, either as a standalone script, or integrated into an existing PsychoPy script. By default, support is included for EyeLink, SMI and Tobii eye trackers. To run as a standalone script, simply run the `run_for_eyelink.py`, `run_for_smi.py`, or `run_for_tobii.py` scripts. To integrate the validator into an existing script, simply call the `run_validation()` function inside `validator.py`.

## Configuration files
The procedure (e.g., position of the validation targets, and duration for which they are shown) is configured using two configuration files.

### `targetPositions.csv`
The `targetPositions.csv` file contains the positions (and associated ID) for the validation targets. It contains the following columns:
|name|description|
| --- | --- |
| `ID` | Unique identifier of the validation target. |
| `x` | Horizontal position of the center of the validation target. |
| `y` | Vertical position of the center of the validation target. |

Validation target coordinates are in PsychoPy's `'pix'` coordinate system. That is, they are in pixels with respect to the center of the screen, with the positive y direction being upward.

## Output data format

## Adding your own eye tracker


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
- Wrapper functions: finally, the package provides two higher-level wrapper functions, `compute_data_quality_from_validation` and `report_data_quality_table`. `compute_data_quality_from_validation` can be used to take a segment of eye tracking formatted according to [the format produced by the validation tool](#output-data-format) and calculate the various data quality metrics for each presented target in this segment. The `report_data_quality_table` then takes one or multiple tables output by `compute_data_quality_from_validation` and provides an overall summary of data quality across targets and (if multiple tables are provided) recordings/participants, along with a textual report that can be directly pasted into a paper.

# Citation
If you use this tool or any of the code in this repository, please cite:<br>
Niehorster, D.C., Nyström, M., Hessels, R.S., Benjamins, J.S., Andersson, R. & Hooge, I.T.C. (in prep). The fundamentals of eye tracking part 7: Data quality

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
    Title = {The fundamentals of eye tracking part 7: Data quality},
    Year = {in prep},
    doi = {}
}
