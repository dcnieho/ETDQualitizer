import pandas as pd
import numpy as np
import typing

from .version import __version__, __url__, __author__, __email__, __description__

N = typing.TypeVar("N", bound=int)

def accuracy(azi: np.ndarray[tuple[N], np.dtype[np.float64]],
             ele: np.ndarray[tuple[N], np.dtype[np.float64]],
             target_azi: float, target_ele: float,
             central_tendency_fun=np.nanmean) -> tuple[float,float,float]:
    """
    Compute Gaze Accuracy

    Calculates the angular offset between gaze and target directions.

    Parameters
    ----------
    azi : np.ndarray
        Gaze azimuth values in degrees.
    ele : np.ndarray
        Gaze elevation values in degrees.
    target_azi : float
        Target azimuth in degrees.
    target_ele : float
        Target elevation in degrees.
    central_tendency_fun : callable, optional
        Function to compute central tendency (default: np.nanmean).

    Returns
    -------
    tuple of float
        A tuple with (offset, offset_azi, offset_ele), representing total, horizontal,
        and vertical offset of gaze from the target (in degrees).

    Examples
    --------
    >>> accuracy(np.array([1, 2]), np.array([1, 2]), 0, 0)
    (2.1211587518891997, 1.4995429717992865, 1.5002284247552145)
    """
    # get unit vectors for gaze and target
    g_x,g_y,g_z = Fick_to_vector(       azi,        ele)
    t_x,t_y,t_z = Fick_to_vector(target_azi, target_ele)
    # calculate angular offset for each sample using dot product
    offsets     = np.arccos(np.dot(np.vstack((g_x,g_y,g_z)).T, np.array([t_x,t_y,t_z])))
    # calculate on-screen orientation so we can decompose offset into x and y
    direction   = np.arctan2(g_y/g_z-t_y/t_z, g_x/g_z-t_x/t_z)  # compute direction on tangent screen (divide by z to project to screen at 1m)
    offsets_2D  = np.degrees(offsets.reshape((-1,1))*np.array([np.cos(direction), np.sin(direction)]).T)
    # calculate mean horizontal and vertical offset
    offset_azi  = central_tendency_fun(offsets_2D[:,0])
    offset_ele  = central_tendency_fun(offsets_2D[:,1])
    # calculate offset of centroid
    return float(np.hypot(offset_azi, offset_ele)), float(offset_azi), float(offset_ele)

def rms_s2s(azi: np.ndarray[tuple[N], np.dtype[np.float64]], ele: np.ndarray[tuple[N], np.dtype[np.float64]], central_tendency_fun=np.nanmean) -> tuple[float,float,float]:
    """
    RMS of Sample-to-Sample Differences

    Computes the root mean square (RMS) of differences between successive gaze samples.

    Parameters
    ----------
    azi : np.ndarray
        Azimuth values in degrees.
    ele : np.ndarray
        Elevation values in degrees.
    central_tendency_fun : callable, optional
        Function to compute central tendency (default: np.nanmean).

    Returns
    -------
    tuple of float
        A tuple with (rms, rms_azi, rms_ele), representing:
        - rms: Total RMS of sample-to-sample distances.
        - rms_azi: RMS of azimuthal component.
        - rms_ele: RMS of elevation component.
        All values are in degrees.

    Examples
    --------
    >>> rms_s2s(np.array([1, 2, 3]), np.array([1, 2, 3]))
    (1.4142135623730951, 1.0, 1.0)
    """
    a_diff = np.diff(azi)**2
    e_diff = np.diff(ele)**2
    # N.B.: cannot simplify to np.hypot(rms_azi, rms_ele)
    # as that is only equivalent when mean() is used as central tendency estimator
    return float(np.sqrt(central_tendency_fun(a_diff + e_diff))), \
           float(np.sqrt(central_tendency_fun(a_diff))), \
           float(np.sqrt(central_tendency_fun(e_diff)))

def std(azi: np.ndarray[tuple[N], np.dtype[np.float64]], ele: np.ndarray[tuple[N], np.dtype[np.float64]]) -> tuple[float,float,float]:
    """
    Standard Deviation of Gaze Samples

    Computes the standard deviation of azimuth and elevation gaze samples.

    Parameters
    ----------
    azi : np.ndarray
        Azimuth values in degrees.
    ele : np.ndarray
        Elevation values in degrees.

    Returns
    -------
    tuple of float
        A tuple containing:
        - std: Total standard deviation (combined azimuth and elevation).
        - std_azi: Standard deviation of azimuth.
        - std_ele: Standard deviation of elevation.
        All values are in degrees.

    Examples
    --------
    >>> std(np.array([1, 2, 3]), np.array([1, 2, 3]))
    (1.4142135623730951, 0.816496580927726, 0.816496580927726)
    """
    std_azi = np.nanstd(azi, ddof=0)
    std_ele = np.nanstd(ele, ddof=0)
    return float(np.hypot(std_azi, std_ele)), \
           float(std_azi), \
           float(std_ele)

def bcea(azi: np.ndarray[tuple[N], np.dtype[np.float64]], ele: np.ndarray[tuple[N], np.dtype[np.float64]], P: float = 0.68) -> tuple[float,float,float,float,float]:
    """
    Bivariate Contour Ellipse Area (BCEA)

    Computes BCEA and ellipse parameters for gaze precision.

    Parameters
    ----------
    azi : np.ndarray
        Azimuth values in degrees.
    ele : np.ndarray
        Elevation values in degrees.
    P : float, optional
        Cumulative probability (default: 0.68).

    Returns
    -------
    tuple of float
        A tuple with:
        - area: BCEA value.
        - orientation: Orientation of the ellipse in degrees.
        - ax1: Length of the major axis.
        - ax2: Length of the minor axis.
        - aspect_ratio: Ratio of major to minor axis.

    Examples
    --------
    >>> bcea(np.random.randn(100), np.random.randn(100))
    (7.415178348580135, -7.406813335934995, 1.1929964622017741, 0.9892420685862205, 1.2059702069754776)
    """
    k = np.log(1./(1-P))    # turn cumulative probability of area under the multivariate normal into scale factor
    azi = np.delete(azi, np.isnan(azi))
    ele = np.delete(ele, np.isnan(ele))
    std_x = np.std(azi, ddof=1)
    std_y = np.std(ele, ddof=1)
    rho = np.corrcoef(azi, ele)[0,1]
    area = 2*k*np.pi*std_x*std_y*np.sqrt(1-rho**2)
    # compute major and minor axis radii, and orientation, of the BCEA ellipse
    d,v = np.linalg.eig(np.cov(azi,ele))
    i = np.argmax(d)
    orientation = np.degrees(np.arctan2(v[1,i], v[0,i]))
    ax1 = np.sqrt(k*d[i])
    ax2 = np.sqrt(k*d[1-i])
    aspect_ratio = max([ax1, ax2])/min([ax1, ax2])
    # sanity check: this (formula for area of ellipse) should
    # closely match directly computed area from above
    # 2*np.pi*ax1*ax2
    return float(area), float(orientation), float(ax1), float(ax2), float(aspect_ratio)

def data_loss(azi: np.ndarray[tuple[N], np.dtype[np.float64]], ele: np.ndarray[tuple[N], np.dtype[np.float64]]):
    """
    Compute Data Loss Percentage

    Calculates the percentage of missing gaze samples.

    Parameters
    ----------
    azi : np.ndarray
        Azimuth values.
    ele : np.ndarray
        Elevation values.

    Returns
    -------
    float
        Percentage of missing samples.

    Examples
    --------
    >>> data_loss(np.array([1, np.nan, 3]), np.array([1, 2, np.nan]))
    66.666
    """
    missing = np.isnan(azi) | np.isnan(ele)
    return np.sum(missing)/missing.size*100

def data_loss_from_expected(azi: np.ndarray[tuple[N], np.dtype[np.float64]], ele: np.ndarray[tuple[N], np.dtype[np.float64]], duration: float, frequency: float):
    """
    Compute Data Loss from Expected Sample Count

    Calculates data loss based on expected number of samples.

    Parameters
    ----------
    azi : np.ndarray
        Azimuth values.
    ele : np.ndarray
        Elevation values.
    duration : float
        Duration in seconds.
    frequency : float
        Sampling frequency in Hz.

    Returns
    -------
    float
        Percentage of data loss.

    Examples
    --------
    >>> data_loss_from_expected(np.array([1, np.nan, 3]), np.array([1, 2, np.nan]), duration=1, frequency=3)
    33.333
    """
    N_valid = np.count_nonzero(~(np.isnan(azi) | np.isnan(ele)))
    return (1-N_valid/(duration*frequency))*100

def effective_frequency(azi: np.ndarray[tuple[N], np.dtype[np.float64]], ele: np.ndarray[tuple[N], np.dtype[np.float64]], duration: float):
    """
    Compute Effective Sampling Frequency

    Calculates effective frequency based on valid samples.

    Parameters
    ----------
    azi : np.ndarray
        Azimuth values.
    ele : np.ndarray
        Elevation values.
    duration : float
        Duration in seconds.

    Returns
    -------
    float
        Effective frequency in Hz.

    Examples
    --------
    >>> effective_frequency(np.array([1, np.nan, 3]), np.array([1, 2, np.nan]), duration=1)
    2.0
    """
    N_valid = np.count_nonzero(~(np.isnan(azi) | np.isnan(ele)))
    return N_valid/duration


def precision_using_moving_window(azi: np.ndarray[tuple[N], np.dtype[np.float64]], ele: np.ndarray[tuple[N], np.dtype[np.float64]],
                                  window_length: int, metric: str, aggregation_fun=np.nanmedian, **kwargs) -> float:
    """
    Precision Using Moving Window

    Computes gaze precision using a moving window and selected metric.

    Parameters
    ----------
    azi : np.ndarray
        Azimuth values.
    ele : np.ndarray
        Elevation values.
    window_length : int
        Window size in samples.
    metric : str
        Precision metric: "RMS-S2S", "STD", or "BCEA".
    aggregation_fun : callable, optional
        Function to aggregate precision values across the windows (default: np.nanmedian).
    **kwargs : dict
        Additional arguments passed to the metric function.

    Returns
    -------
    float
        Aggregated precision value.

    Examples
    --------
    >>> precision_using_moving_window(np.random.randn(100), np.random.randn(100), 10, "STD")
    1.2  # example output
    """

    match metric:
        case 'RMS-S2S':
            fun =  rms_s2s
        case 'STD':
            fun =  std
        case 'BCEA':
            fun =  bcea
        case _:
            raise ValueError(f'metric "{metric}" is not understood')

    # get number of samples in data
    ns  = azi.shape[0]

    if window_length < ns:  # if number of samples in data exceeds window size
        values = np.full((ns-window_length+1,), np.nan)  # pre-allocate
        for p in range(0,ns-window_length+1):
            values[p] = fun(azi[p:p+window_length], ele[p:p+window_length], **kwargs)[0]
        precision = aggregation_fun(values)
    else:
        # if too few samples in data
        precision = np.nan
    return precision


class ScreenConfiguration:
    """
    ScreenConfiguration

    Provides methods for converting between pixel, millimeter, and degree units.

    Parameters
    ----------
    screen_size_x_mm : float
        Screen width in millimeters.
    screen_size_y_mm : float
        Screen height in millimeters.
    screen_res_x_pix : int
        Horizontal screen resolution in pixels.
    screen_res_y_pix : int
        Vertical screen resolution in pixels.
    viewing_distance_mm : float
        Viewing distance in millimeters.

    Examples
    --------
    >>> sc = ScreenConfiguration(500, 300, 1920, 1080, 600)
    >>> sc.pix_to_deg(960, 540)
    """
    def __init__(self,
                 screen_size_x_mm: float, screen_size_y_mm: float,
                 screen_res_x_pix: int  , screen_res_y_pix: int,
                 viewing_distance_mm: float):
        self.screen_size_x_mm   = screen_size_x_mm
        self.screen_size_y_mm   = screen_size_y_mm
        self.screen_res_x_pix   = screen_res_x_pix
        self.screen_res_y_pix   = screen_res_y_pix
        self.viewing_distance_mm= viewing_distance_mm

    def pix_to_mm(self, x: float, y: float) -> tuple[float,float]:
        """
        Convert pixel coordinates to millimeter coordinates on the screen.

        Parameters
        ----------
        x : float
            Horizontal pixel coordinate.
        y : float
            Vertical pixel coordinate.

        Returns
        -------
        tuple of float
            x and y in millimeters.
        """
        x_mm = x/self.screen_res_x_pix*self.screen_size_x_mm
        y_mm = y/self.screen_res_y_pix*self.screen_size_y_mm
        return x_mm, y_mm

    def pix_to_deg(self, x: float, y: float) -> tuple[float,float]:
        """
        Convert pixel coordinates to angular gaze direction in degrees (Fick angles).

        Parameters
        ----------
        x : float
            Horizontal pixel coordinate.
        y : float
            Vertical pixel coordinate.

        Returns
        -------
        tuple of float
            Azimuth and elevation in degrees.
        """
        # N.B.: output is in Fick angles
        x_mm, y_mm = self.pix_to_mm(x, y)
        return self.mm_to_deg(x_mm, y_mm)

    def mm_to_deg(self, x: float, y: float) -> tuple[float,float]:
        """
        Convert millimeter coordinates to angular gaze direction in degrees (Fick angles).

        Parameters
        ----------
        x : float
            Horizontal position in millimeters.
        y : float
            Vertical position in millimeters.

        Returns
        -------
        tuple of float
            Azimuth and elevation in degrees.
        """
        # N.B.: output is in Fick angles
        azi = np.arctan2(x,self.viewing_distance_mm)
        ele = np.arctan2(y,np.hypot(self.viewing_distance_mm,x))
        return np.degrees(azi), np.degrees(ele)

    def mm_to_pix(self, x: float, y: float) -> tuple[float,float]:
        """
        Convert millimeter coordinates to pixel coordinates.

        Parameters
        ----------
        x : float
            Horizontal position in millimeters.
        y : float
            Vertical position in millimeters.

        Returns
        -------
        tuple of float
            x and y in pixels.
        """
        x_pix = x/self.screen_size_x_mm*self.screen_res_x_pix
        y_pix = y/self.screen_size_y_mm*self.screen_res_y_pix
        return x_pix, y_pix

    def deg_to_pix(self, azi: float, ele: float) -> tuple[float,float]:
        """
        Convert angular gaze direction in degrees (Fick angles) to pixel coordinates.

        Parameters
        ----------
        azi : float
            Azimuth in degrees.
        ele : float
            Elevation in degrees.

        Returns
        -------
        tuple of float
            x and y in pixels.
        """
        x_mm, y_mm = self.deg_to_mm(azi, ele)
        return self.mm_to_pix(x_mm, y_mm)

    def deg_to_mm(self, azi: float, ele: float) -> tuple[float,float]:
        """
        Convert angular gaze direction in degrees (Fick angles) to millimeter coordinates.

        Parameters
        ----------
        azi : float
            Azimuth in degrees.
        ele : float
            Elevation in degrees.

        Returns
        -------
        tuple of float
            x and y in millimeters.
        """
        azi = np.radians(azi)
        ele = np.radians(ele)
        x_mm = self.viewing_distance_mm*np.tan(azi)
        y_mm = self.viewing_distance_mm*np.tan(ele)/np.cos(azi)
        return x_mm, y_mm

    def screen_extents(self) -> tuple[float,float]:
        """
        Compute the horizontal and vertical extents of the screen in degrees.

        Returns
        -------
        tuple of float
            Width and height in degrees.
        """
        [x_deg, _] = self.mm_to_deg(self.screen_size_x_mm/2, 0)
        [_, y_deg] = self.mm_to_deg(0, self.screen_size_y_mm/2)
        return x_deg*2, y_deg*2


def Fick_to_vector(azi: np.ndarray[tuple[N], np.dtype[np.float64]]|float, ele: np.ndarray[tuple[N], np.dtype[np.float64]]|float, r: np.ndarray[tuple[N], np.dtype[np.float64]]|float=1.) -> tuple[np.ndarray[tuple[N], np.dtype[np.float64]], np.ndarray[tuple[N], np.dtype[np.float64]], np.ndarray[tuple[N], np.dtype[np.float64]]]:
    """
    Convert Fick Angles to 3D Vector

    Converts azimuth and elevation angles (in degrees) to a 3D unit vector.

    Parameters
    ----------
    azi : float or np.ndarray
        Azimuth angle(s) in degrees.
    ele : float or np.ndarray
        Elevation angle(s) in degrees.
    r : float or np.ndarray, optional
        Radius (default is 1.0).

    Returns
    -------
    tuple of np.ndarray
        A tuple with components (x, y, z) representing the 3D vector.

    Examples
    --------
    >>> Fick_to_vector(30, 10)
    (np.float64(0.49240387650610395), np.float64(0.17364817766693033), np.float64(0.8528685319524433))
    """
    azr = np.radians(azi)
    elr = np.radians(ele)
    r_cos_ele = r*np.cos(elr)

    x = r_cos_ele * np.sin(azr)
    y =         r * np.sin(elr)
    z = r_cos_ele * np.cos(azr)
    return x,y,z

def vector_to_Fick(x: np.ndarray[tuple[N], np.dtype[np.float64]]|float, y: np.ndarray[tuple[N], np.dtype[np.float64]]|float, z: np.ndarray[tuple[N], np.dtype[np.float64]]|float) -> tuple[np.ndarray[tuple[N], np.dtype[np.float64]], np.ndarray[tuple[N], np.dtype[np.float64]]]:
    """
    Convert 3D Vector to Fick Angles

    Converts a 3D vector to azimuth and elevation angles (in degrees).

    Parameters
    ----------
    x : float or np.ndarray
        X component of the vector.
    y : float or np.ndarray
        Y component of the vector.
    z : float or np.ndarray
        Z component of the vector.

    Returns
    -------
    tuple of np.ndarray
        A tuple with components (azi, ele) in degrees.

    Examples
    --------
    >>> vector_to_Fick(0.5, 0.2, 0.8)
    (np.float64(32.005383208083494), np.float64(11.969463124607309))
    """
    return np.degrees(np.arctan2(x,z)), \
           np.degrees(np.arctan2(y,np.hypot(x,z)))


class DataQuality:
    """
    DataQuality

    Provides methods for assessing the quality of gaze data, including accuracy, precision,
    data loss, and effective sampling frequency.

    Notes
    -----
    - Missing data should be coded as NaN, not as special values like (0,0) or (-xres,-yres).
    - Missing samples should not be removed, or RMS calculations will be incorrect.
    - Timestamps should be in seconds.
    - All angular positions are expected to be expressed in Fick angles.

    Parameters
    ----------
    gaze_x : np.ndarray
        Horizontal gaze positions (pixels or degrees).
    gaze_y : np.ndarray
        Vertical gaze positions (pixels or degrees).
    timestamps : np.ndarray
        Vector of timestamps in seconds.
    unit : str
        Unit of gaze data: either "pixels" or "degrees".
    screen : ScreenConfiguration, optional
        Required if unit is "pixels".

    Examples
    --------
    >>> sc = ScreenConfiguration(500, 300, 1920, 1080, 600)
    >>> dq = DataQuality([0, 1, -1], [0, 1, -1], [0, 1, 2], unit="pixels", screen=sc)
    >>> dq.accuracy(0, 0)
    >>> dq.precision_RMS_S2S()
    >>> dq.data_loss()
    """
    def __init__(self,
                 gaze_x     : np.ndarray[tuple[N], np.dtype[np.float64]],
                 gaze_y     : np.ndarray[tuple[N], np.dtype[np.float64]],
                 timestamps : np.ndarray[tuple[N], np.dtype[np.float64]],
                 unit       : str,
                 screen     : ScreenConfiguration|None = None):
        self.timestamps = np.array(timestamps)

        gaze_x = np.array(gaze_x)
        gaze_y = np.array(gaze_y)
        if unit=='pixels':
            if screen is None:
                raise ValueError('If unit is "pixels", a screen configuration must be supplied')
            gaze_x, gaze_y = screen.pix_to_deg(gaze_x, gaze_y)
        elif unit!='degrees':
            raise ValueError('unit should be "pixels" or "degrees"')
        self.azi = gaze_x
        self.ele = gaze_y

    def accuracy(self, target_azi: float, target_ele: float, central_tendency_fun=np.nanmean) -> tuple[float,float,float]:
        """
        Calculates the accuracy of gaze data relative to a known target location.

        Returns
        -------
        tuple of float
            Total, horizontal, and vertical accuracy in degrees.
        """
        # get unit vectors for gaze and target
        return accuracy(self.azi, self.ele, target_azi, target_ele, central_tendency_fun)

    def precision_RMS_S2S(self, central_tendency_fun=np.nanmean) -> tuple[float,float,float]:
        """
        Calculates precision as root mean square of sample-to-sample distances.

        Returns
        -------
        tuple of float
            Total, azimuthal, and elevation RMS precision in degrees.
        """
        return rms_s2s(self.azi, self.ele, central_tendency_fun)

    def precision_STD(self) -> tuple[float,float,float]:
        """
        Calculates precision as standard deviation of gaze positions.

        Returns
        -------
        tuple of float
            Total, azimuthal, and elevation standard deviation in degrees.
        """
        return std(self.azi, self.ele)

    def precision_BCEA(self, P: float = 0.68) -> tuple[float,float,float,float,float]:
        """
        Calculates the Bivariate Contour Ellipse Area (BCEA) and ellipse parameters.

        Parameters
        ----------
        P : float
            Proportion of data to include in the ellipse.

        Returns
        -------
        tuple of float
            BCEA area, orientation, major axis, minor axis, and aspect ratio.
        """
        return bcea(self.azi, self.ele, P)

    def data_loss(self):
        """
        Calculates the proportion of missing data (coded as NaN).

        Returns
        -------
        float
            Percentage of missing samples.
        """
        return data_loss(self.azi, self.ele)

    def data_loss_from_expected(self, frequency):
        """
        Estimates data loss based on expected number of samples.

        Parameters
        ----------
        frequency : float
            Expected sampling frequency in Hz.

        Returns
        -------
        float
            Percentage of missing samples.
        """
        return data_loss_from_expected(self.azi, self.ele, self.get_duration(), frequency)

    def effective_frequency(self):
        """
        Calculates the effective sampling frequency based on timestamps.

        Returns
        -------
        float
            Effective frequency in Hz.
        """
        return effective_frequency(self.azi, self.ele, self.get_duration())

    def get_duration(self) -> float:
        """
        Computes the total duration of the gaze recording, including the last sample.

        Returns
        -------
        float
            Duration in seconds.
        """
        # to get duration right, we need to include duration of last sample
        isi = np.median(np.diff(self.timestamps))
        return self.timestamps[-1]-self.timestamps[0]+isi


    def precision_using_moving_window(self, window_length, metric, aggregation_fun=np.nanmedian, **kwargs) -> float:
        """
        Calculates precision using a moving window approach.

        Parameters
        ----------
        window_length : int
            Length of the moving window in number of samples.
        metric : str
            Precision metric to use ("RMS-S2S", "STD", or "BCEA").
        aggregation_fun : callable
            Function to aggregate windowed precision values.
        **kwargs : dict
            Additional arguments passed to the precision metric function.

        Returns
        -------
        float
            Aggregated precision value.
        """
        return precision_using_moving_window(self.azi, self.ele, window_length, metric, aggregation_fun, **kwargs)


def compute_data_quality_from_validation(gaze               : pd.DataFrame,
                                         unit               : str,
                                         screen             : ScreenConfiguration|None = None,
                                         advanced           : bool = False, # if True, report all metrics. If False, only simple subset
                                         include_data_loss  : bool = False) -> pd.DataFrame:
    # get all targets
    targets         = sorted([t for t in gaze['target_id'].unique() if t!=-1])
    target_locations= np.array([gaze.loc[gaze.index[(gaze['target_id'].values==t).argmax()], ['tar_x','tar_y']] for t in targets])

    # ensure we have target locations in degrees
    if unit=='pixels':
        if screen is None:
            raise ValueError('If unit is "pixels", a screen configuration must be supplied')
        target_locations[:,0], target_locations[:,1] = screen.pix_to_deg(target_locations[:,0], target_locations[:,1])
    elif unit!='degrees':
        raise ValueError('unit should be "pixels" or "degrees"')

    # now, per target, compute data quality metrics
    rows = []
    for e in ('left','right'):
        if f'{e}_x' not in gaze.columns:
            continue
        for i,t_id in enumerate(targets):
            is_target = gaze['target_id'].values==t_id
            dq = DataQuality(gaze[f'{e}_x'][is_target], gaze[f'{e}_y'][is_target], gaze['timestamp'][is_target]/1000, unit, screen) # timestamps are in ms in the file
            row = {'eye': e, 'target_id': t_id}
            for k,v in zip(('offset','offset_x','offset_y'),dq.accuracy(*target_locations[i])):
                row[k] = v
            for k,v in zip(('rms_s2s','rms_s2s_x','rms_s2s_y'),dq.precision_RMS_S2S()):
                row[k] = v
            for k,v in zip(('std','std_x','std_y'),dq.precision_STD()):
                row[k] = v
            for k,v in zip(('bcea','bcea_orientation','bcea_ax1','bcea_ax2','bcea_aspect_ratio'),dq.precision_BCEA()):
                row[k] = v
            if include_data_loss:
                row['data_loss'] = dq.data_loss()
                row['effective_frequency'] = dq.effective_frequency()
            rows.append(row)

    dq_df = pd.DataFrame.from_records(rows).set_index(['eye','target_id'])
    if not advanced:
        dq_df = dq_df.drop(columns=[c for c in dq_df.columns if c not in ('eye', 'target_id', 'offset', 'rms_s2s', 'std', 'bcea', 'data_loss', 'effective_frequency')])
    return dq_df


def report_data_quality_table(dq_table: pd.DataFrame) -> tuple[str,dict[str,pd.DataFrame]]:
    measures = {}
    # average over targets and eyes
    measures['all']  = dq_table.groupby('file').mean()

    # do summary statistics
    m = {}
    m['mean'] = measures['all'].mean()
    m['std']  = measures['all'].std()
    m['min']  = measures['all'].min()
    m['max']  = measures['all'].max()
    measures['summary'] = pd.DataFrame(m).T

    # make text. A little overcomplete, user can trim what they don't want
    # N.B.: do not include data loss/effective frequency, nor bcea. Bcea is
    # niche, user who wants it can get that themselves. Data loss you'd really
    # want to report for all the analyzed data, not just this validation
    # procedure.
    n_target = dq_table.index.get_level_values('target_id').nunique()
    n_subj   = measures['all'].shape[0]
    txt = (
        f"For {n_subj} participants, the average inaccuracy in the data determined from a {n_target}-point validation procedure using ETDQualitizer v{__version__} (Niehorster et al., in prep) "
        f"was {measures['summary'].loc['mean'].offset:.2f}° (SD={measures['summary'].loc['std'].offset:.2f}°, range={measures['summary'].loc['min'].offset:.2f}°--{measures['summary'].loc['max'].offset:.2f}°). "
        f"Average RMS-S2S precision was {measures['summary'].loc['mean'].rms_s2s:.3f}° (SD={measures['summary'].loc['std'].rms_s2s:.3f}°, range={measures['summary'].loc['min'].rms_s2s:.3f}°--{measures['summary'].loc['max'].rms_s2s:.3f}°) "
        f"and STD precision {measures['summary'].loc['mean']['std']:.3f}° (SD={measures['summary'].loc['std']['std']:.3f}°, range={measures['summary'].loc['min']['std']:.3f}°--{measures['summary'].loc['max']['std']:.3f}°)."
    )
    return txt, measures