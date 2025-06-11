import pandas as pd
import numpy as np
import typing

from .version import __version__, __url__, __author__, __email__, __description__

N = typing.TypeVar("N", bound=int)

class ScreenConfiguration:
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
        x_mm = x/self.screen_res_x_pix*self.screen_size_x_mm
        y_mm = y/self.screen_res_y_pix*self.screen_size_y_mm
        return x_mm, y_mm

    def pix_to_deg(self, x: float, y: float) -> tuple[float,float]:
        # N.B.: output is in Fick angles
        x_mm, y_mm = self.pix_to_mm(x, y)
        return self.mm_to_deg(x_mm, y_mm)

    def mm_to_deg(self, x: float, y: float) -> tuple[float,float]:
        # N.B.: output is in Fick angles
        azi = np.arctan2(x,self.viewing_distance_mm)
        ele = np.arctan2(y,np.hypot(self.viewing_distance_mm,x))
        return np.degrees(azi), np.degrees(ele)


def _Fick_to_cartesian(azi: np.ndarray[tuple[N], np.dtype[np.float64]], ele: np.ndarray[tuple[N], np.dtype[np.float64]], r: float=1.) -> tuple[np.ndarray[tuple[N], np.dtype[np.float64]], np.ndarray[tuple[N], np.dtype[np.float64]], np.ndarray[tuple[N], np.dtype[np.float64]]]:
    azi = np.radians(azi)
    ele = np.radians(ele)
    r_cos_ele = r*np.cos(ele)

    x = r_cos_ele * np.sin(azi)
    y =         r * np.sin(ele)
    z = r_cos_ele * np.cos(azi)
    return x,y,z


class DataQuality:
    # N.B: for this module it is assumed that any missing data are not coded with some special value
    # such as (0,0) or (-xres,-yres) but as nan. Missing data should also not be removed, or the RMS
    # calculation would be incorrect.
    #
    # timestamps should be in seconds.
    #
    # all angular positions are expected to be expressed in Fick angles.
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
        self.x = gaze_x
        self.y = gaze_y

    def accuracy(self, target_x_deg: float, target_y_deg: float, central_tendency_fun=np.nanmean) -> tuple[float,float,float]:
        # get unit vectors for gaze and target
        g_x,g_y,g_z = _Fick_to_cartesian(  self.x,       self.y)
        t_x,t_y,t_z = _Fick_to_cartesian(target_x_deg, target_y_deg)
        # calculate angular offset for each sample using dot product
        offsets     = np.arccos(np.dot(np.vstack((g_x,g_y,g_z)).T, np.array([t_x,t_y,t_z])))
        # calculate on-screen orientation so we can decompose offset into x and y
        direction   = np.arctan2(g_y/g_z-t_y/t_z, g_x/g_z-t_x/t_z)  # compute direction on tangent screen (divide by z to project to screen at 1m)
        offsets_2D  = np.degrees(offsets.reshape((-1,1))*np.array([np.cos(direction), np.sin(direction)]).T)
        # calculate mean horizontal and vertical offset
        offset_x    = central_tendency_fun(offsets_2D[:,0])
        offset_y    = central_tendency_fun(offsets_2D[:,1])
        # calculate offset of centroid
        return float(np.hypot(offset_x, offset_y)), float(offset_x), float(offset_y)

    def precision_RMS_S2S(self, central_tendency_fun=np.nanmean) -> tuple[float,float,float]:
        return _RMS_S2S_impl(self.x, self.y, central_tendency_fun)

    def precision_STD(self) -> tuple[float,float,float]:
        return _STD_impl(self.x, self.y)

    def precision_BCEA(self, P: float = 0.6827) -> tuple[float,float,float,float,float]:
        return _BCEA_impl(self.x, self.y, P)

    def data_loss_percentage(self):
        missing = np.isnan(self.x) | np.isnan(self.y)
        return np.sum(missing)/missing.size*100

    def data_loss_percentage_nominal(self, frequency):
        N_valid = np.count_nonzero(~(np.isnan(self.x) | np.isnan(self.y)))
        return (1-N_valid/(self.get_duration()*frequency))*100

    def effective_frequency(self):
        N_valid = np.count_nonzero(~(np.isnan(self.x) | np.isnan(self.y)))
        return N_valid/self.get_duration()

    def get_duration(self) -> float:
        # to get duration right, we need to include duration of last sample
        isi = np.median(np.diff(self.timestamps))
        return self.timestamps[-1]-self.timestamps[0]+isi


    def precision_using_moving_window(self, window_length, metric, aggregation_fun=np.nanmedian, **kwargs) -> float:
        match metric:
            case 'RMS_S2S':
                fun =  _RMS_S2S_impl
            case 'STD':
                fun =  _STD_impl
            case 'BCEA':
                fun =  _BCEA_impl
            case _:
                raise ValueError(f'metric "{metric}" is not understood')

        # get number of samples in data
        ns  = self.x.shape[0]

        if window_length < ns:  # if number of samples in data exceeds window size
            values = np.full((ns-window_length+1,), np.nan)  # pre-allocate
            for p in range(0,ns-window_length+1):
                values[p] = fun(self.x[p:p+window_length], self.y[p:p+window_length], **kwargs)[0]
            precision = aggregation_fun(values)
        else:
            # if too few samples in data
            precision = np.nan
        return precision


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
                row['data_loss'] = dq.data_loss_percentage()
                row['effective_frequency'] = dq.effective_frequency()
            rows.append(row)

    dq_df = pd.DataFrame.from_records(rows).set_index(['eye','target_id'])
    if not advanced:
        dq_df = dq_df.drop(columns=[c for c in dq_df.columns if c not in ('eye', 'target_id', 'offset', 'rms_s2s', 'std', 'bcea', 'data_loss', 'effective_frequency')])
    return dq_df



def _RMS_S2S_impl(x: np.ndarray[tuple[N], np.dtype[np.float64]], y: np.ndarray[tuple[N], np.dtype[np.float64]], central_tendency_fun=np.nanmean) -> tuple[float,float,float]:
    x_diff = np.diff(x)**2
    y_diff = np.diff(y)**2
    # N.B.: cannot simplify to np.hypot(rms_x, rms_y)
    # as that is only equivalent when mean() is used as central tendency estimator
    return float(np.sqrt(central_tendency_fun(x_diff + y_diff))), \
           float(np.sqrt(central_tendency_fun(x_diff))), \
           float(np.sqrt(central_tendency_fun(y_diff)))

def _STD_impl(x: np.ndarray[tuple[N], np.dtype[np.float64]], y: np.ndarray[tuple[N], np.dtype[np.float64]]) -> tuple[float,float,float]:
    std_x = np.nanstd(x, ddof=0)
    std_y = np.nanstd(y, ddof=0)
    return float(np.hypot(std_x, std_y)), \
           float(std_x), \
           float(std_y)

def _BCEA_impl(x: np.ndarray[tuple[N], np.dtype[np.float64]], y: np.ndarray[tuple[N], np.dtype[np.float64]], P: float = 0.68) -> tuple[float,float,float,float,float]:
    k = np.log(1./(1-P))    # turn cumulative probability of area under the multivariate normal into scale factor

    x = np.delete(x, np.isnan(x))
    y = np.delete(y, np.isnan(y))
    std_x = np.std(x, ddof=1)
    std_y = np.std(y, ddof=1)
    rho = np.corrcoef(x, y)[0,1]
    area = 2*k*np.pi*std_x*std_y*np.sqrt(1-rho**2)
    # compute major and minor axis radii, and orientation, of the BCEA ellipse
    d,v = np.linalg.eig(np.cov(x,y))
    i = np.argmax(d)
    orientation = np.degrees(np.arctan2(v[1,i], v[0,i]))
    ax1 = np.sqrt(k*d[i])
    ax2 = np.sqrt(k*d[1-i])
    aspect_ratio = max([ax1, ax2])/min([ax1, ax2])
    # sanity check: this (formula for area of ellipse) should
    # closely match directly computed area from above
    # 2*np.pi*ax1*ax2
    return float(area), float(orientation), float(ax1), float(ax2), float(aspect_ratio)