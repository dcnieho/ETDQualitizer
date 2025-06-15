import pathlib
import numpy as np
import pandas as pd
from collections import defaultdict

import ETDQualitizer
print(f'Using ETDQualitizer version {ETDQualitizer.__version__}')

screen = ETDQualitizer.ScreenConfiguration(
    screen_size_x_mm = 528.0, screen_size_y_mm = 296.9997253417969,
    screen_res_x_pix = 1920, screen_res_y_pix = 1080,
    viewing_distance_mm = 650
)

# per data file, run the analysis
all_dfs: dict[str,pd.DataFrame] = {}
for f in (pathlib.Path(__file__).parent / 'data').glob('*.tsv'):
    print(f'----------\n{f.name}')
    gaze = pd.read_csv(f, sep='\t', dtype=defaultdict(lambda: float, {'target_id': int, 'tar_x': int, 'tar_y': int}))

    # automatically compute data quality measures per target
    dq_df = ETDQualitizer.compute_data_quality_from_validation(gaze, 'pixels', screen, advanced=False, include_data_loss=True)    # include_data_loss for testing, this is probably *not* what you want
    print(dq_df.to_string(float_format='%.4f'))
    all_dfs[f.name] = dq_df

    # manually perform some further data quality computations on data from
    # the whole validation instead of per target
    for e in ('left','right'):
        if f'{e}_x' not in gaze.columns:
            continue
        # and RMS-S2S calculated in two ways over the whole datafile
        dq = ETDQualitizer.DataQuality(gaze[f'{e}_x'].to_numpy(),gaze[f'{e}_y'].to_numpy(),gaze['timestamp'].to_numpy()/1000,'pixels',screen)  # timestamps are in ms in the file

        # determine sampling frequency from the filename (assumes our test
        # files which end in '<xxx>Hz')
        fs = int(f.stem.split('_')[-1].removesuffix('Hz'))
        window_len = int(.2*fs) # 200 ms

        print(f'RMS S2S using median ({e} eye): {dq.precision_RMS_S2S(central_tendency_fun=np.nanmedian)[0]:.4f} deg')
        print(f'RMS S2S using moving window ({e} eye): {dq.precision_using_moving_window(window_len,"RMS_S2S"):.4f} deg')

        # data loss and effective frequency
        print(f'Data loss ({e} eye): {dq.data_loss():.1f}%')
        print(f'Data loss ({e} eye): {dq.data_loss_nominal(fs):.1f}%')
        print(f'Effective frequency ({e} eye): {dq.effective_frequency():.1f} Hz')

all_df = pd.concat(all_dfs, names=['file'])

# make a text one can directly put in a paper. Note that normally all_df
# would contain data for multiple subjects measured under the same
# conditions, not like in this case data from different devices with
# different settings
print('--------\nAll:')
dq_txt, summary_dq = ETDQualitizer.report_data_quality_table(all_df)
print(dq_txt)
summary_dq['all']