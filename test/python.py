import pathlib
import numpy as np
import pandas as pd
from collections import defaultdict

import ETDQualitizer

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

    dq_df = ETDQualitizer.compute_data_quality_from_validation(gaze, 'pixels', screen, advanced=False, include_data_loss=True)    # include_data_loss for testing, this is probably *not* what you want
    print(dq_df.to_string(float_format='%.4f'))
    all_dfs[f.name] = dq_df

    for e in ('left','right'):
        if f'{e}_x' not in gaze.columns:
            continue
        # and RMS S2S calculated in two ways over the whole datafile
        dq = ETDQualitizer.DataQuality(gaze[f'{e}_x'].to_numpy(),gaze[f'{e}_y'].to_numpy(),gaze['timestamp'].to_numpy()/1000,'pixels',screen)  # timestamps are in ms in the file

        fs = int(f.stem.removesuffix('Hz'))
        window_len = int(.2*fs) # 200 ms

        print(f'RMS S2S using median ({e} eye): {dq.precision_RMS_S2S(central_tendency_fun=np.nanmedian)[0]:.4f} deg')
        print(f'RMS S2S using moving window ({e} eye): {dq.precision_using_moving_window(window_len,"RMS_S2S"):.4f} deg')

        # data loss and effective frequency
        print(f'Data loss ({e} eye): {dq.data_loss_percentage():.1f}%')
        print(f'Effective frequency ({e} eye): {dq.effective_frequency():.1f} Hz')

print('--------\nAll:')
all_df = pd.concat(all_dfs, names=['file'])
print(all_df.to_string(float_format='%.4f'))