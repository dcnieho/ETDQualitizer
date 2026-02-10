import json
import traceback
import pathlib
from py_smite import SMITE

from validator import run_validation, open_screen

def main():
    # read protocol setup
    with open(pathlib.Path(__file__).parent.resolve()/"setup.json") as fp:
        config = json.load(fp)

    try:
        # Open window, check
        win = open_screen(config)

        # set up which SMI eye tracker to use
        settings = SMITE.get_defaults('RED')
        settings.autoaccept = 1
        settings.n_cal_points = 5
        tracker = SMITE.Connect(settings)
        tracker.init()
        tracker.calibrate(win)

        file_name = run_validation(win, config, tracker)

    except Exception as e:
        tb_lines = traceback.format_exception(type(e), e, e.__traceback__)
        print("".join(tb_lines))
    finally:
        if 'win' in locals():
            win.close()

if __name__=="__main__":
    main()