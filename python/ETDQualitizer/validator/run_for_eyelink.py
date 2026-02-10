import json
import traceback
import pathlib
from eyelink import eyelink

from validator import run_validation, open_screen

def main():
    # read protocol setup
    with open(pathlib.Path(__file__).parent.resolve()/"setup.json") as fp:
        config = json.load(fp)

    try:
        # Open window, check
        win = open_screen(config)

        # set up the EyeLink eye tracker to use
        settings = eyelink.Settings()
        settings.FILENAME = 'test'
        tracker = eyelink.Connect(settings, use_sample_buffer=True, sample_buffer_length=None)
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