import json
import traceback
import pathlib
from titta import Titta

from validator import run_validation, open_demo_screen

def main():
    # read protocol setup
    with open(pathlib.Path(__file__).parent.resolve()/"setup.json") as fp:
        config = json.load(fp)

    try:
        # Open window, check
        win = open_demo_screen(config)

        # use a Tobii eye tracker
        settings = Titta.get_defaults('Tobii Pro Spectrum')
        tracker = Titta.Connect(settings)
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