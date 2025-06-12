# -*- coding: utf-8 -*-
from psychopy import visual, monitors, core, event
import numpy as np
import pandas as pd
import random
from collections import defaultdict
import pathlib
import abc
import time


class EyeTrackerBase(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def start_recording(self):
        pass
    @abc.abstractmethod
    def send_message(self, msg:str):
        pass
    @abc.abstractmethod
    def stop_recording(self):
        pass
    @abc.abstractmethod
    def save_data(self, file_stem: str, mon: monitors.Monitor) -> str:
        pass

class TobiiEyeTracker(EyeTrackerBase):
    def __init__(self, tracker: 'Tobii.myTobii'):
        super().__init__()
        self.tracker = tracker
        self.start_ts: int = None
        self.stop_ts: int = None
    def start_recording(self):
        self.tracker.start_recording(gaze=True)
        self.start_ts = self.tracker.get_system_time_stamp()
        core.wait(0.5)
    def send_message(self, msg:str):
        self.tracker.send_message(msg)
    def stop_recording(self):
        self.stop_ts = self.tracker.get_system_time_stamp()
        self.tracker.stop_recording(gaze=True)
    def save_data(self, file_stem: str, mon: monitors.Monitor) -> str:
        samples = self.tracker.buffer.peek_time_range('gaze', self.start_ts, self.stop_ts)
        messages= [msg for msg in self.tracker.msg_container if msg[0]>=self.start_ts and msg[0]<=self.stop_ts]

        # parse target info
        msgs = [(m[0],pm) for m in messages if (pm:=parse_target_msg(m[1])) is not None]

        # collect and organize data
        to_include = np.logical_and(samples['system_time_stamp']>=msgs[0][0], samples['system_time_stamp']<=msgs[-1][0])
        # output: [timestamp, left_x, left_y, right_x, right_y, target_id, tar_x, tar_y] (timestamp in ms)
        # get all et data
        n_samp = sum(to_include)
        data = np.full((n_samp,5),np.nan)
        px,py = mon.currentCalib['sizePix']
        for i,(field,fac,sub) in enumerate((('system_time_stamp',0.001,0.),('left_gaze_point_on_display_area_x',px,px/2),('left_gaze_point_on_display_area_y',py,py/2),('right_gaze_point_on_display_area_x',px,px/2),('right_gaze_point_on_display_area_y',py,py/2))):
            data[:,i] = samples[field][to_include]*fac-sub
        # turn into dataframe
        data = pd.DataFrame(data,columns=['timestamp', 'left_x', 'left_y', 'right_x', 'right_y'])
        # add target ID
        target_ids_pos = np.full((n_samp,3),-1,'int32')
        ts = samples['system_time_stamp'][to_include]
        for ti in range(0,len(msgs),2):
            is_target = np.logical_and(ts>=msgs[ti][0], ts<=msgs[ti+1][0])
            target_ids_pos[is_target,:] = msgs[ti][1][2:]
        target_ids_pos[:,2] = -target_ids_pos[:,2]  # positive y direction should be downward, not upwards like only PsychoPy does
        data[['target_id','tar_x','tar_y']] = target_ids_pos

        # done, now store to file. First get filename
        file_name = get_filename(file_stem,'.tsv')
        # store to file
        data.to_csv(file_name, '\t', na_rep='nan', index=False, float_format='%.8f')

        return file_name

class EyeLinkTracker(EyeTrackerBase):
    def __init__(self, tracker: 'eyelink.Connect'):
        super().__init__()
        self.tracker = tracker
        self.start_ts: int = None
        self.stop_ts: int = None
    def start_recording(self):
        self.tracker.start_recording(sendlink=True)
        self.start_ts = self.tracker.get_time_stamp()
        core.wait(0.5)
    def send_message(self, msg:str):
        self.tracker.send_message(msg)
    def stop_recording(self):
        self.stop_ts = self.tracker.get_time_stamp()
        self.tracker.stop_recording()
    def save_data(self, file_stem: str, mon: monitors.Monitor) -> str:
        samples = np.array(self.tracker.buf.peek_time_range(self.start_ts, self.stop_ts))
        # Remove pupil size values
        eye_tracked = self.tracker.settings.EYE_TRACKED.lower()
        is_binocular = eye_tracked=='both'
        to_delete = np.s_[3, 6] if is_binocular else np.s_[3]
        samples = np.delete(samples, to_delete, 1)
        messages= [msg for msg in self.tracker.msg_buffer if msg[0]>=self.start_ts and msg[0]<=self.stop_ts]

        # parse target info
        msgs = [(m[0],pm) for m in messages if (pm:=parse_target_msg(m[1])) is not None]

        # collect and organize data
        to_include = np.logical_and(samples[:, 0]>=msgs[0][0], samples[:, 0]<=msgs[-1][0])
        samples = samples[to_include, :]
        
        # output: [timestamp, left_x, left_y, right_x, right_y, target_id, tar_x, tar_y] (timestamp in ms)
        # get all et data, make relative to screen center
        px,py = mon.currentCalib['sizePix']
        eye_fac = (1,px/2),(1,py/2)
        n_eye = 2 if is_binocular else 1
        for i,(fac,sub) in enumerate(((1,0.),)+n_eye*eye_fac):
            samples[:, i] = samples[:, i] * fac - sub

        # turn into dataframe
        cols = ['timestamp']
        if is_binocular or eye_tracked=='left':
            cols += ['left_x', 'left_y']
        if is_binocular or eye_tracked=='right':
            cols += ['right_x', 'right_y']
        data = pd.DataFrame(samples,columns=cols).sort_values(by='timestamp')
        # add target ID
        n_samp = data.shape[0]
        target_ids_pos = np.full((n_samp,3),-1,'int32')
        ts = data['timestamp'].to_numpy()
        for ti in range(0,len(msgs),2):
            is_target = np.logical_and(ts>=msgs[ti][0], ts<=msgs[ti+1][0])
            target_ids_pos[is_target,:] = msgs[ti][1][2:]
        target_ids_pos[:,2] = -target_ids_pos[:,2]  # positive y direction should be downward, not upwards like only PsychoPy does
        data[['target_id','tar_x','tar_y']] = target_ids_pos

        # done, now store to file. First get filename
        file_name = get_filename(file_stem,'.tsv')
        # store to file
        data.to_csv(file_name, '\t', na_rep='nan', index=False, float_format='%.8f')

        return file_name

class SMIEyeTracker(EyeTrackerBase):
    def __init__(self, tracker: 'SMITE_ET.Connect'):
        super().__init__()
        self.tracker = tracker
        self.start_ts: int = None
        self.stop_ts: int = None
        self.msgs = []
    def start_recording(self):
        self.tracker.start_recording()
        self.tracker.start_buffer(sample_buffer_length=None)
        core.wait(0.5)
        self.start_ts = self.tracker.rawSMI.get_current_time_stamp()
    def send_message(self, msg:str):
        self.tracker.send_message(msg)
        self.msgs.append([self.tracker.rawSMI.get_current_time_stamp(), msg])
    def stop_recording(self):
        self.stop_ts = self.tracker.rawSMI.get_current_time_stamp()
        self.tracker.stop_buffer()
        self.tracker.stop_recording()
    def save_data(self, file_stem: str, mon: monitors.Monitor) -> str:
        samples = self.tracker.peek_buffer_data_time_range(self.start_ts, self.stop_ts)
        samples = np.array([[s.timestamp,s.leftEye.gazeX,s.leftEye.gazeY,s.rightEye.gazeX,s.rightEye.gazeY] for s in samples])
        # replace missing data (coded as 0,0) with nan
        is_miss = np.logical_and(samples[:,1]==0, samples[:,2]==0)
        samples[is_miss,1:3] = np.nan
        is_miss = np.logical_and(samples[:,3]==0, samples[:,4]==0)
        samples[is_miss,3:5] = np.nan
        # get messages
        messages= [msg for msg in self.msgs if msg[0]>=self.start_ts and msg[0]<=self.stop_ts]

        # parse target info (message timestamp us->ms)
        msgs = [(m[0],pm) for m in messages if (pm:=parse_target_msg(m[1])) is not None]

        # collect and organize data
        to_include = np.logical_and(samples[:, 0]>=msgs[0][0], samples[:, 0]<=msgs[-1][0])
        samples = samples[to_include, :]

        # output: [timestamp, left_x, left_y, right_x, right_y, target_id, tar_x, tar_y] (timestamp in ms)
        # get all et data, make relative to screen center
        px,py = mon.currentCalib['sizePix']
        eye_fac = (1,px/2),(1,py/2)
        n_eye = 2
        for i,(fac,sub) in enumerate(((0.001,0.),)+n_eye*eye_fac):
            samples[:, i] = samples[:, i] * fac - sub
        # also scale message timestamp us->ms)
        msgs = [(m[0]/1000,m[1]) for m in msgs]

        # turn into dataframe
        cols = ['timestamp','left_x','left_y','right_x','right_y']
        data = pd.DataFrame(samples,columns=cols).sort_values(by='timestamp')

        # add target ID and store
        return add_target_and_save_data(file_stem, data, msgs)
    
def add_target_and_save_data(file_stem: str, data: pd.DataFrame, msgs: list[tuple[int,list]]) -> str:
    # add target ID
    n_samp = data.shape[0]
    target_ids_pos = np.full((n_samp,3),-1,'int32')
    target_ids_pos[:,2] = 1 # sign will be flipped below
    ts = data['timestamp'].to_numpy()
    for ti in range(0,len(msgs),2):
        is_target = np.logical_and(ts>=msgs[ti][0], ts<=msgs[ti+1][0])
        target_ids_pos[is_target,:] = msgs[ti][1][2:]
    target_ids_pos[:,2] = -target_ids_pos[:,2]  # positive y direction should be downward, not upwards like only PsychoPy does
    data[['target_id','tar_x','tar_y']] = target_ids_pos

    # done, now store to file. First get filename
    file_name = get_filename(file_stem,'.tsv')
    # store to file
    data.to_csv(file_name, sep='\t', na_rep='nan', index=False, float_format='%.8f')

    return file_name

def parse_target_msg(message: str):
    # format: onset/offset target N (x,y)
    if not message.startswith(('onset','offset')):
        return None
    return [t(s) for t,s in zip((str,str,int,int,int),[f for e in message.split() for f in ([e] if not e.startswith('(') else e.strip('()').split(','))])]

def get_eye_tracker_wrapper(tracker):
    track_qual_name = f"{tracker.__module__}.{tracker.__class__.__qualname__}"
    
    if track_qual_name=="titta.Tobii.myTobii":
        return TobiiEyeTracker(tracker)
    elif track_qual_name=="eyelink.eyelink.Connect":
        return EyeLinkTracker(tracker)
    elif track_qual_name=="py_smite.SMITE_ET.Connect":
        return SMIEyeTracker(tracker)
    else:
        raise NotImplementedError(f'Support for a "{track_qual_name}" eye tracker is not implemented')

def get_filename(stem: str, ext: str):
    if stem is None:
        stem = 'ETDQualitizer_'+time.strftime("%Y%m%d-%H_%M_%S")

    files = list(pathlib.Path.cwd().glob(f'*{ext}'))

    i = 1
    file_add = ''
    while True:
        # Go through the files and look for a match
        filename_exists = False
        for f in files:
            # if the file exists
            if stem + file_add == f.stem:
                # append '_i to filename
                file_add = '_' + str(i)
                i += 1
                filename_exists = True
                break

        # If we've gone through all files without
        # a match, we ready!
        if not filename_exists:
            break

    # Add the new extension the the filename
    return stem+file_add+ext

class ABCFixPoint:
    def __init__(self, win, outer_diameter=14, inner_diameter=2,
                 outer_color='black', inner_color='white', units='pix'):
        self.outer_dot      = visual.Circle(win, fillColor = outer_color,
                                            radius = outer_diameter/2,
                                            units = units)
        self.inner_dot      = visual.Circle(win, fillColor = outer_color,
                                            radius = inner_diameter/2,
                                            units = units)
        self.line_vertical  = visual.Rect(win, fillColor=inner_color,
                                          width=inner_diameter, height=outer_diameter,
                                          units = units)
        self.line_horizontal= visual.Rect(win, fillColor=inner_color,
                                          width=outer_diameter, height=inner_diameter,
                                          units = units)

    def set_size(self, size: float):
        self.outer_dot.radius       = size / 2
        self.line_vertical.size     = (self.line_vertical.width, size)
        self.line_horizontal.size   = (size, self.line_horizontal.height)

    def set_pos(self, pos):
        self.outer_dot.pos          = pos
        self.inner_dot.pos          = pos
        self.line_vertical.pos      = pos
        self.line_horizontal.pos    = pos

    def get_pos(self):
        return self.outer_dot.pos

    def get_size(self):
        return self.outer_dot.size

    def draw(self):
        self.outer_dot.draw()
        self.line_vertical.draw()
        self.line_horizontal.draw()
        self.inner_dot.draw()

def read_coord_file(file):
    return pd.read_csv(file, dtype=defaultdict(lambda: np.float32, ID='int32', color='str')).dropna(axis=0, how='all').set_index('ID')

def check_escape(win: visual.Window):
    if 'escape' in event.getKeys():
        win.close()
        core.quit()


def prepare_validation(win: visual.Window, config: dict):
    # get markers and target positions
    file = pathlib.Path(config["targets"]["file"])
    if not file.is_absolute():
        file = pathlib.Path(__file__).parent.resolve() / file
    target_positions = read_coord_file(file)

    # Create a background with circle placeholders where target will appear
    stimList = []
    if (use_file_color:='color' in target_positions.columns) or config["targets"]["placeholder"]["color"]:
        # if color specified in the the target_positions file, use that. Else, fall back on placeholder color
        # specified in config. If neither is specified, do not draw placeholders
        for _, target in target_positions.iterrows():
            circle = visual.Circle(win, radius=config["targets"]["placeholder"]["diameter"]/2, fillColor=target.color if use_file_color else config["targets"]["placeholder"]["color"])
            circle.pos = [target.x, target.y]
            stimList.append(circle)
            circle.draw()

    # Get screenshot of background, so that we can draw unchanging things at once
    background = visual.ImageStim(win, visual.BufferImageStim(win, stim=stimList).image)    # because https://github.com/psychopy/psychopy/issues/840
    return {'background': background, 'target_positions': target_positions}

def show_validation(win: visual.Window, config: dict, refresh_rate: int, task_vars: dict, tracker: EyeTrackerBase):
    # prepare visual objects
    fixation_target = ABCFixPoint(win, outer_color=config["targets"]["look"]["outer_color"], inner_color=config["targets"]["look"]["inner_color"])
    if (have_cue:=config["targets"]["cue"]["color"] is not None):
        cue = visual.Circle(win, radius=config["targets"]["cue"]["diameter"]/2, fillColor=config["targets"]["cue"]["color"])

    # prepare target order
    targets = [t for t in task_vars['target_positions'].index.to_list()]
    if config["targets"]["randomize_order"]:
        random.shuffle(targets)

    # prepare task parameters
    n_shrink_frames = int(config["targets"]["shrink"]["duration"]*refresh_rate)
    shrink_sizes    = np.linspace(config["targets"]["look"]["diameter_max"], config["targets"]["look"]["diameter_min"], n_shrink_frames)
    n_capture_frames= int(config["targets"]["duration"]*refresh_rate)

    # show fixation target sequence
    old_pos = task_vars['target_positions'].loc[targets[0],['x','y']].to_numpy()
    for t_id in targets:
        check_escape(win)

        # Move target to new position
        pos   = task_vars['target_positions'].loc[t_id,['x','y']].tolist()
        dist  = np.sqrt((old_pos[0]-pos[0])**2 + (old_pos[1]-pos[1])**2)
        # Target should move at constant speed regardless of distance to cover, duration contains time to move
        # over width of whole screen. Adjust time to proportion of screen width covered by current move
        move_duration = max(config["targets"]["move"]["min_duration"], config["targets"]["move"]["duration"]*dist/win.size[0])
        n_move_frames = int(move_duration*refresh_rate)
        if config["targets"]["move"]["move_with_acceleration"]:
            accel   = 0 if not dist else dist/(move_duration/2)**2    # solve x=.5*a*t^2 for a, use d/2 for x
            moveVec = [0.,0.] if not dist else [(x-y)/dist for x,y in zip(pos,old_pos)]
            def calc_pos(frac):
                if frac<.5:
                    return [p+m*.5*accel*(frac*move_duration)**2 for p,m in zip(old_pos,moveVec)]
                else:
                    # implement deceleration by accelerating from the other side in backward time
                    return [p-m*.5*accel*((1-frac)*move_duration)**2 for p,m in zip(pos,moveVec)]
            tar_pos = [calc_pos(x) for x in np.linspace(0., 1., n_move_frames)]
        else:
            x_tar_pos = np.linspace(old_pos[0], pos[0], n_move_frames)
            y_tar_pos = np.linspace(old_pos[1], pos[1], n_move_frames)
            tar_pos = [[x,y] for x,y in zip(x_tar_pos,y_tar_pos)]

        # move to next target
        fixation_target.set_size(config["targets"]["move"]["diameter"])
        cue.pos = pos
        for p in tar_pos:
            task_vars['background'].draw()
            if have_cue:
                cue.draw()
            fixation_target.set_pos(p)
            fixation_target.draw()
            win.flip()

        # shrink
        for s in shrink_sizes:
            task_vars['background'].draw()
            fixation_target.set_size(s)
            fixation_target.draw()
            win.flip()

        # fixation period
        msg = f'target {t_id} ({pos[0]:.0f},{pos[1]:.0f})'
        tracker.send_message('onset '+msg)
        for _ in range(n_capture_frames):
            task_vars['background'].draw()
            fixation_target.draw()
            win.flip()
        tracker.send_message('offset '+msg)

        old_pos = pos

def run_validation(win: visual.Window, config: dict, tracker) -> str:
    # get eye tracker wrapper
    tracker = get_eye_tracker_wrapper(tracker)
    # prepare instruction
    textstim = visual.TextStim(win, text="", height=config["instruction_text"]["height"], color=config["instruction_text"]["color"], wrapWidth=9999.)

    # validation instruction
    textstim.text = 'Follow the fixation target around, and keep fixating the center of it\n\n(Press the spacebar to start)'
    textstim.draw()
    win.flip()
    event.waitKeys()

    # prepare validation
    task_vars = prepare_validation(win, config["validation"])

    # run validation
    tracker.start_recording()
    for _ in range(config["validation"]["n_repetitions"]):
        show_validation(win, config["validation"], config["screen"]["refresh_rate"], task_vars, tracker)
    tracker.stop_recording()

    # save the validation data in our common format
    return tracker.save_data(None, win.monitor)


def open_demo_screen(config: dict) -> visual.Window:
    # Open window, check
    mon = monitors.Monitor('temp')
    mon.setWidth(config["screen"]["width"])
    mon.setDistance(config["screen"]["viewing_distance"])
    mon.setSizePix(config["screen"]["resolution"]) 
    win = visual.Window(monitor=mon,
                        fullscr=True,
                        color=config["screen"]["background_color"],
                        screen=config["screen"]["which_monitor"] or 0,
                        units='pix',
                        allowGUI=False,
                        multiSample=True,
                        numSamples=4)
    win.mouseVisible = False
    if not all((x==y for x,y in zip(win.size,config["screen"]["resolution"]))):
        raise RuntimeError(f'expected resolution of {config["screen"]["resolution"]}, but got {win.size}')
    if 1/win.monitorFramePeriod<config["screen"]["refresh_rate"]-2 or 1/win.monitorFramePeriod>config["screen"]["refresh_rate"]+2:  # check within 2 Hz
        raise RuntimeError(f'expected framerate of {config["screen"]["refresh_rate"]}, but got {1/win.monitorFramePeriod}')
    
    return win