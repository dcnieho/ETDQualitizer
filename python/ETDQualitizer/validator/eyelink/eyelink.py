# Pylink wrapper for Psychopy
import time
import threading
import numpy as np
import pylink
from . import psychocal
from psychopy.tools.monitorunittools import deg2pix, pix2deg
from psychopy import event, core
from collections import deque
import copy
import os

class Settings(object):

    def __init__(self):
        ''' Default settings for eye tracker
        '''

        self.MY_MONITOR                  = 'testMonitor' # needs to exists in PsychoPy monitor center
        
        self.FILENAME                    = 'test'
        self.FILEPATH                    = ''
        self.SAMPLE_RATE                 = 1000  # hertz # Always use 1000 - lower rates are filtered/downsampled versions of 1000 Hz data
        
        self.N_CAL_TARGETS               = 9 #  Use 9 point calib. 13 point is for widescreens
        self.PACING_INTERVAL             = 1000 # in ms
        self.CALIBRATION_CORNER_SCALING  = 1
        self.VALIDATION_CORNER_SCALING   = 0.8
        self.CALIBRATION_AREA_PROPORTION = [0.9, 0.9]
        self.VALIDATION_AREA_PROPORTION  = [0.9, 0.9]
        self.CALIBRATION_TARGET_COLOR    = -1
        
        
        self.SCREEN_RES                  = [1920,1080]
        self.SCREEN_WIDTH                = 53 # cm
        self.SCREEN_HEIGHT               = 30 # cm
        self.VIEWING_DIST                = 80  # distance from eye to center of screen (cm)
        self.VIEWING_DIST_TOP_BOTTOM     = None #[830, 810] # distance to top of screen, distance to bottom of screen (mm, use if you want valid parser output)
        
        self.PUPIL_TRACKING_MODE         = 'CENTROID' # CENTROID OR ELLIPSE
        self.PUPIL_SIZE_MODE             = 'AREA' # 'AREA' OR 'DIAMETER'
        self.HEURISTIC_FILTER            = [0, 0]  # Link and file (levels: 0 - off, 1 - normal, 2 - extra) TODO, set different sensitivities for link and file
                                             # [1,2]default EyeLink II or EyeLink 1000 Heuristic filter settings]                                        
        self.SET_HEURISTIC_FILTER        = True  # Activate filter or not (need to be set every time the recording is started!)
        
                                               
        self.FILE_EVENT_FILTER           = 'LEFT,RIGHT,MESSAGE,BUTTON,INPUT'
        self.LINK_EVENT_FILTER           = 'LEFT,RIGHT,FIXATION,SACCADE,BLINK,MESSAGE,BUTTON,INPUT'
        self.LINK_SAMPLE_DATA            = 'LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS,HTARGET'
        self.FILE_SAMPLE_DATA            = 'LEFT,RIGHT,GAZE,GAZERES,AREA,HREF,PUPIL,STATUS,INPUT,HMARKER,HTARGET'
        
        self.ENABLE_SEARCH_LIMITS        = 'OFF' # ON (default) or off
        self.ILLUMINATION_POWER          = 2  # 'elcl_tt_power 2',2 is 75%, 1-100%, 3, 50%
        self.HOST_IP                     = '100.1.1.1' 
        
        
        # Can you set the physical setup?
        self.EL_CONFIGURATION            = 'BTABLER' # 'elcl_select_configuration, options [MTABLER, BTABLER, RTABLER, RBTABLER, AMTABLER, ARTABLER, BTOWER]
        self.EYE_TRACKED                 = 'Both' #'Both' #  BOTH/LEFT/RIGHT. binocular_enabled = YES to track both eyes. active_eye = RIGHT/LEFT to choose eye in monocular mode
        
        self.PATH2CONVERTER = '"C:\\Program Files (x86)\\SR Research\\EyeLink\\EDF_Access_API\\Example\\edf2asc.exe"'
        self.CONVERSION_OPTIONS = '-sg -miss nan -y'        

class Connect(object):
    """
    Provides functions for interacting with the EyeLink via Pylink.

    :param window: Psychopy window object.
    :param edfname: Desired name of the EDF file.
    :param record_raw_data: Set to True is you want to record pupil/CR data
    """

    def __init__(self, settings, record_raw_data=False, 
                   sample_buffer_length=0, 
                   use_sample_buffer=False,
                   read_from_eyelink_buffer=True,
                   event_buffer_length=0):
        """ Init
        Args:
            window: PsychoPy window used to draw in
            edfname: name of edf file
            record_raw_data: active recording of pupil and CR (raw) data
            sample_buffer_length: Store samples in buffer if you need more than the latest
            use_sample_buffer: store data from buffer in RingBuffer
            event_buffer_length: Store events in buffer if you need more than the latest
            use_eyelink_buffer: use getNextData() instead of getNewestSample() 
                the former draws from an internal buffer and should thus miss fewer samples
        """
        
        self.settings = settings
        
        # Init variables 
        self.sres = settings.SCREEN_RES
        self.record_raw_data = record_raw_data
        self.sample_buffer_length = sample_buffer_length
        self.read_from_eyelink_buffer = read_from_eyelink_buffer
        self.use_sample_buffer = use_sample_buffer

        self.msg_buffer = []
        
        # Make pylink accessible
        self.pylink = pylink
        
        # Initiate raw sample buffer (default length 200 samples)
        if self.sample_buffer_length != 0:
            self.buf = RingBuffer(maxlen=sample_buffer_length)
            
        # Initiate event buffer (default length 200 samples)
        if event_buffer_length != 0:
            self.fixdur_buf = RingBuffer(maxlen=event_buffer_length)           
            self.blinkdur_buf = RingBuffer(maxlen=event_buffer_length)        
            self.pupsize_buf = RingBuffer(maxlen=event_buffer_length)            
        
        # Make filename
        self.edfname = settings.FILENAME + '.edf'

        # Initialize connection with eye-tracker
        try:
            self.tracker = pylink.EyeLink(trackeraddress=settings.HOST_IP)
            self.realconnect = True
        except:
            self.tracker = pylink.EyeLink(None)
            self.realconnect = False
            
        # Double check that the hardware and software configruation agree
        print(settings.EL_CONFIGURATION)
        if self.realconnect:
            # Stop trackign if tracker is running
            self.stop_recording()
            
            self.tracker.readRequest("elcl_select_configuration") #  Read config
            time.sleep(.5) #  Wait for answer      
            reply = self.tracker.readReply()
            if settings.EL_CONFIGURATION not in reply:
                print(reply)
                print("Hardward and software config. mismatch")
                self.send_command("elcl_select_configuration = "+ settings.EL_CONFIGURATION)
                time.sleep(1)
                
        # Which eye should be tracked?
        self.select_eye(eye_tracked=settings.EYE_TRACKED)
                
        # Open EDF
        self.tracker.openDataFile(self.edfname)
        pylink.flushGetkeyQueue()
        self.tracker.setOfflineMode()
        
        # Override default settings
        self.set_constants()
                 

        # Keep track of left and right eyes
        self.left_eye  = 0
        self.right_eye = 1
        self.binocular = 2
        self.t_old     = 0  # Init timestamp of sample
        
        # Create self as a thread and initialize it
        self.samplethread = threading.Thread(target = self.sample_thread)
        self.__samplestop = True        
        
        # Create self as a thread and initialize it
        self.thread = threading.Thread(target = self.raw_thread)
        self.__stop = True
        
        # Create self as a thread and initialize it
        self.eventthread = threading.Thread(target = self.event_thread)
        self.__event_stop = True    
                
        # Setup the EDF-file such that it adds 'raw' data
        if self.record_raw_data:
            self.send_command("file_sample_raw_pcr = 0") # Don't write raw data to file... 
            self.send_command("link_sample_raw_pcr = 1" )  # only over link
            self.send_command("raw_pcr_dual_corneal = 1")# Enable tracking of two CR

            self.send_command("inputword_is_window = ON" )
            self.send_command("file_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,HREF,PUPIL,STATUS,INPUT,HMARKER,HTARGET" )
            self.send_command("link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,HREF,PUPIL,STATUS,INPUT,HMARKER,HTARGET" )
        else:
            self.send_command("file_sample_raw_pcr = 0")
            self.send_command("link_sample_raw_pcr = 0" )
            self.send_command("raw_pcr_dual_corneal = 0")

            self.send_command("file_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,HREF,PUPIL,STATUS,INPUT,HMARKER,HTARGET" )
            self.send_command("link_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,HREF,PUPIL,STATUS,INPUT,HMARKER,HTARGET" )

    def select_eye(self, eye_tracked='both'):
        '''
        Select eye to track
        eye - 'both', 'left', or 'right'
        '''
        
        if 'BOTH' in eye_tracked.upper():
            self.send_command("binocular_enabled = YES")
        else:
            self.send_command("binocular_enabled = NO")
            self.send_command("active_eye = "+ eye_tracked.upper())
            
    def calibrate(self, win, record_samples=False):
        """
        Calibrates eye-tracker using psychopy stimuli.

        :param cnum: Number of points to use for calibration. Options are 3, 5,
                     9, 13.
        :type cnum: int
        :param paval: Pacing of calibration, i.e. how long you have to fixate
                      each target.
        :type paval: int

        """

        # Generate custom calibration stimuli
        genv = psychocal.Psychocal(self.settings, self.tracker, win, self.settings.CALIBRATION_TARGET_COLOR)

        if self.realconnect:
            # Set calibration type
            calst = 'HV{}'.format(self.settings.N_CAL_TARGETS)
            self.tracker.setCalibrationType(calst)

            # Set calibraiton pacing
            self.tracker.setAutoCalibrationPacing(self.settings.PACING_INTERVAL)

            # Execute custom calibration display
            pylink.openGraphicsEx(genv)
            
            # Record samples during calibration and validation and store in edffile
            if record_samples:
                if self.record_raw_data:
                    self.tracker.sendCommand('sticky_mode_data_enable DATA = 1 1 1 1')
                else:
                    self.tracker.sendCommand('sticky_mode_data_enable DATA = 1 1 0 0')

            # Calibrate
            self.tracker.doTrackerSetup(self.sres[0], self.sres[1])
            
            # Stop sending sampels
            if record_samples:
                 self.tracker.sendCommand('sticky_mode_data_enable')
                 self.tracker.sendCommand('set_idle_mode') 
                 time.sleep(0.1) # Wait to finish mode transition

                 #  sticky_mode_data_enable is only switched off when there is an actual
                 #  mode change. the above set_idle_mode is a no-op as the tracker is
                 #  already offline at that point. If sticky mode is not switched off
                 #  properly, we end up with junk samples in a small bit of the edf file,
                 #  overwriting part of the next trial
                 self.tracker.sendCommand('setup_menu_mode') 
                 time.sleep(0.1) # Wait to finish mode transition
                 self.tracker.sendCommand('set_idle_mode') 
                 time.sleep(0.1) # Wait to finish mode transition                 
        else:
            genv.dummynote()
            
    def pupil_only(self):
        ''' Set tracker in pupil only mode
        '''
            
        # Active the selection of pupil only mode
        self.tracker.sendCommand("force_corneal_reflection = OFF")      # Default OFF
        self.tracker.sendCommand("allow_pupil_without_cr = ON")         # overwritten in Pupil/CR mode
        self.tracker.sendCommand("elcl_hold_if_no_corneal = OFF")       # Default OFF
        self.tracker.sendCommand("elcl_search_if_no_corneal = OFF")     # Default OFF
        self.tracker.sendCommand("elcl_use_pcr_matching = OFF")         # Default ON
        
        # ..and select it!
        self.tracker.sendCommand("corneal_mode = NO")         # Default ON
         
        

    def set_constants(self):
        """
        Overrides values in final.ini to make sure that the proper settings are used
        Values are imported from constants.py
        """
            
        # Set illumination power
        self.tracker.sendCommand('elcl_tt_power '+str(self.settings.ILLUMINATION_POWER))
        
        # Set display coords for dataviewer
        disptxt = 'DISPLAY_COORDS 0 0 %d %d' % (self.sres[0]-1,self.sres[1]-1)
        self.tracker.sendMessage(disptxt)
        
        scrtxt = 'screen_pixel_coords 0 0 %d %d' % (self.sres[0]-1,self.sres[1]-1)
        self.tracker.sendCommand(scrtxt)     
        
        # Set geometry to be able to use parser output (left, top, right, bottom, in mm)
        scrtxt = 'screen_phys_coords %f %f %f %f' % (-self.settings.SCREEN_WIDTH / 2.0 * 10, self.settings.SCREEN_HEIGHT / 2.0 * 10, 
                                                       self.settings.SCREEN_WIDTH / 2.0 * 10, -self.settings.SCREEN_HEIGHT / 2.0 * 10)
        self.tracker.sendCommand(scrtxt)      
        
        if self.settings.VIEWING_DIST_TOP_BOTTOM:
            scrtxt = 'screen_distance %d %d' % (self.settings.VIEWING_DIST_TOP_BOTTOM[0], 
                                                self.settings.VIEWING_DIST_TOP_BOTTOM[1])
            self.tracker.sendCommand(scrtxt)
                             
        # Set content of edf file
        self.tracker.sendCommand('file_event_filter = ' + self.settings.FILE_EVENT_FILTER)
        self.tracker.sendCommand('link_event_filter = ' + self.settings.LINK_EVENT_FILTER)
        self.tracker.sendCommand('link_sample_data = ' + self.settings.LINK_SAMPLE_DATA)
        self.tracker.sendCommand('file_sample_data = ' + self.settings.FILE_SAMPLE_DATA)        
        
        self.tracker.sendCommand("screen_distance = %d %d" % (self.settings.VIEWING_DIST*10, 
                                                              self.settings.VIEWING_DIST*10))
        self.tracker.sendCommand('screen_phys_coords = %ld %ld %ld %ld' % 
                            (-self.settings.SCREEN_WIDTH/2.0*10,  self.settings.SCREEN_HEIGHT/2.0*10, 
                              self.settings.SCREEN_WIDTH/2.0*10, -self.settings.SCREEN_HEIGHT/2.0*10));
        self.tracker.sendCommand(' '.join(['sample_rate','=',str(self.settings.SAMPLE_RATE)])) 
        self.tracker.sendCommand(' '.join(['pupil_size_diameter','=',self.settings.PUPIL_SIZE_MODE])) 
        
        self.tracker.sendCommand(' '.join(['calibration_corner_scaling','=',str(self.settings.CALIBRATION_CORNER_SCALING)]))   
        self.tracker.sendCommand(' '.join(['validation_corner_scaling','=',str(self.settings.VALIDATION_CORNER_SCALING)]))    
        self.tracker.sendCommand(' '.join(['calibration_area_proportion','=',' '.join([str(i) for i in self.settings.CALIBRATION_AREA_PROPORTION])]))   
        self.tracker.sendCommand(' '.join(['validation_area_proportion','=',' '.join([str(i) for i in self.settings.VALIDATION_AREA_PROPORTION])])) 
        self.tracker.sendCommand('heuristic_filter %d %d' % (self.settings.HEURISTIC_FILTER[0],
                                                             self.settings.HEURISTIC_FILTER[1]))
        
        
        if 'CENTROID' in self.settings.PUPIL_TRACKING_MODE:
            self.tracker.sendCommand("use_ellipse_fitter = NO")   # Use centroid mode per defautl
        else:
            self.tracker.sendCommand("use_ellipse_fitter = YES") 

    def get_time_stamp(self):
        """
        Gets timestamp of latest sample
        """
        
        timestamp = np.nan
        if self.realconnect:
            sample = self.tracker.getNewestSample()
            timestamp = sample.getTime()

        return timestamp
        
                              
    def set_status_host(self, message):
        """
        Sets status message to on host's screen appear while recording.

        :param message: Text object to send, must be < 80 char
        :type message: str
        """
        msg = "record_status_message '{}'".format(message)
        self.tracker.sendCommand(msg)

    def set_trial_id(self, idval=1):
        """
        Sends message that indicates start of trial in EDF.

        :param idval: Values to set TRIALID.
        """

        tid = 'TRIALID {}'.format(idval)
        self.tracker.sendMessage(tid)
        
    def set_trial_result(self, rval=0, scrcol=0):
        """
        Sends trial result to indiciate trial end in EDF and clears screen on
        EyeLink Display.

        :param rval: Value to set for TRIAL_RESULT.
        :type rval: float, str, or int
        :param scrcol: Color to clear screen to. Defaults to black.
        :type scrol: int
        """

        trmsg = 'TRIAL_RESULT {}'.format(rval)
        cscmd = 'clear_screen {}'.format(scrcol)

        self.tracker.sendMessage(trmsg)
        self.tracker.sendCommand(cscmd)        

    def start_recording(self, sendlink=False):
        """
        Starts recording. Waits 50ms to allow eyelink to prepare.

        :param sendlink: Toggle for sending eye data over the link to the
                         display computer during recording.
        :type sendlink: bool
        
        'heuristic_filter 0 0';     % switch off filters
        in eyelink so we get the raw data. Need to set this each time
        recording is started as according to manual it is reset at recording
        stop! It's on per default
        The first parameter 
        """
        if self.settings.SET_HEURISTIC_FILTER:
            cstr = 'heuristic_filter %d %d' % (self.settings.HEURISTIC_FILTER[0],
                                               self.settings.HEURISTIC_FILTER[1])
            self.tracker.sendCommand(cstr)

        self.tracker.sendCommand('set_idle_mode')
        time.sleep(.05)
        
        if self.record_raw_data:
            self.enable_raw(do_enable = True)
            sendlink = True
            
        if sendlink:
            self.tracker.startRecording(1, 1, 1, 1)
        else:
            self.tracker.startRecording(1, 1, 0, 0)
            
        if self.record_raw_data:
            self.enable_realtime_mode()
            self.start_raw()   
            
        if self.use_sample_buffer:
            self.start_sample()

    def stop_recording(self):
        """
        Stops recording.
        """
        if self.record_raw_data:
            self.stop_raw() # Stop thread
            self.disable_realtime_mode()
            
        self.stop_sample()
        self.stop_event_buffer()            
        
        self.tracker.stopRecording()
        
    def el_pix2deg(self, x, y, mon):
        """
        converts from eyelink pixels coords to psychopy deg.
        mon - psychopy monitor
        """

        elx = pix2deg(x-self.sres[0] / 2.0, mon) 
        ely = -(pix2deg(y - (self.sres[1] / 2.0), mon))
        
        return elx,ely

    def draw_AOI_host(self, x, y, size, index, color, name, mon):
        """
        Draws square interest area in EDF and a corresponding filled box on
        eye-tracker display.

        :param x: X coordinate in degrees visual angle for center of check area.
        :type x: float or int
        :param y: Y coordinate in degrees visual angle for center of check area.
        :type y: float or int
        :param size: length of one edge of square in degrees visual angle.
        :type size: float or int
        :param index: number to assign interest area in EDF
        :type index: int
        :param color: color of box drawn on eye-tracker display (0 - 15)
        :type color: int
        :param name: Name interest area in EDF
        :type name: str
        mon - psychopy monitor
        """

        # Convert units to eyelink space
        elx = deg2pix(x, mon) + (self.sres[0] / 2.0)
        ely = -(deg2pix(y, mon) - (self.sres[1] / 2.0))
        elsz = deg2pix(size, mon) / 2.0

        # Make top left / bottom right coordinates for square
        tplf = map(round, [elx - elsz, ely - elsz])
        btrh = map(round, [elx + elsz, ely + elsz])

        # Construct command strings
        flist = [index, name, color] + tplf + btrh
        iamsg = '!V IAREA RECTANGLE {0} {3} {4} {5} {6} {1}'.format(*flist)
        bxmsg = 'draw_filled_box {3} {4} {5} {6} {2}'.format(*flist)

        # Send commands
        self.tracker.sendMessage(iamsg)
        self.tracker.sendCommand(bxmsg)

    def send_var(self, name, value):
        """
        Sends a trial variable to the EDF file.

        :param name: Name of variable.
        :type name: str
        :param value: Value of variable.
        :type value: float, str, or int
        """

        # Make string
        varmsg = '!V TRIAL_VAR {} {}'.format(name, value)

        # Send message
        self.tracker.sendMessage(varmsg)


    def end_experiment(self, spath):
        """
        Closes and transfers the EDF file.
        WARNING - don't retreive a file using PsychoPy. Start exp program from cmd
        Otherwise file transfer can be very slow.

        :param spath: File path of where to save EDF file. Include trailing
                      slash.
        :type spath: str
        """

        # File transfer and cleanup!
        self.tracker.setOfflineMode()
        time.sleep(.5)

        # Generate file path
        fpath = spath + self.edfname

        # Close the file and transfer it to Display PC
        self.tracker.closeDataFile()
        time.sleep(1)
        print(self.edfname, fpath)
        self.tracker.receiveDataFile(self.edfname, fpath)
        self.tracker.close()
        
    def enable_realtime_mode(self):
        '''
        Enables EyeLink realtime mode
        '''
        self.pylink.beginRealTimeMode(100)
        
    def disable_realtime_mode(self):
        '''
        Disables EyeLink realtime mode
        '''
        self.pylink.endRealTimeMode()  
        
    def get_event_from_buffer(self):
        """
        Get the latest blink or fixation event over the link
        The link must have been activated first
        self.tracker.startRecording(0, 0, 1, 1)
        If both eyes are used, the left one is chosen by default
        
        0 - left eye
        1 - right eye
        2 - binocular
        """
        
        # Use these values if nothing else if produced
        event_type = None
        event_prop = []
        
        if self.realconnect:  # only start check loop if real connection
            
            # Keep getting samples until a new sample is found
            timeout = 0.010 # Don't wait for a new samples longer than 10 ms
            t0 = time.time()
            while (time.time() - t0) < timeout:
                data_type = self.tracker.getNextData()
                if data_type == 4: # 4 means ENDBLINK
                    event_type = 'blink'
                    blinkEvent = self.tracker.getFloatData() # Grab sample from buffer
                    event_prop.append(blinkEvent.getEndTime() - blinkEvent.getStartTime())
                elif data_type == 8: #8 means ENDFIX
                    event_type = 'fixation'
                    fixEvent = self.tracker.getFloatData() # Grab sample from buffer
                    event_prop.append(fixEvent.getEndTime() - fixEvent.getStartTime())
                    event_prop.append(fixEvent.getAveragePupilSize())
                elif data_type ==0x3F or data_type == 0:
                    break
                else:
                    continue    
        return event_type, event_prop 
                    
    def start_event_buffer(self):
        
        # First clear buffer from old stuff
        self.fixdur_buf.clear()
        self.pupsize_buf.clear()
        self.blinkdur_buf.clear()
        
        # Starts the thread 'event_thread'
        self.__event_stop = False
        self.eventthread.start()
        
    def event_thread(self):
        # Called by the e.g., eyetracker.event_thread()
        # Continously read events into the ringbuffer 
        
        k = 0
        while True:
            if self.__event_stop:        
                break        
            
            event_type, event_prop  = self.get_event_from_buffer()
            
            # If an actual event happened, write to buffer
            if event_type:
                # Take different action depending of what event it is
                if 'blink' in event_type:
                    self.blinkdur_buf.append(event_prop[0])
                elif 'fixation' in event_type:
                    self.fixdur_buf.append(event_prop[0])
                    self.pupsize_buf.append(event_prop[1])
                    
            if np.mod(k, 10) == 0:
                time.sleep(0.001)
            k += 1
            
    # Stop event thread
    def stop_event_buffer(self):
        self.__event_stop = True   # the run() method monitors this flag         
        
        # Start a new thread (each thread can be started only once)
        self.eventthread = threading.Thread(target = self.event_thread)
        
    def get_sample(self, write_to_edf = False):
        """
        Get the latest gaze sample over the link
        The link must have been activated first
        self.tracker.startRecording(0, 0, 1, 1)
        If both eyes are used, the left one is chosen by default
        
        0 - left eye
        1 - right eye
        2 - binocular
        
        eye_desired = 2   returns len(sample_info) = 4 [lx,ly,rx,ry]
        eye_desired = 0|1 returns len(sample_info) = 2 [Xx,Xy]
        
        write_to_msg = True mans that gaze data are written to the edf-file as messages
        This functionality is implemented when you want to get both raw pupil, cr, 
        and gaze samples; everything cannot be stored in the buffer.
        """
        
        # Use these values if nothing else if produced
        t = -1
        sample_info = None
        
        # Check which eye is being recorded
        eye_used = self.tracker.eyeAvailable()
        if self.realconnect:  # only start check loop if real connection
            # Keep getting samples until a new sample is found
            timeout = 0.010 # Don't wait for a new samples longer than 10 ms
            t0 = time.time()
            while (time.time() - t0) < timeout:
                sample = self.tracker.getNewestSample() # Grab latest sample      
                t = sample.getTime()
                if t != self.t_old and sample != None: 
                    sample_info = self.unpack_sample(t, sample, eye_used, write_to_edf)        
                    break
                        
                self.t_old = t

        return t, sample_info  
        
    def get_sample_from_buffer(self,write_to_edf = False):
        """
        Get the latest gaze sample over the link
        The link must have been activated first
        self.tracker.startRecording(0, 0, 1, 1)
        If both eyes are used, the left one is chosen by default
        
        0 - left eye
        1 - right eye
        2 - binocular
        
        eye_desired = 2   returns len(sample_info) = 4 [lx,ly,rx,ry, pr, pr]
        eye_desired = 0|1 returns len(sample_info) = 2 [Xx,Xy,pX]
        
        write_to_msg = True mans that gaze data are written to the edf-file as messages
        This functionality is implemented when you want to get both raw pupil, cr, 
        and gaze samples; everything cannot be stored in the buffer.
        """
        
        # Use these values if nothing else if produced
        t = -1
        sample_info = None
        
        # Check which eye is being recorded
        eye_used = self.tracker.eyeAvailable()
        if self.realconnect:  # only start check loop if real connection
            # Keep getting samples until a new sample is found
            timeout = 0.010 # Don't wait for a new samples longer than 10 ms
            t0 = time.time()
            while (time.time() - t0) < timeout:
                data_type = self.tracker.getNextData() 
                if data_type == 200: # 200 means it a sample, not an event
                    sample = self.tracker.getFloatData() # Grab sample from buffer
                    t = sample.getTime()
                    if sample != None: 
                        sample_info = self.unpack_sample(t, sample, eye_used, write_to_edf)        
                        break
                elif data_type ==0x3F or data_type == 0:
                    break
                else:
                    continue

        return t, sample_info        
        
    def unpack_sample(self, t, sample, eye_used, write_to_edf):
        """ converts sample into message string """        
        # Extract gaze positions and write as msg
        if eye_used == self.right_eye and sample.isRightSample():
            gaze_position =  sample.getRightEye().getGaze()
            pupil_size = sample.getRightEye().getPupilSize()
            sample_info = list(gaze_position) + [pupil_size]
            if write_to_edf:
                sample_str = ' '.join([str(t),str(gaze_position[0]),
                                    str(gaze_position[1]), str(pupil_size)])
        elif eye_used == self.left_eye and sample.isLeftSample():
            gaze_position = sample.getLeftEye().getGaze()
            pupil_size = sample.getLeftEye().getPupilSize()
            sample_info = list(gaze_position) + [pupil_size]
            if write_to_edf:
                sample_str = ' '.join([str(t),str(gaze_position[0]),
                                    str(gaze_position[1]), str(pupil_size)])
                 
        elif eye_used  == self.binocular and sample.isBinocular():
            r = sample.getRightEye().getGaze()
            pr = sample.getRightEye().getPupilSize()
            l = sample.getLeftEye().getGaze()
            pl = sample.getLeftEye().getPupilSize()        
            sample_info = [l[0],l[1], pl, r[0],r[1], pr]
            if write_to_edf:
                sample_str = ' '.join([str(t), str(sample_info[0]),
                                                str(sample_info[1]),
                                                str(sample_info[2]), 
                                                str(sample_info[3]),
                                                str(sample_info[4]),
                                                str(sample_info[5])])
        # Write to edf?
        if write_to_edf:
            self.send_message(sample_str)      
        
        # Add data to buffer if activated
        if self.sample_buffer_length != 0:
            self.buf.append([t] + sample_info) 
            
        return sample_info
        
    def enable_raw(self,do_enable = True):
        '''
        Enable/disable raw pupil and CR in online sample data over link
        '''
        #switch tracker to idle and give it time to complete mode switch
        self.tracker.setOfflineMode()
        time.sleep(0.050)
        pylink.enablePCRSample(do_enable)
        
    def get_raw_sample(self,write_to_edf = False):
        """
        Get the latest raw sample over the link
        The link must have been activated first
        self.tracker.startRecording(0, 0, 1, 1)
        If both eyes are used, the left one is chosen by default
        
        0 - left eye
        1 - right eye
        2 - binocular
        
        eye_desired = 2   returns len(raw) = 2x
        eye_desired = 0|1 returns len(raw) = x
        
        write_to_msg = True mans that gaze data are written to the edf-file as messages
        This functionality is implemented when you want to get both raw pupil, cr, 
        and gaze samples; everything cannot be stored in the buffer.
        """
        
        # Use these values if nothing else if produced
        t = -1
        raw = None
        
        # Check which eye is being recorded
        eye_used = self.tracker.eyeAvailable()
        if self.realconnect:  # only start check loop if real connection
            
            # Keep getting samples until a new sample is found
            timeout = 0.010 # Don't wait for a new samples longer than 10 ms
            t0 = time.time()
            while (time.time() - t0) < timeout:
                rawsample = self.tracker.getNewestSample() # Grab latest sample
                t = rawsample.getRawSampleTime()
                               
                # if the timestamps between old and new samples differ, it's new
                if t != self.t_old and rawsample != None: 
                    raw = self.unpack_raw_sample(t, rawsample, eye_used, write_to_edf)
                    break
                
                self.t_old = t
                
        return t, raw
        
    def get_raw_sample_from_buffer(self,write_to_edf = False):
        """
        Get the latest raw sample over the link
        The link must have been activsated first
        self.tracker.startRecording(0, 0, 1, 1)
        If both eyes are used, the left one is chosen by default
        
        0 - left eye
        1 - right eye
        2 - binocular
        
        eye_desired = 2   returns len(raw) = 2x
        eye_desired = 0|1 returns len(raw) = x
        
        write_to_msg = True mans that gaze data are written to the edf-file as messages
        This functionality is implemented when you want to get both raw pupil, cr, 
        and gaze samples; everything cannot be stored in the buffer.
        """
        
        # Use these values if nothing else if produced
        t = -1
        raw = None
        
        # Check which eye is being recorded
        eye_used = self.tracker.eyeAvailable()
        if self.realconnect:  # only start check loop if real connection
            # Keep getting samples until a new sample is found
            timeout = 0.010 # Don't wait for a new samples longer than 10 ms
            t0 = time.time()
            while (time.time() - t0) < timeout:
                data_type =  self.tracker.getNextData()
                if data_type == 200:
                    rawsample = self.tracker.getFloatData()

                    # if the timestamps between old and new samples differ, it's new
                    if t != self.t_old and rawsample != None:
                        t = rawsample.getRawSampleTime()
                        raw = self.unpack_raw_sample(t, rawsample, eye_used, write_to_edf)
                        self.t_old = t
                        break
                        
                elif data_type==0x3F or data_type==0: 
                    break
                else:
                    continue
                
        return t, raw        
        
    def unpack_raw_sample(self, t, rawsample, eye_used, write_to_edf):
        """ converts raw sample into message string """
        
        # Extract gaze positions and write as msg
        if eye_used == self.right_eye:
            
            raw = [t,
                   rawsample.getRightrRawPupil()[0], rawsample.getRightrRawPupil()[1],
                   rawsample.getRightPupilArea(),
                   rawsample.getRightPupilDimension()[0], rawsample.getRightPupilDimension()[1],
                   rawsample.getRightRawCr()[0], rawsample.getRightRawCr()[1],
                   rawsample.getRightCrArea(),
                   rawsample.getRightRawCr2()[0], rawsample.getRightRawCr2()[1],
                   rawsample.getRightCrArea2()]                       

            msg = ' '.join(['R',' '.join([str(r) for r in raw])])
            if write_to_edf:
                self.send_message(msg)       
        
        elif eye_used == self.left_eye:
            raw = [t,
                   rawsample.getLeftrRawPupil()[0], rawsample.getLeftrRawPupil()[1],
                   rawsample.getLeftPupilArea(),
                   rawsample.getLeftPupilDimension()[0], rawsample.getLeftPupilDimension()[1],
                   rawsample.getLeftRawCr()[0], rawsample.getLeftRawCr()[1],
                   rawsample.getLeftCrArea(),
                   rawsample.getLeftRawCr2()[0], rawsample.getLeftRawCr2()[1],
                   rawsample.getLeftCrArea2()]                           

            msg = ' '.join(['L',' '.join([str(l) for l in raw])])
            if write_to_edf:
                self.send_message(msg)                               
            
        elif eye_used == self.binocular:                        
            rawL = [t,
                   rawsample.getLeftrRawPupil()[0], rawsample.getLeftrRawPupil()[1],
                   rawsample.getLeftPupilArea(),
                   rawsample.getLeftPupilDimension()[0], rawsample.getLeftPupilDimension()[1],
                   rawsample.getLeftRawCr()[0], rawsample.getLeftRawCr()[1],
                   rawsample.getLeftCrArea(),
                   rawsample.getLeftRawCr2()[0], rawsample.getLeftRawCr2()[1],
                   rawsample.getLeftCrArea2()]                              
            
            rawR = [rawsample.getRightrRawPupil()[0], rawsample.getRightrRawPupil()[1],
                   rawsample.getRightPupilArea(),
                   rawsample.getRightPupilDimension()[0], rawsample.getRightPupilDimension()[1],
                   rawsample.getRightRawCr()[0], rawsample.getRightRawCr()[1],
                   rawsample.getRightCrArea(),
                   rawsample.getRightRawCr2()[0], rawsample.getRightRawCr2()[1],
                   rawsample.getRightCrArea2()]       
            raw =  rawL + rawR
            msg = ' '.join(['L',' '.join([str(lraw) for lraw in rawL]),'R',' '.join([str(rraw) for rraw in rawR])])
            if write_to_edf:
                self.send_message(msg)   

            # Add data to buffer if activated
            if self.sample_buffer_length != 0:
                self.buf.append(raw)
            
        else:
            print('Something went wrong...')
            
        return raw

    def start_sample(self):
        # Starts the thread 'sample_thread'
        self.__samplestop = False
        self.samplethread.start()
        print('sample thread started')
        
    def sample_thread(self):
        # Called by the e.g., eyetracker.start()
        # Continously read raw data into the ringbuffer (convert to deg)
        k = 0
        while True:
            if self.__samplestop:        
                break        
            if self.read_from_eyelink_buffer:
                self.get_sample_from_buffer(write_to_edf = False)
            else:
                self.get_sample(write_to_edf = False)
                
            if np.mod(k, 10) == 0:
                time.sleep(0.001)
            k += 1
        
    # Stop thread
    def stop_sample(self):
        self.__samplestop = True   # the run() method monitors this flag          
        
        # Start a new thread (each thread can be started only once)
        self.samplethread = threading.Thread(target = self.sample_thread)
        
    def start_raw(self):
        # Starts the thread 'raw_thread'
        self.__stop = False
        self.thread.start()
        
    def raw_thread(self):
        # Called by the e.g., eyetracker.start()
        # Continously read raw data into the ringbuffer (convert to deg)
        k = 0
        while True:
            if self.__stop:        
                break        
            if self.read_from_eyelink_buffer:
                self.get_raw_sample_from_buffer(write_to_edf = True)
            else:
                self.get_raw_sample(write_to_edf = True)
                
            if np.mod(k, 10) == 0:
                time.sleep(0.001)
            k += 1
        
    # Stop thread
    def stop_raw(self):
        self.__stop = True   # the run() method monitors this flag          
        
        # Start a new thread (each thread can be started only once)
        self.thread = threading.Thread(target = self.raw_thread)

    def fix_check(self, size, ftime, button, mon):
        """
        Checks that fixation is maintained for certain time.

        :param size: Length of one side of box in degrees visual angle.
        :type size: float or int
        :param ftime: Length of time to check for fixation in seconds.
        :type ftime: int
        :param button: Key to press to recalibrate eye-tracker.
        :type button: char
        mon - psychopy monitor
        """

        # Calculate Fix check borders
        cenX = self.sres[0] / 2.0
        cenY = self.sres[1] / 2.0
        size = deg2pix(size, mon) / 2.0

        xbdr = [cenX - size, cenX + size]
        ybdr = [cenY - size, cenY + size]

        # Set status message & Draw box
        self.set_status_host('Fixation Check')
        bxmsg = 'draw_box {} {} {} {} 1'.format(xbdr[0], ybdr[0], xbdr[1],
                                                ybdr[1])
        self.tracker.sendCommand(bxmsg)

        # Begin recording
        self.tracker.startRecording(0, 0, 1, 1)

        # Check which eye is being recorded
        eye_used = self.tracker.eyeAvailable()
        RIGHT_EYE = 1
        LEFT_EYE = 0

        # Begin polling
        keys = []
        t0 = time.time()
        while self.realconnect:  # only start check loop if real connection
            # Grab latest sample
            sample = self.tracker.getNewestSample()

            # Extract gaze coordinates
            if eye_used == RIGHT_EYE:
                gaze = sample.getRightEye().getGaze()
            else:
                gaze = sample.getLeftEye().getGaze()

            # Are we in the box?
            if xbdr[0] < gaze[0] < xbdr[1] and ybdr[0] < gaze[1] < ybdr[1]:
                # Have we been in the box long enough?
                if (time.time() - t0) > ftime:
                    self.tracker.stopRecording()
                    break
            else:
                # Reset clock if not in box
                t0 = time.time()
                
    def send_message(self, txt):
        """
        Sends a message to the tracker that is recorded in the EDF.

        :param txt: Message to send.
        :type txt: str
        """

        # Send message
        self.tracker.sendMessage(txt)

        # Add message to msg buffer
        self.msg_buffer.append([self.get_time_stamp(), txt])

    def send_command(self, cmd):
        """
        Sends a command to the Eyelink.

        :param cmd: Command to send.
        :type cmd: str
        """

        # Send Command
        self.tracker.sendCommand(cmd)

    def draw_text_host(self, msg):
        """
        Draws text on eye-tracker screen.

        :param msg: Text to draw.
        :type msg: str
        """

        # Figure out center
        x = self.sres[0] / 2

        # Send message
        txt = '"{}"'.format(msg)
        self.tracker.drawText(txt, (x, 50))
       
      
class RingBuffer(object):
    """ A simple ring buffer based on the deque class"""
    def __init__(self, maxlen=200):        
        # Create que with maxlen 
        self.maxlen = maxlen
        self._b = deque(maxlen=maxlen)  

    def clear(self):
        """ Clears buffer """
        return(self._b.clear())
        
    def get_all(self):
        """ Returns all samples from buffer and empties the buffer"""
        lenb = len(self._b)
        return(list([self._b.popleft() for i in range(lenb)]))
        
    def peek(self):
        """ Returns all samples from buffer without emptying the buffer
        First remove an element, then add it again
        """
        b_temp = copy.copy(self._b)
        c = []
        if len(b_temp) > 0:
            for i in range(len(b_temp)):
                c.append(b_temp.pop())

        return(c)

    def peek_time_range(self, t0, t1):
        """ Returns all samples from buffer without emptying the buffer
        First remove an element, then add it again
        """
        b_temp = copy.copy(self._b)
        c = []
        if len(b_temp) > 0:
            for i in range(len(b_temp)):
                sample = b_temp.pop()
                if (sample[0] >= t0) and (sample[0] <= t1):
                    c.append(sample)

        return(c)		
        
    def append(self, L):
        """"Append buffer with the most recent sample (list L)"""
        self._b.append(L)      