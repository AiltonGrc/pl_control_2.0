from PyQt5.QtCore import QThread, pyqtSignal
import numpy as np
import time


class Measurement_1D(QThread):
    device_string = pyqtSignal(str)
    change_device_position = pyqtSignal(float)
    toggle_ui = pyqtSignal(bool)
    update_progress = pyqtSignal(int)
    update_label = pyqtSignal(tuple)
    toggle_button = pyqtSignal(bool)
    device_disconnect = pyqtSignal(str)
    update_position = pyqtSignal(str)
    toggle_trigger = pyqtSignal(bool)
    reset_frames = pyqtSignal(bool)


    def __init__(self,connections,
                 l_start: float, l_end: float, l_step: int, device_name: str, offset: float,starting_pos: float=None):

        super(Measurement_1D, self).__init__()

        self.com = connections

        if starting_pos == None:
            self.starting_pos = l_start
        else:
            self.starting_pos = starting_pos
        self.l_start = l_start
        self.l_end = l_end
        self.l_step = int(l_step)
        self.device_name = device_name
        self.terminate = False	
        self.offset = offset
	
    def run(self) -> None:

        # turn off UI input elements while measuring
        self.toggle_ui.emit(False)
        self.toggle_button.emit(False)

        start = time.time()

        lpositions = np.linspace(self.l_start, self.l_end, self.l_step)
        

        for l, lpos in enumerate(lpositions,1):
            try:
                print(lpos)
                status, curr_pos = self.com.request(self.device_name, "set_pos",lpos)
            except:
                self.terminate = True

            try:
                final = (l) == len(lpositions)
                
                if final:
                    self.com.request("michelson_linear", "trigger_spectrometer", False)
                    
                else:
                    self.com.request("michelson_linear", "trigger_spectrometer", True) 
   
            except:
                self.terminate = True

            elapsed = time.time() - start

            remaining = elapsed / l * len(lpositions) - elapsed

            pos = np.round(float(curr_pos)-float(self.offset),4)

            self.update_label.emit((elapsed, remaining))
            self.update_position.emit(str(pos))

            self.update_progress.emit(int(l)/ len(lpositions) * 100)
            
            if self.terminate:
                status, curr_pos = self.com.request(self.device_name, "set_pos", self.starting_pos)
                time.sleep(0.1)
                self.update_position.emit(str(float(curr_pos) - float(self.offset)))
                time.sleep(0.1)
                self.update_progress.emit(100)
                time.sleep(0.1)
                self.reset_frames.emit(True)
                time.sleep(0.1)
                self.toggle_trigger.emit(False)
                time.sleep(0.1)
                self.toggle_ui.emit(True)
                time.sleep(0.1)
                self.toggle_button.emit(True)
                time.sleep(0.1)
                return
        
        status, curr_pos = self.com.request(self.device_name, "set_pos",self.starting_pos)
        time.sleep(0.1)
        self.update_position.emit(str(float(curr_pos)-float(self.offset)))
        time.sleep(0.1)

        self.update_progress.emit(100)
        time.sleep(0.1)
        self.toggle_ui.emit(True)
        time.sleep(0.1)
        self.toggle_button.emit(True)
        time.sleep(0.1)

        return

