from PyQt5.QtCore import QThread, pyqtSignal
import numpy as np
import time


class MichelsonMeasurement(QThread):
    change_linear_position = pyqtSignal(str)
    change_piezo_position = pyqtSignal(str)
    toggle_ui = pyqtSignal(bool)
    update_progress = pyqtSignal(int)
    update_label = pyqtSignal(tuple)
    toggle_button = pyqtSignal(bool)
    device_disconnect = pyqtSignal(str)
    toggle_trigger = pyqtSignal(bool)
    reset_frames = pyqtSignal(bool)

    def __init__(self, connections,
                 l_start: float, l_end: float, l_step: int,
                 p_start: float, p_end: float, p_step: int,equal_length: float):

        super(MichelsonMeasurement, self).__init__()

        self.com = connections

        self.l_start = l_start
        self.l_end = l_end
        self.l_step = l_step
        self.p_start = p_start
        self.p_end = p_end
        self.p_step = p_step
        self.equal_length = equal_length
        
        self.terminate = False


    def run(self) -> None:

        # turn off UI input elements while measuring
        self.toggle_ui.emit(False)
        self.toggle_button.emit(False)

        start = time.time()

        lpositions = np.linspace(self.l_start, self.l_end, self.l_step)
        ppositions = np.linspace(self.p_start, self.p_end, self.p_step)
        
        for l, lpos in enumerate(lpositions,1):
            try:
                self.com.request("michelson_linear", "set_vel",25)
                time.sleep(0.1)
                status, curr_pos_stage = self.com.request("michelson_linear", "set_pos",lpos)
                curr_pos_stage = float(np.round(curr_pos_stage,4))
                self.change_linear_position.emit(str(curr_pos_stage))
                time.sleep(0.1)
                self.com.request("michelson_linear", "set_vel",0)
            except Exception as e:
                print("Exception in setting linear stage.",e)
                self.terminate = True

            for p, ppos in enumerate(ppositions):
                p += 1
                try:
                    status, curr_pos_piezo = self.com.request("michelson_piezo", "set_pos", ppos)
                    curr_pos_piezo = float(np.round(curr_pos_piezo, 4))
                    self.change_piezo_position.emit(str(curr_pos_piezo))
                    time.sleep(0.1)
                except Exception as e:
                    print("Exception in setting piezo stage.",e)
                    self.terminate = True

                try:
                    final = (l) == len(lpositions) and p == len(ppositions)
                    
                    if final:
                        self.com.request("michelson_linear", "trigger_spectrometer", False)
                        time.sleep(0.1)
                    
                    else:
                        self.com.request("michelson_linear", "trigger_spectrometer", True)                    
                        time.sleep(0.1)
                    
                except Exception as e:
                    print("Exception in setting linear stage trigger.",e)
                    self.device_disconnect.emit("michelson_linear")
                    self.terminate = True

                elapsed = time.time() - start
                remaining = elapsed / max(l * len(ppositions) + p, 1) * len(lpositions) * len(ppositions) - elapsed

                self.update_label.emit((elapsed, remaining))

                self.update_progress.emit(int((l * len(ppositions) + p) / len(lpositions) / len(ppositions) * 100))
            
                if self.terminate:
                    self.com.request("michelson_linear", "set_vel",25)
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
                print("linear stage =  {:.4}, piezo stage= {:.4}".format(lpos,ppos))
                
        self.com.request("michelson_linear", "trigger_spectrometer", True)   
        time.sleep(0.1)


        self.com.request("michelson_linear", "set_vel",25)
        time.sleep(0.1) 
        self.update_progress.emit(100)
        time.sleep(0.1)
        self.toggle_ui.emit(True)
        time.sleep(0.1)
        self.toggle_button.emit(True)
        time.sleep(0.1)
        self.reset_frames.emit(True)
        time.sleep(0.1)
        self.toggle_trigger.emit(False)
        time.sleep(0.1)

