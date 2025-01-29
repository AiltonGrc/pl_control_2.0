from __future__ import annotations

from PyQt5.QtCore import QThread, pyqtSignal
import numpy as np
import time


class Measurement_2D(QThread):
    device_string = pyqtSignal(str)
    change_long_position = pyqtSignal(str)
    change_fast_position = pyqtSignal(str)
    toggle_ui = pyqtSignal(bool)
    update_progress = pyqtSignal(int)
    update_label = pyqtSignal(tuple)
    toggle_button = pyqtSignal(bool)
    device_disconnect = pyqtSignal(str)
    update_pos_slow = pyqtSignal(object)
    update_pos_fast = pyqtSignal(object)
    toggle_trigger = pyqtSignal(bool)
    reset_frames = pyqtSignal(bool)

    def __init__(self, connections,
                 slow_start: float, slow_end: float, slow_step: int, offset_slow: float,
                 fast_start: float, fast_end: float, fast_step: float, offset_fast: float,
                 device_name_slow: str, device_name_fast: str):

        super(Measurement_2D, self).__init__()

        self.com = connections

        self.slow_start = slow_start
        self.slow_end = slow_end
        self.slow_step = int(slow_step)
        self.slow_offset = offset_slow
        self.fast_start = fast_start
        self.fast_end = fast_end
        self.fast_step = int(fast_step)
        self.fast_offset = offset_fast
        self.device_name_slow = device_name_slow
        self.device_name_fast = device_name_fast
        self.terminate = False

    def req_update_pos(self,device: str,pos):
        
        if 'kethley' in device.lower().strip():
            status, (curr_V, curr_I) = self.com.request(device, "apply_voltage", pos)
            curr_pos = (float(np.round(curr_V, 4)), float(np.round(curr_I, 9)))
        else:
            status, curr_pos = self.com.request(device, "set_pos", pos)
            curr_pos = float(np.round(curr_pos, 4))

        return status,curr_pos
    
    def emit_position(self,device:str,axis:str, pos):
        
        if axis == "fast":
            offset = self.fast_offset
        else:
            offset = self.slow_offset
        
        if "kethley" in device:
            pos_to_emit = (float(pos[0]),float(pos[1]))
            print(pos_to_emit)
        else:
            pos_to_emit = str(float(pos) - float(offset))

            
        if axis == "fast":
            self.update_pos_fast.emit(pos_to_emit)
        else:
            self.update_pos_slow.emit(pos_to_emit)
            

    def run(self) -> None:
        
        ppos = self.fast_start
        
        # turn off UI input elements while measuring
        self.toggle_ui.emit(False)
        self.toggle_button.emit(False)

        start = time.time()

        lpositions = np.linspace(self.slow_start, self.slow_end, self.slow_step)
        ppositions = np.linspace(self.fast_start, self.fast_end, self.fast_step)

        for l, lpos in enumerate(lpositions,1):
            try:
                status, curr_pos = self.req_update_pos(self.device_name_slow,lpos)
                self.emit_position(self.device_name_slow, "slow", curr_pos)
            except Exception as e:
                print("Error in slow axis with {}".format(self.device_name_slow),e)
                self.device_disconnect.emit(self.device_name_slow)
                self.terminate = True

            for p, ppos in enumerate(ppositions):
                p += 1
                try:
                    status, curr_pos2 = self.req_update_pos(self.device_name_fast,ppos)
                    self.emit_position(self.device_name_fast, "fast", curr_pos2)
                except Exception as e:
                    print("Error in fast axis with {}".format(self.device_name_fast),e)
                    self.device_disconnect.emit(self.device_name_fast)
                    self.terminate = True

                try:
                    final = (l) == len(lpositions) and p == len(ppositions)
                    
                    if final:
                        self.com.request("michelson_linear", "trigger_spectrometer", False)
                    
                    else:
                        self.com.request("michelson_linear", "trigger_spectrometer", True)                    
                    
                except:
                    self.device_disconnect.emit(self.device_name_fast)
                    self.terminate = True

                elapsed = time.time() - start
                remaining = elapsed / max(l * len(ppositions) + p, 1) * len(lpositions) * len(ppositions) - elapsed

                self.update_label.emit((elapsed, remaining))

                self.update_progress.emit(int((l * len(ppositions) + p) / len(lpositions) / len(ppositions) * 100))
                
                if self.terminate:
                    status, curr_pos = self.req_update_pos(self.device_name_fast,self.fast_start)
                    time.sleep(0.1)
                    self.emit_position(self.device_name_fast, "fast", curr_pos)
                    time.sleep(0.1)
                    status, curr_pos = self.req_update_pos(self.device_name_slow,self.slow_start)
                    time.sleep(0.1)
                    self.emit_position(self.device_name_slow, "slow", curr_pos)
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
                print(ppos)

        status, curr_pos = self.req_update_pos(self.device_name_fast,self.fast_start)
        time.sleep(0.1)
        self.emit_position(self.device_name_fast, "fast", curr_pos)
        time.sleep(0.1)
        status, curr_pos = self.req_update_pos(self.device_name_slow,self.slow_start)
        time.sleep(0.1)
        self.emit_position(self.device_name_slow, "slow", curr_pos)
        time.sleep(0.1)
        self.update_progress.emit(0)
        time.sleep(0.1)
        self.toggle_ui.emit(True)
        time.sleep(0.1)
        self.toggle_button.emit(True)
        time.sleep(0.1)
        return

