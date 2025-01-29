from PyQt5.QtCore import QThread, pyqtSignal
import numpy as np
import time


class Measurement_IV_notrig(QThread):
    device_string = pyqtSignal(str)
    change_device_position = pyqtSignal(float)
    toggle_ui = pyqtSignal(bool)
    update_progress = pyqtSignal(int)
    update_label = pyqtSignal(tuple)
    toggle_button = pyqtSignal(bool)
    device_disconnect = pyqtSignal(str)
    update_position = pyqtSignal(tuple)
    return_iv_data = pyqtSignal(tuple)
    process_done = pyqtSignal(bool)

    def __init__(self,connections,
                 l_start: float, l_end: float, l_step: int, device_name: str):

        super(Measurement_IV_notrig, self).__init__()

        self.com = connections

        self.l_start = l_start
        self.l_end = l_end
        self.l_step = int(l_step)
        self.device_name = device_name
        self.terminate = False	
        self.curr_V = l_start
        self.curr_I = 0
        self.i_v_list = [[],[]] #Voltage and current data
        
	
    def run(self) -> None:

        # turn off UI input elements while measuring
        self.toggle_ui.emit(False)
        self.toggle_button.emit(False)

        start = time.time()
        
        lpositions = np.linspace(self.l_start, self.l_end, self.l_step)
        

        for l, lpos in enumerate(lpositions,1):
            try:
                
                status, (self.curr_V,self.curr_I) = self.com.request(self.device_name, "apply_voltage" ,lpos)                 
            except:
                self.device_disconnect.emit(self.device_name)
                self.terminate = True
                
                self.i_v_list[0].append(self.curr_V)
                self.i_v_list[1].append(self.curr_I)
                
            elapsed = time.time() - start

            remaining = elapsed / l * len(lpositions) - elapsed

            self.update_label.emit((elapsed, remaining))
            self.update_position.emit((self.curr_V,self.curr_I))
            self.return_iv_data.emit((self.curr_V,self.curr_I))

            self.update_progress.emit(int(l)/ len(lpositions) * 100)
            
            if self.terminate:
                status, (self.curr_V, self.curr_I) = self.com.request(self.device_name, "apply_voltage", self.l_start)
                time.sleep(0.1)
                self.update_position.emit((self.curr_V, self.curr_I))
                time.sleep(0.1)
                self.update_progress.emit(100)
                time.sleep(0.1)
                self.toggle_ui.emit(True)
                time.sleep(0.1)
                self.toggle_button.emit(True)
                time.sleep(0.1)
                return

            # print(lpos,self.curr_V,self.curr_I)

        status, (self.curr_V, self.curr_I) = self.com.request(self.device_name, "apply_voltage", self.l_start)
        time.sleep(0.1)
        self.update_position.emit((self.curr_V, self.curr_I))
        time.sleep(0.1)
        self.process_done.emit(True)
        time.sleep(0.1)
        self.update_progress.emit(100)
        time.sleep(0.1)
        self.toggle_ui.emit(True)
        time.sleep(0.1)
        self.toggle_button.emit(True)
        time.sleep(0.1)


