from PyQt5.QtCore import QThread, pyqtSignal
import numpy as np
import time


class Test_thread(QThread):
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

    def __init__(self,connections,device_name: str):

        super(Test_thread, self).__init__()

        self.com = connections
        self.device_name = device_name
        
    def run(self) -> None:

        # turn off UI input elements while measurin

        start = time.time()

        self.toggle_ui.emit(False)
        self.toggle_button.emit(False)
        self.offset = 0
        self.l_start = int(40)
        self.l_end = int(60)
        self.l_step = int((self.l_end-self.l_start)/30)
        self.terminate = False

        lpositions = np.linspace(self.l_start, self.l_end, 20)
        
        
        i=0
        
        while i<10: 
        
            for l, lpos in enumerate(lpositions,1):
                
                try:
                    print(lpos)
                    status, curr_pos = self.com.request(self.device_name, "set_pos",lpos)
                except Exception as e:
                    print("Exception: ",e)
                    self.device_disconnect.emit(self.device_name)
                    self.terminate = True
    
                elapsed = time.time() - start
    
                remaining = elapsed / l * len(lpositions) - elapsed
    
                self.update_label.emit((elapsed, remaining))
                self.update_position.emit(str(float(curr_pos)-float(self.offset)))
    
                self.update_progress.emit(int(l)/ len(lpositions) * 100)
                
                print("Elapsed time: ", elapsed)
                print("Thread Name: ", QThread.currentThread())
                
                if self.terminate:
                    print("Terminating Thread")
                    self.update_progress.emit(0)
                    self.toggle_ui.emit(True)
                    self.toggle_button.emit(True)
                    return
                
            lpositions = np.flip(lpositions)
            print(lpositions)
            
            i+=1

        self.update_progress.emit(100)
        time.sleep(0.1)
        self.toggle_ui.emit(True)
        time.sleep(0.1)
        self.toggle_button.emit(True)
        time.sleep(0.1)

