from PyQt5.QtCore import QThread, pyqtSignal
import time
import numpy as np

class LoopLinear(QThread):
    change_linear_position = pyqtSignal(str)
    toggle_ui = pyqtSignal(bool)
    toggle_button = pyqtSignal(bool)
    device_disconnect = pyqtSignal(str)

    def __init__(self, connections, iterations: int, start_pos: float, end_pos: float):
        super(LoopLinear, self).__init__()

        self.com = connections
        self.iterations = iterations
        self.start_pos = start_pos
        self.end_pos = end_pos
        
        self.terminate = False

    def run(self) -> None:
        # turn off UI input elements while measuring
        self.toggle_ui.emit(False)
        self.toggle_button.emit(False)

        for _ in range(self.iterations):
            try:
                status,pos = self.com.request("michelson_linear", "set_pos", self.start_pos)
                self.change_linear_position.emit(str(pos))
            except:
                self.device_disconnect.emit("Linear Stage")
                self.terminate = True
            
            if self.terminate: break

            try:
                status,pos = self.com.request("michelson_linear", "set_pos", self.end_pos)
                self.change_linear_position.emit(str(pos))
            except:
                self.device_disconnect.emit("Linear Stage")
                self.terminate = True
            
            if self.terminate: break
        
        self.toggle_ui.emit(True)
        self.toggle_button.emit(True)



class LoopPiezo(QThread):
    change_piezo_position = pyqtSignal(str)
    toggle_ui = pyqtSignal(bool)
    toggle_button = pyqtSignal(bool)
    device_disconnect = pyqtSignal(str)

    def __init__(self, connections, iterations: int, start_pos: float, end_pos: float, steps: int):
        super(LoopPiezo, self).__init__()

        self.com = connections
        self.iterations = iterations
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.steps = steps
        
        self.terminate = False

    def run(self):

        # turn off UI input elements while measuring
        self.toggle_ui.emit(False)
        self.toggle_button.emit(False)

        for _ in range(self.iterations):
            for i in range(self.steps):
                try:
                    pos = self.start_pos + i * (self.end_pos - self.start_pos) / self.steps
                    status,pos = self.com.request("michelson_piezo", "set_pos", pos)
                    self.change_piezo_position.emit(str(pos))
                except:
                    self.device_disconnect.emit("Piezo Stage")
                    self.terminate = True
                if self.terminate: break
                time.sleep(0.3)

            try:
                status,pos = self.com.request("michelson_piezo", "set_pos", self.end_pos)
                self.change_piezo_position.emit(str(pos))
            except:
                self.device_disconnect.emit("Piezo Stage")
                self.terminate = True

            for i in reversed(range(self.steps)):
                time.sleep(0.3)
                try:
                    pos = self.start_pos + i * (self.end_pos - self.start_pos) / self.steps
                    status,pos = self.com.request("michelson_piezo", "set_pos", pos)
                    self.change_piezo_position.emit(str(pos))
                except:
                    self.device_disconnect.emit("Piezo Stage")
                    self.terminate = True
                if self.terminate: break
            if self.terminate: break

        self.toggle_ui.emit(True)
        self.toggle_button.emit(True)