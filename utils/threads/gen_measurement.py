from PyQt5.QtCore import QThread, pyqtSignal
import numpy as np
import time


class MichelsonMeasurement(QThread):
    change_linear_position = pyqtSignal(float)
    change_piezo_position = pyqtSignal(float)
    toggle_ui = pyqtSignal(bool)
    update_progress = pyqtSignal(int)
    update_label = pyqtSignal(tuple)
    toggle_button = pyqtSignal(bool)
    device_disconnect = pyqtSignal(str)

    def __init__(self, device_manager: DeviceManager,
                 l_start: float, l_end: float, l_step: int,
                 p_start: float, p_end: float, p_step: float):

        super(MichelsonMeasurement, self).__init__()

        self.device_manager = device_manager

        self.l_start = l_start
        self.l_end = l_end
        self.l_step = l_step
        self.p_start = p_start
        self.p_end = p_end
        self.p_step = p_step
        
        self.terminate = False


    def run(self) -> None:

        # turn off UI input elements while measuring
        self.toggle_ui.emit(False)
        self.toggle_button.emit(False)

        start = time.time()

        lpositions = np.linspace(self.l_start, self.l_end, self.l_step)
        ppositions = np.linspace(self.p_start, self.p_end, self.p_step)


        for l, lpos in enumerate(lpositions):
            try:
                self.device_manager.goto_linear(lpos)
                self.change_linear_position.emit(self.device_manager.get_linear_pos())
            except:
                self.device_disconnect.emit("Linear Stage")
                self.terminate = True

            for p, ppos in enumerate(ppositions):
                p += 1
                try:
                    self.device_manager.goto_piezo(ppos)
                    self.change_piezo_position.emit(self.device_manager.get_piezo_pos())
                except:
                    self.device_disconnect.emit("Piezo Stage")
                    self.terminate = True

                try:
                    final = (l+1) == len(lpositions) and p == len(ppositions)
                    self.device_manager.trig(low_handshake=(not final))
                except:
                    self.device_disconnect.emit("Linear Stage")
                    self.terminate = True

                elapsed = time.time() - start
                remaining = elapsed / max(l * len(ppositions) + p, 1) * len(lpositions) * len(ppositions) - elapsed

                self.update_label.emit((elapsed, remaining))

                self.update_progress.emit(int((l * len(ppositions) + p) / len(lpositions) / len(ppositions) * 100))
                
                if self.terminate:
                    self.update_progress.emit(0)
                    self.toggle_ui.emit(True)
                    self.toggle_button.emit(True)
                    return
                print(ppos)

        self.update_progress.emit(100)
        self.toggle_ui.emit(True)
        self.toggle_button.emit(True)

