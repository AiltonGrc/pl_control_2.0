from PyQt5.QtCore import QThread, pyqtSignal
import numpy as np
import time, traceback


class XYMeasurement_fast(QThread):
    device_string = pyqtSignal(str)
    change_long_position = pyqtSignal(float)
    change_fast_position = pyqtSignal(float)
    toggle_ui = pyqtSignal(bool)
    update_progress = pyqtSignal(int)
    update_label = pyqtSignal(tuple)
    toggle_button = pyqtSignal(bool)
    device_disconnect = pyqtSignal(str)
    update_pos1 = pyqtSignal(str)
    update_pos2 = pyqtSignal(str)
    toggle_trigger = pyqtSignal(bool)
    reset_frames = pyqtSignal(bool)

    def __init__(self,connections,positions_array,device_name_xy: str,starting_pos:list,frames:list):

        super(XYMeasurement_fast, self).__init__()

        self.com = connections

        self.device_name = device_name_xy
        self.terminate = False
        self.starting_pos = starting_pos
        self.nodes = positions_array
        self.frames = frames
    
    def run(self) -> None:
        
        
        # turn off UI input elements while measuring
        self.toggle_ui.emit(False)
        self.toggle_button.emit(False)

        start = time.time()
        
        print(len(self.nodes),self.nodes)

        curr_posx = curr_posy = None
        trigger_now = False

        for node in np.arange(int(len(self.nodes))):
            try:
                
                (xpos,ypos) = self.nodes[node][0],self.nodes[node][1]
                
                # print("target: ", xpos,ypos, type(xpos) , type(ypos) )
                
                xpos = np.round(float(xpos),4)
                ypos = np.round(float(ypos),4)
                
                # print("target rounded: ", xpos,ypos, type(xpos) , type(ypos) )
                
                self.com.request(self.device_name, "set_pos", f"{xpos};{ypos}")

                try:
                    delta_x = abs(xpos - curr_posx)
                    delta_y = abs(ypos - curr_posy)

                    if abs(delta_x - delta_y) >= 1e-3:
                        trigger_now = True


                except Exception as e:
                    print(e)
                    trigger_now = True


                # print(status)
                
                # print("acquired: ", curr_posx,curr_posy, type(curr_posx) , type(curr_posy) )
                
                curr_posx,curr_posy = str(np.round(xpos,4)),str(np.round(ypos,4))
                
                # print("acquired rounded: ", curr_posx,curr_posy, type(curr_posx) , type(curr_posy) )
                
                self.update_pos1.emit(curr_posx)
                time.sleep(0.1)
                self.update_pos2.emit(curr_posy)
            except Exception as err:
                traceback.print_tb(err.__traceback__)
                self.device_disconnect.emit(self.device_name)
                self.terminate = True


            if trigger_now:
                try:
                    for _ in np.arange(int(self.frames[1])):
                        self.com.request("michelson_linear", "trigger_spectrometer", False)
                        time.sleep(float(self.frames[2])+50e-3)

                except:
                    self.terminate = True

            elapsed = time.time() - start
            try:
                remaining = remaining = elapsed / (node+1) * (len(self.nodes)+1) - elapsed
            except:
                remaining = 100
            
            self.update_label.emit((elapsed, remaining))

            self.update_progress.emit(int(node)/ len(self.nodes) * 100)
            
            if self.terminate:
                try:
                    status, (curr_posx, curr_posy) = self.com.request(self.device_name_xy, "set_pos", f"{self.starting_pos[0]};{self.starting_pos[1]}")
                    curr_posx, curr_posy = np.round(curr_posx, 4), np.round(curr_posy, 4)
                    time.sleep(0.1)
                    self.update_pos1.emit(str(curr_posx))
                    time.sleep(0.1)
                    self.update_pos2.emit(str(curr_posy))
                    time.sleep(0.1)
                except:
                    pass

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
        
        try:
            time.sleep(0.1)
            self.com.request("michelson_linear", "trigger_spectrometer", True) 
        except:
            pass
            
        try:
            status, (curr_posx,curr_posy) = self.com.request(self.device_name_xy, "set_pos", f"{self.starting_pos[0]};{self.starting_pos[1]}")
            curr_posx,curr_posy = np.round(curr_posx,4),np.round(curr_posy,4)
            time.sleep(0.1)
            
            self.update_pos1.emit(str(curr_posx))
            time.sleep(0.1)
            self.update_pos2.emit(str(curr_posy))
            time.sleep(0.1)
        except:
            pass

        time.sleep(0.1)
        self.update_progress.emit(100)
        time.sleep(0.1)
        self.toggle_ui.emit(True)
        time.sleep(0.1)
        self.toggle_button.emit(True)
        time.sleep(0.1)

        return

