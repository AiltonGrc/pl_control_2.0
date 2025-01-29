from PyQt5.QtCore import QThread, pyqtSignal
import numpy as np
import time
from la_com_request_parser.com_request_parser import ReqParser
from __future__ import annotations

class Movement_request(QThread):
    dev_return = pyqtSignal(list)
    movement_done = pyqtSignal(bool)


    def __init__(self,connections, dev:str, action: str , pos: float | None):

        super(Movement_request, self).__init__()

        self.com = connections
        self.action = action
        self.device = dev
        self.pos = pos
	
    def run(self) -> None:

        status, pos = False, float(0)

        if self.action == "set_pos":
            status, curr_pos = self.com.request(self.device, "set_pos",self.pos)
        elif self.action == "get_pos":
            status, curr_pos = self.com.request(self.device, "get_pos")

        time.sleep(0.1)
        self.dev_return.emit([status,curr_pos])

        return

