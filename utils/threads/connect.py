from PyQt5.QtCore import QThread, pyqtSignal
from la_com_request_parser.com_request_parser import ReqParser


class Connect(QThread):
    connected = pyqtSignal(bool)
    message = pyqtSignal(str)
    connect_message = pyqtSignal(str)

    def __init__(self, com: ReqParser, devices):
        super(Connect, self).__init__()

        self.com = com
        self.devices = devices

    def run(self):
        status,status_,_ = True,False,None
        for dev in self.devices:
            try:
                self.message.emit(f"Connecting {dev}")
                status_, _ = self.com.request(dev, "connect")
                if status_:
                    self.connect_message.emit("{} Connected.".format(dev))
                else:
                    self.connect_message.emit("{} Not Connected due to {}.".format(dev,_))
                self.connected.emit(status_)
            except:
                self.connect_message.emit("{} Not Connected due to {}.".format(dev, _))
                self.connected.emit(status_)

