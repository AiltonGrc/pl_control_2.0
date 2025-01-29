from PyQt5.QtCore import QTimer
from PyQt5.QtWidgets import QDialog,QComboBox,QCheckBox,QDialogButtonBox,QFormLayout
import time

class StartupDialog(QDialog):

    def __init__(self):
        super().__init__()

        self.piezo = QComboBox(self)
        self.spectrometer_connect = QCheckBox(self)
        self.cw_3900s = QCheckBox(self)
        self.power_meter = QCheckBox(self)
        self.kethley_connect = QCheckBox(self)
        self.debug_mode = QCheckBox(self)

        self.piezo.addItem("PL Piezo")
        self.piezo.addItem("Imaging Piezo")
        self.piezo.addItem("None")
        self.spectrometer_connect.setChecked(True)
        self.kethley_connect.setChecked(False)
        self.cw_3900s.setChecked(False)
        self.power_meter.setChecked(False)
        self.debug_mode.setChecked(False)

        self.tmrClose = QTimer(self)
        self.tmrClose.setInterval(30000)
        self.tmrClose.timeout.connect(self.tmrCloseTimeout)

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok , self)

        layout = QFormLayout(self)
        layout.addRow("Piezo: ", self.piezo)
        layout.addRow("Connect to Spectrometer: ", self.spectrometer_connect)
        layout.addRow("Connect to kethley: ", self.kethley_connect)
        layout.addRow("CW 3900s Motor: ", self.cw_3900s)
        layout.addRow("Power Meter: ", self.power_meter)
        layout.addRow("Debug Mode: ", self.debug_mode)
        layout.addWidget(buttonBox)

        print(self.spectrometer_connect.isChecked())

        self.tmrClose.start()

        buttonBox.accepted.connect(self.getInputs)
        buttonBox.rejected.connect(self.reject)

    def getInputs(self):
        self.tmrClose.stop()
        time.sleep(0.1)
        self.accept()

    def tmrCloseTimeout(self):
        print("Timed Out")
        self.tmrClose.stop()
        time.sleep(0.1)
        self.reject()