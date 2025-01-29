from PyQt5.QtCore import QThread, pyqtSignal
from PyQt5.QtWidgets import QDialog,QComboBox,QCheckBox,QDialogButtonBox,QFormLayout
from scipy import interpolate
import numpy as np
import time

class InputDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.piezo = QComboBox(self)
        self.spectrometer_connect = QCheckBox(self)

        self.piezo.addItem("PL Piezo")
        self.piezo.addItem("Imaging Piezo")
        self.spectrometer_connect.setChecked(True)

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok , self);

        layout = QFormLayout(self)
        layout.addRow("Piezo: ", self.piezo)
        layout.addRow("Connect to Spectrometer: ", self.spectrometer_connect)
        layout.addWidget(buttonBox)

        print(self.spectrometer_connect.isChecked())

        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)

    def getInputs(self):
        return (self.piezo.currentText(), self.spectrometer_connect.isChecked())

    def close(self):
        self.reject()

class show_dialogue_box(QThread):
    piezo_name_request = pyqtSignal(str)
    spectrometer_connect_state = pyqtSignal(bool)

    dialog = InputDialog()
    dialog.exec()

    def kill(self):
        piezo, spectrometer_request = dialog.getInputs()
        dialog.close()
        self.piezo_name_request.emit(piezo)
        self.spectrometer_connect_state.emit(spectrometer_request)



