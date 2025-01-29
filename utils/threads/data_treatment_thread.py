from PyQt5.QtCore import QThread, pyqtSignal
import os

from .PL_data_treatment.data_fit import Fss_fit,Pow_fit,IV_fit

class Data_treatment(QThread):
    result_dict_signal = pyqtSignal(dict,str)
    message = pyqtSignal(str)

    def __init__(self, file: str, plot: bool = False, type='fss'):

        super(Data_treatment, self).__init__()

        self.file = file
        self.plot = plot
        self.wl = []
        self.intensity = []
        self.fit_result_exciton = {}
        self.fit_result_meas = {}
        self.result_dict = {}

    def run(self):

        print('thread running')

        if self.file == 'test':
            meas_axis = 'test'
            self.result_dict = Fss_fit(root=self.file, filename=self.file,flag_plot=self.plot)

        elif not os.path.isfile(self.file) or self is None:
            print("Error in filename")
            return

        elif 'polarization' in self.file.lower().strip():
            meas_axis = 'polarization'
            self.result_dict = Fss_fit(root=None,filename=self.file,flag_plot=self.plot)

        elif 'power_map' in self.file.lower().strip():
            meas_axis = 'power_map'
            self.result_dict = Pow_fit(root=None, filename=self.file, flag_plot=self.plot)

        elif 'voltage' in self.file.lower().strip():
            meas_axis = 'voltage'
            self.result_dict = IV_fit(root=None, filename=self.file, flag_plot=self.plot)

        else:
            print("No fitting procedure setup!")
            return

        print('thread done')
        self.result_dict_signal.emit(self.result_dict,meas_axis)

        return
