# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 16:43:09 2022

@author: Ailton Garcia/AK124573
"""

# TODO include PS4 controller to move at a fixed rate which the motor/connection can keep up
# TODO make movement on a separate thread


from __future__ import annotations

import numpy as np
import regex as re
from datetime import date, datetime
from pathlib import Path
import csv,os,sys,json,time
import matplotlib.pyplot as plt

from PyQt5 import QtWidgets, QtCore
from PyQt5.QtWidgets import QFileDialog, QDialog, QMessageBox, \
    QTableWidgetItem, QRadioButton

from PyQt5.QtGui import QIcon

from utils.ui.PL_GUI import Ui_MainWindow
from utils.ui.q_double_range_validator import QDoubleRangeValidator
from utils.ui.DialogWindow import StartupDialog

from utils.threads.michelson_measurement import MichelsonMeasurement
from utils.threads.meas_1d import Measurement_1D
from utils.threads.meas_pol_sup import Measurement_pol_sup
from utils.threads.meas_IV import Measurement_IV
from utils.threads.meas_IV_notrig import Measurement_IV_notrig
from utils.threads.meas_2d import Measurement_2D
from utils.threads.xy_meas import XYMeasurement
from utils.threads.fast_xy_meas import XYMeasurement_fast
from utils.threads.mesh_measurement import Meshed_xy_measurement
from utils.threads.fast_xy_meas_rad import rad_limited_map
from utils.threads.loop import LoopLinear, LoopPiezo
from utils.threads.test_thread import Test_thread

from utils.threads.movement_thread import Movement_request

from utils.threads.data_treatment_thread import Data_treatment

from utils.device.ps4_controller import PS4_controller

from utils.ErrorLogger import ErrorLogger

from utils.threads.connect import Connect

from la_com_request_parser.com_request_parser import ReqParser


class Pl_Control:  # is parenthesis needed here?

    def __init__(self):
        app = QtWidgets.QApplication(sys.argv)
        self.threadpool = QtCore.QThreadPool()
        self.main_window = QtWidgets.QMainWindow()
        # layout = QVBoxLayout() #commented out on 29.03.23
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.main_window)

        # disconnect devices upon closing
        app.aboutToQuit.connect(self.disconnect)
        # check boxes
        self.check_text()

        self.logo = "PL_icon.png"
        self.main_window.setWindowIcon(QIcon(os.path.join("utils", "ui", self.logo)))
        self.main_window.raise_()

        # PlaceHolder variables
        self.connect_thread = None
        self.measurement_thread = None
        self.loop_thread = None
        self.com = None
        self.white_led_status = False
        self.blue_led_status = False
        self.red_led_status = False
        self.flip_mirror_state = False
        self.thr_bus_state = False
        self.periscope_piezo_state = False
        self.michelson_piezo_state = False
        self.camera_piezo_state = False
        self.spectrometer_connected = False
        self.dialog = StartupDialog()
        self.path = str(self.ui.base_directory.text()).strip()
        self.thread_frames = None

        # Devices names from device server
        self.piezo_name = "focus_piezo"  # "imaging_piezo"
        self.lin_stage_name = "michelson_linear"
        self.piezo_stage_name = "michelson_piezo"
        self.xy_stage_name = "cryo_xy"
        self.spectrometer_name = "magnetolab_spectrometer"
        self.laser_power_name = "laser_power_hwp"
        # self.thr_lab_piezo_name = ['magneto_michelson_in', 'magneto_michelson_out']
        self.fss_name = "fss_hwp"
        self.smu_name = "kethley_magneto"  # last on the list
        self.thr_bus_controller = "tl_elliptec"
        self.camera_piezo_name = "pl_bs"
        self.periscope_piezo_name = "periscope_piezo"
        self.michelson_piezo_name = "magneto_michelson_elliptec"
        self.handheld = 'ps4_controller'
        self.cw_name = '3900s_motor'
        self.cw_hwp_name = "3900s_hwp"
        self.arduino_name = "arduino_led_controller"
        self.pwr_meter_name = "magnetolab_powermeter"
        self.curr_path = "spec"
        #TODO: FORCE TO SET SPEC PATH AT START

        # XY stage has to be 3rd on the list, z-piezo has to be first and spectrometer has to be last
        self.device_name_list = [self.piezo_stage_name, self.lin_stage_name, self.xy_stage_name,
                                 self.laser_power_name, self.fss_name,  self.thr_bus_controller]#self.arduino_name,

        # Dialog to select piezo, if spectrometer and if kethley
        self.startup_dialog()

        # Store positions
        self.curr_x, self.curr_y = 25, 25
        self.start_x, self.start_y = 25, 25
        self.curr_z = 0
        self.curr_fss = 0
        self.curr_laser = 0
        self.curr_V = 0
        self.curr_I = 0
        self.curr_cw = [0, 0]
        self.curr_cw_hwp = 0
        self.cw_limits = [0,25,600,1050]
        self.curr_linear_pos = 0
        self.curr_mich_piezo_pos = 25
        self.temp_holder = ''
        self.running_measurement = None
        self.fss_offset = get_float(self.ui.pol_offset.text())
        self.laser_offset = get_float(self.ui.laser_offset.text())
        self.cw_laser_offset = get_float(self.ui.cw_laser_offset.text())
        self.michelson_equal_length = get_float(self.ui.equal_path_length.text())
        self.wl_calibration_file = os.path.join("utils", "cw_calibration_data.txt")
        self.wl_calibration = [[],[]]

        self.import_cw_calib(self.wl_calibration_file)

        # Starting graphs
        self.iv_data, = self.ui.iv_curve.canvas.ax.plot(
            [], [], '-b', label='I-V Curve')
        self.i_v_list = [[], []]
        self.iv_data_neg, = self.ui.iv_curve.canvas.ax.plot(
            [], [], '-r', label='I-V Curve')
        self.i_v_list_neg = [[], []]

        self.intensity_data, = self.ui.treated_data.canvas.ax.plot(
            [], [], '-b', label='Intenity')
        self.intensity_fit_data, = self.ui.treated_data.canvas.ax.plot(
            [], [], '--r', label='Peak Fit')

        self.experiment_data, = self.ui.treated_data_2.canvas.ax.plot(
            [], [], '-b', label='Experiment Data')
        self.experiment_data_fit, = self.ui.treated_data_2.canvas.ax.plot(
            [], [], '--r', label='Experiment Data Fit')

        # virtual devices and processes
        self.logger = ErrorLogger("Error_logfile")

        # creating threads to avoid simultaneous requests
        self.main_thread = QtCore.QThread.currentThread()
        self.device_thread = QtCore.QThread()
        self.measurement_thread = QtCore.QThread()
        self.thread_queue = list()
        self.ps4_controller_thread = QtCore.QThread()
        self.treatment_thread = QtCore.QThread()
        self.movement_thread = QtCore.QThread()

        # self.connect_ps_controller()

        # connect devices
        try:
            self.com = self.connect()
        except Exception as e:
            self.show_message("Error in connection!", e)

        print('connection', self.com)

        # if not noconnect:
        #     self.connect()
        if self.com is not None:
            try:
                status, _ = self.com.request("michelson_linear", "set_whitelight", False)
                status, _ = self.com.request("michelson_linear", "set_red_led", False)
                status, _ = self.com.request("michelson_linear", "set_blue_led", False)
            except Exception as e:
                self.show_message("Error in setting LED status due to {}".fomat(e),e)

        # minor startup things
        self.update_fss_step()
        self.laser_step()
        self.volt_frames()
        self.space_map_frames()
        self.set_label_filename()
        self.standard_settings()
        self.update_offsets()

        self.piezos_positions("pl",forced=True)

        self.connect_ps_controller()

        # validators for inputs
        self.ui.linear_start.setValidator(
            QDoubleRangeValidator(bottom=0., top=305., decimals=1))
        self.ui.linear_end.setValidator(
            QDoubleRangeValidator(bottom=0., top=305., decimals=1))
        self.ui.linear_goto.setValidator(
            QDoubleRangeValidator(bottom=0., top=305., decimals=1))
        self.ui.piezo_start.setValidator(
            QDoubleRangeValidator(bottom=0., top=50., decimals=4))
        self.ui.piezo_end.setValidator(
            QDoubleRangeValidator(bottom=0., top=50., decimals=4))
        self.ui.piezo_goto.setValidator(
            QDoubleRangeValidator(bottom=0., top=50., decimals=4))

        # connect button of controller to function
        try:
            self.ps4_controller_thread.button.connect(self.ps4_action)
        except:
            self.ps4_controller_thread = None

        # buttons michelson
        self.ui.linear_loop.clicked.connect(self.loop_linear)
        self.ui.piezo_loop.clicked.connect(self.loop_piezo)
        self.ui.save_config.clicked.connect(self.save_config)
        self.ui.load_config.clicked.connect(self.load_config)
        if self.spectrometer_name in self.device_name_list:
            self.ui.retry_connection.clicked.connect(lambda: (self.connect(self.device_name_list[1:])))
        else:

            self.ui.retry_connection.clicked.connect(lambda: (self.connect(self.device_name_list[:])))
        self.ui.go_button.clicked.connect(self.goto_linear)
        self.ui.go_button_2.clicked.connect(self.goto_piezo)
        self.ui.update_equal_path_length.clicked.connect(self.update_equal_length_path_value)
        self.ui.measurement_button.clicked.connect(self.michelson_measurement)
        self.ui.apply_button.clicked.connect(self.standard_settings)

        # =============================================================================
        #         # buttons movement
        # =============================================================================
        self.ui.goto_xy.clicked.connect(lambda: (self.make_move(0.1, 'xy', 1, goto=[get_float(self.ui.target_x.text()),
                                                                                    get_float(
                                                                                        self.ui.target_y.text())])))
        self.ui.goto_z.clicked.connect(lambda: (self.make_move(0.1, 'z', 1, goto=get_float(self.ui.target_z.text()))))

        # X axis
        self.ui.xp_x01.clicked.connect(lambda: (self.make_move(0.1, 'x', -1)))
        self.ui.xp_x1.clicked.connect(lambda: (self.make_move(1, 'x', -1)))
        self.ui.xp_x10.clicked.connect(lambda: (self.make_move(10, 'x', -1)))
        self.ui.xm_x01.clicked.connect(lambda: (self.make_move(0.1, 'x', 1)))
        self.ui.xm_x1.clicked.connect(lambda: (self.make_move(1, 'x', 1)))
        self.ui.xm_x10.clicked.connect(lambda: (self.make_move(10, 'x', 1)))
        # Y axis
        self.ui.yp_x01.clicked.connect(lambda: (self.make_move(0.1, 'y', 1)))
        self.ui.yp_x1.clicked.connect(lambda: (self.make_move(1, 'y', 1)))
        self.ui.yp_x10.clicked.connect(lambda: (self.make_move(10, 'y', 1)))
        self.ui.ym_x01.clicked.connect(lambda: (self.make_move(0.1, 'y', -1)))
        self.ui.ym_x1.clicked.connect(lambda: (self.make_move(1, 'y', -1)))
        self.ui.ym_x10.clicked.connect(lambda: (self.make_move(10, 'y', -1)))
        # Z axis
        self.ui.zp_x01.clicked.connect(lambda: (self.make_move(0.1, 'z', 1)))
        self.ui.zp_x1.clicked.connect(lambda: (self.make_move(1, 'z', 1)))
        self.ui.zp_x10.clicked.connect(lambda: (self.make_move(10, 'z', 1)))
        self.ui.zm_x01.clicked.connect(lambda: (self.make_move(0.1, 'z', -1)))
        self.ui.zm_x1.clicked.connect(lambda: (self.make_move(1, 'z', -1)))
        self.ui.zm_x10.clicked.connect(lambda: (self.make_move(10, 'z', -1)))
        self.ui.home_z_axis.clicked.connect(
            lambda: (self.make_move(0.1, 'z', -1, goto=50)))
        # laser_hwp
        self.ui.hwp_laser_p05.clicked.connect(lambda: (self.make_move(0.5, 'laser', 1)))
        self.ui.hwp_laser_p1.clicked.connect(lambda: (self.make_move(1, 'laser', 1)))
        self.ui.hwp_laser_p2.clicked.connect(lambda: (self.make_move(2, 'laser', 1)))
        self.ui.hwp_laser_p5.clicked.connect(lambda: (self.make_move(5, 'laser', 1)))
        self.ui.hwp_laser_m05.clicked.connect(lambda: (self.make_move(0.5, 'laser', -1)))
        self.ui.hwp_laser_m1.clicked.connect(lambda: (self.make_move(1, 'laser', -1)))
        self.ui.hwp_laser_m2.clicked.connect(lambda: (self.make_move(2, 'laser', -1)))
        self.ui.hwp_laser_m5.clicked.connect(lambda: (self.make_move(5, 'laser', -1)))
        # fss_hwp
        self.ui.fss_p05.clicked.connect(lambda: (self.make_move(0.5, 'fss', 1)))
        self.ui.fss_p1.clicked.connect(lambda: (self.make_move(1, 'fss', 1)))
        self.ui.fss_p2.clicked.connect(lambda: (self.make_move(2, 'fss', 1)))
        self.ui.fss_p5.clicked.connect(lambda: (self.make_move(5, 'fss', 1)))
        self.ui.fss_m05.clicked.connect(lambda: (self.make_move(0.5, 'fss', -1)))
        self.ui.fss_m1.clicked.connect(lambda: (self.make_move(1, 'fss', -1)))
        self.ui.fss_m2.clicked.connect(lambda: (self.make_move(2, 'fss', -1)))
        self.ui.fss_m5.clicked.connect(lambda: (self.make_move(5, 'fss', -1)))

        # move to angles
        # fss
        self.ui.fss_goto_1.clicked.connect(lambda: (self.make_move(1, 'fss', -1,
                                                                   goto=get_float(self.ui.fss_angle_1.text()))))

        self.ui.fss_goto_2.clicked.connect(lambda: (self.make_move(1, 'fss', -1,
                                                                   goto=get_float(self.ui.fss_angle_2.text()))))

        self.ui.fss_goto_3.clicked.connect(lambda: (self.make_move(1, 'fss', -1,
                                                                   goto=get_float(self.ui.fss_angle_3.text()))))

        self.ui.fss_goto_4.clicked.connect(lambda: (self.make_move(1, 'fss', -1,
                                                                   goto=(get_float(self.ui.fss_angle_4.text()) / 2))))

        # laser
        self.ui.hwp_laser_goto_1.clicked.connect(lambda: (self.make_move(1, 'laser', -1,
                                                                         goto=get_float(
                                                                             self.ui.pow_map_angle_1.text()))))

        self.ui.hwp_laser_goto_2.clicked.connect(lambda: (self.make_move(1, 'laser', -1,
                                                                         goto=get_float(
                                                                             self.ui.pow_map_angle_2.text()))))

        self.ui.hwp_laser_goto_3.clicked.connect(lambda: (self.make_move(1, 'laser', -1,
                                                                         goto=get_float(
                                                                             self.ui.pow_map_angle_3.text()))))
        self.ui.hwp_laser_goto_4.clicked.connect(lambda: (self.make_move(1, 'laser', -1,
                                                                         goto=get_float(
                                                                             self.ui.pow_map_angle_4.text()))))


        #CW Laser


        self.ui.cw_hwp_laser_n.clicked.connect(lambda: (self.make_move(get_float(
                                                                             self.ui.cw_hwp_rel_n.text()), 'cw_hwp', -1)))
        self.ui.cw_hwp_laser_p.clicked.connect(lambda: (self.make_move(get_float(
                                                                             self.ui.cw_hwp_rel_p.text()), 'cw_hwp', 1)))
        self.ui.cw_hwp_goto1.clicked.connect(lambda: (self.make_move(1, 'cw_hwp', -1,
                                                                         goto=get_float(
                                                                             self.ui.cw_hwp_1.text()))))
        self.ui.cw_hwp_goto2.clicked.connect(lambda: (self.make_move(1, 'cw_hwp', -1,
                                                                     goto=get_float(
                                                                         self.ui.cw_hwp_2.text()))))

        #'cw_laser', 'cw_laser_wl'


        #absolute move in wl
        self.ui.cw_wl_goto.clicked.connect(lambda: (self.make_move(1, 'cw_laser_wl', -1,
                                                                      goto=get_float(
                                                                          self.ui.cw_wl_setp.text()))))
        self.ui.cw_wl_goto_2.clicked.connect(lambda: (self.make_move(1, 'cw_laser_wl', -1,
                                                                   goto=get_float(
                                                                       self.ui.cw_wl_setp_2.text()))))
        #relative move in wl
        self.ui.cw_laser_wl_m.clicked.connect(lambda: (self.make_move(get_float(self.ui.cw_wl_m_label.text()
                                                                                ), 'cw_laser_wl', -1)))
        self.ui.cw_laser_wl_p.clicked.connect(lambda: (self.make_move(get_float(self.ui.cw_wl_p_label.text()
                                                                                ), 'cw_laser_wl', 1)))
        #absolute move in motor pos
        self.ui.cw_motor_goto.clicked.connect(lambda: (self.make_move(1, 'cw_laser', -1,
                                                                   goto=get_float(
                                                                       self.ui.cw_motor_setp.text()))))
        self.ui.cw_motor_goto_2.clicked.connect(lambda: (self.make_move(1, 'cw_laser', -1,
                                                                     goto=get_float(
                                                                         self.ui.cw_motor_setp_2.text()))))

        #relative move in motor pos
        self.ui.cw_motor_m.clicked.connect(lambda: (self.make_move(get_float(self.ui.cw_motor_m_label.text()
                                                                                ), 'cw_laser', -1, )))
        self.ui.cw_motor_p.clicked.connect(lambda: (self.make_move(get_float(self.ui.cw_motor_p_label.text()
                                                                                ), 'cw_laser', 1, )))


        # store positions
        self.ui.add_position.clicked.connect(lambda: (self.store_position('stored_positions')))
        self.ui.remove_position.clicked.connect(lambda: (self.remove_position('stored_positions')))
        self.ui.go_to_position.clicked.connect(lambda: (self.goto_store_position('stored_positions')))

        self.ui.add_position_2.clicked.connect(lambda: (self.store_position('plane_positions')))
        self.ui.remove_position_2.clicked.connect(lambda: (self.remove_position('plane_positions')))
        self.ui.go_to_position_2.clicked.connect(lambda: (self.goto_store_position('plane_positions')))

        # Load/Save pos to file
        self.ui.load_pos.clicked.connect(lambda: (self.import_export_pos(True)))
        self.ui.save_pos.clicked.connect(lambda: (self.import_export_pos(False)))

        # =============================================================================
        # measurement buttons
        # =============================================================================

        # space_map
        self.ui.space_map_button.clicked.connect(self.space_map)
        self.ui.line_scan_button.clicked.connect(self.line_scan)
        self.ui.check_spacemap_points.clicked.connect(self.check_spacemap_points)

        # buttons fss
        self.ui.fss_button.clicked.connect(self.fss_measurement)

        # buttons power_map
        self.ui.power_map_button.clicked.connect(lambda: (self.power_map(laser="Green")))
        self.ui.power_map_button.clicked.connect(lambda: (self.power_map(laser="3900s")))

        #PLE button
        self.ui.ple_button.clicked.connect(self.ple_meas)
        self.ui.cw_pwr_read.clicked.connect(self.meas_cw_power)

        # buttons 2d_measurement
        self.ui.measure_2d_button.clicked.connect(self.setup_2d_meas)

        # button volt sweep
        self.ui.volt_sweep_button.clicked.connect(self.voltage_sweep)
        self.ui.set_v.clicked.connect(self.set_voltage)
        self.ui.meas_IV.clicked.connect(self.meas_iv_curve)
        self.ui.double_v_sweep.clicked.connect(self.meas_double_IV_curve)

        # Abort measurement
        self.ui.abort_measurent.clicked.connect(self.abort_current_measurement)

        # =============================================================================
        #         # buttons camera/led
        # =============================================================================
        self.ui.w_LED_button.clicked.connect(self.white_led_on_off)
        self.ui.b_LED_button.clicked.connect(self.blue_led_on_off)
        self.ui.r_LED_button.clicked.connect(self.red_led_on_off)
        self.ui.flip_mirror_button.clicked.connect(self.flip_mirror)  # replace with self.flip_mirror, test function for now
        self.ui.michelson_piezos_button.clicked.connect(lambda: (self.piezos_positions("michelson")))
        self.ui.pl_piezos_button.clicked.connect(lambda: (self.piezos_positions("pl")))
        self.ui.camera_piezo_button.clicked.connect(lambda: (self.piezos_positions("imaging")))
        self.ui.update_settings.clicked.connect(self.update_standard_settings)

        self.ui.periscope_jog_left.clicked.connect(lambda: (self.jog_piezo("periscope", "backward")))
        self.ui.periscope_jog_right.clicked.connect(lambda: (self.jog_piezo("periscope", "forward")))
        self.ui.michelson_jog_left.clicked.connect(lambda: (self.jog_piezo("michelson", "backward")))
        self.ui.michelson_jog_right.clicked.connect(lambda: (self.jog_piezo("michelson", "forward")))
        self.ui.camera_jog_left.clicked.connect(lambda: (self.jog_piezo("camera", "backward")))
        self.ui.camera_jog_right.clicked.connect(lambda: (self.jog_piezo("camera", "forward")))

        self.ui.update_periscope_in.clicked.connect(lambda: (self.update_piezo_pos("periscope", "in")))
        self.ui.update_periscope_out.clicked.connect(lambda: (self.update_piezo_pos("periscope", "out")))
        self.ui.update_michelson_in.clicked.connect(lambda: (self.update_piezo_pos("periscope", "in")))
        self.ui.update_michelson_out.clicked.connect(lambda: (self.update_piezo_pos("michelson", "out")))
        self.ui.update_camera_in.clicked.connect(lambda: (self.update_piezo_pos("camera", "in")))
        self.ui.update_camera_out.clicked.connect(lambda: (self.update_piezo_pos("camera", "out")))

        self.ui.update_offsets_button.clicked.connect(self.update_offsets)
        self.ui.updatefilename_button.clicked.connect(self.update_filename_spectrometer)
        self.ui.home_cw.clicked.connect(self.home_cw)
        self.ui.update_cw_calib.clicked.connect(lambda: (self.import_cw_calib(self.wl_calibration_file)))

        # =============================================================================
        #         update labels
        # =============================================================================

        # update values of self.ui.fss_step.text()
        self.ui.fss_frames.textChanged.connect(self.update_fss_step)
        self.ui.fss_i.textChanged.connect(self.update_fss_step)
        self.ui.fss_end.textChanged.connect(self.update_fss_step)

        #update values of ple step
        self.ui.ple_frames.textChanged.connect(self.ple_frames)
        self.ui.ple_start.textChanged.connect(self.ple_frames)
        self.ui.ple_end.textChanged.connect(self.ple_frames)

        # update values of self.ui.laser_step.text()
        self.ui.pow_map_frames.textChanged.connect(self.laser_step)
        self.ui.power_map_i.textChanged.connect(self.laser_step)
        self.ui.power_map_f.textChanged.connect(self.laser_step)
        self.ui.cw_pow_map_frames.textChanged.connect(self.laser_step)
        self.ui.cw_power_map_i.textChanged.connect(self.laser_step)
        self.ui.cw_power_map_f.textChanged.connect(self.laser_step)

        # update values of self.ui.volt_frames.text()
        self.ui.volt_max.textChanged.connect(self.volt_frames)
        self.ui.volt_min.textChanged.connect(self.volt_frames)
        self.ui.volt_step.textChanged.connect(self.volt_frames)

        # update values of self.ui.space_map_frames.text()
        self.ui.space_map_colums.textChanged.connect(self.space_map_frames)
        self.ui.space_map_lines.textChanged.connect(self.space_map_frames)

        # update values of 2dmeasurement
        self.ui.slow_axis_start.editingFinished.connect(self.meas_2d_frames)
        self.ui.slow_axis_steps.editingFinished.connect(self.meas_2d_frames)
        self.ui.slow_axis_end.editingFinished.connect(self.meas_2d_frames)
        self.ui.fast_axis_start.editingFinished.connect(self.meas_2d_frames)
        self.ui.fast_axis_steps.editingFinished.connect(self.meas_2d_frames)
        self.ui.fast_axis_end.editingFinished.connect(self.meas_2d_frames)

        # text edits
        self.ui.linear_start.editingFinished.connect(self.check_text)
        self.ui.linear_start.returnPressed.connect(self.check_text)
        self.ui.linear_start.editingFinished.connect(self.set_label_filename)
        self.ui.linear_start.returnPressed.connect(self.set_label_filename)
        self.ui.linear_end.editingFinished.connect(self.check_text)
        self.ui.linear_end.returnPressed.connect(self.check_text)
        self.ui.linear_end.editingFinished.connect(self.set_label_filename)
        self.ui.linear_end.returnPressed.connect(self.set_label_filename)
        self.ui.piezo_start.editingFinished.connect(self.check_text)
        self.ui.piezo_start.returnPressed.connect(self.check_text)
        self.ui.piezo_start.editingFinished.connect(self.set_label_filename)
        self.ui.piezo_start.returnPressed.connect(self.set_label_filename)
        self.ui.piezo_end.editingFinished.connect(self.check_text)
        self.ui.piezo_end.returnPressed.connect(self.check_text)
        self.ui.piezo_end.editingFinished.connect(self.set_label_filename)
        self.ui.piezo_end.returnPressed.connect(self.set_label_filename)

        # change filename
        self.ui.sample_name.editingFinished.connect(self.set_label_filename)
        self.ui.sample_name_2.editingFinished.connect(self.set_label_filename)
        self.ui.temperature.editingFinished.connect(self.set_label_filename)
        self.ui.filter_name.editingFinished.connect(self.set_label_filename)
        self.ui.slit.editingFinished.connect(self.set_label_filename)
        self.ui.exposure.editingFinished.connect(self.set_label_filename)
        self.ui.grating.currentIndexChanged.connect(self.set_label_filename)
        self.ui.center_wl.editingFinished.connect(self.set_label_filename)
        self.ui.start_wl.editingFinished.connect(self.set_label_filename)
        self.ui.end_wl.editingFinished.connect(self.set_label_filename)
        self.ui.cw_wl.editingFinished.connect(self.set_label_filename)
        self.ui.sources_group.buttonClicked.connect(self.set_label_filename)
        self.ui.measurements_group.buttonClicked.connect(self.set_label_filename)
        self.ui.space_map_colums.editingFinished.connect(self.set_label_filename)
        self.ui.space_map_lines.editingFinished.connect(self.set_label_filename)
        self.ui.column_step.editingFinished.connect(self.set_label_filename)
        self.ui.line_step.editingFinished.connect(self.set_label_filename)

        # steps
        self.ui.linear_step.valueChanged.connect(self.check_text)
        self.ui.linear_step.valueChanged.connect(self.set_label_filename)
        self.ui.piezo_step.valueChanged.connect(self.check_text)
        self.ui.piezo_step.valueChanged.connect(self.set_label_filename)

        # linear vel/acc/dec
        self.ui.vel_spin.valueChanged.connect(self.set_linear_velocity)
        self.ui.acc_spin.valueChanged.connect(self.set_linear_acceleration)
        self.ui.dec_spin.valueChanged.connect(self.set_linear_deceleration)

        # copy to clipboard
        self.ui.copy_filename.setIcon(QIcon(os.path.join("ui", "copy.png")))
        self.ui.copy_filename.clicked.connect(self.copy_filename)
        self.ui.copy_number_of_spectra.setIcon(
            QIcon(os.path.join("ui", "copy.png")))
        self.ui.copy_number_of_spectra.clicked.connect(
            self.copy_number_of_spectra)

        self.clipboard = app.clipboard()

        # show ui
        self.main_window.show()
        sys.exit(app.exec_())

    def __del__(self):
        """
        destructor makes sure that devices are disconnected and made available again
        set LED states
        :return: No return
        """

        try:
            self.com.request("michelson_linear", "set_red_led", False)
            self.com.request("michelson_linear", "set_blue_led", False)
            self.com.request("michelson_linear", "set_white_led", False)
            for dev in self.device_name_list:
                self.com.request(dev, "disconnect")
        except AttributeError:
            # devices were not initialized in the first place, cannot be disconnected, no need to interfere
            pass
        try:
            self.update_standard_settings()
        except Exception as e:
            print("Error in updating Settings.", e)

    #     Main control functions

    def ps4_action(self, button: str, value: float):

        """
        Gets button and value from PS controller

        :param button: button name (left,right,up,down L1,R1,LT,RT)
        :params value: float of multiplier for movement

        """

        multiplier, axis, direction = None, None, None

        if button == 'left':
            axis = 'x'
            multiplier = value
            direction = 1
        elif button == 'right':
            axis = 'x'
            multiplier = value
            direction = -1
        elif button == 'up':
            axis = 'y'
            multiplier = value
            direction = 1
        elif button == 'down':
            axis = 'y'
            multiplier = value
            direction = -1
        elif button == 'L1':
            axis = 'fss'
            multiplier = 1
            direction = -1
        elif button == 'R1':
            axis = 'fss'
            multiplier = 1
            direction = 1
        elif button == 'LT':
            axis = 'laser'
            multiplier = 1
            direction = -1
        elif button == 'RT':
            axis = 'laser'
            multiplier = 1
            direction = 1
        elif button == 'z-up':
            axis = 'z'
            multiplier = 2
            direction = 1
        elif button == 'z-down':
            axis = 'z'
            multiplier = 2
            direction = -1

        if None in (multiplier, axis, direction):
            pass
        else:
            try:
                self.make_move(multiplier,axis,direction)
            except Exception as e:
                self.show_message("Error in PS4 action: ", e)

    def connect_ps_controller(self):
        try:
            self.ps4_controller_thread = PS4_controller()
            self.ps4_controller_thread.start()
        except:
            self.ps4_controller_thread = None

        return

    def test_function(self):

        print("test Function")
        print(self.com.request(self.pwr_meter_name,"read"))

        return

    def treat_last_data(self, file: str):

        self.treatment_thread = Data_treatment(file=file, plot=False, type=self.running_measurement)
        self.treatment_thread.start()
        self.treatment_thread.result_dict_signal.connect(self.update_treated_data)

    def home_cw(self):
        status, _ = self.com.request(self.cw_name, "home")

    def wl_regression(self, target: float, file=None):

        try:
            if file is not None or self.wl_calibration is None:
                self.wl_calibration = self.import_cw_calib(file)


            if target <=self.cw_limits[0] or target>=self.cw_limits[1]:
                if target<=self.cw_limits[2] or target>=self.cw_limits[3]:
                    raise Exception("Motor position higher than set limit (25mm).")
                else:
                    pos = self.wl_calibration[0][0] * target ** 3 + self.wl_calibration[0][1] * target ** 2 + \
                          self.wl_calibration[0][2] * target ** 1 + self.wl_calibration[0][3]
            elif target<=self.cw_limits[2] or target>=self.cw_limits[3]:
                if target <= self.cw_limits[0] or target >= self.cw_limits[1]:
                    raise Exception("Motor position higher than set limit (25mm).")
                else:
                    pos = self.wl_calibration[1][0] * target ** 3 + self.wl_calibration[1][1] * target ** 2 + \
                          self.wl_calibration[1][2] * target ** 1 + self.wl_calibration[1][3]

            if pos >=self.cw_limits[1] or pos <= self.cw_limits[0]:
                if pos <=self.cw_limits[2] or pos>=self.cw_limits[3]:
                    pos = "Error"
                    raise Exception("Motor position higher than set limit (0.1 - 24.9mm).1")
            elif pos <=self.cw_limits[2] or pos>=self.cw_limits[3]:
                if pos >=self.cw_limits[1] or pos <= self.cw_limits[0]:
                    pos = "Error"
                    raise Exception("Motor position higher than set limit (0.1 - 24.9mm).2")

        except Exception as e:
            pos = "Error"
            print(e)
            self.show_message("Could not find motor position from motor fit or pos>25mm")

        return pos


    def import_cw_calib(self, file=None):


        try:
            if file is not None:
                file = file
            else:
                file = self.wl_calibration_file
        except Exception as e:
            self.show_message("Could not open calibration file.",e)

        wl = []
        motor_pos = []

        try:
            with open(file,'r') as f:
                next(f)
                for line in f:
                    line.strip()
                    columns=line.split(";")
                    wl.append(float(columns[0]))
                    motor_pos.append(float(columns[1]))
        except Exception as e:
            self.show_message("Could not open file.",e)

        try:
            self.wl_calibration[0] = np.polyfit(wl,motor_pos,3)
            self.wl_calibration[1] = np.polyfit(motor_pos,wl,3)
        except Exception as e:
            self.show_message("could not make polynomial fit for data.",e)

        return

    def update_treated_data(self, result_dict, experiment_axis_name='test'):

        if result_dict is None: return

        try:

            spectra = int(np.round(len(result_dict['df'].intensity) / 2, 0))

            self.intensity_data.set_data(result_dict['df'].wavelength, result_dict['df'].intensity[spectra])
            self.intensity_fit_data.set_data(result_dict['peak_x_coord'], result_dict['exciton_fit'])
            self.experiment_data.set_data(result_dict['second_x_axis'], result_dict['center_ev'])
            # ax1.plot(angles_rad,resid(result_fss.init_vals,angles),'-g',label='start')
            self.experiment_data_fit.set_data(result_dict['second_x_axis'], result_dict['fit_data'])

            if result_dict['dev_FSS'] is not None:
                self.ui.treated_data_2.canvas.ax.text(.7, .7, "{:.3f} ± {:.3f} µeV".format(result_dict['FSS'][0],
                                                                                           result_dict['dev_FSS'])
                                                      , transform=self.ui.treated_data_2.canvas.ax.transAxes)
            else:
                self.ui.treated_data_2.canvas.ax.text(.7, .7, "{:.3f} µeV".format(result_dict['FSS'][0])
                                                      , transform=self.ui.treated_data_2.canvas.ax.transAxes)

            self.ui.treated_data.canvas.ax.set_xlim(
                [min(result_dict['df'].wavelength), max(result_dict['df'].wavelength)])
            self.ui.treated_data.canvas.ax.set_ylim(
                [min(result_dict['df'].intensity[spectra]), max(result_dict['df'].intensity[spectra]) * 1.25])
            self.ui.treated_data.canvas.ax.set_xlabel('Wavelength (nm)')
            self.ui.treated_data.canvas.ax.set_ylabel('Intensity (arb. un.)')
            self.ui.treated_data.canvas.ax.legend(loc=2)
            self.ui.treated_data.canvas.fig.tight_layout()

            self.ui.treated_data_2.canvas.ax.set_xlim(
                [min(result_dict['second_x_axis']), max(result_dict['second_x_axis'])])
            min_fss = np.mean(result_dict['center_ev']) - 1.5 * result_dict['FSS'][0] * 1e-6
            max_fss = np.mean(result_dict['center_ev']) + 1.5 * result_dict['FSS'][0] * 1e-6
            self.ui.treated_data_2.canvas.ax.set_ylim([min_fss, max_fss])
            # [min(result_dict['center_ev'])*.85, max(result_dict['center_ev']) * 1.25])
            self.ui.treated_data_2.canvas.ax.set_xlabel(experiment_axis_name)
            self.ui.treated_data_2.canvas.ax.set_ylabel('Energy (ueV)')
            self.ui.treated_data_2.canvas.ax.legend(loc=2)
            self.ui.treated_data_2.canvas.fig.tight_layout()
            self.ui.treated_data_2.canvas.draw()
            self.ui.treated_data.canvas.draw()

        except Exception as e:
            self.show_message("Could not update Data graph", e)
            print(e)

        return

    def show_message(self, message: str, error: Exception = None):

        """
        Display error message in UI and log it

        :param message: Error message to be displayed on UI as string
        :type message: str
        :param error: Exception error to be logged
        :type error: Exception, optional
        :return: No return
        """

        self.ui.statusbar.setStyleSheet("color: red;")
        self.ui.statusbar.showMessage(message, 3000)

        today = datetime.today()

        year = today.strftime("%Y")
        month = today.strftime("%m")
        day = today.strftime("%d")

        hour = today.strftime("%H")
        minute = today.strftime("%M")
        second = today.strftime("%S")

        filename = "Logg_{}_{}_{}.txt".format(year, month, day)
        absolute_path = Path.cwd()
        filename = os.path.join(absolute_path, "logfiles", filename)

        if error is not None:

            try:
                if os.path.isdir(os.path.join(absolute_path, "logfiles")):
                    pass
                else:
                    os.mkdir(os.path.join(absolute_path, "logfiles"))
            except Exception as e:
                self.ui.statusbar.showMessage("Logfiles folder does not exist and could not be created", 3000)
                return

            try:

                if not os.path.isfile(filename):
                    with open(filename, "w+") as f:
                        f.write("Error log: \n")

                with open(filename, "a") as f:
                    f.write("{}:{}:{} : \nError Message: {} \nError: {}\n"
                            .format(hour, minute, second, message, error))

            except Exception as e:
                print("Error in logging error, :(. I try to log again ", e)

                try:

                    if not os.path.isfile(filename):
                        with open(filename, "w+") as f:
                            f.write("Error log: \n")

                    with open(filename, "a") as f:
                        f.write("{}:{}:{} : \nError Message: {} \nError: {}\n"
                                .format(hour, minute, second, message, error))
                except Exception as e:
                    print("NOPE!", e)

    def check_text(self) -> bool:

        """

        Check to set that labels are set properly for number of spectras for michelson, michelson linear piezo step and michelson piezo step

        :return: No return
        """

        checked = True

        self.ui.label_spectrum.setText(
            f"{(self.ui.linear_step.value() + 1) * (self.ui.piezo_step.value() + 1):,} spectra to be taken".replace(",",
                                                                                                                    " "))

        if self.ui.force_equal_path_2.isChecked():
            if get_float(self.ui.linear_step.text()) < 4:
                self.ui.linear_stepsize.setText("Needs at least 4 steps!")

            else:
                try:
                    equal_path = get_float(self.ui.equal_path_length.text())
                    end_linear = get_float(self.ui.linear_end.text())
                    max_delay = get_float(self.ui.max_delay.text())
                    total_steps = int(get_float(self.ui.linear_step.text()))
                    step_size = get_float(self.ui.linear_stepsize.text())
                    start = get_float(self.ui.linear_start.text())

                    if total_steps > 10:
                        before_steps = 2
                    elif total_steps > 20:
                        before_steps = 3
                    else:
                        before_steps = 1

                    if np.isnan([end_linear, start, step_size]).any():
                        start = 0.0
                        end_linear = 91.0
                        total_steps = 5
                        before_steps = 1

                    if before_steps < 1: before_steps = 1

                    after_steps = total_steps - before_steps

                    step_size = (end_linear - equal_path) / after_steps

                    start = equal_path - (before_steps * step_size)

                    michelson_pos_error = start < 0 or end_linear > 300

                    while michelson_pos_error:

                        step_size = abs(step_size - .5)

                        start = equal_path - before_steps * step_size

                        michelson_pos_error = start < 0 or end_linear > 300

                        if step_size < 0.1: break

                    print(end_linear, start, step_size)

                    if np.isnan([end_linear, start, step_size]).any():
                        start = 0.0
                        end_linear = 91.0
                        step_size = 22.75

                        self.ui.linear_end.setText("91")
                        self.ui.linear_start.setText("0")
                        self.ui.linear_stepsize.setText("22.75")
                        self.ui.linear_step.setProperty("value", int(4))

                    else:
                        end_linear = equal_path + (after_steps - 1) * step_size

                    print('end linear: ', end_linear, type(end_linear))

                    self.ui.linear_end.setText(f"{end_linear:.4}")
                    self.ui.linear_stepsize.setText(f"{step_size:.4f}")
                    if step_size < 0.01:
                        self.ui.linear_start.setText("Error! Start pos <0 mm!")
                    else:
                        self.ui.linear_start.setText(f"{start:.4f}")

                except Exception as e:
                    self.show_message("Error in Manual setup of Michelson!", e)
                    print(e)
                    start_linear = 0
                    end_linear = 91
                    step_size = 22.75
                    self.ui.linear_end.setText(f"{end_linear:.4}")
                    self.ui.linear_start.setText(f"{start_linear:.4f}")
                    self.ui.linear_stepsize.setText(f"{step_size:.4f}")
                    self.ui.linear_step.setProperty("value", 4)


        elif self.ui.force_equal_path.isChecked():

            step_size = get_float(self.ui.linear_stepsize.text())
            start = get_float(self.ui.linear_start.text())
            equal_path = get_float(self.ui.equal_path_length.text())
            end_linear = get_float(self.ui.linear_end.text())
            max_delay = get_float(self.ui.max_delay.text())
            if end_linear > max_delay * 1.5e2: linear_end = max_delay * 1.5e2

            total_steps = int(get_float(self.ui.linear_step.text()))

            try:
                if get_float(self.ui.linear_step.text()) < 4:
                    self.ui.linear_stepsize.setText("Needs at least 4 steps!")

                else:

                    if total_steps > 10:
                        before_steps = 2
                    elif total_steps > 20:
                        before_steps = 3
                    else:
                        before_steps = 1

                    if np.isnan([end_linear, start, step_size]).any():
                        start = 0.0
                        end_linear = 91.0
                        total_steps = 5
                        before_steps = 1

                    if before_steps < 1: before_steps = 1

                    after_steps = total_steps - before_steps

                    step_size = (end_linear - equal_path) / after_steps

                    start = equal_path - (before_steps * step_size)

                    michelson_pos_error = start < 0

                    while michelson_pos_error:

                        step_size = abs(step_size - .5)

                        start = equal_path - before_steps * step_size

                        michelson_pos_error = start < 0

                        if step_size < 0.01: break

                    end_linear = equal_path + (after_steps - 1) * step_size

                    self.ui.linear_end.setText(f"{end_linear:.4}")
                    self.ui.linear_stepsize.setText(f"{step_size:.4f}")
                    if step_size < 0.01:
                        self.ui.linear_start.setText("Error! Start pos <0 mm!")
                    else:
                        self.ui.linear_start.setText(f"{start:.4f}")

            except ValueError:
                self.ui.linear_stepsize.setText("Error!")
                checked = False

            try:
                self.ui.piezo_stepsize.setText(
                    f"{(get_float(self.ui.piezo_end.text()) - get_float(self.ui.piezo_start.text())) / (self.ui.piezo_step.value()):.4f}")
            except ValueError:
                self.ui.piezo_stepsize.setText("Error")
                checked = False
        else:

            try:
                self.ui.linear_stepsize.setText(
                    f"{(get_float(self.ui.linear_end.text()) - get_float(self.ui.linear_start.text())) / (self.ui.linear_step.value()):.4f}")
            except ValueError:
                self.ui.linear_stepsize.setText("")
                checked = False

            try:
                self.ui.piezo_stepsize.setText(
                    f"{(get_float(self.ui.piezo_end.text()) - get_float(self.ui.piezo_start.text())) / (self.ui.piezo_step.value()):.4f}")
            except ValueError:
                self.ui.piezo_stepsize.setText("")
                checked = False

        return checked

    def copy_filename(self):

        """
        Copy the filename set into the label to the clipboard

        :return: No return
        """

        self.clipboard.setText(self.ui.label_filename.text().replace(
            ".spe", "").replace("filename", ""))

    def copy_number_of_spectra(self):

        """
        Copy the michelson number of spectras set into the label to the clipboard

        :return: No return
        """

        self.clipboard.setText(
            "".join(self.ui.label_spectrum.text().split(" ")[:-4]))

    def show_error_msg(self, device: str):

        """
        Display Error windows message and disconnect the device who failed.
        Also sends message on GUI and logs the error.

        :param device: Device who failed
        :return: None
        """

        QMessageBox.critical(
            self.main_window, "Device disconnected", f"{device} disconnected")
        self.show_message("Error in {}. Disconnecting.".format(device))
        self.disconnect(dev_list=[device])

    def startup_dialog(self):

        """
        Opens a dialog box at software startup to select z-focus piezo, and if spectrometer/kethley should be connected


        :return: No return
        """

        self.dialog.exec()

        self.piezo_name = "None"

        if self.dialog.piezo.currentText() == "None":
            self.piezo_name = "None"
        elif self.dialog.piezo.currentText() == "PL Piezo":
            self.piezo_name = "focus_piezo"

        if self.piezo_name != "None": self.device_name_list.insert(0, self.piezo_name)

        self.ui.label_trigger_5.setText(self.dialog.piezo.currentText() + ": ")

        if self.dialog.spectrometer_connect.isChecked(): self.device_name_list.insert(0, self.spectrometer_name)
        if self.dialog.kethley_connect.isChecked(): self.device_name_list.append(self.smu_name)
        if self.dialog.cw_3900s.isChecked():
            self.device_name_list.append(self.cw_name)
            self.device_name_list.append(self.cw_hwp_name)
        if self.dialog.power_meter.isChecked(): self.device_name_list.append(self.pwr_meter_name)
        if self.dialog.debug_mode.isChecked(): self.device_name_list = []

        print(self.device_name_list)

    # Connection/disconnect functions

    def connect(self, dev_list: list = None):

        """

        Create a communication instance with the device servers required (defined by requested devices)

        :params dev_list: list of devices to connect
        :type dev_list: list
        :return: self.com, instance of communication variable with device server
        :return: self.com, instance of communication variable with device server
        """

        device_name_list = self.device_name_list

        if dev_list is not None and type(dev_list) is not bool:
            device_name_list = dev_list

        # Easy communication interface
        if self.com is None:
            print("starting connection thread")
            try:
                self.com = ReqParser("pl_control", device_name_list)
            except Exception as e:
                self.show_message("Disconnecting devices", e)
        else:
            self.show_message("Disconnecting devices")

            for dev in device_name_list:
                try:
                    self.com.request(dev, "disconnect")
                except Exception as e:
                    self.com = None
                    self.show_message("Disconnecting devices", e)

            print("starting connection thread")

            connect_tries=0
            while self.com == None and connect_tries<5:
                self.com = ReqParser("pl_control", device_name_list)
                # self.show_message("Error in connection.", e)
                connect_tries+=1

            if connect_tries>=5:
                self.show_message("Error in connection.", e)

        try:
            self.connect_thread = Connect(self.com, device_name_list)
            self.connect_thread.message.connect(self.show_message)
            # self.connect_thread.connected.connect(self.connected)
            self.connect_thread.connect_message.connect(self.show_message)
            self.connect_thread.start()
            self.connect_thread.finished.connect(lambda: self.connected(True))
            print("reached connection")

        except Exception as e:

            self.show_message("Error in connection.", e)

        return self.com

        # self.update_current_positions("startup")

    def connected(self, status: bool):

        """
            Check the connection status for each device after connected thread is finished (or on request)


        :param status: Unused by now
        :return: No return
        """

        if not status:
            self.show_message("Error in connected? {}".format(str(status)))

        try:
            linear_connected, self.curr_linear_pos = self.com.request(self.lin_stage_name, "get_pos")
        except Exception as e:
            self.show_message("Michelson linear stage not connected.", e)
            linear_connected = False

        try:
            piezo_connected, self.curr_mich_piezo_pos = self.com.request(self.piezo_stage_name, "get_pos")
        except Exception as e:
            self.show_message("Michelson Piezo not connected.", e)
            piezo_connected = False

        try:
            cryo_connected, (self.curr_x, self.curr_y) = self.com.request(self.xy_stage_name, "get_pos")
        except Exception as e:
            self.show_message("XY Stage not connected.", e)
            cryo_connected = False

        try:
            focus_connected, (self.curr_z) = self.com.request(self.piezo_name, "get_pos")
        except Exception as e:
            self.show_message("Z-Piezo Not connected.", e)
            focus_connected = False

        try:
            fss_connected, self.curr_fss = self.com.request(self.fss_name, "get_pos")
        except Exception as e:
            self.show_message("Collection HWP not connected.", e)
            fss_connected = False

        try:
            laser_connected, self.curr_laser = self.com.request(self.laser_power_name, "get_pos")
        except Exception as e:
            self.show_message("Laser HWP not connected.", e)
            laser_connected = False

        if self.smu_name in self.device_name_list:
            try:
                kethley_connected, (self.curr_V, self.curr_I) = self.com.request(self.smu_name, "apply_voltage", 0)
                self.ui.curr_V.setText(str(np.round(get_float(self.curr_V), 4)))
                self.ui.curr_I.setText(str(np.round(get_float(self.curr_I)*1e6, 4)))
            except Exception as e:
                self.show_message("Kethley not connected", e)
                kethley_connected = False
        else:
            kethley_connected = False
            self.ui.curr_V.setText(str("Not connected"))
            self.ui.curr_I.setText(str("Not connected"))

        if self.cw_name in self.device_name_list:
            try:
                cw_3900s_connected, self.curr_cw[0] = self.com.request(self.cw_name, "set_pos", 0.1)
            except Exception as e:
                self.show_message("General Error message", e)
                cw_3900s_connected = False
        else:
            cw_3900s_connected = False

        if self.spectrometer_name in self.device_name_list:
            try:
                self.spectrometer_connected = self.com.request(self.spectrometer_name, "test")
            except Exception as e:
                self.show_message("General Error message", e)
                self.spectrometer_connected = False
        else:
            self.spectrometer_connected = False

        waveplate_velocity = 20  # deg/s maximum/standard 20deg/s

        self.toggle_all(True)

        try:
            self.thr_bus_state = True

        except Exception as e:
            self.show_message(str("thr elliptec not connected: " + str(self.thr_bus_state)), e)
            self.thr_bus_state = False

        if self.thr_bus_state:
            try:
                self.camera_piezo_state, _ = self.com.request(self.thr_bus_controller, "get_pos",
                                                              {self.camera_piezo_name: None})

            except Exception as e:
                self.show_message("Could not get pos from elliptec piezo for camera piezo.", e)
                self.camera_piezo_state = False

            try:
                self.michelson_piezo_state,_ = self.com.request(self.thr_bus_controller, "get_pos",
                                                              {self.michelson_piezo_name: None})
            except Exception as e:
                self.show_message("Could not get pos from elliptec piezo for michelson piezo", e)
                self.camera_piezo_state = False

        if cryo_connected:
            self.ui.status_cryo.setText("Connected")
            self.ui.xy_movement.setEnabled(True)
            self.ui.space_map_group.setEnabled(True)
            self.ui.line_scan_group.setEnabled(True)
            self.ui.stored_positions_3.setEnabled(True)
        else:
            self.ui.status_cryo.setText("Not connected")

        if focus_connected:
            self.ui.status_focus_piezo.setText("Connected")
            self.ui.z_movement.setEnabled(True)
        else:
            self.ui.status_focus_piezo.setText("Not connected")

        if cryo_connected and focus_connected:
            self.ui.stored_positions_2.setEnabled(True)
            self.ui.fit_plane.setEnabled(True)
            self.ui.fast_map.setEnabled(True)
            self.ui.snake_like.setEnabled(True)

        if fss_connected:
            self.ui.status_fss.setText("Connected")
            if not status:
                self.ui.status_fss.setText("Not connected")
            self.ui.polarization_control.setEnabled(True)
        else:
            self.ui.status_fss.setText("Not connected")

        if laser_connected:
            self.ui.status_laser.setText("Connected")
            if not status:
                self.ui.status_laser.setText("Not connected")
            self.ui.power_mapping.setEnabled(True)
        else:
            self.ui.status_laser.setText("Not connected")

        if kethley_connected:
            self.ui.status_kethley.setText("Connected")
            self.ui.kethley_box.setEnabled(True)
        else:
            self.ui.status_kethley.setText("Not connected")

        if cw_3900s_connected:
            self.ui.status_3900s.setText("Connected")
            self.ui.cw_3900s_tab.setEnabled(True)
        else:
            self.ui.status_3900s.setText("Not connected")

        if self.spectrometer_connected:
            experiment = str(self.ui.experiment_name.text()).strip()
            self.com.request(self.spectrometer_name, "load_settings", experiment)
            time.sleep(0.1)
            self.ui.status_spectrometer.setText("Connected")
        else:
            self.ui.status_spectrometer.setText("Not connected")

        if self.camera_piezo_state:
            self.ui.status_camera_piezo.setText("Connected")
            self.ui.camera_bs_box.setEnabled(True)
        else:
            self.ui.status_camera_piezo.setText("Not connected")

        if self.michelson_piezo_state:
            self.ui.status_piezo_michelson.setText("Connected")
            self.ui.michelson_miror_box.setEnabled(True)
        else:
            self.ui.status_piezo_michelson.setText("Not connected")

        if linear_connected:
            self.ui.linear_cur_position.setText(f"{self.curr_linear_pos:.4} mm")
            self.ui.status_linear.setText("Connected")
            self.ui.status_trigger.setText("Connected")
            self.ui.devices_box.setEnabled(True)
            self.ui.linear_goto.returnPressed.connect(self.goto_linear)

            # set to current values
            self.set_linear_velocity()
            self.set_linear_acceleration()
            self.set_linear_deceleration()
        else:
            self.ui.go_button.setEnabled(False)
            self.ui.linear_loop.setEnabled(False)
            self.ui.status_linear.setText("Not connected")
            self.ui.status_trigger.setText("Not connected")

        if piezo_connected:
            self.ui.piezo_cur_position.setText(f"{self.curr_mich_piezo_pos:.4} µm")
            self.ui.status_piezo.setText("Connected")
            self.ui.piezo_goto.returnPressed.connect(self.goto_piezo)
        else:
            self.ui.go_button_2.setEnabled(False)
            self.ui.piezo_loop.setEnabled(False)
            self.ui.status_piezo.setText("Not connected")

        if linear_connected and piezo_connected:
            self.ui.label_connect_devices.setText("")
        else:
            self.ui.label_connect_devices.setText("need to connect all devices!")
            self.ui.measurement_button.setEnabled(False)

        self.update_current_positions(axis=None)

    def disconnect(self, dev_list: list = None):

        """

        Disconnect each device in the device list

        :params dev_list: list of devices, default = None
        :type dev_list: list
        :return: No return
        """

        if dev_list is None:
            dev_list = self.device_name_list
        

        
        try:
            if self.com is not None:
                info = None
                i = 0
                while len(dev_list)>0 and i<50:
                    for dev in dev_list:
                        status = self.com.request(dev, "disconnect")
                        if len(status)>1:
                            status,info = status[0],status[1]
                        if not status:
                            self.show_message("Could not disconect from {}, due to {}".format(dev,info))
                        else:
                            dev_list.pop(dev_list.index(dev))
                            i+=1
                        if i==50:
                            raise Exception("Too many tries to disconnet! {} devices where not disconnected.".format(dev_list))
                    # self.toggle_measurement_ui(False)
                self.show_message("Devices disconnected.")
        except Exception as e:
            self.show_message("Error in disconnection.",e)

        return

    # =============================================================================
    # Movement/position storage functions
    # =============================================================================

    def update_text_axis(self, axis: str):

        """
        Update text axis in the UI

        :param axis: Axis to be updated
        :type axis: str
        :return: No return
        """

        if axis in ('x', 'y', 'xy'):

            if self.ui.invert_xy.isChecked():
                if self.curr_x is str:
                    self.ui.curr_y.setText('Error')
                else:
                    self.curr_x = np.round(self.curr_x, 4)
                    self.ui.curr_y.setText(str(self.curr_x))

                if self.curr_y is str:
                    self.ui.curr_x.setText('Error')
                else:
                    self.curr_y = np.round(self.curr_y, 4)
                    self.ui.curr_x.setText(str(self.curr_y))

            else:
                if self.curr_x is str:
                    self.ui.curr_x.setText('Error')
                else:
                    self.curr_x = np.round(self.curr_x, 4)
                    self.ui.curr_x.setText(str(self.curr_x))

                if self.curr_y is str:
                    self.ui.curr_y.setText('Error')
                else:
                    self.curr_y = np.round(self.curr_y, 4)
                    self.ui.curr_y.setText(str(self.curr_y))
        elif axis == 'xyz':

            if self.ui.invert_xy.isChecked():
                if self.curr_x is str:
                    self.ui.curr_y.setText('Error')
                else:
                    self.curr_x = np.round(self.curr_x, 4)
                    self.ui.curr_y.setText(str(self.curr_x))

                if self.curr_y is str:
                    self.ui.curr_x.setText('Error')
                else:
                    self.curr_y = np.round(self.curr_y, 4)
                    self.ui.curr_x.setText(str(self.curr_y))

            else:
                if self.curr_x is str:
                    self.ui.curr_x.setText('Error')
                else:
                    self.curr_x = np.round(self.curr_x, 4)
                    self.ui.curr_x.setText(str(self.curr_x))

                if self.curr_y is str:
                    self.ui.curr_y.setText('Error')
                else:
                    self.curr_y = np.round(self.curr_y, 4)
                    self.ui.curr_y.setText(str(self.curr_y))
            if self.curr_z is str:
                self.ui.curr_z.setText('Error')
            else:
                self.curr_z = np.round(float(self.curr_z), 4)
                self.ui.curr_z.setText(str(self.curr_z))
        elif axis == 'z':
            if self.curr_z is str:
                self.ui.curr_z.setText('Error')
            else:
                self.curr_z = np.round(float(self.curr_z), 4)
                self.ui.curr_z.setText(str(self.curr_z))
        elif axis == 'fss':
            if self.curr_fss is str:
                self.ui.curr_fss.setText('Error')
            else:
                self.curr_fss = np.round(float(self.curr_fss), 4)
                self.ui.curr_fss.setText(str(np.round(float(self.curr_fss) - float(self.fss_offset), 4)))
        elif axis == 'laser':
            if self.curr_laser is str:
                self.ui.curr_laser.setText('Error')
            else:
                self.curr_laser = np.round(float(self.curr_laser), 4)
                self.ui.curr_laser.setText(str(np.round(float(self.curr_laser) - float(self.laser_offset), 4)))
        elif axis == 'cw_hwp':
            try:
                self.curr_cw_hwp = np.round(float(self.curr_cw_hwp), 4)
                self.ui.curr_cw_hwp.setText(str(np.round(float(self.curr_cw_hwp) - float(self.cw_laser_offset), 4)))
            except Exception as e:
                self.show_message("Error in setting hwp angle.",e)
                self.ui.curr_cw_hwp.setText('Error')


        elif axis in ('cw_laser', 'cw_laser_wl'):

            if self.curr_cw[0] is str:
                self.ui.curr_cw_pos.setText('Error')
                return
            else:
                self.ui.curr_cw_pos.setText(str(round(float(self.curr_cw[0]), 4)))
            if self.curr_cw[1] is str :
                self.ui.curr_cw_wl.setText('Error')
                return

            try:
                self.curr_cw[1] = self.wl_regression(self.curr_cw[0])
            except Exception as e:
                self.ui.curr_cw_wl.setText(str("Error!"))
                return
            if self.curr_cw[1] is None or self.curr_cw[1] == "Error":
                self.ui.curr_cw_wl.setText(str("Error!"))
                return
            try:
                self.curr_cw[0] = np.round(float(self.curr_cw[0]), 4)
                self.curr_cw[1] = np.round(float(self.curr_cw[1]), 4)
                self.ui.curr_cw_pos.setText(str(self.curr_cw[0]))
                self.ui.curr_cw_wl.setText(str(self.curr_cw[1]))
            except Exception as e:
                self.ui.curr_cw_pos.setText('Error')
                self.ui.curr_cw_wl.setText('Error')
                self.show_message("Could not set pos/wl for cw laser!", e)

        return

    def update_current_positions(self, axis: str = None):

        """
        Update current positions in the UI according to the given axis.
        Default axis value is None which request the update of all axis.

        :param axis: Axis to be updated or None for all axis.
        :return: No return
        """

        if isinstance(type(axis), type(None)):
            axis = None

        if axis is None or isinstance(type(axis), type(None)) or axis == "startup":

            try:
                status, (self.curr_x, self.curr_y) = self.move_request(self.xy_stage_name, "get_pos")
                self.update_text_axis('xy')
            except Exception as e:
                self.show_message("Error in getting XY motors position, pass", e)
                pass
            try:
                status, self.curr_z = self.move_request(self.piezo_name, "get_pos")
                self.update_text_axis('z')
            except Exception as e:
                self.show_message("Error in getting Z-piezo position, pass", e)
                pass
            try:
                status, self.curr_fss = self.move_request(self.fss_name, "get_pos")
                self.update_text_axis('fss')
            except Exception as e:
                self.show_message("Error in getting collection HWP position, pass", e)
                pass
            try:
                status, self.curr_laser = self.move_request(self.laser_power_name, "get_pos")
                self.update_text_axis('laser')
            except Exception as e:
                self.show_message("Error in getting laser HWP position, pass", e)
                pass
            try:
                status, self.curr_cw[0] = self.move_request(self.cw_name, "set_pos",0.1)
                self.update_text_axis('cw_laser')
            except Exception as e:
                self.show_message("Error in getting cw motor position, pass", e)

        else:
            if axis in ('x', 'y'):
                try:
                    status, (self.curr_x, self.curr_y) = self.move_request(self.xy_stage_name, "get_pos")
                    self.update_text_axis(axis)
                except Exception as e:
                    self.show_message("Error in getting XY motors position, pass", e)
                    pass
            if axis == "xyz":
                try:
                    status, (self.curr_x, self.curr_y) = self.move_request(self.xy_stage_name, "get_pos")
                    self.update_text_axis("xy")
                except Exception as e:
                    self.show_message("Error in getting XY motors position, pass", e)
                    pass
                try:
                    status, self.curr_z = self.move_request(self.piezo_name, "get_pos")
                    self.update_text_axis("z")
                except Exception as e:
                    self.show_message("Error in getting Z-piezo position, pass", e)
                    pass
            if axis == "z":
                try:
                    status, self.curr_z = self.move_request(self.piezo_name, "get_pos")
                    self.update_text_axis(axis)
                except Exception as e:
                    self.show_message("Error in getting Z-piezo position, pass", e)
                    pass
            elif axis == 'fss':
                try:
                    status, self.curr_fss = self.move_request(self.fss_name, "get_pos")
                    self.update_text_axis(axis)
                except Exception as e:
                    self.show_message("Error in getting collection HWP position, pass", e)
                    pass
            elif axis == 'laser':
                try:
                    status, self.curr_laser = self.move_request(self.laser_power_name, "get_pos")
                    self.update_text_axis(axis)
                except Exception as e:
                    self.show_message("Error in getting laser HWP position, pass", e)
                    pass
            elif axis == 'kethley':
                try:
                    self.update_iv_label((self.curr_V, self.curr_I))
                except Exception as e:
                    self.show_message("Error in getting Kethley Voltage/Current, pass", e)
                    pass
            elif axis == 'cw_hwp':
                try:
                    status, self.curr_cw_hwp = self.move_request(self.cw_hwp_name, "get_pos")
                    self.update_text_axis('cw_hwp')
                except Exception as e:
                    self.show_message("Error in getting cw hwp position, pass", e)

            elif axis in ('cw_laser', 'cw_laser_wl'):
                try:
                    status, self.curr_cw[0] = self.move_request(self.cw_name, "get_pos")
                    self.update_text_axis('cw_laser')
                except Exception as e:
                    self.show_message("Error in getting cw motor position, pass", e)

        print('X: ', self.curr_x, 'Y: ', self.curr_y, 'Z: ', self.curr_z, 'fss: ', self.curr_fss, 'laser: ',
              self.curr_laser)

        return

    def update_offsets(self):

        """
        Get offsets for fss and laser from UI values

        :return: No return
        """

        self.fss_offset = get_float(self.ui.pol_offset.text())
        self.laser_offset = get_float(self.ui.laser_offset.text())
        self.cw_laser_offset = get_float(self.ui.cw_laser_offset.text())
        return

    # Store and remove positions from tables

    def import_export_pos(self, load_from_file: bool):

        today = date.today()

        year = today.strftime("%Y")
        month = today.strftime("%m")
        day = today.strftime("%d")

        open_path = os.path.join(self.path, year, month, str(day))  # , str(self.ui.sample_name.text()).strip())

        path = file_dialog('Choose positions File Location', load_from_file, fmt='csv', open_path=open_path)

        if path is not None:
            try:
                if load_from_file:
                    self.ui.stored_positions.clear()
                    with open(path) as csvfile:
                        reader = csv.reader(csvfile)
                        header = next(reader)
                        self.ui.stored_positions.setColumnCount(len(header))
                        self.ui.stored_positions.setHorizontalHeaderLabels(header)
                        for row, values in enumerate(reader):
                            self.ui.stored_positions.insertRow(row)
                            for column, value in enumerate(values):
                                self.ui.stored_positions.setItem(
                                    row, column, QtWidgets.QTableWidgetItem(value))


                else:
                    columns = range(self.ui.stored_positions.columnCount())
                    header = [self.ui.stored_positions.horizontalHeaderItem(column).text()
                              for column in columns]
                    with open(path, 'w') as csvfile:
                        writer = csv.writer(
                            csvfile, dialect='excel', lineterminator='\n')
                        writer.writerow(header)
                        for row in range(self.ui.stored_positions.rowCount()):
                            try:
                                writer.writerow(
                                    self.ui.stored_positions.item(row, column).text()
                                    for column in columns)
                            except:
                                pass

            except Exception as e:
                self.show_message("Error in import/exporting positions", e)

    def store_position(self, table: str):

        """
        Store current position in a table

        :param table: Table to be stored, stored_positions for simple movements,
        or plane_position for meshed XYZ plane fit
        :type table: str
        :return: No return

        """

        status, (self.curr_x, self.curr_y) = self.com.request(self.xy_stage_name, "get_pos")
        status, self.curr_z = self.com.request(self.piezo_name, "get_pos")

        self.curr_x = str(np.round(self.curr_x, 4))
        self.curr_y = str(np.round(self.curr_y, 4))
        self.curr_z = str(np.round(self.curr_z, 4))


        if table == 'stored_positions':
            table = self.ui.stored_positions
        elif table == 'plane_positions':
            table = self.ui.plane_positions

        rowposition = table.rowCount()

        if len(table.selectedIndexes()) > 0:
            if table.selectedIndexes() not in (None, 'Nan', '', 'None', 'Null'):
                rowposition = table.selectedIndexes()[0].row()

        if table == self.ui.stored_positions:
            table.insertRow(rowposition)
            table.setItem(rowposition, 1, QTableWidgetItem(self.curr_x))
            table.setItem(rowposition, 2, QTableWidgetItem(self.curr_y))
            table.setItem(rowposition, 3, QTableWidgetItem(self.curr_z))
        elif table == self.ui.plane_positions:
            table.insertRow(rowposition)
            table.setItem(rowposition, 0, QTableWidgetItem(self.curr_x))
            table.setItem(rowposition, 1, QTableWidgetItem(self.curr_y))
            table.setItem(rowposition, 2, QTableWidgetItem(self.curr_z))

        return

    def remove_position(self, table: str):

        """
        Remove current position from the table

        :param table: table from where the position is to be removed. stored_positions or plane_positions
        :type table: str
        :return: No return
        """

        if table == 'stored_positions':
            table = self.ui.stored_positions
        elif table == 'plane_positions':
            table = self.ui.plane_positions

        irows = []

        if table.selectedIndexes() not in (None, 'Nan', '', 'None', 'Null'):
            for obj in table.selectedIndexes():
                irow = obj.row()
                if irow not in irows:
                    irows.append(irow)

        irows.sort()

        for irow in reversed(irows):
            if table.verticalHeaderItem(irow):
                pass
            else:
                table.removeRow(irow)

        return

    def goto_store_position(self, table: str):

        """
        Remove current position from the table

        :param table: table from where the position is to go is taken. stored_positions or plane_positions
        :type table: str
        :return: No return
        """

        if table == 'stored_positions':
            table = self.ui.stored_positions
        elif table == 'plane_positions':
            table = self.ui.plane_positions

        rowposition = None

        if len(table.selectedIndexes()) > 0:
            if table.selectedIndexes() not in (None, 'Nan', '', 'None', 'Null'):
                # for obj in self.ui.custom_structure.selectedIndexes():
                rowposition = table.selectedIndexes()[0].row()
        else:
            return

        if rowposition is None:
            self.show_message("No position selected")
            return

        if not table.verticalHeaderItem(rowposition):
            if table == self.ui.plane_positions:
                target_x = get_float(table.item(rowposition, 0).data(0))
                target_y = get_float(table.item(rowposition, 1).data(0))
                target_z = get_float(table.item(rowposition, 2).data(0))
            elif table == self.ui.stored_positions:
                target_x = get_float(table.item(rowposition, 1).data(0))
                target_y = get_float(table.item(rowposition, 2).data(0))
                target_z = get_float(table.item(rowposition, 3).data(0))

            self.make_move(0, 'xyz', 0, goto=[target_x, target_y, target_z])

        return

    def make_move(self, multiplier: float, axis: str, direction: float, goto: float | list | tuple = None):

        """

        Move motors/stages to a relative position (if goto==None)
        final position = current+(multiplier*direction)*1 (relative un)
        or to an absolute position by using the goto variable

        :param multiplier: value by which the change is multiplied, standard value is 1*multiplier
        :param axis: axis to be moved, x,y,z,xyz, fss, laser,
        :param direction: direction in which the movement is made, either 1 or -1
        :param goto: goto absolute position, float, (x,y), or (x,y,z)
        :return: No return
        :type multiplier: float
        :type axis: str
        :type direction: float
        :type goto: float,list or tuple
        """

        target = 0

        # start_time = time.time()

        if self.ui.invert_xy.isChecked():
            if axis == 'y':
                axis = 'x'
            elif axis == 'x':
                axis = 'y'

        if self.ui.invert_x.isChecked():
            if axis == 'x':
                multiplier = -1 * multiplier
        if self.ui.invert_y.isChecked():
            if axis == 'y':
                multiplier = -1 * multiplier
        if axis == "laser" and goto is not None:
            goto = float(goto) + float(self.laser_offset)
        if axis == "fss" and goto is not None:
            goto = float(goto) + float(self.fss_offset)
        if axis == "cw_hwp" and goto is not None:
            goto = float(goto) + float(self.cw_laser_offset)

        xy_min_limit = 0
        xy_max_limit = 50
        z_min_limit = 0
        z_max_limit = 100
        cw_limit = 25
        cw_wl_limit = (600, 1100)

        direction = int(direction)
        multiplier = float(multiplier)

        # check for safe positions
        if goto is not None:
            if axis in ('x', 'y'):
                if goto < xy_min_limit or goto > xy_max_limit:
                    return
            elif axis == 'xy':
                if isinstance(goto, float) or len(goto) < 2:
                    self.show_message("Not enough coordinates!")
                    return
                elif goto[0] < xy_min_limit or goto[0] > xy_max_limit or \
                        goto[1] < xy_min_limit or goto[1] > xy_max_limit:
                    self.show_message("Tried to move XY stage out of range X:{},Y:{}".format(goto[0], goto[1]))
                    return
            elif axis == 'xyz':
                if isinstance(goto, float) or len(goto) < 3:
                    self.show_message("Not enough coordinates!")
                    return
                elif goto[0] < xy_min_limit or goto[0] > xy_max_limit or \
                        goto[1] < xy_min_limit or goto[1] > xy_max_limit or \
                        goto[2] < z_min_limit or goto[2] > z_max_limit:
                    self.show_message("Tried to move XY stage out of range X:{},Y:{},Z:" \
                                      .format(goto[0], goto[1], goto[2]))
                    return
            elif axis == 'z':
                if goto < z_min_limit or goto > z_max_limit:
                    return

            elif axis == 'cw_laser':
                if abs(goto) > cw_limit:
                    return
            elif axis == 'cw_laser_wl':
                if abs(goto) < cw_wl_limit[0] or abs(goto) > cw_wl_limit[1]:
                    return

            elif axis == 'cw_hwp':
                if abs(goto) < 0 or abs(goto) > 360:
                    goto = goto % 360

        else:
            if axis in ('x', 'y'):
                if axis == 'x':
                    target = float(self.curr_x) + get_float(
                        self.ui.step_size_xy.text()) * multiplier * direction * 1e-3
                elif axis == 'y':
                    target = float(self.curr_y) + get_float(
                        self.ui.step_size_xy.text()) * multiplier * direction * 1e-3

                if target < xy_min_limit or target > xy_max_limit:
                    self.show_message(
                        str("Tentative to move XY stage above maximum range: " + str(axis) + ": " + str(target)))
                    return

            elif axis == 'z':
                target_z = float(self.curr_z) + get_float(self.ui.step_size_z.text()) * multiplier
                if target_z < z_min_limit or target_z > z_max_limit:
                    return

            elif axis == 'fss':

                status, self.curr_fss = self.move_request(self.fss_name, "get_pos")

                if not status:
                    return

            elif axis == 'laser':

                status, self.curr_laser = self.move_request(self.laser_power_name, "get_pos")

                if not status:
                    return

            elif axis in ('cw_laser', 'cw_laser_wl'):
                status, self.curr_cw[0] = self.move_request(self.cw_name, "get_pos")

                if not status:
                    return

            elif axis == 'cw_hwp':
                status, self.curr_cw_hwp = self.move_request(self.cw_hwp_name, "get_pos")

        # actual movement happens below here
        if goto is not None:
            if axis == 'x':
                goto = goto

                """status, _=self.com.request(
                    self.xy_stage_name, "set_pos", f"{goto};{ypos}")"""

                status, (self.curr_x, self.curr_y) = self.move_request(self.xy_stage_name, "set_pos",
                                                                      f"{goto};{self.curr_y}")

                if not status:
                    return

            elif axis == 'y':

                status, (self.curr_x, self.curr_y) = self.move_request(self.xy_stage_name, "set_pos",
                                                                      f"{self.curr_x};{goto}")

            elif axis == 'xy':
                if isinstance(goto, float) or len(goto) < 2:
                    self.show_message("Not enough coordinates!")
                    return

                goto_x = goto[0]
                goto_y = goto[1]

                status, (self.curr_x, self.curr_y) = self.move_request(self.xy_stage_name, "set_pos",
                                                                      f"{goto_x};{goto_y}")

            elif axis == 'xyz':
                if isinstance(goto, float) or len(goto) < 3:
                    self.show_message("Not enough coordinates!")
                    return

                goto_x = goto[0]
                goto_y = goto[1]
                goto_z = goto[2]

                status, (self.curr_x, self.curr_y) = self.move_request(self.xy_stage_name, "set_pos",
                                                                      f"{goto_x};{goto_y}")
                status, self.curr_z = self.move_request(self.piezo_name, "set_pos", goto_z)

            elif axis == 'z':
                status, self.curr_z = self.move_request(self.piezo_name, "set_pos", goto)

            elif axis == 'laser':
                status, self.curr_laser = self.move_request(self.laser_power_name, "set_pos", goto)
            elif axis == 'fss':
                status, self.curr_fss = self.move_request(self.fss_name, "set_pos", goto)
            elif axis in ('cw_laser', 'cw_laser_wl'):

                if axis == 'cw_laser_wl':
                    self.curr_cw[1] = goto
                    try:
                        goto = self.wl_regression(self.curr_cw[1])
                    except:
                        goto = self.curr_cw[1]

                status, self.curr_cw[0] = self.move_request(self.cw_name, "set_pos", goto)
            elif axis == 'cw_hwp':
                status, self.curr_cw_hwp = self.move_request(self.cw_hwp_name, "set_pos", goto)


        # if goto == None
        elif axis in ('x', 'y'):

            step = get_float(self.ui.step_size_xy.text()) * multiplier * direction * 1e-3

            if axis == 'x':
                xpos = float(self.curr_x) + step
                ypos = float(self.curr_y)
            else:
                xpos = float(self.curr_x)
                ypos = float(self.curr_y) + step

            for pos in (xpos, ypos):
                if pos <= xy_min_limit:
                    pos = xy_min_limit
                elif pos >= xy_max_limit:
                    pos = xy_max_limit

            status, (self.curr_x, self.curr_y) = self.move_request(self.xy_stage_name, "set_pos", f"{xpos};{ypos}")

            if not status:
                return

        elif axis == 'z':

            step = get_float(self.ui.step_size_z.text()) * multiplier
            zpos = float(self.curr_z) + direction * step

            if zpos <= z_min_limit:
                zpos = z_min_limit
            elif zpos >= z_max_limit:
                zpos = z_max_limit

            status, self.curr_z = self.move_request(self.piezo_name, "set_pos", zpos)

            if not status:
                return

        elif axis == 'fss':

            step = multiplier
            pos = float(self.curr_fss) + direction * get_float(step)
            status, self.curr_fss = self.move_request(self.fss_name, "set_pos", pos)

            if not status:
                return

        elif axis == 'laser':

            step = multiplier
            pos = float(self.curr_laser) + direction * get_float(step)
            status, self.curr_laser = self.move_request(self.laser_power_name, "set_pos", pos)

            if not status:
                return

        elif axis == 'cw_hwp':
            step = multiplier
            pos = float(self.curr_cw_hwp) + direction * get_float(step)
            print(pos, step, direction, self.curr_cw_hwp)
            status, self.curr_cw_hwp = self.move_request(self.cw_hwp_name, "set_pos", pos)

        elif axis in ('cw_laser', 'cw_laser_wl'):

            if axis == 'cw_laser_wl':
                step = multiplier
                try:
                    pos = self.wl_regression(get_float(self.curr_cw[1]) + direction * get_float(step))
                except Exception as e:
                    pos = float(self.curr_cw[1])
                    self.show_message("Could not find motor pos.", e)
                if pos is None or pos == "Error":
                    pos = float(self.curr_cw[0])
                else:
                    status, self.curr_cw[0] = self.move_request(self.cw_name, "set_pos", pos)

            else:
                step = multiplier
                pos = float(self.curr_cw[0]) + direction * get_float(step)
                status, self.curr_cw[0] = self.move_request(self.cw_name, "set_pos", pos)

            if not status:
                return

        # end_time = time.time() - start_time
        # self.show_message("movement of {} took {} s".format(axis,end_time))

        self.update_text_axis(axis)

        return

    def move_request (self, device_name,action, pos=None):

        """
        Starts measurement thread with the Measurement_1D function

        Updates position after each spectrum

        When finished, returns to initial positions and resets spectrometer to default

        :return: No return
        """

        status, pos = None,None

        def handle_result(dev_return_list):

            status, pos = dev_return_list

            if not status:
                self.show_message("Error in movement status: {}, pos {}".format(dev_return_list[0],dev_return_list[1]),Exception("Generic exception"))

            return

        if self.check_text():
            try:
                if self.movement_thread.isRunning():
                    self.show_message("Movement Thread Busy")
                    return
            except Exception as e:
                self.show_message("Error in Movement with {}".format(str(device_name)), e)

            curr_return = curr_status = None

            self.movement_thread = Movement_request(self.com,device_name,action, pos)

            self.measurement_thread.dev_return.connect(handle_result)

            self.measurement_thread.start()

        return status, pos

    def set_voltage(self):

        """
        Sets Kethley voltage to value from GUI and writes current and voltage to GUI

        :return: No return
        """

        target_voltage = get_float(self.ui.target_V.text())

        status, (self.curr_V, self.curr_I) = self.com.request(self.smu_name, "apply_voltage", target_voltage)

        self.ui.curr_V.setText(str(self.curr_V))

        display_current = np.round((self.curr_I * 1e6),4)

        self.ui.curr_I.setText(str(display_current))

        if not status:
            return

        return

    # IV data handling

    def update_iv_data(self, data: tuple):
        """
        appends data given into self.i_v_list, data is in the format voltage,current

        :param data: Tuple with data to be stored in self.i_v_list in the format voltage,current
        :return: No return
        """

        # print(data)
        self.i_v_list[0].append(data[0])
        self.i_v_list[1].append(data[1])

        if len(self.i_v_list) > get_float(self.ui.volt_frames.text()):
            self.i_v_list_neg[0].append(data[0])
            self.i_v_list_neg[1].append(data[1])

        return

    def update_iv_graph(self):

        """
        Updates graph with final IV data, I current is multiplied by 1e6 to plot in uA
        Exports graph to png file with same name as the spe file
        Exports IV data to txt file with same name as spe file

        Option for semilog plot via checkbox, threshold is given in GUI

        :return: No return
        """

        self.path = str(self.ui.base_directory.text()).strip()

        i_v_list = self.i_v_list
        i_v_list_neg = self.i_v_list_neg

        multiplied_ydata = [ydata * 1e6 for ydata in i_v_list[1]]

        try:
            self.iv_data.set_data(i_v_list[0], multiplied_ydata)
        except Exception as e:
            self.show_message("Could not create IV data from dataset", e)

        if len(self.i_v_list_neg[0]) > 1:

            multiplied_ydata_neg = [ydata * 1e6 for ydata in i_v_list_neg[1]]

            try:
                self.iv_data_neg.set_data(i_v_list_neg[0], multiplied_ydata_neg)

            except Exception as e:
                self.show_message("Could not create IV data from dataset, 'negative' side", e)
        else:
            self.iv_data_neg.set_visible(False)

        try:
            self.ui.iv_curve.canvas.ax.set_xlabel('Voltage (V)')
            if self.ui.semilog_iv_check.isChecked():
                threshold = get_float(self.ui.threshold_current.text())
                self.ui.iv_curve.canvas.ax.set_yscale('symlog', linthresh=threshold)
                self.ui.iv_curve.canvas.ax.set_ylabel('Current (uA) (semilog)')
            else:
                self.ui.iv_curve.canvas.ax.set_yscale("linear")
                self.ui.iv_curve.canvas.ax.set_ylabel('Current (uA)')

            self.ui.iv_curve.canvas.ax.set_xlim([min(i_v_list[0]),
                                                 max(i_v_list[0])])
            self.ui.iv_curve.canvas.ax.set_ylim([min(multiplied_ydata) * 1.25, max(multiplied_ydata) * 1.25])

            if len(self.i_v_list_neg[0]) > 1:
                self.ui.iv_curve.canvas.ax.set_xlim([min(min(i_v_list[0]), min(i_v_list_neg[0])),
                                                 max(max(i_v_list[0]), max(i_v_list_neg[0]))])
                self.ui.iv_curve.canvas.ax.set_ylim([min(min(multiplied_ydata_neg), min(multiplied_ydata)),
                                                    max(max(multiplied_ydata_neg), max(multiplied_ydata))])

            self.ui.iv_curve.canvas.ax.legend(loc=2)
            self.ui.iv_curve.canvas.fig.tight_layout()
            self.ui.iv_curve.canvas.draw()
        except Exception as e:
            self.show_message("Could not create IV graph", e)

        base_name = self.ui.label_filename.text().replace(".spe", "")

        filename = str(base_name) + '_IV_trace_' + str(i_v_list[0][0]) + str(i_v_list[0][-1]) + '.png'

        today = date.today()

        year = today.strftime("%Y")
        month = today.strftime("%m")
        day = today.strftime("%d")

        path = os.path.join(self.path, year, month, str(day), str(self.ui.sample_name.text()).strip())

        try:
            pix = QtWidgets.QWidget.grab(self.ui.iv_curve)
            pix.save(os.path.join(path, filename), "PNG")

        except Exception as e:
            self.show_message("Could not create png file of IV curve.", e)
            self.ui.iv_curve.canvas.ax.set_title('NOT SAVED')
            self.ui.iv_curve.canvas.fig.tight_layout()
            self.ui.iv_curve.canvas.draw()

        filename = os.path.join(path, filename.replace(".png", ".txt"))
        # if os.path.isfile(filename): pass
        # else: open(filename,'x')

        try:
            with open(filename, "w+") as f:
                f.write("Current (uA);Voltage (V)\n")
                for line in np.arange(len(self.i_v_list[0])):
                    f.write("{};{}\n".format(i_v_list[1][line], i_v_list[0][line]))
                if len(self.i_v_list_neg[0]) > 1:
                    for line in np.arange(len(self.i_v_list_neg[0])):
                        f.write("{};{}\n".format(i_v_list[1][line], i_v_list[0][line]))
        except Exception as e:
            self.show_message("Could not create txt file of IV curve", e)

        return

    # =============================================================================
    # Turn digital devices ON/OFF
    # =============================================================================

    def white_led_on_off(self):

        """
        Sets white LED ON/OFF, defined by current state of LED

        :return: No return
        """

        if self.white_led_status:
            status, _ = self.com.request("michelson_linear", "set_whitelight", False)
            self.white_led_status = False
        else:
            status, _ = self.com.request("michelson_linear", "set_whitelight", True)
            self.white_led_status = True

            if not status:
                return

        return

    def blue_led_on_off(self):

        """
        Sets blue LED ON/OFF, defined by current state of LED

        :return: No return
        """

        if self.blue_led_status:
            status, _ = self.com.request("michelson_linear", "set_blue_led", False)
            self.blue_led_status = False
        else:
            status, _ = self.com.request("michelson_linear", "set_blue_led", True)
            self.blue_led_status = True

            if not status:
                return

        return

    def red_led_on_off(self):

        """
        Sets red LED ON/OFF, defined by current state of LED

        :return: No return
        """

        if self.red_led_status:
            status, _ = self.com.request("michelson_linear", "set_red_led", False)
            self.red_led_status = False
        else:
            status, _ = self.com.request("michelson_linear", "set_red_led", True)
            self.red_led_status = True

            if not status:
                return

        return

    def flip_mirror(self):  # flip mirror button function, currently used for testing

        try:
            status, _ = self.com.request("michelson_linear", "set_blue_led", True)
            status, _ = self.com.request("michelson_linear", "set_blue_led", False)

            if not status:
                raise Exception("Status while flipping Mirror!")
        except Exception as e:
            self.show_message("Error in flipping Mirror.",e)

        return

    # Elliptec devices control

    def piezos_positions(self, target: str = None, forced=False):

        """
        Define position for elliptec piezos depending on the type of measurement to be done.

        :param target: pl, imaging,michelson
        :type target: str
        :return: No return
        """

        if target is None:
            return

        if target == "pl":
            try:
                if forced:
                    try:
                        self.com.request(self.thr_bus_controller, "set_pos",
                                         {self.michelson_piezo_name: get_float(self.ui.michelson_mirror_out.text())})
                    except: pass

                    time.sleep(3)

                    try:
                        self.com.request(self.thr_bus_controller, "set_pos",
                                         {self.camera_piezo_name: get_float(self.ui.camera_mirror_out.text())})
                    except: pass
                    time.sleep(3)

                if "michelson" in self.curr_path:
                    self.flip_mirror()
                    time.sleep(0.1)

                    status,_ = self.com.request(self.thr_bus_controller, "get_pos",
                                    {self.michelson_piezo_name:None})

                    if (abs(_ - get_float(self.ui.michelson_mirror_out.text())) < 1e-5):
                        pass
                    else:

                        self.com.request(self.thr_bus_controller, "set_pos",
                                        {self.michelson_piezo_name: get_float(self.ui.michelson_mirror_out.text())})
                        time.sleep(0.5)

                if "imaging" in self.curr_path:

                    status,_ =  self.com.request(self.thr_bus_controller, "get_pos",
                                    {self.camera_piezo_name:None})

                    if abs(_ -get_float(self.ui.camera_mirror_out.text())) < 1e-3:
                        pass
                    else:

                        self.com.request(self.thr_bus_controller, "set_pos",
                                        {self.camera_piezo_name: get_float(self.ui.camera_mirror_out.text())})
                        time.sleep(0.5)

                    time.sleep(0.1)

                self.curr_path = "spec"
                print(self.curr_path)

            except Exception as e:
                self.show_message("Could not move Piezos", e)

            return

        if target == "imaging":

            if "imaging" in self.curr_path:
                try:
                    status, _ = self.com.request(self.thr_bus_controller, "get_pos",
                                                 {self.camera_piezo_name: None})
                    time.sleep(1)

                    if abs(_ - get_float(self.ui.camera_mirror_out.text())) < 1e-3:
                        pass
                    else:

                        self.com.request(self.thr_bus_controller, "set_pos",
                                         {self.camera_piezo_name: get_float(self.ui.camera_mirror_out.text())})
                        time.sleep(1)

                    self.curr_path = self.curr_path.replace(" imaging","")
                    print(self.curr_path)

                except Exception as e:
                    self.show_message("Could not move Piezos", e)
            else:
                try:
                    status,_ =  self.com.request(self.thr_bus_controller, "get_pos",
                                    {self.camera_piezo_name:None})
                    time.sleep(1)

                    if abs(_ -get_float(self.ui.camera_mirror_in.text())) < 1e-4:
                        pass
                    else:

                        self.com.request(self.thr_bus_controller, "set_pos",
                                        {self.camera_piezo_name: get_float(self.ui.camera_mirror_in.text())})
                        time.sleep(1)

                    self.curr_path += " imaging"
                    print(self.curr_path)

                except Exception as e:
                    self.show_message("Could not move Piezos", e)
                return

        if target == "michelson":
            try:
                if "spec" in self.curr_path :
                    self.flip_mirror()
                    time.sleep(1)

                    status,_ =  self.com.request(self.thr_bus_controller, "get_pos",
                                    {self.michelson_piezo_name:None})

                    print(status,_)

                    if abs(_ -get_float(self.ui.michelson_mirror_in.text())) < 1e-5:
                        pass
                    else:

                        self.com.request(self.thr_bus_controller, "set_pos",
                                        {self.michelson_piezo_name: get_float(self.ui.michelson_mirror_in.text())})
                        time.sleep(1)

                if "imaging" in self.curr_path:

                    status,_ =  self.com.request(self.thr_bus_controller, "get_pos",
                                    {self.camera_piezo_name:None})

                    if abs(_ -get_float(self.ui.camera_mirror_out.text())) < 1e-5:
                        pass
                    else:

                        self.com.request(self.thr_bus_controller, "set_pos",
                                        {self.camera_piezo_name: get_float(self.ui.camera_mirror_out.text())})
                        time.sleep(3)

                self.curr_path = "michelson"
                print(self.curr_path)

            except Exception as e:
                self.show_message("Could not move Piezos", e)

            return
            

    def update_piezo_pos(self, piezo: str, state: str):
        """

        Update the in/ou positions of the elliptec piezos within the GUI.
        Request position from DCS for the piezo and "print" it into the GUI for future use.

        :param piezo: name of the piezo to update on the GUI
        :type piezo: str
        :param state: in or out, if the piezo is in the optical path or out of the optical path
        :type state: str
        :return: No return
        """
        pos = None

        if piezo == "periscope":

            if not self.periscope_piezo_state:
                return

            try:

                status, pos = self.com.request(self.thr_bus_controller, "get_pos", {self.periscope_piezo_name: None})

                if status:

                    if state == 'in':
                        self.ui.periscope_mirror_in.setText(str(pos))
                    elif state == 'out':
                        self.ui.periscope_mirror_out.setText(str(pos))

                return
            except Exception as e:
                self.show_message("General Error while updating periscope elliptec piezo pos {}".format(pos), e)
                self.periscope_piezo_state = False
                return

        elif piezo == "michelson":

            if not self.michelson_piezo_state:
                return

            try:

                status, pos = self.com.request(self.thr_bus_controller, "get_pos", {self.michelson_piezo_name: None})

                if status:

                    if state == 'in':
                        self.ui.michelson_mirror_in.setText(str(pos))
                    elif state == 'out':
                        self.ui.michelson_mirror_out.setText(str(pos))

                return
            except Exception as e:
                self.show_message("General Error while updating michelson elliptec piezo pos {}".format(pos), e)
                self.michelson_piezo_state = False
                return

        elif piezo == "camera":

            if not self.camera_piezo_state:
                return

            try:

                status, pos = self.com.request(self.thr_bus_controller, "get_pos", {self.camera_piezo_name: None})

                if status:

                    if state == 'in':
                        self.ui.camera_mirror_in.setText(str(pos))
                    elif state == 'out':
                        self.ui.camera_mirror_out.setText(str(pos))

                return
            except Exception as e:
                self.show_message("General Error while updating camera elliptec piezo pos {}".format(pos), e)
                self.camera_piezo_state = False
                return

    def jog_piezo(self, piezo: str, direction: str):

        """
        Set elliptec piezo jog distance from GUI value for the equivalent piezo and move it in the desired direction.
        Step size is given in um

        :param piezo: Piezo name: periscope, michelson or camera.
        :type piezo: str
        :param direction: forward or backward
        :type direction: str
        :return: No return
        """

        piezo_name = piezo
        jog = 1
        pos = None

        if piezo == "periscope":
            jog = get_float(self.ui.periscope_mirror_jog.text())
            piezo_name = self.periscope_piezo_name

            if not self.periscope_piezo_state:
                return

        elif piezo == "michelson":

            jog = get_float(self.ui.michelson_mirror_jog.text())
            piezo_name = self.michelson_piezo_name

            if not self.michelson_piezo_state:
                return

        elif piezo == "camera":

            jog = get_float(
                self.ui.camera_mirror_jog.text())  # float(self.ui.camera_mirror_jog.text().replace(",","."))
            piezo_name = self.camera_piezo_name

            if not self.camera_piezo_state:
                return

        try:

            status, pos = self.com.request(self.thr_bus_controller, "set_jog", {piezo_name: jog})
            time.sleep(.1)
            status, pos = self.com.request(self.thr_bus_controller, "jog", {piezo_name: direction})

            if not status:
                return

        except Exception as e:
            self.show_message("General Error while jogging {} elliptec piezo pos {}".format(piezo_name, pos), e)

        return

    # =============================================================================
    # Set frames/steps
    # =============================================================================

    def volt_frames(self):

        """

        Check number of frames needed for the desired kethley measurement based on voltages (min, max) in GUI and voltage step
        and add write it in the GUI

        :return: No return
        """

        v_max = get_float(self.ui.volt_max.text())
        v_min = get_float(self.ui.volt_min.text())
        v_ste = get_float(self.ui.volt_step.text())

        try:

            frames = (v_max - v_min) / v_ste

            frames = str(int(frames) + 1)

            self.ui.volt_frames.setText(frames)
        except Exception as e:
            self.show_message("Error in setting frames for voltage.", e)

        return

    def update_fss_step(self):

        """

        Check the step size of polarization measurement based on initial and end angle and number of frames in the GUI and write result to GUI

        :return: No return
        """

        frames = get_float(self.ui.fss_frames.text()) - 1
        starting_angle = get_float(self.ui.fss_i.text())
        end_angle = get_float(self.ui.fss_end.text())

        try:
            step = (end_angle - starting_angle) / frames

            step = round(step, 3)
            self.ui.fss_step.setText(str(step))
        except Exception as e:
            self.show_message("Error in setting frames for Polarization.", e)

        return

    def laser_step(self):

        """

        Check the step size of power dependant measurement based on initial and end angle and number of frames in the GUI and write result to GUI

        :return: No return
        """

        frames = get_float(self.ui.pow_map_frames.text()) - 1
        starting_angle = get_float(self.ui.power_map_i.text())
        end_angle = get_float(self.ui.power_map_f.text())

        cw_frames = get_float(self.ui.pow_map_frames.text()) - 1
        cw_starting_angle = get_float(self.ui.cw_power_map_i.text())
        cw_end_angle = get_float(self.ui.cw_power_map_f.text())

        try:

            step = (end_angle - starting_angle) / frames

            step = round(step, 3)

            self.ui.step_laser.setText(str(step))
        except Exception as e:
            self.show_message("Error in setting laser power frames.", e)

        try:

            step = (cw_end_angle - cw_starting_angle) / cw_frames

            step = round(step, 3)

            self.ui.cw_laser_step.setText(str(step))
        except Exception as e:
            self.show_message("Error in setting cw laser power frames.", e)

    def ple_frames(self):

        ple_frames = get_float(self.ui.ple_frames.text()) - 1
        ple_starting_wl = get_float(self.ui.ple_start.text())
        ple_end_wl = get_float(self.ui.ple_end.text())

        try:

            step = (ple_end_wl - ple_starting_wl) / ple_frames

            step = round(step, 3)

            self.ui.ple_step.setText(str(step))
        except Exception as e:
            self.show_message("Error in setting cw laser power frames.", e)

    def meas_2d_frames(self):
        """
        Calculates number of frames for 2d measurement

        :return: No return
        """
        try:
            start_slow = get_float(self.ui.slow_axis_start.text())
            end_slow = get_float(self.ui.slow_axis_end.text())
            step_slow = get_float(self.ui.slow_axis_steps.text())
            slow_frames = np.round((np.abs(end_slow - start_slow) / step_slow))
            slow_frames = int(np.round((np.abs(end_slow - start_slow) / step_slow), 0)) + 1

            self.ui.slow_axis_frames.setText(str(slow_frames))

            start_fast = get_float(self.ui.fast_axis_start.text())
            end_fast = get_float(self.ui.fast_axis_end.text())
            step_fast = get_float(self.ui.fast_axis_steps.text())
            fast_frames = int(np.round((np.abs(end_fast - start_fast) / step_fast), 0)) + 1

            self.ui.fast_axis_frames.setText(str(fast_frames))


        except Exception as e:
            self.show_message("Error in setting 2d frames.", e)
            print("Error in setting 2d frames.", e)

        return

    def space_map_frames(self):

        """
        Calculates number of frames needed for space mapping based on number of rows and columns from the GUI and write the result in the GUI


        :return: No return
        """

        ncol = get_float(self.ui.space_map_colums.text())
        nlines = get_float(self.ui.space_map_lines.text())

        try:
            frames = ncol * nlines

            frames = str(int(frames))

            self.ui.space_map_frames.setText(frames)
        except Exception as e:
            self.show_message("Error in setting space map frames.", e)

        return

    def change_frames(self, frames: int | bool):

        """
        Change the number of frames that te spectrometer should prepare to acquire

        :param frames: Number of desired frames
        :type frames: int or bool
        :return: No return
        """

        if self.spectrometer_connected:
            if self.running_measurement is None:
                if type(frames) is bool and frames is True:
                    self.com.request(self.spectrometer_name, "stop")
                    time.sleep(0.1)
                    try:
                        self.com.request(self.spectrometer_name, "frames_to_save", "1")
                    except Exception as e:
                        self.show_message("Error in setting number of frames into spectrometer.", e)
                        return

                else:
                    frames = int(np.round(frames, 0))
                    try:
                        self.com.request(self.spectrometer_name, "frames_to_save", frames)
                    except Exception as e:
                        self.show_message("Error in setting number of frames into spectrometer.", e)

    # Set filename label from experiment settings

    def set_label_filename(self):

        """
        Set the appropriate filename for the measurement into the GUI label based on the GUI info and selected measurement type.
        The filename is set according to the requirements given in the latest XRSP documentation.

        :return: No return
        """

        sample_name = str(self.ui.sample_name.text()).strip()
        sample_name2 = str(self.ui.sample_name_2.text()).strip()
        sample_name3 = str(self.ui.sample_name_3.text()).strip()
        temperature = str(self.ui.temperature.text()).strip()
        filter_name = str(self.ui.filter_name.text()).strip()
        slit = str(self.ui.slit.text()).strip()
        exposure = str(self.ui.exposure.text()).strip()
        grating = str(self.ui.grating.currentText()).strip()  # combo_box
        center_wl = str(self.ui.center_wl.text()).strip()
        start_wl = str(self.ui.start_wl.text()).strip()
        end_wl = str(self.ui.end_wl.text()).strip()
        cw_wl = str(self.ui.cw_wl.text()).strip()
        sufix_name = ""

        try:
            status, info = self.com.request(self.spectrometer_name, "get_acq_info")
            
            if info[0] == None:
                pass
            else:
                exposure = info[0]

            if info[1] == None:
                pass
            else:
                center_wl = info[1]

            if info[2] == None:
                pass
            else:
                grating = info[2]

            if info[3] == None:
                pass
            else:
                slit = info[3]

            try:
                self.ui.slit.setText(str(slit))
                self.ui.exposure.setText(str(exposure))
                # self.ui.grating.setText(str(grating)) #TODO: How to do this properly
                self.ui.center_wl.setText(str(center_wl))

            except Exception as e:
                self.show_message("Could not set spectrometer settings to GUI.", e)

        except Exception as e:
            self.show_message("Error in getting acquisition info from spectrometer!", e)

        source = find_radio_button(self.ui.sources_group.buttons())
        source.strip().replace(' ', '').replace(" ", "")

        if slit in ("manual", "side"): slit+="_"
        else: slit+="µmSlit_"

        if source == "Green":
            source = "GreenLaser"
        elif source == 'CW':
            source = 'CW_' + str(cw_wl)

        # _QD1_50umslit_1800grating_greenlaser_FEL600_1s_7K#Pol(deg)_0_360

        standard_name = sample_name

        exposure_units = 'ms'

        try:
            if int(exposure) < 10:
                exposure_units = 's'
        except:
            pass

        if sample_name2.strip() not in (None, ''):
            standard_name += '_' + sample_name2 + '_' + slit + grating + 'grating_' + source + \
                             '_' + filter_name + '_' + exposure + exposure_units + '_'
        else:
            standard_name += '_' + slit + grating + 'grating_' + source + \
                             '_' + filter_name + '_' + exposure + 'ms_' +"center_wl_"+ center_wl + 'nm_' 

        if temperature.lower() == 'rt':
            standard_name += 'RT'
        else:
            standard_name += temperature + 'K'

        # standard_pl,step_and_glue,space_map,fss_measure,power_map,michelson_int

        measurement = find_radio_button(self.ui.measurements_group.buttons())

        if measurement.lower() == 'polarization': measurement = "FSS"

        if measurement == 'Standard':
            sufix_name = ''
        elif measurement == 'Step and Glue':
            sufix_name = start_wl + '_' + end_wl

        elif measurement == 'Line Scan':

            axis = str(find_radio_button(self.ui.line_sca_group.buttons()))

            step = str(self.ui.line_scan_step.text())
            spectras = str(self.ui.line_scan_points.text())

            sufix_name = center_wl + 'nm_' + '#line_scan_'+axis+'#(µm)_' + spectras + '_' + step

        elif measurement == 'Space map':

            ncol = str(self.ui.space_map_colums.text())
            nlines = str(self.ui.space_map_lines.text())
            step_x = str(self.ui.column_step.text())
            step_y = str(self.ui.line_step.text())

            if self.ui.snake_like.isChecked():
                sufix_name+="_SnakeLike_#spacemap#(µm)_-" + ncol + '_' + step_x

            else:
                sufix_name += "#spacemap#(µm)_" + ncol + '_' + step_x

            if nlines != '0' and step_y != step_x:
                sufix_name += '_' + step_y

        elif measurement == 'FSS':
            starting_angle = str(2 * get_float(self.ui.fss_i.text()))
            end_angle = str(2 * get_float(self.ui.fss_end.text()))

            sufix_name = '#Pol(deg)_' + starting_angle + '_' + end_angle

        elif measurement == 'Power Map':

            starting_angle = str(self.ui.power_map_i.text())
            end_angle = str(self.ui.power_map_f.text())
            min_power = str(self.ui.minpower.text())
            max_power = str(self.ui.maxpower.text())

            sufix_name = '#power_map#(deg)_' + starting_angle + '_' + \
                         end_angle #+ '_' + min_power + '_' + max_power

        elif measurement == 'Voltage Sweep':

            volt_min = str(self.ui.volt_min.text())
            volt_max = str(self.ui.volt_max.text())

            sufix_name = '#Voltage(V)_' + volt_min + '_' + volt_max

        elif measurement == 'Michelson':

            piezo_steps = get_float(self.ui.piezo_step.value()) + 1
            piezo_step_size = get_float(self.ui.piezo_stepsize.text())
            linear_step_size = int(get_float(self.ui.linear_stepsize.text()) * 1000)

            sufix_name = '#michelson#' + str(piezo_steps) + '_' + str(piezo_step_size) + '_' + str(linear_step_size)

        elif measurement == '2D Measurement':
            slow_axis = find_radio_button(self.ui.slow_axis_radio.buttons())
            fast_axis = find_radio_button(self.ui.fast_axis_radio.buttons())

            slow_axis_unit = fast_axis_unit = ""

            slow_axis_start = get_float(self.ui.slow_axis_start.text())
            slow_axis_end = get_float(self.ui.slow_axis_end.text())

            fast_axis_start = get_float(self.ui.fast_axis_start.text())
            fast_axis_end = get_float(self.ui.fast_axis_end.text())

            if slow_axis.lower() in ('polarization', 'power map'):
                slow_axis_unit = '(deg)'
                slow_axis_start *= 2
                slow_axis_end *= 2
            elif slow_axis.lower() == 'focus':
                slow_axis_unit = '(um)'
            elif slow_axis.lower() == 'volage':
                slow_axis_unit = '(V)'

            if fast_axis.lower() in ('polarization', 'power map'):
                fast_axis_unit = '(deg)'
                fast_axis_start *= 2
                fast_axis_end *= 2
            elif fast_axis.lower() == 'focus':
                fast_axis_unit = '(um)'
            elif fast_axis.lower() == 'volage':
                fast_axis_unit = '(V)'

            sufix_name = '#' + slow_axis + slow_axis_unit + str(slow_axis_start) + '_' + str(slow_axis_end) + \
                         '#' + fast_axis + fast_axis_unit + str(fast_axis_start) + '_' + str(fast_axis_end)

        if sufix_name != '':
            sufix_name = '_' + sufix_name

        name_noextension = standard_name + sufix_name

        if sample_name3.strip() not in (None, ''):
            name_noextension += sample_name3

        complete_name = name_noextension + '.spe'

        self.ui.label_filename.setText(complete_name)

    # Spectrometer functions

    def update_filename_spectrometer(self):

        """
        Sends the filename from the GUI label to lightfield for the next acquired spectra
        Name convention follows the XRSP documentation.

        :return: No return
        """

        self.set_label_filename()
        time.sleep(0.2)

        sample_name = str(self.ui.sample_name.text()).strip()
        filename = self.ui.label_filename.text().replace('.spe', "").strip()

        today = date.today()

        year = today.strftime("%Y")
        month = today.strftime("%m")
        day = today.strftime("%d")

        target_dir = os.path.join(self.path, year, month, str(day))

        if sample_name.strip() not in (None, ''):
            target_dir = os.path.join(target_dir, sample_name.strip())

        self.path = str(self.ui.base_directory.text()).strip()

        try:
            Path(target_dir).mkdir(parents=True, exist_ok=True)
        except Exception as e:
            self.show_message("Error in creating folder for experiment/sample", e)

        if self.running_measurement is None:

            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)

                try:
                    self.com.request(self.spectrometer_name, "change_folder", str(target_dir))
                except Exception as e:
                    self.show_message("Error in setting folder in spectrometer.", e)

                time.sleep(0.1)

                self.com.request(self.spectrometer_name, "change_filename",
                                 filename)

    def toggle_trigger(self, state: bool):

        """
        Toggle trigger state of spectrometer between no response and readout per trigger
        True -> readout per trigger
        False -> No response

        :param state: Readout per tigger (triggering ON, state=True) or No response (triggering OFF, state=False)
        :type: bool
        :return: No return
        """

        trigger = "noresponse"

        if self.spectrometer_connected:
            if self.running_measurement is None:
                if state is True:
                    trigger = "readout"
                elif state is False:
                    trigger = "noresponse"

                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)

                try:
                    self.com.request(self.spectrometer_name, "trigger_settings", trigger)
                except Exception as e:
                    self.show_message("Error in setting spectrometer trigger.", e)
                    return

    # =============================================================================
    # State UI toggling
    # =============================================================================

    def toggle_all(self, state: bool):

        """
        Toggle all elements of GUI to the specified state.

        :param state: Final UI state
        :type state: bool
        :return: No return
        """

        if state:
            meas_axis = self.running_measurement
        if not state:
            meas_axis = 'None'

        self.ui.measurement_button.setEnabled(state)
        self.ui.abort_measurent.setEnabled(not state)
        self.ui.linear_loop.setEnabled(state)
        self.ui.piezo_loop.setEnabled(state)
        self.toggle_ui(state)

        if state:
            filename = self.ui.label_filename.text() + '.spe'

            today = date.today()

            year = today.strftime("%Y")
            month = today.strftime("%m")
            day = today.strftime("%d")

            self.path = str(self.ui.base_directory.text()).strip()
            sample_name = str(self.ui.sample_name.text()).strip()

            target_dir = os.path.join(self.path, year, month, str(day), sample_name.strip())

            filename = os.path.join(target_dir, filename)

            if os.path.isfile(filename) and self.ui.treat_data_checkbox.isChecked():
                try:
                    self.treat_last_data(filename)
                except Exception as e:
                    self.show_message("Error in data treatment!", e)

    def toggle_measurement(self, state: bool):

        """
        Toggle GUI elements to state, and afterwards toggle the button of the selected measurement to not state
         in order to stop the measurement if needed

        :param state: state of the GUI
        :type state: bool
        :return: No return
        """

        self.ui.linear_loop.setEnabled(state)
        self.ui.piezo_loop.setEnabled(state)
        self.toggle_ui(state)

        self.ui.space_map_group.setEnabled(state)
        self.ui.line_scan_group.setEnabled(state)
        self.ui.common_laser_tab.setEnabled(state)
        self.ui.cw_group_box.setEnabled(state)
        self.ui.cw_powermap_box.setEnabled(state)
        self.ui.ple_box.setEnabled(state)
        self.ui.kethley_box.setEnabled(state)
        self.ui.polarization_control.setEnabled(state)

        if not state:

            if self.running_measurement == 'space_map':
                self.ui.space_map_button.setEnabled(not state)
            if self.running_measurement == 'line_scan':
                self.ui.line_scan_button.setEnabled(not state)
            if self.running_measurement == 'power_map':
                self.ui.power_map_button.setEnabled(not state)
            if self.running_measurement == 'power_map_3900s':
                self.ui.cw_power_map_button.setEnabled(not state)
            if self.running_measurement == 'kethley':
                self.ui.volt_sweep_button.setEnabled(not state)
            if self.running_measurement == 'fss':
                self.ui.fss_button.setEnabled(not state)
            if self.running_measurement == 'ple':
                self.ui.ple_button.setEnabled(not state)
            if self.running_measurement == '2d_meas':
                self.ui.measure_2d_button.setEnabled(not state)

    def toggle_linear_loop(self, state: bool):

        """
        Toggle elements of the michelson linear stage to state during loop.

        :param state: Inverse of state of loop
        :type state: bool
        :return: No return
        """

        self.ui.measurement_button.setEnabled(state)
        self.toggle_ui(state)
        self.ui.linear_loop.setEnabled(not state)

    def toggle_piezo_loop(self, state: bool):

        """
        Toggle elements of the michelson piezo to state during loop.

        :param state: Inverse of state of loop
        :type state: bool
        :return: No return
        """

        self.ui.measurement_button.setEnabled(state)
        self.toggle_ui(state)
        self.ui.piezo_loop.setEnabled(not state)

    def toggle_abort_button(self, completed: bool):

        """
        Set the abort button in the GUI to the inverse of the state given

        :param completed: Inverse of the state to be set, it translates to if the measurement was completed or not
        :return: No return
        """

        self.ui.abort_measurent.setEnabled(not completed)

    def toggle_ui(self, state: bool):

        """
        Set the GUI elements to stat
        If there was a previous ongoing measurement request the spectrometer to stop any ongoing acquisition
        return it to the initial state and updates current positions related to the previous measurement


        :param state: State the set GUI
        :return: No return
        """

        # buttons
        self.ui.go_button.setEnabled(state)
        self.ui.go_button_2.setEnabled(state)
        self.ui.save_config.setEnabled(state)
        self.ui.load_config.setEnabled(state)
        self.ui.retry_connection.setEnabled(state)
        self.ui.check_spacemap_points.setEnabled(state)
        # inputs
        self.ui.linear_start.setEnabled(state)
        self.ui.linear_end.setEnabled(state)
        self.ui.linear_step.setEnabled(state)
        self.ui.linear_goto.setEnabled(state)
        self.ui.linear_spin.setEnabled(state)
        self.ui.piezo_start.setEnabled(state)
        self.ui.piezo_end.setEnabled(state)
        self.ui.piezo_step.setEnabled(state)
        self.ui.piezo_goto.setEnabled(state)
        self.ui.piezo_spin.setEnabled(state)
        self.ui.vel_spin.setEnabled(state)
        self.ui.acc_spin.setEnabled(state)
        self.ui.dec_spin.setEnabled(state)
        self.ui.trigger_linear.setEnabled(state)
        self.ui.trigger_piezo.setEnabled(state)
        # standard tab
        self.ui.apply_button.setEnabled(state)
        self.ui.quality.setEnabled(state)
        self.ui.linewidth.setEnabled(state)
        self.ui.xy_movement.setEnabled(state)
        self.ui.space_map_group.setEnabled(state)
        self.ui.line_scan_group.setEnabled(state)
        self.ui.stored_positions_3.setEnabled(state)
        self.ui.z_movement.setEnabled(state)
        self.ui.fit_plane.setEnabled(state)
        # self.ui.fast_map.setEnabled(state)
        self.ui.snake_like.setEnabled(state)
        self.ui.power_mapping.setEnabled(state)
        self.ui.experiment_settings_box.setEnabled(state)
        self.ui.kethley_box.setEnabled(state)
        self.ui.devices_box.setEnabled(state)
        self.ui.polarization_control.setEnabled(state)
        self.ui.cw_group_box.setEnabled(state)
        self.ui.ple_tab.setEnabled(state)
        self.ui.cw_meas_tabs.setEnabled(state)
        self.ui.cw_powemap_tab.setEnabled(state)

        self.ui.stored_positions_2.setEnabled(state)
        self.ui.load_pos.setEnabled(state)
        self.ui.add_position.setEnabled(state)

        if state and self.running_measurement is not None:

            self.com.request(self.spectrometer_name, "stop")
            time.sleep(0.1)

            self.com.request(self.xy_stage_name, 'set_acc', 25)

            if self.running_measurement == 'space_map':
                self.update_current_positions('xyz')
            elif self.running_measurement == 'line_scan':
                self.update_current_positions('xyz')
            elif self.running_measurement == 'power_map':
                self.update_current_positions('laser')
            elif self.running_measurement == 'kethley':
                self.update_current_positions('kethley')
            elif self.running_measurement == 'fss':
                self.update_current_positions('fss')
            elif self.running_measurement == 'ple':
                self.update_current_positions('cw_laser')
            else:
                self.update_current_positions(None)
            self.running_measurement = None
            self.ui.standard_pl.setChecked(True)
            self.update_filename_spectrometer()
            time.sleep(0.1)

            self.change_frames(1)
            time.sleep(0.1)
            self.toggle_trigger(False)
            time.sleep(0.1)

    def toggle_linear_loop_button(self, completed: bool):
        """
        Changes the text of  Michelson linear loop button according to state

        :param completed: If True set text to Loop, else set it to End
        :type completed: bool
        :return: No return
        """

        if completed:
            self.ui.linear_loop.setText("Loop")
        else:
            self.ui.linear_loop.setText("End")

    def toggle_piezo_loop_button(self, completed):
        """
        Changes the text of  Michelson piezo loop button according to state

        :param completed: If True set text to Loop, else set it to End
        :type completed: bool
        :return: No return
        """

        if completed:
            self.ui.piezo_loop.setText("Loop")
        else:
            self.ui.piezo_loop.setText("End")

    # =============================================================================
    # Load/store settings in .json
    # =============================================================================

    def standard_settings(self):

        """
        Loads a set of standard setting from settings.json file
        Each of the defined settings is then written into the appropriate place in the GUI
        If error occurs it is shown and logged

        :return: No return
        """

        filename = "settings.json"

        cwd = os.getcwd()

        filename = os.path.join(cwd,"configs", filename)

        with open(filename, "r") as f:
            settings = json.load(f)
        # just for info, if this is no longer the case, one has to change that we still hit the equal_path_length
        try:
            self.michelson_equal_length = equal_path_length = settings["equal_path_length"]  # mm
        except Exception as e:
            self.show_message("Error in setting equal length path.", e)
            self.michelson_equal_length = equal_path_length = 45.5
        finally:
            self.ui.equal_path_length.setText(str(self.michelson_equal_length))

        linewidth = self.ui.linewidth.currentIndex()
        quality = self.ui.quality.currentIndex()

        # steps are total images - 1 as an image before changing for a step is taken
        steps = [5, 10, 30, 60][quality]
        width = [100, 50, 20][linewidth]

        min_linear_stage = float(settings["first_point"])
        max_linear_stage = float(settings["last_point"])

        if self.ui.quality.currentIndex() == 0:
            first_point = equal_path_length - min(width / 2, equal_path_length)
            last_point = equal_path_length + 2 * min(width / 2, equal_path_length)

        else:
            first_point = max(0, equal_path_length - width)
            last_point = equal_path_length + width

            stepsize = ((last_point - first_point) * 10) // steps / 10

            # print(stepsize, equal_path_length % stepsize)

            first_point = first_point + round(equal_path_length % stepsize, 2)
            last_point = first_point + steps * stepsize

        if first_point < min_linear_stage:
            first_point = min_linear_stage
        if last_point > max_linear_stage:
            last_point = max_linear_stage

        try:
            self.ui.base_directory.setText(settings["base_directory"])
        except Exception as e:
            self.show_message("Error in reading json settings.", e)
            self.ui.base_directory.setText("I:\\public\\NANOSCALE SEMICONDUCTOR GROUP\\1. DATA\\SMALL-LAB")
        try:
            self.ui.experiment_name.setText(settings["experiment_name"])
        except Exception as e:
            self.show_message("Error in reading json settings.", e)
            self.ui.experiment_name.setText("2022_07_standard_experiment")
        try:
            self.ui.periscope_mirror_in.setText(str(settings["periscope_mirror_in"]))
        except Exception as e:
            self.show_message("Error in reading json settings.", e)
            self.ui.periscope_mirror_in.setText("0.0")
        try:
            self.ui.periscope_mirror_out.setText(str(settings["periscope_mirror_out"]))
        except Exception as e:
            self.show_message("Error in reading json settings.", e)
            self.ui.periscope_mirror_out.setText("58.0")
        try:
            self.ui.michelson_mirror_in.setText(str(settings["michelson_mirror_in"]))
        except Exception as e:
            self.show_message("Error in reading json settings.", e)
            self.ui.michelson_mirror_in.setText("0.0")
        try:
            self.ui.michelson_mirror_out.setText(str(settings["michelson_mirror_out"]))
        except Exception as e:
            self.show_message("Error in reading json settings.", e)
            self.ui.michelson_mirror_out.setText("60.0")
        try:
            self.ui.camera_mirror_in.setText(str(settings["camera_mirror_in"]))
        except Exception as e:
            self.show_message("Error in reading json settings.", e)
            self.ui.camera_mirror_in.setText("0.0")
        try:
            self.ui.camera_mirror_out.setText(str(settings["camera_mirror_out"]))
        except Exception as e:
            self.show_message("Error in reading json settings.", e)
            self.ui.camera_mirror_out.setText("60.0")
        try:
            self.ui.linear_start.setText(str(first_point))
        except Exception as e:
            self.show_message("Error in reading json settings.", e)
            self.ui.linear_start.setText("3.30")
        try:
            self.ui.linear_end.setText(str(last_point))
        except Exception as e:
            self.show_message("Error in reading json settings.", e)
            self.ui.linear_end.setText("180.0")
        try:
            self.ui.pol_offset.setText(str(settings["pol_offset"]))
        except Exception as e:
            self.show_message("Error in reading json settings.", e)
            self.ui.pol_offset.setText("25.0")
        try:
            self.ui.laser_offset.setText(str(settings["laser_offset"]))
        except Exception as e:
            self.show_message("Error in reading json settings.", e)
            self.ui.laser_offset.setText("69.49")

        # set to one and to correct value to trigger change event
        self.ui.linear_step.setValue(1)
        self.ui.linear_step.setValue(steps)

    def update_standard_settings(self):
        standard_settings = {"equal_path_length"   : get_float(self.ui.equal_path_length.text()),
                             "base_directory"      : self.ui.base_directory.text(),
                             "experiment_name"     : self.ui.experiment_name.text(),
                             "periscope_mirror_in" : self.ui.periscope_mirror_in.text(),
                             "periscope_mirror_out": self.ui.periscope_mirror_out.text(),
                             "michelson_mirror_in" : self.ui.michelson_mirror_in.text(),
                             "michelson_mirror_out": self.ui.michelson_mirror_out.text(),
                             "camera_mirror_in"    : self.ui.camera_mirror_in.text(),
                             "camera_mirror_out"   : self.ui.camera_mirror_out.text(),
                             "first_point"         : self.ui.linear_start.text(),
                             "last_point"          : self.ui.linear_end.text(),
                             "pol_offset"          : get_float(self.ui.pol_offset.text()),
                             "laser_offset"        : get_float(self.ui.laser_offset.text()),
                             "cw_laser_offset"     : get_float(self.ui.cw_laser_offset.text())
                             }


        filename = "settings_PL.json"

        cwd = os.getcwd()

        filename = os.path.join(cwd,"configs", filename)

        with open(filename, "w") as f:
            settings = f
            settings.write(json.dumps(standard_settings))
            settings.close()

        return

    def save_config(self):

        """
        Write current michelson settings to a json file

        :return: No return
        """

        if not self.check_text():
            return
        path = file_dialog('Choose Config File Location', False)
        if path is not None:
            with open(path, "w") as f:
                f.write(
                    "linear_start(mm);linear_end(mm);linear_step;linear_stepsize(mm);linear_goto(mm);")
                f.write(
                    "piezo_start(um);piezo_end(um);piezo_step;piezo_stepsize(um);piezo_goto(um);")
                f.write("linear_vel(mm/s);linear_acc(mm/s^2);linear_dec(mm/s^2)\n")
                f.write(self.ui.linear_start.text() + ";")
                f.write(self.ui.linear_end.text() + ";")
                f.write(f"{self.ui.linear_step.value()};")
                f.write(self.ui.linear_stepsize.text() + ";")
                f.write(self.ui.linear_goto.text() + ";")
                f.write(self.ui.piezo_start.text() + ";")
                f.write(self.ui.piezo_end.text() + ";")
                f.write(f"{self.ui.piezo_step.value()};")
                f.write(self.ui.piezo_stepsize.text() + ";")
                f.write(self.ui.piezo_goto.text() + ";")
                f.write(f"{self.ui.vel_spin.value()};")
                f.write(f"{self.ui.acc_spin.value()};")
                f.write(f"{self.ui.dec_spin.value()}")

    def load_config(self):
        """
        Load michelson settings from a json file, and updates velocity, acceleration and deceleration of michelson linear stage.

        :return: No return
        """

        path = file_dialog('Choose Config File Location', True)
        if path is not None:
            with open(path, "r") as f:
                for i, data in enumerate(f):
                    if i == 0:
                        continue
                    data = data.split(";")
                    self.ui.linear_start.setText(data[0])
                    self.ui.linear_end.setText(data[1])
                    self.ui.linear_step.setValue(int(data[2]))
                    self.ui.linear_stepsize.setText(data[3])
                    self.ui.linear_goto.setText(data[4])
                    self.ui.piezo_start.setText(data[5])
                    self.ui.piezo_end.setText(data[6])
                    self.ui.piezo_step.setValue(int(data[7]))
                    self.ui.piezo_stepsize.setText(data[8])
                    self.ui.piezo_goto.setText(data[9])
                    self.ui.vel_spin.setValue(int(data[10]))
                    self.ui.acc_spin.setValue(int(data[11]))
                    self.ui.dec_spin.setValue(int(data[12]))

        self.check_text()
        self.set_linear_velocity()
        self.set_linear_acceleration()
        self.set_linear_deceleration()

    # =============================================================================
    # Set devices configs for michelson piezo/linear
    # =============================================================================

    def set_linear_velocity(self):
        """
        Set velocity of michelson linear stage according to value in the GUI

        :return: No return
        """

        if self.check_text():
            try:
                self.com.request("michelson_linear", "vel", self.ui.vel_spin.value())
            except Exception as e:
                self.show_message("Error while setting liner stage velocity.", e)

    def set_linear_acceleration(self):
        """
        Set acceleration of michelson linear stage according to value in the GUI

        :return: No return
        """

        if self.check_text():
            try:
                self.com.request("michelson_linear", "acc", self.ui.vel_spin.value())
            except Exception as e:
                self.show_message("Error while setting liner stage acceleration.", e)

    def set_linear_deceleration(self):
        """
        Set deceleration of michelson linear stage according to value in the GUI

        :return: No return
        """

        if self.check_text():
            try:
                self.com.request("michelson_linear", "dec", self.ui.vel_spin.value())

            except Exception as e:
                self.show_message("Error in setting linear Stage deceleration.", e)

    def goto_linear(self):
        """
        Moves michelson linear stage to the position in the GUI, if the trigger option is selected also sends a trigger to spectrometer.

        :return: No return
        """

        if self.check_text():
            target_point = get_float(self.ui.linear_goto.text())
            try:
                status, self.curr_linear_pos = self.com.request("michelson_linear", "set_pos", target_point)
                self.ui.linear_cur_position.setText(str(self.curr_linear_pos))
            except Exception as e:
                self.show_message("Error in moving linear Stage", e)
            if self.ui.trigger_linear.isChecked():
                try:
                    self.com.request("michelson_linear", "vel", self.ui.vel_spin.value())
                    self.com.request("michelson_linear", "trig")
                except Exception as e:
                    self.show_message("Error in setting linear Stage trigger.", e)

    def goto_piezo(self):
        """
        Moves michelson piezo stage to the position in the GUI, if the trigger option is selected also sends
        a trigger to spectrometer.

        :return: No return
        """

        if self.check_text():
            target_point = get_float(self.ui.piezo_goto.text())
            try:
                status, self.curr_mich_piezo_pos = self.com.request("michelson_piezo", "set_pos", target_point)
                self.ui.piezo_cur_position.setText(str(self.curr_mich_piezo_pos))
            except Exception as e:
                self.show_message("Error in moving Michelson Piezo Stage", e)
            if self.ui.trigger_piezo.isChecked():
                try:
                    self.com.request("michelson_linear", "trig")
                except Exception as e:
                    self.show_message("Error in triggering with linear Stage", e)

    def update_cw_wl_label(self,motor_pos):
        """updates cw wl and motor pos from motor pos input
        :params motor_pos: float -> float of motor position
        :return:
        """
        self.curr_cw[0] = get_float(motor_pos)
        self.ui.curr_cw_pos.setText(str(np.round(self.curr_cw[0], 4)))
        try:
            self.curr_cw[1] = self.wl_regression(self.curr_cw[0])
            self.ui.curr_cw_wl.setText(str(np.round(self.curr_cw[1], 4)))
        except Exception as e:
            self.show_message("Could not perform wavelenght regression from motor pos.",e)
            self.ui.curr_cw_wl.setText("Error!")

    def update_iv_label(self, meas: tuple):
        """
        Updates IV labels with the floats in the given tuple, (Voltage,Current)

        :param meas: tuple -> (float,float), (Voltage,Current)
        :return:
        """

        (self.curr_V, self.curr_I) = meas
        current = self.curr_I * 1e6
        self.ui.curr_V.setText(str(np.round(self.curr_V, 4)))
        self.ui.curr_I.setText(str(np.round(current, 4)))

    def update_label(self, times: tuple):
        """
        Updates the elapsed, remaining time with the given times

        :param times: tuple-> (float,float), (elapsed, remaining)
        :return: No return
        """

        elapsed, remaining = times
        self.ui.time_label.setText(
            f"Elapsed / Remaining Time: {elapsed // 3600:02.0f}:{(elapsed % 3600) // 60:02.0f}:{elapsed % 60:02.0f}\
             /  {remaining // 3600:02.0f}:{(remaining % 3600) // 60:02.0f}:{remaining % 60:02.0f}")

    # =============================================================================
    # Starting measurement threads
    # =============================================================================

    def loop_linear(self):

        """
        Starts thread to loop michelson linear stage by the number of times give in the GUI
        between the values given by linear start/end in the GUI.

        If the motor is already looping stops it.

        """

        if self.check_text():
            try:
                if self.loop_thread.isRunning():
                    print("quit linear loop")
                    self.loop_thread.terminate = True
                    return
            except Exception as e:
                self.show_message("Error in looping Linear stage", e)

            if self.ui.trigger_linear.isChecked():

                self.loop_thread = Measurement_1D(self.com,
                                                  get_float(self.ui.linear_start.text()),
                                                  get_float(self.ui.linear_end.text()),
                                                  self.ui.linear_step.value() + 1,
                                                  self.lin_stage_name, 0)

                self.loop_thread.toggle_ui.connect(lambda: self.toggle_linear_loop)
                self.loop_thread.update_position.connect(self.ui.linear_cur_position.setText)
                self.loop_thread.update_progress.connect(self.ui.progressBar.setValue)
                self.loop_thread.toggle_button.connect(lambda: self.toggle_linear_loop_button)
                self.loop_thread.update_label.connect(self.update_label)
                self.loop_thread.device_disconnect.connect(self.show_error_msg)

                self.loop_thread.start()

            else:
                self.loop_thread = LoopLinear(self.com,
                                              self.ui.linear_spin.value(),
                                              get_float(
                                                  self.ui.linear_start.text()),
                                              get_float(self.ui.linear_end.text()))

                self.loop_thread.toggle_ui.connect(lambda: self.toggle_linear_loop)
                self.loop_thread.change_linear_position.connect(self.ui.linear_cur_position.setText)
                self.loop_thread.toggle_button.connect(lambda: self.toggle_linear_loop_button)
                self.loop_thread.device_disconnect.connect(self.show_error_msg)

                self.loop_thread.start()

    def loop_piezo(self):

        """
        Starts thread to loop michelson piezo stage by the number of times give in the GUI
        between the values given by linear start/end in the GUI.

        If the piezo is already looping stops it.

        """

        if self.check_text():

            try:
                if self.loop_thread.isRunning():
                    print("quit piezo loop")
                    self.loop_thread.terminate = True
                    return
            except Exception as e:
                self.show_message("Error in looping Michelson Piezo stage", e)

            if self.ui.trigger_piezo.isChecked():

                self.loop_thread = LoopPiezo(self.com,
                                             self.ui.piezo_spin.value(),
                                             get_float(self.ui.piezo_start.text()),
                                             get_float(self.ui.piezo_end.text()),
                                             self.ui.piezo_step.value())

                self.loop_thread.toggle_ui.connect(lambda: self.toggle_piezo_loop)
                self.loop_thread.change_piezo_position.connect(self.ui.piezo_cur_position.setText)
                self.loop_thread.toggle_button.connect(lambda: self.toggle_piezo_loop_button)
                self.loop_thread.device_disconnect.connect(self.show_error_msg)

                self.loop_thread.start()

            else:

                self.loop_thread = LoopPiezo(self.com,
                                             self.ui.piezo_spin.value(),
                                             get_float(self.ui.piezo_start.text()),
                                             get_float(self.ui.piezo_end.text()),
                                             self.ui.piezo_step.value())

                self.loop_thread.toggle_ui.connect(lambda: self.toggle_piezo_loop)
                self.loop_thread.change_piezo_position.connect(self.ui.piezo_cur_position.setText)
                self.loop_thread.toggle_button.connect(lambda: self.toggle_piezo_loop_button)
                self.loop_thread.device_disconnect.connect(self.show_error_msg)

                self.loop_thread.start()

    def measure_michelson(self):

        """
        Starts measurement thread with the MichelsonMeasurement function, takes values from GUI
        Sends number of spectras, filename, location and trigger setting to spectrometer through lightfield

        Updates position after each spectrum

        When finished, returns to initial positions and resets spectrometer to default

        :return: No return
        """

        if self.check_text():
            try:
                if self.measurement_thread.isRunning():
                    print("quit measurement")
                    self.measurement_thread.terminate = True

                    return
            except Exception as e:
                self.show_message("Error in Michelson Measurement.", e)

            self.measurement_thread = MichelsonMeasurement(self.com,
                                                           get_float(
                                                               self.ui.linear_start.text()),
                                                           get_float(
                                                               self.ui.linear_end.text()),
                                                           self.ui.linear_step.value() + 1,
                                                           get_float(
                                                               self.ui.piezo_start.text()),
                                                           get_float(
                                                               self.ui.piezo_end.text()),
                                                           self.ui.piezo_step.value() + 1,
                                                           self.michelson_equal_length)

            self.measurement_thread.toggle_ui.connect(self.toggle_measurement)
            self.measurement_thread.update_progress.connect(self.ui.progressBar.setValue)
            self.measurement_thread.change_linear_position.connect(self.ui.linear_cur_position.setText)
            self.measurement_thread.change_piezo_position.connect(self.ui.piezo_cur_position.setText)
            self.measurement_thread.toggle_button.connect(self.toggle_abort_button)
            self.measurement_thread.update_label.connect(self.update_label)
            self.measurement_thread.device_disconnect.connect(self.show_error_msg)

            self.measurement_thread.start()

    def measure_pol_sup(self, device_name: str, start: float, end: float, step: float, offset: float, starting_pos: float=None
                        , start_pos_sup: float= None):
        """
        Starts measurement thread with the Measurement_1D function

        Updates position after each spectrum

        When finished, returns to initial positions and resets spectrometer to default

        :return: No return
        """

        if self.check_text():
            try:
                if self.measurement_thread.isRunning():
                    print("quit measurement")
                    self.measurement_thread.terminate = True

                    return
            except Exception as e:
                self.show_message("Error in Measurement with {}".format(str(device_name)), e)

            excitation_laser = find_radio_button(self.ui.sources_group.buttons())
            if ("Green", "HeNe") in excitation_laser:
                pol_device_name = self.laser_power_name
                pol_sup_offset = self.curr_laser
            elif excitation_laser == "cw":
                pol_device_name = self.cw_hwp_name
                pol_sup_offset = self.cw_laser_offset
            else:
                self.show_message("Excitation source not found!")
                return

            self.measurement_thread = Measurement_pol_sup(self.com,
                                                     get_float(start),
                                                     get_float(end),
                                                     int(step),
                                                     device_name, offset,starting_pos,
                                                          pol_device_name,start_pos_sup,pol_sup_offset)

            self.measurement_thread.toggle_ui.connect(self.toggle_measurement)

            self.measurement_thread.update_position.connect(self.ui.curr_fss.setText)

            if device_name == self.laser_power_name:
                self.measurement_thread.update_pol_position.connect(self.ui.curr_laser.setText)
            elif device_name == self.cw_name:
                self.measurement_thread.update_pol_position.connect(self.ui.curr_cw_hwp.setText)


            self.measurement_thread.update_progress.connect(self.ui.progressBar.setValue)
            self.measurement_thread.toggle_button.connect(self.toggle_abort_button)
            self.measurement_thread.update_label.connect(self.update_label)
            self.measurement_thread.device_disconnect.connect(self.show_error_msg)

            self.measurement_thread.start()

    def measure_1d(self, device_name: str, start: float, end: float, step: float, offset: float, starting_pos: float=None):
        """
        Starts measurement thread with the Measurement_1D function

        Updates position after each spectrum

        When finished, returns to initial positions and resets spectrometer to default

        :return: No return
        """

        if self.check_text():
            try:
                if self.measurement_thread.isRunning():
                    print("quit measurement")
                    self.measurement_thread.terminate = True

                    return
            except Exception as e:
                self.show_message("Error in Measurement with {}".format(str(device_name)), e)

            self.measurement_thread = Measurement_1D(self.com,
                                                     get_float(start),
                                                     get_float(end),
                                                     int(step),
                                                     device_name, offset,starting_pos)

            self.measurement_thread.toggle_ui.connect(self.toggle_measurement)
            if device_name == self.piezo_name:
                self.measurement_thread.update_position.connect(self.ui.curr_z.setText)
            elif device_name == self.laser_power_name:
                self.measurement_thread.update_position.connect(self.ui.curr_laser.setText)
            elif device_name == self.fss_name:
                self.measurement_thread.update_position.connect(self.ui.curr_fss.setText)
            elif device_name == self.cw_name:
                self.measurement_thread.update_position.connect(self.update_cw_wl_label)
            # elif device_name == self.cw_hwp_name: #TODO: uncomment once we have hwp for cw
            #     self.measurement_thread.update_position.connect(self.ui.curr_cw_hwp.setText)

            self.measurement_thread.update_progress.connect(self.ui.progressBar.setValue)
            self.measurement_thread.toggle_button.connect(self.toggle_abort_button)
            self.measurement_thread.update_label.connect(self.update_label)
            self.measurement_thread.device_disconnect.connect(self.show_error_msg)

            self.measurement_thread.start()

    def measure_2d(self, device_name_slow: str, start_slow: float, end_slow: float, step_slow: float,
                   offset_slow: float,
                   device_name_fast: str, start_fast: float, end_fast: float, step_fast: float, offset_fast: float):
        """
        Starts measurement thread with the Measurement_2D function

        Updates position after each spectrum

        When finished, returns to initial positions and resets spectrometer to default

        :return: No return
        """

        if self.check_text():
            try:
                if self.measurement_thread.isRunning():
                    print("quit measurement")
                    self.measurement_thread.terminate = True

                    return
            except Exception as e:
                self.show_message(
                    "Error in Measurement with {} abd {}".format(str(device_name_slow), str(device_name_fast)), e)

            self.measurement_thread = Measurement_2D(self.com,
                                                     get_float(start_slow),
                                                     get_float(end_slow),
                                                     int(step_slow), offset_slow,
                                                     get_float(start_fast),
                                                     get_float(end_fast),
                                                     int(step_fast),
                                                     offset_fast,
                                                     device_name_slow,
                                                     device_name_fast
                                                     )

            self.measurement_thread.toggle_ui.connect(self.toggle_measurement)

            try:
                if device_name_slow == self.piezo_name:
                    self.measurement_thread.update_pos_slow.connect(self.ui.curr_z.setText)
                elif device_name_slow == self.laser_power_name:
                    self.measurement_thread.update_pos_slow.connect(self.ui.curr_laser.setText)
                elif device_name_slow == self.fss_name:
                    self.measurement_thread.update_pos_slow.connect(self.ui.curr_fss.setText)
                elif device_name_slow == self.smu_name:
                    self.measurement_thread.update_pos_slow.connect(self.update_iv_label)

                if device_name_slow in (self.piezo_name, self.laser_power_name, self.fss_name):
                    self.measurement_thread.update_pos_slow.connect(self.ui.curr_slow_1.setText)
                    self.ui.curr_slow_2.setText('None')
                else:
                    self.ui.curr_slow_1.setText(self.ui.curr_V.text())
                    self.ui.curr_slow_2.setText(self.ui.curr_I.text())

                if device_name_fast == self.piezo_name:
                    self.measurement_thread.update_pos_fast.connect(self.ui.curr_z.setText)
                elif device_name_fast == self.laser_power_name:
                    self.measurement_thread.update_pos_fast.connect(self.ui.curr_laser.setText)
                elif device_name_fast == self.fss_name:
                    self.measurement_thread.update_pos_fast.connect(self.ui.curr_fss.setText)
                elif device_name_fast == self.smu_name:
                    self.measurement_thread.update_pos_fast.connect(self.update_iv_label)

                if device_name_fast in (self.piezo_name, self.laser_power_name, self.fss_name):
                    self.measurement_thread.update_pos_fast.connect(self.ui.curr_fast_1.setText)
                    self.ui.curr_fast_2.setText('None')
                else:
                    self.ui.curr_fast.setText(self.ui.curr_V.text())
                    self.ui.curr_fast_2.setText(self.ui.curr_I.text())
            except Exception as e:
                self.show_message("Error in updating curr pos of 2d meas.", e)
                print("Error in updating curr pos of 2d meas.", e)

            self.measurement_thread.update_progress.connect(self.ui.progressBar.setValue)
            self.measurement_thread.toggle_button.connect(self.toggle_abort_button)
            self.measurement_thread.update_label.connect(self.update_label)
            self.measurement_thread.device_disconnect.connect(self.show_error_msg)

            self.measurement_thread.start()

    def volt_sweep(self, device_name: str, start: float, end: float, step: float):
        """
        Starts measurement thread with the Measurement_IV function
        Sends number of spectras, filename, location and trigger setting to spectrometer through lightfield

        Updates position after each spectrum

        When finished, returns to initial positions and resets spectrometer to default

        :return: No return
        """

        if self.check_text():
            try:
                if self.measurement_thread.isRunning():
                    print("quit measurement")
                    self.measurement_thread.terminate = True
                    return
            except Exception as e:
                self.show_message("Error in Voltage Sweep.", e)

            # self.change_frames(step)
            # self.toggle_trigger(True)

            self.measurement_thread = Measurement_IV(self.com,
                                                     get_float(start),
                                                     get_float(end),
                                                     int(step),
                                                     device_name)

            self.measurement_thread.toggle_ui.connect(self.toggle_measurement)
            self.measurement_thread.update_progress.connect(self.ui.progressBar.setValue)
            self.measurement_thread.update_position.connect(self.update_iv_label)
            self.measurement_thread.toggle_button.connect(self.toggle_abort_button)
            self.measurement_thread.update_label.connect(self.update_label)
            self.measurement_thread.device_disconnect.connect(self.show_error_msg)

            self.measurement_thread.start()

    def no_trigger_iv(self, device_name: str, start: float, end: float, frames: int, double: tuple = False):

        """
        Starts measurement thread with the Measure_IV_notrig function
        DOES NOT trigger spectrometer or sends any filename

        Updates position after each spectrum

        When finished, returns to initial positions and resets spectrometer to default

        :return: No return
        """

        self.i_v_list = [[], []]  # Voltage and current data
        self.i_v_list_neg = [[], []]  # Voltage and current data

        if self.check_text():
            try:
                if self.measurement_thread.isRunning():
                    print("quit measurement")
                    self.measurement_thread.terminate = True
                    return
            except Exception as e:
                self.show_message("Error in Measurement of IV curve", e)

            self.measurement_thread = Measurement_IV_notrig(self.com,
                                                            start,
                                                            end,
                                                            int(frames),
                                                            device_name)

            self.measurement_thread.toggle_ui.connect(self.toggle_measurement)
            self.measurement_thread.update_progress.connect(self.ui.progressBar.setValue)
            self.measurement_thread.update_position.connect(self.update_iv_label)
            self.measurement_thread.toggle_button.connect(self.toggle_abort_button)
            self.measurement_thread.update_label.connect(self.update_label)
            self.measurement_thread.device_disconnect.connect(self.show_error_msg)
            self.measurement_thread.return_iv_data.connect(self.update_iv_data)
            self.measurement_thread.process_done.connect(self.update_iv_graph)

            self.measurement_thread.start()

        if double:

            time.sleep(1)

            if self.check_text():
                try:
                    if self.measurement_thread.isRunning():
                        print("quit measurement")
                        self.measurement_thread.terminate = True
                        return
                except Exception as e:
                    self.show_message("Error in Measurement of IV curve", e)

                self.measurement_thread = Measurement_IV_notrig(self.com,
                                                                end,
                                                                start,
                                                                int(frames),
                                                                device_name)

                self.measurement_thread.toggle_ui.connect(self.toggle_measurement)
                self.measurement_thread.update_progress.connect(self.ui.progressBar.setValue)
                self.measurement_thread.update_position.connect(self.update_iv_label)
                self.measurement_thread.toggle_button.connect(self.toggle_abort_button)
                self.measurement_thread.update_label.connect(self.update_label)
                self.measurement_thread.device_disconnect.connect(self.show_error_msg)
                self.measurement_thread.return_iv_data.connect(self.update_iv_data)
                self.measurement_thread.process_done.connect(self.update_iv_graph)

                self.measurement_thread.start()

    def measure_xy(self, device_name_xy: str,
                   positions_array:list,starting_pos:list = None):
        """
        Starts measurement thread with the XYMeasurement function

        Z position is not fitted to the XY plane

        Updates position after each spectrum

        When finished, returns to initial positions and resets spectrometer to default

        :return: No return
        """

        if self.check_text():
            try:
                if self.measurement_thread.isRunning():
                    print("quit measurement")
                    self.measurement_thread.terminate = True

                    return
            except Exception as e:
                self.show_message("Error in space map measurement.", e)

            if self.ui.fast_map.isChecked():
                self.measurement_thread = XYMeasurement_fast(self.com,
                                                        positions_array,
                                                        device_name_xy, starting_pos,self.thread_frames)
            else:
                self.measurement_thread = XYMeasurement(self.com,
                                                        positions_array,
                                                        device_name_xy,starting_pos)

            self.measurement_thread.toggle_ui.connect(self.toggle_measurement)
            self.measurement_thread.update_progress.connect(self.ui.progressBar.setValue)
            self.measurement_thread.toggle_button.connect(self.toggle_abort_button)
            self.measurement_thread.update_label.connect(self.update_label)
            self.measurement_thread.device_disconnect.connect(self.show_error_msg)

            self.measurement_thread.update_pos1.connect(self.ui.curr_x.setText)
            self.measurement_thread.update_pos2.connect(self.ui.curr_y.setText)

            self.measurement_thread.start()

    def measure_meshed_xy(self, device_name_xy: str, device_name_z: str,
                          positions_array:list,starting_pos:list = None):

        """
        Starts measurement thread with the Meshed_xy_measurement function, takes values from GUI

        Z-position IS fitted to a plane in this function

        Updates position after each spectrum

        When finished, returns to initial positions and resets spectrometer to default

        :return: No return
        """

        if self.check_text():
            try:
                if self.measurement_thread.isRunning():
                    print("quit measurement")
                    self.measurement_thread.terminate = True

                    return
            except Exception as e:
                self.show_message("Error in Meshed space map measurement.", e)

            if self.ui.fast_map.isChecked():

                rad = get_float(self.ui.map_rad.text())

                self.measurement_thread = rad_limited_map(self.com,
                                                                positions_array,
                                                                device_name_xy, device_name_z,
                                                          starting_pos,rad)
            else:
                self.measurement_thread = Meshed_xy_measurement(self.com,
                                                                positions_array,
                                                                device_name_xy, device_name_z,starting_pos)

            self.measurement_thread.toggle_ui.connect(self.toggle_measurement)
            self.measurement_thread.update_progress.connect(self.ui.progressBar.setValue)
            self.measurement_thread.toggle_button.connect(self.toggle_abort_button)
            self.measurement_thread.update_label.connect(self.update_label)
            self.measurement_thread.device_disconnect.connect(self.show_error_msg)
            self.measurement_thread.update_pos1.connect(self.ui.curr_x.setText)
            self.measurement_thread.update_pos2.connect(self.ui.curr_y.setText)
            self.measurement_thread.update_pos3.connect(self.ui.curr_z.setText)

            self.measurement_thread.start()

    # =============================================================================
    #   Setup measurement and start/abort thread
    # =============================================================================

    def abort_current_measurement(self):

        """
        Aborts current measurement, kills the measurement thread, request current positions of stages and write it to GUI.
        Resets spectrometer to standard settings.

        :return: No return
        """

        if self.running_measurement == 'space_map':
            self.space_map()
        elif self.running_measurement == 'line_scan':
            self.line_scan()
        elif self.running_measurement == "power_map":
            self.power_map(laser="Green")
        elif self.running_measurement == "power_map_3900s":
            self.power_map(laser="3900s")
        elif self.running_measurement == 'kethley':
            self.voltage_sweep()
        elif self.running_measurement == 'fss':
            self.fss_measurement()
        elif self.running_measurement == 'michelson':
            self.measure_michelson()
        elif self.running_measurement == 'ple':
            self.ple_meas()

        self.toggle_ui(True)

    def fss_measurement(self):

        """
        Send settings to function to start FSS measurement (gets info from GUI) or stops measurement,
        kills thread and restore initial positions if thread is running.

        Sends number of spectras, filename, location and trigger setting to spectrometer through lightfield or resets it.


        :return: No return
        """

        if self.running_measurement not in ('fss', None):
            return

        frames = get_float(self.ui.fss_frames.text())
        frames = int(np.round(frames, 0))
        starting_angle = get_float(self.ui.fss_i.text()) + self.fss_offset
        end_angle = get_float(self.ui.fss_end.text()) + self.fss_offset
        end_angle = get_float(self.ui.fss_end.text()) + self.fss_offset

        if self.running_measurement == 'fss':
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)

        elif self.running_measurement is None:
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.2)
                self.change_frames(frames)
                time.sleep(0.1)
                self.toggle_trigger(True)
                time.sleep(0.1)
                self.ui.fss_measure.setChecked(True)
                time.sleep(0.1)
                self.update_filename_spectrometer()
                time.sleep(0.1)
                self.com.request(self.spectrometer_name, "acquire_file")
                self.running_measurement = 'fss'
                print("Starting experiment: ", self.running_measurement)
        if self.ui.pol_sup_check.isChecked():
            self.measure_pol_sup(self.fss_name, starting_angle, end_angle, frames, self.fss_offset, self.curr_fss
            , self.curr_cw_hwp)
        else:
            self.measure_1d(self.fss_name, starting_angle, end_angle, frames, self.fss_offset,self.curr_fss)

    def update_equal_length_path_value(self):

        linear_connected, self.curr_linear_pos = self.com.request(self.lin_stage_name, "get_pos")
        self.michelson_equal_length = self.curr_linear_pos
        self.ui.equal_path_length.setText(str(self.michelson_equal_length))

        return

    def michelson_measurement(self):

        """
        Calls function to start michelson measurement or to kill the thread.
        restore initial positions if thread is running.

        Sends number of spectras, filename, location and trigger setting to spectrometer through lightfield
        or resets it.


        :return: No return
        """

        frames = (self.ui.linear_step.value() + 1) * (self.ui.piezo_step.value() + 1)


        if self.running_measurement not in ('michelson', None): return

        if self.running_measurement == 'michelson':
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")


        elif self.running_measurement is None:
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)
                self.change_frames(frames)
                time.sleep(0.1)
                self.toggle_trigger(True)
                time.sleep(0.1)
                self.ui.michelson_int.setChecked(True)
                time.sleep(0.1)
                self.update_filename_spectrometer()
                time.sleep(0.1)
                self.com.request(self.spectrometer_name, "acquire_file")
                self.running_measurement = 'michelson'
                print("Starting experiment: ", self.running_measurement)

        self.measure_michelson()

    def meas_cw_power(self):

        frames = int(np.round(get_float(self.ui.ple_frames.text()), 0))
        start_wl = get_float(self.ui.ple_start.text())
        end_wl = get_float(self.ui.ple_end.text())

        start_motor_pos = end_motor_pos = self.curr_cw[0]

        start_motor_pos = self.wl_regression(start_wl)
        end_motor_pos = self.wl_regression(end_wl)

        motor_pos_list = np.linspace(start_motor_pos,end_motor_pos,frames)

        pwr_readings = []

        if self.pwr_meter_name not in self.device_name_list: return
        if self.cw_name not in self.device_name_list: return

        for pos in motor_pos_list:
            try:
                status, curr_motor_pos = self.com.request(self.cw_name, "set_pos",pos)
                self.update_cw_wl_label(curr_motor_pos)

                status, pwr_read = self.com.request(self.pwr_meter_name, "read")
                pwr_readings.append([self.wl_regression(curr_motor_pos),pwr_read])
            except Exception as e:
                self.show_message("Error in reading laser power due to ",e)

        today = datetime.today()

        year = today.strftime("%Y")
        month = today.strftime("%m")
        day = today.strftime("%d")

        hour = today.strftime("%H")
        minute = today.strftime("%M")
        second = today.strftime("%S")

        base_path = os.path.join(self.path, year, month, str(day), str(self.ui.sample_name.text()).strip())
        filename = "cw_power_read_{}_{}_{}_{.2f}nm_to_{.2f}_millenia_power_2750mW.txt".format(day,month,year,start_wl,end_wl)

        try:
            with open(os.path.join(base_path, filename), "w+") as f:
                f.write("Wavelength(nm);Power Read (mW)")
                for item in pwr_readings:
                    f.write("{};{}".format(item[0],item[1]))

        except:
            self.show_message("Could not save power reading to file")
            print(pwr_readings)

    def ple_meas(self):
        """
        Send settings to function to start PLE measurement (gets info from GUI) or stops measurement,
        kills thread and restore initial positions if thread is running.

        Sends number of spectras, filename, location and trigger setting to spectrometer through lightfield
        or resets it.


        :return: No return
        """

        frames = int(np.round(get_float(self.ui.ple_frames.text()), 0))
        start_wl = get_float(self.ui.ple_start.text())
        end_wl = get_float(self.ui.ple_end.text())

        start_pos = end_pos = self.curr_cw[0]

        if self.running_measurement not in ('ple', None): return

        if self.running_measurement == 'ple':
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)

        elif self.running_measurement is None:

            try:
                start_pos = self.wl_regression(end_wl)
                end_pos = self.wl_regression(start_wl)
                print(start_pos,end_pos)
            except Exception as e:
                self.show_message("Could not find motor positions for PLE.", e)
                return

            self.com.request(self.cw_name,"set_pos",start_pos)

            time.sleep(2)

            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)
                self.change_frames(frames)
                time.sleep(0.1)
                self.toggle_trigger(True)
                time.sleep(0.1)
                self.ui.ple_rb.setChecked(True)
                time.sleep(0.1)
                self.update_filename_spectrometer()
                time.sleep(0.1)
                self.com.request(self.spectrometer_name, "acquire_file")
                self.running_measurement = 'ple'
                print("Starting experiment: ", self.running_measurement)

        self.measure_1d(self.cw_name, start_pos, end_pos, frames, 0,self.curr_cw[0])

    def power_map(self,laser:str = None):

        """
        Send settings to function to start power mapping measurement (gets info from GUI) or stops measurement,
        kills thread and restore initial positions if thread is running.

        Sends number of spectras, filename, location and trigger setting to spectrometer through lightfield or resets it.
        :params laser: Name of the laser used or None for Green/HeNe

        :return: No return
        """

        if laser == "Green":
            frames = int(np.round(get_float(self.ui.pow_map_frames.text()), 0))
            starting_angle = get_float(self.ui.power_map_i.text()) + self.laser_offset
            end_angle = get_float(self.ui.power_map_f.text()) + self.laser_offset
        elif laser == "3900s":
            return #TODO remove once cw has a hwp
            frames = int(np.round(get_float(self.ui.cw_pow_map_frames.text()), 0))
            starting_angle = get_float(self.ui.cw_power_map_i.text()) + self.cw_laser_offset
            end_angle = get_float(self.ui.cw_power_map_f.text()) + self.cw_laser_offset


        if self.running_measurement not in ('power_map','power_map_3900s', None):
            return


        if self.running_measurement is not None and 'power_map' in self.running_measurement:
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)

        elif self.running_measurement is None:
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)
                self.change_frames(frames)
                time.sleep(0.1)
                self.toggle_trigger(True)
                time.sleep(0.1)
                self.ui.power_map.setChecked(True)
                time.sleep(0.1)
                self.update_filename_spectrometer()
                time.sleep(0.1)
                self.com.request(self.spectrometer_name, "acquire_file")
                if laser == "Green":
                    self.running_measurement = 'power_map'
                elif laser == '3900s':
                    self.running_measurement = 'power_map_3900s'
                print("Starting experiment: ", self.running_measurement)

        if laser == "Green":
            self.measure_1d(self.laser_power_name, starting_angle, end_angle, frames, self.laser_offset,self.curr_laser)
        elif laser == '3900s':
            #self.measure_1d(self.cw_hwp_name, starting_angle, end_angle, frames, self.cw_laser_offset,self.curr_cw_hwp)
            pass #TODO remove once cw has a hwp
    def setup_2d_meas(self):
        """
        Gets values from UI and setup 2d measurement

        :return: No return
        """

        def find_device_name(button_group):
            button = find_radio_button(button_group.buttons())

            if button is None:
                return 'break'

            elif button.lower().strip() == 'polarization':
                device = self.fss_name
            elif button.lower().strip() == 'laser power':
                device = self.laser_power_name
            elif button.lower().strip() == 'focus':
                device = self.piezo_name
            elif button.lower().strip() == 'voltage':
                device = self.smu_name
            elif button.lower().strip() == 'ple':
                device = self.cw_name
            elif button.lower().strip() == '3900s Power':
                device = self.cw_hwp_name

            return device

        device_name_slow = find_device_name(self.ui.slow_axis_radio)
        device_name_fast = find_device_name(self.ui.fast_axis_radio)

        if device_name_slow == self.laser_power_name:
            offset_slow = self.laser_offset
        elif device_name_slow == self.fss_name:
            offset_slow = self.fss_offset
        elif device_name_slow == self.cw_hwp_name:
            offset_slow = self.cw_laser_offset
        else:
            offset_slow = 0

        if device_name_fast == self.laser_power_name:
            offset_fast = self.laser_offset
        elif device_name_fast == self.fss_name:
            offset_fast = self.fss_offset
        elif device_name_fast == self.cw_hwp_name:
            offset_fast = self.cw_laser_offset
        else:
            offset_fast = 0

        start_slow = get_float(self.ui.slow_axis_start.text()) + offset_slow
        slow_frames = get_float(self.ui.slow_axis_frames.text())
        end_slow = get_float(self.ui.slow_axis_end.text()) + offset_slow
        step_slow = get_float(self.ui.slow_axis_steps.text())

        if device_name_slow == self.cw_name:
            try:
                start_pos = self.wl_regression(start_slow)
                end_pos = self.wl_regression(end_slow)
            except Exception as e:
                self.show_message("Could not find motor positions for PLE.", e)
                return

        start_fast = get_float(self.ui.fast_axis_start.text()) + offset_fast
        end_fast = get_float(self.ui.fast_axis_end.text())
        step_fast = get_float(self.ui.fast_axis_steps.text()) + offset_fast
        fast_frames = get_float(self.ui.fast_axis_frames.text())

        if device_name_fast == self.cw_name:
            try:
                start_pos = self.wl_regression(start_fast)
                end_pos = self.wl_regression(end_fast)
            except Exception as e:
                self.show_message("Could not find motor positions for PLE.", e)
                return

        frames = fast_frames * slow_frames

        if device_name_slow is None or device_name_fast is None:
            return

        if self.running_measurement not in ('2d_meas', None): return

        if self.running_measurement == '2d_meas':
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)

        elif self.running_measurement is None:
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)
                self.change_frames(frames)
                time.sleep(0.1)
                self.toggle_trigger(True)
                time.sleep(0.1)
                self.ui.meas_2d.setChecked(True)
                time.sleep(0.1)
                self.update_filename_spectrometer()
                time.sleep(0.1)
                self.com.request(self.spectrometer_name, "acquire_file")
                self.running_measurement = '2d_meas'
                print("Starting experiment: ", self.running_measurement)

        self.measure_2d(device_name_slow, start_slow, end_slow, slow_frames, offset_slow,
                        device_name_fast, start_fast, end_fast, fast_frames, offset_fast)

    def meas_double_IV_curve(self):
        def meas_iv_curve(self):
            """
            Send settings to function to start double iv curve (gets info from GUI) or stops measurement,
            kills thread and restore initial positions if thread is running.

            Spectrometer is not involved in this measurement

            :return: No return
            """

            v_max = get_float(self.ui.volt_max.text())
            v_min = get_float(self.ui.volt_min.text())
            frames = get_float(self.ui.volt_frames.text())
            frames += 1

            frames = int(np.round(frames, 0))

            self.no_trigger_iv(self.smu_name, v_min, v_max, frames, True)

    def meas_iv_curve(self):

        """
        Send settings to function to start iv curve (gets info from GUI) or stops measurement,
        kills thread and restore initial positions if thread is running.

        Spectrometer is not involved in this measurement

        :return: No return
        """

        v_max = get_float(self.ui.volt_max.text())
        v_min = get_float(self.ui.volt_min.text())
        frames = get_float(self.ui.volt_frames.text())
        frames += 1

        frames = int(np.round(frames, 0))

        self.no_trigger_iv(self.smu_name, v_min, v_max, frames)

    def voltage_sweep(self):

        """
        Send settings to function to start voltage sweep measurement (gets info from GUI) or stops measurement,
        kills thread and restore initial positions if thread is running.

        Sends number of spectras, filename, location and trigger setting to spectrometer through lightfield
        or resets it.


        :return: No return
        """

        v_max = get_float(self.ui.volt_max.text())
        v_min = get_float(self.ui.volt_min.text())

        frames = get_float(self.ui.volt_frames.text())
        frames = int(np.round(frames, 0))

        if self.running_measurement not in ('kethley', None): return

        if self.running_measurement == 'kethley':
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)

        elif self.running_measurement is None:
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)
                self.change_frames(frames)
                time.sleep(0.1)
                self.toggle_trigger(True)
                time.sleep(0.1)
                self.ui.volt_sweep.setChecked(True)
                time.sleep(0.1)
                self.update_filename_spectrometer()
                time.sleep(0.1)
                self.com.request(self.spectrometer_name, "acquire_file")
                time.sleep(0.1)
                self.running_measurement = 'kethley'
                print("Starting experiment: ", self.running_measurement)

        self.volt_sweep(self.smu_name, v_min, v_max, frames)

    def line_scan(self):

        """
        Send settings to function to start space map (gets info from GUI) or stops measurement,
        kills thread and restore initial positions if thread is running.

        Sends number of spectras, filename, location and trigger setting to spectrometer through lightfield
        or resets it.

        if plane fit checkbox is checked will fit a plane equation to the three given points and extrapolate focus positions

        :return: No return
        """

        npoints = get_float(self.ui.line_scan_points.text())
        frames = int(npoints)

        self.com.request(self.xy_stage_name, 'set_acc', 30)

        if self.running_measurement not in ('line_scan', None): return

        if self.running_measurement == 'line_scan':
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)

                if self.ui.fit_plane.isChecked():
                    positions_array = [(self.start_x, self.start_y, 50)]

                    self.measure_meshed_xy(self.xy_stage_name, self.piezo_name, positions_array,[self.start_x,self.start_y])
                else:
                    self.measure_xy(self.xy_stage_name, [self.start_x, self.start_y],[self.start_x,self.start_y])

        elif self.running_measurement is None:
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)

                time.sleep(0.1)
                self.change_frames(frames)
                time.sleep(0.1)
                self.toggle_trigger(True)
                time.sleep(0.1)
                self.ui.line_scan.setChecked(True)
                time.sleep(0.1)
                self.update_filename_spectrometer()
                time.sleep(0.1)
                self.com.request(self.spectrometer_name, "acquire_file")
                self.running_measurement = 'line_scan'
                print("Starting experiment: ", self.running_measurement)

                positions_array, start_x, start_y, start_z = self.check_spacemap_points(line = True)

                self.make_move(1, 'xy', -1, goto=[start_x, start_y])

                if self.ui.line_scan_z.isChecked():
                    self.measure_1d(self.piezo_name, start_z, positions_array, frames - 1, 0, self.curr_z)

                else:
                    try:
                        if self.ui.fit_plane.isChecked():
                            self.measure_meshed_xy(self.xy_stage_name, self.piezo_name, positions_array,
                                                   [self.start_x,self.start_y,self.start_z])

                        else:
                            self.measure_xy(self.xy_stage_name, positions_array,[self.start_x,self.start_y])

                    except Exception as e:
                        self.show_message('Error in device Handling or quit measurement', e)

                        if self.ui.fit_plane.isChecked():
                            positions_array = [(self.start_x, self.start_y, self.curr_z)]

                            self.measure_meshed_xy(self.xy_stage_name, self.piezo_name, positions_array,[self.start_x,self.start_y,self.start_z])
                        else:
                            self.measure_xy(self.xy_stage_name, [self.start_x, self.start_y],[self.start_x,self.start_y])

    # noinspection PyTypeChecker

    def check_spacemap_points(self, line=False):
        ncol = get_float(self.ui.space_map_colums.text())
        nlines = get_float(self.ui.space_map_lines.text())
        step_x = get_float(self.ui.column_step.text())
        step_y = get_float(self.ui.line_step.text())

        def z_pos(x_coord, y_coord, eq_factors):

            eq_factors = [float(_) for _ in eq_factors]

            try:
                z_coord = (eq_factors[0] * x_coord + eq_factors[1] * y_coord + eq_factors[3]) / (-1 * eq_factors[2])
            except Exception as error:
                self.show_message("Error in calculating z_coord positions for meshed measurement.", error)
                z_coord = (0 * x_coord + 0 * y_coord + float(self.curr_z)) / 1

            if z_coord > 97e-3:
                return 97
            elif z_coord < 5e-3:
                return 5
            else:
                z_coord = 1e3 * z_coord
            return z_coord

        status, (self.start_x, self.start_y) = self.com.request(self.xy_stage_name, "get_pos")
        status, self.start_z = self.com.request(self.piezo_name, "get_pos")
        starting_z = self.start_z

        # assume measurement start at center and move to one edge
        if not line:
            lat_size = nlines * step_x * 1e-3
            hor_size = ncol * step_y * 1e-3

            start_x = self.start_x - lat_size / 2
            start_y = self.start_y - hor_size / 2
            end_x = self.start_x + lat_size / 2
            end_y = self.start_y + hor_size / 2

        else:
            start_x = end_x = self.start_x
            start_y = end_y = self.start_y
            start_z = end_z = self.start_z
            nlines = ncol = 1
            npoints = get_float(self.ui.line_scan_points.text())
            step = get_float(self.ui.line_scan_step.text())
            line_size = npoints * step * 1e-3

            if self.ui.line_scan_x.isChecked():
                start_x = self.start_x - line_size / 2
                end_x = self.start_x + line_size / 2
                frames = ncol
            elif self.ui.line_scan_y.isChecked():
                start_y = self.start_y - line_size / 2
                end_y = self.start_y + line_size / 2
                frames = nlines
            elif self.ui.line_scan_z.isChecked():
                if float(line_size) > 100 - self.start_z:
                    line_size = 100 - self.start_z
                start_z = self.start_z - line_size / 2
                end_z = self.start_z + line_size / 2

                return end_z, start_x, start_y, start_z

        try:

            x_array = np.linspace(start_x, end_x, int(ncol))
            y_array = np.linspace(start_y, end_y, int(nlines))
            z_array = np.array([])
            positions_array = []

            if self.ui.fit_plane.isChecked():

                status, starting_z = self.com.request(self.piezo_name, "get_pos")

                points = [[], [], []]

                def equation_plane(p1, p2, p3):

                    x1, y1, z1 = p1[0], p1[1], p1[2] * 1e-3
                    x2, y2, z2 = p2[0], p2[1], p2[2] * 1e-3
                    x3, y3, z3 = p3[0], p3[1], p3[2] * 1e-3

                    a1 = x2 - x1
                    b1 = y2 - y1
                    c1 = z2 - z1
                    a2 = x3 - x1
                    b2 = y3 - y1
                    c2 = z3 - z1
                    a = b1 * c2 - b2 * c1
                    b = a2 * c1 - a1 * c2
                    c = a1 * b2 - b1 * a2
                    d = (- a * x1 - b * y1 - c * z1)

                    params = (a, b, c, d)

                    return params

                rowcount = self.ui.plane_positions.rowCount()

                for _ in np.arange(rowcount):

                    if not self.ui.plane_positions.verticalHeaderItem(_):
                        points[_].append(float(self.ui.plane_positions.item(_, 0).data(0)))
                        points[_].append(float(self.ui.plane_positions.item(_, 1).data(0)))
                        points[_].append(float(self.ui.plane_positions.item(_, 2).data(0)))

                try:
                    eq_factors = equation_plane(points[0], points[1], points[2])
                except Exception as e:
                    self.show_message("Error in plane equation.", e)
                    eq_factors = [0, 0, -1, self.curr_z]

                if self.ui.snake_like.isChecked():

                    for y in y_array:
                        for x in x_array:
                            z = z_pos(x, y, eq_factors)
                            positions_array.append((x, y, z))
                            z_array = np.append(z_array, z)
                        x_array = np.flip(x_array)

                else:
                    for y in y_array:
                        for x in x_array:
                            z = z_pos(x, y, eq_factors)
                            positions_array.append((x, y, z))
                            z_array = np.append(z_array, z)

            else:

                if self.ui.snake_like.isChecked():
                    for y in y_array:
                        for x in x_array:
                            positions_array.append((x, y))
                            x_array = np.flip(x_array)
                            z_array = np.append(z_array, starting_z)

                else:
                    for y in y_array:
                        for x in x_array:
                            positions_array.append((x, y))
                            z_array = np.append(z_array, starting_z)

            # Set labels
            min_x, max_x = min(x_array), max(x_array)
            min_y, max_y = min(y_array), max(y_array)
            if len(z_array) > 2:
                min_z, max_z = min(z_array), max(z_array)
            else:
                min_z = max_z = "unchanged"

            # Create plot mesh:
            z_array_plot = np.array([])
            positions_array_plot = []

            # For plotting plane
            if self.ui.fit_plane.isChecked():
                for y in y_array:
                    for x in x_array:
                        z = z_pos(x, y, eq_factors)
                        positions_array_plot.append((x, y, z))
                        z_array_plot = np.append(z_array_plot, z)
            else:
                for y in y_array:
                    for x in x_array:
                        positions_array_plot.append((x, y, starting_z))
                        z_array_plot = np.append(z_array_plot, starting_z)

            try:
                self.ui.start_x.setText(str(round(min_x, 4)))
                self.ui.start_y.setText(str(round(min_y, 4)))
                self.ui.start_z.setText(str(round(min_z, 4)))
                self.ui.end_x.setText(str(round(max_x, 4)))
                self.ui.end_y.setText(str(round(max_y, 4)))
                self.ui.end_z.setText(str(round(max_z, 4)))

            except:
                self.ui.start_x.setText("ERROR")
                self.ui.start_y.setText("ERROR")
                self.ui.start_z.setText("ERROR")
                self.ui.end_x.setText("ERROR")
                self.ui.end_y.setText("ERROR")
                self.ui.end_z.setText("ERROR")

            X, Y = np.meshgrid(x_array, y_array)
            Z = z_array_plot.reshape(len(x_array), len(y_array))

            self.set_label_filename()

            today = date.today()

            year = today.strftime("%Y")
            month = today.strftime("%m")
            day = today.strftime("%d")

            path = os.path.join(self.path, year, month, str(day), str(self.ui.sample_name.text()).strip())

            filename = str(self.ui.label_filename.text()).replace(".spe", "_plane_fit.png")

            filename = os.path.join(path, filename)

            plane_fig, ax_plane = plt.subplots(subplot_kw={'projection': '3d'})
            ax_plane.plot_surface(X, Y, Z, rstride=1, cstride=1,
                                  cmap='viridis', edgecolor='none')
            ax_plane.set_title('Wafe plane fit')
            plane_fig.savefig(filename)

            return positions_array, start_x, start_y, starting_z

        except Exception as e:
            self.show_message("Error in space map points", e)

    def space_map(self):

        """
        Send settings to function to start space map (gets info from GUI) or stops measurement,
        kills thread and restore initial positions if thread is running.

        Sends number of spectras, filename, location and trigger setting to spectrometer through lightfield
        or resets it.

        if plane fit checkbox is checked will fit a plane equation to the three given points and extrapolate focus positions

        :return: No return
        """

        ncol = get_float(self.ui.space_map_colums.text())
        nlines = get_float(self.ui.space_map_lines.text())
        frames = ncol*nlines

        self.com.request(self.xy_stage_name, 'set_acc', 30)

        if self.running_measurement not in ('space_map', None): return

        if self.running_measurement == 'space_map':
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)

                if self.ui.fit_plane.isChecked():
                    positions_array = [(self.start_x, self.start_y, 50)]

                    self.measure_meshed_xy(self.xy_stage_name, self.piezo_name, positions_array,[self.start_x,self.start_y])
                else:
                    self.measure_xy(self.xy_stage_name, [self.start_x, self.start_y],[self.start_x,self.start_y])

        elif self.running_measurement is None:
            if self.spectrometer_connected:
                self.com.request(self.spectrometer_name, "stop")
                time.sleep(0.1)

                status, (self.start_x, self.start_y) = self.com.request(self.xy_stage_name, "get_pos")

                # assume measurement start at center and move to one edge

                time.sleep(0.1)
                self.change_frames(int(frames))
                time.sleep(0.1)
                self.toggle_trigger(True)
                time.sleep(0.1)
                self.ui.space_map.setChecked(True)
                time.sleep(0.1)
                self.update_filename_spectrometer()
                time.sleep(0.1)
                self.com.request(self.spectrometer_name, "acquire_file")
                self.running_measurement = 'space_map'
                print("Starting experiment: ", self.running_measurement)

                try:

                    positions_array, start_x, start_y, starting_z = self.check_spacemap_points()

                    self.make_move(1, 'xy', -1, goto=[start_x, start_y])

                    if self.ui.fit_plane.isChecked():

                        self.measure_meshed_xy(self.xy_stage_name, self.piezo_name, positions_array,[self.start_x,self.start_y,starting_z])

                    else:
                        self.measure_xy(self.xy_stage_name, positions_array,
                                        [self.start_x,self.start_y])

                except Exception as e:
                    self.show_message('Error in device Handling or quit measurement', e)

                    if self.ui.fit_plane.isChecked():
                        positions_array = [(self.start_x, self.start_y, self.curr_z)]

                        self.measure_meshed_xy(self.xy_stage_name, self.piezo_name, positions_array,[self.start_x,self.start_y,starting_z])
                    else:
                        self.measure_xy(self.xy_stage_name, [self.start_x, self.start_y],[self.start_x,self.start_y])


def get_float(variable: str | float):
    """
    Converts string or float to float, usually to change "," to "."

    :param variable: variable to be converted
    :return: float of variable
    """

    if isinstance(variable, float): return variable
    elif isinstance(variable, np.floating): return variable

    negative = False

    if "-" in str(variable): negative = True

    variable = str(variable)
    variable = variable.strip().replace(",", ".").replace(",", ".").replace(",", ".")
    variable = re.sub("[^0-9,\.]", "", variable)
    variable = variable.strip().replace(",", ".").replace(",", ".").replace(",", ".")

    if "-" in variable:pass
    elif negative: variable = "-"+variable

    try:
        float_variable = float(variable)
        return float_variable

    except Exception as e:
        print("Error in converting to {} float".format(variable), e)
        try:
            float_variable = float(variable.replace("-",""))
        except Exception as e:
            print("Error in converting to {} float".format(variable), e)
            return float(0.01)

def find_radio_button(group):
    """
    Returns the label of the selected radio button within the GUI group

    :param group: group of radio buttons from GUI
    :return: text from selected radiobutton in the given group
    """

    box_elements = group

    checked_on_rb = None

    radio_buttons = [
        elem for elem in box_elements if isinstance(elem, QRadioButton)]

    for rb in radio_buttons:
        if rb.isChecked():
            checked_on_rb = rb.text()

    return checked_on_rb


def file_dialog(title: str, for_open: bool = True, fmt: str = 'conf', is_folder: bool = False, open_path: str = None):
    """
    Dialog to open prompt and choose file

    :param title: Window title
    :param for_open: File is going ot be open
    :param fmt: File format
    :param is_folder: Only choose directory
    :type title: str
    :type for_open: bool
    :type fmt: str
    :type is_folder: bool
    :return:
    """

    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog
    options |= QFileDialog.DontUseCustomDirectoryIcons
    dialog = QFileDialog()
    dialog.setOptions(options)
    dialog.setWindowTitle(title)
    dialog.setFilter(dialog.filter() | QtCore.QDir.Hidden)

    # ARE WE TALKING ABOUT FILES OR FOLDERS
    if is_folder:
        dialog.setFileMode(QFileDialog.DirectoryOnly)
    else:
        dialog.setFileMode(QFileDialog.AnyFile)
    # OPENING OR SAVING
    dialog.setAcceptMode(QFileDialog.AcceptOpen) if for_open else dialog.setAcceptMode(
        QFileDialog.AcceptSave)

    # SET FORMAT, IF SPECIFIED
    if fmt != '' and is_folder is False:
        dialog.setDefaultSuffix(fmt)
        dialog.setNameFilters([f'{fmt} (*.{fmt})'])

    # SET THE STARTING DIRECTORY
    if open_path is None:
        dialog.setDirectory(os.path.join(os.getcwd(), "configs"))
    else:
        dialog.setDirectory(open_path)

    if dialog.exec_() == QDialog.Accepted:
        path = dialog.selectedFiles()[0]  # returns a list
        return path
    else:
        return None

if __name__ == "__main__":
    # parser=argparse.ArgumentParser()
    # # Do not connect to device by default
    # parser.add_argument("--noconnect", action="store_true")
    # args=parser.parse_args()

    Pl_Control()