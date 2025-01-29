import numpy as np
import os
import sys
from pathlib import Path
from datetime import datetime
import cv2
from PIL import Image
import math
from skspatial.objects import Plane, Points
from source.sample_structure.sample import Sample
from PyQt5.QtCore import QObject, pyqtSignal
from source.devices.Camera import Camera
from source.devices.Cryo_XY_Stage import CryoStage
from source.devices.Imaging_Piezo import ImagingPiezo
from source.devices.Linear_Stage import LinearStage
from source.devices.Waveplate import Waveplate
from source.devices.Spectrometer import Spectrometer
from source.qd_detection.model.wl_regression import WLRegression
import torch
from source.utils.utils import standardize


from source.utils.qd_detection import get_qd_positions_from_image
from source.utils.marker_detection import find_marker_origin_conventional, square_image

class SampleScanner(QObject):
    error_signal = pyqtSignal(tuple)
    device_available_signal = pyqtSignal(str, bool)
    progress_signal = pyqtSignal(int)
    triangpos_parameter_signal = pyqtSignal(tuple)
    exciton_wl_signal = pyqtSignal(int)
    finished_signal = pyqtSignal()

    def __init__(self, cam: Camera, xy_stage: CryoStage, imaging_piezo: ImagingPiezo, lin_stage: LinearStage,
                 spectrometer: Spectrometer, waveplate: Waveplate, wl_regression, image_path):
        super().__init__()
        self.name = "SampleScanner"
        self.image_path = image_path
        self.device_available = True

        self.sample = None

        # // improve this by using com.devices?
        self.camera = cam
        self.imaging_piezo = imaging_piezo
        self.lin_stage = lin_stage
        self.xy_stage = xy_stage
        self.spectrometer = spectrometer
        self.waveplate = waveplate
        self.wl_regressor = wl_regression

    def add_sample(self, sample: Sample):
        self.device_available_signal.emit(self.name, False)
        print("In adding Sample to Scanner")
        self.sample = sample
        self.device_available_signal.emit(self.name, True)
        self.finished_signal.emit()
        return

    def __calc_image_metrics__(self, image):
        # calculate histogram, mean & var
        hist = cv2.calcHist(images=[image], channels=[0], mask=None, histSize=[self.imaging_piezo.bins],
                            ranges=[0, 2 ** 16])
        mean = 0
        for j in range(self.imaging_piezo.bins):
            mean = mean + (j + 1) * hist[j] / np.sum(hist)
        var = 0
        for j in range(self.imaging_piezo.bins):
            var = var + ((j + 1) - mean) ** 2 * hist[j] / np.sum(hist)
        std = np.sqrt(var)

        return std

    def autofocus(self, standalone_call=True, keep_img=False):
        """
        standalone_call = True, if the function is called outside of the scanning / saving functions
        """

        # TODO: Make it faster? I remember Lars said it was slow. Maybe this can be done by measuring images in larger steps
        # #TODO: around the center (current pos) then analyzing the images and moving to smaller steps
        self.__set_signal_block_status__(True)

        if standalone_call:
            self.device_available_signal.emit(self.name, False)

        if keep_img:
            # create temporary folder for images, if necessary or desired
            path = Path.cwd()
            temp_path = Path(str(path) + "/temp_focus")
            os.makedirs(temp_path, exist_ok=True)
            os.chdir(temp_path)
            with open('piezo_pos.txt', 'w') as f:
                f.write("")

        # go through all piezo steps & capture images
        start_pos = 0
        optimal_piezo = 0
        optimal_std = 0

        # optimal piezo position is found in 3 different magnitudes (steps of 10, 1, 0.1)
        for magnitude in range(3):
            end_pos = min((self.imaging_piezo.max_piezo_steps // (10 ** magnitude)) + start_pos, self.imaging_piezo.max_piezo_steps)
            step_size = self.imaging_piezo.max_piezo_steps / (10 ** (magnitude + 1))
            print(f"Start: {start_pos}, End: {end_pos}, Step: {step_size}")
            piezo_values = []
            # break_flag = False
            for i, piezo_pos in enumerate(np.arange(start_pos, end_pos + 1, step_size)):
                # move piezo to next position
                status = self.imaging_piezo.set_piezo_pos(piezo_pos)
                if not status:
                    # error signal has been emitted by piezo already
                    self.__set_signal_block_status__(False)
                    if standalone_call:
                        self.device_available_signal.emit(self.name, True)
                        self.finished_signal.emit()
                    return None, None

                # capture image
                original_image = self.camera.capture_image()
                if original_image is None:
                    # error signal has been emitted by camera already
                    self.__set_signal_block_status__(False)
                    if standalone_call:
                        self.device_available_signal.emit(self.name, True)
                        self.finished_signal.emit()
                    return None, None

                if keep_img:
                    cv2.imwrite(f"{magnitude}_{str(piezo_pos)}_um_focus.png", original_image.astype(np.uint16))

                std = self.__calc_image_metrics__(original_image)

                # save standard deviation
                piezo_values.append([float(piezo_pos), float(std)])

                # finding maximum value
                # if  not break_flag and
                if i >= 3 and (piezo_values[i - 3][1] > piezo_values[i - 2][1]) \
                        and (piezo_values[i - 3][1] > piezo_values[i - 1][1]) \
                        and (piezo_values[i - 3][1] > piezo_values[i][1]):
                    start_pos = max((piezo_values[i - 3][0]) - (step_size / 2), 0)
                    optimal_piezo = piezo_values[i - 3][0]
                    optimal_std = piezo_values[i - 3][1]
                    # break_flag = True
                    break
                elif piezo_pos == end_pos:
                    start_pos = max((piezo_values[i][0]) - (step_size / 2), 0)
                    optimal_piezo = piezo_values[i][0]
                    optimal_std = piezo_values[i - 3][1]
                    break

        self.__set_signal_block_status__(False)

        # set piezo to optimal position
        status = self.imaging_piezo.set_piezo_pos(optimal_piezo)
        if not status:
            # error signal has been emitted by piezo already
            if standalone_call:
                self.device_available_signal.emit(self.name, True)
                self.finished_signal.emit()
            return None, None

        if standalone_call:
            self.device_available_signal.emit(self.name, True)
            self.finished_signal.emit()
        return optimal_piezo, optimal_std

    def convert_px_to_abs_coordinates(self, px_coord_array: np.ndarray, reference_abs_pos, reference_px_coords, mm_px_scaling):
        if px_coord_array.shape[1] != 2 or len(px_coord_array) < 1:
            info = "List of coordinates must be of form (n, 2) [y, x]"
            self.error_signal.emit((self.name, False, info))
            return None

        # rotation does not need to be considered, because coordinates will be transformed later
        abs_pos_array = np.zeros(px_coord_array.shape)
        # pixel coordinates = [y, x] --> absolute coordinates = [x, y]
        for i, px_coords in enumerate(px_coord_array):
            # y coords
            if px_coords[0] > reference_px_coords[0]:
                abs_pos_array[i, 1] = round(reference_abs_pos[1] - (px_coords[0] - reference_px_coords[0]) * mm_px_scaling, 4)
            else:
                abs_pos_array[i, 1] = round(reference_abs_pos[1] + (reference_px_coords[0] - px_coords[0]) * mm_px_scaling, 4)

            # x coords
            if px_coords[1] > reference_px_coords[1]:
                abs_pos_array[i, 0] = round(reference_abs_pos[0] + (px_coords[1] - reference_px_coords[1]) * mm_px_scaling, 4)
            else:
                abs_pos_array[i, 0] = round(reference_abs_pos[0] - (reference_px_coords[1] - px_coords[1]) * mm_px_scaling, 4)

        return abs_pos_array

    def get_mm_px_scaling(self, marker_centers_px_pos):
        median_marker_px_distance = np.median(np.asarray([
            math.hypot(marker_centers_px_pos[1, 1] - marker_centers_px_pos[0, 1],
                       marker_centers_px_pos[1, 0] - marker_centers_px_pos[0, 0]),
            math.hypot(marker_centers_px_pos[3, 1] - marker_centers_px_pos[2, 1],
                       marker_centers_px_pos[3, 0] - marker_centers_px_pos[2, 0]),
            math.hypot(marker_centers_px_pos[2, 1] - marker_centers_px_pos[0, 1],
                       marker_centers_px_pos[2, 0] - marker_centers_px_pos[0, 0]),
            math.hypot(marker_centers_px_pos[3, 1] - marker_centers_px_pos[1, 1],
                       marker_centers_px_pos[3, 0] - marker_centers_px_pos[1, 0])]))

        mm_px_scaling = (self.sample.mf_size / median_marker_px_distance) / 1000
        return mm_px_scaling

    def save_triangulation_parameters(self, point_idx):
        # point_idx: 0 = bottom-left, 1 = bottom-right, 2 = top-left
        self.device_available_signal.emit(self.name, False)

        # "manual" focusing
        original_image = self.camera.capture_image()
        if original_image is None:
            self.device_available_signal.emit(self.name, True)
            self.finished_signal.emit()
            return
        optimal_std = self.__calc_image_metrics__(original_image)
        optimal_piezo_pos = self.imaging_piezo.get_piezo_pos(standalone_call=False)

        """# autofocus to determine z-coordiante
        optimal_piezo_pos, optimal_std = self.autofocus(standalone_call=False)
        if optimal_piezo_pos is None or optimal_std is None:
            # error caught by device
            self.device_available_signal.emit(self.name, True)
            self.finished_signal.emit()
            print("Autofocus failed")
            return
        
        # capture imagewith optimal properties again
        original_image = self.camera.capture_image()
        if original_image is None:
            self.device_available_signal.emit(self.name, True)
            self.finished_signal.emit()
            return
        print("Autofocus done successfully")"""

        # acquiring x- and y-coordinates from stage
        stage_pos_x, stage_pos_y = self.xy_stage.get_stage_position(standalone_call=False)
        if stage_pos_x is None or stage_pos_y is None:
            # error caught by device
            self.device_available_signal.emit(self.name, True)
            self.finished_signal.emit()
            return

        # // TODO: apply for rotated images
        # Heavy rework! --> needs function structure, e.g. for access to QDs

        cropped_img = square_image(original_image)
        marker_centers_px_pos = find_marker_origin_conventional(cropped_img)

        if marker_centers_px_pos is False:
            self.error_signal.emit((self.name, False, f"Marker centers could not be detected!"))
            self.device_available_signal.emit(self.name, True)
            self.finished_signal.emit()
            return

        mm_px_scaling = self.get_mm_px_scaling(marker_centers_px_pos)

        reference_px_coords = np.asarray([cropped_img.shape[0]//2, cropped_img.shape[1]//2])
        reference_abs_pos = np.asarray([stage_pos_x, stage_pos_y])
        # rotation does not need to be considered, because coordinates will be transformed later
        abs_pos_array = self.convert_px_to_abs_coordinates(marker_centers_px_pos,
                                                           reference_abs_pos=reference_abs_pos,
                                                           reference_px_coords=reference_px_coords,
                                                           mm_px_scaling=mm_px_scaling)

        # taking the relevant markers (bl: point_idx = 0, marker_idx = 2; br: point_idx = 1, marker_idx = 3; tl: point_idx = 2, marker_idx = 0
        if point_idx == 2:
            triang_coordinates = abs_pos_array[0]
        elif point_idx == 1:
            triang_coordinates = abs_pos_array[3]
        else:
            triang_coordinates = abs_pos_array[2]

        self.triangpos_parameter_signal.emit((point_idx, triang_coordinates, optimal_piezo_pos, optimal_std, stage_pos_x, stage_pos_y))
        self.device_available_signal.emit(self.name, True)
        self.finished_signal.emit()
        return

    def calc_triangulation_parameters(self):
        # check if triangulation has been performed, make sure all 3 points are non-zero
        if self.sample.triangulation_points[0].sum() == 0 or \
                self.sample.triangulation_points[1].sum() == 0 or \
                self.sample.triangulation_points[2].sum() == 0:
            self.error_signal.emit((self.name, False, "Sample triangulation has not been completed yet!"))
            return False

        ### calculations (µm for sample, mm for stage)
        # rotation matrix in case sample is not perpendicular
        # first index: 0 = bottom-left, 1 = bottom-right, 2 = top-left

        # determine length of the sample via x ([0]) and y coordinates ([1])
        self.sample.length_x = math.hypot(
            self.sample.triangulation_points[1, 0] - self.sample.triangulation_points[0, 0],
            self.sample.triangulation_points[1, 1] - self.sample.triangulation_points[0, 1]) * 1000
        self.sample.length_y = math.hypot(
            self.sample.triangulation_points[2, 0] - self.sample.triangulation_points[0, 0],
            self.sample.triangulation_points[2, 1] - self.sample.triangulation_points[0, 1]) * 1000
        sin_theta = (self.sample.triangulation_points[1, 1] - self.sample.triangulation_points[0, 1]) / (self.sample.length_x / 1000)
        cos_theta = (self.sample.triangulation_points[1, 0] - self.sample.triangulation_points[0, 0]) / (self.sample.length_x / 1000)

        self.sample.rotation_matrix = np.asarray([[cos_theta, -sin_theta], [sin_theta, cos_theta]])

        self.sample.step_deviation_factor_x = self.sample.length_x / self.sample.actual_sample_size
        self.sample.step_deviation_factor_y = self.sample.length_y / self.sample.actual_sample_size

        # calculation of the plane of the piezo position (to determine optimal piezo focus)
        focus_points = Points([self.sample.triangulation_points[0, :3],
                               self.sample.triangulation_points[1, :3],
                               self.sample.triangulation_points[2, :3]])

        focus_plane = Plane.best_fit(focus_points)

        # cartesian focus plane parameters
        a_x, b_y, c_z, d = focus_plane.cartesian()
        self.sample.piezo_pos_parameters = np.asarray([a_x, b_y, c_z, d])

        # calculation regarding to the std of histograms --> not necessary, worked in tests just fine

        # save position values & angles in txt
        with open('triangulation_params.txt', 'w') as f:
            f.write("Sample properties: \n")
            f.write(f"Width (x-axis): {self.sample.length_x},\t height (y-axis): {self.sample.length_y}\n")
            f.write(f"Focus properties: a_x: {a_x}, b_y: {b_y}, c_z: {c_z}, d: {d}\n")
            f.write("Triangulation position (X, Y) \t values\n")
            for triang_point in self.sample.triangulation_points:
                f.write(f"{triang_point[0]}, {triang_point[1]}\n")
                f.write(f"\tFocus (stage position): {triang_point[4]}, {triang_point[5]}\n")
                f.write(f"\tFocus (piezo position): {triang_point[2]}\n")
                f.write(f"\tStandard deviation of histograms (piezo metric): {triang_point[3]}\n")
                f.write("\n")

        return True

    def __set_signal_block_status__(self, state):
        self.imaging_piezo.blockSignals(state)
        self.camera.blockSignals(state)
        self.xy_stage.blockSignals(state)

    def __finish_action__(self):
        self.device_available_signal.emit(self.name, True)
        self.finished_signal.emit()

    def scanning_sample(self, save_positions_in_file=False):
        self.device_available_signal.emit(self.name, False)

        if not self.calc_triangulation_parameters():
            # error signal has been emitted by function
            self.device_available_signal.emit(self.name, True)
            self.finished_signal.emit()
            return False

        processed_mf_counter = 0

        """# set LED states
        if not self.lin_stage.set_led_state("red", True):
            # error signal already sent by michelson linear stage
            self.device_available_signal.emit(self.name, True)
            return False
        if not self.lin_stage.set_led_state("blue", True):
            # error signal already sent by michelson linear stage
            self.device_available_signal.emit(self.name, True)
            return False"""

        os.makedirs(self.image_path, exist_ok=True)
        os.chdir(self.image_path)
        if save_positions_in_file:
            with open('mf_positions.txt', 'w') as f:
                f.write("Sample properties: \n")
                f.write(f"Cells: {self.sample.n_cells}; Marker fields per cell: {self.sample.n_marker_fields}\n")
                f.write("Cell_ID, MF_ID; Cell no. x, cell no. y, mf no. x, mf no. y \n")
                f.write("\tmetrics: (x, y, z, std):\n")

        # scanning the sample
        # moving to starting point
        if not self.xy_stage.move_stage_absolute(dest_x=self.sample.triangulation_points[0, 4], dest_y=self.sample.triangulation_points[0, 5]):
            # error signal already sent by stage
            self.device_available_signal.emit(self.name, True)
            self.finished_signal.emit()
            return False

        # setting piezo position:
        self.imaging_piezo.set_piezo_pos(self.sample.triangulation_points[0, 2])

        max_std_deviation = 2  # value that the focus might deviate from the plane before autofocus is required

        self.__set_signal_block_status__(True)  # making sure that not all camera images are displayed (would keep GUI too busy)

        for cell_y in range(int(np.sqrt(self.sample.n_cells_y))):
            # follow cells in S-form to minimize required travel distance
            if cell_y % 2 == 0:
                cell_order_x = range(int(self.sample.n_cells_x), 1)
                mf_order_x = range(0, int(np.sqrt(self.sample.n_marker_fields)), 1)
                direction_x = 1
            else:
                cell_order_x = range(int(self.sample.n_cells_x) - 1, -1, -1)
                mf_order_x = range(int(self.sample.n_marker_fields_x) - 1, -1, -1)
                direction_x = -1

            for n_cx, cell_x in enumerate(cell_order_x):
                # follow marker fields in S-form to minimize required travel distance
                for n_mfx, mf_x in enumerate(mf_order_x):
                    if mf_x % 2 == 0:
                        # going from left to right
                        mf_order_y = range(0, int(self.sample.n_marker_fields_y)-1, 1)
                        direction_factor_mf_y = 1 * direction_x
                    else:
                        mf_order_y = range(int(self.sample.n_marker_fields_y) - 1, -1, -1)
                        direction_factor_mf_y = -1 * direction_x

                    # reset marker field rows when going backwards on sample
                    if cell_y % 2 != 0:
                        mf_order_y = reversed(mf_order_y)

                    for n_mfy, mf_y in enumerate(mf_order_y):
                        # current cell and MF identifier
                        cell_id = self.sample.oriented_cells_layout[cell_y, cell_x]
                        mf_id = self.sample.cell[cell_id].oriented_mfs_layout[mf_y, mf_x]

                        self.progress_signal.emit(processed_mf_counter)
                        processed_mf_counter += 1

                        # get current (position) values
                        curr_x_pos, curr_y_pos = self.xy_stage.get_stage_position()
                        if curr_x_pos is None or curr_y_pos is None:
                            # error signal already sent by stage
                            self.__set_signal_block_status__(False)
                            self.device_available_signal.emit(self.name, True)
                            self.finished_signal.emit()
                            return False

                        # Focusing by using cartesian form of plane and calculating point on plane
                        a_x, b_y, c_z, d = self.sample.piezo_pos_parameters
                        piezo_pos = -1*(a_x*curr_x_pos + b_y*curr_y_pos + d)/c_z

                        print("Calculated piezo position, while scanning: ", piezo_pos)

                        if not self.imaging_piezo.set_piezo_pos(piezo_pos):
                            # error signal already sent by piezo
                            self.__set_signal_block_status__(False)
                            self.device_available_signal.emit(self.name, True)
                            self.finished_signal.emit()
                            return False

                        # capture image
                        original_image = self.camera.capture_image()
                        if original_image is None:
                            info = f"Camera failed at Cell: {cell_id}, Marker Field: {mf_id}. Continues with next image"
                            self.error_signal.emit((self.name, True, info))
                            continue

                        std = self.__calc_image_metrics__(original_image)

                        # check if std in range
                        """predicted_image_std = curr_x_pos * self.sample.focus_std_parameters[0] + \
                                              curr_y_pos * self.sample.focus_std_parameters[1] + self.sample.focus_std_parameters[2]

                        # // optimize autofocus around current position
                        # // optimize determination of autofocus value
                        # 2 here is just a rough estimate of the std of the std
                        if not ((predicted_image_std + max_std_deviation) > std > (predicted_image_std - max_std_deviation)):
                            print(f"Autofocus during scan would have applied, because {std} outside of +-2 {predicted_image_std}")

                            # deviation too broad, refocus
                            print("Predicted std should have been: ", predicted_image_std)
                            print("Current std: ", std)
                            keep_img = False
                            piezo_pos, optimal_std = self.autofocus(standalone_call=False, keep_img=keep_img)
                            if piezo_pos is None or optimal_std is None:
                                # error signal already sent by devices
                                self.__set_signal_block_status__(False)
                                self.device_available_signal.emit(self.name, True)
                                self.finished_signal.emit()
                                return False

                            if keep_img:
                                # images are kept separately; setting the path back to where it was before
                                os.chdir(self.image_path)
                            print("Optimal std after focus: ", optimal_std)
                            print(f"Moving Piezo to {piezo_pos}")
                            if not self.imaging_piezo.set_piezo_pos(piezo_pos):
                                # error signal already sent by piezo
                                self.__set_signal_block_status__(False)
                                self.device_available_signal.emit(self.name, True)
                                self.finished_signal.emit()
                                return False

                            # capture new image with optimized focus
                            original_image = self.camera.capture_image()
                            if original_image is None:
                                # error signal already sent by devices
                                self.__set_signal_block_status__(False)
                                self.device_available_signal.emit(self.name, True)
                                self.finished_signal.emit()
                                return False"""

                        # save all data in respective mf_field etc.
                        self.sample.cell[cell_id].marker_field[mf_id].camera_pos[:] = curr_x_pos, curr_y_pos
                        self.sample.cell[cell_id].marker_field[mf_id].focus_piezo_pos = piezo_pos

                        print(f"Current Cell: {cell_id}, current mf: {mf_id}")

                        # save positions:
                        # done here and not at the end to safe results in case of break
                        if save_positions_in_file:
                            try:
                                with open('mf_positions.txt', 'a') as f:
                                    f.write(f"\n{cell_id}, {mf_id};\t{cell_y}, {cell_x},\t{mf_x}, {mf_y}\n")
                                    f.write(f"\t{curr_x_pos}, {curr_y_pos}\t-\t{piezo_pos}, {std}\n")
                            except FileNotFoundError:
                                info = f"File to store marker positions not found."
                                self.error_signal.emit((self.name, True, info))
                                self.__set_signal_block_status__(False)
                                self.device_available_signal.emit(self.name, True)
                                self.finished_signal.emit()
                                return False

                        # save image
                        cv2.imwrite(f"{self.sample.name}_{cell_id}-{mf_id.lower()}.png", original_image.astype(np.uint16))

                        # stage needs to be moved only n-1 times, but action needs to be executed n times
                        if n_mfy < int(np.sqrt(self.sample.n_marker_fields)) - 1:
                            # move stage to next marker field in y direction
                            rel_mf_distance_y = direction_factor_mf_y * self.sample.mf_size * self.sample.step_deviation_factor_y
                            if not self.xy_stage.move_stage_relative(0, rel_mf_distance_y, self.sample.rotation_matrix):
                                # error signal already sent by stage
                                self.__set_signal_block_status__(False)
                                self.device_available_signal.emit(self.name, True)
                                self.finished_signal.emit()
                                return False

                    # move stage to next marker field in x direction
                    if n_mfx < int(np.sqrt(self.sample.n_marker_fields)) - 1:
                        rel_mf_distance_x = direction_x * self.sample.mf_size * self.sample.step_deviation_factor_x
                        if not self.xy_stage.move_stage_relative(rel_mf_distance_x, 0, self.sample.rotation_matrix):
                            # error signal already sent by stage
                            self.__set_signal_block_status__(False)
                            self.device_available_signal.emit(self.name, True)
                            self.finished_signal.emit()
                            return False

                # move stage to next cell in x direction
                if n_cx < int(np.sqrt(self.sample.n_cells)) - 1:
                    rel_cell_distance_x = direction_x * (self.sample.inter_cell_distance + self.sample.mf_size) * self.sample.step_deviation_factor_x
                    if not self.xy_stage.move_stage_relative(rel_cell_distance_x, 0, self.sample.rotation_matrix):
                        # error signal already sent by stage
                        self.__set_signal_block_status__(False)
                        self.device_available_signal.emit(self.name, True)
                        self.finished_signal.emit()
                        return False

            # move stage to next cell in y direction (upwards)

            rel_cell_distance_y = (np.sqrt(self.sample.n_marker_fields) * self.sample.mf_size + self.sample.inter_cell_distance) * self.sample.step_deviation_factor_y
            if not self.xy_stage.move_stage_relative(0, rel_cell_distance_y, self.sample.rotation_matrix):
                # error signal already sent by stage
                self.__set_signal_block_status__(False)
                self.device_available_signal.emit(self.name, True)
                self.finished_signal.emit()
                return False

        self.__set_signal_block_status__(False)
        self.device_available_signal.emit(self.name, True)
        self.finished_signal.emit()
        return True


    def capture_image_series(self, move_stage=False, n_images=30):
        self.__set_signal_block_status__(True)
        self.device_available_signal.emit(self.name, False)

        if move_stage:
            stage_step_samples = np.arange(-1, 1.1, 0.1)  # stage steps between -1 and 1 µm
            base_x_pos, base_y_pos = self.xy_stage.get_stage_position()
            if base_x_pos is None or base_y_pos is None:
                # error signal already sent by stage
                self.__set_signal_block_status__(False)
                self.__finish_action__()
                return False

        dir_path = os.path.join(str(Path(__file__).parent.parent), "resources", "MPE_analysis_",
                                datetime.now().strftime("%Y%m%d_%H-%M-%S"))
        os.makedirs(dir_path, exist_ok=True)
        cwd = Path.cwd()
        os.chdir(dir_path)

        for step in range(n_images):
            original_image = self.camera.capture_image()
            if original_image is None:
                info = f"Camera failed"
                os.chdir(cwd)
                self.error_signal.emit((self.name, False, info))
                self.__set_signal_block_status__(False)
                self.__finish_action__()
                return False

            # save image
            # // ToDO: include device parameters once everything is hooked up
            cv2.imwrite(f"{step}.png", original_image.astype(np.uint16))

            if move_stage:
                curr_x_pos, curr_y_pos = self.xy_stage.get_stage_position()
                if curr_x_pos is None or curr_y_pos is None:
                    # error signal already sent by stage
                    self.__set_signal_block_status__(False)
                    self.__finish_action__()
                    return False

                sample_indices = np.random.randint(0, len(stage_step_samples), 2)

                x_dist = round(stage_step_samples[sample_indices[0]], 4)  # values between 0 and 0.9[µm]
                y_dist = round(stage_step_samples[sample_indices[1]], 4)

                if round(abs(base_x_pos - (curr_x_pos + x_dist / 1000)), 4) > 0.001:  # making sure that stage does not leave image
                    x_dist = 0
                if round(abs(base_y_pos - (curr_y_pos + y_dist / 1000)), 4) > 0.001:  # making sure that stage does not leave image
                    y_dist = 0

                # print(f"\tMoving stage by: {x_dist}, {y_dist} [µm]")

                if not self.xy_stage.move_stage_relative(x_dist, y_dist):
                    # error signal already sent by stage
                    os.chdir(cwd)
                    self.__set_signal_block_status__(False)
                    self.__finish_action__()
                    return False

        os.chdir(cwd)
        if not self.xy_stage.move_stage_absolute(base_x_pos, base_y_pos):
            # error signal already sent by stage
            self.__set_signal_block_status__(False)
            self.__finish_action__()
            return False

        self.__set_signal_block_status__(False)
        self.__finish_action__()
        return True

    def scan_designated_cells(self, cell_list: list, standalone_call=True):
        if self.sample is None:
            self.error_signal.emit((self.name, False, "Cannot go to cell before sample has been created!"))
            return False

        if standalone_call:
            self.device_available_signal.emit(False)
        for cell_name in cell_list:
            cell_name = cell_name.upper()
            if cell_name not in self.sample.basic_cells_layout:
                self.error_signal.emit((self.name, False, f"Cell {cell_name} does not exist"))
                if standalone_call:
                    self.__finish_action__()
                return False

            # // ToDO: Continue here

        if standalone_call:
            self.__finish_action__()

    def qd_spectroscopy(self, n_splits_fss, cell_name, marker_field_name, use_only_search_pos=False, qd_detection_thresh=1.26, standalone_call=True):
        print("In MarkerScanner QD spectroscopy")

        # get current image
        image = self.camera.capture_image()
        if image is None:
            # error signal has been emitted by camera already
            self.__set_signal_block_status__(False)
            if standalone_call:
                self.__finish_action__()
            return False

        # find markers
        image = square_image(image)
        marker_centers_px_pos = find_marker_origin_conventional(image)

        if marker_centers_px_pos is False:
            self.error_signal.emit((self.name, False, f"Marker centers could not be detected!"))
            self.__finish_action__()
            return False
        marker_min_y = marker_centers_px_pos[:, 0].min()
        marker_max_y = marker_centers_px_pos[:, 0].max()
        marker_min_x = marker_centers_px_pos[:, 1].min()
        marker_max_x = marker_centers_px_pos[:, 1].max()

        stage_pos_x, stage_pos_y = self.xy_stage.get_stage_position(standalone_call=False)
        if stage_pos_x is None or stage_pos_y is None:
            self.__finish_action__()
            return False

        # find QDs, parameters found via trial and error (function returns peak as Y, X)
        peak_list, sigmas = get_qd_positions_from_image(image, area=17, thresh=(
                    np.mean(image) * qd_detection_thresh))

        # only consider QDs within Marker field
        peak_list = peak_list[np.where(np.logical_and(np.logical_and(marker_min_y < peak_list[:, 0],
                                                                     peak_list[:, 0] < marker_max_y),
                                                      np.logical_and(marker_min_x < peak_list[:, 1],
                                                                     peak_list[:, 1] < marker_max_x)))]

        # // ToDo: reduce number of QDs as in semi_automatic mode (via Flag)
        if use_only_search_pos:
            print("Reuse stuff from QD Detector")

        # get stage coordinates from the image
        mm_px_scaling = self.get_mm_px_scaling(marker_centers_px_pos)  # // Needs to use x / y

        reference_px_coords = np.asarray([image.shape[0] // 2, image.shape[1] // 2])
        reference_abs_pos = np.asarray([stage_pos_x, stage_pos_y])
        # rotation does not need to be considered, because coordinates will be transformed later
        abs_marker_center_pos_array = self.convert_px_to_abs_coordinates(marker_centers_px_pos,
                                                           reference_abs_pos=reference_abs_pos,
                                                           reference_px_coords=reference_px_coords,
                                                           mm_px_scaling=mm_px_scaling)

        abs_qd_pos_array = self.convert_px_to_abs_coordinates(peak_list,
                                                              reference_abs_pos=reference_abs_pos,
                                                              reference_px_coords=reference_px_coords,
                                                              mm_px_scaling=mm_px_scaling)

        if abs_marker_center_pos_array is None or abs_qd_pos_array is None:
            self.error_signal.emit((self.name, False, f"Stage positions could not be inferred!"))
            self.__finish_action__()
            return False

        # turn off IR LED for spectroscopy
        if not self.lin_stage.set_led_state("red", False):
            self.__finish_action__()
            return False

        # // ToDo: Move thr_elliptec linear stage in new setup
        print("Linear stage has not been resolved")

        target_dir = os.path.join(self.image_path, self.sample.name, cell_name, marker_field_name)
        Path(target_dir).mkdir(parents=True, exist_ok=True)
        os.chdir(target_dir)

        # create output structure for
        excitons_wl_array = np.zeros((len(abs_qd_pos_array), ))
        exciton_txt_filename = f"exciton_wavelength_{self.sample.name}_{cell_name}_{marker_field_name}.txt"
        with open(exciton_txt_filename, 'w') as f:
            f.write(f"Sample:\t{self.sample.name}\n")
            f.write(f"Cell:\t{cell_name}\n")
            f.write(f"Marker Field:\t{marker_field_name}\n")
            f.write(f"QD:\texciton wavelength\n\n")

        # iterate over all qd positions
        for n_qd, qd_pos in enumerate(abs_qd_pos_array):
            if not self.xy_stage.move_stage_absolute(qd_pos[0], qd_pos[1]):
                self.__finish_action__()
                return False

            # #Pol(deg)_0_180 = starting angle_final angle for FSS
            filename_noextension = f"{self.sample.name}_{cell_name}_{marker_field_name}_{n_qd}#Pol(deg)_0_360 "
            filename = filename_noextension + '.spe'

            # Copied from measurement-py from Kai, I suggest changing this for "acquire" (?) to keep the spe files
            spectrometer_status, spectrum = self.spectrometer.acquire_spectrum(standalone_call=False)
            if not spectrometer_status:
                self.__finish_action__()
                return False
            wl = self.wl_regressor(torch.FloatTensor(standardize(np.asarray(spectrum[1])))).detach().cpu().item()
            idx = int(len(spectrum[0]) * wl)
            # not on edge to enable average later on
            exciton_wavelength = max(min(idx, len(spectrum[0]) - 2), 1)

            # publish results
            self.exciton_wl_signal.emit(int(exciton_wavelength))
            excitons_wl_array[n_qd] = exciton_wavelength
            with open(exciton_txt_filename, 'a') as f:
                f.write(f"{n_qd}\t{exciton_wavelength}\n")

            # spectrometer settings for FSS measurement according to meas1d.py from Ailton
            spectrometer_status = self.spectrometer.stop(standalonecall=False)
            if spectrometer_status is False:
                self.__finish_action__()
                return False

            spectrometer_status = self.spectrometer.frames_to_save(frames=n_splits_fss, standalonecall=False)
            if spectrometer_status is False:
                self.__finish_action__()
                return False

            spectrometer_status = self.spectrometer.trigger_settings("readout", standalone_call=False)
            if spectrometer_status is False:
                self.__finish_action__()
                return False

            spectrometer_status = self.spectrometer.change_folder(folder_name=str(target_dir), standalonecall=False)
            if spectrometer_status is False:
                self.__finish_action__()
                return False

            spectrometer_status = self.spectrometer.change_filename(file_name=filename_noextension, standalonecall=False)
            if spectrometer_status is False:
                self.__finish_action__()
                return False

            spectrometer_status = self.spectrometer.acquire_file(standalone_call=False)
            if spectrometer_status is False:
                self.__finish_action__()
                return False

            # going into Fine Structure Splitting (FSS) (Kai)
            for split_no, j in enumerate(range(n_splits_fss), 1):
                fss_pos = (j * 180 / n_splits_fss)
                wp_status = self.waveplate.set_position(fss_pos)
                if wp_status is False:
                    self.__finish_action__()
                    return False

                if split_no == len(n_splits_fss):  # send falling edge to signal end of spectroscopy measurement
                    lin_stage_status = self.lin_stage.trigger_spectrometer(False, standalonecall=False)
                else:
                    lin_stage_status = self.lin_stage.trigger_spectrometer(True, standalonecall=False)

                if lin_stage_status is False:
                    self.__finish_action__()
                    return False

        #might need to add an extra trigger signal to finish the measurement.

        # // ToDo: move thr_elliptec linear stage out of the visual way

        # Spectroscopy + FSS for all QDs in image has been done
        status = self.xy_stage.move_stage_absolute(stage_pos_x, stage_pos_y)  # return to center of MF for scanning routine
        if not status:
            self.error_signal.emit((self.name, True, "Stage could not move back to MF center, will screw up scanning routine"))

        _led_status = self.lin_stage.set_led_state("red", True)
        self.__finish_action__()
        return status

