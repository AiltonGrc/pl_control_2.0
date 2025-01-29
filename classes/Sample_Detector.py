import os
import numpy as np
import string
from glob import glob
import skimage.io
import cv2
from pathlib import Path, PurePath
from source.sample_structure.sample import Sample
from source.utils.marker_detection import find_marker_origin_conventional, square_image
from source.utils.qd_detection import get_qd_positions_from_image
from PyQt5.QtCore import QObject, pyqtSignal
import shutil
import math
import csv
from tqdm import tqdm
import time
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from source.utils.qd_detector import get_image_qd


# Define the Gaussian function
def gauss(x, A, mu, sigma, offset):
    return offset + A * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))

# Definition of log-normal // --> should be moved to utils
def log_normal(x, A, mu, sigma):
    return A / (x * sigma * np.sqrt(2 * np.pi)) * np.exp(- ((np.log(x) - mu)**2 / (2 * sigma**2) ))

class SampleDetector(QObject):
    progress_signal = pyqtSignal(int)
    error_signal = pyqtSignal(tuple)
    found_files_signal = pyqtSignal(int)
    device_available_signal = pyqtSignal(str, bool)
    finished_signal = pyqtSignal()

    def __init__(self, image_path, output_path):
        super().__init__()
        self.image = None
        self.img_height = 0
        self.img_width = 0
        self.image_path = image_path
        self.output_path = output_path
        self.config_file_path = Path(__file__).parent.parent / "resources" / "settings" / "semi_manual_qd_detection.conf"
        self.file_list = None
        self.sample = None
        self.color = (2 ** 16, 2 ** 16, 2 ** 16)
        self.mf_type_threshold = 20000

        self.name = "SampleProcessor"

        # limiting the detected QD output by pre-determining coordiantes (in µm)
        self.qd_search_positions = dict()
        self.qd_analysis_parameter = dict()
        self.selected_qd_dict = dict()

    def add_sample(self, sample: Sample):
        self.device_available_signal.emit(self.name, False)
        print("In adding sample to Detector")
        self.sample = sample
        self.device_available_signal.emit(self.name, True)
        self.finished_signal.emit()
        return

    def read_image(self, filename, read_fits=False):
        # os.chdir(self.image_path)
        if read_fits:
            try:
                img = skimage.io.imread(filename, as_gray=True, plugin='fits')
            except OSError:
                info = f"Image {filename} corrupted"
                self.error_signal.emit((self.name, True, info))
                return False
            except FileNotFoundError:
                info = f"Image {filename} not found. Please check path"
                self.error_signal.emit((self.name, True, info))
                return False

            # skimage.io flips the image along its horizontal axis --> needs to be revereted here
            img = np.flipud(img).astype(np.uint16)
        else:
            img = cv2.imread(filename, -1)
            if img is None:
                info = f"Image {filename} corrupted"
                self.error_signal.emit((self.name, True, info))
                return False

        self.img_height = img.shape[0]
        self.img_width = img.shape[1]

        # checking image properties
        min_size = 600
        if self.img_width < min_size or self.img_height < min_size:
            info = f"Image of invalid format. Needs to be at least {min_size} px."
            self.error_signal.emit((self.name, True, info))
            return False
        if img.var() <= 0:
            info = f"Image of invalid format. Needs to be at least {min_size} px."
            self.error_signal.emit((self.name, True, info))
            return False

        if img.mean() < self.mf_type_threshold:
            self.color = (0, 0, 0)

        self.image = img

        return True

    def read_config_file(self):
        current_path = os.getcwd()

        config_directory = self.config_file_path.parent
        config_file_name = self.config_file_path.name

        os.chdir(config_directory)
        qd_analysis_dict = dict()
        qd_name_list = list(string.ascii_lowercase)[:11]  # limited to 10 points

        with open(config_file_name, "r") as f:
            for line in f:
                # exclude comment lines
                if line.strip().startswith("#") or line == '' or line == "\n":
                    continue
                # exclude inline comments
                comment_idx = line.find('#')
                if comment_idx != -1:
                    line = line[:comment_idx]
                kv = line.split("=")

                try:
                    config_key = kv[0].strip()
                    config_value = kv[1].strip()
                except ValueError:
                    self.gui_error_signal(("GUI", False, "Reading in the config-file failed due to wrong input."))
                    return False

                # read lines and store values accordingly
                if config_key in qd_name_list:
                    self.qd_search_positions[config_key] = np.asarray(list(map(int, config_value.split(","))))
                # check numbers
                elif '"' in config_value:  # string
                    qd_analysis_dict[config_key] = config_value.replace('"', '')
                elif any(char.isdigit() for char in config_value):
                    if "." in config_value:
                        qd_analysis_dict[config_key] = float(config_value)
                    else:  # must be int
                        qd_analysis_dict[config_key] = int(config_value)
                # check boolean
                elif config_value.lower() == "true":
                    qd_analysis_dict[config_key] = True
                elif config_value.lower() == "false":
                    qd_analysis_dict[config_key] = False

        self.qd_analysis_parameter = qd_analysis_dict

        os.chdir(current_path)
        return True

    def create_file_list(self, path, read_fits=False):
        # used to grant GUI information for the progressbar
        # setting progress bar "busy" while looking for files
        self.found_files_signal.emit(0)

        if read_fits:
            self.file_list = sorted(glob(os.path.join(path, '**', '*.fits'), recursive=True))
        else:
            self.file_list = sorted(glob(os.path.join(path, '**', '*.png'), recursive=True))

        n_files = len(self.file_list)

        if n_files == 0:
            info = "No images to analyze. Check image path.or start scanning first!"
            self.error_signal.emit(("Detector", False, info))
            # set progress bar as available again
            self.found_files_signal.emit(100)
            return False

        self.found_files_signal.emit(n_files)

        return True

    def save_information_to_file(self, filename, object_name, coordinates: np.ndarray):
        path = os.path.join(os.getcwd(), filename)
        if not os.path.exists(path):
            with open(filename, 'w') as f:
                f.write(f"Sample: {self.sample.name} \n")
                f.write("Information are given as relative distances to the last origin position in the following form:\n")
                f.write("Object\t-\t x-coordinate [um], y-coordinate [um]\n")
                f.write("Reference points are stored with their respective absolute coordinates.\n")

        with open(filename, 'a') as f:
            for coord in coordinates:
                if "origin" in object_name:
                    f.write("\n")
                    f.write(f"Reference origin: {object_name}\n")
                    f.write(f"\tAbsolute stage position: {coord[0]} um, {coord[1]} um\n")
                else:
                    f.write(f"{object_name}\t-\t {coord[0]} um, {coord[1]} um\n")

    def __convert_px_to_abs_coordinates__(self, px_coord_array: np.ndarray, reference_abs_pos, reference_px_coords, um_px_scale_y, um_px_scale_x):
        if px_coord_array.shape[1] != 2 or len(px_coord_array) < 1:
            info = "List of coordinates must be of form (n, 2) [y, x]"
            self.error_signal.emit(self.name, False, info)
            return None

        # rotation does not need to be considered, because coordinates will be transformed later
        abs_pos_array = np.zeros(px_coord_array.shape)
        # pixel coordinates = [y, x] --> absolute coordinates = [x, y]
        for i, px_coords in enumerate(px_coord_array):
            # y coords
            if px_coords[0] > reference_px_coords[0]:
                abs_pos_array[i, 1] = reference_abs_pos[1] - (px_coords[0] - reference_px_coords[0]) * um_px_scale_y
            else:
                abs_pos_array[i, 1] = reference_abs_pos[1] + (reference_px_coords[0] - px_coords[0]) * um_px_scale_y

            # x coords
            if px_coords[1] > reference_px_coords[1]:
                abs_pos_array[i, 0] = reference_abs_pos[0] + (px_coords[1] - reference_px_coords[1]) * um_px_scale_x
            else:
                abs_pos_array[i, 0] = reference_abs_pos[0] - (reference_px_coords[1] - px_coords[1]) * um_px_scale_x

        return abs_pos_array

    def __find_in_layout__(self, obj_id: str, layout_array: np.ndarray):
        key_location_y, key_location_x = np.where(layout_array == obj_id)

        if key_location_y.shape[0] == 0:
            # in case that ID cannot be found!
            print("Object ID was not found in layout")
            return np.asarray([0, 0])

        return np.asarray([key_location_x[0], key_location_y[0]])

    def detect_dots_from_manual_images(self, marker_type_dark=True):
        self.device_available_signal.emit(self.name, False)

        # read relevanbt information from config file
        if not self.read_config_file():
            # error signal already sent
            self.device_available_signal.emit(self.name, True)
            self.finished_signal.emit()
            return

        # define frequently used or long-named variables
        sample_name = self.qd_analysis_parameter["sample_name"]
        marker_center_origin_idx = self.qd_analysis_parameter["marker_center_origin_idx"]
        mf_size = self.qd_analysis_parameter["marker_field_size"]
        cell_to_cell_distance = self.qd_analysis_parameter["inter_cell_distance"] + np.sqrt(self.qd_analysis_parameter["n_marker_fields"]) * mf_size
        cell_origin_offset = np.asarray([self.qd_analysis_parameter["cell_origin_offset"]]*2)
        sample_origin_offset = np.asarray([self.qd_analysis_parameter["sample_origin_offset"]]*2)
        keep_images = self.qd_analysis_parameter["keep_images"]
        save_additional_info = self.qd_analysis_parameter["save_additional_info"]
        cell_directory_list = sorted(glob(os.path.join(self.qd_analysis_parameter["search_directory_path"], "*", "")))

        if self.sample is None:
            self.sample = Sample(name=sample_name,
                                 n_cells=self.qd_analysis_parameter["n_cells"],
                                 n_mfs=self.qd_analysis_parameter["n_marker_fields"],
                                 mf_size=mf_size,
                                 inter_cell_distance=self.qd_analysis_parameter["inter_cell_distance"])

        # making sure there actually are directories
        if len(cell_directory_list) == 0:
            self.error_signal.emit((self.name, False, "Could not find any cell-subdirectories in path"))
            self.device_available_signal.emit(self.name, True)
            self.finished_signal.emit()
            return

        for cell_directory in tqdm(cell_directory_list):
            cell_directory_name = str(Path(cell_directory).name)

            # ignore all directories that do not belong to the sample
            if sample_name not in cell_directory_name:
                continue

            mf_directories_list = sorted(glob(os.path.join(cell_directory, '**', '*.fits'), recursive=True))
            if len(mf_directories_list) == 0:
                continue

            cell_id = cell_directory_name.split("_")[-1].strip()
            csv_file_name = f"{sample_name}_{cell_id}.csv"

            os.chdir(cell_directory)
            with open(csv_file_name, 'w', newline='') as csv_f:
                csv_qd_writer = csv.writer(csv_f, delimiter=";")
                if save_additional_info:
                    csv_qd_writer.writerow(
                        ["Sample name", "Cell-ID", "MF-ID", "QD-identifier", "x coordinate", "y coordinate",
                         "x distance to cell origin", "y distance to cell origin", "x distance to sample origin",
                         "y distance to sample origin"])
                else:
                    csv_qd_writer.writerow(
                        ["Sample name", "Cell-ID", "MF-ID", "QD-identifier", "x coordinate", "y coordinate"])

            for i, mf_image_path in enumerate(mf_directories_list):
                self.progress_signal.emit(i)

                filename = mf_image_path.split("\\")[-1]

                # get marker field identification from path:
                mf_id = filename.split("-")[1].split("_")[0].strip().upper()

                read_state = self.read_image(filename, read_fits=True)
                if not read_state:
                    info = f"Marker Field {mf_id} of cell {cell_id} skipped, because read_image failed"
                    self.error_signal.emit((self.name, True, info))
                    # shutil.move(mf_image_path, os.path.join(manual_detection_path, filename)) # //
                    continue

                # image cannot be rotated --> screws up dimensions
                # self.image = np.rot90(self.image, k=rotation_steps, axes=(1, 0)) # rotates image by 90° k times to the right: ^ > v <

                if self.qd_analysis_parameter["flip_image_lr"]:
                    self.image = np.fliplr(self.image).astype(np.uint16)

                # find marker centers
                marker_centers_px_pos = find_marker_origin_conventional(self.image, marker_type_dark=marker_type_dark, show_images=False)  # at 0° rotation: tl = 0, tr = 1, bl = 2, br = 3

                # MF image could not be detected properly, meaning less then 4 corners were detected
                if marker_centers_px_pos is False:
                    info = f"Marker Field {mf_id} of cell {cell_id} skipped, because marker detection failed"
                    # shutil.move(mf_image_path, os.path.join(manual_detection_path, filename)) # //
                    self.error_signal.emit((self.name, True, info))
                    continue

                if keep_images:
                    detected_image = (self.image.copy()/256).astype(np.uint8)
                    # draw lines between marker centers
                    cv2.line(detected_image,
                             (marker_centers_px_pos[0, 1], marker_centers_px_pos[0, 0]),
                             (marker_centers_px_pos[1, 1], marker_centers_px_pos[1, 0]),
                             (256, 256, 256), 1, cv2.LINE_AA)
                    cv2.line(detected_image,
                             (marker_centers_px_pos[0, 1], marker_centers_px_pos[0, 0]),
                             (marker_centers_px_pos[2, 1], marker_centers_px_pos[2, 0]),
                             (256, 256, 256), 1, cv2.LINE_AA)
                    cv2.line(detected_image,
                             (marker_centers_px_pos[3, 1], marker_centers_px_pos[3, 0]),
                             (marker_centers_px_pos[1, 1], marker_centers_px_pos[1, 0]),
                             (256, 256, 256), 1, cv2.LINE_AA)
                    cv2.line(detected_image,
                             (marker_centers_px_pos[3, 1], marker_centers_px_pos[3, 0]),
                             (marker_centers_px_pos[2, 1], marker_centers_px_pos[2, 0]),
                             (256, 256, 256), 1, cv2.LINE_AA)

                # determine the reference marker (origin)
                marker_field_origin_px_pos = marker_centers_px_pos[marker_center_origin_idx]

                # find QDs, parameters found via trial and error (function returns peak as Y, X)
                peak_list, sigmas = get_qd_positions_from_image(self.image, area=17, thresh=(np.mean(self.image)*self.qd_analysis_parameter["qd_detection_threshold"]))

                # determine relative px-position: absolute allows calculation without differentiation of positive & negative values
                qd_rel_px_pos = np.absolute(marker_field_origin_px_pos - peak_list)  # returns [y_image, x_image]

                # determining the scaling factor how [µm/px]
                # Assumption that no axis-specific distortion appears; taking the median to avoid outliers
                # // ToDO: replace with best distance?
                median_marker_px_distance = np.median(np.asarray([
                    math.hypot(marker_centers_px_pos[1, 1] - marker_centers_px_pos[0, 1],
                               marker_centers_px_pos[1, 0] - marker_centers_px_pos[0, 0]),
                    math.hypot(marker_centers_px_pos[3, 1] - marker_centers_px_pos[2, 1],
                               marker_centers_px_pos[3, 0] - marker_centers_px_pos[2, 0]),
                    math.hypot(marker_centers_px_pos[2, 1] - marker_centers_px_pos[0, 1],
                               marker_centers_px_pos[2, 0] - marker_centers_px_pos[0, 0]),
                    math.hypot(marker_centers_px_pos[3, 1] - marker_centers_px_pos[1, 1],
                               marker_centers_px_pos[3, 0] - marker_centers_px_pos[1, 0])]))

                um_px_scaling = mf_size/median_marker_px_distance

                # calculate physical distances of QDs relative to marker field origin [µm]
                # convert px-data [y, x] to physical coordiantes [x, y] --> swap axes
                # if the image is rotated by 90° or 270°, axes need to be swapped again --> cancels out previous step
                if marker_center_origin_idx == 1 or marker_center_origin_idx == 2:
                    # coordinate system [x, y] lines up with image, but image axes [y, x] need to be transformed [x, y]
                    qd_relative_um_distance = qd_rel_px_pos[:, ::-1] * um_px_scaling
                else:
                    qd_relative_um_distance = qd_rel_px_pos * um_px_scaling

                # find mf and cell in layout for calculation of origin coordinates (additional information)
                if self.sample is not None:
                    cell_idx_in_sample = self.__find_in_layout__(cell_id, self.sample.basic_cells_layout)
                    mf_idx_in_cell = self.__find_in_layout__(mf_id, self.sample.cell[cell_id].markerfields_layout)

                # marker selection: trying to find the closest points to a given array of approximate coordinates:
                for qd_num, (qd_key, qd_search_pos) in enumerate(self.qd_search_positions.items()):
                    if keep_images:
                        cv2.drawMarker(img=detected_image,
                                       position=(int(qd_search_pos[1]/um_px_scaling+marker_field_origin_px_pos[1]),
                                                 int(qd_search_pos[0]/um_px_scaling+marker_field_origin_px_pos[0])),
                                       color=(256, 256, 256),
                                       markerSize=15)

                    # determine squared euclidean distance between points
                    squared_distance_array = np.sum((qd_relative_um_distance - qd_search_pos) ** 2, axis=1)
                    # find index of element closest to point that is not too close to previous points
                    closest_found_qd_idx = squared_distance_array.argmin()
                    # for closest_found_qd_idx in np.argsort(squared_distance_array):
                    closest_qd_um_distance = qd_relative_um_distance[closest_found_qd_idx]
                    closest_qd_px_pos = peak_list[closest_found_qd_idx]  # Peak list = Y, X !

                    # making sure that the search range is not exceeded
                    if round(np.sqrt(np.sum((closest_qd_um_distance - qd_search_pos) ** 2)), 1) \
                            <= self.qd_analysis_parameter["qd_search_range"]:
                        # print(closest_qd_um_distance)

                        # determining all relative positions:
                        qd_cell_origin_distance = closest_qd_um_distance + mf_size * mf_idx_in_cell + cell_origin_offset
                        qd_sample_origin_distance = qd_cell_origin_distance - cell_origin_offset + cell_to_cell_distance * cell_idx_in_sample + sample_origin_offset

                        if save_additional_info:
                            self.selected_qd_dict[qd_key] = np.asarray([closest_qd_um_distance, closest_qd_px_pos,
                                                                        qd_cell_origin_distance, qd_sample_origin_distance])
                        else:
                            self.selected_qd_dict[qd_key] = np.asarray([closest_qd_um_distance, closest_qd_px_pos])

                        if keep_images:
                            cv2.drawMarker(img=detected_image,
                                           position=(int(closest_qd_px_pos[1]), int(closest_qd_px_pos[0])),
                                           color=(0, 0, 0),
                                           markerSize=20)
                    else:
                        # print("Distance is: ", round(np.sqrt(np.sum((closest_qd_um_distance - qd_search_pos) ** 2)), 1))
                        if save_additional_info:
                            self.selected_qd_dict[qd_key] = np.asarray([["-", "-"]] * 4)
                        else:
                            self.selected_qd_dict[qd_key] = np.asarray([["-", "-"]] * 2)
                        if self.qd_analysis_parameter["mark_discarded_qd_positions"]:
                            cv2.drawMarker(img=detected_image,
                                           position=(int(closest_qd_px_pos[1]), int(closest_qd_px_pos[0])),
                                           color=(0, 0, 0),
                                           markerSize=10)

                with open(csv_file_name, 'a', newline='') as csv_f:
                    csv_qd_writer = csv.writer(csv_f, delimiter=";")
                    for qd_key, qd_positions in self.selected_qd_dict.items():
                        csv_qd_writer.writerow([sample_name, cell_id, mf_id, qd_key, qd_positions[0, 0], qd_positions[0, 1]])

                if keep_images:
                    cv2.imwrite(f"{sample_name}_{cell_id}-{mf_id.lower()}_detections.png", detected_image)

        self.device_available_signal.emit(self.name, True)
        self.finished_signal.emit()
        return True

    def calc_measurement_precision_error(self, data, sample_name="noname", fit_lim=250, read_fits=False, marker_type_dark=True, save_as_text=True, plot_results=False, verbose=False):
        t1 = time.perf_counter()
        # sample_name = "SA656_3R2"  # // what to do with this?

        if isinstance(data, (str, PurePath)):
            os.chdir(data)
        else:
            print(f"Var 'data' is {type(data)}, but needs to be {str} or {PurePath}")
            print("Calculating MPE from arrays has not been implemented yet")
            return

        # find image files --> this should yield a list of 30 images
        if read_fits:
            file_suffix = "fits"
        else:
            file_suffix = "png"

        mf_image_list = sorted(glob(os.path.join(data, '**', f'*.{file_suffix}'), recursive=True))
        if len(mf_image_list) == 0:
            info = f"No {file_suffix}-images found at {data}"
            print(info)
            self.error_signal.emit((self.name, False, info))
            return -1

        qd_array = None

        if verbose:
            file_generator = enumerate(tqdm(mf_image_list))
        else:
            file_generator = enumerate(mf_image_list)

        # saving files
        save_path = os.path.join(os.getcwd(), f"{sample_name}_saves")
        os.makedirs(save_path, exist_ok=True)
        if False:
            os.chdir(save_path)
            with open("Marker_poses_mv3.txt", 'w') as f:
                f.write("Marker coords V3: \n\n")
            os.chdir(data)

        for n_img, mf_file in file_generator:
            filename = mf_file.split("\\")[-1]
            if verbose:
                print("Time stats: ", filename)

            read_state = self.read_image(filename, read_fits=read_fits)
            if not read_state:
                info = f"Image {filename} skipped, because read_image failed"
                self.error_signal.emit((self.name, True, info))
                continue

            # find marker centers
            self.image = square_image(self.image)
            marker_centers_px_pos = find_marker_origin_conventional(self.image, marker_type_dark=marker_type_dark,
                                                                    show_images=False)  # at 0° rotation: tl = 0, tr = 1, bl = 2, br = 3

            # MF image could not be detected properly, meaning less then 4 corners were detected
            if marker_centers_px_pos is False:
                info = f"Marker Field {filename} skipped, because marker detection failed"
                self.error_signal.emit((self.name, True, info))
                continue

            if False:
                os.chdir(save_path)
                with open("Marker_poses_mv3.txt", 'a') as f:
                    for marker_entry in marker_centers_px_pos:
                        f.write(str(marker_entry))
                        f.write("\n")
                    f.write("\n")
                os.chdir(data)

            # determine the reference marker (origin)
            marker_center_origin_idx = 2  # // replace with read-in parameters
            marker_field_origin_px_pos = marker_centers_px_pos[marker_center_origin_idx]

            marker_min_y = marker_centers_px_pos[:, 0].min()
            marker_max_y = marker_centers_px_pos[:, 0].max()
            marker_min_x = marker_centers_px_pos[:, 1].min()
            marker_max_x = marker_centers_px_pos[:, 1].max()

            # determining the scaling factor how [µm/px]
            # Assumption that no axis-specific distortion appears; taking the median to avoid outliers
            # // ToDO: replace with best distance?
            marker_distances = np.asarray([
                math.hypot(marker_centers_px_pos[1, 1] - marker_centers_px_pos[0, 1],
                           marker_centers_px_pos[1, 0] - marker_centers_px_pos[0, 0]),
                math.hypot(marker_centers_px_pos[3, 1] - marker_centers_px_pos[2, 1],
                           marker_centers_px_pos[3, 0] - marker_centers_px_pos[2, 0]),
                math.hypot(marker_centers_px_pos[2, 1] - marker_centers_px_pos[0, 1],
                           marker_centers_px_pos[2, 0] - marker_centers_px_pos[0, 0]),
                math.hypot(marker_centers_px_pos[3, 1] - marker_centers_px_pos[1, 1],
                           marker_centers_px_pos[3, 0] - marker_centers_px_pos[1, 0])])

            median_marker_px_distance = np.median(marker_distances)
            marker_px_distance_x = np.mean(marker_distances[:2])
            marker_px_distance_y = np.mean(marker_distances[2:])
            um_px_x_scale = 50 / marker_px_distance_x
            um_px_y_scale = 50 / marker_px_distance_y

            # find QDs, parameters found via trial and error (function returns peak as Y, X)
            # peak_list, sigmas = get_qd_positions_from_image(self.image, gauss_fit_iterations=3, area=17, thresh=(
            #             np.mean(self.image) * 1.26))  # // needs to be a parameter

            peak_list, sigmas = get_image_qd(self.image, marker_centers_px_pos, area=10)

            # Consider only QDs within marker field!!!
            peak_list = peak_list[np.where(np.logical_and(np.logical_and(marker_min_y < peak_list[:, 0],
                                                                         peak_list[:, 0] < marker_max_y),
                                                          np.logical_and(marker_min_x < peak_list[:, 1],
                                                                         peak_list[:, 1] < marker_max_x)))]
            if n_img == 0:
                initial_qds = peak_list.shape[0]

            # determine relative px-position: absolute allows calculation without differentiation of positive & negative values
            qd_rel_px_pos = np.absolute(marker_field_origin_px_pos - peak_list)  # returns [y_image, x_image]

            qd_relative_um_distance = np.zeros_like(qd_rel_px_pos, dtype=np.float64)
            # calculate physical distances of QDs relative to marker field origin [µm]
            # convert px-data [y, x] to physical coordiantes [x, y] --> swap axes
            # if the image is rotated by 90° or 270°, axes need to be swapped again --> cancels out previous step
            if marker_center_origin_idx == 1 or marker_center_origin_idx == 2:
                # coordinate system [x, y] lines up with image, but image axes [y, x] need to be transformed [x, y]
                # qd_relative_um_distance = qd_rel_px_pos[:, ::-1] * um_px_scaling
                qd_relative_um_distance[:, 0] = qd_rel_px_pos[:, 1] * um_px_x_scale
                qd_relative_um_distance[:, 1] = qd_rel_px_pos[:, 0] * um_px_y_scale

            else:
                # qd_relative_um_distance = qd_rel_px_pos * um_px_scaling
                qd_relative_um_distance[:, 0] = qd_rel_px_pos[:, 0] * um_px_y_scale
                qd_relative_um_distance[:, 1] = qd_rel_px_pos[:, 1] * um_px_x_scale

            # // ToDO: keep parameters as class --> should be read at the beginning / maybe differentiate between rotation & basic

            # iterating through the relative coordinates of all QDs in one image (for comparability)
            for qd in qd_relative_um_distance:
                # the first image determines the inital size of the qd_array
                if n_img == 0:
                    # values are stored as an array: [QD_number, image_number, (x, y)]
                    qd_array = np.zeros((len(qd_relative_um_distance), len(mf_image_list), 2))
                    qd_array[:, n_img, :] = qd_relative_um_distance
                else:
                    # find the closest QD in the initial QD list (hence 0)
                    squared_distance_array = np.sum((qd_array[:, 0, :] - qd) ** 2, axis=1)
                    qd_idx = squared_distance_array.argmin()

                    # print(f"{qd} is closest to {qd_array[qd_idx, 0, :]}")

                    # if distance is smaller than threshold, consider point to be found again
                    if np.sqrt(squared_distance_array.min()) <= 1.5:  # // replace with variable
                        qd_array[qd_idx, n_img, :] = qd
                    else:
                        # point did not match, is appended (in the QD column of the first image asfor easier referencing)
                        new_value = np.zeros((1, len(mf_image_list), 2))
                        new_value[0, 0, :] = qd
                        qd_array = np.append(qd_array, new_value, axis=0)
            t7 = time.perf_counter()
            #print(f"\tReading image: {t3-t2}")
            #print(f"\tMarker detection: {t4-t3}")
            #print(f"\tPeak finder: {t5-t4}")
            #print(f"\tKalddaradatsch: {t6-t5}")
            #print(f"\tCalculations: {t7-t6}")

        # final assessment
        if qd_array is None:
            info = f"Not enough QDs detected"
            print(info)
            self.error_signal.emit((self.name, False, info))
            return -1

        qd_value_array = []

        # data_mean = qd_array.sum(1) / (qd_array != 0).sum(1)

        for q, qd in enumerate(qd_array):
            found_qds = qd[np.where(qd.sum(1) != 0)]
            # only take variances into account that are statistically relevant
            if found_qds.shape[0] < int(len(mf_image_list)*2/3) - 1:  # QD must be detected in at least 2/3s of all images
                 continue
            # qd_mean = np.mean(found_qds, axis=0)
            # qd_std = np.std(found_qds, axis=0)  # creates biased estimator (population, not sample)
            qd_covars = np.cov(np.transpose(found_qds))
            if qd_covars.any() is np.NaN:
                print("Condition met")
                continue
            qd_value_array.append([np.sqrt(qd_covars[0, 0]), np.sqrt(qd_covars[1, 1]), found_qds.shape[0]])


        qd_value_array = np.asarray(qd_value_array)

        t8 = time.perf_counter()

        if plot_results:
            print(f"Total calculation time: {t8-t1}")
            print("Initial QDs: ", initial_qds)
            print("Final QDs: ", qd_array.shape[0])

        radial_std = np.sqrt(qd_value_array[:,0]**2 + qd_value_array[:,1]**2)
        radial_freq_array = np.stack((radial_std, qd_value_array[:, 2]), axis=1)

        if save_as_text:
            os.chdir(save_path)
            np.savetxt(f'{sample_name}_qd_stds_list.txt', qd_value_array, delimiter=',')
            np.savetxt(f'{sample_name}_qd_radial_stds_counts_list.txt', radial_freq_array, delimiter=',')
            np.save("plot_data.npy", radial_freq_array)

        rad_std_nm = np.round(radial_freq_array[:, 0] * 1000)
        binwidth = 1  # binning in 1 nm
        bins = range(0, int(rad_std_nm.max()) + binwidth, binwidth)

        hist, edges = np.histogram(rad_std_nm, bins=bins)

        bin_centers = edges[:-2] - binwidth / 2
        mpe_value = None

        # init_param = [100, np.log(hist.argmax()), 0.0025]
        try:
            lognormal_parameters, cov = curve_fit(log_normal, bin_centers[:fit_lim], hist[:fit_lim])
            amp, mean, std = lognormal_parameters
            if plot_results:
                print("Curve parameters: ", amp, mean, std)
        except RuntimeError:
            # Curve fit failed, do not add point as center
            print("fit failed")
            mpe_value = -1

        if mpe_value is None:
            pred_vals = log_normal(bin_centers, amp, mean, std)
            mpe_value = pred_vals.argmax()
            if plot_results:
                print(mpe_value)

        if plot_results:
            if mpe_value != -1:
                plt.plot(bin_centers, pred_vals, 'red')
            plt.bar(bin_centers, hist, width=1)
            plt.xlim([0, fit_lim])
            plt.xlabel(r"Radial standard deviation $\sigma$ of position [nm]")
            plt.ylabel("Frequency of position deviation")
            plt.title("Distribution of radial standard deviations")
            if save_as_text:
                os.chdir(save_path)
                plt.savefig(f"{sample_name}_radial_std_histogram.png")
                np.savetxt(f"{sample_name}_MPE.txt", np.asarray([mpe_value, ]))
                os.chdir(data)
            plt.show()

        return mpe_value

    def mpe_from_data(self, path):
        os.chdir(path)

        radial_freq_array = np.load("plot_data.npy")

        """Add here the graph path etc. """

    def calculate_directory_mpe(self, dir_search_path, sample_name=None):
        cell_directory_list = sorted(glob(os.path.join(dir_search_path, "*", "")))

        # making sure there actually are directories
        if len(cell_directory_list) == 0:
            # self.error_signal.emit((self.name, False, "Could not find any cell-subdirectories in path"))
            # self.device_available_signal.emit(self.name, True)
            # self.finished_signal.emit()
            return

        mpe_array = np.zeros((len(cell_directory_list), ))

        os.chdir(dir_search_path)
        text_file_name = "MPE_stat_list.txt"
        with open(text_file_name, 'w') as f:
            f.write("Measurement precision estimate\tFile\n")

        for idx, cell_directory in enumerate(tqdm(cell_directory_list)):
            if sample_name is None:
                sample_name = str(Path(cell_directory.name))

            mpe = self.calc_measurement_precision_error(data=cell_directory, sample_name=sample_name, fit_lim=600,
                                                        read_fits=False, save_as_text=True, plot_results=True)
            print(f"Mpe for {Path(cell_directory).name}:\t{mpe}")
            mpe_array[idx] = mpe

            os.chdir(dir_search_path)
            with open(text_file_name, 'a') as f:
                f.write(f"{mpe}\t{cell_directory}\n")

        return mpe_array


    def detect_dots(self):
        if self.sample is None:
            return

        self.device_available_signal.emit(self.name, False)

        if self.file_list is None:
            state = self.create_file_list(self.image_path)
            if not state:
                # error signal has been sent by function
                self.device_available_signal.emit(self.name, True)
                self.finished_signal.emit()
                return
        os.chdir(self.image_path)

        # differentiating between auto-detectable images and images that need to be detected by hand or are not useful
        valid_image_path = os.path.join(self.image_path, "auto_detection")
        os.makedirs(valid_image_path, exist_ok=True)
        manual_detection_path = os.path.join(self.image_path, "manual_detection")
        os.makedirs(manual_detection_path, exist_ok=True)

        # create rotation matrix to transform coordinates according to sample orientation
        rot_angle = (4 - self.sample.n_rot_steps) / 2 * np.pi
        rotation_matrix = np.asarray([[np.cos(rot_angle), -np.sin(rot_angle)], [np.sin(rot_angle), np.cos(rot_angle)]])

        marker_origin_index = 2  # due to the structure of the marker detector, the bottom-left cross has the index 2

        # detect all pixel coordinates and save them in the sample object
        for i, filepath in enumerate(self.file_list):

            self.progress_signal.emit(i)

            filename = filepath.split("\\")[-1]
            # determine image_id based on name cell-x_cell-y_mf-x_mf_y.png
            mf_information = filename.split(".")[0].split("_")  # cuts off the ending (.png) and then split
            cell_id = mf_information[0]
            mf_id = mf_information[1]

            # no cell origin detected, cell is skipped
            if not self.sample.cell[cell_id].detection_succssful:
                # no detection, mark MF
                self.sample.cell[cell_id].marker_field[mf_id].detection_successful = False
                shutil.move(filepath, os.path.join(manual_detection_path, filename))
                self.error_signal.emit(self.name, True, f"Marker Field {mf_id} skipped, because cell origin of {cell_id} was not detectable")
                continue

            status = self.read_image(filename)
            if not status:
                info = f"Marker Field {mf_id} of cell {cell_id} skipped, because read_image failed"
                self.sample.cell[cell_id].marker_field[mf_id].detection_successful = False
                if mf_id == self.sample.origin_mf_id:
                    info = f"Marker Field {mf_id} (origin) of cell {cell_id} could not be read. Entire cell will be skipped!"
                    self.sample.cell[cell_id].detection_succssful = False
                self.error_signal.emit(self.name, True, info)
                shutil.move(filepath, os.path.join(manual_detection_path, filename))
                continue

            marker_centers_px_pos = find_marker_origin_conventional(self.image)

            # MF image could not be detected properly, meaning less then 4 corners were detected
            if marker_centers_px_pos is False:
                # MF is useless - skip
                info = f"Crosses of Marker Field {mf_id} of cell {cell_id} could not be detected. "
                if mf_id == self.sample.origin_mf_id:
                    info += f"Cell origin could not be detected. All MFs are skipped"
                    self.sample.cell[cell_id].detection_succssful = False
                else:
                    info += "MF is skipped."
                self.sample.cell[cell_id].marker_field[mf_id].detection_successful = False
                self.error_signal.emit(self.name, True, info)
                continue

            # Marker detection function yields 4 pairs of y,x pairs, (top-left, top-right, bottom-left, ...)
            # since sample is rotated, the values need to be assigned correctly
            marker_idx_array = np.arange(4).reshape(2, 2)
            marker_idx_array = np.rot90(marker_idx_array, k=self.sample.n_rot_steps, axes=(1, 0))
            for idx, marker_center in zip(np.nditer(marker_idx_array, order='C'), marker_centers_px_pos):
                # now the marker centers order follows the same system as the marker detection (at 0° rotation)
                # (top-left, top-right, bottom-left, ...)
                self.sample.cell[cell_id].marker_field[mf_id].marker_centers_px[idx, :] = marker_center

            # calculate image scaling factor (um/px)
            # idea: markers are always mf_size apart, the mean of the differences of pixel positions can then be used as a scaling factor
            # marker detector always yields the marker centers in the same order --> differentiation between x and y possible
            # careful: image might still be rotated --> consideration not necessary here yet because we work on pixel level
            self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_y = (self.sample.mf_size/(marker_centers_px_pos[1,0] - marker_centers_px_pos[0,0]) +
                                                                                    self.sample.mf_size/(marker_centers_px_pos[3, 0] - marker_centers_px_pos[2,0])) \
                                                                                   / 2
            self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_x = (self.sample.mf_size/(marker_centers_px_pos[2, 1] - marker_centers_px_pos[0,1]) +
                                                                                    self.sample.mf_size/(marker_centers_px_pos[3, 1] - marker_centers_px_pos[1,1])) \
                                                                                   / 2

            # determine cell origin
            # rotated to enable easier, relative coordiante calculations and revisiting
            if mf_id == self.sample.origin_mf_id:
                origin_coords = self.__convert_px_to_abs_coordinates__(
                    px_coord_array=self.sample.cell[cell_id].marker_field[mf_id].marker_centers_px[2],
                    reference_abs_pos=self.sample.cell[cell_id].marker_field[mf_id].camera_pos,
                    reference_px_coords=np.asarray([self.img_height, self.img_width]),
                    um_px_scale_y=self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_y,
                    um_px_scale_x=self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_x
                )
                if origin_coords is None:
                    self.device_available_signal.emit(self.name, True)
                    self.finished_signal.emit()
                    return
                self.sample.cell[cell_id].cell_origin[:] = origin_coords @ rotation_matrix


            # detect QDs
            peak_list, sigmas = get_qd_positions_from_image(self.image)

            # add suitable QDs into dict
            # forbidden area:
            distance_from_markers = 40
            y_min = marker_centers_px_pos[0, 1] + distance_from_markers
            y_max = marker_centers_px_pos[3, 1] - distance_from_markers
            x_min = marker_centers_px_pos[0, 0] + distance_from_markers
            x_max = marker_centers_px_pos[3, 0] - distance_from_markers


            qd_counter = 0  # chosen instead of enumerate to guarantee consecutive rising numbers even if some QDs are neglected
            for y, x in peak_list:
                # dot is not in "forbidden" area
                if (y_min < y < y_max) and (x_min < x < x_max):
                    self.sample.cell[cell_id].marker_field[mf_id].addQD(f"{cell_id}_{mf_id}_{qd_counter}", x, y)
                    qd_counter += 1

            # move images into directory for final calculations
            shutil.move(filepath, os.path.join(valid_image_path, filename))

        # calculating relative positions between objects
        _ = self.create_file_list(path=valid_image_path)
        os.chdir(valid_image_path)

        for i, filepath in enumerate(self.file_list):
            filename = filepath.split("\\")[-1]
            # determine image_id based on name cell-x_cell-y_mf-x_mf_y.png
            mf_information = filename.split(".")[0].split("_")  # cuts off the ending (.png) and then split
            cell_id = mf_information[0]
            mf_id = mf_information[1]

            self.progress_signal.emit(i)

            # making sure that calculation make sense
            if not self.sample.cell[cell_id].marker_field[mf_id].detection_successful:
                info = f"Detection impossible, because MF {mf_id} of cell {cell_id} not properly imaged"
                self.error_signal.emit(self.name, True, info)
                continue

            # convert marker centers
            for i, marker_px_coord in enumerate(self.sample.cell[cell_id].marker_field[mf_id].marker_centers_px):
                marker_abs_pos = self.__convert_px_to_abs_coordinates__(
                px_coord_array=marker_px_coord,
                reference_abs_pos=self.sample.cell[cell_id].marker_field[mf_id].camera_pos,
                reference_px_coords=np.asarray([self.img_height, self.img_width]),
                um_px_scale_y=self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_y,
                um_px_scale_x=self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_x
                )
                if marker_abs_pos is None:
                    self.device_available_signal.emit(self.name, True)
                    self.finished_signal.emit()
                    return

                marker_pos_rotated = rotation_matrix @ marker_abs_pos
                if i == marker_origin_index:
                    self.sample.cell[cell_id].marker_field[mf_id].origin_absolute_position[:] = marker_pos_rotated

                    # save in file
                    if mf_id == self.sample.origin_mf_id:
                        self.save_information_to_file(filename=f"{self.sample.name}_{cell_id}.txt",
                                                      object_name=f"Cell origin {cell_id}",
                                                      coordinates=marker_pos_rotated)

                    self.save_information_to_file(filename=f"{self.sample.name}_{cell_id}.txt",
                                                  object_name=f"Marker origin {mf_id}",
                                                  coordinates=marker_pos_rotated)

                rel_marker_distance = marker_pos_rotated - self.sample.cell[cell_id].cell_origin
                self.sample.cell[cell_id].marker_field[mf_id].marker_pos_relative_cell[i, :] = rel_marker_distance

                # save in file
                self.save_information_to_file(filename=f"{self.sample.name}_{cell_id}.txt",
                                              object_name=f"Marker center {mf_id}-{i}",
                                              coordinates=rel_marker_distance)

            # QDs (relative to MF, relative to cell)
            for qd_key, qd_obj in self.sample.cell[cell_id].marker_field[mf_id].QD.itmes():
                # absolute position (unrotated)

                abs_stage_pos_coords = self.__convert_px_to_abs_coordinates__(
                    px_coord_array=np.asarray([qd_obj.ypos, qd_obj.xpos]),
                    reference_abs_pos=self.sample.cell[cell_id].marker_field[mf_id].camera_pos,
                    reference_px_coords=np.asarray([self.img_height, self.img_width]),
                    um_px_scale_y=self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_y,
                    um_px_scale_x=self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_x
                )
                if abs_stage_pos_coords is None:
                    self.device_available_signal.emit(self.name, True)
                    self.finished_signal.emit()
                    return
                # transfer coordinates into QD
                qd_obj.absolute_stage_pos[:] = abs_stage_pos_coords

                qd_rotated_pos = rotation_matrix @ qd_obj.absolute_stage_pos
                qd_obj.distance_relative_mf[:] = qd_rotated_pos - self.sample.cell[cell_id].marker_field[mf_id].origin_absolute_position
                qd_obj.distance_relative_cell[:] = qd_rotated_pos - self.sample.cell[cell_id].cell_origin

                # save in file
                self.save_information_to_file(filename=f"{self.sample.name}_{cell_id}.txt",
                                              object_name=f"QD {qd_key} relative to MF",
                                              coordinates=qd_obj.distance_relative_mf)

        self.device_available_signal.emit(self.name, True)
        self.finished_signal.emit()
        return True

    def detect_dots_stage_pos(self, calc_absolute_stage_pos=True):
        # ToDO: Heavy rework!!!
        self.device_available_signal.emit(self.name, False)
        if self.file_list is None:
            file_list_state = self.create_file_list(self.image_path)
            if not file_list_state:
                # error signal has been sent by function
                self.device_available_signal.emit(self.name, True)
                self.finished_signal.emit()
                return
        os.chdir(self.image_path)

        # // /*

        # differentiating between auto-detectable images and images that need to be detected by hand or are not useful
        valid_image_path = os.path.join(self.image_path, "auto_detection")
        os.makedirs(valid_image_path, exist_ok=True)
        manual_detection_path = os.path.join(self.image_path, "manual_detection")
        os.makedirs(manual_detection_path, exist_ok=True)

        # create rotation matrix to transform coordinates according to sample orientation
        # (so that the detected coordinates of the markers are all at 0° rotation)
        rot_angle = (4 - self.sample_n_rot_steps) / 2 * np.pi  # //
        rotation_matrix = np.asarray([[np.cos(rot_angle), -np.sin(rot_angle)], [np.sin(rot_angle), np.cos(rot_angle)]])

        marker_origin_index = 2  # due to the structure of the marker detector, the bottom-left cross has the index 2

        # detect all pixel coordinates and save them in the sample object
        for i, filepath in enumerate(self.file_list):

            self.progress_signal.emit(i)

            filename = filepath.split("\\")[-1]
            # determine image_id based on name: "cell-x_cell-y_mf-x_mf_y.png"
            mf_information = filename.split(".")[0].strip.split("_")  # cuts off the ending (.png) and then split
            cell_id = mf_information[0]
            mf_id = mf_information[1]

            read_state = self.read_image(filename)
            if not read_state:
                info = f"Marker Field {mf_id} of cell {cell_id} skipped, because read_image failed"
                # self.sample.cell[cell_id].marker_field[mf_id].detection_successful = False # //
                self.error_signal.emit((self.name, True, info))
                shutil.move(filepath, os.path.join(manual_detection_path, filename))
                continue

            # detect marker centers via Hough transform & Gaussian fit
            # image is not rotated to make sure the actual stage positions do match reality
            marker_centers_px_pos = find_marker_origin_conventional(
                self.image)  # at 0° rotation: tl = 0, tr = 1, bl = 2, br = 3

            # MF image could not be detected properly, meaning less then 4 corners were detected
            if marker_centers_px_pos is False:
                # MF is useless - skip
                info = f"Crosses of Marker Field {mf_id} of cell {cell_id} could not be detected. MF is skipped."
                # self.sample.cell[cell_id].marker_field[mf_id].detection_successful = False #//
                shutil.move(filepath, os.path.join(manual_detection_path, filename))
                self.error_signal.emit((self.name, True, info))
                continue

            # Marker detection function yields 4 pairs of y,x pairs, (top-left, top-right, bottom-left, ...)
            # since sample is rotated, the values need to be assigned correctly
            marker_idx_array = np.arange(4).reshape(2, 2)
            marker_idx_array = np.rot90(marker_idx_array, k=self.sample_n_rot_steps, axes=(1, 0))
            # going through detected markers (which might have been rotated) and the respectively rotated indices & matching them
            for idx, marker_center in zip(np.nditer(marker_idx_array, order='C'), marker_centers_px_pos):
                # now the marker centers order follows the same system as the marker detection (at 0° rotation)
                # (top-left, top-right, bottom-left, ...)
                # self.sample.cell[cell_id].marker_field[mf_id].marker_centers_px[idx, :] = marker_center //
                self.center_px_pos[i, idx, :] = marker_center

            # calculate image scaling factor (um/px)
            # idea: markers are always mf_size apart, the mean of the differences of pixel positions can then be used as a scaling factor
            # marker detector always yields the marker centers in the same order --> differentiation between x and y possible
            # careful: image might still be rotated --> consideration not necessary here yet because we work on pixel level
            # // /*
            self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_x = (self.sample.mf_size / (
                        marker_centers_px_pos[1, 1] - marker_centers_px_pos[0, 1]) +
                                                                                    self.sample.mf_size / (
                                                                                                marker_centers_px_pos[
                                                                                                    3, 1] -
                                                                                                marker_centers_px_pos[
                                                                                                    2, 1])) \
                                                                                   / 2
            self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_y = (self.sample.mf_size / (
                        marker_centers_px_pos[2, 0] - marker_centers_px_pos[0, 0]) +
                                                                                    self.sample.mf_size / (
                                                                                                marker_centers_px_pos[
                                                                                                    3, 0] -
                                                                                                marker_centers_px_pos[
                                                                                                    1, 0])) \
                                                                                   / 2
            self.px_scaling_x[i] = (self.sample.mf_size / (marker_centers_px_pos[1, 1] - marker_centers_px_pos[0, 1]) +
                                    self.sample.mf_size / (marker_centers_px_pos[3, 1] - marker_centers_px_pos[2, 1])) \
                                   / 2
            self.px_scaling_y[i] = (self.sample.mf_size / (marker_centers_px_pos[2, 0] - marker_centers_px_pos[0, 0]) +
                                    self.sample.mf_size / (marker_centers_px_pos[3, 0] - marker_centers_px_pos[1, 0])) \
                                   / 2

            # absolute calculations (not necessary / possible when completely offline) # //
            if calc_absolute_stage_pos:
                # determine cell origin
                # rotated to enable easier, relative coordiante calculations and revisiting
                if mf_id == self.sample.origin_mf_id:
                    origin_coords = self.__convert_px_to_abs_coordinates__(
                        px_coord_array=self.sample.cell[cell_id].marker_field[mf_id].marker_centers_px[
                            marker_origin_index],
                        reference_abs_pos=self.sample.cell[cell_id].marker_field[mf_id].camera_pos,
                        reference_px_coords=np.asarray([self.img_height, self.img_width]),
                        um_px_scale_y=self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_y,
                        um_px_scale_x=self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_x
                    )
                    if origin_coords is None:
                        self.device_available_signal.emit(self.name, True)
                        self.finished_signal.emit()
                        return

                    # rotate marker coordinates so that --> // why?
                    self.sample.cell[cell_id].cell_origin[:] = origin_coords @ rotation_matrix

            # detect QDs
            peak_list, sigmas = get_qd_positions_from_image(self.image)

            # add suitable QDs into dict
            # forbidden area:
            distance_from_markers = 40
            y_min = marker_centers_px_pos[0, 1] + distance_from_markers
            y_max = marker_centers_px_pos[3, 1] - distance_from_markers
            x_min = marker_centers_px_pos[0, 0] + distance_from_markers
            x_max = marker_centers_px_pos[3, 0] - distance_from_markers

            qd_counter = 0  # chosen instead of enumerate to guarantee consecutive rising numbers even if some QDs are neglected
            for y, x in peak_list:
                # dot is not in "forbidden" area
                if (y_min < y < y_max) and (x_min < x < x_max):
                    self.sample.cell[cell_id].marker_field[mf_id].addQD(f"{cell_id}_{mf_id}_{qd_counter}", x, y)
                    qd_counter += 1

            # move images into directory for final calculations
            shutil.move(filepath, os.path.join(valid_image_path, filename))

        # calculating relative positions between objects
        _ = self.create_file_list(path=valid_image_path)
        os.chdir(valid_image_path)

        for i, filepath in enumerate(self.file_list):
            filename = filepath.split("\\")[-1]
            # determine image_id based on name cell-x_cell-y_mf-x_mf_y.png
            mf_information = filename.split(".")[0].split("_")  # cuts off the ending (.png) and then split
            cell_id = mf_information[0]
            mf_id = mf_information[1]

            self.progress_signal.emit(i)

            # making sure that calculation make sense
            if not self.sample.cell[cell_id].marker_field[mf_id].detection_successful:
                info = f"Detection impossible, because MF {mf_id} of cell {cell_id} not properly imaged"
                self.error_signal.emit((self.name, True, info))
                continue

            # convert marker centers
            for i, marker_px_coord in enumerate(self.sample.cell[cell_id].marker_field[mf_id].marker_centers_px):
                marker_abs_pos = self.__convert_px_to_abs_coordinates__(
                    px_coord_array=marker_px_coord,
                    reference_abs_pos=self.sample.cell[cell_id].marker_field[mf_id].camera_pos,
                    reference_px_coords=np.asarray([self.img_height, self.img_width]),
                    um_px_scale_y=self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_y,
                    um_px_scale_x=self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_x
                )
                if marker_abs_pos is None:
                    self.device_available_signal.emit(self.name, True)
                    self.finished_signal.emit()
                    return

                marker_pos_rotated = rotation_matrix @ marker_abs_pos
                if i == marker_origin_index:
                    self.sample.cell[cell_id].marker_field[mf_id].origin_absolute_position[:] = marker_pos_rotated

                    # save in file
                    if mf_id == self.sample.origin_mf_id:
                        self.save_information_to_file(filename=f"{self.sample.name}_{cell_id}.txt",
                                                      object_name=f"Cell origin {cell_id}",
                                                      coordinates=marker_pos_rotated)

                    self.save_information_to_file(filename=f"{self.sample.name}_{cell_id}.txt",
                                                  object_name=f"Marker origin {mf_id}",
                                                  coordinates=marker_pos_rotated)

                rel_marker_distance = marker_pos_rotated - self.sample.cell[cell_id].cell_origin
                self.sample.cell[cell_id].marker_field[mf_id].marker_pos_relative_cell[i, :] = rel_marker_distance

                # save in file
                self.save_information_to_file(filename=f"{self.sample.name}_{cell_id}.txt",
                                              object_name=f"Marker center {mf_id}-{i}",
                                              coordinates=rel_marker_distance)

            # QDs (relative to MF, relative to cell)
            for qd_key, qd_obj in self.sample.cell[cell_id].marker_field[mf_id].QD.itmes():
                # absolute position (unrotated)

                abs_stage_pos_coords = self.__convert_px_to_abs_coordinates__(
                    px_coord_array=np.asarray([qd_obj.ypos, qd_obj.xpos]),
                    reference_abs_pos=self.sample.cell[cell_id].marker_field[mf_id].camera_pos,
                    reference_px_coords=np.asarray([self.img_height, self.img_width]),
                    um_px_scale_y=self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_y,
                    um_px_scale_x=self.sample.cell[cell_id].marker_field[mf_id].um_px_scaling_factor_x
                )
                if abs_stage_pos_coords is None:
                    self.device_available_signal.emit(self.name, True)
                    self.finished_signal.emit()
                    return
                # transfer coordinates into QD
                qd_obj.absolute_stage_pos[:] = abs_stage_pos_coords

                qd_rotated_pos = rotation_matrix @ qd_obj.absolute_stage_pos
                qd_obj.distance_relative_mf[:] = qd_rotated_pos - self.sample.cell[cell_id].marker_field[
                    mf_id].origin_absolute_position
                qd_obj.distance_relative_cell[:] = qd_rotated_pos - self.sample.cell[cell_id].cell_origin

                # save in file
                self.save_information_to_file(filename=f"{self.sample.name}_{cell_id}.txt",
                                              object_name=f"QD {qd_key} relative to MF",
                                              coordinates=qd_obj.distance_relative_mf)

        self.device_available_signal.emit(self.name, True)
        self.finished_signal.emit()
        return True