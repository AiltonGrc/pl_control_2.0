from source.sample_structure.quantumdot import QuantumDot
from collections import OrderedDict
import numpy as np

class MarkerField:
    def __init__(self,name, markerSize=50):
        self.name = name
        self.markerSize = markerSize
        self.originPosition = None
        self.pictureMarker = None
        self.markerPosition = None  # coordinates at the sample

        self.detection_successful = True
        self.camera_pos = np.zeros((1, 2))  # absolute stage coordinates for each MarkerField
        self.focus_piezo_pos = None

        self.marker_centers_px = np.zeros((4, 2))  # pixel position in the image of the maker centers
        self.marker_pos_relative_cell = np.zeros((4, 2))  # location of the four marker centers relative to cell origin
        self.origin_absolute_position = np.zeros((1, 2))  # location of the chosen origin (marker center bottom-left) in absolute stage coordiantes, rotated
        self.um_px_scaling_factor_x = None  # [um / px]
        self.um_px_scaling_factor_y = None  # [um / px]
        self.QD = dict()

    def addQD(self, name, x_coord, y_coord):
        self.QD[name] = QuantumDot(name, x_coord, y_coord)
    