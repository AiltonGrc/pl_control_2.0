# this class has all important data of
# the sample in it
from source.sample_structure.markerField import MarkerField
from collections import OrderedDict
from util.source.sample_structure.cell import Cell
import numpy as np
import string

class Sample:
    def __init__(self, name, n_cells_x=5,n_cells_y=5, bottom_left_cell="AA", n_mfs_x=6,n_mfs_y=6, mf_size_x=50,mf_size_y=50,
                 inter_cell_distance_x=200,inter_cell_distance_y=200, origin_mf_id="A1"):
        self.name = name
        self.n_cells = n_cells_x*n_cells_y
        self.n_cells_x = n_cells_x
        self.n_cells_y = n_cells_y
        self.n_marker_fields = n_mfs_x*n_mfs_y
        self.basic_cells_layout = self.createCellsLayout()
        self.n_rot_steps, self.oriented_cells_layout = self.createOrientedCellsLayout(bottom_left_cell)
        self.sample_orientation = self.n_rot_steps / 2 * np.pi  # mathematically positive definition (left turns of sample)
        self.n_mf_x = n_mfs_x
        self.n_mf_y = n_mfs_y
        self.n_mf = n_mfs_y*n_mfs_x
        self.mf_size_x = mf_size_x
        self.mf_size_y = mf_size_y
        self.cell = self.createCellsDict(self.n_mf)
        self.inter_cell_distance_x = inter_cell_distance_x
        self.inter_cell_distance_y = inter_cell_distance_y
        self.cell_origin = None
        self.actual_sample_size_x = self.n_cells_x*self.n_mf_x * mf_size_x +(self.n_cells-1)*inter_cell_distance_x
        self.actual_sample_size_y = self.n_cells_y*self.n_mf_y * mf_size_y + (self.n_cells_y-1)*inter_cell_distance_y

        # sample properties that will be calculated
        self.length_x = None
        self.length_y = None
        self.step_deviation_factor_x = 1
        self.step_deviation_factor_y = 1
        self.rotation_matrix = np.asarray([[np.cos(0), -np.sin(0)], [np.sin(0), np.cos(0)]])

        self.triangulation_points = np.zeros((3, 6))  # triang_x, triang_y, piezo, std, stage_x, stage_y
        self.piezo_pos_parameters = np.zeros((3,))
        self.focus_std_parameters = np.zeros((3,))

        self.cellDict = OrderedDict()
        self.selectedCellName = None
        self.origin_mf_id = origin_mf_id

    def addCell(self, cellName):
        self.cellDict[cellName] = Cell(cellName)

    def getCell(self,cellName):
        return self.cellDict[cellName]
    
    def removeCell(self,cellName):
        if self.selectedCellName is cellName:
            self.selectedCellName = None
        del self.cellDict[cellName]
        print('deleted cell')
        print(self.cellDict)
    
    def chooseSelectedCell(self,cellName):
        self.selectedCellName = cellName
        return self.returnSelectedCell()
        
    def returnSelectedCell(self):
        if self.selectedCellName is not None:
            return self.cellDict[self.selectedCellName]

    def createCellsLayout(self):
        letters = np.asarray(list(string.ascii_uppercase), dtype=str)

        cells = letters[:int(np.sqrt(self.n_cells))]
        cells_layout = np.zeros((len(cells), len(cells)), dtype=object)

        for row_idx, cell_letter in enumerate(cells):
            for col_idx, cell_no in enumerate(cells):
                cells_layout[row_idx, col_idx] = cell_letter + cell_no

        return cells_layout

    def createOrientedCellsLayout(self, bl_cell):
        # find location of cell (key) id in basic layout
        key_location_y, key_location_x = np.where(self.basic_cells_layout == bl_cell)

        if key_location_y.shape[0] == 0:
            # in case that ID cannot be found!
            print("Cell ID was not found in layout")
            return 0, self.basic_cells_layout

        # rotate arrays counter-clockwise
        if key_location_y[0] != 0 and key_location_x[0] == 0:
            n_rotation_steps = 1
        elif key_location_y[0] != 0 and key_location_x[0] != 0:
            n_rotation_steps = 2
        elif key_location_y[0] == 0 and key_location_x[0] != 0:
            n_rotation_steps = 3
        else:
            n_rotation_steps = 0

        # rotating the layout counter-clockwise n times by 90°
        oriented_cells_layout = np.rot90(self.basic_cells_layout, k=n_rotation_steps, axes=(1, 0))

        return n_rotation_steps, oriented_cells_layout

    def createCellsDict(self, n_mfs):
        cell_dict = dict()
        for row in self.oriented_cells_layout:
            for cell_id in row:
                cell_dict[cell_id] = Cell(name=cell_id, n_mf=n_mfs, n_rotations=self.n_rot_steps)

        return cell_dict
