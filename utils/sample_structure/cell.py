#this is the template for the write field
from source.sample_structure.markerField import MarkerField
from collections import OrderedDict
import numpy as np
import string

class Cell:
    def __init__(self, name, n_mf_x=6, n_mf_y=6, n_rotations=0):
        self.name = name
        self.n_mf = n_mf_x*n_mf_y
        self.n_mf_x = n_mf_x
        self.n_mf_y = n_mf_y
        self.originX = None
        self.originY = None
        self.markerFieldDict = OrderedDict()
        self.selectedMarkerFieldName = None
    
        self.markerfields_layout = self.createMarkerFieldsLayout()
        self.oriented_mfs_layout = np.rot90(self.markerfields_layout, k=n_rotations, axes=(1, 0))
        self.marker_field = self.createMarkerFieldsDict()
        self.detection_successful = True
        self.cell_origin = np.zeros((1, 2))  # absolute position [um], rotated
        
    def addMarkerField(self, markerFieldName,markerSize):
        self.markerFieldDict[markerFieldName] = MarkerField(markerFieldName,markerSize)

    def getMarkerField(self,markerFieldName):
        return self.markerFieldDict[markerFieldName]
    
    def removeMarkerField(self,markerFieldName):
        if self.selectedMarkerFieldName is markerFieldName:
            self.selectedMarkerFieldName = None
        del self.markerFieldDict[markerFieldName]
        print('deleted field')
        print(self.markerFieldDict)
        
    def chooseSelectedMarkerField(self,markerfieldName):
        self.selectedMarkerFieldName = markerfieldName
        return self.markerFieldDict[markerfieldName]
    
    def returnSelectedMarkerField(self):
        if self.selectedMarkerFieldName is not None:
            return self.markerFieldDict[self.selectedMarkerFieldName]

    def createMarkerFieldsLayout(self):
        letters = np.asarray(list(string.ascii_uppercase), dtype=str)

        mfs_rows = letters[:int((self.n_mf_y))]
        mfs_columns = list(range(1, int((self.n_mf_x)) + 1))

        mfs_layout = np.zeros((len(mfs_rows), len(mfs_columns)), dtype=object)

        for row_idx, mf_letter in enumerate(mfs_rows):
            for col_idx, mf_no in enumerate(mfs_columns):
                mfs_layout[row_idx, col_idx] = mf_letter + str(mf_no)

        return mfs_layout

    def createMarkerFieldsDict(self):
        mf_dict = dict()
        for row in self.oriented_mfs_layout:
            for mf_id in row:
                mf_dict[mf_id] = MarkerField(name=mf_id)

        return mf_dict
