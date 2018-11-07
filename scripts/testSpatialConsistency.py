# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 14:56:38 2018

@author: onerva

A small script for testing the parallelization of consistency calculations (for
a single subject first)
"""
import numpy as np
import matplotlib.pylab as plt
from scipy import io
import cPickle as pickle

import functions
import onion_parameters as params

allVoxelTsPath = params.testSubjectFolders[0] + params.ROIVoxelTsFileName
ROIInfoFile = params.ROIInfoFile

allVoxelTs = io.loadmat(allVoxelTsPath)['roi_voxel_data'][0]['roi_voxel_ts'][0]

_,_,voxelCoordinates,ROIMaps = functions.readROICentroids(ROIInfoFile,readVoxels=True,fixCentroids=True)

ROIIndices = []

for ROIMap in ROIMaps:
    indices = np.zeros(len(ROIMap),dtype=int)
    for i, voxel in enumerate(ROIMap):
        indices[i] = np.where((voxelCoordinates == voxel).all(axis=1)==1)[0][0]
    ROIIndices.append(indices)

testROIs = ROIIndices[0:20]
spatialConsistencies = functions.calculateSpatialConsistencyInParallel(testROIs,allVoxelTs)
