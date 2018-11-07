# -*- coding: utf-8 -*-
"""
Created on Tue Nov  6 14:31:45 2018

@author: onerva

A script for validating the consistency-based optimization of ROIs by comparing
the consistency of optimized ROIs to the consistency of ROIs of the Brainnetome
parcellation.

In future, this script should calculate the consistency distributions form pooled
data of all subjects.
"""
import numpy as np
import matplotlib.pylab as plt
from scipy import io
import cPickle as pickle

import functions
import onion_parameters as params

subjects = params.testSubjectFolders
originalROIInfoPath = params.originalROIInfoFile
optimizedROIInfoFile = params.optimizedROIInfoFile
allVoxelTsFile = params.ROIVoxelTsFileName
optimizedConsistencySaveName = params.optimizedSpatialConsistencySaveName
originalConsistencySavePath = params.originalSpatialConsistencySavePath
figureSavePath = params.spatialConsistencyValidationPath

fig = plt.figure()
ax = fig.add_subplot(111)

for i, subject in enumerate(subjects):
    optimizedROIInfoPath = subject + optimizedROIInfoFile
    allVoxelTsPath = subject + allVoxelTsFile
    
    allVoxelTs = allVoxelTs = io.loadmat(allVoxelTsPath)['roi_voxel_data'][0]['roi_voxel_ts'][0]
    _,_,voxelCoordinates,ROIMaps = functions.readROICentroids(optimizedROIInfoPath,readVoxels=True)
    
    ROIIndices = []
    for ROIMap in ROIMaps:
        indices = np.zeros(len(ROIMap),dtype=int)
        for i, voxel in enumerate(ROIMap):
            indices[i] = np.where((voxelCoordinates == voxel).all(axis=1)==1)[0][0]
        ROIIndices.append(indices)
      
    spatialConsistencies = functions.calculateSpatialConsistencyInParallel(ROIIndices,allVoxelTs)
    spatialConsistencyData = {'spatialConsistencies':spatialConsistencies,'type':'optimized'}
    savePath = subject + optimizedConsistencySaveName
    with open(savePath, 'wb') as f:
        pickle.dump(spatialConsistencyData, f, -1)
    consistencyDistribution,binCenters = functions.getDistribution(spatialConsistencies,params.nConsistencyBins)
    
    if i == 0:
        ax.plot(binCenters,consistencyDistribution,color=params.optimizedColor,alpha=params.optimizedAlpha,label='Optimized ROIs')
    else:
        ax.plot(binCenters,consistencyDistribution,color=params.optimizedColor,alpha=params.optimizedAlpha)
    
_,_,voxelCoordinates,ROIMaps = functions.readROICentroids(originalROIInfoPath,readVoxels=True)

ROIIndices = []
for ROIMap in ROIMaps:
    indices = np.zeros(len(ROIMap),dtype=int)
    for i, voxel in enumerate(ROIMap):
        indices[i] = np.where((voxelCoordinates == voxel).all(axis=1)==1)[0][0]
    ROIIndices.append(indices)
  
spatialConsistencies = functions.calculateSpatialConsistencyInParallel(ROIIndices,allVoxelTs)
spatialConsistencyData = {'spaitalConsistencies':spatialConsistencies,'type':'original Brainnetome'}
with open(originalConsistencySavePath, 'wb') as f:
        pickle.dump(spatialConsistencyData, f, -1)

consistencyDistribution,binCenters = functions.getDistribution(spatialConsistencies)

ax.plot(binCenters,consistencyDistribution,color=params.optimizedColor,alpha=params.optimizedAlpha,label='Original ROIs')

ax.set_xlabel('Spatial consistency')
ax.set_ylabel('PDF')
ax.legend()

plt.tight_layout()
plt.savefig(figureSavePath,format='pdf',bbox_inches='tight')
    
