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

#import sys
#try:
#    sys.path.insert(0,'/home/onerva/consistency-as-onion')
#except:
#    sys.path.insert(0,'/home/onerva/projects/consistency-as-onion')   

import functions
import onion_parameters as params

subjects = [params.testSubjectFolders[0]]
originalROIInfoPath = params.originalROIInfoFile
optimizedROIInfoFile = params.optimizedROIInfoFile
allVoxelTsFile = params.ROIVoxelTsFileName
optimizedConsistencySaveName = params.optimizedSpatialConsistencySaveName
originalConsistencySavePath = params.originalSpatialConsistencySavePath
optimizedCorrelationSaveName = params.optimizedCorrelationSaveName
originalCorrelationSavePath = params.originalCorrelationSavePath
figureSavePath = params.spatialConsistencyValidationPath

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

originalCentroids,_,originalVoxelCoordinates,originalROIMaps = functions.readROICentroids(originalROIInfoPath,readVoxels=True,fixCentroids=True)
centroidIndices = np.zeros(len(originalCentroids),dtype=int)
for i, originalCentroid in enumerate(originalCentroids):
    centroidIndices[i] = np.where((originalVoxelCoordinates==originalCentroid).all(axis=1)==1)[0][0]

for i, subject in enumerate(subjects):
    optimizedROIInfoPath = subject + optimizedROIInfoFile
    allVoxelTsPath = subject + allVoxelTsFile
    
    allVoxelTs = io.loadmat(allVoxelTsPath)['roi_voxel_data'][0]['roi_voxel_ts'][0]
    _,_,voxelCoordinates,ROIMaps = functions.readROICentroids(optimizedROIInfoPath,readVoxels=True)
    
    ROIIndices = []
    for ROIMap in ROIMaps:
        indices = np.zeros(len(ROIMap),dtype=int)
        if len(ROIMap.shape) == 1:
            indices[0] = np.where((voxelCoordinates == ROIMap).all(axis=1)==1)[0][0]
        else:
            for j, voxel in enumerate(ROIMap):
                #print 'something'
                indices[j] = np.where((voxelCoordinates == voxel).all(axis=1)==1)[0][0]
        ROIIndices.append(indices)
        
    spatialConsistencies = functions.calculateSpatialConsistencyInParallel(ROIIndices,allVoxelTs)
    spatialConsistencyData = {'spatialConsistencies':spatialConsistencies,'type':'optimized'}
    savePath = subject + optimizedConsistencySaveName
    with open(savePath, 'wb') as f:
        pickle.dump(spatialConsistencyData, f, -1)
    consistencyDistribution,consistencyBinCenters = functions.getDistribution(spatialConsistencies,params.nConsistencyBins)
    
    correlationsToCentroid = functions.calculateCorrelationsToCentroidInParallel(ROIIndices,allVoxelTs,centroidIndices)
    correlationData = {'correlationsToCentroid':correlationsToCentroid,'type':'optimized'}
    savePath = subject + optimizedCorrelationSaveName
    with open(savePath, 'wb') as f:
        pickle.dump(correlationData, f, -1)
    correlationDistribution,correlationBinCenters = functions.getDistribution(correlationsToCentroid,params.nConsistencyBins)
    
    if i == 0:
        ax1.plot(consistencyBinCenters,consistencyDistribution,color=params.optimizedColor,alpha=params.optimizedAlpha,label='Optimized ROIs')
        ax2.plot(correlationBinCenters,correlationDistribution,color=params.optimizedColor,alpha=params.optimizedAlpha,label='Optimized ROIs')
    else:
        ax1.plot(consistencyBinCenters,consistencyDistribution,color=params.optimizedColor,alpha=params.optimizedAlpha)
        ax2.plot(correlationBinCenters,correlationDistribution,color=params.optimizedColor,alpha=params.optimizedAlpha)   
#_,_,voxelCoordinates,ROIMaps = functions.readROICentroids(originalROIInfoPath,readVoxels=True)

ROIIndices = []
for ROIMap in originalROIMaps:
    indices = np.zeros(len(ROIMap),dtype=int)
    for i, voxel in enumerate(ROIMap):
        indices[i] = np.where((originalVoxelCoordinates == voxel).all(axis=1)==1)[0][0]
    ROIIndices.append(indices)
  
spatialConsistencies = functions.calculateSpatialConsistencyInParallel(ROIIndices,allVoxelTs)
spatialConsistencyData = {'spaitalConsistencies':spatialConsistencies,'type':'original Brainnetome'}
with open(originalConsistencySavePath, 'wb') as f:
        pickle.dump(spatialConsistencyData, f, -1)
consistencyDistribution,consistencyBinCenters = functions.getDistribution(spatialConsistencies,params.nConsistencyBins)

correlationsToCentroid = functions.calculateCorrelationsToCentroidInParallel(ROIIndices,allVoxelTs,centroidIndices)
correlationData = {'correlationsToCentroid':correlationsToCentroid,'type':'original Brainnetome'}
savePath = subject + optimizedCorrelationSaveName
with open(savePath, 'wb') as f:
    pickle.dump(correlationData, f, -1)
correlationDistribution,correlationBinCenters = functions.getDistribution(correlationsToCentroid,params.nConsistencyBins)

ax1.plot(consistencyBinCenters,consistencyDistribution,color=params.originalColor,alpha=params.originalAlpha,label='Original ROIs')
ax2.plot(correlationBinCenters,correlationDistribution,color=params.originalColor,alpha=params.originalAlpha,label='Original ROIs')

ax1.set_xlabel('Spatial consistency')
ax1.set_ylabel('PDF')
ax1.legend()

ax2.set_xlabel('Correllation to ROI centroid')
ax2.set_ylabel('PDF')
ax2.legend()

plt.tight_layout()
plt.savefig(figureSavePath,format='pdf',bbox_inches='tight')
    
