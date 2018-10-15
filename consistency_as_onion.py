#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 15:07:57 2018

@author: onerva

A script for analysing the dependency of spatial consistency on ROI size. The
idea is to grow ROIs as spherical onions starting from set centroids and to see, how consistency
decreases with the increasing radius.
"""
import numpy as np

import functions
import onion_parameters as params

subjectFolders = params.testSubjectFolders
nSubjects = len(subjectFolders)
roiInfoFile = 'atlases/brainnetome/brainnetome_thr25_4mm_rois.mat'
resolution = 4
distanceMatrixPath = params.distanceMatrixPath

radiusRange = np.arange(1,20,2)
nRadia = len(radiusRange)

centroids, _, voxelCoords = functions.readROICentroids(roiInfoFile,readVoxels=True)
nROIs = len(centroids)

consistencies = np.zeros((nSubjects,nRadia,nROIs))

for j, radius in enumerate(radiusRange):
    if j == 0: # calculating the distance matrix only once
        roiInfo = functions.defineSphericalROIs(centroids,voxelCoords,radius,resolution=resolution,save=True,savePath=distanceMatrixPath)
    else:
        roiInfo = functions.defineSphericalROIs(centroids,voxelCoords,radius,resolution=resolution,distanceMatrixPath=distanceMatrixPath)
    voxelIndices = roiInfo['ROIVoxels']
    for i, subject in enumerate(subjectFolders):
        voxelTsFilePath = subject + params.ROIVoxelTsFileName
        for k, ROI in enumerate(voxelIndices):
            consistency = functions.calculateSpatialConsistency(ROI,voxelTsFilePath)
            consistencies[i,j,k] = consistency

meanConsistencies = np.mean(consistencies,axis=(0,2)) # average & std over subjects and ROIs
stdConsistencies = np.std(consistencies,axis=(0,2))
            
            
            
        


