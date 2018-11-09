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
import matplotlib.pylab as plt
from scipy import io
import cPickle as pickle

import sys
try:
    sys.path.insert(0,'/home/onerva/consistency-as-onion')
except:
    sys.path.insert(0,'/home/onerva/projects/consistency-as-onion')   

import functions
import onion_parameters as params

visualizeOnly = False
consistencySavePath = params.consistencyVsRadiusPath
figureSavePath = '/media/onerva/KINGSTON/test-data/outcome/consistency-vs-radius.pdf'

if not visualizeOnly:
    subjectFolders = params.testSubjectFolders
    nSubjects = len(subjectFolders)
    roiInfoFile = '/media/onerva/KINGSTON/test-data/group_roi_mask-30-4mm_with_subcortl_and_cerebellum.mat'
    resolution = 4
    distanceMatrixPath = params.distanceMatrixPath
    
    radiusRange = np.arange(1,20,2)
    nRadia = len(radiusRange)
    
    centroids, _, voxelCoords, _ = functions.readROICentroids(roiInfoFile,readVoxels=True)
    nROIs = len(centroids)
    
    consistencies = np.zeros((nSubjects,nRadia,nROIs))
    
    allVoxelTs = []
    
    
    for i, subject in enumerate(subjectFolders): # reading the data here and not inside the loop to save time
        voxelTsFilePath = subject + params.ROIVoxelTsFileName
        allVoxelTs.append(io.loadmat(voxelTsFilePath)['roi_voxel_data'][0]['roi_voxel_ts'][0])
    
    for j, radius in enumerate(radiusRange):
        if j == 0: # calculating the distance matrix only once
            roiInfo = functions.defineSphericalROIs(centroids,voxelCoords,radius,resolution=resolution,save=True,savePath=distanceMatrixPath)
        else:
            roiInfo = functions.defineSphericalROIs(centroids,voxelCoords,radius,resolution=resolution,distanceMatrixPath=distanceMatrixPath)
        voxelIndices = roiInfo['ROIVoxels']
        for i, subject in enumerate(allVoxelTs):
            for k, ROI in enumerate(voxelIndices):
                cfg = {'allROITs':subject}
                params = (cfg,ROI)
                consistency = functions.calculateSpatialConsistency(params)
                consistencies[i,j,k] = consistency
                
    consistencyData = {'radia':radiusRange,'subjects':subjectFolders,'consistencies':consistencies}   
    with open(consistencySavePath, 'wb') as f:
        pickle.dump(consistencyData, f, -1)
      
else:
    f = open(consistencySavePath, "rb")
    consistencyData = pickle.load(f)
    f.close()
    consistencies = consistencyData['consistencies']
    radiusRange = consistencyData['radia']
    
meanConsistencies = np.mean(consistencies,axis=(0,2)) # average & std over subjects and ROIs
stdConsistencies = np.std(consistencies,axis=(0,2))
    
fig = plt.figure()
ax = fig.add_subplot(111)
#plt.plot(radiusRange,meanConsistencies,ls='',marker='.',color='k')
plt.errorbar(radiusRange,meanConsistencies,yerr=stdConsistencies,ls='',marker='.',color='k',ecolor='k')
plt.plot([radiusRange[0],radiusRange[-1]],[0,0],ls='--',color='k',alpha=0.7)
ax.set_xlabel('Radius (in voxels)')
ax.set_ylabel('Mean spatial consistency')
plt.tight_layout()
plt.savefig(figureSavePath,format='pdf',bbox_inches='tight')


            
            
            
        


