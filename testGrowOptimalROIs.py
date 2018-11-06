# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 12:17:59 2018

@author: onerva

A script for testing functions.growOprimalROIs
"""
import functions
import onion_parameters as params

import numpy as np
from scipy import io


allVoxelTsPath = params.testSubjectFolders[0] + params.ROIVoxelTsFileName
ROIInfoFile = params.originalROIInfoFile
savePath = '/media/onerva/KINGSTON/test-data/optimized-rois-test.nii'

ROICentroids,_,voxelCoordinates,_ = functions.readROICentroids(ROIInfoFile,readVoxels=True,fixCentroids=True)
allVoxelTs = io.loadmat(allVoxelTsPath)['roi_voxel_data'][0]['roi_voxel_ts'][0]

#voxelCoordinates = voxelCoordinates[0:1000,:]

cfg = {}
cfg['ROICentroids'] = ROICentroids
cfg['voxelCoordinates'] = voxelCoordinates
cfg['names'] = ''
cfg['allVoxelTs'] = allVoxelTs

ROIInfo = functions.growOptimizedROIs(cfg)

functions.createNii(ROIInfo, savePath, imgSize=[45,54,45], affine=np.eye(4))

