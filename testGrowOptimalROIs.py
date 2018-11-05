# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 12:17:59 2018

@author: onerva

A script for testing functions.growOprimalROIs
"""
import functions
import onion_parameters as params

from scipy import io

allVoxelTsPath = params.testSubjectFolders[0] + params.ROIVoxelTsFileName
ROIInfoFile = params.ROIInfoFile

ROICentroids,_,voxelCoordinates,_ = functions.readROICentroids(ROIInfoFile,readVoxels=True,fixCentroids=True)
allVoxelTs = io.loadmat(allVoxelTsPath)['roi_voxel_data'][0]['roi_voxel_ts'][0]

cfg = {}
cfg['ROICentroids'] = ROICentroids
cfg['voxelCoordinates'] = voxelCoordinates
cfg['names'] = ''
cfg['allVoxelTs'] = allVoxelTs

ROIInfo = functions.growOptimizedROIs(cfg)


