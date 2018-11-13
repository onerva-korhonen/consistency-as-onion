# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 12:17:59 2018

@author: onerva

A script for testing functions.growOprimalROIs
"""
import numpy as np
from scipy import io
import cPickle as pickle
import matplotlib.pylab as plt

import os.path
import sys
if os.path.exists('/home/onerva/consistency-as-onion'):
    sys.path.insert(0,'/home/onerva/consistency-as-onion')
else:
    sys.path.insert(0,'/home/onerva/projects/consistency-as-onion')   

import functions
import onion_parameters as params

allVoxelTsPath = params.testSubjectFolders[0] + params.ROIVoxelTsFileName
ROIInfoFile = params.originalROIInfoFile
niiSavePath = params.testSubjectFolders[0] + '/optimized-rois-test-spatialConsistency.nii'
pickleSavePath = params.testSubjectFolders[0] + '/optimized-rois-test-spatialConsistency.pkl'

ROICentroids,_,voxelCoordinates,_ = functions.readROICentroids(ROIInfoFile,readVoxels=True,fixCentroids=True)
allVoxelTs = io.loadmat(allVoxelTsPath)['roi_voxel_data'][0]['roi_voxel_ts'][0]

cfg = {}
cfg['ROICentroids'] = ROICentroids
cfg['voxelCoordinates'] = voxelCoordinates
cfg['names'] = ''
cfg['allVoxelTs'] = allVoxelTs

cfg['threshold'] = -1.0
cfg['targetFunction'] = 'spatialConsistency'

ROIInfo, selectedMeasures = functions.growOptimizedROIs(cfg)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(selectedMeasures)
ax.set_xlabel('Iteration step')
ax.set_ylabel('Maximal similarity index')
plt.tight_layout()
plt.savefig(params.maximalMeasureFigurePath,format='pdf',bbox_inches='tight')


functions.createNii(ROIInfo, niiSavePath, imgSize=[45,54,45], affine=np.eye(4))
with open(pickleSavePath, 'wb') as f:
        pickle.dump(ROIInfo, f, -1)

