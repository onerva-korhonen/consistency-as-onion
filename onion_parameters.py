# -*- coding: utf-8 -*-
"""
Created on Tue May  8 13:40:42 2018

@author: aokorhon

A storage file for the parameters used in the "Consistency as onion" project.
"""
# Data locations and other file paths:

testSubjectFolders = ['/media/onerva/KINGSTON/test-data/010/',
                      '/media/onerva/KINGSTON/test-data/011/']

# Input file names:

ROIVoxelTsFileName = 'roi_voxel_ts_all_rois4mm_FWHM0.mat'
ROIVoxelTsInfoFileName = 'roi_voxel_ts_all_rois4mm_WHM0_info.mat'

originalROIInfoFile = '/media/onerva/KINGSTON/test-data/group_roi_mask-30-4mm_with_subcortl_and_cerebellum.mat'
optimizedROIInfoFile = 'optimized-rois-test.nii' # note: this is not a full path; this is the name of the file that is located in each subject's folder (as optimization is done separately for each subject)

# Paths for saving

distanceMatrixPath = '/media/onerva/KINGSTON/test-data/outcome/distance-matrix-brainnetome-4mm.pkl'
consistencyVsRadiusPath = '/media/onerva/KINGSTON/test-data/outcome/consistency-over-radia-brainnetome-4mm-test.pkl'

spatialConsistencyValidationPath = '/media/onerva/KINGSTON/test-data/outcome/spatial-consisistency-validation.pdf'

# Parallelization
nCPUs = 5

# Distributions and binning
nConsistencyBins = 50

# Visualization
optimizedColor = 'r'
optimizedAlpha = 0.9
originalColor = 'k'
originalAlpha = 0.5

