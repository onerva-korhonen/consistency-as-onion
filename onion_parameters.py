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
optimizedROIInfoFile = 'optimized-rois-test-spatialConsistency.pkl' # note: this is not a full path; this is the name of the file that is located in each subject's folder (as optimization is done separately for each subject)

# Paths for saving

originalSpatialConsistencySavePath = '/media/onerva/KINGSTON/test-data/spatial-consistency-original.pkl'
optimizedSpatialConsistencySaveName = 'spatial-consistency-optimized-spatialConsistency.pkl'
originalCorrelationSavePath = '/media/onerva/KINGSTON/test-data/correlation-to-centroid-original.pkl'
optimizedCorrelationSaveName = 'correlation-to-centroid-optimized-spatialConsistency.pkl'

distanceMatrixPath = '/media/onerva/KINGSTON/test-data/outcome/distance-matrix-brainnetome-4mm.pkl'
consistencyVsRadiusPath = '/media/onerva/KINGSTON/test-data/outcome/consistency-over-radia-brainnetome-4mm-test.pkl'

spatialConsistencyValidationPath = '/media/onerva/KINGSTON/test-data/outcome/spatial-consisistency-validation-weighted-mean-consistency.pdf'
sizeSavePath = '/media/onerva/KINGSTON/test-data/outcome/spatial-consistency-validation-local-weighted-consistency-sizes.pdf'
maximalMeasureFigurePath = '/media/onerva/KINGSTON/test-data/outcome/maximal-measure-spatialConsistency.pdf'

# Parallelization
nCPUs = 5

# Distributions and binning
nConsistencyBins = 50
nSizeBins = 50

# Visualization
optimizedColor = 'r'
optimizedAlpha = 0.9
originalColor = 'k'
originalAlpha = 0.5
inROILs = '-'
betweenROILs = '--'
