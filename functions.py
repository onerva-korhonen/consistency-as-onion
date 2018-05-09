# -*- coding: utf-8 -*-
"""
Created on Tue May  8 13:27:05 2018

@author: aokorhon

Functions for analysing the relationship between ROI's size and spatial consistency. Call either from
a frontend script or interactively.
"""
import numpy as np
import cPickle as pickle

from scipy import io


def readROICentroids(ROIInfoFile, readVoxels=False):
    """
    Reads ROI data, in particular the coordinates of ROI centroids. An option
    for reading coordinates of all voxels exist
    
    Parameters:
    -----------
        ROIInfoFile: str, path to the file that contains information about ROIs
        readVoxels: bool, will voxel coordinates be read? (default=False)
        
    Returns:
    --------
        ROICentroids: nROIs x 3 np.array, coordinates of the centroids of the ROIs
        voxelCoordinates: nVoxels x 3 np.array, coordinates of all voxels
    """
    infoData = io.loadmat(ROIInfoFile)
    ROIInfo = infoData['rois'][0]
    nROIs = ROIInfo.shape[0] # number of ROIs
    ROIcentroids = np.zeros((nROIs,3))
    voxelCoordinates = []
    for i, ROI in enumerate(ROIInfo):
        ROIcentroids[i,:] = ROI['centroid'][0]
        if readVoxels:
            voxelCoordinates.extend(list(ROI['map']))
    voxelCoordinates = np.array(voxelCoordinates)
    return ROIcentroids, voxelCoordinates
    
def getDistanceMatrix(ROICentroids, voxelCoordinates, save=False, savePath=''):
    """
    Calculates the centroid-to-voxel Euclidean distance between all ROI-voxel
    pairs. If the resolution (and number of voxels) is high, consider setting
    save=True in order to avoid unnecessarely repeating timely calculations.
    
    Parameters:
    -----------
    ROICentroids: nROIs x 3 np.array, coordinates of the centroids of the ROIs
    voxelCoordinates: nVoxels x 3 np.array, coordinates of all voxels
    save: bool, will the distance matrix be saved in a file? (default=False)
    savePath: str, path for saving the distance matrix (default='')
        
    Returns:
    --------
    distanceMatrix: nROIs x nVoxels np.array, matrix of distances of each voxel from each ROI centroid
    """
    nROIs = ROICentroids.shape[0]
    nVoxels = voxelCoordinates.shape[0]
    distanceMatrix = np.zeros((nROIs, nVoxels))
    for i, centroid in enumerate(ROICentroids):
        distanceMatrix[i,:] = np.sqrt(np.sum((voxelCoordinates-centroid)**2,axis=1))
    if save:
        distanceData = {}
        distanceData['distence_matrix'] = distanceMatrix
        with open(savePath, 'wb') as f:
            pickle.dump(distanceData, f, -1)
    return distanceMatrix
    
def calculateSpatialConsistency(voxelIndices, voxelTsFilePath, type='pearson c', fTransform=False):
    """
    Calculates the spatial consistency of a chunk of voxels. By default,
    spatial consistency is defined as the mean Pearson correlation coefficient
    between the voxel time series.
    
    Parameters:
    -----------
    voxelIdices: np.array, indices of voxels; these indices should refer to voxels' 
            locations in the file containing voxel time series; nete that the chunk
            must contain more than one voxel
    voxelTsFilePath: str, path to a file that contains the voxel time series;
            the file should contain a dictionary with a key 'roi_voxel_data' (and
            possible additional keys), value assigned to this key is a structured
            np.array with a field name 'roi_voxel_ts' (and possible additional 
            fields), this field contains voxel time series
    type: str, definition of spatial consistency to be used; default:
          'pearson c' (mean Pearson correlation coefficient), other options:
          TODO: add other consistency options
    fTransform: bool, are the correlations Fisher f transformed before averaging
                when type = 'pearson c' (default=False)
          
    Returns:
    --------
    spatialConsistency: dbl, value of spatial consistency
    """
    try:
        if np.amax(voxelIndices.shape) == 1:
            raise ValueError('voxelIndices must contain  more than one element')
        
        allVoxelTs = io.loadmat(voxelTsFilePath)['roi_voxel_data'][0]['roi_voxel_ts'][0]
        voxelTs = allVoxelTs[voxelIndices,:]
        if type == 'pearson c':
            correlations = np.corrcoef(voxelTs)
            correlations = correlations[np.tril_indices(voxelTs.shape[0],k=-1)] # keeping only the lower triangle, diagonal is discarded
            if fTransform:
                correlations = np.arctanh(correlations)
                spatialConsistency = np.tanh(np.mean(correlations))
            else:
                spatialConsistency = np.mean(correlations)
        else:
            spatialConsistency = 0
        return spatialConsistency
    except ValueError as error:
        print(error.message)
        
        
    
    
