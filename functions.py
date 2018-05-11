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
    
def getDistanceMatrix(ROICentroids, voxelCoords, save=False, savePath=''):
    """
    Calculates the centroid-to-voxel Euclidean distance between all ROI-voxel
    pairs. If the resolution (and number of voxels) is high, consider setting
    save=True in order to avoid unnecessarely repeating timely calculations.
    
    Parameters:
    -----------
    ROICentroids: nROIs x 3 np.array, coordinates of the centroids of the ROIs.
                  This can be a ROI centroid from an atlas but also any other
                  (arbitrary) point.
    voxelCoords: nVoxels x 3 np.array, coordinates of all voxels
    save: bool, will the distance matrix be saved in a file? (default=False)
    savePath: str, path for saving the distance matrix (default='')
        
    Returns:
    --------
    distanceMatrix: nROIs x nVoxels np.array, matrix of distances of each voxel from each ROI centroid
    """
    nROIs = ROICentroids.shape[0]
    nVoxels = voxelCoords.shape[0]
    distanceMatrix = np.zeros((nROIs, nVoxels))
    for i, centroid in enumerate(ROICentroids):
        distanceMatrix[i,:] = np.sqrt(np.sum((voxelCoords-centroid)**2,axis=1))
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
        if np.amax(voxelIndices.shape) == 0:
            raise ValueError('voxelIndices is empty, cannot calculate consistency')
            
        if np.amax(voxelIndices.shape) == 1:
            spatialConsistency = 1. # a single voxel is always fully consistent
        else:
        
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
        print error.message
        
def defineSphericalROIs(ROICentroids, voxelCoords, radius, resolution=4.0, names='', save=False, savePath=''):
    """
    Constructs a set of (approximately) spherical ROIs with a given radius located
    around a given set of centroid points.
    
    Parameters:
    -----------
    ROICentroids: nROIs x 3 np.array, coordinates of the centroids of the ROIs.
                  This can be a ROI centroid from an atlas but also any other
                  (arbitrary) point.
    voxelCoords: nVoxels x 3 np.array, coordinates of all voxels. Must be
                      in the same space with ROICentoids
    radius: dbl, radius of the sphere measured by number of voxels. A sphere with
            radius = 1 contains the centroid and it's six closest neighbors. If 
            radius < 1, the sphere will contain only the centroid.
    resolution: dbl, length of the voxel edge in mm. Resolution gives the length
                of radius in mm. Default = 4.
    names: list of strs, names of the ROIs, can be e.g. the anatomical name associated with
           the centroid. Default = ''.
    save: bool, will the distance matrix between ROI centroids and voxels
          be saved in a file? (default=False). Consider saving if resolution
          is very high and number of voxels very large.
    savePath: str, path for saving the distance matrix (default='')
                
    Returns:
    --------
    sphericalROIs: dic, contains:
                  ROICentroids: nROIs x 3 np.array, coordinates of the centroids
                  ROIMaps: list of ROISizes x 3 np.arrays, coordinates of voxels
                           belonging to each ROI. len(ROIMaps) = nROIs.
                  ROIVoxels: list of ROISizes x 1 np.array, indices of the voxels belonging
                             to each ROI. These indices refer to the columns of
                             the distance matrix, as well as to the rows of the
                             voxel time series file. len(ROIVoxels) = nROIs.
                  ROISizes: nROIs x 1 np.array of ints. Sizes of ROIs defined as
                            number of voxels in the sphere
                  ROINames: list of strs, name of the spherical ROI. If no name is given as
                            input parameter for a ROI, this is set to ''. 
                            len(ROINames) = NROIs.
    """
    nROIs = np.amax(ROICentroids.shape)    
    
    distanceMatrix = getDistanceMatrix(ROICentroids, voxelCoords, save=save, savePath=savePath)
    physicalRadius = radius*resolution
    
    ROIMaps = []
    ROIVoxels = []
    ROISizes = np.zeros((nROIs, 1))
    
    for i, (distances, centroid) in enumerate(zip(distanceMatrix, ROICentroids)):
        if radius < 1: # In this case, the sphere contains only the centroid
            voxelIndices = np.where(distances == 0)
            ROIMap = centroid
            ROISize = 1
        else:
            voxelIndices = np.where(distances <= physicalRadius)[0]
            ROIMap = voxelCoords[voxelIndices,:]
            ROISize = np.amax(voxelIndices.shape)
        
        ROIMaps.append(ROIMap)
        ROIVoxels.append(voxelIndices)
        ROISizes[i] = ROISize
        
    sphericalROIs = {'ROICentroids':ROICentroids,'ROIMaps':ROIMaps,'ROIVoxels':ROIVoxels,'ROISizes':ROISizes,'ROINames':names}
    
    return sphericalROIs
        
        
    
        
    
  
        
    
    
