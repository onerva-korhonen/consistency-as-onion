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


def readROICentroids(ROIInfoFile, readVoxels=False, fixCentroids=False):
    """
    Reads ROI data, in particular the coordinates of ROI centroids. An option
    for reading coordinates of all voxels exist. 
    
    Parameters:
    -----------
    ROIInfoFile: str, path to the file that contains information about ROIs
    readVoxels: bool, will voxel coordinates be read? (default=False)
    fixCentroids: bool, if fixCentroids=True, the function will return as centroid 
                  the one of ROI's voxels that is closest to the actual centroid.
        
    Returns:
    --------
    ROICentroids: nROIs x 3 np.array, coordinates (in voxels) of the centroids of the ROIs
    ROIMNICentroids: nROIs x 3 np.array, coordinates (in mm) of the centroids of the ROIs
    voxelCoordinates: nVoxels x 3 np.array, coordinates (in voxels) of all voxels
    """
    infoData = io.loadmat(ROIInfoFile)
    ROIInfo = infoData['rois'][0]
    nROIs = ROIInfo.shape[0] # number of ROIs
    ROICentroids = np.zeros((nROIs,3))
    ROIMNICentroids = np.zeros((nROIs,3))
    voxelCoordinates = []
    for i, ROI in enumerate(ROIInfo):
        centroid = ROI['centroid'][0]
        if fixCentroids:
            ROIMap = ROI['map']
            distances = np.zeros(ROIMap.shape[0])
            for j, voxel in enumerate(ROIMap):
                distances[j] = np.sqrt(np.sum((voxel-centroid)**2,axis=1))
            centroid = ROIMap[np.where(distances==np.amin(distances))]
        ROICentroids[i,:] = centroid
        ROIMNICentroids[i,:] = ROI['centroidMNI'][0]
        if readVoxels:
            voxelCoordinates.extend(list(ROI['map']))
    voxelCoordinates = np.array(voxelCoordinates)
    return ROICentroids, ROIMNICentroids, voxelCoordinates
    
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
    if len(ROICentroids.shape) == 1:
        nROIs = 1
    else:
        nROIs = ROICentroids.shape[0]
    nVoxels = voxelCoords.shape[0]
    distanceMatrix = np.zeros((nROIs, nVoxels))
    if nROIs == 1:
        distanceMatrix[0,:] = np.sqrt(np.sum((voxelCoords-ROICentroids)**2,axis=1))
    else:
        for i, centroid in enumerate(ROICentroids):
            distanceMatrix[i,:] = np.sqrt(np.sum((voxelCoords-centroid)**2,axis=1))
    if save:
        distanceData = {}
        distanceData['distanceMatrix'] = distanceMatrix
        with open(savePath, 'wb') as f:
            pickle.dump(distanceData, f, -1)
    return distanceMatrix
    
def calculateSpatialConsistency(voxelIndices, allVoxelTs, type='pearson c', fTransform=False):
    """
    Calculates the spatial consistency of a chunk of voxels. By default,
    spatial consistency is defined as the mean Pearson correlation coefficient
    between the voxel time series.
    
    Parameters:
    -----------
    voxelIdices: np.array, indices of voxels; these indices should refer to voxels' 
            locations in the file containing voxel time series; note that the chunk
            must contain at least one voxel
    allVoxelTs: structured np.array with a field name 'roi_voxel_ts' (and possible additional 
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
# NOTE: What follows is ugly... At the moment, I'll accept empty voxel groups (and
# return consistency = 0). However, it's probably that the try-except structure will be reactivated
# later so let's keep it here as a comment...
#    try: 
#        if np.amax(voxelIndices.shape) == 0:
#            raise ValueError('voxelIndices is empty, cannot calculate consistency')
#            
#        if np.amax(voxelIndices.shape) == 1:
#            spatialConsistency = 1. # a single voxel is always fully consistent
#        else:
#            allVoxelTs = io.loadmat(voxelTsFilePath)['roi_voxel_data'][0]['roi_voxel_ts'][0]
#            voxelTs = allVoxelTs[voxelIndices,:]
#            if type == 'pearson c':
#                correlations = np.corrcoef(voxelTs)
#                correlations = correlations[np.tril_indices(voxelTs.shape[0],k=-1)] # keeping only the lower triangle, diagonal is discarded
#                if fTransform:
#                    correlations = np.arctanh(correlations)
#                    spatialConsistency = np.tanh(np.mean(correlations))
#                else:
#                    spatialConsistency = np.mean(correlations)
#            else:
#                spatialConsistency = 0
#        return spatialConsistency
#    except ValueError as error:
#        print error.message
    if np.amax(voxelIndices.shape) == 0:
        spatialConsistency = 0
        print "Detected an empty ROI, set consistency to 0."
    elif np.amax(voxelIndices.shape) == 1:
        spatialConsistency = 1. # a single voxel is always fully consistent
    else: 
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
    if np.isnan(spatialConsistency):
        print 'nan detected!'
    return spatialConsistency
        
def defineSphericalROIs(ROICentroids, voxelCoords, radius, resolution=4.0, names='', distanceMatrixPath='', save=False, savePath=''):
    """
    Constructs a set of (approximately) spherical ROIs with a given radius located
    around a given set of centroid points.
    
    Parameters:
    -----------
    ROICentroids: nROIs x 3 np.array, coordinates of the centroids of the ROIs.
                  This can be a ROI centroid from an atlas but also any other
                  (arbitrary) point.
    voxelCoords: nVoxels x 3 np.array, coordinates of all voxels. Must be
                 in the same space with ROICentoids. Coordinates must be ordered by ROIs:
                 first all voxels of ROI 1, then all voxels of ROI 2, etc.
    radius: dbl, radius of the sphere measured by number of voxels. A sphere with
            radius = 1 contains the centroid and it's six closest neighbors. If 
            radius < 1, the sphere will contain only the centroid.
    resolution: dbl, length of the voxel edge in mm. Resolution gives the length
                of radius in mm. Default = 4.
    names: list of strs, names of the ROIs, can be e.g. the anatomical name associated with
           the centroid. Default = ''.
    distanceMatrixPath: str, path of an existing distance matrix. If distanceMatrixPath is given,
                        the existing matrix will be used instead of calculating a new one.
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
    if len(ROICentroids.shape) == 1:
        nROIs = 1
    else:
        nROIs = ROICentroids.shape[0]    
    
    if distanceMatrixPath == '':
        distanceMatrix = getDistanceMatrix(ROICentroids, voxelCoords, save=save, savePath=savePath)
    else:
         f = open(distanceMatrixPath, "rb")
         data = pickle.load(f)
         f.close()
         distanceMatrix = data['distanceMatrix']
    
    ROIMaps = []
    ROIVoxels = []
    ROISizes = np.zeros((nROIs, 1))
    
    for i, (distances, centroid) in enumerate(zip(distanceMatrix, ROICentroids)):
        if radius < 1: # In this case, the sphere contains only the centroid
            voxelIndices = np.where(distances == 0)
            ROIMap = centroid
            ROISize = 1
        else:
            voxelIndices = np.where(distances <= radius)[0]
            ROIMap = voxelCoords[voxelIndices,:]
            ROISize = np.amax(voxelIndices.shape)
        
        ROIMaps.append(ROIMap)
        ROIVoxels.append(voxelIndices)
        ROISizes[i] = ROISize
        
    sphericalROIs = {'ROICentroids':ROICentroids,'ROIMaps':ROIMaps,'ROIVoxels':ROIVoxels,'ROISizes':ROISizes,'ROINames':names}
    
    return sphericalROIs
    
def findNeighbors(voxelCoords, resolution=1, allVoxels=[]):
    """
    Returns the 6 closest neighbors (the ones sharing a face) of a voxel.
    
    Parameters:
    -----------
    voxelCoords: 1x3 np.array, coordinates of a voxel (either in voxels or in mm)
    resolution: double, distance between voxels if coordinates are given in mm;
                if coordinates are given in voxels, use the default value 1 (voxels
                are 1 voxel away from each other).
    allVoxels: iterable, coordinates of all acceptable voxels. If allVoxels is given,
               only neighbors in allVoxels are returned (default: []).
                
    Returns:
    --------    
    neighbors: 6x6 np.array, coordinates of the closest neighbors of the voxel
    """
    x = voxelCoords[0]
    y = voxelCoords[1]
    z = voxelCoords[2]    
    
    neighbors = np.array([[x+resolution,y,z],
                         [x-resolution,y,z],
                         [x,y+resolution,z],
                         [x,y-resolution,z],
                         [x,y,z+resolution],
                         [x,y,z-resolution]])
                         
    if not len(allVoxels) == 0:
        accepted = np.zeros(neighbors.shape[0])
        for i, neighbor in enumerate(neighbors):
            accepted[i] = np.any((np.array(allVoxels) == neighbor).all(axis=1))
        neighbors = neighbors[np.where(accepted)]
                         
    return neighbors
    
def findROIlessVoxels(voxelCoordinates,ROIInfo):
    """
    Returns the indices and coordinates of voxels that do not belong to any ROI.
    
    Parameters:
    -----------
    voxelCoordinates: nVoxels x 3 np.array, coordinates (in voxels) of all voxels
    ROIInfo: dic, contains:
             ROIMaps: list of ROISizes x 3 np.arrays, coordinates of voxels
                           belonging to each ROI. len(ROIMaps) = nROIs.
             ROIVoxels: list of ROISizes x 1 np.array, indices of the voxels belonging
                     to each ROI. These indices refer to the columns of
                     the distance matrix, as well as to the rows of the
                     voxel time series file. len(ROIVoxels) = nROIs.
    Returns:
    --------
    ROIlessVoxels: dic, contains:
                   ROIlessIndices: NROIless x 1 np.array, indices of ROIless voxels in voxelCoordinates
                   ROIlessMap: NROIless x 3 np.array, coordinates of the ROIless voxels
    """
    ROIMaps = ROIInfo['ROIMaps']
    for i, ROI in enumerate(ROIMaps):
        if i == 0:
            inROIVoxels = ROI
        else:
            inROIVoxels = np.concatenate((inROIVoxels,ROI),axis=0)
    ROIlessIndices = []
    ROIlessMap = []
    for i, voxel in enumerate(voxelCoordinates):
        if not np.any((inROIVoxels == voxel).all(axis=1)): # voxel is not found in any ROI map
            ROIlessIndices.append(i)
            ROIlessMap.append(voxel)
    ROIlessIndices = np.array(ROIlessIndices)
    ROIlessMap = np.array(ROIlessMap)
    ROIlessVoxels = {'ROIlessIndices':ROIlessIndices,'ROIlessMap':ROIlessMap}
    return ROIlessVoxels
    
def findROIlessNeighbors(ROIIndex,voxelCoordinates,ROIInfo):
    """
    Finds the neighboring voxels of a ROI that do not belong to any ROI.
    
    Parameters:
    -----------
    ROIIndex: int, index of the ROI in the lists of ROIInfo (see below)
    voxelCoordinates: nVoxels x 3 np.array, coordinates (in voxels) of all voxels
    ROIInfo: dic, contains:
             ROIMaps: list of ROISizes x 3 np.arrays, coordinates of voxels
                           belonging to each ROI. len(ROIMaps) = nROIs.
             ROIVoxels: list of ROISizes x 1 np.array, indices of the voxels belonging
                     to each ROI. These indices refer to the columns of
                     the distance matrix, as well as to the rows of the
                     voxel time series file. len(ROIVoxels) = nROIs.
    Returns:
    --------
    ROIlessNeighbors: dic, contains:
                      ROIlessIndices: list, indices of ROIless neighbor voxels in voxelCoordinates
                      ROIlessMap: NNeighbors x 3 np.array, coordinates of ROIless neighbor voxels
    """
    ROIMap = ROIInfo['ROIMaps'][ROIIndex]
    for i, voxel in enumerate(ROIMap):
       neighbors = findNeighbors(voxel,allVoxels=voxelCoordinates)
       if i == 0:
           ROINeighbors = neighbors
       else:
           ROINeighbors = np.concatenate((ROINeighbors,neighbors),axis=0)
    ROINeighbors = np.unique(ROINeighbors,axis=0) # removing dublicates
    ROIlessNeighbors = findROIlessVoxels(ROINeighbors,ROIInfo)
    return ROIlessNeighbors
    
    


   
        
    
        
    
  
        
    
    
