# -*- coding: utf-8 -*-
"""
Created on Tue May  8 13:27:05 2018

@author: aokorhon

Functions for analysing the relationship between ROI's size and spatial consistency. Call either from
a frontend script or interactively.
"""
import numpy as np
import cPickle as pickle
import string

from scipy import io
from scipy.stats import binned_statistic
from concurrent.futures import ProcessPoolExecutor as Pool

import pymnet # the multilayer module by Mikko Kivelä


def readROICentroids(ROIInfoFile, readVoxels=False, fixCentroids=False, resolution=4):
    """
    Reads ROI data, in particular the coordinates of ROI centroids. An option
    for reading coordinates of all voxels exist. 
    
    Parameters:
    -----------
    ROIInfoFile: str, path to the file that contains information about ROIs
    readVoxels: bool, will voxel coordinates be read? (default=False)
    fixCentroids: bool, if fixCentroids=True, the function will return as centroid 
                  the one of ROI's voxels that is closest to the actual centroid.
    resolution: int, the physical distance (in mm) between two points in the voxel space
                (used for transforming centroids into the MNI space if the ROIInfoFile
                doesn't contain the MNI coordinates; default: 4).
        
    Returns:
    --------
    ROICentroids: nROIs x 3 np.array, coordinates (in voxels) of the centroids of the ROIs
    ROIMNICentroids: nROIs x 3 np.array, coordinates (in mm) of the centroids of the ROIs
    voxelCoordinates: nVoxels x 3 np.array, coordinates (in voxels) of all voxels
    ROIMaps: list of ROISize x 3 np.arrays, coordinates (in voxels) of voxels belonging to each ROI
    """
    if ROIInfoFile[-4::] == '.mat':
        infoData = io.loadmat(ROIInfoFile)
        ROIInfo = infoData['rois'][0]
        nROIs = ROIInfo.shape[0] # number of ROIs
        ROICentroids = np.zeros((nROIs,3),dtype=int)
        ROIMNICentroids = np.zeros((nROIs,3))
        voxelCoordinates = []
        ROIMaps = []
        for i, ROI in enumerate(ROIInfo):
            centroid = np.array(ROI['centroid'][0]) - np.array([1,1,1]) # correcting for the indexing difference between Matlab and Spyder
            if fixCentroids:
                ROIMap = ROI['map'] - np.ones(ROI['map'].shape,dtype=int)
                distances = np.zeros(ROIMap.shape[0])
                for j, voxel in enumerate(ROIMap):
                    distances[j] = np.sqrt(np.sum((voxel-centroid)**2))
                centroid = ROIMap[np.where(distances==np.amin(distances))[0][0]] # if multiple voxels are at the same distance from the centroid, the first one is picked
            ROICentroids[i,:] = centroid
            ROIMNICentroids[i,:] = ROI['centroidMNI'][0]
            if readVoxels:
                voxelCoordinates.extend(list(ROI['map'] - np.ones(ROI['map'].shape,dtype=int)))
                ROIMaps.append(ROI['map']-np.ones(ROI['map'].shape,dtype=int))
        voxelCoordinates = np.array(voxelCoordinates)
    # TODO: check the reading of pickle files when the file structure is more established    
    elif ROIInfoFile[-4::] == '.pkl':
        f = open(ROIInfoFile, "rb")
        ROIInfo = pickle.load(f)
        f.close()
        ROIMaps = ROIInfo['ROIMaps']
        ROICentroids = np.zeros((len(ROIMaps),3),dtype=int)
        ROIMNICentroids = np.zeros((len(ROIMaps),3),dtype=int)
        for i, ROIMap in enumerate(ROIMaps):
            centroid = getCentroid(ROIMap)
            ROICentroids[i,:] = centroid
            ROIMNICentroids[i,:] = spaceToMNI(centroid,resolution)
        voxelCoordinates = []
        if readVoxels:
            for ROIMap in ROIMaps:
                if len(ROIMap.shape) == 1:
                    voxelCoordinates.append(ROIMap)
                else:
                    voxelCoordinates.extend(ROIMap)
            voxelCoordinates = np.array(voxelCoordinates)
    else:
        print 'Unknown file format; accepted formats are .mat and .pkl'
        
    return ROICentroids, ROIMNICentroids, voxelCoordinates, ROIMaps
    
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
    
def calculateSpatialConsistency(params):
    """
    Calculates the spatial consistency of a chunk of voxels. By default,
    spatial consistency is defined as the mean Pearson correlation coefficient
    between the voxel time series.
    
    Parameters:
    -----------
    params: tuple, containing:
    
        cfg: dict, containing:
            allVoxelTs: structured np.array with a field name 'roi_voxel_ts' (and possible additional 
                    fields), this field contains voxel time series
            consistencyType: str, definition of spatial consistency to be used; default:
                  'pearson c' (mean Pearson correlation coefficient), other options:
                  TODO: add other consistency options
            fTransform: bool, are the correlations Fisher f transformed before averaging
                        when consistencyType = 'pearson c' (default=False)
                
        voxelIndices: np.array, indices of voxels; these indices should refer to voxels' 
                locations in the file containing voxel time series; note that the chunk
                must contain at least one voxel
          
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
#            if consistencyType == 'pearson c':
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
    cfg = params[0]
    allVoxelTs = cfg['allVoxelTs']
    consistencyType = cfg['consistencyType']
    fTransform = cfg['fTransform']
    voxelIndices = params[1]
    
    #print 'Calculating spatial consistency'
    
    if np.amax(voxelIndices.shape) == 0:
        spatialConsistency = 0
        print "Detected an empty ROI, set consistency to 0."
    elif np.amax(voxelIndices.shape) == 1:
        spatialConsistency = 1. # a single voxel is always fully consistent
    else: 
        voxelTs = allVoxelTs[voxelIndices,:]
        if consistencyType == 'pearson c':
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
    
def calculateSpatialConsistencyInParallel(voxelIndices,allVoxelTs,consistencyType='pearson c', fTransform=False, nCPUs=5):
    """
    A wrapper function for calculating the spatial consistency in parallel across
    ROIs.
    
    Parameters:
    -----------
    voxelIndices: list of np.arrays, each array containing indices of voxels of one ROI; 
               these indices should refer to voxels' 
               locations in the file containing voxel time series; note that the chunk
               must contain at least one voxel
    allVoxelTs: structured np.array with a field name 'roi_voxel_ts' (and possible additional 
                fields), this field contains voxel time series
    consistencyType: str, definition of spatial consistency to be used; default:
          'pearson c' (mean Pearson correlation coefficient), other options:
    fTransform: bool, are the correlations Fisher f transformed before averaging
                when consistencyType = 'pearson c' (default=False)
    nCPUs = int, number of CPUs to be used for the parallel computing (default = 5)
    
    
    Returns:
    --------
    spatialConsistencies: list of doubles, spatial consistencies of the ROIs defined
                          by voxelIndices
    """
    cfg = {'allVoxelTs':allVoxelTs,'consistencyType':consistencyType,'fTransform':fTransform}
    paramSpace = [(cfg,voxelInd) for voxelInd in voxelIndices]
    pool = Pool(max_workers = nCPUs)
    spatialConsistencies = list(pool.map(calculateSpatialConsistency,paramSpace,chunksize=1))
    return spatialConsistencies

def calculateCorrelationToCentroid(params):
    """
    Calculates the mean correlation between ROI centroid and all voxel time
    series in the ROI.
    
    Parameters:
    -----------
    params: tuple, containing:
    
        cfg: dict, containing:
            allVoxelTs: structured np.array with a field name 'roi_voxel_ts' (and possible additional 
                    fields), this field contains voxel time series
                
        individualInfo: dict, containing:
            voxelIndices: np.array, indices of voxels; these indices should refer to voxels' 
                    locations in the file containing voxel time series; note that the chunk
                    must contain at least one voxel
            centroidIndex: int, index of the ROI centroid in the allVoxelTs array
            
    Returns:
    --------
    correlationToCentroid: float, mean correlation between centroid and voxel ts.
    """
    cfg = params[0]
    allVoxelTs = cfg['allVoxelTs']
    individualInfo = params[1]
    voxelIndices = individualInfo['voxelIndices']
    centroidIndex = individualInfo['centroidIndex']
    
    print 'Calculating correlation to centroid'
    
    if np.amax(voxelIndices.shape) == 0:
        correlationToCentroid = 0
        print "Detected an empty ROI, set correlation to centroid to 0."
    elif np.amax(voxelIndices.shape) == 1:
        correlationToCentroid = 1. # a single voxel is always fully consistent
    else: 
        voxelTs = allVoxelTs[voxelIndices,:]
        centroidTs = allVoxelTs[centroidIndex,:]
        correlations = []
        for voxel in voxelTs:
            correlations.append(np.corrcoef(voxel,centroidTs)[0][1])
        correlationToCentroid = np.mean(correlations)
    if np.isnan(correlationToCentroid):
        print 'nan detected!'
    return correlationToCentroid

def calculateCorrelationsToCentroidInParallel(voxelIndices,allVoxelTs,centroidIndices,nCPUs=5):
    """
    A wrapper function for calculating the correlation to ROI centroid in parallel across
    ROIs.
    
    Parameters:
    -----------
    voxelIndices: list of np.arrays, each array containing indices of voxels of one ROI; 
               these indices should refer to voxels' 
               locations in the file containing voxel time series; note that the chunk
               must contain at least one voxel
    allVoxelTs: structured np.array with a field name 'roi_voxel_ts' (and possible additional 
                fields), this field contains voxel time series
    centroidIndices: np.array of indices of ROI centroids in allVoxelTs
    nCPUs = int, number of CPUs to be used for the parallel computing (default = 5)
    
    
    Returns:
    --------
    correlationsToCentroid: list of doubles, correlations of the voxel ts defined
                          by voxelIndices to the ROICentroids
    """
    cfg = {'allVoxelTs':allVoxelTs}
    paramSpace = [(cfg,{'voxelIndices':voxelInd,'centroidIndex':centroidIndex}) for voxelInd,centroidIndex in zip(voxelIndices,centroidIndices)]
    pool = Pool(max_workers = nCPUs)
    correlationsToCentroid = list(pool.map(calculateCorrelationToCentroid,paramSpace,chunksize=1))
    return correlationsToCentroid
        
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
        if len(ROI.shape) == 1:
            ROI = np.array([ROI]) # adding an extra dimension to enable concatenation
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
    Returns:
    --------
    ROIlessNeighbors: dic, contains:
                      ROIlessIndices: NNeihbors x 1 np.array, indices of ROIless neighbor voxels in voxelCoordinates
                      ROIlessMap: NNeighbors x 3 np.array, coordinates of ROIless neighbor voxels
    """
    ROIMap = ROIInfo['ROIMaps'][ROIIndex]
    if len(ROIMap.shape) == 1:
        ROIMap = np.array([ROIMap]) # adding an outermost dimension to enable proper indexing later on
    for i, voxel in enumerate(ROIMap):
        neighbors = findNeighbors(voxel,allVoxels=voxelCoordinates)
        if i == 0:
            ROINeighbors = neighbors
        else:
            ROINeighbors = np.concatenate((ROINeighbors,neighbors),axis=0)
    if len(ROINeighbors) > 0:
        #print 'Found ' + str(len(ROINeighbors)) + 'ROIless neighbors'
        ROINeighbors = np.unique(ROINeighbors,axis=0) # removing dublicates
        ROIlessNeighbors = findROIlessVoxels(ROINeighbors,ROIInfo) 
        ROIlessMap = ROIlessNeighbors['ROIlessMap']
        ROIlessIndices = np.zeros(ROIlessMap.shape[0],dtype=int) # indices in the list of neighbors
        for i, voxel in enumerate(ROIlessMap):
            ROIlessIndices[i] = np.where((voxelCoordinates==voxel).all(axis=1)==1)[0][0] # finding indices in the voxelCoordinates array (indexing assumes that a voxel is present in the voxelCoordinates only once)
    else:
        ROIlessMap = np.array([])
        ROIlessIndices = np.array([])
    ROIlessNeighbors = {'ROIlessIndices':ROIlessIndices,'ROIlessMap':ROIlessMap}
    return ROIlessNeighbors

def updateROI(ROIIndex,voxelCoordinates,ROIInfo):
    """
    Updates the ROI by adding to it all of its neighbors that don't belong to
    any ROI yet.
    
    Parameters:
    -----------
    ROIIndex: int, index of the ROI in the lists of ROIInfo (see below)
    candidateVoxels: nVoxels x 3 np.array, coordinates (in voxels) of all voxels
    ROIInfo: dict, contains:
             ROIMaps: list of ROISizes x 3 np.arrays, coordinates of voxels
                           belonging to each ROI. len(ROIMaps) = nROIs.
             ROIVoxels: list of ROISizes x 1 np.arrays, indices of the voxels belonging
                     to each ROI. These indices refer to the columns of
                     the distance matrix, as well as to the rows of the
                     voxel time series file. len(ROIVoxels) = nROIs.
             ROISizes: nROIs x 1 np.array of ints. Sizes of ROIs defined as
                     number of voxels in the sphere
            
    Returns:
    --------
    ROIInfo: dict, input ROIInfo with lists updated with the new voxels of the ROI
    """
    ROIlessNeighbors = findROIlessNeighbors(ROIIndex,voxelCoordinates,ROIInfo)
    if max(ROIlessNeighbors['ROIlessMap'].shape) > 0:
        ROIMap = ROIInfo['ROIMaps'][ROIIndex]
        if len(ROIMap.shape) == 1: # adding an outermost dimension for successful concatenation later on
            ROIMap = np.array([ROIMap])
        ROIMap = np.concatenate((ROIMap,ROIlessNeighbors['ROIlessMap']),axis=0)
        ROIVoxels = ROIInfo['ROIVoxels'][ROIIndex]
        ROIVoxels = np.concatenate((ROIVoxels,ROIlessNeighbors['ROIlessIndices']),axis=0)
        ROIInfo['ROIMaps'][ROIIndex] = ROIMap
        ROIInfo['ROIVoxels'][ROIIndex] = ROIVoxels
        ROIInfo['ROISizes'][ROIIndex] = len(ROIVoxels)
    return ROIInfo

def addVoxel(ROIIndex, voxelIndex, ROIInfo, voxelCoordinates):
    """
    Adds the given voxel to the given ROI by updating the ROIInfo dictionary
    accordingly.
    
    Parameters:
    -----------
    ROIIndex: int, index of the ROI to be updated in ROIInfo['ROIMaps'] etc.
    voxelIndex: int, index of the voxel to be added to the ROI in voxelCoordinates
    ROIInfo: dict, containing:
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
                            
    Returns:
    --------
    ROIInfo: dict, original ROIInfo with the given ROI updated
    """
    ROIMap = ROIInfo['ROIMaps'][ROIIndex]
    ROIVoxels = ROIInfo['ROIVoxels'][ROIIndex]
    voxel = np.array([voxelCoordinates[voxelIndex,:]])
    if len(ROIMap.shape) == 1: # adding an outermost dimension for successful concatenation later on
            ROIMap = np.array([ROIMap])
    ROIMap = np.concatenate((ROIMap,voxel),axis=0)
    ROIVoxels = np.concatenate((ROIVoxels,np.array([voxelIndex])),axis=0)
    ROIInfo['ROIMaps'][ROIIndex] = ROIMap
    ROIInfo['ROIVoxels'][ROIIndex] = ROIVoxels
    ROIInfo['ROISizes'][ROIIndex] = len(ROIVoxels)
    return ROIInfo

def growROIs(ROICentroids,voxelCoordinates,names=''):
    """
    Divides the voxels to ROIs by growing each ROI spherewise by its ROIless
    neighbors. There is a first-come-first-served principle: if a voxel is neighbor
    of several ROIs, its added to the first one. Voxels that have no neighbors
    are excluded from the new ROI map.
    
    Parameters:
    -----------
    ROICentroids: nROIs x 3 np.array, coordinates of the centroids of the ROIs.
                  This can be a ROI centroid from an atlas but also any other
                  (arbitrary) point.
    voxelCoords: nVoxels x 3 np.array, coordinates of all voxels. Must be
                 in the same space with ROICentoids. Coordinates must be ordered by ROIs:
                 first all voxels of ROI 1, then all voxels of ROI 2, etc.
    names: list of strs, names of the ROIs, can be e.g. the anatomical name associated with
           the centroid. Default = ''.
    
    Returns:
    --------
    ROIInfo: dict, contains:
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
    if not isinstance(ROICentroids,list):
        ROICentroids = [np.array(centroid) for centroid in ROICentroids.tolist()] # ugly hack for ensuring that ROIInfo['ROIMaps'] is a list of arrays as expected later on
    nROIs = len(ROICentroids)
    ROIMaps = [centroid for centroid in ROICentroids]
    ROIInfo = {'ROICentroids':ROICentroids,'ROIMaps':ROIMaps,'ROISizes':np.ones(nROIs),'ROINames':names}
    ROIVoxels = []
    for centroid in ROICentroids:
        ROIVoxels.append(np.where((voxelCoordinates==centroid).all(axis=1)==1)[0]) # finding indices of centroids in voxelCoordinates
    ROIInfo['ROIVoxels'] = ROIVoxels
    
    nROIless = len(findROIlessVoxels(voxelCoordinates,ROIInfo)['ROIlessIndices'])
    
    while nROIless>0:
        
        for ROIIndex in range(nROIs):
            ROIInfo = updateROI(ROIIndex,voxelCoordinates,ROIInfo)
            
        print str(nROIless) + ' ROIless voxels found'
        nROIlessPrevious = nROIless
        nROIless = len(findROIlessVoxels(voxelCoordinates,ROIInfo)['ROIlessIndices'])

        if nROIless == nROIlessPrevious: # handling the case where the mask contains neighborless voxels that can't be added to any ROI
            ROIlessVoxels = findROIlessVoxels(voxelCoordinates,ROIInfo)['ROIlessMap']
            noNeighbors = np.zeros(nROIless)
            for i, voxel in enumerate(ROIlessVoxels):
                neighbors = findNeighbors(voxel,allVoxels=voxelCoordinates)
                if len(neighbors) == 0:
                    noNeighbors[i] = True
            if all(noNeighbors):
                break
            
    return ROIInfo

    
def growOptimizedROIs(cfg):
    """
    Starting from given centroids, grows a set of ROIs optimized in terms of
    spatial consistency. Optimization is based on a priority queue system: at each
    step, a measure is calculated for each centroid-voxel pair and the voxel
    with the highest measure value is added to the ROI in question.
    
    Possible measures are:
        1) the correlation between the ROI centroid time series and voxel time series
        2) the spatial consistency (mean Pearson correlation coefficient) of the voxels
           already in the ROI and the candidate voxel
    
    Parameters:
    -----------
    cfg: dict, contains:
         ROICentroids: nROIs x 3 np.array, coordinates of the centroids of the ROIs.
                  This can be a ROI centroid from an atlas but also any other
                  (arbitrary) point.
         voxelCoordinates: nVoxels x 3 np.array, coordinates of all voxels. Must be
                  in the same space with ROICentoids. Coordinates must be ordered by ROIs:
                  first all voxels of ROI 1, then all voxels of ROI 2, etc.
         names: list of strs, names of the ROIs, can be e.g. the anatomical name associated with
                  the centroid. Default = ''.
         allVoxelTs: nVoxels x nTime np.array, time series of all voxels
         threshold: float, the lowest centroid-voxel correlation that leads to adding a voxel (default=-1)
         targetFunction: str, measure that will be optimized (options: correlationWithCentroid, spatialConsistency)
         consistencyType: str, definition of spatial consistency to be used if 
                          targetFunction == 'spatialConsistency' (default: 'pearson c' (mean Pearson correlation coefficient))
         fTransform: bool, should Fisher Z transform be applied if targetFunction == 'spatialConsistency' 
                     (default=False)
    
    Returns:
    --------
    ROIInfo: dict, contains:
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
    # Setting up: reading parameters
    ROICentroids = cfg['ROICentroids']
    voxelCoordinates = cfg['voxelCoordinates']
    allVoxelTs = cfg['allVoxelTs']
    if 'threshold' in cfg.keys():
        threshold = cfg['threshold']
    else:
        threshold = -1
    targetFunction = cfg['targetFunction']
    if targetFunction == 'spatialConsistency':
        if 'consistencyType' in cfg.keys():
            consistencyType = cfg['consistencyType']
        else:
            consistencyType = 'pearson c'
        if 'fTransform' in cfg.keys():
            fTransform = cfg['fTransform']
        else:
            fTransform = False
    
    nROIs = len(ROICentroids)
    nTime = allVoxelTs.shape[1]
    
    # Setting up: defining initial priority queues, priority measures (centroid-voxel correlations) and candidate voxels to be added per ROI
    ROIMaps = [np.array(centroid) for centroid in ROICentroids]
    ROIVoxels = [np.array(np.where((voxelCoordinates==centroid).all(axis=1)==1)[0]) for centroid in ROICentroids]
    ROIInfo = {'ROIMaps':ROIMaps,'ROIVoxels':ROIVoxels,'ROISizes':np.ones(nROIs,dtype=int),'ROINames':cfg['names']}
    priorityQueues = [findROIlessNeighbors(i,voxelCoordinates,{'ROIMaps':ROIMaps})['ROIlessIndices'].tolist() for i in range(nROIs)] # priority queues change so it's better to keep them as lists
    centroidTs = np.zeros((nROIs,nTime))
    for i, ROIIndex in enumerate(ROIVoxels):
        centroidTs[i,:] = allVoxelTs[ROIIndex[0],:]

    additionCandidates = np.zeros(nROIs,dtype=int)
    maximalMeasures = np.zeros(nROIs)
    
    for i,(priorityQueue,centroid,ROI) in enumerate(zip(priorityQueues,centroidTs,ROIInfo['ROIVoxels'])):
        if targetFunction == 'correlationWithCentroid':
            priorityMeasures = [np.corrcoef(centroid,allVoxelTs[priorityIndex])[0][1] for priorityIndex in priorityQueue]
        elif targetFunction == 'spatialConsistency':
            priorityMeasures = np.zeros(len(priorityQueue))
            for j, voxel in enumerate(priorityQueue):
                voxelIndices = np.concatenate((ROI,np.array([voxel])))
                priorityMeasures[j] = calculateSpatialConsistency(({'allVoxelTs':allVoxelTs,'consistencyType':consistencyType,'fTransform':fTransform},voxelIndices))
        additionCandidates[i] = priorityQueue[np.argmax(priorityMeasures)]
        maximalMeasures[i] = np.amax(priorityMeasures)
        
    nInQueue = sum([len(priorityQueue) for priorityQueue in priorityQueues])
    selectedMeasures = []
    
    # Actual optimization takes place inside the while loop:    
    while nInQueue>0 and np.amax(maximalMeasures) >= threshold :
        print str(nInQueue) + ' voxels in priority queues'
        totalROISize = sum(len(ROIVox) for ROIVox in ROIInfo['ROIVoxels'])
        print str(totalROISize) + ' voxels in ROIs'

        # Selecting the ROI to be updated and voxel to be added to that ROI (based on the highest centroid-voxel correlation)
        ROIToUpdate = np.argmax(maximalMeasures)
        voxelToAdd = additionCandidates[ROIToUpdate]
        ROIInfo = addVoxel(ROIToUpdate,voxelToAdd,ROIInfo,voxelCoordinates)
        selectedMeasures.append(np.amax(maximalMeasures))
        
        # Updating priority queues: adding the ROIless neighbors of the updated ROI
        neighbors = findROIlessVoxels(findNeighbors(voxelCoordinates[voxelToAdd,:],allVoxels=voxelCoordinates),ROIInfo)['ROIlessMap']
        ROIlessIndices = [np.where((voxelCoordinates == neighbor).all(axis=1)==1)[0][0] for neighbor in neighbors]
        for ROIlessIndex in ROIlessIndices:
            if not ROIlessIndex in priorityQueues[ROIToUpdate]:
                priorityQueues[ROIToUpdate].append(ROIlessIndex)
        
        # Updating priority queues: removing the added voxel from all priority queues and updating the candidate voxels and maximal measures of the queues that changed        
        for i, priorityQueue in enumerate(priorityQueues):
            if voxelToAdd in priorityQueue:
                priorityQueue.remove(voxelToAdd) # voxel can belong to more than one priority que; let's remove it from all of them
                if len(priorityQueue) > 0:
                    if targetFunction == 'correlationWithCentroid':
                        priorityMeasures = [np.corrcoef(centroidTs[i],allVoxelTs[priorityIndex])[0][1] for priorityIndex in priorityQueue]
                    elif targetFunction == 'spatialConsistency':    
                        priorityMeasures = np.zeros(len(priorityQueue))
                        for j, voxel in enumerate(priorityQueue):
                            voxelIndices = np.concatenate((ROIInfo['ROIVoxels'][i],np.array([voxel])))
                            priorityMeasures[j] = calculateSpatialConsistency(({'allVoxelTs':allVoxelTs,'consistencyType':consistencyType,'fTransform':fTransform},voxelIndices))
                    additionCandidates[i] = priorityQueue[np.argmax(priorityMeasures)]
                    maximalMeasures[i] = np.amax(priorityMeasures)
                else:
                    maximalMeasures[i] = -1
        
        nInQueue = sum([len(priorityQueue) for priorityQueue in priorityQueues])

    return ROIInfo, selectedMeasures
    
# Multilayer construction:

def constructMultilayer(nLayers,layers=[],nodes=[],edges=[]):
    """
    Using the pymnet module by Mikko Kivelä (https://bitbucket.org/bolozna/multilayer-networks-library/overview),
    creates a multilayer network object.
    
    Parameters:
    -----------
    nLayers: int, number of layers in the network
    layers: list of strings, names of the layers (default = [], in which case the
            layers will be automatically named with letters in alphabetical order)
    nodes: list of objects, nodes of the network. All nodes are present at all layers.
           (default = [], in which case no nodes are added; they can be added later on.)
    edges: list of tuples, edges of the network. Each edge should be a 
           ((source node, source layer), (target node, target layer) or
           ((source node, source layer), (target node, target layer), weight) tuple (in the
           previous case, weight is set to 1) (default = [], in which case no nodes are added; 
           they can be added later on.)
           
    Returns:
    --------
    mnet: a pymnet multilayer network object
    """
    mnet = pymnet.MultilayerNetwork(aspects=1)
    if len(layers) == 0:
        letters = list(string.ascii_lowercase)
        nLetters = len(letters)
        if nLayers > nLetters:
            layers = []
            n = -1
            for i in range(0,nLayers):
                if np.remainder(i,nLetters) == 0:
                    n = n + 1
                layers.append(letters[np.remainder(i,nLetters)] + str(n))             
        else:
            layers = letters[0:nLayers]
    for layer in layers:
        mnet.add_layer(layer)
    if len(nodes) > 0:
        for node in nodes:
            mnet.add_node(node)
    if len(edges) > 0:
        for edge in edges:
            sourceNode = edge[0][0]
            sourceLayer = edge[0][1]
            targetNode = edge[1][0]
            targetLayer = edge[1][1]
            if len(edge) == 3:
                weight = edge[2]
            else:
                weight = 1
            mnet[sourceNode,sourceLayer][targetNode,targetLayer] = weight
    return mnet
       
    

    
# Accessories
    
def getDistribution(data, nBins):
    """
    Calculates the PDF of the given data
    
    Parameters:
    -----------
    data: a container of data points, e.g. list or np.array
    nBins: int, number of bins used to calculate the distribution
    
    Returns:
    --------
    pdf: np.array, PDF of the data
    binCenters: np.array, points where pdf has been calculated
    """
    count, binEdges, _ = binned_statistic(data, data, statistic='count', bins=nBins)
    pdf = count/float(np.sum(count))
    binCenters = 0.5*(binEdges[:-1]+binEdges[1:])
    
    return pdf, binCenters
    
def getCentroid(coordinates):
    """
    Calculates the centroid (rounded mean) of a given set of coordinates.
    
    Parameters:
    -----------
    coordinates: nRows x 3 np.array, points for calculating the centroid
    
    Returns:
    --------
    centroid: 1 x 3 np.array, centroid of the coordinates
    """
    centroid = np.round(np.mean(coordinates,axis=0)).astype(int)
    return centroid

def spaceToMNI(coordinate,resolution):
    """
    Transform a set of array indices into the MNI space.
    
    Parameters:
    -----------
    coordinate: 1 x 3 np.array, the indices to be transformed
    resolution: int, the physical distance (in mm) between the elements of the array
    
    Returns:
    --------
    MNICoordinate: 1 x 3 np.array, the corresponding coordinates in MNI space
    
    TODO: move this function of nifti_tools
    """
    origin = np.array([90, -126, -72])
    m = np.array([-1*resolution,resolution,resolution])
    MNICoordinate = coordinate * m + origin
    return MNICoordinate





   
        
    
        
    
  
        
    
    
