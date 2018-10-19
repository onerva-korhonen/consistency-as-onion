# -*- coding: utf-8 -*-
"""
Created on Fri May 11 14:18:58 2018

@author: aokorhon

A script for testing functions of the consistency as onion project
"""
import functions

import numpy as np

testDefSphericalROIs = False
testFindROIlessVoxels = True
testFindROIlessNeighbors = True
testUpdateROI = True

# testing defSphericalROIs
if testDefSphericalROIs: 
    ROIPath = 'atlases/brainnetome/brainnetome_thr25_4mm_rois.mat'
    centroids, _, voxelCoords = functions.readROICentroids(ROIPath, readVoxels=True)
    
    # testcase 1: ROIs with radius < 1 should contain only ROI centroids
    spheres = functions.defineSphericalROIs(centroids, voxelCoords, 0)
    sphereMaps = []
    for sphereMap in spheres['ROIMaps']:
        sphereMaps.append(sphereMap)
    sphereMaps = np.array(sphereMaps)
    
    diff = np.sum(np.abs(sphereMaps - centroids))
    
    if diff == 0:
        print "defSphericalROIs: testcase 1 OK"
        
    # testcase 2: small example space
    space = np.zeros((4,4,4))
    w1, w2, w3 = np.where(space==0)
    voxels = np.zeros((64,3))
    for i, (x, y, z) in enumerate(zip(w1,w2,w3)):
        voxels[i,:] = [x,y,z]
    centroids = np.array([[2,2,2],[1,2,3]])
    map1 = np.array([[2,2,2],[1,2,2],[2,1,2],[2,2,1],[3,2,2],[2,3,2],[2,2,3]])
    map1 = np.sort(map1, axis=0)
    map2 = np.array([[1,2,3],[0,2,3],[1,1,3],[1,2,2],[2,2,3],[1,3,3]]) #[1,2,4] is outside of the space and thus not included
    map2 = np.sort(map2, axis=0)
    
    spheres = functions.defineSphericalROIs(centroids, voxels, 1, resolution=1)
    diff = 0
    for sphereMap,trueMap in zip(spheres['ROIMaps'],[map1,map2]):
        sphereMap = np.sort(sphereMap, axis=0) # sorts the array columnwise; distorts the original rows but gives same output for all arrays with the same contents (independent of the row order)
        diff = diff + np.sum(np.abs(sphereMap - trueMap))
    
    if diff == 0:
        print "defSphericalROIs: testcase 2 OK"
        
    # testcase 3: changing sphere size:
    map3 = np.array([[2,2,2],[1,2,2],[2,1,2],[2,2,1],[3,2,2],[2,3,2],[2,2,3], # centroid + at distance 1
                     [0,2,2],[2,0,2],[2,2,0],[1,1,2],[2,1,1],[1,2,1], # at distance 2
                     [3,3,2],[2,3,3],[3,2,3],[1,3,2],[3,1,2],[2,1,3],[2,3,1],[1,2,3],[3,2,1], # at distance sqrt(2)
                     [1,1,1],[1,3,1],[1,1,3],[3,1,1],[3,3,1],[3,1,3],[1,3,3],[3,3,3]]) #at distance sqrt(3)
    map3 = np.sort(map3, axis=0)    

    spheres = functions.defineSphericalROIs(centroids[0], voxels, 2, resolution=1)
    sphereMap = np.sort(spheres['ROIMaps'][0], axis=0)
    diff = np.sum(np.abs(sphereMap - map3))
    
    if diff == 0:
        print 'defSphericalROIs: testcase 3 OK'
        
    # testcase 4: changing resolution
    # NOTE: I've removed the old testcase 4 on 13 Oct 2018: It assumed that voxel coordinates
    # are saved in millimeters while they are actually saved in voxels. So, the testcase didn't
    # match the reality and lead to wrong conclusions.
        
# testing getROIlessVoxels
        
if testFindROIlessVoxels:
    # testcase 1: a very small space with two small ROIs
    space = np.zeros((2,2,2))
    w1, w2, w3 = np.where(space==0)
    voxels = np.zeros((8,3))
    for i, (x, y, z) in enumerate(zip(w1,w2,w3)):
        voxels[i,:] = [x,y,z]
    
    ROI1 = np.array([[0,0,0],[0,1,0]])
    ROI2 = np.array([[0,0,1],[1,0,0],[0,1,1]])
    ROIMaps = [ROI1,ROI2]
    ROIVoxels = [np.array([0,1]),np.array([2,3,4])]
    ROIInfo = {'ROIMaps':ROIMaps,'ROIVoxels':ROIVoxels}
    
    trueROIless = np.array([[1,0,1],[1,1,0],[1,1,1]])
    trueROIlessIndices = [5,6,7]
    
    ROIlessVoxels = functions.findROIlessVoxels(voxels,ROIInfo)
    testROIless = ROIlessVoxels['ROIlessMap']
    testROIlessIndices = ROIlessVoxels['ROIlessIndices']
    
    mapDif = np.sum(np.abs(trueROIless - testROIless))
    indDif = np.sum(np.abs(np.array(trueROIlessIndices)-np.array(testROIlessIndices)))
    
    if mapDif == 0:
        print 'findROIlessVoxels: maps of ROIless voxels OK'
    if indDif == 0:
        print 'findROIlessVoxels: indices of ROIless voxels OK'
        
if testFindROIlessNeighbors:
    # testcase 1: a very small space with two small ROIs    
    space = np.zeros((2,2,2))
    w1, w2, w3 = np.where(space==0)
    voxels = np.zeros((8,3))
    for i, (x, y, z) in enumerate(zip(w1,w2,w3)):
        voxels[i,:] = [x,y,z]
        
    ROI1 = np.array([[0,0,0],[0,1,0],[1,1,0],[1,0,0]])
    ROI2 = np.array([[1,0,1]])
    ROIMaps = [ROI1,ROI2]
    ROIVoxels = [np.array([0,2,6,4]),np.array([5])]
    ROIInfo = {'ROIMaps':ROIMaps,'ROIVoxels':ROIVoxels}
    
    trueROIless1 = np.array([[0,0,1],[0,1,1],[1,1,1],])
    trueROIlessIndices1 = np.array([1,3,7])
    trueROIless2 = np.array([[0,0,1],[1,1,1]])
    trueROIlessIndices2 = np.array([1,7])
    
    testROIlessNeighbors1 = functions.findROIlessNeighbors(0,voxels,ROIInfo)
    testROIless1 = testROIlessNeighbors1['ROIlessMap']
    testROIlessIndices1 = testROIlessNeighbors1['ROIlessIndices']
    
    mapDif1 = np.sum(np.abs(trueROIless1 - testROIless1))
    indDif1 = np.sum(np.abs(trueROIlessIndices1 - testROIlessIndices1))
    
    testROIlessNeighbors2 = functions.findROIlessNeighbors(1,voxels,ROIInfo)
    testROIless2 = testROIlessNeighbors2['ROIlessMap']
    testROIlessIndices2 = testROIlessNeighbors2['ROIlessIndices']
    
    mapDif2 = np.sum(np.abs(trueROIless2 - testROIless2))
    indDif2 = np.sum(np.abs(trueROIlessIndices2 - testROIlessIndices2))
    
    if max(mapDif1,mapDif2) == 0:
        print 'findROIlessNeighbors: maps of ROIless neighbors OK'
        
    if max(indDif1,indDif2) == 0:
        print 'findROIlessNeighbors: indices of ROIless neighbors OK'

if testUpdateROI:
    # testcase 1: one small ROI in a small space
    space = np.zeros((2,2,2))
    w1, w2, w3 = np.where(space==0)
    voxels = np.zeros((8,3))
    for i, (x, y, z) in enumerate(zip(w1,w2,w3)):
        voxels[i,:] = [x,y,z]
    ROIMaps = [np.array([[0,0,0]])]
    ROIVoxels = [np.array([0])]
    ROIInfo = {'ROIMaps':ROIMaps,'ROIVoxels':ROIVoxels,'ROISizes':np.array([1])}
    
    trueROIMaps = [np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])]
    trueROIVoxels = [np.array([0,4,2,1])]
    
    testROIInfo = functions.updateROI(0,voxels,ROIInfo)
    testMaps = testROIInfo['ROIMaps']
    testROIVoxels = testROIInfo['ROIVoxels']
    
    mapCorrect = []
    indicesCorrect = []
    
    for ROI, testROI, ROIVoxel, testROIVoxel in zip(trueROIMaps,testMaps,trueROIVoxels,testROIVoxels):
        if not (len(ROI) == len(testROI)):
            mapCorrect.append(False)
        else:
            tempCorrect = np.zeros(ROI.shape[0])
            for i,testVoxel in enumerate(testROI):
                tempCorrect[i] = np.any((ROI == testVoxel).all(axis=1))
            mapCorrect.append(all(tempCorrect))
        if not (len(ROIVoxel) == len(testROIVoxel)):
            indicesCorrect.append(False)
        else:
            tempCorrect = np.zeros(len(ROIVoxel))
            for i, index in enumerate(testROIVoxel):
                tempCorrect[i] = np.any(ROIVoxel == index)
            indicesCorrect.append(all(tempCorrect))
    
    if all(mapCorrect):
        print 'updateROI: testcase 1: ROI map OK'
        
    if all(indicesCorrect):
        print 'updateROI: testcase 1: ROI voxel indices OK'
        
    if testROIInfo['ROISizes'] == np.array([4]):
        print 'updateROI: testcase 1: ROI sizes OK'
        
    # testcase 2: two small ROIs in a small space:
    ROI1 = np.array([[0,0,0]])
    ROI2 = np.array([[1,0,1]])
    ROIMaps = [ROI1,ROI2]
    ROIVoxels = [np.array([0]),np.array([5])]
    ROIInfo = {'ROIMaps':ROIMaps,'ROIVoxels':ROIVoxels,'ROISizes':np.array([1,1])}
    
    trueROI1Map = np.array([[0,0,0],[1,0,0],[0,0,1],[0,1,0]])
    trueROI1Voxels = np.array([0,1,2,4])
    trueROIMaps = [trueROI1Map,np.array([[1,0,1]])]
    trueROIVoxels = [trueROI1Voxels,np.array([5])]
    
    testROIInfo = functions.updateROI(0,voxels,ROIInfo)
    testROIMaps = testROIInfo['ROIMaps']
    testROIVoxels = testROIInfo['ROIVoxels']
    testROISizes = testROIInfo['ROISizes']
    
    mapCorrect = []
    indicesCorrect = []
    
    for ROI, testROI, ROIVoxel, testROIVoxel in zip(trueROIMaps,testMaps,trueROIVoxels,testROIVoxels):
        if not (len(ROI) == len(testROI)):
            mapCorrect.append(False)
        else:
            tempCorrect = np.zeros(ROI.shape[0])
            for i,testVoxel in enumerate(testROI):
                tempCorrect[i] = np.any((ROI == testVoxel).all(axis=1))
            mapCorrect.append(all(tempCorrect))
        if not (len(ROIVoxel) == len(testROIVoxel)):
            indicesCorrect.append(False)
        else:
            tempCorrect = np.zeros(len(ROIVoxel))
            for i, index in enumerate(testROIVoxel):
                tempCorrect[i] = np.any(ROIVoxel == index)
            indicesCorrect.append(all(tempCorrect))
            
    if all(mapCorrect) and len(testROIMaps) == len(trueROIMaps):
        print 'updateROI: testcase 2: ROI maps OK'
        
    if all(indicesCorrect) and len(testROIVoxels) == len(trueROIVoxels):
        print 'updateROI: testcase 2: ROI voxel indices OK'
        
    if np.all(testROISizes == np.array([7,1])):
        print 'updateROI: testcase 2: ROI sizes OK'
    
    
    
    
    



