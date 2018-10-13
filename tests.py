# -*- coding: utf-8 -*-
"""
Created on Fri May 11 14:18:58 2018

@author: aokorhon

A script for testing functions of the consistency as onion project
"""
import functions

import numpy as np

testDefSphericalROIs = True

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
    
    



