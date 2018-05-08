# -*- coding: utf-8 -*-
"""
Created on Tue May  8 13:27:05 2018

@author: aokorhon

Functions for analysing the relationship between ROI's size and spatial consistency. Call either from
a frontend script or interactively.
"""
import numpy as np
from scipy import io


def readROICentroids(ROIInfoFile):
    """
    Reads ROI data, in particular the coordinates of ROI centroids.
    
    Parameters:
    -----------
    ROIInfoFile: str, path to the file that contains information about ROIs
        
    Returns:
    --------
    ROICentroids: Nx3 np.array, coordinates of the centroids of the ROIs
    """
    infoData = io.loadmat(ROIInfoFile)
    ROIInfo = infoData['rois'][0].T # transpose makes looping easier
    nROIs = ROIInfo.shape[0] # number of ROIs
    ROIcentroids = np.zeros((nROIs,3))
    for i, ROI in enumerate(ROIInfo):
        ROIcentroids[i,:] = ROI['centroid'][0]
    return ROIcentroids
        
    
