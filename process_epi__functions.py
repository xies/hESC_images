#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 13:26:58 2017

@author: xies@stanford

"""
import numpy as np
from skimage import io, filters, morphology, measure, feature
from scipy.ndimage import distance_transform_edt
import pandas as pd

"""

"""
def segment(image,Nerode=10,min_size=300):
    
    # Otsu threshold
    global_thresh = filters.threshold_otsu(image)
    mask = image > global_thresh
    
    # Remove salt & pepper noise
    mask = morphology.closing(mask)
    mask = morphology.opening(mask)
    
    # Remove objects and holes that are too small
    mask = morphology.remove_small_objects(mask,min_size)
    mask = morphology.remove_small_holes(mask,min_size)
    
    # Iteratively erode the raw mask for Nerode times to obtain more reliable
    # foreground markers
    foreground = mask
    for i in xrange(Nerode):
        foreground = morphology.binary_erosion(foreground)
    
    # Euclidean distance transform to separate 'touching' foregrounds
    distTransform = distance_transform_edt(foreground)
    watershedImg = -distTransform / distTransform.max() + 1
    # Use local maixma to find seeds and label them
    seeds = feature.peak_local_max(watershedImg, min_distance=20,indices=False)
    seeds = measure.label(~seeds)
    
    # Watershed segment
    nuclei = morphology.watershed( watershedImg, seeds)
    nuclei[~mask]=0
    # Label 1 is actually the background -> Need to set to 0
    nuclei[nuclei == 1] = 0
    
    return nuclei

    
"""
Get geometric or intensity statistics from segmented objects

Parameters

"""
def measure_props(super_idx,nuclei,dapi,sox,nanog):

    nucProps = measure.regionprops(nuclei, dapi) #
    soxProps = measure.regionprops(nuclei, sox)
    nanogProps = measure.regionprops(nuclei, nanog)

    # Measure the following properties
    columns = ('position','nuclear_area','nuclear_perimeter',
    'mean_dapi','mean_sox','mean_nanog')
    
    object_properties = []
    indices = []
    
    for [i,d] in enumerate(nucProps):
        object_properties.append([super_idx,d.area,d.perimeter,d.mean_intensity,
                                  soxProps[i].mean_intensity,nanogProps[i].mean_intensity]
                                  )
        indices.append(d.label + super_idx*10000)
        #tentatively using +10000 to aovoid index-clashing
        
    #if not len(indices):
    #    object_properties = pd.DataFrame([],index=[])
    indices = pd.Index(indices, name='objectID')
    object_properties = pd.DataFrame(object_properties, index=indices, columns=columns)
    
    return object_properties
    
