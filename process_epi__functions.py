#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 13:26:58 2017

@author: mimi
"""
import numpy as np
from skimage import io, filters, morphology, measure, feature
from scipy.ndimage import distance_transform_edt
import pandas as pd

def segment(image,Nerode=10,min_size=300):
    
    global_thresh = filters.threshold_otsu(image)
    mask = image > global_thresh
    mask = morphology.closing(mask)
    mask = morphology.opening(mask)
    mask = morphology.remove_small_objects(mask,min_size)
    
    foreground = mask
    for i in xrange(Nerode):
        foreground = morphology.binary_erosion(foreground)
    
    distTransform = distance_transform_edt(foreground)
    watershedImg = -distTransform / distTransform.max() + 1
    seeds = feature.peak_local_max(watershedImg, min_distance=20,indices=False)
    seeds = measure.label(~seeds)
    
    nuclei = morphology.watershed( watershedImg, seeds)
    nuclei[~mask]=0
    nuclei[nuclei == 1] = 0
    return nuclei

    
"""
Get object statistics
"""
def measure_props(super_idx,nuclei,dapi,sox,nanog):

    nucProps = measure.regionprops(nuclei, dapi) #
    soxProps = measure.regionprops(nuclei, sox)
    nanogProps = measure.regionprops(nuclei, nanog)
    #rbProps = measure.regionprops(labels,rb)
    columns = ('position','nuclear_area','mean_dapi','mean_sox','mean_nanog')
    
    object_properties = []
    indices = []
    
    for [i,d] in enumerate(nucProps):
        area = d.area
        object_properties.append([super_idx,area,d.mean_intensity,
                                  soxProps[i].mean_intensity,nanogProps[i].mean_intensity]
                                  )
        indices.append(d.label + super_idx*10000)
        
    #if not len(indices):
    #    object_properties = pd.DataFrame([],index=[])
    indices = pd.Index(indices, name='objectID')
    object_properties = pd.DataFrame(object_properties, index=indices, columns=columns)
    
    return object_properties
    
