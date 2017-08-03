#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 11:34:23 2017

@author: mimi
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from skimage import io, filters, morphology, measure, util, feature
from scipy.ndimage import distance_transform_edt
import pandas as pd
import seaborn as sns

um_per_px = 0.0000023

which_replicate = '1'
conditions = ['DMSO','400nM_PD','9uM_RO']
c = 2

all_props = pd.DataFrame()
for frame in range(11):
    frame += 1
    base = "".join(['/Users/mimi/Box Sync/hESCs/Images/7_30_2017_H9_hESC_1day_inhibitors_DAPI_Sox2488_Nanog594/',
            'Rep', which_replicate, '/', conditions[c], '/',str(frame), '/'])
    
    filename = ''.join([base,str(frame),'_01.tif'])
    im = io.imread(filename)
    dapi = filters.gaussian(im,sigma=1)
    
    filename = ''.join([base,str(frame),'_02.tif'])
    im = io.imread(filename)
    sox = filters.gaussian(im,sigma=1)
    
    filename = ''.join([base,str(frame),'_03.tif'])
    im = io.imread(filename)
    nanog = filters.gaussian(im,sigma=1)

    nuclei = segment(dapi,Nerode=11,min_size=300)
    nuclei = morphology.remove_small_objects(nuclei,min_size=300)
    properties = measure_props(frame,nuclei,dapi,sox,nanog)
    
    #Filter by obj size
    I = (properties['nuclear_area'] > 2000) | (properties['nuclear_area'] < 300)
    idx = properties.index[I] - 10000 * frame
    
    notrightsize = np.isnan(nuclei)
    for i in idx:
        notrightsize = notrightsize | (nuclei == i)
    
    good_labels = nuclei.copy(); good_labels[notrightsize] = 0
    nucs = properties[~I]
    
    sp.misc.imsave(''.join([base,'nuclei.tif']),good_labels)
    
    all_props = pd.concat([all_props,nucs])
    print "Done with ", frame


