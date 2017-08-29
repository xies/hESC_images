#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 18:09:47 2017

@author: mimi
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 11:34:23 2017

@author: mimi

@todo: 1) isoperimetric threshold

"""

# Core modules
import numpy as np
import scipy as sp
import pandas as pd

# Image modules
from skimage import io
from libtiff import TIFFimage

# Plotting modules
import matplotlib.pyplot as plt
import seaborn as sns


um_per_px = 0.0000023

# Object size threshold in pixels
min_size = 300
max_size = 2000

# Isoperimetric threshold
ip_threshold = 17

# Which experiment set to process
conditions = ['DMSO','500nM_PD','1uM_Roscov']
c = 1

# Go through all positions available for the replicate + condition
all_props = pd.DataFrame()
for position in range(10):
    position += 1
    
    # Filename base
    base = "".join(['/Users/mimi/Box Sync/hESCs/Images/8_10_2017_H9_hESC_22hr_inhibitors_DAPI_Sox2488_Nanog594/',
            conditions[c], '/',str(position), '/'])
    
    # DAPI channel is _01
    filename = ''.join([base,str(position),'.tif'])
    im = io.imread(filename)
    im = filters.gaussian(im,sigma=1)
    dapi = im[...,1]
    sox = im[...,2]
    nanog = im[...,3]

    # Segment the DAPI channel
    nuclei = segment(dapi,Nerode=10,min_size=300)
    nuclei = morphology.remove_small_objects(nuclei,min_size=300) #pre-filter
    
    # Get object properties
    properties = measure_props(position,nuclei,dapi,sox,nanog)
    
    #Filter by obj size again
    # 'I' are the 'bad' objects
    I = (properties['nuclear_area'] > max_size) | (properties['nuclear_area'] < min_size)
    
    #Filter by perimeter/area ratio
    L2 = properties['nuclear_perimeter'] ** 2
    A = properties['nuclear_area']
    I = I | (L2 / A > ip_threshold).as_matrix()
    
    # Find the proper index in the labeled image
    idx = properties.index[I] - 10000 * position
    # Filter bad idx from labeled image
    notrightsize = np.isnan(nuclei)
    for i in idx:
        notrightsize = notrightsize | (nuclei == i)
    
    # Save only good nuclei in image and properties table
    good_labels = nuclei.copy(); good_labels[notrightsize] = 0
    nucs = properties[~I]
    
    # Save the segmentation as 32-bit TIFF
    # NB: only use the libtiff package since we need 32-bit
    tiff = TIFFimage(good_labels, description='')
    tiff.write_file(''.join((base, 'nuclei.tif')),
                    compression='none')
    # DEBUG ONLY: also save the raw segmentation
#    tiff = TIFFimage(nuclei, description='')
#    tiff.write_file(''.join((base, 'raw_nuclei.tif')),
#                    compression='none')   

    # Save the DataFrame as CSV
    nucs.to_csv(''.join((base,'quantification.csv')), sep='\t')
    
    all_props = pd.concat([all_props,nucs])
    print "Done with ", position


