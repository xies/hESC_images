#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 20:54:47 2017

@author: mimi
"""

import numpy as np
from skimage import io, morphology
import pandas as pd
import seaborn as sb
from os import path

mtesr = pd.DataFrame()
for pos in xrange(6):
    filename = ''.join( ('/Users/mimi/Box Sync/hESCs/Images/',
                         '10_13_2017 H9 mTeSR1 DE anti-SOX17 anti-FoxA2/DE_DMSO_1/Pos',
                         str(pos), '/img_channel001_position00',str(pos),
                         '_time000000000_z000_Object Predictions.tiff') )
    mask = io.imread(filename)
    mask[mask == 2] = 0
    imsave(''.join( (path.dirname(filename),'/mask.tif') ) ,mask)

    
for pos in xrange(3):
    df = pd.DataFrame()
    filename = ''.join( ('/Users/mimi/Box Sync/hESCs/Images/',
                         '10_13_2017 H9 mTeSR1 DE anti-SOX17 anti-FoxA2/mTeSR1_DMSO_1/Pos',
                         str(pos), '/mask.tif') )
    mask = io.imread(filename)
    labels = morphology.label(mask)
    
    filename = ''.join( ('/Users/mimi/Box Sync/hESCs/Images/',
                     '10_13_2017 H9 mTeSR1 DE anti-SOX17 anti-FoxA2/mTeSR1_DMSO_1/Pos',
                     str(pos), '/img_channel002_position00',str(pos),
                     '_time000000000_z000.tif') )
    sox17 = io.imread(filename)
    
    filename = ''.join( ('/Users/mimi/Box Sync/hESCs/Images/',
                     '10_13_2017 H9 mTeSR1 DE anti-SOX17 anti-FoxA2/mTeSR1_DMSO_1/Pos',
                     str(pos), '/img_channel003_position00',str(pos),
                     '_time000000000_z000.tif') )
    fox2a = io.imread(filename)
    
    all_labels = np.unique(labels)
    mean_sox = np.zeros(all_labels.shape)
    mean_fox = np.zeros(all_labels.shape)
    
    for (i,l) in enumerate(all_labels):
        this_nucleus = labels == l
        mean_sox[i] = mean(sox17[this_nucleus])
        mean_fox[i] = mean(fox2a[this_nucleus])
        this_datum = pd.DataFrame()

    df['Sox17'] = mean_sox
    df['Fox2A'] = mean_fox
    df['Condition'] = 'mTeSR1'
    df['Position'] = pos
    
    mtesr = pd.concat((mtesr,df))
    print 'Done with ', pos
        
        
    