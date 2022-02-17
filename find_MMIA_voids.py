#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 11:35:58 2022

@author: jaimemorandominguez
"""

import numpy as np
import pickle
import matplotlib.pyplot as plt
import math

def split_MMIA_events(event, split_window):
    
    derivative = np.diff(event[:,1])/np.diff(event[:,0]) # dy/dx
    non_linear_positions = np.zeros((len(derivative),2))
    
    for i in range(len(derivative)):
        non_linear_positions[i,0] = event[i,0]
        if derivative[i] == 0:
            non_linear_positions[i,1] = 0.5
    
    in_event_pos = []

    for i in range(len(non_linear_positions)):
        if non_linear_positions[i,1] != 0.5:
            in_event_pos.append(i)
    
    non_linear_positions = np.delete(non_linear_positions, in_event_pos, axis=0)
    
    # Miro saltos sin datos de más de split_window samples
    
    event_borders = []
    event_borders.append(0)
    
    for i in range(1,len(non_linear_positions)):
        
        current_pos_in_signal = np.where(event[:,0] <= non_linear_positions[i,0])[0][-1]
        prev_pos_in_signal = np.where(event[:,0] <= non_linear_positions[i-1,0])[0][-1]
        
        if (current_pos_in_signal - prev_pos_in_signal) >= split_window:

            event_borders.append(prev_pos_in_signal)
            event_borders.append(current_pos_in_signal)
    event_borders.append(-1)

    return event_borders






ssd_path = '/Users/jaimemorandominguez/Desktop/test_descarga_GLM_2'
xcorr_bin = ssd_path + '/xcorr_bin'

# Uploading signals
f = open(xcorr_bin + '/20200607_signals.pckl', 'rb')
[GLM_xcorr, MMIA_xcorr, GLM_xcorr_norm, MMIA_xcorr_norm] = pickle.load(f)
f.close()

event_borders = split_MMIA_events(MMIA_xcorr_norm[1], 3000)

plt.figure()
plt.plot(MMIA_xcorr_norm[1][:,0], MMIA_xcorr_norm[1][:,1])
#plt.scatter(non_linear_positions[:,0], non_linear_positions[:,1], color='r')
plt.show()

if len(event_borders) != 1:
    print(event_borders)
    events_num = int(len(event_borders)/2)

    new_event_list = [None] * events_num

    for i in range(events_num):
        my_start = event_borders[2*i]
        my_end = event_borders[2*i+1]
        my_event = MMIA_xcorr_norm[1][my_start:my_end,:]
        new_event_list[i] = my_event
    plt.figure()
    plt.plot(MMIA_xcorr_norm[1][:,0], MMIA_xcorr_norm[1][:,1])
    for i in range(events_num):
        #plt.figure()
        plt.plot(new_event_list[i][:,0], new_event_list[i][:,1])
        #plt.show()