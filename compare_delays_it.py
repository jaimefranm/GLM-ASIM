#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 20:53:09 2022

@author: jaimemorandominguez
"""

import os
import pickle
import matplotlib.pyplot as plt

def merge_days(matches, bin_path, filename_end_to_merge):
    # La idea de la funcion es generar una matrix de datos independiente de
    # los datos de la matrix (delays, numero de muestras...)
    
    # Para ello los datos deben ser guardados por dia en binarios de antemano!
    
    # Generation of the blank datatable
    matrix = [None] * len(matches)
    
    # If there is just one possible name type in bin_path
    if filename_end_to_merge == 'no_type':
        with os.scandir(bin_path) as files:
            files = [file.name for file in files if file.is_file() and file.name.endswith('.pckl')]
    else: # Multiple namings in bin_path
        with os.scandir(bin_path) as files:
            files = [file.name for file in files if file.is_file() and file.name.endswith(filename_end_to_merge + '.pckl')]
    
    for i in range(len(files)):
        
        # Get positioning in 'matrix'
        day = files[i][0:8]
        pos_in_matches = matches.index(day)
        
        # Upload contents of the binary into 'matrix'
        if filename_end_to_merge == 'no_type':
            f = open(bin_path + '/' + day + '.pckl', 'rb')
            day_info_vector = pickle.load(f)
            f.close()
        else: # Multiple namings in bin_path
            f = open(bin_path + '/' + day + '_' + filename_end_to_merge + '.pckl', 'rb')
            day_info_vector = pickle.load(f)
            f.close()
        
        matrix[pos_in_matches] = day_info_vector
    
    return matrix


general_variables_path = '/Users/jaimemorandominguez/Desktop/delays_it/general_variables_bin'
statistics_bin_it = '/Users/jaimemorandominguez/Desktop/delays_it/stats_bin_it'
statistics_bin = '/Volumes/Jaime_F_HD/mmia_2020/stats_bin'

f = open(general_variables_path+'/matches.pckl', 'rb')
matches = pickle.load(f)
f.close()

delays_it_mat = merge_days(matches, statistics_bin_it, 'delays')
delays_mat = merge_days(matches, statistics_bin, 'delays')

delays_it = []
delays = []

it_disparado_counter = 0
disparado_counter = 0

for i in range(len(delays_mat)):
    for j in range(len(delays_mat[i])):
        if type(delays_it_mat[i][j]) == int: 
            delays_it.append(delays_it_mat[i][j])
            if abs(delays_it_mat[i][j]) >= 10000:
                it_disparado_counter = it_disparado_counter+1
        if type(delays_mat[i][j]) == int:
            delays.append(delays_mat[i][j])
            if abs(delays_mat[i][j]) >= 10000:
                disparado_counter = disparado_counter+1

plt.figure()
plt.hist(delays_it, bins = 200, rwidth=0.85, range=[-100000, -8000])
plt.title('Iterative')
plt.xlabel('Delay [s]')
plt.ylabel('Frequency')
plt.grid('on')
plt.show()
# Clear the current axes
plt.cla() 
# Clear the current figure
plt.clf() 
# Closes all the figure windows
plt.close('all')

# Delays histogram for all delays in range
plt.figure()
plt.hist(delays, bins = 200, rwidth=0.85, range=[-100000, -8000])
plt.title('No iterative')
plt.xlabel('Delay [s]')
plt.ylabel('Frequency')
plt.grid('on')
plt.show()
# Clear the current axes
plt.cla() 
# Clear the current figure
plt.clf() 
# Closes all the figure windows
plt.close('all')