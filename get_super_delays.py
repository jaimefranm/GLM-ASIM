#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 16 19:07:25 2022

@author: jaimemorandominguez
"""

import pickle
import os


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




pick_path = '/Volumes/Jaime_F_HD/mmia_2020/mmia_20/'
general_variables_path = '/Volumes/Jaime_F_HD/mmia_2020/general_variables_bin'
statistics_bin = '/Users/jaimemorandominguez/Desktop/delays_it/stats_bin'


f = open(general_variables_path+'/matches.pckl', 'rb')
matches = pickle.load(f)
f.close()

# Import filenames bla bla bla

f = open(general_variables_path+'/event_filenames.pckl', 'rb')
event_filenames = pickle.load(f)
f.close()

# Create delays matrix

delays_mat = merge_days(matches, statistics_bin, 'delays')

# Copy those filenames into a new folder in Desktop

super_delay_list = []

for i in range(len(event_filenames)):
    for j in range(len(event_filenames[i])):
        if (type(delays_mat[i][j]) == int) and (delays_mat[i][j] >= 8000):
            for k in range(len(event_filenames[i][j])):
                super_delay_list.append(event_filenames[i][j][k])


os.system('mkdir ~/Desktop/super_delays')
for i in range(len(super_delay_list)):
    os.system('cp '+pick_path+super_delay_list[i]+' ~/Desktop/super_delays')









