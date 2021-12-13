#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Queda por hacer:
    1. Mirar bien cómo hacer una función del código de descarga de Jesús
    2. Implementar la función
    3. Mover cada archivo descargado directamente a un nuevo directorio por trigger
    4. Ya se podrá sacar un .txt por cada trigger

@author: jaimemorandominguez
"""

import numpy as np
import pandas as pd
import os
import datetime
import matlab.engine
import scipy.io as sio
import re

#path_to_mmia_files = '/media/lrg/012EE6107EB7CB6B/mmia_20'
#ssd_path = '/media/lrg/012EE6107EB7CB6B'
path_to_mmia_files = '/Users/jaimemorandominguez/Desktop/Final/MMIA_archivos/cdf'
ssd_path = '/Users/jaimemorandominguez/Desktop/test_descarga_GLM'
trigger_length = 4





def get_MMIA_triggers(path_to_mmia_files, trigger_length):
    
    print('Generating triggers from MMIA files...\n')
    ######### Ordering MMIA files by name in ascending time #########
    with os.scandir(path_to_mmia_files) as files:
        files = [file.name for file in files if file.is_file() and file.name.endswith('.cdf')]

    mmia_files_data = np.zeros((len(files),2))

    for i in range(len(files)):
        mmia_files_data[i,0] = int(i)
    
        time = files[i][50:69]
        year = time[0:4]
        day_of_year = datetime.datetime.strptime(time[0:10], "%Y-%m-%d").strftime("%j")
        hour = time[11:19]
        hour = float(hour[0:2])*3600 + float(hour[3:5])*60 + float(hour[6:8])
        mmia_files_data[i,1] = int(year)*10e7 + int(day_of_year)*10e4 + hour
    
        mmia_files_df = pd.DataFrame(data=mmia_files_data, columns=["index", "time"])
        mmia_files_df.sort_values(by=['time'], inplace=True)
        mmia_data = mmia_files_df.to_numpy()
    
    mmia_files = [] # List of MMIA file names in order of ascending time
    matches = []
    
    for i in range(len(files)):
        
        mmia_files.append([])
        
        # Generating matches vector (vector of existing dates)
        day_to_matches = datetime.datetime.strptime(str(mmia_data[i,1])[0:7], "%Y%j").strftime("%Y%m%d")
        if matches.count(day_to_matches) == 0:
            matches.append(day_to_matches)
        

    # Reorder MMIA files vector sorting by ascending time
    for i in range(len(files)):
        # Find every new position
        index = np.where(mmia_data[:,0] == i)[0][0]
        mmia_files[index] = files[i]
        
    ######### End of ordering #########
    
    # Trigger characterization
    
    trigger_info = []
    for i in range(len(matches)):
        trigger_info.append([])
    
    day = datetime.datetime.strptime(str(mmia_data[0,1])[0:7], "%Y%j").strftime("%Y%m%d")
    
    for i in range(len(mmia_files)):
        
        if i == 0:  # First file
            pivoting_time = mmia_data[0,1]
            file_list = []
            file_list.append(mmia_files[0])
            
        else:   #All of the rest of files
        
            if mmia_data[i,1] - pivoting_time <= trigger_length:   # Same trigger
                file_list.append(mmia_files[i])
                day = datetime.datetime.strptime(str(mmia_data[i,1])[0:7], "%Y%j").strftime("%Y%m%d")
            
            else:   # New trigger
                trigger_info[matches.index(day)].append(file_list)
                pivoting_time = mmia_data[i,1]
                file_list = []
                file_list.append(mmia_files[i])
                day = datetime.datetime.strptime(str(mmia_data[i,1])[0:7], "%Y%j").strftime("%Y%m%d")
                
    print('Trigger generation done!\n')
    return [matches, trigger_info]
        
def create_MMIA_trigger_directories(matches, trigger_info, path_to_mmia_files, ssd_path):
    
    print('Copying MMIA files into trigger directories...\n')
    os.system('mkdir ' + ssd_path + '/mmia_dirs')
    for i in range(len(trigger_info)):
        for j in range(len(trigger_info[i])):
            # Creating the trigger directory
            os.system('mkdir ' + ssd_path + '/mmia_dirs/' + matches[i] + '_' + str(j))
            for k in range(len(trigger_info[i][j])):
                file_name = trigger_info[i][j][k]
                # Moving the current file to its trigger directory
                os.system('cp ' + path_to_mmia_files + '/' + file_name + ' ' + ssd_path + '/mmia_dirs/' + matches[i] + '_' + str(j))
    
    print('MMIA classification done!\n')
    
def extract_trigger_info(ssd_path, trigger_info):
    
    mmia_mat_files_path = ssd_path + '/mmia_mat'
    
    os.system('mkdir ' + mmia_mat_files_path)
    
    trigger_limits = [None] * len(trigger_info)
    mmia_raw = [None] * len(trigger_info)
    
    for i in range(len(trigger_info)):
        
        # Creating the template for trigger data and info
        triggers = [None] * len(trigger_info[i])
        trigger_limits[i] = triggers
        mmia_raw[i] = triggers
        
        # Extracting data for every trigger
        for j in range(len(trigger_info[i])):
            
            print('Starting the MatLab engine and extracting data from .cdf files for day %d, trigger %d / %d...' % (int(matches[i]), j, len(trigger_info[i])))
            eng = matlab.engine.start_matlab()
            path = ssd_path + '/mmia_dirs/' + matches[i] + '_' + str(j) + '/'
            eng.workspace['str'] = path
            eng.MMIA_symplified_v5(nargout=0)
            eng.quit()
            wd = os.getcwd()
            os.system('mv '+wd+'/MMIA_data.mat ' + mmia_mat_files_path+'/'+matches[i]+'_' + str(j) +'_data.mat')
            os.system('mv '+wd+'/MMIA_space_time.mat ' + mmia_mat_files_path+'/'+matches[i]+'_' + str(j) +'_info.mat')
            
            # Filling mmia raw data and trigger info variables
            data_mat = sio.loadmat(mmia_mat_files_path+'/'+matches[i]+'_' + str(j) +'_data.mat')
            info_mat = sio.loadmat(mmia_mat_files_path+'/'+matches[i]+'_' + str(j) +'_info.mat')
            current_data = data_mat.get('MMIA_all')
            current_info = info_mat.get('space_time')
            mmia_raw[i][j] = current_data
            trigger_limits[i][j] = current_info
    
    return [mmia_raw, trigger_limits]
        
[matches, trigger_info] = get_MMIA_triggers(path_to_mmia_files, trigger_length)

create_MMIA_trigger_directories(matches, trigger_info, path_to_mmia_files, ssd_path)
        
[mmia_raw, trigger_limits] = extract_trigger_info(ssd_path, trigger_info)
        
