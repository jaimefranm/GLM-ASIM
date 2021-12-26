#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Queda por hacer:
    1. Mirar bien cómo hacer una función del código de descarga de Jesús
    2. Implementar la función
    3. Mover cada archivo descargado directamente a un nuevo directorio por trigger
    4. Ya se podrá sacar un .txt por cada trigger
    5. Cuidado con poder subir datos de los .mat a Python directamente sin
        tener que extraer los .mat (función aparte que sólo lea el directorio
        y suba los datos de data e info)

@author: jaimemorandominguez
"""

import numpy as np
import pandas as pd
import os
import datetime
import matlab.engine
import scipy.io as sio
import re
from google.cloud import storage
import math

# PC UPC
#path_to_mmia_files = '/media/lrg/012EE6107EB7CB6B/mmia_20'
#path_to_mmia_files = '/home/lrg/Desktop/test_cdf'
#ssd_path = '/media/lrg/012EE6107EB7CB6B'

path_to_mmia_files = '/media/lrg/mmia_20'
ssd_path = '/media/lrg'

# LOCAL mac

# For testing
#path_to_mmia_files = '/Users/jaimemorandominguez/Desktop/Final/MMIA_archivos/cdf'
#path_to_mmia_files = '/Users/jaimemorandominguez/Desktop/test_cdf'
#ssd_path = '/Users/jaimemorandominguez/Desktop/test_descarga_GLM'

# For doing something
#path_to_mmia_files = '/Volumes/Jaime_F_HD/mmia_2020/mmia_20'
#ssd_path = '/Volumes/Jaime_F_HD/mmia_2020'

trigger_length = 2 # [s]
pre_extracted_MMIA = True
pre_downloaded_GLM = False


def get_MMIA_triggers(path_to_mmia_files, trigger_length):

    print('Generating triggers from MMIA files...\n')
    ######### Ordering MMIA files by name in ascending time #########

    # Moving away all files with size 0
    with os.scandir(path_to_mmia_files) as files:
        files = [file.name for file in files if file.is_file() and file.name.endswith('.cdf')]

    files_size = [
        (f,os.stat(os.path.join(path_to_mmia_files, f)).st_size)
        for f in files
    ]
    path_to_no_size = path_to_mmia_files + '/no_size'
    os.system('mkdir '+ path_to_no_size)

    for i in range(len(files_size)):
        if files_size[i][1] == 0:
            os.system('mv ' + path_to_mmia_files + '/' + files_size[i][0] + ' ' + path_to_no_size)


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
                if i == (len(mmia_files)-1):
                    trigger_info[matches.index(day)].append(file_list)

            else:   # New trigger
                trigger_info[matches.index(day)].append(file_list)
                pivoting_time = mmia_data[i,1]
                file_list = []
                file_list.append(mmia_files[i])
                day = datetime.datetime.strptime(str(mmia_data[i,1])[0:7], "%Y%j").strftime("%Y%m%d")

    print('Trigger generation done!\n')
    return [matches, trigger_info]

def create_MMIA_trigger_directories(matches, trigger_filenames, path_to_mmia_files, ssd_path):

    print('Copying MMIA files into trigger directories...\n')
    os.system('mkdir ' + ssd_path + '/mmia_dirs')
    for i in range(len(trigger_filenames)):
        for j in range(len(trigger_filenames[i])):
            # Creating the trigger directory
            os.system('mkdir ' + ssd_path + '/mmia_dirs/' + matches[i] + '_' + str(j))
            for k in range(len(trigger_filenames[i][j])):
                file_name = trigger_filenames[i][j][k]
                # Moving the current file to its trigger directory
                os.system('cp ' + path_to_mmia_files + '/' + file_name + ' ' + ssd_path + '/mmia_dirs/' + matches[i] + '_' + str(j))

    print('MMIA classification done!\n')

def extract_trigger_info(ssd_path, trigger_filenames, matches):

    mmia_mat_files_path = ssd_path + '/mmia_mat'

    os.system('mkdir ' + mmia_mat_files_path)

    trigger_limits = [None] * len(trigger_filenames)
    mmia_raw = [None] * len(trigger_filenames)

    for i in range(150,len(trigger_filenames)): # ------------------------------- CAMBIAAAAAAAR

        # Creating the template for trigger data and info
        triggers1 = [None] * len(trigger_filenames[i])
        triggers2 = [None] * len(trigger_filenames[i])
        trigger_limits[i] = triggers1
        mmia_raw[i] = triggers2

        # Extracting data for every trigger
        for j in range(len(trigger_filenames[i])):

            print('Starting the MatLab engine and extracting data from .cdf files for day %d (%d / %d), trigger %d / %d...' % (int(matches[i]), i+1, len(matches), j, len(trigger_filenames[i])))
            eng = matlab.engine.start_matlab()
            path = ssd_path + '/mmia_dirs/' + matches[i] + '_' + str(j) + '/'
            eng.workspace['str'] = path
            eng.MMIA_symplified_v5(nargout=0)
            eng.quit()
            wd = os.getcwd()

            with os.scandir(wd) as files:
                files = [file.name for file in files if file.is_file() and file.name.endswith('.mat')]

            if len(files) != 0:

                os.system('mv '+wd+'/MMIA_data.mat ' + mmia_mat_files_path+'/'+matches[i]+'_' + str(j) +'_data.mat')
                os.system('mv '+wd+'/MMIA_space_time.mat ' + mmia_mat_files_path+'/'+matches[i]+'_' + str(j) +'_info.mat')

                # Filling mmia raw data and trigger info variables
                data_mat = sio.loadmat(mmia_mat_files_path+'/'+matches[i]+'_' + str(j) +'_data.mat')
                info_mat = sio.loadmat(mmia_mat_files_path+'/'+matches[i]+'_' + str(j) +'_info.mat')
                current_data = data_mat.get('MMIA_all')
                current_info = info_mat.get('space_time')
                mmia_raw[i][j] = current_data
                trigger_limits[i][j] = current_info
            else:
                print('No MMIA data could be extracted for day %s trigger %d' % (matches[i], j))
        print(' ')  # For separating dates

    return [mmia_raw, trigger_limits]

def download_GLM(ssd_path, trigger_filenames, mmia_raw):

    print("Downloading GLM's .nc files from Google Cloud Storage...")

    path_to_downloaded_glm_nc = ssd_path + '/GLM_downloaded_nc_files'

    os.system('mkdir ' + path_to_downloaded_glm_nc)

    for i in range(len(trigger_filenames)):
        for j in range(len(trigger_filenames[i])):
            if type(mmia_raw[i][j]) == np.ndarray:
                current_trigger_download_path = path_to_downloaded_glm_nc + '/' + matches[i] + '_' + str(j)
                os.system('mkdir ' + current_trigger_download_path)

                date = trigger_filenames[i][j][0][50:60]
                year = date[0:4]
                day_year = datetime.datetime.strptime(date, "%Y-%m-%d").strftime("%j")

                download_glm_from_google(ssd_path, year, day_year, mmia_raw[i][j][0,0], mmia_raw[i][j][-1,0])

                with os.scandir(ssd_path) as files:
                    files = [file.name for file in files if file.is_file() and file.name.endswith('.nc')]

                for k in range(len(files)):
                    os.system('mv ' + ssd_path + '/' + files[k] + ' ' + current_trigger_download_path)

    print('Download done! Your .nc files can be accessed per trigger at %s' % path_to_downloaded_glm_nc)

def download_glm_from_google(ssd_path, year, day_year, t_ini, t_end):

    hour_ini = str(int(t_ini // 3600))
    if len(hour_ini) == 1:
        hour_ini = '0'+hour_ini
    min_ini = str(int((t_ini/3600 - t_ini//3600) * 60))
    seg_ini = str(int(((t_ini/3600 - t_ini//3600) * 60 - int((t_ini/3600 - t_ini//3600) * 60)) * 60))
    mseg_ini = '000'


    hour_end = str(int(t_end // 3600))
    if len(hour_end) == 1:
        hour_end = '0'+hour_end
    min_end = str(int((t_end/3600 - t_end//3600) * 60))
    seg_end = str(int(((t_end/3600 - t_end//3600) * 60 - int((t_end/3600 - t_end//3600) * 60)) * 60))
    mseg_end = '999'



    # Instantiates a client
    storage_client = storage.Client()
    # Get GCS bucket
    bucket_name='gcp-public-data-goes-16'
    # GLM-L2-LCFA/ + AÑO + DAY OF YEAR + HORA (sólo hora entera)
    prefix='GLM-L2-LCFA/'+year+'/'+day_year+'/'+hour_ini # Objeto
    delimiter=','
    # Get blobs in bucket (including all subdirectories)
    bucket = storage_client.get_bucket(bucket_name)
    blobs = storage_client.list_blobs(bucket_name, prefix=prefix, delimiter=delimiter)
    # print("Blobs:")
    glm_name=''
    for blob in blobs:
        glm_name +=(blob.name+'\n') # Ficheros .nc GLM en "String" para la hora completa de búsqueda

    glm_name_list = glm_name.splitlines() # "Nombre de los ficheros String pero en list"

    # String list permite identificar por cada "línea" como si se tratase de un array.

    # Find la fecha inicial de los ficheros.
    date_files= re.findall(r'_s(.*)_e\S', glm_name)

    # Convierte la fecha de string a entero o flotante
    date_files_float = list(map(int, date_files))

    # Fecha de búsqueda, primero en string para concatenar y posterior en int
    # AÑO + DIA DEL AÑO + HORA + MINUTO + mseg
    find_ini = year + day_year + hour_ini + min_ini + mseg_ini
    find_end = year + day_year + hour_end + min_end + mseg_end

    int_find_ini=int(find_ini)
    int_find_end=int(find_end)

    # Posiciones de los ficheros según las fechas de búsqueda
    date_in = np.array([date_files_float])
    pos_in=(np.where(np.logical_and(date_in>=int_find_ini, date_in<=int_find_end)))

    # Seleciona los ficheros de ls búsqueda y los almacena en destination_file_name
    for x in range(len(pos_in[1])):
        # file=bucket_name+'/'+glm_name_list[pos_in[1][x]]
        file=glm_name_list[pos_in[1][x]]
        # HAY UN BUG CUANDO LOS ALMACENO CON NOMBRES DIFERENTES AL ORIGINAL
        destination_file_name=ssd_path+'/'+str(glm_name_list[pos_in[1][x]]).replace("/","_")
        blob=bucket.blob(file)
        blob.download_to_filename(destination_file_name)
        print("Blob {} downloaded to {}.".format(file, destination_file_name))

def upload_MMIA_mats(ssd_path, trigger_filenames, matches):

    mmia_mat_files_path = ssd_path + '/mmia_mat'

    with os.scandir(mmia_mat_files_path) as data_files:
        data_files = [file.name for file in data_files if file.is_file() and file.name.endswith('data.mat')]

    with os.scandir(mmia_mat_files_path) as info_files:
        info_files = [file.name for file in info_files if file.is_file() and file.name.endswith('info.mat')]

    trigger_limits = [None] * len(trigger_filenames)
    mmia_raw = [None] * len(trigger_filenames)

    for i in range(len(trigger_filenames)):

        # Creating the template for trigger data and info
        triggers1 = [None] * len(trigger_filenames[i])
        triggers2 = [None] * len(trigger_filenames[i])
        trigger_limits[i] = triggers1
        mmia_raw[i] = triggers2

    for i in range(len(data_files)):

        # Filling mmia raw data and trigger info variables

        data_date = data_files[i][0:8]
        
        if len(data_files[i]) == 19: # Trigger number is 1 digit
            data_trigger = data_files[i][9:10]
        else:
            data_trigger = data_files[i][9:11]

        data_mat = sio.loadmat(mmia_mat_files_path+'/'+data_date+'_' + data_trigger +'_data.mat')
        info_mat = sio.loadmat(mmia_mat_files_path+'/'+data_date+'_' + data_trigger +'_info.mat')
        current_data = data_mat.get('MMIA_all')
        current_info = info_mat.get('space_time')
        mmia_raw[matches.index(data_date)][int(data_trigger)] = current_data
        trigger_limits[matches.index(data_date)][int(data_trigger)] = current_info

    return [mmia_raw, trigger_limits]


[matches, trigger_filenames] = get_MMIA_triggers(path_to_mmia_files, trigger_length)

if pre_extracted_MMIA == False:
    
    create_MMIA_trigger_directories(matches, trigger_filenames, path_to_mmia_files, ssd_path)

    [mmia_raw, trigger_limits] = extract_trigger_info(ssd_path, trigger_filenames, matches)

else:
    [mmia_raw, trigger_limits] = upload_MMIA_mats(ssd_path, trigger_filenames, matches)
    
if pre_downloaded_GLM == False:
    download_GLM(ssd_path, trigger_filenames, mmia_raw)
