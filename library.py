'''
Library of necessary functions to run 'main.py' to download GLM data, condition
MMIA and GLM data and compare it

Please note every 'trigger' notation is an 'event' if existing

Morán Domínguez, Jaime Francisco
jaime.francisco.moran@upc.edu
'''

import matplotlib.pyplot as plt
import numpy as np
import math
import os
import re
import netCDF4 as nc
from netCDF4 import Dataset
import numpy.matlib
import scipy.io as sio
import datetime
import pandas as pd
from scipy.signal import find_peaks
from scipy.signal import correlate
from google.cloud import storage
import pickle
import library as TFG


# MMIA data extraction, upload and conditioning

def get_MMIA_events(path_to_mmia_files, event_length):
    '''
    This function gets the MMIA .cdf files and separates them in events according
    to their time specified in the name of the .cdf file.

    Parameters
    ----------
    path_to_mmia_files : str
        Path to where all MMIA .cdf files are located.
    event_length : int, double
        Maximum length (in seconds) to consider two files to be part of the same
        event.

    Returns
    -------
    matches : list
        List containing all dates with existing MMIA data in the form YYMMDD.
    event_info : list
        List containing the names of all .cdf files that correspond to a single
        event per position.

    '''

    print('Generating events from MMIA files...')
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

    event_info = []
    for i in range(len(matches)):
        event_info.append([])

    day = datetime.datetime.strptime(str(mmia_data[0,1])[0:7], "%Y%j").strftime("%Y%m%d")

    for i in range(len(mmia_files)):

        if i == 0:  # First file
            pivoting_time = mmia_data[0,1]
            file_list = []
            file_list.append(mmia_files[0])

        else:   #All of the rest of files

            if mmia_data[i,1] - pivoting_time <= event_length:   # Same event
                file_list.append(mmia_files[i])
                day = datetime.datetime.strptime(str(mmia_data[i,1])[0:7], "%Y%j").strftime("%Y%m%d")
                if i == (len(mmia_files)-1):
                    event_info[matches.index(day)].append(file_list)

            else:   # New trigger
                event_info[matches.index(day)].append(file_list)
                pivoting_time = mmia_data[i,1]
                file_list = []
                file_list.append(mmia_files[i])
                day = datetime.datetime.strptime(str(mmia_data[i,1])[0:7], "%Y%j").strftime("%Y%m%d")

    print('Event generation done!\n')
    return [matches, event_info]

def create_MMIA_event_directories(matches, event_filenames, path_to_mmia_files, ssd_path):
    '''
    This function creates a directory in 'ssd_path/mmia_dirs' for every generated
    MMIA event and drops the corresponding .cdf files in there.

    Parameters
    ----------
    matches : list
        List containing all dates with existing MMIA data in the form YYMMDD.
    event_filenames : list
        List containing the names of all .cdf files that correspond to a single
        event per position.
    path_to_mmia_files : str
        Path to where all MMIA .cdf files are located.
    ssd_path : str
        Path to where all outputted data will be stored.

    Returns
    -------
    None.

    '''
    
    print('Copying MMIA files into event directories...')
    os.system('mkdir ' + ssd_path + '/mmia_dirs')
    for i in range(len(event_filenames)):
        for j in range(len(event_filenames[i])):
            # Creating the trigger directory
            os.system('mkdir ' + ssd_path + '/mmia_dirs/' + matches[i] + '_' + str(j))
            for k in range(len(event_filenames[i][j])):
                file_name = event_filenames[i][j][k]
                # Moving the current file to its trigger directory
                os.system('cp ' + path_to_mmia_files + '/' + file_name + ' ' + ssd_path + '/mmia_dirs/' + matches[i] + '_' + str(j))

    print('MMIA classification done!\n')

def extract_MMIA_event_info(path_to_mmia_dirs, mmia_mats_files_path, matlab_path):
    '''
    Function that extracts MMIA data from their .cdf files per event using the 
    MATLAB script 'MMIA_extraction.m' developed by Jesús López.

    Parameters
    ----------
    path_to_mmia_dirs : str
        Path to where all pre-created MMIA event directories are located
        ('ssd_path/mmia_dirs')
    mmia_mats_files_path : str
        Path to where all resulting .mat files from the extraction will be
        located ('ssd_path/mmia_mat')
    matlab_path : str
        Path to where the MATLAB binary is located.

    Returns
    -------
    None.

    '''
    
    # Extracting data for every trigger

    print('Starting MatLab and extracting data from .cdf files (this process can take a while)...')
    
    # Writing those paths MATLAB will use in a .txt
    os.system("echo '" + path_to_mmia_dirs + "'/ > mmia_dirs_path.txt")
    os.system("echo '" + mmia_mats_files_path + "'/ > mmia_mats_path.txt")
    
    # Executing MATLAB
    my_string = '\"' + matlab_path + '\" -nojvm -nodisplay -nosplash -nodesktop -r \"MMIA_extraction; quit;\"'
    os.system(my_string)
    
    # Removing .txt's with paths (no longer used)
    os.system('rm mmia_dirs_path.txt mmia_mats_path.txt')
    
    with os.scandir(mmia_mats_files_path) as files:
        files = [file.name for file in files if file.is_file() and file.name.endswith('.mat')]

    if len(files) == 0:
        print('No MMIA data could be extracted for any day!')

def upload_MMIA_mats(ssd_path, event_filenames, matches, current_day):
    '''
    This function uploads all MMIA .mat files that correspond to a same day
    as a list in order to further analyze the signals.

    Parameters
    ----------
    ssd_path : str
        Path to where all outputted data will be stored.
    event_filenames : list
        List containing the names of all .cdf files that correspond to a single
        event per position.
    matches : list
        List containing all dates with existing MMIA data in the form YYMMDD.
    current_day : int
        Current iteration day as the current index of 'matches'.

    Returns
    -------
    mmia_raw : list
        List of n positions representing n events for a single day with MMIA
        data. Contains time vector and 3 photometer signals.
    event_limits : list
        List containing min and max times, latitudes and longitudes for every
        event.

    '''

    mmia_mat_files_path = ssd_path + '/mmia_mat'

    with os.scandir(mmia_mat_files_path) as data_files:
        data_files = [file.name for file in data_files if file.is_file() and file.name.endswith('data.mat')]

    with os.scandir(mmia_mat_files_path) as info_files:
        info_files = [file.name for file in info_files if file.is_file() and file.name.endswith('info.mat')]

    event_limits = [None] * len(event_filenames[current_day])
    mmia_raw = [None] * len(event_filenames[current_day])

    for i in range(len(data_files)):

        # Filling mmia raw data and trigger info variables

        data_date = data_files[i][0:8]
        
        if data_date == matches[current_day]:
        
            if len(data_files[i]) == 19: # Event number is 1 digit
                data_trigger = data_files[i][9:10]
            elif len(data_files[i]) == 20:  # Event number is 2 digits
                data_trigger = data_files[i][9:11]
            elif len(data_files[i]) == 21:  # Event number is 3 digits (not expected but sometimes)
                data_trigger = data_files[i][9:12]

            data_mat = sio.loadmat(mmia_mat_files_path+'/'+data_date+'_' + data_trigger +'_data.mat')
            info_mat = sio.loadmat(mmia_mat_files_path+'/'+data_date+'_' + data_trigger +'_info.mat')
            current_data = data_mat.get('MMIA_all')
            current_info = info_mat.get('space_time')
            mmia_raw[int(data_trigger)] = current_data
            event_limits[int(data_trigger)] = current_info

    return [mmia_raw, event_limits]

def split_MMIA_events(mmia_raw, event_limits, event_filenames_on_day, split_window):
    '''
    This function detects events that actually correspond to 2 or more events.
    This can happen when more than 1 event fits inside 'event_length' seconds
    (see function 'get_MMIA_events()')

    Parameters
    ----------
    mmia_raw : list
        List of n positions representing n events for a single day with MMIA
        data. Contains time vector and 3 photometer signals.
    event_limits : list
        List containing min and max times, latitudes and longitudes for every
        event.
    event_filenames_on_day : list
        Variable 'event_filenames' at position 'day' of the daily iterations of 
        'matches'. That is, the names of the .cdf files that conform all events
        for a single day ordered by event.
    split_window : int, double
        Minimum time of signal separation to consider as two different events.

    Returns
    -------
    Updated 'mmia_raw', 'event_limits' and 'event_filenames_on_day' variables
    accounting for new events that may have been created (new positions).

    '''
    
    print('Looking for splittable events...')
    for i in range(len(mmia_raw)):
        
        if type(mmia_raw[i]) == np.ndarray:
        
            event_times = mmia_raw[i][:,0]
            
            split_positions = []
            
            for j in range(1, len(event_times)):
                if (event_times[j] - event_times[j-1]) >= split_window:
                    split_positions.append(j)
            
            if len(split_positions) != 0:   # If there were time jumps
                
                if len(split_positions) == 1:
                    
                    current_mmia_raw_copy = mmia_raw[i]
                    mmia_raw[i] = mmia_raw[i][0:split_positions[0]-2,:]
                    mmia_raw.append(current_mmia_raw_copy[split_positions[0]+2:-1,:])
                    event_limits.append(event_limits[i])
                    event_filenames_on_day.append(event_filenames_on_day[i])
                
                else:
                    
                    current_mmia_raw_copy = mmia_raw[i]
                    split_positions.append(-1)
            
                    for j in range(len(split_positions)):

                        if j == 0:
                            mmia_raw[i] = mmia_raw[i][0:split_positions[0]-2,:]
                        elif j == (len(split_positions)-1):
                            mmia_raw.append(current_mmia_raw_copy[split_positions[j-1]+2:-1,:])
                            event_limits.append(event_limits[i])
                            event_filenames_on_day.append(event_filenames_on_day[i])
                        else:
                            mmia_raw.append(current_mmia_raw_copy[split_positions[j-1]+2:split_positions[j]-2,:])
                            event_limits.append(event_limits[i])
                            event_filenames_on_day.append(event_filenames_on_day[i])
        
    print('Done!\n')
    return [mmia_raw, event_limits, event_filenames_on_day]

def merge_days(matches, bin_path, filename_end_to_merge):
    '''
    Auxiliary function to create a matrix of any daily-list type variable, where
    every position corresponds to a single day.

    Parameters
    ----------
    matches : list
        List containing all dates with existing MMIA data in the form YYMMDD.
    bin_path : str
        Path to where the binaries for the selected variable to merge are located.
    filename_end_to_merge : str
        Possible values = ['no_type', str], for binaries special nomenclature.

    Returns
    -------
    matrix : list
        List of daily lists where every position corresponds to the selected
        variable as lists of daily events.

    '''
    
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

def condition_MMIA_data(MMIA_data, matches, show_plots, mmia_threshold, current_day):
    '''
    This functions takes 'MMIA_data', a list of MMIA tables of information
    and applies a filter in 777.4nm photometer information vector to reduce
    noise.
    --> Revision: no longer filtering MMIA since introduction of TCE conversion.
    It also plots every signal with and without the filter applied
    if 'show_plots' is True.

    Parameters
    ----------
    mmia_data : list
        List of n positions representing n events for a single day with MMIA
        data. Contains time vector and 3 photometer signals.
    matches : list
        List containing all dates with existing MMIA data in the form YYMMDD.
    show_plots : bool
        Boolean variable for plotting throughout the program.
    mmia_threshold : double
        Maximum value of MMIA signal to be considered noise. If all signal is
        below this value, will be considered as just noise.
    current_day : int
        Current iteration day as the current index of 'matches'.

    Returns
    -------
    MMIA_filtered : list
        A list of MMIA tables of information with a filter applied
        and only regarding time and 777.4nm photometer information
    '''

    # Creation of the new list of lists of MMIA data with Moving Average
    MMIA_filtered = [None] * len(MMIA_data)

    for j in range(len(MMIA_data)): # For every day of MMIA data

        if type(MMIA_data[j]) == np.ndarray:
            
            print('Conditioning MMIA signal, date %d snippet %d / %d...' % (int(matches[current_day]), j, len(MMIA_data)))
            
            current_data=np.zeros((len(MMIA_data[j]),2))
            current_data[:,0] = MMIA_data[j][:,0]
            current_data[:,1] = MMIA_data[j][:,3] # 4th column as 777.4 nm

            if show_plots == 1:
                plt.figure(figsize=(9, 3))
                figname = matches[current_day] + '_' + str(j) + '.pdf'
                plt.plot(current_data[:,0],current_data[:,1],linewidth=0.5, color='r')
                plt.savefig(figname)
                # Clear the current axes
                plt.cla() 
                # Clear the current figure
                plt.clf() 
                # Closes all the figure windows
                plt.close('all')
            
            if (current_data[:,1] < mmia_threshold).all() == True:
                MMIA_filtered[j] = None
                print('MMIA data for day %d snippet %d was just noise!' % (int(matches[current_day]), j))
                print(' ')
                
            else:

                # Assuring continuity in MMIA_filtered timesteps
                MMIA_filtered[j] = TFG.fit_vector_in_MMIA_timesteps(current_data, int(matches[current_day]), j, show_plots, True)
                
                if show_plots == 1:
                    # MMIA representation with filter and with rectified time
                    plt.figure()
                    plt.plot(MMIA_filtered[j][:,0],MMIA_filtered[j][:,1],linewidth=0.5, color='b')
                    plt.title("Untreated (red) and filtered (blue) MMIA 777.4nm photometer detections of day %d snippet %d" % (int(matches[current_day]), j))
                    plt.xlabel('Time [s]')
                    plt.grid('on')
                    plt.ylabel(r"Irradiance $\left[\dfrac{\mu W}{m^2}\right]$")
                    plt.legend(['Untreated signal', 'Filtered signal'])
                    plt.show()
                    # Clear the current axes
                    plt.cla() 
                    # Clear the current figure
                    plt.clf() 
                    # Closes all the figure windows
                    plt.close('all')
                    
            # Freeing memory
            del current_data
                
        else:
            print('There is no MMIA data for day %d, snippet %d' % (int(matches[current_day]), j))
            print(' ')

    print('Done!')
    print(' ')
    return MMIA_filtered

def fit_vector_in_MMIA_timesteps(GLM_int_data, day, event, show_plots, is_MMIA):
    '''
    This function takes a time and signal pair of vectors and accommodates it
    into MMIA timesteps of 0.00001s. Inexistent values in between are filled
    by simple linear regression.

    Parameters
    ----------
    GLM_int_data : list
        Pair of time-signal vectors of an event. IMPORANT: NOT necessarily GLM data.
    day : str
        'matches[day]', where 'day' is the current daily iteration of 'matches'.
    event : int
        Index of the current event inside the daily vector of data.
    show_plots : bool
        Boolean variable for showing plots all through the program.
    is_MMIA : bool
        Boolean variable for sepparating GLM expansion from MMIA expansion.

    Returns
    -------
    GLM_current_data : list
        Pair of time-signal vectors of an event like "GLM_int_data" input, but with
        accomodation to 0.00001s timesteps done. IMPORANT: NOT necessarily GLM data.
    '''

    # Expanding event to fit missing MMIA timesteps to cross-correlate data
    if is_MMIA == 0:
        print('Fitting GLM data in MMIA timesteps date %d event %d...' % (day,event))
    else:
        print('Completing MMIA data in MMIA timesteps date %d event %d...' % (day,event))

    new_length = 1          # New length of the timestep-wise matrix
    acumulated_voids = 0    # Number of non-existing timesteps up to current line
    void_info = np.zeros((len(GLM_int_data),4)) # Matrix of special info for each line:
        # 1st column: .txt row number
        # 2nd column: void timesteps after that row until next existing timestep
        # 3rd column: Accumulated void timesteps before current line
        # 4th column: Differential energy between existing timesteps ([i]-[i-1])

    # Updating new_length value to make a new table with 1st dimension being new_length

    GLM_int_data[0,0] = round(GLM_int_data[0,0],5)     # Rounding to MMIA period

    for j in range(1,len(GLM_int_data)):

        GLM_int_data[j,0] = round(GLM_int_data[j,0],5) # Rounding to MMIA period
        void_info[j][0] = j    # Filling first void_info column

        if GLM_int_data[j,0] == GLM_int_data[j-1,0] + 0.00001:  # Exactly one timestep ahead
            new_length = new_length + 1
            void_info[j][2] = acumulated_voids

        elif GLM_int_data[j,0] < GLM_int_data[j-1,0] + 0.00001: # Less than a whole timestep (sometimes occur)
            new_length = new_length + 1
            void_info[j][2] = acumulated_voids

        else:   # There are missing timesteps in between current and last row
            void_timesteps = round((GLM_int_data[j,0] - GLM_int_data[j-1,0])/0.00001) - 1
            new_length = new_length + 1 + void_timesteps
            void_info[j-1][1] = void_timesteps
            acumulated_voids = acumulated_voids + void_timesteps
            void_info[j][2] = acumulated_voids
            void_info[j-1][3] = GLM_int_data[j,1] - GLM_int_data[j-1,1]

    # Filling the new time-wise matrix

    GLM_current_data = np.zeros((new_length,2)) # New matrix with void lines for non-existing timesteps

    for j in range(0,len(GLM_int_data)):
        new_j = int(j + void_info[j,2])    # Row position in the new matrix
        GLM_current_data[new_j,:] = GLM_int_data[j,:] # Filling rows with existing data

        if void_info[j,1] != 0:   # Lines with non-existing timesteps afterwards
            counter = 1       # Adds 0.00001s and a linear energy fraction
            for k in range(new_j+1, new_j+1+int(void_info[j,1])):
                GLM_current_data[k,0] = GLM_int_data[j,0] + counter * 0.00001
                GLM_current_data[k,1] = GLM_int_data[j][1] + counter * (void_info[j][3]/void_info[j][1])
                counter = counter + 1

    if show_plots == 1 and is_MMIA == 0:
        # GLM time representation at MMIA sample rate
        plt.figure()
        plt.plot(GLM_current_data[:,0])
        #plt.title('GLM Time vector of day %d event %d with 0.00001s period' % (day, event))
        plt.xlabel('Samples')
        plt.ylabel('Time [s]')
        plt.grid('on')
        plt.show()
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')

        # Radiance vs time graph representation
        plt.figure()
        plt.plot(GLM_current_data[:,0],GLM_current_data[:,1], linewidth=0.5, color='black')
        plt.grid('on')
        #plt.title('GLM signal of day %d event %d with MMIA sample rate (0.00001s)' % (day, event))
        plt.xlabel('Time (second of the day) [s]')
        plt.ylabel('Radiance [J]')
        plt.show()
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')

    print('Date %d event %d fit' % (day, event))
    print(' ')

    return GLM_current_data


# GLM data download, extraction, upload and conditioning

def download_GLM(ssd_path, event_filenames, MMIA_filtered, matches, current_day):
    '''
    This function sets the variables for downloading GLM .nc files directly
    from the Google Cloud Storage servers.

    Parameters
    ----------
    ssd_path : str
        Path to where all outputted data will be stored.
    event_filenames : list
        List containing the names of all .cdf files that correspond to a single
        event per position.
    MMIA_filtered : list
        A list of MMIA tables of information with a filter applied
        and only regarding time and 777.4nm photometer information
    matches : list
        List containing all dates with existing MMIA data in the form YYMMDD.
    current_day : int
        Current iteration day as the current index of 'matches'.

    Returns
    -------
    None. Downloaded GLM .nc files in 'ssd_path/glm_downl_nc_files'.

    '''

    print("Downloading GLM's .nc files from Google Cloud Storage for day %s...\n" % matches[current_day])

    path_to_downloaded_glm_nc = ssd_path + '/glm_downl_nc_files'
    
    if current_day == 0:
        os.system('mkdir ' + path_to_downloaded_glm_nc)

    for j in range(len(event_filenames[current_day])):
        if type(MMIA_filtered[j]) == np.ndarray:
            current_trigger_download_path = path_to_downloaded_glm_nc + '/' + matches[current_day] + '_' + str(j)
            os.system('mkdir ' + current_trigger_download_path)

            date = event_filenames[current_day][j][0][50:60]
            year = date[0:4]
            day_year = datetime.datetime.strptime(date, "%Y-%m-%d").strftime("%j")

            download_glm_from_google(ssd_path, year, day_year, MMIA_filtered[j][0,0], MMIA_filtered[j][-1,0])

            with os.scandir(ssd_path) as files:
                files = [file.name for file in files if file.is_file() and file.name.endswith('.nc')]

            for k in range(len(files)):
                os.system('mv ' + ssd_path + '/' + files[k] + ' ' + current_trigger_download_path)    

    print('Download done! Your .nc files for day %s can be accessed per trigger at %s' % (matches[current_day], path_to_downloaded_glm_nc))
    print(' ')

def download_glm_from_google(ssd_path, year, day_year, t_ini, t_end):
    '''
    This function actually downloads GLM .nc files from Google Cloud Services
    servers. It downloads the whole second of data (i.e. 3 files) as the minimum
    time.
    Modified from original by Jesús López.

    Parameters
    ----------
    ssd_path : str
        Path to where all outputted data will be stored.
    year : str
        Year of the data willing to download.
    day_year : str
        Day of the year of the data willing to download (i.e. 1-365).
    t_ini : double
        First time of the corresponding MMIA signal for that event.
    t_end : double
        Last time of the corresponding MMIA signal for that event.

    Returns
    -------
    None. Downloaded GLM .nc files in 'ssd_path/glm_downl_nc_files'.

    '''

    hour_ini = str(int(t_ini // 3600))
    if len(hour_ini) == 1:
        hour_ini = '0'+hour_ini
    min_ini = "%02d" % ((t_ini/3600 - t_ini//3600) * 60)
    seg_ini = str(int(((t_ini/3600 - t_ini//3600) * 60 - int((t_ini/3600 - t_ini//3600) * 60)) * 60))
    mseg_ini = '000'


    hour_end = str(int(t_end // 3600))
    if len(hour_end) == 1:
        hour_end = '0'+hour_end
    #min_end = str(int((t_end/3600 - t_end//3600) * 60))
    min_end = "%02d" % ((t_end/3600 - t_end//3600) * 60)
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
        print("Blob {} downloaded to {}".format(file, destination_file_name))

def extract_GLM(dir_path, output_path, event_limits, matches, MMIA_filtered, angle_margin, cropping_margin, current_day):
    '''
    This function calls every directory with .nc files and extracts the data
    of all the files in it via GLM_processing function.

    Parameters
    ----------
    dir_path : string
        Path to the directory where the daily-ordered directories with
        ordered separated .nc files are located.
    output_path : string
        Path to the existing directory where the resulting daily .txt files
        will be located.
    event_limits : list
        List containing min and max times, latitudes and longitudes for every
        event.
    matches : list
        List containing all dates with existing MMIA data in the form YYMMDD.
    MMIA_filtered : list
        A list of MMIA tables of information with a filter applied
        and only regarding time and 777.4nm photometer information
    angle_margin : float
        Plus of latitude and longitude angle for extracting GLM data with
        respect to MMIA's data.
    cropping_margin : float
        Plus of time before and after MMIA event times for extracting GLM data.
    current_day : int
        Current iteration day as the current index of 'matches'.

    Returns
    -------
    .txt files for every event with GLM data prepared to be analyzed.

    '''

    print('Extracting data from GLM .nc files into event .txt...')
    print(' ')
    
    if current_day == 0:
        # Creating a directory to store all resulting .txt
        os.system('mkdir ' + output_path)
    
    for j in range(len(event_limits)):  # Analyzing each directory's .nc files

        if type(MMIA_filtered[j]) == np.ndarray and type(event_limits[j]) == np.ndarray:

            print('Extracting GLM data for date %d event %d...' % (int(matches[current_day]), j))
            
            trigger_name = matches[current_day] + '_' + str(j)
            trigger_path = dir_path + '/' + trigger_name
            
            with os.scandir(trigger_path) as files:
                files = [file.name for file in files if file.is_file() and file.name.endswith('.nc')]
            
            if len(files) == 0:
                print('No GLM .nc files could be downloaded for day %s, event %d\n' % (matches[current_day], j))
            else:

                min_lat = event_limits[j][0,0] - angle_margin
                max_lat = event_limits[j][0,1] + angle_margin

                min_lon = event_limits[j][0,2] - angle_margin
                max_lon = event_limits[j][0,3] + angle_margin

                start_time = MMIA_filtered[j][0,0] - cropping_margin
                end_time = MMIA_filtered[j][-1,0] + cropping_margin

                #start_time = event_limits[j][0,4] - cropping_margin
                #end_time = event_limits[j][0,5] + cropping_margin

                TFG.GLM_processing(dir_path+'/'+trigger_name+'/', output_path, trigger_name, min_lat, max_lat, min_lon, max_lon, start_time, end_time)
                
                print(' ')
                print('Date %s event %d done\n' % (matches[current_day], j))
        else:
            print('GLM data for date %d event %d will not be extracted due to lack of MMIA data\n' % (int(matches[current_day]), j))

    print('Your processed .txt files for day %s can be accessed at %s' % (matches[current_day], output_path))
    print(' ')

def GLM_processing(read_path, save_path, name, min_lat, max_lat, min_lon, max_lon, begin, end):
    '''
    This function takes a directory with GLM .nc files and extracts their
    important data inside space and time windows into a .txt file.
    Modified from original by Jesús López.

    Parameters
    ----------
    read_path : str
        Path to the directory with GLM .nc files for a single event
    save_path : str
        Path to the existing directory where the resulting daily .txt files
        will be located.
    name : str
        Name for the generated .txt file (default is event date and number)
    min_lat : float
        Minimum latitude to take data from.
    max_lat : float
        Maximum latitude to take data from.
    min_lon : float
        Minimum longitude to take data from.
    max_lon : float
        Maximum longitude to take data from.
    begin : float
        Start second of the day to take data from (default MMIA signal starting
                                                   time minus cropping_margin)
    end : float
        Final second of the day to take data from (default MMIA signal starting
                                                   time plus cropping_margin)

    Returns
    -------
    .txt files for every event with GLM data prepared to be analyzed.

    '''

    length=len([f for f in os.listdir(read_path) if os.path.isfile(os.path.join(read_path, f))])
    files=os.listdir(read_path); 
    #%%
    for ind in range(0,length,1):

            file1 = read_path+files[ind]
            path = files[ind]
            total = len(os.listdir(read_path))
                
                ############ FIX UNSIGNED ATTRIBUTE ##############
            print('GLM file '+str(ind+1)+'/'+str(total))
            print('GLM file '+ file1)    
                
    #%%        
            ############ FLASH ID ##############
            # CUANDO SE ALMACENA LOS FICHEROS CON NOMBRES ORIGINALES        
            Start = (path[path.find("_s")+1:path.find("_e")])
            Year = int(Start[1:5])
            Day = int(Start[5:8])
            Seconds = int(Start[8:10])*3600+int(Start[10:12])*60+int(Start[12:14])
            
            
            data = datetime.datetime(Year, 1, 1) + datetime.timedelta(Day - 1)
            data_int = int(data.strftime("%Y%m%d"))

            g16glm = Dataset(file1,'r')       
            event_time = g16glm.variables['event_time_offset'][:]
            
            a = len(event_time)
            if a != 0:         
                event_time_2=np.matlib.zeros((a, 1))
                k2 = 0
                if event_time[a-1]>20: # Cada fichero cotiene hasta 20 segundos?
                    for y in event_time:
                        if y>30000:
                            y = y-131072
                        event_time_2[k2] = y/1000
                        k2 = k2 +1
                else:    
                    event_time_2 = event_time # Si no hay más eventos fuera del rango de los 20 segundos, no hay necesidad ser ajustados. 
            
            time_seconds = Seconds + event_time_2 # Segundos del día.
            
            # DF de los datos! 
            # Inicio con DF temporal de las coordenadas de los flashes
            
            df_latlon_flash=pd.DataFrame(g16glm.variables['flash_id'][:],columns=['flash_id'])
            df_latlon_flash['flash_lat']=g16glm.variables['flash_lat'][:]
            df_latlon_flash['flash_lon']=g16glm.variables['flash_lon'][:]
            df_latlon_flash['flash_energy']=g16glm.variables['flash_energy'][:]

            # DF de los flashes para hacer la búsqueda y asignar a los eventos el ID del flash! 
                
            df_glm_flash=pd.DataFrame(g16glm.variables['group_id'][:],columns=['group_id'])
            df_glm_flash['group_pfid']=g16glm.variables['group_parent_flash_id'][:]
            
            df_glm_flash['flash_lat']=df_glm_flash['group_pfid'].map(df_latlon_flash.set_index('flash_id')['flash_lat'])
            df_glm_flash['flash_lon']=df_glm_flash['group_pfid'].map(df_latlon_flash.set_index('flash_id')['flash_lon'])
            df_glm_flash['flash_energy']=df_glm_flash['group_pfid'].map(df_latlon_flash.set_index('flash_id')['flash_energy'])
            
            # DF de eventos
                
            df_glm=pd.DataFrame(time_seconds,columns=['seconds_day'])
            df_glm['event_lat']=g16glm.variables['event_lat'][:]
            df_glm['event_lon']=g16glm.variables['event_lon'][:]    
            ##### BÚSQUEDA DIRECTA DE LOS ID DE FLASHES DE SUS CORRESPONDIENTES EVENTOS
            df_glm['event_pgid']=g16glm.variables['event_parent_group_id'][:]
            df_glm['flash_id']=df_glm['event_pgid'].map(df_glm_flash.set_index('group_id')['group_pfid'])
            df_glm['flash_lat']=df_glm['event_pgid'].map(df_glm_flash.set_index('group_id')['flash_lat'])
            df_glm['flash_lon']=df_glm['event_pgid'].map(df_glm_flash.set_index('group_id')['flash_lon'])
            df_glm['event_energy']=g16glm.variables['event_energy'][:]     
            df_glm['flash_energy']=df_glm['event_pgid'].map(df_glm_flash.set_index('group_id')['flash_energy'])
            
            # Elimino la columna del event_parent_group_id
            df_glm.drop("event_pgid", axis=1, inplace=True)
            
            df_glm_in = df_glm[(df_glm['seconds_day'] >= begin) & (df_glm['seconds_day'] <= end) &
                               (df_glm['event_lat'] >= min_lat) & (df_glm['event_lat'] <= max_lat) &
                               (df_glm['event_lon'] >= min_lon) & (df_glm['event_lon'] <= max_lon)]
            
            if ind==0:
                
                df_export=df_glm_in; 
            
            else:
                df_export=pd.concat([df_export,df_glm_in])
            
            g16glm.close()

    df_export.to_csv(save_path+'/'+name+'.txt', sep=' ', index=False,header=False)

def unify_GLM_data(output_path, MMIA_filtered, matches, current_day):
    '''
    This function gets all the GLM's extracted data .txt files from the
    directory output_path and creates and returns list GLM_raw_data as a variable.

    Parameters
    ----------
    output_path : string
        Path to the existing directory where the resulting daily .txt files
        will be located.
    MMIA_filtered : list
        A list of MMIA tables of information with a filter applied
        and only regarding time and 777.4nm photometer information
    matches : list
        List containing all dates with existing MMIA data in the form YYMMDD.
    current_day : int
        Current iteration day as the current index of 'matches'.

    Returns
    -------
    GLM_raw_data : list
        A list with all the data found in the .txt files, ordered by events for
        a single day.
    '''

    print('Uploading GLM data from .txt files...')
    # Creation of a new list of daily GLM events
    GLM_raw_data = [None] * len(MMIA_filtered)

    with os.scandir(output_path) as files:
        files = [file.name for file in files if file.is_file() and file.name.endswith('.txt')]
    size = len(files)
    if size == 0:
        print('No GLM .txt files found!')

    column_subset = ['Time', 'Event_lat', 'Event_lon', 'Event_ID', 'Flash_lat', 'Flash_lon', 'Event_radiance', 'Flash_radiance']

    for i in range(size): # For every .txt file (for every event)

        day = files[i][0:8]
        
        if day == matches[current_day]:
            
            if len(files[i]) == 14: # Trigger number is 1 digit
                event = int(files[i][9:10])
            elif len(files[i]) == 15:   # Trigger number is 2 digits
                event = int(files[i][9:11])
            elif len(files[i]) == 16:   # Trigger number is 3 digits (not expected)
                event = int(files[i][9:12])

            current_path = output_path + '/' + files[i]

            # Uploading the GLM's .txt using Pandas to sort it by time
            #current_data = pd.read_csv(current_path, names=column_subset, sep='\s+')
            current_data = pd.read_csv(current_path, names=column_subset, sep=' ')
            # Sorting data by time
            current_data = current_data.sort_values(by='Time')
            # Translating Pandas Dataframe to Numpy Matrix for easy data access
            current_data = current_data.to_numpy()
            # Appending current day to GLM_raw_data
            GLM_raw_data[event] = current_data
            # Freeing memory
            del current_data

    print('Done!')
    print(' ')

    return GLM_raw_data

def condition_GLM_data(GLM_total_raw_data, matches, show_plots, current_day):
    '''
    This function takes all the extracted data from GLM .txt files and prepares
    it for further analysis, as well as filtering those events with not
    enough data.

    Parameters
    ----------
    GLM_total_raw_data : list
        A list with all the data found in the .txt files, ordered by events for
        a single day.
    matches : list
        List containing all dates with existing MMIA data in the form YYMMDD.
    show_plots : bool
        Boolean variable for plotting throughout the program.
    current_day : int
        Current iteration day as the current index of 'matches'.

    Returns
    -------
    GLM_data : list
        A list with all the conditioned data found in GLM_total_raw_data,
        ordered by events for a single day.
    '''

    # Creating a new set of data
    GLM_data = [None] * len(GLM_total_raw_data)

    # Integrating and extending GLM vectors by date
    for j in range(len(GLM_total_raw_data)):   # For every event with GLM data
        
        # If the resulting .txt/.csv has 0 or 1 line of data
        if type(GLM_total_raw_data[j]) == np.ndarray and len(GLM_total_raw_data[j]) <= 1:
            print('GLM detection for day %d event %d is void!' % (int(matches[current_day]), j))
            print(' ')
            GLM_total_raw_data[j] = None

        # If not, check if only one timestep (still not enough data)
        elif type(GLM_total_raw_data[j]) == np.ndarray and len(GLM_total_raw_data[j]) != 0:
            just_one_timestep = 1
            check_pos = 1
            while just_one_timestep == 1:
                # If last position and same as timestep before in the .txt
                if check_pos == (len(GLM_total_raw_data[j])-1) and GLM_total_raw_data[j][check_pos-1,0] == GLM_total_raw_data[j][check_pos,0]:
                    just_one_timestep = 2
                else:
                    # Different timestep in the .txt as in line before
                    if GLM_total_raw_data[j][check_pos-1,0] != GLM_total_raw_data[j][check_pos,0]:
                        just_one_timestep = 0
                    else: # Same timestep in the .txt as in line before
                        check_pos = check_pos + 1

        if type(GLM_total_raw_data[j]) == np.ndarray and just_one_timestep == 2:
            print('GLM detection for day %d event %d contains only 1 timestep and will not be compared' % (int(matches[current_day]), j))
            print(' ')
            GLM_total_raw_data[j] = None

        if type(GLM_total_raw_data[j]) == np.ndarray:
            
            # Showing non-inflated time vector
            if show_plots == 1:
                plt.figure()
                plt.plot(GLM_total_raw_data[j][:,0])
                plt.xlabel('Samples')
                plt.ylabel('Time [s]')
                plt.title('Original GLM Time VS Samples for date %d event %d' % (int(matches[current_day]), j))
                plt.grid('on')
                plt.show()
                # Clear the current axes
                plt.cla() 
                # Clear the current figure
                plt.clf() 
                # Closes all the figure windows
                plt.close('all')


            GLM_data[j] = GLM_total_raw_data[j][:,[0,6]] # Time and instrument signal for event radiance
            
            # Check for too short snippet vectors
            if len(GLM_data[j])<=5:
                print('Data for day %s event %d is too poor, only %d samples. This event will be omitted.' % (matches[current_day], j, len(GLM_data[j])))
                GLM_data[j] = None

            if show_plots == 1 and type(GLM_data[j]) == np.ndarray:
                # Plotting lineality in GLM time vector with GLM sampling rate
                plt.figure()
                plt.plot(GLM_data[j][:,0])
                plt.title('GLM Time vector of day %d event %d with 0.002s period' % (int(matches[current_day]), j))
                plt.xlabel('Samples')
                plt.ylabel('Time [s]')
                plt.grid('on')
                plt.show()
                # Clear the current axes
                plt.cla() 
                # Clear the current figure
                plt.clf() 
                # Closes all the figure windows
                plt.close('all')

                # Integrated GLM radiance vs time graph representation
                plt.figure()
                plt.plot(GLM_data[j][:,0],GLM_data[j][:,1], linewidth=0.5, color='black')
                plt.grid('on')
                plt.title('GLM signal of day %d event %d with GLM sample rate (0.002s)' % (int(matches[current_day]), j))
                plt.xlabel('Time [s]')
                plt.ylabel('Radiance [J]')
                plt.show()
                # Clear the current axes
                plt.cla() 
                # Clear the current figure
                plt.clf() 
                # Closes all the figure windows
                plt.close('all')
                
    print('Conditioning of GLM vectors done!')
    print(' ')

    return GLM_data


# Conversion to Top Cloud Energy (TCE)

def top_cloud_energy(GLM_data, MMIA_filtered, current_day, show_plots, tce_figures_path, glm_pix_size, cropping_margin):
    '''
    This function takes all "raw" GLM and MMIA instrument data and applies a
    conversion to Top Cloud Energy to the signals, as per Van der Velde et al,
    2020.

    Parameters
    ----------
    GLM_data : list
        A list with all the conditioned data found in GLM_total_raw_data,
        ordered by events for a single day.
    MMIA_filtered : list
        A list of MMIA tables of information with a filter applied
        and only regarding time and 777.4nm photometer information
    current_day : int
        Current iteration day as the current index of 'matches'.
    show_plots : bool
        Boolean variable for plotting throughout the program.
    tce_figures_path : str
        Path to the directory where all TCE figures will be stored.
    glm_pix_size : float
        Size of the projection of a GLM pixel at the observation area [km^2]
    cropping_margin : float
        Plus of time before and after MMIA event times for extracting GLM data.

    Returns
    -------
    glm_tce : list
        List of TCE-converted GLM events for a single day
    mmia_tce : list
        List of TCE-converted MMIA events for a single day

    '''
    
    glm_tce = [None] * len(GLM_data)
    mmia_tce = [None] * len(MMIA_filtered)
    
    for i in range(len(GLM_data)): # For every event
    
        if type(GLM_data[i]) == np.ndarray and type(MMIA_filtered[i]) == np.ndarray:
            print('Converting instrumental signals to Top Cloud Energy for day %s event %d / %d' % (current_day, i, len(GLM_data)))
            
            # GLM
            GLM_cloud_E = GLM_data[i]
            GLM_cloud_E[:,1] = GLM_data[i][:,1] * 6611570247.933885 * glm_pix_size*1e6 #[J]
            glm_tce[i] = integrate_signal_002(GLM_cloud_E, True, MMIA_filtered[i][0,0]-cropping_margin, MMIA_filtered[i][-1,0]+cropping_margin)
            del GLM_cloud_E


            # MMIA
            # Computing the integral over MMIA signal
            MMIA_cloud_E = integrate_signal_002(MMIA_filtered[i], False, MMIA_filtered[i][0,0]-cropping_margin, MMIA_filtered[i][-1,0]+cropping_margin) # [micro J/m^2]
            MMIA_cloud_E[:,1] = MMIA_cloud_E[:,1]*1e-6*(math.pi)*(400e3**2) #(Van der Velde et al 2020), [J]
            mmia_tce[i] = MMIA_cloud_E
            del MMIA_cloud_E

            if show_plots == True:
                plt.figure(figsize=(10, 6))
                figure_name = current_day + '_' + str(i)
                plt.plot(glm_tce[i][:,0], glm_tce[i][:,1], color='black', linewidth=0.5)
                plt.plot(mmia_tce[i][:,0], mmia_tce[i][:,1], color='red', linewidth=0.5)
                plt.yscale('log')
                plt.legend(['GLM','MMIA'])
                plt.title('GLM (black) and MMIA (red) correlated signals converted to Top Cloud Energy for day %s event %d' % (current_day, i))
                plt.xlabel('Time [s]')
                plt.ylabel('Top Cloud Energy [J]')
                plt.grid('on')
                plt.savefig(tce_figures_path + '/' + figure_name + '_both.pdf')
                #plt.show()
                # Clear the current axes
                plt.cla() 
                # Clear the current figure
                plt.clf() 
                # Closes all the figure windows
                plt.close('all')
                
                plt.figure(figsize=(10, 6))
                figure_name = current_day + '_' + str(i)
                plt.plot(glm_tce[i][:,0], glm_tce[i][:,1], color='black', linewidth=0.5)
                plt.title('GLM correlated signal converted to Top Cloud Energy for day %s event %d' % (current_day, i))
                plt.xlabel('Time [s]')
                plt.ylabel('Top Cloud Energy [J]')
                plt.grid('on')
                plt.savefig(tce_figures_path + '/' + figure_name + '_glm.pdf')
                #plt.show()
                # Clear the current axes
                plt.cla() 
                # Clear the current figure
                plt.clf() 
                # Closes all the figure windows
                plt.close('all')
                
                plt.figure(figsize=(10, 6))
                figure_name = current_day + '_' + str(i)
                plt.plot(mmia_tce[i][:,0], mmia_tce[i][:,1], color='red', linewidth=0.5)
                plt.title('MMIA correlated signal converted to Top Cloud Energy for day %s event %d' % (current_day, i))
                plt.xlabel('Time [s]')
                plt.ylabel('Top Cloud Energy [J]')
                plt.grid('on')
                plt.savefig(tce_figures_path + '/' + figure_name + '_mmia.pdf')
                #plt.show()
                # Clear the current axes
                plt.cla() 
                # Clear the current figure
                plt.clf() 
                # Closes all the figure windows
                plt.close('all')
            
        else:
            print('Signals for day %s event %d could not be correlated, so no conversion is possible' % (current_day, i))
    
    return [glm_tce, mmia_tce]

def integrate_signal_002(event, isGLM, begin, end):
    '''
    This function takes a pair of signal-time vectors and integrates the vector
    in 0.002s windows as for help in TCE conversion.

    Parameters
    ----------
    event : np matrix
        n x 2 matrix where the first vector represents signal values and the
        second represents time values for a GLM/MMIA event
    isGLM : bool
        Bool for selecting integration method for GLM or MMIA
    begin : float
        Starting integration time value.
    end : float
        Ending integration time value.

    Returns
    -------
    int_data : TYPE
        m x 2 matrix where the first vector represents TCE-converted signal
        values and the second represents time values for a GLM/MMIA event

    '''

    new_length = math.ceil((end - begin) / 0.002)
    int_data = np.zeros((new_length, 2))

    for k in range(new_length): # For every sample accounting zeros at GLM rate
        
        # Fill time instance
        int_data[k,0] = round(begin + k*0.002, 3)
        
        # Create windows of 0.002s and integrate their content
        t_min = int_data[k,0]
        t_max = t_min + 0.002
        
        # Check positions in the signal vector which are inside current window
        window_indx = np.where(np.logical_and(event[:,0] >= t_min, event[:,0] < t_max))[0]
        
        if len(window_indx) == 0:
            int_data[k,1] = 1e-11
        
        else: # Data inside this window exists
            
            if isGLM == True: # Event represents GLM signal
                # Simply add values inside window
                int_data[k,1] = sum(event[window_indx,1])
                
            else: # Event represents MMIA signal
                # Trapezoid integration
                int_data[k,1] = np.trapz(event[window_indx,1], x=event[window_indx,0])
                
    return int_data


# Cross-Correlation of GLM and MMIA

def cross_correlate_GLM_MMIA(GLM_snippets, MMIA_snippets, GLM_norm, MMIA_norm, matches, show_plots, current_day, xcorr_figures, GLM_ordered_outputs, tfg, ssd_path):
    '''
    This function takes all GLM and MMIA event data for a single day and
    iterates a cross-correlation method for computing MMIA delay and synchronize
    both signals with GLM as a reference. Usual number of iterations is 1 for
    similar signals.

    Parameters
    ----------
    GLM_snippets : list
        List of TCE-converted GLM events for a single day
    MMIA_snippets : list
        List of TCE-converted MMIA events for a single day
    GLM_norm : list
        List of normalized TCE-converted GLM events for a single day
    MMIA_norm : list
        List of normalized TCE-converted MMIA events for a single day
    matches : list
        List containing all dates with existing MMIA data in the form YYMMDD.
    show_plots : bool
        Boolean variable for plotting throughout the program.
    current_day : int
        Current iteration day as the current index of 'matches'.
    xcorr_figures : str
        Path to the directory where all cross-correlation figures will be stored
    GLM_ordered_outputs : str
        Path to the existing directory where the resulting daily .txt files
        will be located.
    tfg : bool
        Boolean variable for defining TFG events (outputs GLM .txt's that could
                                                  be cross-correlated)
    ssd_path : str
        Path to where all outputted data will be stored.

    Returns
    -------
    GLM_xcorr : list
        List of TCE-converted cross-correlated GLM events for a single day
    MMIA_xcorr : list
        List of TCE-converted cross-correlated MMIA events for a single day
    GLM_xcorr_norm : list
        List of normalized TCE-converted cross-correlated GLM events for a
        single day
    MMIA_xcorr_norm : list
        List of normalized TCE-converted cross-correlated MMIA events for a
        single day
    delays : list
        List of MMIA delays for every event. Positive: MMIA anticipated GLM

    '''

    print('Starting cross-correlation of events and syncronization of signals...')
    print(' ')

    # Creation of new lists for daily GLM and MMIA cross-correlated data
    GLM_xcorr = [None] * len(GLM_snippets)
    MMIA_xcorr = [None] * len(MMIA_snippets)
    GLM_xcorr_norm = [None] * len(GLM_snippets)
    MMIA_xcorr_norm = [None] * len(MMIA_snippets)
    delays = [None] * len(GLM_snippets)
    
    if current_day == 0:
        os.system('mkdir ' + xcorr_figures)

    for j in range(len(GLM_snippets)):
        
        # If there's no info for this snippet due to lack of .nc or .cdf files

        if type(GLM_snippets[j]) == np.ndarray and type(MMIA_snippets[j]) == np.ndarray:
            
            GLM_starts_after_MMIA_ends = (GLM_snippets[j][0,0] >= MMIA_snippets[j][-1,0])
            MMIA_starts_after_GLM_ends = (MMIA_snippets[j][0,0] >= GLM_snippets[j][-1,0])
            
            no_overlap_conditions = [GLM_starts_after_MMIA_ends, MMIA_starts_after_GLM_ends]
            
        if type(GLM_snippets[j]) == np.ndarray and type(MMIA_snippets[j]) == np.ndarray and any(no_overlap_conditions) == False:
            
            current_GLM = GLM_norm[j]
            current_MMIA = MMIA_norm[j]
            
            if show_plots == 1:
                # Plotting and saving non cross-correlated nor syncronized GLM and MMIA normalized signals
                figure_name = matches[current_day] + '_' + str(j)
                plt.figure(figsize=(12, 4))
                plt.plot(current_MMIA[:,0], current_MMIA[:,1], color = 'r', linewidth = 0.5)
                plt.plot(current_GLM[:,0], current_GLM[:,1], color = 'black', linewidth = 0.5)
                plt.legend(['MMIA','GLM'])
                plt.title('GLM (black) and MMIA (red) non-correlated normalized signals for day %d event %d' % (int(matches[current_day]), j))
                plt.xlabel('Time [s]')
                plt.ylabel('Normalized energy')
                plt.grid('on')
                plt.savefig(xcorr_figures + '/' + figure_name + '_no_xc.pdf')
                #plt.show()
                # Clear the current axes
                plt.cla() 
                # Clear the current figure
                plt.clf() 
                # Closes all the figure windows
                plt.close('all')


            # Calculation of delay in samples of GLM with respect to MMIA
            
            #delay = int(TFG.signal_delay(current_GLM, current_MMIA, show_plots, int(matches[current_day]), j))
            
            # Assuring cross-correlation is inside accuracy in time by cropping signals if necessary
            
            delay = 100000  # Very high delay value in order to enter the while loop
            prev_min_delay = delay
            delay_crop = 25 # In samples, or 0.05s at 0.002s in frequency
            max_viable_delay = 40 # In samples, or 0.08s at 0.002s in frequency
            counter = 0
            max_it = 1000
            
            it_current_GLM = current_GLM
            it_current_MMIA = current_MMIA

            while abs(delay) >= max_viable_delay and counter < max_it and len(it_current_GLM) > delay_crop:
                
                delay = int(TFG.signal_delay(it_current_GLM, it_current_MMIA, show_plots, int(matches[current_day]), j))
                
                if counter == 0 or abs(delay) < abs(prev_min_delay):
                    prev_min_delay = delay
                
                if counter != 0 and delay >= prev_min_delay:
                    delay = prev_min_delay
                    counter = max_it # Cut while loop
                
                counter = counter+1
                
                if abs(delay) >= max_viable_delay and counter == 1:
                    
                    # Crop GLM and MMIA first
                    try: # Try to find where in GLM the MMIA signal starts
                        MMIA_start_in_GLM_pos = np.where(current_GLM[:,0] <= current_MMIA[0,0])[0][-1]
                        GLM_start_in_MMIA_pos = 0
                    except IndexError: # MMIA signal actually starts before GLM signal
                        MMIA_start_in_GLM_pos = 0
                        GLM_start_in_MMIA_pos = np.where(current_MMIA[:,0] <= current_GLM[0,0])[0][-1]
                    
                    try: # Try to find where in GLM the MMIA signal ends
                        MMIA_end_in_GLM_pos = np.where(current_GLM[:,0] <= current_MMIA[-1,0])[0][-1]
                        GLM_end_in_MMIA_pos = -1
                    except IndexError: # MMIA signal actually ends after GLM signal
                        MMIA_end_in_GLM_pos = -1
                        GLM_end_in_MMIA_pos = np.where(current_MMIA[:,0] <= current_GLM[-1,0])[0][-1]
                    
                    del it_current_GLM
                    del it_current_MMIA
                    it_current_GLM = current_GLM[MMIA_start_in_GLM_pos:MMIA_end_in_GLM_pos,:]
                    it_current_MMIA = current_MMIA[GLM_start_in_MMIA_pos:GLM_end_in_MMIA_pos,:]
                    
                if abs(delay) >= max_viable_delay and counter > 1:
                    
                    new_GLM_length = len(it_current_GLM) - delay_crop
                    new_current_GLM = np.zeros((new_GLM_length,2))
                    
                    if delay < 0:      # MMIA delays too far from time accuracy
                        new_current_GLM[:,:] = it_current_GLM[0:new_GLM_length, :]
                        del it_current_GLM
                        it_current_GLM = new_current_GLM
                        del new_current_GLM
                        
                    elif delay > 0:    # MMIA anticipates too far from time accuracy
                        new_current_GLM[:,:] = it_current_GLM[delay_crop:len(it_current_GLM), :]
                        del it_current_GLM
                        it_current_GLM = new_current_GLM
                        del new_current_GLM
                        
            del it_current_GLM
            del it_current_MMIA
            
            delays[j] = delay

            GLM_xc = current_GLM     # Normalized!
            del current_GLM
            MMIA_xc = np.zeros((len(current_MMIA),2))

            for k in range(len(current_MMIA)):
                if delay != 0: # There is delay
                    # Adjust Normalized vector
                    MMIA_xc[k,0] = current_MMIA[k,0] + delay*0.002
                    MMIA_xc[k,1] = current_MMIA[k,1]
                    # Adjust original vector
                    MMIA_snippets[j][k,0] = MMIA_snippets[j][k,0] + delay*0.002
                else: # delay==0 so no delay at all
                    MMIA_xc[k,:] = current_MMIA[k,:]
            del current_MMIA

            if tfg == True:
                txt_path = GLM_ordered_outputs + '/' + matches[current_day] + '_' + str(j) + '.txt'
                os.system('cp ' + txt_path + ' ' + ssd_path +'/GLM_TGF_match_txt')

            if show_plots == 1:
                # Plotting cross-correlated and syncronized GLM and MMIA normalized signals
                figure_name = matches[current_day] + '_' + str(j)
                plt.figure(figsize=(12, 4))
                plt.plot(MMIA_xc[:,0], MMIA_xc[:,1], color = 'r', linewidth = 0.5)
                plt.plot(GLM_xc[:,0], GLM_xc[:,1], color = 'black', linewidth = 0.5)
                plt.legend(['MMIA','GLM'])
                plt.title('GLM (black) and MMIA (red) correlated normalized signals for day %d event %d' % (int(matches[current_day]), j))
                plt.xlabel('Time [s]')
                plt.ylabel('Normalized energy')
                plt.grid('on')
                plt.savefig(xcorr_figures + '/' + figure_name + '_xc.pdf')
                #plt.show()
                # Clear the current axes
                plt.cla() 
                # Clear the current figure
                plt.clf() 
                # Closes all the figure windows
                plt.close('all')
                
            
            # Storing non-normalized correlated data
            GLM_xcorr[j] = GLM_snippets[j]
            MMIA_xcorr[j] = MMIA_snippets[j]
            
            # Storing normalized correlated data
            GLM_xcorr_norm[j] = GLM_xc
            del GLM_xc
            MMIA_xcorr_norm[j] = MMIA_xc
            del MMIA_xc

            print('Date %d event %d / %d cross-correlated and aligned!' % (int(matches[current_day]), j, len(GLM_xcorr_norm)))
            
        elif type(GLM_snippets[j]) == np.ndarray and type(MMIA_snippets[j]) != np.ndarray:
            print('Date %d event %d / %d was pre-avoided for lack of MMIA data' % (int(matches[current_day]), j, len(GLM_xcorr_norm)))
            
        elif type(GLM_snippets[j]) != np.ndarray and type(MMIA_snippets[j]) == np.ndarray:
            print('Date %d event %d / %d was pre-avoided for lack of GLM data' % (int(matches[current_day]), j, len(GLM_xcorr_norm)))
            
        elif type(GLM_snippets[j]) != np.ndarray and type(MMIA_snippets[j]) != np.ndarray:
            print('Date %d event %d / %d was pre-avoided for lack of GLM and MMIA data' % (int(matches[current_day]), j, len(GLM_xcorr_norm)))
        elif any(no_overlap_conditions):
            print('GLM and MMIA signals for day %s event %d / %d do not overlap at all' % (matches[current_day], j, len(GLM_xcorr_norm)))
    print(' ')
    print('All events checked!')
    print(' ')
    return [GLM_xcorr, MMIA_xcorr, GLM_xcorr_norm, MMIA_xcorr_norm, delays]

def signal_delay(data1, data2, show_plots, day, snip):
    '''
    This function determines the delay in samples of two signals using a
    cross-correlation method.

    Parameters
    ----------
    data1 : Array
        First signal to be analyzed. In this case, GLM signal.
    data2 : Array
        Second signal to be analyzed. In this case, MMIA signal.
    show_plots : bool
        Boolean variable for showing plots all through the program.
    day : int
        Day of the snip to be cross-correlated in shape YearMonthDay
    snip : int
        Position of snip in day "day"

    Returns
    -------
    real_delay_samples : The number of samples that data1 is shifted with
        respect to data2.
    '''

    xcorr_factors = correlate(data1[:,1], data2[:,1], mode='full', method = 'auto')

    len_x = len(data1)+len(data2)-1
    x = np.empty(len_x)

    for i in range(len_x):
        if (len_x % 2) == 0: # Even number
            x[i] = (i - (len_x/2))
        if (len_x % 2) != 0: # Odd number
            x[i] = (i - (len_x/2 - 0.5))
    
    if show_plots == 1:
        plt.subplot(2, 1, 1)
        plt.plot(data2[:,1],'-r', linewidth = 0.5)
        plt.plot(data1[:,1],'-k', linewidth = 0.5)
        plt.title('Non-correlated GLM (black) and MMIA (red) signals, day %d snippet %d' % (day, snip))
        plt.ylabel('Normalized Energy')
        plt.xlabel('Vector samples')
        plt.legend(['MMIA signal', 'GLM signal'])
        plt.grid('on')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')

        plt.subplot(2, 1, 2)
        plt.plot(x, xcorr_factors, '-b', linewidth = 0.5)
        plt.xlabel('Diff. Samples')
        plt.ylabel('Correlation Factor')
        plt.grid('on')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')

    max_factor_pos = np.where(xcorr_factors == max(xcorr_factors))[0][0]

    if ((len(data1)+len(data2)) % 2 == 0): # len(x) is Odd
        delay_samples = x[max_factor_pos]+(len(data1)-len(data2))/2

    if ((len(data1)+len(data2)) % 2 != 0): # len(x) is Even
        delay_samples = x[max_factor_pos]+(len(data1)-len(data2))/2 + 0.5
    
    # Delay samples accounting actual positioning due to time:
    if data1[0,0] > data2[0,0]: # GLM vector starts later
        pos_MMIA_start_GLM = np.where(data2[:,0] <= data1[0,0])[0][-1]
        real_delay_samples = delay_samples + pos_MMIA_start_GLM
    elif data1[0,0] < data2[0,0]: # GLM vector starts earlier
        pos_GLM_start_MMIA = np.where(data1[:,0] <= data2[0,0])[0][-1]
        real_delay_samples = delay_samples - pos_GLM_start_MMIA
    else:
        real_delay_samples = delay_samples

    return (real_delay_samples)

def normalize(vector):
    '''
    This function takes a signal vector and normalizes it, i.e. makes it fit
    in range 0 to 1 according to its absolute maximum.

    Parameters
    ----------
    vector : array
        Vector to normalize

    Returns
    -------
    None. Same vector but normalized.
    '''
    
    vector_max = max(vector)
    vector_norm = [None]*len(vector)
    for i in range(len(vector)):
        vector_norm[i] = vector[i]/vector_max
    return vector_norm

def signal_data_to_mat(mmia_raw, GLM_xcorr, MMIA_xcorr, delays, current_day, save_path):
    
    # Correct MMIA time on original signals
    for j in range(len(mmia_raw)):
        if type(delays[j]) == int: # Avoid None type delays (no xcorr was possible)
            mmia_raw[j][:,0] = mmia_raw[j][:,0] + delays[j] * 0.002
    
    # Convert variables to final expected output
    delays_t = np.array(delays)
    for j in range(len(delays_t)):
        if type(delays_t[j]) == np.int64: # Avoid None type delays (no xcorr was possible)
            delays_t[j] = delays_t[j] * 0.002
    
    # Save variables into a .mat per event
    for j in range(len(GLM_xcorr)):
        if type(GLM_xcorr[j]) == np.ndarray:
            # Create a dictionary of the variables to save
            vars_to_save = {}
            vars_to_save['corr_mmia_all'] = np.array(mmia_raw[j],dtype=object)
            vars_to_save['GLM_xcorr'] = np.array(GLM_xcorr[j],dtype=object)
            vars_to_save['MMIA_xcorr'] = np.array(MMIA_xcorr[j],dtype=object)
            vars_to_save['delay_t'] = delays_t[j]
            # Save the library into a MATLAB .mat binary file
            sio.savemat(save_path + '/' + current_day + '_' + str(j) + '.mat', vars_to_save)


# Peak characterization and matching

def get_GLM_MMIA_peaks(GLM_xcorr, MMIA_xcorr, GLM_xcorr_norm, MMIA_xcorr_norm, matches, show_plots, current_day, peaks_path, glm_min_peak_num, mmia_min_peak_num):
    '''
    This function gets the cross-correlated vector snippets from GLM and MMIA
    and finds their indexes for every prominent peak in their signals.
    It returns a list of lists of index vectors for every snippet.

    Parameters
    ----------
    GLM_xcorr : list
        List of daily lists of synchronized GLM data.
    MMIA_xcorr : list
        List of daily lists of synchronized MMIA data.
    matches : list
        List of existing GLM and MMIA dates
    show_plots : bool
        Boolean variable for plotting. Not ploting makes the program faster.

    Returns
    -------
    GLM_peaks : list
        List of daily lists with vectors of GLM_xcorr indexes for peaks in the
        signal.
    GLM_peaks : list
        List of daily lists with vectors of MMIA_xcorr indexes for peaks in the
        signal.
    '''
    
    print('Starting peak detection...\n')
    
    GLM_peaks = [None] * len(GLM_xcorr)
    
    MMIA_peaks = [None] * len(MMIA_xcorr)

    for j in range(len(GLM_xcorr)):
            
        if type(GLM_xcorr[j]) == np.ndarray and type(MMIA_xcorr[j]) == np.ndarray:
            
            print('Finding peaks in GLM and MMIA cross-correlated signals for date %s event %d / %d' % (matches[current_day], j, len(GLM_xcorr)))
            
            # Cropping in order to have the same time to compare
            
            # Not overlapping conditions
            GLM_left_cond = GLM_xcorr[j][-1,0]<=MMIA_xcorr[j][0,0]
            GLM_right_cond = GLM_xcorr[j][0,0]>=MMIA_xcorr[j][-1,0]
            
            if GLM_left_cond == True or GLM_right_cond == True:
                print('Correlated events for date %s event %d do not overlap at all' % (matches[current_day], j))
            else:
            
                # Finding the starting position
                GLM_first = 0
                if GLM_xcorr[j][0,0] < MMIA_xcorr[j][0,0]: # GLM starts first
                    start_pos = np.where(GLM_xcorr[j][:,0] <= MMIA_xcorr[j][0,0])[0][-1]
                    GLM_first = 1
                elif GLM_xcorr[j][0,0] > MMIA_xcorr[j][0,0]: # MMIA starts first
                    start_pos = np.where(MMIA_xcorr[j][:,0] <= GLM_xcorr[j][0,0])[0][-1]
                else: # Both start at the sime timestep
                    start_pos = 0
            
                # Finding the end position
                GLM_last = 0
                if GLM_xcorr[j][-1,0] < MMIA_xcorr[j][-1,0]: # GLM ends first
                    end_pos = np.where(MMIA_xcorr[j][:,0] <= GLM_xcorr[j][-1,0])[0][-1]
                elif GLM_xcorr[j][-1,0] > MMIA_xcorr[j][-1,0]: # MMIA ends first
                    end_pos = np.where(GLM_xcorr[j][:,0] <= MMIA_xcorr[j][-1,0])[0][-1]
                    GLM_last = 1
                else: # Both end at the sime timestep
                    end_pos = -1
            
                # Cropping vectors accordingly
                if GLM_first == 1 and GLM_last == 1:
                    GLM_vector = GLM_xcorr[j][start_pos:end_pos,1]
                    GLM_time_vector = GLM_xcorr[j][start_pos:end_pos,0]
                    MMIA_vector = MMIA_xcorr[j][:,1]
                    MMIA_time_vector = MMIA_xcorr[j][:,0]
                
                elif GLM_first == 1 and GLM_last != 1:
                    GLM_vector = GLM_xcorr[j][start_pos:-1,1]
                    GLM_time_vector = GLM_xcorr[j][start_pos:-1,0]
                    MMIA_vector = MMIA_xcorr[j][0:end_pos,1]
                    MMIA_time_vector = MMIA_xcorr[j][0:end_pos,0]
                
                elif GLM_first != 1 and GLM_last == 1:
                    GLM_vector = GLM_xcorr[j][0:end_pos,1]
                    GLM_time_vector = GLM_xcorr[j][0:end_pos,0]
                    MMIA_vector = MMIA_xcorr[j][start_pos:-1,1]
                    MMIA_time_vector = MMIA_xcorr[j][start_pos:-1,0]
                
                elif GLM_first != 1 and GLM_last != 1:
                    GLM_vector = GLM_xcorr[j][:,1]
                    GLM_time_vector = GLM_xcorr[j][:,0]
                    MMIA_vector = MMIA_xcorr[j][start_pos:end_pos,1]
                    MMIA_time_vector = MMIA_xcorr[j][start_pos:end_pos,0]
                
                
                # Cropping GLM_xcorr and MMIA x_corr
                GLM_xcorr_new_snippet = np.zeros((len(GLM_vector),2))
                GLM_xcorr_new_snippet[:,0] = GLM_time_vector
                GLM_xcorr_new_snippet[:,1] = GLM_vector
                GLM_xcorr[j] = GLM_xcorr_new_snippet
                
                MMIA_xcorr_new_snippet = np.zeros((len(MMIA_vector),2))
                MMIA_xcorr_new_snippet[:,0] = MMIA_time_vector
                MMIA_xcorr_new_snippet[:,1] = MMIA_vector
                MMIA_xcorr[j] = MMIA_xcorr_new_snippet
                
                # Cropping GLM_xcorr_norm and MMIA_xcorr_norm
                GLM_norm_vecs = np.zeros((len(GLM_xcorr[j]),2))
                GLM_norm_vecs[:,0] = GLM_xcorr[j][:,0]
                GLM_norm_vecs[:,1] = TFG.normalize(GLM_xcorr[j][:,1])
                GLM_xcorr_norm[j] = GLM_norm_vecs
                
                MMIA_norm_vecs = np.zeros((len(MMIA_xcorr[j]),2))
                MMIA_norm_vecs[:,0] = MMIA_xcorr[j][:,0]
                MMIA_norm_vecs[:,1] = TFG.normalize(MMIA_xcorr[j][:,1])
                MMIA_xcorr_norm[j] = MMIA_norm_vecs
                
            
                # Calculating indexes of peaks in GLM signal
                #GLM_peak_vec, _ = find_peaks(GLM_vector, prominence = 0.3e-14, rel_height = 20)
                # Jesus method
                GLM_noise_level = np.percentile(GLM_vector,50, axis=0)
                GLM_peak_vec, _ = find_peaks(GLM_vector, rel_height = 100, height = GLM_noise_level*1.5, prominence = 0.4)


                # Calculating indexes of peaks in MMIA signal
                #MMIA_noise_level = np.percentile(MMIA_vector,90, axis=0)
                #MMIA_peak_vec, _ = find_peaks(MMIA_vector, rel_height = 100, height = MMIA_noise_level, prominence = 0.4, distance=2)
                # Jesus method
                MMIA_noise_level = np.percentile(MMIA_vector,50, axis=0)
                MMIA_peak_vec, _ = find_peaks(MMIA_vector, rel_height = 100, height = MMIA_noise_level*1.5, prominence = 0.4)
                
                
                # Deleting those triggers with only 2 or less peaks (no meaningful sense)
                    # and assigning values
                if len(GLM_peak_vec) >= glm_min_peak_num and len(MMIA_peak_vec) >= mmia_min_peak_num:
                    GLM_peaks[j] = GLM_peak_vec
                    MMIA_peaks[j] = MMIA_peak_vec
                else:
                    print('GLM or MMIA vector for day %s event %d had less than minimum input peak number, so no reliable results can be extracted' % (matches[current_day], j))
                
                
                if show_plots == 1:
                    
                    # GLM peaks
                    plt.figure()
                    figure_name = matches[current_day] + '_' + str(j)
                    plt.plot(GLM_vector, color = 'black', linewidth=0.5)
                    plt.plot(GLM_peak_vec, GLM_vector[GLM_peak_vec], "*", color='b')
                    plt.title('GLM peaks on day %d, event %d' % (int(matches[current_day]), j))
                    plt.xlabel('Samples')
                    plt.ylabel('Top Cloud Energy [J]')
                    plt.grid('on')
                    #plt.show()
                    plt.savefig(peaks_path + '/' + figure_name + '_glm.pdf')
                    # Clear the current axes
                    plt.cla() 
                    # Clear the current figure
                    plt.clf() 
                    # Closes all the figure windows
                    plt.close('all')

                    # MMIA_peaks
                    plt.figure()
                    figure_name = matches[current_day] + '_' + str(j)
                    plt.plot(MMIA_vector, color = 'r', linewidth=0.5)
                    plt.plot(MMIA_peak_vec, MMIA_vector[MMIA_peak_vec], "*", color='b')
                    plt.title('MMIA peaks on day %d, event %d' % (int(matches[current_day]), j))
                    plt.xlabel('Samples')
                    plt.ylabel('Top Cloud Energy [J]')
                    plt.grid('on')
                    #plt.show()
                    plt.savefig(peaks_path + '/' + figure_name + '_mmia.pdf')
                    # Clear the current axes
                    plt.cla() 
                    # Clear the current figure
                    plt.clf() 
                    # Closes all the figure windows
                    plt.close('all')
                    
                    #GLM_vec = GLM_xcorr_norm[j][:,1]
                    #MMIA_vec = MMIA_xcorr_norm[j][:,1]
                    
                    # GLM and MMIA peaks
                    plt.figure()
                    figure_name = matches[current_day] + '_' + str(j)
                    plt.plot(GLM_vector, color = 'black', linewidth=0.5)
                    plt.plot(GLM_peak_vec, GLM_vector[GLM_peak_vec], "*", color='gold')
                    plt.plot(MMIA_vector, color = 'r', linewidth=0.5)
                    plt.plot(MMIA_peak_vec, MMIA_vector[MMIA_peak_vec], "*", color='b')
                    plt.title('GLM (black-yellow) and MMIA (red-blue) peaks on day %d, event %d' % (int(matches[current_day]), j))
                    plt.xlabel('Samples')
                    plt.ylabel('Top Cloud Energy [J]')
                    plt.grid('on')
                    #plt.show()
                    plt.savefig(peaks_path + '/' + figure_name + '_both.pdf')
                    # Clear the current axes
                    plt.cla() 
                    # Clear the current figure
                    plt.clf() 
                    # Closes all the figure windows
                    plt.close('all')

        else:
            print('Date %s event %d was not cross correlated' % (matches[current_day], j))
    
    print(' ')            
    print('Done!')
    print(' ')
    
    return [GLM_peaks, MMIA_peaks]

def get_peak_matches(GLM_xcorr, MMIA_xcorr, GLM_peaks, MMIA_peaks, show_plots, matches, current_day, match_figs_path):
    
    print('Starting peak matching...\n')
    
    common_peaks = [None] * len(GLM_xcorr)
    # Time and signal values for matching peaks for outputting to .mat
    common_peaks_values = [None] * len(GLM_xcorr)
    
    for j in range(len(GLM_xcorr)): # For every event in current day
            
        # If data was successfully correlated and with 3 or more peaks (condition imposed before)
        if type(GLM_peaks[j]) == np.ndarray and type(MMIA_peaks[j]) == np.ndarray:
            
            print('Matching peaks for date %s event %d / %d...' % (matches[current_day], j, len(GLM_xcorr)))
            
            # Assuring all the arrays are of int type
            
            GLM_peaks[j] = GLM_peaks[j].astype(int)
            MMIA_peaks[j] = MMIA_peaks[j].astype(int)
            
            # Checking GLM peaks in MMIA
    
            GLM_peaks_in_MMIA = np.zeros((len(GLM_peaks[j]),2), dtype=int)
            # Structure:
                # 1st column: Boolean for matching GLM peak in MMIA
                # 2nd column: Position of the matching peak in MMIA_peaks[j]
            
            glm_peak_values = GLM_xcorr[j][GLM_peaks[j],1]  # Signal value of the peak
            GLM_peaks_copy = np.zeros((len(GLM_peaks[j])))  # Copying peaks position vectors
            GLM_peaks_copy[:] = GLM_peaks[j][:]
            GLM_peaks_copy = GLM_peaks_copy.astype(int)
            MMIA_peaks_copy = np.zeros((len(MMIA_peaks[j])))
            MMIA_peaks_copy[:] = MMIA_peaks[j][:]
            MMIA_peaks_copy = MMIA_peaks_copy.astype(int)
            
            for k in range(len(GLM_peaks[j])):
                
                # Position of the maximum peak position computation
                
                glm_peak_values_copy = GLM_xcorr[j][GLM_peaks_copy,1]
                
                max_pos = np.where(glm_peak_values == max(glm_peak_values_copy))[0][0]
                
                # Max peak is next max peak for next iteration
                GLM_peaks_copy = np.delete(GLM_peaks_copy, np.where(GLM_peaks_copy == GLM_peaks[j][max_pos]))
                
                # Starting at the maximum peak
                current_peak_pos = GLM_peaks[j][max_pos]
                
                window = np.linspace(current_peak_pos-1,current_peak_pos+1, 3, dtype=int) # Searching in a window of 3 samples
                # Precision of +- 0.002s --> 0.004s window
                # 0.004 [s] / 0.002 [s/sample] = 2 [samples window]
                # 2 samples window --> ** +-1 samples **
                
                
                # Check if any value of the window exists in MMIA_peaks[j]
                bool_peaks_in_mmia = np.isin(MMIA_peaks_copy, window)
                
                if bool_peaks_in_mmia.any(): # If there's at least one 'True'
                
                    # See where in MMIA_peaks_copy is a peak in the window
                    coincidences = np.where(bool_peaks_in_mmia == True)[0]
                
                    # See the values of MMIA_xcorr in the peaks inside the window
                    in_window_peak_values = MMIA_xcorr[j][MMIA_peaks_copy[coincidences],1]
                
                    # Get the highest MMIA peak position in coincidences in that window
                    max_mmia_pos_in_coincidences = np.where(in_window_peak_values == max(in_window_peak_values))[0][0]
                
                    # Selection of the matching peak
                    desired_peak_pos = MMIA_peaks_copy[coincidences[max_mmia_pos_in_coincidences]]
                    
                    # GLM_peaks_in_MMIA compleition
                    GLM_peaks_in_MMIA[max_pos,0] = 1
                    GLM_peaks_in_MMIA[max_pos,1] = np.where(MMIA_peaks[j] == desired_peak_pos)[0][0]
                    
                    # Delete desired_peak_pos from MMIA_peaks_copy to not repeat it
                    MMIA_peaks_copy = np.delete(MMIA_peaks_copy, np.where(MMIA_peaks_copy == MMIA_peaks[j][GLM_peaks_in_MMIA[max_pos,1]]))


            matching_GLM_peaks_pos_pos = [] # Positions in GLM_peaks[j] with matching peaks
            matching_MMIA_peaks_pos_pos = [] # Positions in MMIA_peaks[j] with matching peaks
        
            for k in range(len(GLM_peaks_in_MMIA)):
                
                if GLM_peaks_in_MMIA[k,0] == 1:
                    
                    matching_GLM_peaks_pos_pos.append(k)
                    matching_MMIA_peaks_pos_pos.append(GLM_peaks_in_MMIA[k,1])
                
            common_peaks[j] = [matching_GLM_peaks_pos_pos, matching_MMIA_peaks_pos_pos]
            
            common_peaks_values[j] = [GLM_xcorr[j][GLM_peaks[j][matching_GLM_peaks_pos_pos],0],     # GLM time for matching peaks
                                      GLM_xcorr[j][GLM_peaks[j][matching_GLM_peaks_pos_pos],1],     # GLM signal value for matching peaks
                                      MMIA_xcorr[j][MMIA_peaks[j][matching_MMIA_peaks_pos_pos],0],  # MMIA time for matching peaks
                                      MMIA_xcorr[j][MMIA_peaks[j][matching_MMIA_peaks_pos_pos],1]]  # MMIA signal value for matching peaks


            if show_plots == True:
                
                figure_name = matches[current_day] + '_' + str(j)
                # GLM signal with common peaks
                plt.figure()
                plt.plot(GLM_xcorr[j][:,1], color = 'black', linewidth=0.5)
                plt.plot(GLM_peaks[j], GLM_xcorr[j][GLM_peaks[j],1], "*", color='b')
                plt.plot(GLM_peaks[j][matching_GLM_peaks_pos_pos], GLM_xcorr[j][GLM_peaks[j][matching_GLM_peaks_pos_pos],1], "*", color='gold')
                plt.title('GLM peaks on day %s, event %d' % (matches[current_day], j))
                plt.xlabel('Samples')
                plt.ylabel('Top Cloud Energy [J]')
                plt.grid('on')
                #plt.show()
                plt.savefig(match_figs_path + '/' + figure_name + '_GLM_match.pdf')
                # Clear the current axes
                plt.cla() 
                # Clear the current figure
                plt.clf() 
                # Closes all the figure windows
                plt.close('all')

                plt.figure()
                plt.plot(MMIA_xcorr[j][:,1], color = 'r', linewidth=0.5)
                plt.plot(MMIA_peaks[j], MMIA_xcorr[j][MMIA_peaks[j],1], "*", color='b')
                plt.plot(MMIA_peaks[j][matching_MMIA_peaks_pos_pos], MMIA_xcorr[j][MMIA_peaks[j][matching_MMIA_peaks_pos_pos],1], "*", color='gold')
                plt.title('MMIA peaks on day %s, event %d' % (matches[current_day], j))
                plt.xlabel('Samples')
                plt.ylabel('Top Cloud Energy [J]')
                plt.grid('on')
                #plt.show()
                plt.savefig(match_figs_path + '/' + figure_name + '_MMIA_match.pdf')
                # Clear the current axes
                plt.cla() 
                # Clear the current figure
                plt.clf() 
                # Closes all the figure windows
                plt.close('all')
                
                plt.figure(figsize=(12, 4))
                plt.plot(MMIA_xcorr[j][:,1], color = 'r', linewidth=0.5)
                plt.plot(MMIA_peaks[j], MMIA_xcorr[j][MMIA_peaks[j],1], "*", color='lime')
                plt.plot(MMIA_peaks[j][matching_MMIA_peaks_pos_pos], MMIA_xcorr[j][MMIA_peaks[j][matching_MMIA_peaks_pos_pos],1], "*", color='gold')
                plt.plot(GLM_xcorr[j][:,1], color = 'black', linewidth=0.5)
                plt.plot(GLM_peaks[j], GLM_xcorr[j][GLM_peaks[j],1], "*", color='b')
                plt.plot(GLM_peaks[j][matching_GLM_peaks_pos_pos], GLM_xcorr[j][GLM_peaks[j][matching_GLM_peaks_pos_pos],1], "*", color='darkorange')
                plt.title('GLM (black) and MMIA (red) peaks and common peaks on day %s, event %d' % (matches[current_day], j))
                plt.xlabel('Samples')
                plt.ylabel('Top Cloud Energy [J]')
                plt.legend(['MMIA corr norm signal', 'MMIA peaks', 'MMIA matching peaks', 'GLM corr norm signal', 'GLM peaks', 'GLM matching peaks'])
                plt.grid('on')
                #plt.show()
                plt.savefig(match_figs_path + '/' + figure_name + '_both_match.pdf')
                # Clear the current axes
                plt.cla() 
                # Clear the current figure
                plt.clf() 
                # Closes all the figure windows
                plt.close('all')
                
        else:
            print('Date %s event %d was not suitable for study' % (matches[current_day], j))

    print(' ')
    print('Done!')
    print(' ')
    
    return [common_peaks, common_peaks_values]


# General statistics outputting

def get_ministats(GLM_xcorr, MMIA_xcorr):
    
    # Generation of lists
    GLM_avg = [None] * len(GLM_xcorr)
    MMIA_avg = [None] * len(MMIA_xcorr)
    GLM_std = [None] * len(GLM_xcorr)
    MMIA_std = [None] * len(MMIA_xcorr)
    
    for j in range(len(GLM_xcorr)):
        if type(GLM_xcorr[j]) == np.ndarray:
            GLM_avg[j] = np.average(GLM_xcorr[j][:,1])
            MMIA_avg[j] = np.average(GLM_xcorr[j][:,1])
            GLM_std[j] = np.std(GLM_xcorr[j][:,1])
            MMIA_std[j] = np.std(GLM_xcorr[j][:,1])
    
    return [GLM_avg, MMIA_avg, GLM_std, MMIA_std]

def study_delays(statistics_bin, show_plots, statistics_figures_path, matches, ssd_path, outputting_to_mat):
    '''
    This function computes the average delay of GLM with respect to MMIA,
    both attending absolute values and real values.
    '''
    print('Computing delay statistics...')
    
    # Generating necessary global matrices
    
    delays = TFG.merge_days(matches, statistics_bin, 'delays')
    glm_avg = TFG.merge_days(matches, statistics_bin, 'glm_avg')
    mmia_avg = TFG.merge_days(matches, statistics_bin, 'mmia_avg')
    glm_std = TFG.merge_days(matches, statistics_bin, 'glm_std')
    mmia_std = TFG.merge_days(matches, statistics_bin, 'mmia_std')

    # Initializing variables
    valid_trigger_counter = 0
    delay_list = []
    MMIA_delays_list = []
    GLM_delays_list = []
    no_delays_list = []
    trigger_counter = 0
    
    for i in range(len(delays)):
        for j in range(len(delays[i])):
            trigger_counter = trigger_counter + 1
            if type(delays[i][j]) == int:
                valid_trigger_counter = valid_trigger_counter + 1
                delay_list.append(delays[i][j])
                if delays[i][j] < 0:
                    MMIA_delays_list.append(delays[i][j])
                elif delays[i][j] > 0:
                    GLM_delays_list.append(delays[i][j])
                else:
                    no_delays_list.append(delays[i][j])

    all_delays = np.array(delay_list)
    if len(all_delays) != 0:
        avg_all = np.average(all_delays)
        std_all = np.std(all_delays)
    
    MMIA_delays = np.array(MMIA_delays_list)
    if len(MMIA_delays) != 0:
        avg_negative = np.average(MMIA_delays)
        std_negative = np.std(MMIA_delays)
        
    # (GLM delays is really MMIA sends early)
    
    GLM_delays = np.array(GLM_delays_list)
    if len(GLM_delays) != 0:
        avg_positive = np.average(GLM_delays)
        std_positive = np.std(GLM_delays)
        
    no_delays = np.array(no_delays_list)
    
    valid_triggers = len(MMIA_delays) + len(GLM_delays) + len(no_delays)
    
    # Checking relationship between delay and signal average energy
    
    GLM_delay_signal = np.zeros((len(GLM_delays),3))
    MMIA_delay_signal = np.zeros((len(MMIA_delays),3))
    GLM_counter = 0
    MMIA_counter = 0
    
    # Creating a delays vector just for histogram plotting
    delays_s_vec = []
    
    for i in range(len(delays)):
        for j in range(len(delays[i])):
            if type(delays[i][j]) == int:
                
                delays_s_vec.append(delays[i][j]*0.002)
                
                if delays[i][j] > 0: # MMIA signal anticipates
                    GLM_delay_signal[GLM_counter,0] = delays[i][j] # Delay for that snippet
                    GLM_delay_signal[GLM_counter,1] = glm_avg[i][j] # Average GLM energy for that snippet
                    GLM_delay_signal[GLM_counter,2] = glm_std[i][j] # Std deviation of GLM energy for that snippet
                    GLM_counter = GLM_counter + 1
                elif delays[i][j] < 0: # MMIA signal delays
                    MMIA_delay_signal[MMIA_counter,0] = delays[i][j] # Delay for that snippet
                    MMIA_delay_signal[MMIA_counter,1] = mmia_avg[i][j] # Average MMIA energy for that snippet
                    MMIA_delay_signal[MMIA_counter,2] = mmia_std[i][j] # Std deviation of MMIA energy for that snippet
                    MMIA_counter = MMIA_counter + 1
    
    print('Done!')

    # Adding delays_s_vec to dictionary to outputting .mat
    outputting_to_mat['all_delays_t'] = delays_s_vec

    print('Writing results into ' + ssd_path + '/RESULTS.txt...')
    # Saving results into a .txt
    f =  open(ssd_path + '/RESULTS.txt', 'w')
    f.write('******** RESULTS ********')
    f.write('\n')
    f.write('\n')
    f.write('\n')
    f.write('* EVENT INFO:')
    f.write('\n')
    f.write('    --> Total number of events processed: %d\n' % trigger_counter)
    f.write('    --> Number of valid events: %d\n' % valid_trigger_counter)
    f.write('    --> Percentage of valid over total events: %s%%\n' % str(format(valid_trigger_counter/trigger_counter * 100, '.3f')))
    f.write('\n')
    f.write('\n')
    f.write('* DELAY INFO:')
    f.write('\n')
    f.write('    --> Average delay in samples (accounting for positive and negative values): %s +- %s\n' % (str(format(avg_all, '.3f')), str(format(std_all, '.3f'))))
    f.write('    --> Average delay in seconds (accounting for positive and negative values): %s +- %s\n' % (str(format(avg_all*0.002, '.3f')), str(format(std_all*0.002, '.3f'))))
    f.write('    --> Number of MMIA delays: %d (%s%% over valid events)\n' % (len(MMIA_delays), str(format(len(MMIA_delays)/valid_triggers*100, '.3f'))))
    if 'avg_negative' in locals():
        f.write('       --> Average MMIA delay in samples: %s +- %s\n' % (str(format(avg_negative, '.3f')), str(format(std_negative, '.3f'))))
        f.write('       --> Average MMIA delay in seconds: %s +- %s\n' % (str(format(avg_negative*0.002, '.3f')), str(format(std_negative*0.002, '.3f'))))
    else:
        f.write('       --> Average MMIA delay in samples: -\n')
        f.write('       --> Average MMIA delay in seconds: -\n')
        
    f.write('    --> Number of MMIA anticipations: %d (%s%% over valid events)\n' % (len(GLM_delays), str(format(len(GLM_delays)/valid_triggers*100, '.3f'))))
    if 'avg_positive' in locals():
        f.write('       --> Average MMIA anticipation in samples: %s +- %s\n' % (str(format(avg_positive, '.3f')), str(format(std_positive, '.3f'))))
        f.write('       --> Average MMIA anticipation in seconds: %s +- %s\n' % (str(format(avg_positive*0.002, '.3f')), str(format(std_positive*0.002, '.3f'))))
    else:
        f.write('       --> Average MMIA anticipation in samples: -\n')
        f.write('       --> Average MMIA anticipation in seconds: -\n')
    f.write('    --> Number of no delays: %d (%s%%)\n' % (len(no_delays), str(format(len(no_delays)/valid_triggers*100, '.3f'))))
    f.close()
    print('Done!\n')
    
    if show_plots == True:
        # Good snippets vs total snippets
        labels = ['Valid snippets', 'Non valid snippets']
        colors = ['mediumaquamarine','lightcoral']
        sizes = [valid_triggers, trigger_counter-valid_triggers]
        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, labels=labels, autopct='%1.2f%%', shadow=False, startangle=90, colors=colors)
        ax1.axis('equal')
        #plt.show()
        plt.savefig(statistics_figures_path + '/good_vs_total_triggers.pdf')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')
        
        # MMIA delays vs GLM delays vs No delays
        labels = ['MMIA delayed', 'MMIA anticipated', 'No delay']
        colors = ['darksalmon','slategrey','black']
        sizes = [len(MMIA_delays), len(GLM_delays), len(no_delays)]
        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, labels=labels, autopct='%1.2f%%', shadow=False, startangle=90, colors=colors)
        ax1.axis('equal')
        #plt.show()
        plt.savefig(statistics_figures_path + '/MMIA_vs_GLM_vs_no_delays.pdf')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')
        
        # Averaging in delays
        if ('avg_positive' in locals()) and ('avg_negative' in locals()):
            fig = plt.figure(dpi=800)
            ax = fig.add_axes([0,0,1,1])
            langs = ['Avg. Delay', 'Avg. MMIA Delay', 'Avg. MMIA Anticipation']
            data = [avg_all*0.002, avg_negative*0.002, avg_positive*0.002]
            plt.grid('on')
            colors = ['tab:blue', 'darksalmon', 'slategrey']
            ax.set_ylabel('Delay [s]')
            ax.bar(langs,data, color=colors, yerr=[std_all*0.002,std_negative*0.002,std_positive*0.002])
            #plt.show()
            plt.savefig(statistics_figures_path + '/avg_delays.pdf')
            # Clear the current axes
            plt.cla() 
            # Clear the current figure
            plt.clf() 
            # Closes all the figure windows
            plt.close('all')
        
        # GLM delay vs GLM avg energy
        plt.figure()
        plt.scatter(GLM_delay_signal[:,1]*0.002, GLM_delay_signal[:,0]*0.002, color='black', marker='x')
        plt.xlabel('GLM event average energy [J]')
        plt.ylabel('Trigger delay [s]')
        plt.grid('on')
        #plt.show()
        plt.savefig(statistics_figures_path + '/glm_delay_vs_glm_avg_energy.pdf')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')
        
        # MMIA delay vs MMIA avg energy
        plt.figure()
        plt.scatter(MMIA_delay_signal[:,1]*0.002, MMIA_delay_signal[:,0]*0.002, color='r', marker='x')
        plt.xlabel('MMIA event average energy [J]')
        plt.ylabel('Trigger delay [s]')
        plt.grid('on')
        #plt.show()
        plt.savefig(statistics_figures_path + '/mmia_delay_vs_mmia_avg_energy.pdf')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')
        
        # GLM delay vs GLM std deviation in energy
        plt.figure()
        plt.scatter(GLM_delay_signal[:,2]*0.002, GLM_delay_signal[:,0]*0.002, color='black', marker='x')
        plt.xlabel("GLM event's energy stantard deviation [J]")
        plt.ylabel('Trigger delay [s]')
        plt.grid('on')
        #plt.show()
        plt.savefig(statistics_figures_path + '/glm_delay_vs_glm_std_energy.pdf')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')
        
        # MMIA delay vs MMIA std deviation in energy
        plt.figure()
        plt.scatter(MMIA_delay_signal[:,2]*0.002, MMIA_delay_signal[:,0]*0.002, color='r', marker='x')
        plt.xlabel("MMIA event's energy stantard deviation [J]")
        plt.ylabel('Trigger delay [s]')
        plt.grid('on')
        #plt.show()
        plt.savefig(statistics_figures_path + '/mmia_delay_vs_mmia_std_energy.pdf')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')
        
        # Delays histogram for all delays
        
        plt.figure()
        plt.hist(delays_s_vec, bins = 50, rwidth=0.85)
        plt.title('All Delay Histogram')
        plt.xlabel('Delay [s]')
        plt.ylabel('Frequency')
        plt.grid('on')
        #plt.show()
        plt.savefig(statistics_figures_path + '/all_delay_histogram.pdf')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')
        
        # Delays histogram for all delays in range
        plt.figure()
        plt.hist(delays_s_vec, bins = 50, range = [-0.05, 0.05], rwidth=0.85)
        plt.title('Delay Histogram (up to 0.05s of absolute delay)')
        plt.xlabel('Delay [s]')
        plt.ylabel('Frequency')
        plt.grid('on')
        #plt.show()
        plt.savefig(statistics_figures_path + '/inrange_delay_histogram.pdf')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')
        
        delays_s_vec = np.array(delays_s_vec)
        delays_s_vec = abs(delays_s_vec)
        count, bins_count = np.histogram(delays_s_vec, bins=100)
        pdf = count / sum(count)    # Range 0-1
        cdf = np.cumsum(pdf)
        
        # Cumulative Curve for all delays
        plt.figure()
        plt.plot(bins_count[1:], cdf, label="CDF")
        plt.title('Cumulative Curve')
        plt.xlabel('Delay [s]')
        plt.ylabel('Cumulative Frequency')
        plt.grid('on')
        #plt.show()
        plt.savefig(statistics_figures_path + '/cumulative_all.pdf')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')
        
        count, bins_count = np.histogram(delays_s_vec, bins=100, range=[0,1])
        pdf = count / sum(count)    # Relative range 0-1
        cdf = np.cumsum(pdf)
        
        # Cumulative Curve up to 1s for all delays
        plt.figure()
        plt.plot(bins_count[1:], cdf, label="CDF")
        plt.title('Cumulative Curve (up to 1s of absolute delay)')
        plt.xlabel('Delay [s]')
        plt.ylabel('Cumulative Frequency')
        plt.grid('on')
        #plt.show()
        plt.savefig(statistics_figures_path + '/inrange_cumulative.pdf')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')

def more_statistics(peaks_bin, matches, ssd_path, outputting_to_mat):
    
    print('Computing peak statistics...')
    
    # Generating necessary global matrices
    GLM_peaks = TFG.merge_days(matches, peaks_bin, 'glm')
    MMIA_peaks = TFG.merge_days(matches, peaks_bin, 'mmia')
    matching_peaks = TFG.merge_days(matches, peaks_bin, 'matching')
    matching_peaks_vals = TFG.merge_days(matches, peaks_bin, 'matching_vals')
    
    # Setting variables to study
    counter = 0
    sum_rel_GLM = 0 # en mmia/todos
    sum_rel_MMIA = 0 # en glm/todos
    sum_GLM_peaks = 0 # total peaks
    sum_matching_peaks = 0
    sum_MMIA_peaks = 0
    matching_peak_time_signal_matrix = []
    
    for i in range(len(GLM_peaks)):         # For every day in 'matches'
        for j in range(len(GLM_peaks[i])):  # For every event inside that day
            if type(matching_peaks[i][j]) == list: # If there are matching peaks
                
                counter = counter + 1       # Number of events with matching peaks
                
                # Matched vs Total Peaks for every instrument
                sum_rel_GLM = sum_rel_GLM + len(matching_peaks[i][j][0])/len(GLM_peaks[i][j])
                sum_rel_MMIA = sum_rel_MMIA + len(matching_peaks[i][j][1])/len(MMIA_peaks[i][j])
                
                # Total peaks of every instrument
                sum_GLM_peaks = sum_GLM_peaks + len(GLM_peaks[i][j])
                sum_MMIA_peaks = sum_MMIA_peaks + len(MMIA_peaks[i][j])
                
                # Total matched peaks
                sum_matching_peaks = sum_matching_peaks + len(matching_peaks[i][j][0])

                # Add every matching peak to the outputting matrix
                for k in range(len(matching_peaks_vals[i][j][0])): # For every peak
                    matching_peak_time_signal_matrix.append([matching_peaks_vals[i][j][0][k],   # GLM matching peak time
                                                             matching_peaks_vals[i][j][1][k],   # GLM matching peak signal value
                                                             matching_peaks_vals[i][j][2][k],   # MMIA matching peak time
                                                             matching_peaks_vals[i][j][3][k]])  # MMIA matching peak signal value
    
    avg_GLM_rel = sum_rel_GLM/counter   # Average matched vs total GLM peaks per event
    avg_MMIA_rel = sum_rel_MMIA/counter # Average matched vs total MMIA peaks per event
    
    avg_GLM_peaks = sum_GLM_peaks/counter   # Average found GLM peaks per event
    avg_MMIA_peaks = sum_MMIA_peaks/counter # Average found MMIA peaks per event
    
    avg_matching = sum_matching_peaks/counter # Average matched peaks per event
    
    print('Done!')
    print('Adding peak statistics to ' + ssd_path + '/RESULTS.txt...')
    
    # Writing on the outputting file
    f =  open(ssd_path + '/RESULTS.txt', 'a')
    f.write('\n')
    f.write('\n')
    f.write('* PEAKS INFO:')
    f.write('\n')
    f.write('    --> Average GLM peaks found per event: %s\n' % str(format(avg_GLM_peaks, '.3f')))
    f.write('    --> Average MMIA peaks found per event: %s\n' % str(format(avg_MMIA_peaks, '.3f')))
    f.write('    --> Average matched peaks found per event: %s\n' % str(format(avg_matching, '.3f')))
    f.write('    --> Average percentage of GLM matched peaks over GLM peaks found per event: %s%%\n' % str(format(avg_GLM_rel*100, '.3f')))
    #f.write('    --> Average GLM matched peaks over GLM peaks found per event: %s\n' % str(format(avg_GLM_rel/avg_GLM_peaks, '.3f')))
    f.write('    --> Average percentage of MMIA matched peaks over MMIA peaks found per event: %s%%\n' % str(format(avg_MMIA_rel*100, '.3f')))
    #f.write('    --> Average MMIA matched peaks over MMIA peaks found per event: %s\n' % str(format(avg_MMIA_rel/avg_MMIA_peaks, '.3f')))
    f.close()
    print(' ')
    print('Done! Your final results can be accessed at ' + ssd_path + '/RESULTS.txt\n')
    print('Computing outputting variables for external statistics...')

    # Computing and adding peak values and time distribution to outputting .mat
    matching_peaks_per_time = np.zeros((12,1))
    for i in range(len(matching_peak_time_signal_matrix)):  # For every matching peak
        hour = matching_peak_time_signal_matrix[i][0]/3600  # Taking GLM time as reference
        if hour >= 0 and hour < 2:
            matching_peaks_per_time[0] = matching_peaks_per_time[0] + 1
        elif hour >= 2 and hour < 4:
            matching_peaks_per_time[1] = matching_peaks_per_time[1] + 1
        elif hour >= 4 and hour < 6:
            matching_peaks_per_time[2] = matching_peaks_per_time[2] + 1
        elif hour >= 6 and hour < 8:
            matching_peaks_per_time[3] = matching_peaks_per_time[3] + 1
        elif hour >= 8 and hour < 10:
            matching_peaks_per_time[4] = matching_peaks_per_time[4] + 1
        elif hour >= 10 and hour < 12:
            matching_peaks_per_time[5] = matching_peaks_per_time[5] + 1
        elif hour >= 12 and hour < 14:
            matching_peaks_per_time[6] = matching_peaks_per_time[6] + 1
        elif hour >= 14 and hour < 16:
            matching_peaks_per_time[7] = matching_peaks_per_time[7] + 1
        elif hour >= 16 and hour < 18:
            matching_peaks_per_time[8] = matching_peaks_per_time[8] + 1
        elif hour >= 18 and hour < 20:
            matching_peaks_per_time[9] = matching_peaks_per_time[9] + 1
        elif hour >= 20 and hour < 22:
            matching_peaks_per_time[10] = matching_peaks_per_time[10] + 1
        else:
            matching_peaks_per_time[11] = matching_peaks_per_time[11] + 1
    
    outputting_to_mat['glm_mmia_matching_time_signal'] = np.array(matching_peak_time_signal_matrix,dtype=object)
    outputting_to_mat['matching_time_distribution'] = np.array(matching_peaks_per_time,dtype=object)

    # Computing and adding GLM matching peaks over GLM peaks and same for MMIA per event
    peak_relation = np.zeros((counter, 2))
    pos = 0
    for i in range(len(GLM_peaks)):         # For every day in 'matches'
        for j in range(len(GLM_peaks[i])):  # For every event inside that day
            if type(matching_peaks[i][j]) == list: # If there are matching peaks
                # Matching peaks VS GLM detected peaks
                peak_relation[pos,0] = len(matching_peaks[i][j][0])/len(GLM_peaks[i][j])
                # Matching peaks VS MMIA detected peaks
                peak_relation[pos,1] = len(matching_peaks[i][j][1])/len(MMIA_peaks[i][j])
                pos = pos+1

    outputting_to_mat['peak_relations'] = np.array(peak_relation,dtype=object)
