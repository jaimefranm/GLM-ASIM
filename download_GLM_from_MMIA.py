#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script gets a directory containing MMIA snippets and automatically
downloads corresponding GLM .nc files from AWS servers.



Notes:
    1. No need to daily order GLM files in main_TFG.py
    2. Matches will no longer be useful (matches will just be the days list)
        after all GLM files have been downloaded
    3. Get linet timing will be done by MMIA data directly (no need for LINET here)
    4. Important to extract ISS location and corrected time from MatLab
    5. 


@author: Jaime F. Moran Dominguez (jaime.francisco.moran@upc.edu)
"""
import os
import scipy.io as sio

# Path where MMIA snippet info is located
MMIA_info_path = '/home/lrg/Desktop/MMIA_test/snippet_info'

# Path where GLM snippet datafiles will be located
GLM_ordered_dir = '/home/lrg/Desktop/MMIA_test/GLM_snippets'


def get_GLM_datafiles_from_AWS(MMIA_info_path, GLM_ordered_dir):
    
    print('Starting GLM download from AWS for every MMIA snippet...')
    print(' ')
    # Getting all snippets in the directory
    with os.scandir(MMIA_info_path) as snippets:
        snippets = [file.name for file in snippets if file.is_dir() and file.name.endswith('.mat')]
    
    for i in range(len(snippets)):
        current_path = MMIA_info_path + '/' + snippets[i]
        mat = sio.loadmat(current_path)
        snippet_info = mat.get('MMIA_all')
        # Snippet info structure is:
            # 1st value : Minimum latitude of ISS
            # 2nd value : Maximum latitude of ISS
            # 3rd value : Minimum longitude of ISS
            # 4th value : Maximum longitude of ISS
            # 5th value : Minimum data corrected time
            # 6th value : Maximum data corrected time
        
        # Conversion from MMIA data to DD-MM-YYYY-HH:mm
        
        
        
        os.system("python goesaws.py -i 'glm' --start '" + 09-01-2019-16:00 + "' --end '" + 09-01-2019-16:30 + "' --dl -o " + GLM_ordered_dir)