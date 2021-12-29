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
import library as TFG

# PC UPC
#path_to_mmia_files = '/media/lrg/012EE6107EB7CB6B/mmia_20'
#path_to_mmia_files = '/home/lrg/Desktop/test_cdf'
#ssd_path = '/media/lrg/012EE6107EB7CB6B'

#path_to_mmia_files = '/media/lrg/mmia_20'
#ssd_path = '/media/lrg'

# LOCAL mac

# For testing
#path_to_mmia_files = '/Users/jaimemorandominguez/Desktop/Final/MMIA_archivos/cdf'
#path_to_mmia_files = '/Users/jaimemorandominguez/Desktop/test_cdf'
#ssd_path = '/Users/jaimemorandominguez/Desktop/test_descarga_GLM'

# For doing something
path_to_mmia_files = '/Volumes/Jaime_F_HD/mmia_2020/mmia_20'
ssd_path = '/Volumes/Jaime_F_HD/mmia_2020'

trigger_length = 2 # [s]
pre_extracted_MMIA = True
pre_downloaded_GLM = False





[matches, trigger_filenames] = TFG.get_MMIA_triggers(path_to_mmia_files, trigger_length)

if pre_extracted_MMIA == False:
    
    TFG.create_MMIA_trigger_directories(matches, trigger_filenames, path_to_mmia_files, ssd_path)

    [mmia_raw, trigger_limits] = TFG.extract_trigger_info(ssd_path, trigger_filenames, matches)

else:
    [mmia_raw, trigger_limits] = TFG.upload_MMIA_mats(ssd_path, trigger_filenames, matches)
    
if pre_downloaded_GLM == False:
    TFG.download_GLM(ssd_path, trigger_filenames, mmia_raw)
