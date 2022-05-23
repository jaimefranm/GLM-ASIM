# GLM-ASIM

## Contents
1. [Overview](#overview)
2. [Setup](#setup)
   1. [Python](#python)
   2. [MATLAB](#matlab)
   3. [General](#general)
3. [Fast Revision](#revision)
	1. [Global structure](#structure)
	2. [Input Data](#input)
	3. [Output Data](#output)
	4. [First Execution Checklist](#checklist)
	5. [Flow Diagram for `main.py`](#flowchart)

## Overview <a name="overview"></a>
This program compares GLM to MMIA signal based on input MMIA triggers as \textit{.cdf} files. As main points, this program automatically gets:

1. Given MMIA triggers ordered by event in different directories
2. Downloaded corresponding GLM files from Google Cloud servers as \textit{.nc} files
3. Extracted MMIA and GLM data for those files as \textit{.mat} and \textit{.txt} files, respectively
4. Top Cloud Energy - converted GLM and MMIA signals
5. Cross-correlated GLM-MMIA signals (and thus, MMIA delay with respect to GLM)
6. Peaks found in both signals
7. Matching peaks between both signals (GLM peaks present in MMIA signal)
8. Some statistics about MMIA delays and peak-matching

## Setup and Dependencies <a name="setup"></a>
This script is mostly written in Python although some parts require MATLAB to be installed (with a valid license). Next steps are presented in groups as for Python-specific setup, MATLAB-specific setup and general setup.
>It is important to note that this script was designed to be run on a macOS or Linux-like machine

### Python <a name="python"></a>
As this script is mainly written in Python, a Python interpreter needs to be installed ([Python 3.9.8](https://www.python.org/downloads/release/python-398/) was used in the development). Some packages are needed and their installation is shown below. Apart from those packages, no more Python-specific setup needs to be done in order to run the code properly.

```python
# Numpy
pip install numpy

# Matplotlib (for plotting)
pip install matplotlib

# Scipy
pip install scipy

# Pandas
pip install pandas

# Pickle (for binaries storage)
pip install pickleshare

# netCDF4 for reading GLM .nc files
pip install netCDF4

# Google storage packages
pip install gsutil
pip install google
pip install --upgrade google-api-python-client
pip install google.cloud.storage
```

### MATLAB <a name="matlab"></a>
As some part of the process requires a MATLAB run (specifically MMIA data outputting from *.cdf* files into *.mat* files for further processing), this software needs to be installed and properly activated. For that extraction to be done, a MATLAB patch for *.cdf* files must be downloaded, found [here](https://cdf.gsfc.nasa.gov/html/matlab_cdf_patch.html) for macOS, Linux or Windows systems. Apart from that patch, no more setup is needed for MATLAB issues.

### General <a name="general"></a>
In order for this script to run correctly, some previous minimal setup has to be considered.

First of all, at least 3 different directories should exist previous to the execution of the script. Those directories can be wherever in the computer's file system, separate or nested into one another, but must exist as three different folders for optimal functioning:

1. *Main Directory*: This directory should include both `main.py` and `library.py` files, as well as `MMIA_extraction.m` file (all three explained in [Sec.3](#revision)) and the downloaded MATLAB CDF patch directory explained before.
2. *MMIA Files Directory*: Where all MMIA *.cdf* files to compare are located (all together).
3. *Results Directory*: Where all results will be outputted.

>Important to notice that this directory's name should NOT include strings "*_s*", "*_e*" or "*_c*", as it may interfere with GLM filename system for *start*, *end* and *creation* times of the file.}

Another important step in order to download GLM *.nc* files is to [download the key for the GLM bucket](https://console.cloud.google.com/iam-admin/serviceaccounts?project=balma-280811&supportedpurview=project), store the JSON key wherever in the filesystem and define its path as an environment variable with the name `GOOGLE_APPLICATION_CREDENTIALS`.

> Example: Add line `export GOOGLE_APPLICATION_CREDENTIALS=<path-to-json-key>` to file `~/.bashrc` in Linux systems (or `~/.zshrc` on macOS)

## Fast Revision <a name="revision"></a>
### Global Structure <a name="structure"></a>
The global structure of the program consists on an executable `main.py` script which serves as a complete "list of tasks", being the data-input file and enumerating all steps by calling their respective functions. Those functions are found in `library.py` in order to keep the main script as clean as possible, being this second file just definitions of functions and having no need to enter it (unless willing to understand or change a process, of course).

A third MATLAB script called `MMIA_extraction.m` is needed for MMIA data extraction from *.cdf* to *.mat* files. In order to set it up and make it work with the main python script, it is needed to add the path to the MATLAB CDF patch at the beginning of the *.m* file (clearly marked inside the script).

### Input Data <a name="input"></a>
User input data section on `main.py` looks like the following code extract, and should be the only part of the script with user interaction.

```python
'''
###########################################################
##                     USER INPUT DATA                   ##
###########################################################
'''


### GENERAL ###

# Boolean variable for setting everything for the first execution
first_execution = True

# Boolean variable for generating plots
show_plots = False

# Boolean variable for pre-cross-correlated data
pre_xc = False

# Boolean variable for pre-converted to top cloud energy
pre_tce = True

# Boolean variable for pre-detected peaks
pre_detected_peaks = True

# Boolean variable for pre-studied peaks
pre_studied = False

# Boolean variable for just outputting results
just_results = False

# Boolean variable for pre-oredered events in directories
pre_event_directories = True

# Boolean variable for deleting non-important directories at the end of execution
delete_non_important_directories = False

# Path to Hard Disk (with all MMIA files and where to store all files)
#ssd_path = '/path/to/outputting/directory'

# Path where MMIA's .cdf files are located
#MMIA_files_path = '/path/to/mmia/files'

# Path to MATLAB executable
matlab_path = '/usr/local/MATLAB/R2021b/bin/matlab'


### GLM ###

# Time in seconds to analyze GLM before and after MMIA's time snippet
cropping_margin = 0.1

# GLM pixel size [km^2]
glm_pix_size = 8*8  # Colombia
#glm_pix_size = 78   # USA

# Plus of angle in latitude and longitude to snip GLM data
GLM_radius = 400 # [km]
angle_margin = GLM_radius / 111.11 # or a given value in degrees

# Minumum number of peaks to consider a valid event for GLM
glm_min_peak_num = 3

# Boolean variable for downloading GLM .nc files from Google Cloud Storage
pre_downloaded_GLM = True

# Boolean variable for pre-extracted files
pre_extracted_GLM = True

# Boolean variable for integrating GLM signals if not pre-done
pre_conditioned_GLM = True

# Boolean variable for defining TFG events (outputs GLM .txt's that could be cross-correlated)
tgf = True


### MMIA ###

# Boolean variable for pre-extracted files
pre_extracted_MMIA = True

# Boolean variable for conditioning MMIA data if not done before
pre_conditioned_MMIA = True

# Maximum length in seconds of each event
event_length = 2 # [s]

# Minimum time to consider as two sepparate events
split_window = 2000*0.00001

# Threshold for MMIA signal
mmia_threshold = 1.75   # [micro W / m^2]

# Minumum number of peaks to consider a valid event for MMIA
mmia_min_peak_num = 3


'''
###########################################################
##                END OF USER INPUT DATA                 ##
###########################################################
'''
```
The main variables for setting an execution are, first and foremost:

* `ssd_path`: String defining the path to the directory where all results will be outputted into
* `MMIA_files_path`: String defining the path to the directory containing all MMIA *.cdf* files
* `matlab_path`: String defining the path to MATLAB binary

While there are some other specific variables both for GLM and MMIA or for the functioning of the program being:

* `cropping_margin`: Time in seconds to crop GLM signal before and after MMIA signal's starting and ending time (for better match and cross-correlation).
* `glm_pix_size`: Pixel size for GLM, important for TCE conversion [km$^2$]
* `GLM_radius`: Radius in [km] to search for GLM signals with respect to MMIA coordinates (`angle_margin` is calculated after this value or can be directly given)
* `(glm/mmia)_min_peak_num`: Minimum number of peaks to find for GLM or MMIA, respectively. If found less (but not equal to) in any of them, the event is discarded.
* `event_length`: Maximum time for a MMIA event in order to classify its *.cdf* files in events [s]
* `split_window`: Minimum time with no signal inside a MMIA event to be considered two different events (because of the `event_length` definition)
* `mmia_threshold`: Maximum value of MMIA signal to be considered noise [microW/m^2]
* `show_plots`: Boolean variable for outputting plots throughout the whole program
* `delete_non_important_directories`: Boolean variable for deleting all non-necessary directories used in the execution of the program at its conclusion if not desired (saves storage)
* `tgf`: Boolean variable for stating that all MMIA *.cdf* files represent TGF's, outputting all generated GLM *.txt* files that could cross-correlate with MMIA

The other boolean variables are explained due to the contruction of the `main.py` scipt. The whole document consists of a series of steps as stated previously, which store their results just after finishing that step (storing MMIA conditioned data, then GLM conditioned data, then cross-correlated data, ...) so *if the previous steps were made correctly*, some parts of the code can be deactivated and necessary data will be uploaded from binaries, making the processing of data way faster if one step is in testing stage or a single little change was made. Those main sections and their corresponding boolean variable are:

* Order MMIA *.cdf* files in daily events (`pre_event_directories`)
* Extract MMIA data from *.cdf* files to *.mat* files (`pre_extracted_MMIA`)

And then, for every day with existing MMIA data:


* MMIA data conditioning (`pre_conditioned_MMIA`)
* GLM *.nc* files download (`pre_downloaded_GLM`)
* GLM data conditioning (`pre_conditioned_GLM`)
* Top Cloud Energy Conversion (`pre_tce`)
* GLM-MMIA Cross Correlation (`pre_xc`)
* GLM-MMIA peak detection and peak matching (`pre_detected_peaks`)
* Output results (`pre_studied`)


And some extra special variables:

* `first_execution`: Bool variable for executing the whole program (no need to set all previous boolean variables to False)
* `just_results`: Bool variable for jumping all processing and just computing results. Important to note that all previous steps should be correctly done, being handful if all processing was pre-done and more variables want to be extracted


### Output Data <a name="output"></a>
At the conclusion of the execution, some useful information remains storaged in the outputting directory. Among the most useful, and depending on outputting boolean variables:

* GLM downloaded *.nc* files for their corresponding ordered MMIA *.cdf* files
* GLM and MMIA extracted data (in *.txt* and *.mat* file formats, respectively)
* GLM and MMIA conditioned, TCE-converted signals, as *pickle* binaries
* GLM and MMIA cross-correlated signals, as *pickle* binaries, and their representations pre- and post- cross correlation
* MMIA delays for every event
* All GLM and MMIA peaks as well as their matching peaks and their representation in plots, per event
* Some basic statistics (including the most important variables for further statistical studies in a *.mat* format)

Giving special attention to the output *.mat* file called `vars_for_stats.mat`, where all important results can be found, its fields represent:

* `all_delays_t`: A vector containing all MMIA delays in seconds for every correlated event (MMIA was not noise-only and MMIA and GLM files contained data). *Positive values show GLM coming after MMIA, and negative values show GLM coming before MMIA signal*.
* `glm_mmia_matching_time_signal`: Four vectors (1st and 2nd for GLM, 3rd and 4th for MMIA) containing time (in seconds) and the signal peak value (in joules, after TCE conversion) of all matching peaks. Every row represents a different matching peak with time and signal values for GLM and the same for MMIA.
* `matching_time_distribution`: Twelve position vector counting the number of matching peaks for every 2-hour time segment of the day.
* `peak_relations`: Two vectors, where every row represents a correlated event with a minimum of *glm_min_peak_num* peaks for GLM and `mmia_min_peak_num` peaks for MMIA. First column represents the relation between matching peaks over all GLM peaks for that row-event, and the second column represents the same over all MMIA peaks for that event.

### First Execution Checklist <a name="checklist"></a>

As a matter of summary for making the script work properly for the first time:


1. Clone this repository (will be *main directory*)
2. Download [MATLAB CDF Patch](https://cdf.gsfc.nasa.gov/html/matlab_cdf_patch.html) and move it to *main directory* (with `main.py`, `library.py`, ...)
3. Enter `MMIA_extraction.m` and add path to the patch
4. Install all Python dependencies
5. [Download GLM JSON key](https://console.cloud.google.com/iam-admin/serviceaccounts?project=balma-280811&supportedpurview=project) and add it as an environment variable
6. Create a directory with all desirable MMIA *.cdf* files
7. Create a directory where all results will be outputted
8. Enter `main.py` and define paths to directories and to MATLAB binary, as well as all input variables
9. Execute `main.py`

### Flow Diagram for `main.py` <a name="flowchart"></a>
![Basic flow diagram for the `main.py` script](https://github.com/jaimefranm/GLM-ASIM/blob/master/code_chart.pdf)
