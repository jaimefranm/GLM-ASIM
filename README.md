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
-- In Progress --
### Input Data <a name="input"></a>
-- In Progress --
### Output Data <a name="output"></a>
-- In Progress --
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