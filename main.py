import library as TFG
import numpy as np
import pickle
import os

# Just for plot presentation in LaTeX Style (slows the program)
#plt.rc('font', **{'family': 'serif', 'serif': ['latin modern roman']})
'''

# TODO: Pasar a top cloud energy

'''

'''
###################################################
##               USER INPUT DATA                 ##
###################################################
'''


### GENERAL ###

# Boolean variable for setting everything for the first execution
first_execution = False

# Boolean variable for generating plots
show_plots = False

# Boolean variable for pre-cross-correlated data
pre_xc = False

# Boolean variable for pre-detected peaks
pre_detected_peaks = False

# Boolean variable for pre-studied peaks
pre_studied = False

# Boolean variable for just outputting results
just_results = True

# Boolean variable for pre-oredered events in directories
pre_event_directories = True

# Path to Hard Disk (with all MMIA files and where to store all files)
ssd_path = '/Volumes/Jaime_F_HD/mmia_2020'
#ssd_path = '/Users/jaimemorandominguez/Desktop/test_descarga_GLM'
#ssd_path = '/media/lrg'

# Path where MMIA's .cdf files are located
MMIA_files_path = '/Volumes/Jaime_F_HD/mmia_2020/mmia_20'
#MMIA_files_path = '/Users/jaimemorandominguez/Desktop/test_cdf'
#MMIA_files_path = '/media/lrg/mmia_20'

# Path to MATLAB executable
matlab_path = '/Applications/MATLAB_R2021b.app/bin/matlab'
#matlab_path = '/usr/local/bin/matlab'


### GLM ###

# Time in seconds to analyze GLM before and after MMIA's time snippet
cropping_margin = 0.05

# Plus of angle in latitude and longitude to snip GLM data
GLM_radius = 400 # [km]
angle_margin = GLM_radius / 111.11 # or a given value in degrees [º]

# Boolean variable for downloading GLM .nc files from Google Cloud Storage
pre_downloaded_GLM = True

# Boolean variable for pre-extracted files
pre_extracted_GLM = True

# Boolean variable for integrating GLM signals if not pre-done
pre_integrated_GLM = False


### MMIA ###

# Boolean variable for pre-extracted files
pre_extracted_MMIA = False

# Boolean variable for conditioning MMIA data if not done before
pre_conditioned_MMIA = False

# Maximum length in seconds of each event
event_length = 2 # [s]

# Threshold for MMIA signal
mmia_threshold = 1.75   # [micro W / m^2]


'''
###################################################
##           END OF USER INPUT DATA              ##
###################################################
'''


# New directories paths needed in the future

general_variables_path = ssd_path + '/general_variables_bin'
xcorr_bin = ssd_path + '/xcorr_bin'
xcorr_figures_path = ssd_path + '/xcorr_figures'
peaks_bin = ssd_path + '/peaks_bin'
peaks_figures_path = ssd_path + '/peaks_figures'
statistics_bin = ssd_path + '/stats_bin'
statistics_figures_path = ssd_path + '/stats_figures'

GLM_ordered_dir = ssd_path + '/glm_downl_nc_files'
GLM_ordered_outputs = ssd_path + '/glm_txt'
GLM_integrated_bin = ssd_path + '/glm_int_bin'

MMIA_mats_path = ssd_path + '/mmia_mat'
MMIA_filtered_bin = ssd_path + '/mmia_filt_bin'
path_to_mmia_dirs = ssd_path + '/mmia_dirs'
mmia_mats_files_path = ssd_path + '/mmia_mat'



###### MMIA'S EVENT CHARACTERIZATION ######

if first_execution == True:
    
    pre_event_directories = False
    pre_downloaded_GLM = False
    pre_extracted_GLM = False
    pre_integrated_GLM = False
    pre_extracted_MMIA = False
    pre_conditioned_MMIA = False
    pre_xc = False
    pre_detected_peaks = False
    pre_studied = False

    [matches, event_filenames] = TFG.get_MMIA_events(MMIA_files_path, event_length)
    
    os.system('mkdir ' + general_variables_path)
    
    # Saving 'matches' variable into a binary
    f = open(general_variables_path+'/matches.pckl', 'wb')
    pickle.dump(matches, f)
    f.close()
    
    # Saving 'event_filenames' variable into a binary
    f = open(general_variables_path+'/event_filenames.pckl', 'wb')
    pickle.dump(event_filenames, f)
    f.close()
    
else:
    
    # Uploading 'matches' from binary
    f = open(general_variables_path+'/matches.pckl', 'rb')
    matches = pickle.load(f)
    f.close()
    
    if just_results == False:
    
        # Uploading 'event_filenames' from binary
        f = open(general_variables_path+'/event_filenames.pckl', 'rb')
        event_filenames = pickle.load(f)
        f.close()


if just_results == False:
    
    ###### MMIA'S DATA ORDERING AND EXTRACTION ######
    
    # MMIA file ordering
    if pre_event_directories == False:
        
        TFG.create_MMIA_event_directories(matches, event_filenames, MMIA_files_path, ssd_path)
    
    # MMIA data extraction
    if pre_extracted_MMIA == False:
        os.system('mkdir ' + mmia_mats_files_path)
        TFG.extract_event_info(path_to_mmia_dirs, mmia_mats_files_path, matlab_path)


    #########################################################################################################
    #  From this point every step is made for every day with existing MMIA data until outputting statistics
    #########################################################################################################


    for day in range(len(matches)):
        print(' ')
        print('******************************')
        print('DAY %s, %d of %d' % (matches[day], day+1, len(matches)))
        print('******************************')
        print(' ')
        
        ###### MMIA'S DATA UPLOAD AND CONDITIONING ######

        # Uploading MMIA data and info
        print('All MMIA data was pre-extracted, uploading from %s...' % MMIA_mats_path)
        [mmia_raw, event_limits] = TFG.upload_MMIA_mats(ssd_path, event_filenames, matches, day)
        print('Done!\n')


        if pre_conditioned_MMIA == False:
            # Conditioning MMIA data for further analysis
            MMIA_filtered = TFG.condition_MMIA_data(mmia_raw, matches, show_plots, mmia_threshold, day)
            
            # Saving MMIA filtered data into binary
            print('Saving MMIA conditioned data for day %s...' % matches[day])
            if day == 0:
                os.system('mkdir ' + MMIA_filtered_bin)
            f = open(MMIA_filtered_bin + '/' + matches[day] + '.pckl', 'wb')
            pickle.dump(MMIA_filtered, f)
            f.close()
            print('Done!\n')
        else:
            print('MMIA data was pre-conditioned. Uploading from %s/%s.pckl...' % (MMIA_filtered_bin, matches[day]))
            f = open(MMIA_filtered_bin + '/' + matches[day] + '.pckl', 'rb')
            MMIA_filtered = pickle.load(f)
            f.close()
            print('Done!\n')
        del mmia_raw


        ########### GLM'S DATA DOWNLOAD, EXTRACTION, UPLOAD AND CONDITIONING ###########

        # Downloading GLM data from Google Cloud Services
        if pre_downloaded_GLM == False:
            
            TFG.download_GLM(ssd_path, event_filenames, MMIA_filtered, matches, day)

        # Extracting GLM data into event .txt files
        if pre_extracted_GLM == False:

            TFG.extract_GLM(GLM_ordered_dir, GLM_ordered_outputs, event_limits, matches, MMIA_filtered, angle_margin, cropping_margin, day)
        del event_limits

        # Uploading and integrating GLM signal

        if pre_integrated_GLM == False:
            
            # Unifying all data in a structure of matrices
            GLM_raw_data = TFG.unify_GLM_data(GLM_ordered_outputs, MMIA_filtered, matches, day)

            # Conditioning GLM data for further analysis
            GLM_data = TFG.condition_GLM_data(GLM_raw_data, matches, show_plots, day)
            
            del GLM_raw_data
            
            print('Saving GLM integrated data for day %s...' % matches[day])
            if day == 0:
                os.system('mkdir ' + GLM_integrated_bin)
            f = open(GLM_integrated_bin + '/' + matches[day] + '.pckl', 'wb')
            pickle.dump(GLM_data, f)
            f.close()
            print('Done!\n')
        else:
            print('GLM data was pre-integrated. Uploading from %s/%s.pckl...' % (GLM_integrated_bin, matches[day]))
            f = open(GLM_integrated_bin + '/' + matches[day] + '.pckl', 'rb')
            GLM_data = pickle.load(f)
            f.close()
            print('Done!\n')



        ########### CROSS-CORRELATION ###########

        if pre_xc == False:
            
            # Normalizing GLM data to cross-correlate with MMIA data
            print('Normalizing GLM data for day %s...' % matches[day])
            GLM_norm = [None] * len(GLM_data)

            for j in range(len(GLM_data)):
                if type(GLM_data[j]) == np.ndarray:
                    snippet = np.zeros((len(GLM_data[j]),2))
                    GLM_norm[j] = snippet
                    GLM_norm[j][:,0] = GLM_data[j][:,0]
                    GLM_norm[j][:,1] = TFG.normalize(GLM_data[j][:,1])
            print('Done!')
            print(' ')
            
            # Normalizing MMIA data to cross-correlate with GLM data
            print('Normalizing MMIA data for day %s...' % matches[day])
            MMIA_norm = [None] * len(MMIA_filtered)

            for j in range(len(MMIA_filtered)):
                if type(MMIA_filtered[j]) == np.ndarray:
                    snippet = np.zeros((len(MMIA_filtered[j]),2))
                    MMIA_norm[j] = snippet
                    MMIA_norm[j][:,0] = MMIA_filtered[j][:,0]
                    MMIA_norm[j][:,1] = TFG.normalize(MMIA_filtered[j][:,1])
            print('Done!')
            print(' ')

            # Cross-correlating snippets
            show_plots = True
            
            [GLM_xcorr, MMIA_xcorr, GLM_xcorr_norm, MMIA_xcorr_norm, delays] = TFG.cross_correlate_GLM_MMIA(GLM_data, MMIA_filtered, GLM_norm, MMIA_norm, matches, show_plots, day, xcorr_figures_path)
            
            # Calculating some values for further use
            [GLM_avg, MMIA_avg, GLM_std, MMIA_std] = TFG.get_ministats(GLM_xcorr, MMIA_xcorr)
            
            show_plots = False
            
            del GLM_data
            del MMIA_filtered
            del GLM_norm
            del MMIA_norm
            
            # Saving cross-correlated data
            print('Saving cross-correlated data for day %s...' % matches[day])
            if day == 0:
                os.system('mkdir ' + xcorr_bin)
                os.system('mkdir ' + statistics_bin)
            
            # Saving correlated signals
            f = open(xcorr_bin + '/' + matches[day] + '_signals.pckl', 'wb')
            pickle.dump([GLM_xcorr, MMIA_xcorr, GLM_xcorr_norm, MMIA_xcorr_norm], f)
            f.close()
            
            # Saving results in different binaries to make easier further use
            
            # Saving delays
            f = open(statistics_bin + '/' + matches[day] + '_delays.pckl', 'wb')
            pickle.dump(delays, f)
            f.close()
            
            # Saving GLM_avg
            f = open(statistics_bin + '/' + matches[day] + '_glm_avg.pckl', 'wb')
            pickle.dump(GLM_avg, f)
            f.close()
            
            # Saving MMIA_avg
            f = open(statistics_bin + '/' + matches[day] + '_mmia_avg.pckl', 'wb')
            pickle.dump(MMIA_avg, f)
            f.close()
            
            # Saving GLM_std
            f = open(statistics_bin + '/' + matches[day] + '_glm_std.pckl', 'wb')
            pickle.dump(GLM_std, f)
            f.close()
            
            # Saving MMIA_std
            f = open(statistics_bin + '/' + matches[day] + '_mmia_std.pckl', 'wb')
            pickle.dump(MMIA_std, f)
            f.close()
            
            print('Done!\n')
            
            del delays
            del GLM_avg
            del MMIA_avg
            del GLM_std
            del MMIA_std
            
        else:
            print('GLM and MMIA data for day %s was pre-correlated. Uploading from %s/%s.pckl...' % (matches[day], xcorr_bin, matches[day]))
            
            # Uploading signals
            f = open(xcorr_bin + '/' + matches[day] + '_signals.pckl', 'rb')
            [GLM_xcorr, MMIA_xcorr, GLM_xcorr_norm, MMIA_xcorr_norm] = pickle.load(f)
            f.close()
            


        ########### PEAK DETECTION AND COMPARISON ###########

        if pre_detected_peaks == False:
            
            # Creating directories for figures
            
            peaks_path = peaks_figures_path + '/all_peaks_figures'
            match_figs_path = peaks_figures_path + '/matching_peaks'
            if day == 0:
                os.system('mkdir ' + peaks_figures_path)
                os.system('mkdir ' + peaks_path)
                os.system('mkdir ' + match_figs_path)
            
            # Getting peaks from cross-correlated signals
            show_plots = True
            [GLM_peaks, MMIA_peaks] = TFG.get_GLM_MMIA_peaks(GLM_xcorr, MMIA_xcorr, GLM_xcorr_norm, MMIA_xcorr_norm, matches, show_plots, day, peaks_path)

            
            # Getting matching peaks
            
            matching_peaks = TFG.get_peak_matches(GLM_xcorr, GLM_xcorr_norm, MMIA_xcorr, MMIA_xcorr_norm, GLM_peaks, MMIA_peaks, show_plots, matches, day, match_figs_path)
            show_plots = False
            
            print('Saving peak positions for day %s...' % matches[day])
            if day == 0:
                os.system('mkdir ' + peaks_bin)
            
            # Saving results in different binaries to make easier further use
            
            # Saving GLM_peaks
            f = open(peaks_bin + '/' + matches[day] + '_glm.pckl', 'wb')
            pickle.dump(GLM_peaks, f)
            f.close()
            
            # Saving MMIA_peaks
            f = open(peaks_bin + '/' + matches[day] + '_mmia.pckl', 'wb')
            pickle.dump(MMIA_peaks, f)
            f.close()
            
            # Saving matching_peaks
            f = open(peaks_bin + '/' + matches[day] + '_matching.pckl', 'wb')
            pickle.dump(matching_peaks, f)
            f.close()
            
            print('Done!')
            print(' ')
            
            del GLM_peaks
            del MMIA_peaks
            del matching_peaks
        
        del GLM_xcorr
        del GLM_xcorr_norm
        del MMIA_xcorr
        del MMIA_xcorr_norm
        
    del event_filenames



########### OUTPUTTING VALUABLE STATS ###########

if pre_studied == False:

    # Getting delay statistics

    # Making necessary directories
    os.system('mkdir ' + statistics_figures_path)

    show_plots = True
    TFG.study_delays(statistics_bin, show_plots, statistics_figures_path, matches, ssd_path)

    TFG.more_statistics(peaks_bin, matches, ssd_path)
    show_plots = False
