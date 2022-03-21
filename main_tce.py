import library_tce as TFG
import numpy as np
import pickle
import os

# Just for plot presentation in LaTeX Style (slows the program)
#plt.rc('font', **{'family': 'serif', 'serif': ['latin modern roman']})

# TODO: Empezar con el TCE conversion al principio del programa
# TODO: Comprobar Pixel size y Factor división angle_margin para USA


'''
###########################################################
##                     USER INPUT DATA                   ##
###########################################################
'''


### GENERAL ###

# Boolean variable for setting everything for the first execution
first_execution = False

# Boolean variable for generating plots
show_plots = False

# Boolean variable for pre-cross-correlated data
pre_xc = True

# Boolean variable for pre-converted to top cloud energy
pre_tce = False

# Boolean variable for pre-detected peaks
pre_detected_peaks = True

# Boolean variable for pre-studied peaks
pre_studied = True

# Boolean variable for just outputting results
just_results = False

# Boolean variable for pre-oredered events in directories
pre_event_directories = True

# Path to Hard Disk (with all MMIA files and where to store all files)
#ssd_path = '/Volumes/Jaime_F_HD/mmia_2020'
ssd_path = '/Users/jaimemorandominguez/Desktop/special_tests/van_der_velde_results'
#ssd_path = '/home/lrg/Desktop/TCEpreXCORR'
#ssd_path = '/home/lrg/Desktop/USA'

# Path where MMIA's .cdf files are located
#MMIA_files_path = '/Volumes/Jaime_F_HD/mmia_2020/mmia_20'
MMIA_files_path = '/Users/jaimemorandominguez/Desktop/special_tests/van_der_velde/mmia_cdf'
#MMIA_files_path = '/media/lrg/colombia_2020/mmia_20'
#MMIA_files_path = '/media/lrg/mmia_triggers_usa'

# Path to MATLAB executable
matlab_path = '/Applications/MATLAB_R2021b.app/bin/matlab'
#matlab_path = '/usr/local/MATLAB/R2021b/bin/matlab'


### GLM ###

# Time in seconds to analyze GLM before and after MMIA's time snippet
cropping_margin = 0.5

# GLM pixel size [km²]
glm_pix_size = 8*8  # Colombia
#glm_pix_size = 78   # USA

# Plus of angle in latitude and longitude to snip GLM data
GLM_radius = 400 # [km]
angle_margin = GLM_radius / 111.11 # or a given value in degrees

# Boolean variable for downloading GLM .nc files from Google Cloud Storage
pre_downloaded_GLM = True

# Boolean variable for pre-extracted files
pre_extracted_GLM = True

# Boolean variable for integrating GLM signals if not pre-done
pre_conditioned_GLM = True


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


'''
###########################################################
##                END OF USER INPUT DATA                 ##
###########################################################
'''


# New directories paths needed in the future

general_variables_path = ssd_path + '/general_variables_bin'
xcorr_bin = ssd_path + '/xcorr_bin'
xcorr_figures_path = ssd_path + '/xcorr_figures'
tce_bin = ssd_path + '/tce_bin'
tce_figures_path = ssd_path + '/tce_figures'
peaks_bin = ssd_path + '/peaks_bin'
peaks_figures_path = ssd_path + '/peaks_figures'
statistics_bin = ssd_path + '/stats_bin'
statistics_figures_path = ssd_path + '/stats_figures'

GLM_ordered_dir = ssd_path + '/glm_downl_nc_files'
GLM_ordered_outputs = ssd_path + '/glm_txt'
GLM_conditioned_bin = ssd_path + '/glm_cond_bin'

MMIA_mats_path = ssd_path + '/mmia_mat'
MMIA_filtered_bin = ssd_path + '/mmia_filt_bin'
path_to_mmia_dirs = ssd_path + '/mmia_dirs'
mmia_mats_files_path = ssd_path + '/mmia_mat'



###### MMIA'S EVENT CHARACTERIZATION ######

if first_execution == True:
    
    pre_event_directories = False
    pre_downloaded_GLM = False
    pre_extracted_GLM = False
    pre_conditioned_GLM = False
    pre_extracted_MMIA = False
    pre_conditioned_MMIA = False
    pre_tce = False
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
        TFG.extract_MMIA_event_info(path_to_mmia_dirs, mmia_mats_files_path, matlab_path)


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

            # Splitting those signals where more than one event was found
            [mmia_raw, event_limits, event_filenames[day]] = TFG.split_MMIA_events(mmia_raw, event_limits, event_filenames[day], split_window)

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

        if pre_conditioned_GLM == False:
            
            # Unifying all data in a structure of matrices
            GLM_raw_data = TFG.unify_GLM_data(GLM_ordered_outputs, MMIA_filtered, matches, day)

            # Conditioning GLM data for further analysis
            GLM_data = TFG.condition_GLM_data(GLM_raw_data, matches, show_plots, day)
            
            del GLM_raw_data
            
            print('Saving GLM conditioned data for day %s...' % matches[day])
            if day == 0:
                os.system('mkdir ' + GLM_conditioned_bin)
            f = open(GLM_conditioned_bin + '/' + matches[day] + '.pckl', 'wb')
            pickle.dump(GLM_data, f)
            f.close()
            print('Done!\n')
        else:
            print('GLM data was pre-conditioned. Uploading from %s/%s.pckl...' % (GLM_conditioned_bin, matches[day]))
            f = open(GLM_conditioned_bin + '/' + matches[day] + '.pckl', 'rb')
            GLM_data = pickle.load(f)
            f.close()
            print('Done!\n')


        ########### CONVERSION TO TOP CLOUD ENERGY ###########
        
        if pre_tce == False:

            if day == 0:
                os.system('mkdir ' + tce_bin)
                os.system('mkdir ' + tce_figures_path)
            
            # Convert GLM and MMIA data to Top Cloud Energy data
            show_plots = True
            print('Converting GLM and MMIA conditioned instrumental data into Top Cloud Energy for day %s\n' % matches[day])
            [glm_tce, mmia_tce] = TFG.top_cloud_energy(GLM_data, MMIA_filtered, matches[day], show_plots, tce_figures_path, glm_pix_size)
            print('Done!\n')
            show_plots = False
            
            # Saving GLM and MMIA Top Cloud Energy data
            print('Saving TCE data for day %s...\n' % matches[day])
            f = open(tce_bin + '/' + matches[day] + '.pckl', 'wb')
            pickle.dump([glm_tce, mmia_tce], f)
            f.close()

        else:
            print('Top Cloud Energy values were pre-calculated. Uploading from %s/%s.pckl...' % (tce_bin, matches[day]))

            # Uploading TCE signals
            f = open(tce_bin + '/' + matches[day] + '.pckl', 'rb')
            [glm_tce, mmia_tce] = pickle.load(f)
            f.close()

        #del GLM_data
        del MMIA_filtered



        ########### CROSS-CORRELATION ###########

        if pre_xc == False:

            # Normalizing GLM data to cross-correlate with MMIA data
            print('Normalizing GLM data for day %s...' % matches[day])
            GLM_norm = [None] * len(glm_tce)

            for j in range(len(glm_tce)):
                if type(glm_tce[j]) == np.ndarray:
                    event = np.zeros((len(glm_tce[j]),2))
                    GLM_norm[j] = event
                    GLM_norm[j][:,0] = glm_tce[j][:,0]
                    GLM_norm[j][:,1] = TFG.normalize(glm_tce[j][:,1])
            print('Done!')
            print(' ')

            # Normalizing MMIA data to cross-correlate with GLM data
            print('Normalizing MMIA data for day %s...' % matches[day])
            MMIA_norm = [None] * len(mmia_tce)

            for j in range(len(mmia_tce)):
                if type(mmia_tce[j]) == np.ndarray:
                    event = np.zeros((len(mmia_tce[j]),2))
                    MMIA_norm[j] = event
                    MMIA_norm[j][:,0] = mmia_tce[j][:,0]
                    MMIA_norm[j][:,1] = TFG.normalize(mmia_tce[j][:,1])
            print('Done!')
            print(' ')

            # Cross-correlating snippets
            show_plots = True

            [GLM_xcorr, MMIA_xcorr, GLM_xcorr_norm, MMIA_xcorr_norm, delays] = TFG.cross_correlate_GLM_MMIA(glm_tce, mmia_tce, GLM_norm, MMIA_norm, matches, show_plots, day, xcorr_figures_path)

            # Calculating some values for further use
            [GLM_avg, MMIA_avg, GLM_std, MMIA_std] = TFG.get_ministats(GLM_xcorr, MMIA_xcorr)

            show_plots = False

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
            print('GLM and MMIA data for day %s was pre-correlated. Uploading from %s/%s.pckl...\n' % (matches[day], xcorr_bin, matches[day]))

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

        else:
            print('GLM and MMIA peaks were pre-found for day %s. Uploading from %s/%s.pckl...\n' % (matches[day], peaks_bin, matches[day]))

            # Uploading GLM peaks
            f = open(peaks_bin + '/' + matches[day] + '_glm.pckl', 'rb')
            GLM_peaks = pickle.load(f)
            f.close()

            # Uploading MMIA peaks
            f = open(peaks_bin + '/' + matches[day] + '_mmia.pckl', 'rb')
            MMIA_peaks = pickle.load(f)
            f.close()

            # Uploading GLM peaks
            f = open(peaks_bin + '/' + matches[day] + '_matching.pckl', 'rb')
            matching_peaks = pickle.load(f)
            f.close()


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
