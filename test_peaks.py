import library as TFGimport numpy as npimport matplotlib.pyplot as pltimport pickleprint('Uploading cross-correlated data...')f = open('GLM_MMIA_xcorr_data.pckl', 'rb')[matches, GLM_xcorr, MMIA_xcorr, GLM_xcorr_norm, MMIA_xcorr_norm] = pickle.load(f)f.close()print('Done!\n')show_plots = True# Getting peaks from cross-correlated signals[GLM_peaks, MMIA_peaks] = TFG.get_GLM_MMIA_peaks(GLM_xcorr, MMIA_xcorr, GLM_xcorr_norm, MMIA_xcorr_norm, matches, show_plots)f = open('peaks_data.pckl', 'wb')pickle.dump([GLM_peaks, MMIA_peaks], f)f.close()