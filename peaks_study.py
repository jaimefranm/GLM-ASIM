import pickleimport matplotlib.pyplot as pltimport numpy as npimport TFG_library_swap as TFGplt.rc('font', **{'family': 'serif', 'serif': ['latin modern roman']})def get_peak_statistics(GLM_xcorr, MMIA_xcorr, GLM_peaks, MMIA_peaks, show_plots, matches):        common_peaks = [None] * len(GLM_xcorr)    for i in range(len(GLM_xcorr)):        snip = [None] * len(GLM_xcorr[i])        common_peaks[i] = snip        for i in range(len(GLM_xcorr)):        for j in range(len(GLM_xcorr[i])):                        print('Checking date %s snippet %d...' % (matches[i], j))            if 3<9:#if i!=3 and j!=19:                 if type(GLM_peaks[i][j]) == np.ndarray and type(MMIA_peaks[i][j]) == np.ndarray:                    if len(GLM_peaks[i][j]) >= 1 and len(MMIA_peaks[i][j]) >= 1:                                            # Checking GLM peaks in MMIA                                    GLM_peaks_in_MMIA = np.zeros((len(GLM_peaks[i][j]),2))                                    for k in range(len(GLM_peaks[i][j])):                            current_peak_pos = GLM_peaks[i][j][k]                            if current_peak_pos in MMIA_peaks[i][j][:]: # If the peak exists in the same position                                GLM_peaks_in_MMIA[k][0] = 1                                GLM_peaks_in_MMIA[k][1] = int(np.where(MMIA_peaks[i][j][:] == current_peak_pos)[0][0])                            else: # If the peak does not exist in the same position or does not exist                                counter = 1                                while counter <= 500: # Searching in a range of 1000 samples                                    pos_left = current_peak_pos - counter                                    pos_right = current_peak_pos + counter                                                            if pos_left >= MMIA_peaks[i][j][0] and pos_right <= MMIA_peaks[i][j][-1] and GLM_peaks_in_MMIA[k][0] != 1:                                                                    if ((MMIA_peaks[i][j][:] == pos_left).any() == True) and ((MMIA_peaks[i][j][:] == pos_right).any() == False):                                            GLM_peaks_in_MMIA[k][0] = 1                                            GLM_peaks_in_MMIA[k][1] = int(np.where(MMIA_peaks[i][j][:] == pos_left)[0][0])                                                                            elif ((MMIA_peaks[i][j][:] == pos_left).any() == False) and ((MMIA_peaks[i][j][:] == pos_right).any() == True):                                            GLM_peaks_in_MMIA[k][0] = 1                                            GLM_peaks_in_MMIA[k][1] = int(np.where(MMIA_peaks[i][j][:] == pos_right)[0][0])                                                                            elif ((MMIA_peaks[i][j][:] == pos_left).any() == True) and ((MMIA_peaks[i][j][:] == pos_right).any() == True):                                            GLM_peaks_in_MMIA[k][0] = 1                                            possible_sides = [pos_left, pos_right]                                            winner = possible_sides[int(round(np.random.random()))]                                            GLM_peaks_in_MMIA[k][1] = int(np.where(MMIA_peaks[i][j][:] == winner)[0][0])                                            # If both are false, do nothing (counter = counter+1)                                                                elif pos_left < MMIA_peaks[i][j][0] and pos_right <= MMIA_peaks[i][j][-1] and GLM_peaks_in_MMIA[k][0] != 1:                                        # pos_left can't be True                                        if pos_right in MMIA_peaks[i][j][:]:                                            GLM_peaks_in_MMIA[k][0] = 1                                            GLM_peaks_in_MMIA[k][1] = int(np.where(MMIA_peaks[i][j][:] == pos_right)[0][0])                                                                    elif pos_left >= MMIA_peaks[i][j][0] and pos_right > MMIA_peaks[i][j][-1] and GLM_peaks_in_MMIA[k][0] != 1:                                        # pos_right can't be True                                        if pos_left in MMIA_peaks[i][j][:]:                                            GLM_peaks_in_MMIA[k][0] = 1                                            GLM_peaks_in_MMIA[k][1] = int(np.where(MMIA_peaks[i][j][:] == pos_left)[0][0])                                                            counter = counter+1                                                        GLM_real_peaks_pos_pos_in_MMIA = []                        MMIA_real_peaks_pos_pos_in_GLM = []                                        for k in range(len(GLM_peaks_in_MMIA)):                            if GLM_peaks_in_MMIA[k,0] != 0:                                GLM_real_peaks_pos_pos_in_MMIA.append(k)                            if GLM_peaks_in_MMIA[k,0] != 0:                                MMIA_real_peaks_pos_pos_in_GLM.append(int(GLM_peaks_in_MMIA[k,1]))                                                common_peaks[i][j] = [GLM_real_peaks_pos_pos_in_MMIA, MMIA_real_peaks_pos_pos_in_GLM]                                                if show_plots == 1:                            # GLM singal with common peaks                            plt.figure()                            plt.plot(GLM_xcorr[i][j][:,1], color = 'black', linewidth=0.5)                            plt.plot(GLM_peaks[i][j], GLM_xcorr[i][j][GLM_peaks[i][j],1], "*", color='b')                            plt.plot(GLM_peaks[i][j][GLM_real_peaks_pos_pos_in_MMIA], GLM_xcorr[i][j][GLM_peaks[i][j][GLM_real_peaks_pos_pos_in_MMIA],1], "*", color='gold')                            #plt.title('GLM peaks on day %d, snippet %d' % (int(matches[i]), j))                            plt.xlabel('Samples')                            plt.ylabel('Radiance [J]')                            plt.grid('on')                            plt.show()                            plt.figure()                            plt.plot(MMIA_xcorr[i][j][:,1], color = 'r', linewidth=0.5)                            plt.plot(MMIA_peaks[i][j], MMIA_xcorr[i][j][MMIA_peaks[i][j],1], "*", color='gold')                            plt.plot(MMIA_peaks[i][j][MMIA_real_peaks_pos_pos_in_GLM], MMIA_xcorr[i][j][GLM_peaks[i][j][GLM_real_peaks_pos_pos_in_MMIA],1], "*", color='b')                            #plt.title('MMIA peaks on day %d, snippet %d' % (int(matches[i]), j))                            plt.xlabel('Samples')                            plt.ylabel(r"Irradiance $\left[\dfrac{\mu W}{m^2}\right]$")                            plt.grid('on')                            plt.show()    return common_peaksdef more_statistics(GLM_peaks, MMIA_peaks, GLM_xcorr, MMIA_xcorr, common_peaks):        counter = 0    sum_rel_GLM = 0 # en mmia/todos    sum_rel_MMIA = 0 # en glm/todos    sum_GLM_peaks = 0 # total peaks    sum_GLM_MMIA_peaks = 0    sum_MMIA_peaks = 0    sum_MMIA_GLM_peaks = 0        for i in range(len(GLM_peaks)):        for j in range(len(GLM_peaks[i])):            if 3<9:                if type(common_peaks[i][j]) == list:                    counter = counter+1                    sum_rel_GLM = sum_rel_GLM + len(common_peaks[i][j][1])/len(GLM_peaks[i][j])                    sum_rel_MMIA = sum_rel_MMIA + len(common_peaks[i][j][1])/len(MMIA_peaks[i][j])                                        sum_GLM_peaks = sum_GLM_peaks + len(GLM_peaks[i][j])                    sum_GLM_MMIA_peaks = sum_GLM_MMIA_peaks + len(common_peaks[i][j][1])                                        sum_MMIA_peaks = sum_MMIA_peaks + len(MMIA_peaks[i][j])                    sum_MMIA_GLM_peaks = sum_MMIA_GLM_peaks + len(common_peaks[i][j][1])                        avg_GLM_rel = sum_rel_GLM/counter    avg_MMIA_rel = sum_rel_MMIA/counter        avg_GLM_peaks = sum_GLM_peaks/counter    avg_GLM_MMIA = sum_GLM_MMIA_peaks/counter        avg_MMIA_peaks = sum_MMIA_peaks/counter    avg_MMIA_GLM = sum_MMIA_GLM_peaks/counter        return [avg_GLM_rel, avg_MMIA_rel, avg_GLM_peaks, avg_GLM_MMIA, avg_MMIA_peaks, avg_MMIA_GLM]show_plots = 1print('Uploading data...')f = open('xcorr_data.pckl', 'rb')[GLM_xcorr, MMIA_xcorr, GLM_norm, MMIA_norm, GLM_peaks, MMIA_peaks, matches, delays] = pickle.load(f)f.close()print('Done!')'''f = open('data.pckl', 'rb')[common_peaks] = pickle.load(f)f.close()print('Done!')'''# Comparing peaks of cross-correlated signals#[GLM_peaks, MMIA_peaks] = TFG.get_GLM_MMIA_peaks(GLM_xcorr, MMIA_xcorr, matches, show_plots)show_plots = 0# Getting delay statistics#[total_snippets, avg_all, std_all, avg_MMIA_delay, std_MMIA_delay, avg_GLM_delay, std_GLM_delay, MMIA_delays, GLM_delays, no_delays] = TFG.study_delays(delays, GLM_xcorr, MMIA_xcorr, show_plots)# Getting peak statisticscommon_peaks = get_peak_statistics(GLM_xcorr, MMIA_xcorr, GLM_peaks, MMIA_peaks, show_plots, matches)#[avg_GLM_rel, avg_MMIA_rel, avg_GLM_peaks, avg_GLM_MMIA, avg_MMIA_peaks, avg_MMIA_GLM] = more_statistics(GLM_peaks, MMIA_peaks, GLM_xcorr, MMIA_xcorr, common_peaks)