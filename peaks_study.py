import pickle
import matplotlib.pyplot as plt
import numpy as np
import library as TFG

#plt.rc('font', **{'family': 'serif', 'serif': ['latin modern roman']})

def study_delays(statistics_bin, show_plots, statistics_figures_path, matches, ssd_path):
    '''
    This function computes the average delay of GLM with respect to MMIA,
    both attending absolute values and real values.
    '''
    
    # Generating necessary global matrices
    
    delays = TFG.merge_days(matches, statistics_bin, 'delays')
    glm_avg = TFG.merge_days(matches, statistics_bin, 'glm_avg')
    mmia_avg = TFG.merge_days(matches, statistics_bin, 'mmia_avg')
    glm_std = TFG.merge_days(matches, statistics_bin, 'glm_std')
    mmia_std = TFG.merge_days(matches, statistics_bin, 'mmia_std')
    
    
    valid_trigger_counter = 0
    delay_list = []
    MMIA_delays_list = []
    GLM_delays_list = []
    no_delays_list = []
    trigger_counter = 0
    
    
    print('Computing average delay...')
    
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

    
    all_delays = np.zeros((len(delay_list),1))
    all_delays[:,0] = delay_list[:]
    if len(all_delays) != 0:
        avg_all = np.average(all_delays)
        std_all = np.std(all_delays)
    
    
    MMIA_delays = np.zeros((len(MMIA_delays_list),1))
    MMIA_delays[:,0] = MMIA_delays_list[:]
    if len(MMIA_delays) != 0:
        avg_negative = np.average(MMIA_delays)
        std_negative = np.std(MMIA_delays)
        
    # GLM delays is really MMIA sends early
    
    GLM_delays = np.zeros((len(GLM_delays_list),1))
    GLM_delays[:,0] = GLM_delays_list[:]
    if len(GLM_delays) != 0:
        avg_positive = np.average(GLM_delays)
        std_positive = np.std(GLM_delays)
    

    no_delays = np.zeros((len(no_delays_list),1))
    no_delays[:,0] = no_delays_list[:]
    
    valid_triggers = len(MMIA_delays) + len(GLM_delays) + len(no_delays)
    
    # Checking relationship between delay and signal average energy
    
    GLM_delay_signal = np.zeros((len(GLM_delays),3))
    MMIA_delay_signal = np.zeros((len(MMIA_delays),3))
    GLM_counter = 0
    MMIA_counter = 0
    
    for i in range(len(delays)):
        for j in range(len(delays[i])):
            if type(delays[i][j]) == int:
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
            
            
    # Saving results into a .txt
    f =  open(ssd_path + 'RESULTS.txt', 'w')
    f.write('******** RESULTS ********')
    f.write('\n')
    f.write('\n')
    f.write('\n')
    f.write('* TRIGGER INFO:')
    f.write('\n')
    f.write('    --> Total number of triggers processed: %d\n' % trigger_counter)
    f.write('    --> Number of valid triggers: %d AND %d\n' % (valid_trigger_counter, valid_triggers))
    f.write('    --> Percentage of valid over total triggers: %d%%\n' % valid_trigger_counter/trigger_counter * 100)
    f.write('\n')
    f.write('\n')
    f.write('* DELAY INFO:')
    f.write('\n')
    f.write('    --> Average delay in samples: %d +- %d\n' % (avg_all, std_all))
    f.write('    --> Average delay in seconds: %d +- %d\n' % (avg_all*0.00001, std_all*0.00001))
    f.write('    --> Number of MMIA delays: %d (%d%%)\n' % (len(MMIA_delays), len(MMIA_delays)/valid_triggers*100))
    f.write('    --> Average MMIA delay in samples: %d +- %d\n' % (avg_negative, std_negative))
    f.write('    --> Average MMIA delay in seconds: %d +- %d\n' % (avg_negative*0.00001, std_negative*0.00001))
    f.write('    --> Number of MMIA anticipations: %d (%d%%)\n' % (len(GLM_delays), len(GLM_delays)/valid_triggers*100))
    f.write('    --> Average MMIA anticipation in samples: %d +- %d\n' % (avg_positive, std_positive))
    f.write('    --> Average MMIA anticipation in seconds: %d +- %d\n' % (avg_positive*0.00001, std_positive*0.00001))
    f.write('    --> Number of no delays: %d (%d%%)\n' % (len(no_delays), len(no_delays)/valid_triggers*100))
    f.close()
    
    
    if show_plots == 1:
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
        
        # Average in delays
        fig = plt.figure(dpi=800)
        ax = fig.add_axes([0,0,1,1])
        langs = ['Avg. Delay', 'Avg. MMIA Delay', 'Avg. MMIA Anticipation']
        data = [avg_all*0.00001, avg_negative*0.00001, avg_positive*0.00001]
        plt.grid('on')
        colors = ['tab:blue', 'darksalmon', 'slategrey']
        ax.set_ylabel('Delay [s]')
        ax.bar(langs,data, color=colors, yerr=[std_all*0.00001,std_negative*0.00001,std_positive*0.00001])
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
        plt.scatter(GLM_delay_signal[:,1]*0.00001, GLM_delay_signal[:,0]*0.00001, color='black', marker='x')
        plt.xlabel('GLM snippet average radiance [J]')
        plt.ylabel('Snippet delay [s]')
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
        plt.scatter(MMIA_delay_signal[:,1]*0.00001, MMIA_delay_signal[:,0]*0.00001, color='r', marker='x')
        plt.xlabel(r'MMIA snippet average irradiance $\left[\dfrac{\mu W}{m^2}\right]$')
        plt.ylabel('Snippet delay [s]')
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
        plt.scatter(GLM_delay_signal[:,2]*0.00001, GLM_delay_signal[:,0]*0.00001, color='black', marker='x')
        plt.xlabel("GLM snippet's radiance stantard deviation [J]")
        plt.ylabel('Snippet delay [s]')
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
        plt.scatter(MMIA_delay_signal[:,2]*0.00001, MMIA_delay_signal[:,0]*0.00001, color='r', marker='x')
        plt.xlabel(r"MMIA snippet's irradiance stantard deviation $\left[\dfrac{\mu W}{m^2}\right]$")
        plt.ylabel('Snippet delay [s]')
        plt.grid('on')
        #plt.show()
        plt.savefig(statistics_figures_path + '/mmia_delay_vs_mmia_std_energy.pdf')
        # Clear the current axes
        plt.cla() 
        # Clear the current figure
        plt.clf() 
        # Closes all the figure windows
        plt.close('all')
                     
    print('Done')
    print(' ')

def more_statistics(peaks_bin, matches):
    
    # Generating necessary global matrices
    GLM_peaks = TFG.merge_days(matches, peaks_bin, 'glm')
    MMIA_peaks = TFG.merge_days(matches, peaks_bin, 'mmia')
    matching_peaks = TFG.merge_days(matches, peaks_bin, 'matching')
    
    
    counter = 0
    sum_rel_GLM = 0 # en mmia/todos
    sum_rel_MMIA = 0 # en glm/todos
    sum_GLM_peaks = 0 # total peaks
    sum_GLM_MMIA_peaks = 0
    sum_MMIA_peaks = 0
    sum_MMIA_GLM_peaks = 0
    
    for i in range(len(GLM_peaks)):
        for j in range(len(GLM_peaks[i])):
            if type(matching_peaks[i][j]) == list:
                counter = counter+1
                sum_rel_GLM = sum_rel_GLM + len(matching_peaks[i][j][1])/len(GLM_peaks[i][j])
                sum_rel_MMIA = sum_rel_MMIA + len(matching_peaks[i][j][1])/len(MMIA_peaks[i][j])
                
                sum_GLM_peaks = sum_GLM_peaks + len(GLM_peaks[i][j])
                sum_GLM_MMIA_peaks = sum_GLM_MMIA_peaks + len(matching_peaks[i][j][1])
                
                sum_MMIA_peaks = sum_MMIA_peaks + len(MMIA_peaks[i][j])
                sum_MMIA_GLM_peaks = sum_MMIA_GLM_peaks + len(matching_peaks[i][j][1])
                    
    avg_GLM_rel = sum_rel_GLM/counter
    avg_MMIA_rel = sum_rel_MMIA/counter
    
    avg_GLM_peaks = sum_GLM_peaks/counter
    avg_GLM_MMIA = sum_GLM_MMIA_peaks/counter
    
    avg_MMIA_peaks = sum_MMIA_peaks/counter
    avg_MMIA_GLM = sum_MMIA_GLM_peaks/counter
    
    f =  open(ssd_path + 'RESULTS.txt', 'a')
    avg_GLM_rel, avg_MMIA_rel, avg_GLM_peaks, avg_GLM_MMIA, avg_MMIA_peaks, avg_MMIA_GLM
    f.write('\n')
    f.write('\n')
    f.write('* PEAKS INFO:')
    f.write('\n')
    f.write('    --> Average delay in samples:')
    f.close()
    
    return [avg_GLM_rel, avg_MMIA_rel, avg_GLM_peaks, avg_GLM_MMIA, avg_MMIA_peaks, avg_MMIA_GLM]



ssd_path = ssd_path = '/Users/jaimemorandominguez/Desktop/test_descarga_GLM'
general_variables_path = ssd_path + '/general_variables_bin'
xcorr_bin = xcorr_bin = ssd_path + '/xcorr_bin'

show_plots = 1

f = open(general_variables_path+'/matches.pckl', 'rb')
matches = pickle.load(f)
f.close()

delays = TFG.merge_days(matches, xcorr_bin, 'delays')

'''
# Getting delay statistics
[total_snippets, avg_all, std_all, avg_MMIA_delay, std_MMIA_delay, avg_GLM_delay, std_GLM_delay, MMIA_delays, GLM_delays, no_delays] = TFG.study_delays(delays, GLM_xcorr, MMIA_xcorr, show_plots)

[avg_GLM_rel, avg_MMIA_rel, avg_GLM_peaks, avg_GLM_MMIA, avg_MMIA_peaks, avg_MMIA_GLM] = more_statistics(GLM_peaks, MMIA_peaks, GLM_xcorr, MMIA_xcorr, matching_peaks)
'''