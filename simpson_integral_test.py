import math
import numpy as np

def simpson_integral(signal, t1, t_end):
    
    dt_bin = 0.002 # Intervalos las muestras cada 2 ms
    cnt = 1
    int_signal = [None] * 2
    int_signal[0] = []
    int_signal[1] = []
    
    while t1 <= t_end:
        
        n = 1
        t_bin = range(t1, t1+dt_bin, dt_bin/n) # Subintervalos de division de la muestra
        int_signal = np.zeros((len(t_bin), 2))
        
        for i in range(len(t_bin)-1):

            f_in_signal = np.where(signal[:,0] >= t_bin[i] and signal[:,0] < t_bin[i+1])

            if len(f_in_signal) == 0:
                
                int_signal[0].append(t_bin[0])
                int_signal[1].append(1e-11)
                #int_signal[cnt,0] = t_bin[0] # Tiempo de la integral como el primedio
                #int_signal[cnt,1] = 1e-11
                
            else:
                
                int_signal[0].append(t_bin[1])
                int_signal[1].append(np.trapz(signal[f_in_signal,4])*1e-5)
                #int_signal[cnt,0] = t_bin[1] # Van der velde et al 2020
                #int_signal[cnt,1] = np.trapz(signal[f_in_signal,4])*1e-5 # 777 nm uJ m-2 in a 10 us exposure integrated with trapezium rule over bin of 2 ms
            
            cnt = cnt+1

        t1 = t1 + dt_bin
        
        b = np.array(int_signal)

        int_final_signal = np.zeros((b.shape[1], 2))
        int_final_signal[:,0] = b[0,:]
        int_final_signal[:,1] = b[1,:]
        
    return int_final_signal


def top_cloud_energy(GLM_xcorr, MMIA_xcorr):
    
    glm_tce = [None] * len(GLM_xcorr)
    mmia_tce = [None] * len(MMIA_xcorr)
    
    for i in range(len(GLM_xcorr)): # For every event
        if type(GLM_xcorr[i]) == np.ndarray:
            
            # GLM
            glm_pix_size = 8*8 # [km^2]
            GLM_cloud_E = GLM_xcorr[i]
            GLM_cloud_E[:,1] = 6.612 * (GLM_cloud_E[:,1]*1e15) * glm_pix_size
            glm_tce[i] = GLM_cloud_E
            del GLM_cloud_E

            # MMIA
            t1 = MMIA_xcorr[i][0,0]
            t_end = MMIA_xcorr[i][0,-1]
            #n = 10000 # (????)
            
            # Computing a Simpson integral over MMIA signal
            MMIA_cloud_E = simpson_integral(MMIA_xcorr[i], t1, t_end)

            MMIA_cloud_E[:,1] = MMIA_cloud_E[:,1]*1e-6*math.pi*400e3^2 #(Van der velde et al 2020)
            mmia_tce[i] = MMIA_cloud_E
            del MMIA_cloud_E
    
    return [glm_tce, mmia_tce]











ssd_path = '/Volumes/Jaime_F_HD/mmia_2020'
xcorr_bin = ssd_path + '/xcorr_bin'

