import math

def simpson_integral(MMIA_signal,t1,t_end, n):
    
    dt_bin=0.002 # Intervalos las muestras cada 2 ms
    cnt=1
    while t1<=t_end:
        t_bin=t1:dt_bin/n:t1+dt_bin  # Subintervalos de division de la muestra
        for i in range(len(t_bin)-1):

            f_in_MMIA = find(MMIA_(:,1)>=t_bin(i) & MMIA_(:,1)<(t_bin(i+1)))

            if isempty(f_in_MMIA):
                int_MMIA(cnt,1)=(t_bin(1)) % # Tiempo de la integral como el primedio
                int_MMIA(cnt,2)=1e-11
                int_MMIA(cnt,3)=1e-11
            else:
                int_MMIA(cnt,1)=(t_bin(1)); % Van der velde et al 2020
    %             int_MMIA(cnt,2)=sum(MMIA_(f_in_MMIA,4))*1e-5
    %             int_MMIA(cnt,3)=sum(MMIA_(f_in_MMIA,2))*1e-5
                int_MMIA(cnt,2)=trapz(MMIA_(f_in_MMIA,4))*1e-5 # 777 nm uJ m-2 in a 10 us exposure integrated with trapezium rule over bin of 2 ms.
                int_MMIA(cnt,3)=trapz(MMIA_(f_in_MMIA,2))*1e-5 # 337 nm
                int_MMIA(cnt,4)=trapz(MMIA_(f_in_MMIA,3))*1e-5 # 180 nm
            
            cnt=cnt+1

        t1=t1+dt_bin



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
            n = 10000 # (????)
            
            # Computing a Simpson integral over MMIA signal
            MMIA_cloud_E = simpson_integral(MMIA_xcorr[i], t1, t_end, n)

            MMIA_cloud_E[:,1] = MMIA_cloud_E[:,1]*1e-6*math.pi*400e3^2 #(Van der velde et al 2020)
            mmia_tce[i] = MMIA_cloud_E
            del MMIA_cloud_E
    
    return [glm_tce, mmia_tce]











ssd_path = '/Volumes/Jaime_F_HD/mmia_2020'
xcorr_bin = ssd_path + '/xcorr_bin'

