% INTEGRAL DE SIMPSON PARA MMIA, LIS Y GLM CADA 2 ms
%
function [int_LIS,int_GLM,int_MMIA]=integral_GLM_LIS_MMIA(GLM_in,LIS_in,MMIA_,t1,t_end)

dt_bin=0.002; % Intervalos las muestras cada 2 ms
cnt=1;
while t1<=t_end
    n=1;  % N�mero de subintervalos
    t_bin = t1 : dt_bin/n : t1+dt_bin;  % Subintervalos de divisi�n de la muestra
    for i=1:length(t_bin)-1

        f_in_MMIA=find(MMIA_(:,1)>=t_bin(i) & MMIA_(:,1)<(t_bin(i+1)));

        if isempty(f_in_MMIA)
            int_MMIA(cnt,1)=(t_bin(1)); % # Tiempo de la integral como el primedio
            int_MMIA(cnt,2)=1e-11;
            int_MMIA(cnt,3)=1e-11;
        else
            int_MMIA(cnt,1)=(t_bin(1)); % Van der velde et al 2020
%             int_MMIA(cnt,2)=sum(MMIA_(f_in_MMIA,4))*1e-5;
%             int_MMIA(cnt,3)=sum(MMIA_(f_in_MMIA,2))*1e-5;
            int_MMIA(cnt,2)=trapz(MMIA_(f_in_MMIA,4))*1e-5; % 777 nm uJ m-2 in a 10 us exposure integrated with trapezium rule over bin of 2 ms.
            int_MMIA(cnt,3)=trapz(MMIA_(f_in_MMIA,2))*1e-5; % 337 nm
            int_MMIA(cnt,4)=trapz(MMIA_(f_in_MMIA,3))*1e-5; % 180 nm
        end
        
        cnt=cnt+1;
    end
    t1=t1+dt_bin;
end
end


%% GLM
glm_pix_size = 8*8;% (8kmx8km)
GLM_cloud_E = GLM_in; % Datos GLM en la regi�n  
GLM_cloud_E(:,end) = 6.612 * (GLM_cloud_E(:,end)*1e15) * glm_pix_size; 

%% MMIA
[int_LIS,int_GLM,int_MMIA]=integral_GLM_LIS_MMIA(GLM_cloud_E,LIS_cloud_E,MMIA_all,t1,t_end);

MMIA_cloud_E = int_MMIA;
MMIA_cloud_E(:,2:3) = MMIA_cloud_E(:,2:3)*1e-6*pi*400e3^2;  %(Van der velde et al 2020)



