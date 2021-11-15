% SET cdf MATLAB patch
if isunix
    %addpath '/usr/local/MATLAB/matlab_cdf380_patch-64'
    addpath '/Users/jaimemorandominguez/Desktop/TFG/MMIA/cdf_patch/matlab_cdf380_patch-64'
else
    addpath 'C:\matlab_cdf370_patch'
end

%%%%%%%%%%% Si se quiere probar este script por separado, indicar aquí
%%%%%%%%%%% abajo el path a los archivos .cdf en la variable 'str'

%str='/Users/jaimemorandominguez/Desktop/Pruebas_Python/MMIA_archivos/MMIA_dairy/20200610/';
%str='/Users/jaimemorandominguez/Desktop/TFG/MMIA/Entrada_MMIA/181122/descarga_2021-04-25/';
%str = '/Users/jaimemorandominguez/Desktop/Problems/MMIA_ERROR/20200712/';

%str='C:\Users\Jesús\Downloads\20200924-20210915T214936Z-001\20200924\';
%str = '/Users/jaimemorandominguez/Desktop/Final/MMIA_archivos/MMIA_dairy/20200807_7/';
%str = '/home/lrg/Desktop/20200222_0/';

tresh_frame=100; % UMBRAL DE TRESHOLD PARA CADA FOTÓMETRO

folder_name=(str);
files=dir(fullfile(folder_name,'*.cdf'));
zz=1;
n_files=length(files);

% Predefino las matrices de los fotones.
PHOT1Data_all_tmp=[];
PHOT2Data_all_tmp=[];
PHOT3Data_all_tmp=[];
t_vector_tmp=[];
t_frame_L1=[];

while zz<=n_files
    dirandfile=[str,files(zz).name];
    only_chu=0;     
    INFO=spdfcdfinfo(dirandfile);
    L=spdfcdfinfo(dirandfile);
    DATA=spdfcdfread(dirandfile);
    
    [a,numvars]=size(DATA);
    
    for nv=1:numvars
        eval([L.Variables{nv} '= double(DATA{nv});' ])
    end
    
    % Definition of max and min latitudes and longitudes

    if zz == 1
        min_lat = latitude(1);
        max_lat = latitude(1);
        min_lon = longitude(1);
        max_lon = longitude(1);
    else
        if latitude(1) < min_lat
            min_lat = latitude(1);
        end
        if latitude(1) > max_lat
            max_lat = latitude(1);
        end
        if longitude(1) < min_lon
            min_lon = longitude(1);
        end
        if longitude(1) > max_lon
            max_lon = longitude(1);
        end
    end

    % Corrected time
    t_corrected_l1=(str2num(datestr(corrected_datetime_level1(1,1),'SS.FFF'))-str2num(datestr(frame_time_phot(1,1),'SS.FFF'))) +  str2num(datestr(frame_time_phot(1,:),'SS.FFF'))
    hour = datestr(raw_datetime(1),'HH:MM:SS.FFF')(1:2)
    minute = datestr(raw_datetime(1),'HH:MM:SS.FFF')(4:5)
    


    if zz == 1
        min_t = t_corrected_l1(1);
        max_t = t_corrected_l1(end);
    else
        if t_corrected_l1(1) < min_t
            min_t = t_corrected_l1(1);
        end
        if latitude(1) > max_t
            max_t = t_corrected_l1(end);
        end
    end

    %Frame en el que hay un trigger con MXGS
    frame=1;%find(MXGSTrigger);
    
    %Extraemos las curvas de los fotometros para 1 frame
    
    format long
    % CORRECCIÓN DEL TIEMPO frame_time_phot según el que se espcifica en
    % corrected_datetime_level1.
    t_corrected_l1=(str2num(datestr(corrected_datetime_level1(1,1),'SS.FFF'))-str2num(datestr(frame_time_phot(1,1),'SS.FFF'))) +  str2num(datestr(frame_time_phot(1,:),'SS.FFF'));
    
    % Difenrencia de tiempo entre "frames"
    t_btw_corrected_frame=t_corrected_l1(2:end,1) - t_corrected_l1(1:end-1,1);
    if (t_btw_corrected_frame(:)~=0)
        
        PHOT1Data_all = reshape(PHOT1_photon_flux',[],1);
        PHOT1Data_all=PHOT1Data_all';
        PHOT2Data_all = reshape(PHOT2_photon_flux',[],1);
        PHOT2Data_all=PHOT2Data_all';
        PHOT3Data_all = reshape(PHOT3_photon_flux',[],1);
        PHOT3Data_all=PHOT3Data_all';
        
        tresh_phot1=find(PHOT1Data_all>=tresh_frame);
        tresh_phot2=find(PHOT2Data_all>=tresh_frame);
        tresh_phot3=find(PHOT3Data_all>=tresh_frame);
        
        frames_tresh=find(tresh_phot1(1,2:end)-tresh_phot1(1,1:end-1)>2);

        % MODIFICACIÓN DEL 2021-07-28
        % SI ELIMINO LOS FAKE FRAMES, EL SAMPLE RATE SE ALTERA. INSERTO 0
        PHOT1Data_all(tresh_phot3)=[];
        PHOT2Data_all(tresh_phot3)=[];
        PHOT3Data_all(tresh_phot3)=[];
%         PHOT1Data_all(tresh_phot3)=0;
%         PHOT2Data_all(tresh_phot3)=0;
%         PHOT3Data_all(tresh_phot3)=0;
        
        PHOT1Data_all_tmp=[PHOT1Data_all_tmp,PHOT1Data_all];
        PHOT2Data_all_tmp=[PHOT2Data_all_tmp,PHOT2Data_all];
        PHOT3Data_all_tmp=[PHOT3Data_all_tmp,PHOT3Data_all];
        h_asim_corr=str2num(datestr(corrected_datetime_level1(1),'HH'))*3600;
        m_asim_corr=str2num(datestr(corrected_datetime_level1(1),'MM'))*60;
        s_asim_corr=str2num(datestr(corrected_datetime_level1(1),'SS.FFF'));
        t_ini_asim_corr(zz)=h_asim_corr+m_asim_corr+s_asim_corr;
        sample_r=1e-5;
        
        trigger_length=length(PHOT1Data_all(:)); %%DEBE SER ESTE FOTÓMETRO YA QUE SE ELIMINA LOS FAKES FRAMES
        clear t_asim_corr % DEBO ELIMINAR EL VECTOR DE TIEMPO CREADO ANTERIORMENTE. 
        
        t_asim_corr(1:trigger_length,1)=0;
        t_asim_corr(1)=t_ini_asim_corr(zz);
        
        for i=2:trigger_length
            t_asim_corr(i,1)= t_asim_corr(i-1,1)+sample_r;
        end

        
        % DEBO MIRAR EL t_corrected_level 1 para el siguiente fichero ya
        % que debido al sample rate, es posible que se traslape con el
        % t_corrected del siguiente fichero. 
%         if zz>=2 && (t_asim_end_prev-t_ini_asim_corr(zz)<=0.02)  
%             disp('Traslapo')
%             t_asim_corr(1)=t_asim_end_prev; 
%             for i=2:trigger_length
%                 t_asim_corr(i,1)= t_asim_corr(i-1,1)+sample_r;
%             end
%         else
%             for i=2:trigger_length
%                 t_asim_corr(i,1)= t_asim_corr(i-1,1)+sample_r;
%             end
%             t_asim_end_prev=t_asim_corr(end,1); 
%         end
%         
%         
        t_vectorL1=[t_vector_tmp;t_asim_corr];
        t_vector_tmp=t_vectorL1;
        
        % Implementado para separa los tiempo de cada frame. 
        [frame_size]=size(PHOT1_photon_flux);
        pos_frames=1:floor(length(t_asim_corr)/frame_size(1)):length(t_asim_corr);
        t_frame=t_asim_corr(pos_frames);
        t_frame_L1=[t_frame_L1;t_frame];
%         t_frame_tmp=t_frame_L1;
        % GRÁFICO DE LOS CHU
        CHU1Data=CHU1_photon_flux;
        CHU2Data=CHU2_photon_flux;
        CHU1_pixel_longitude;
        CHU1_pixel_latitude;
        CHU2_pixel_latitude;
        CHU2_pixel_longitude;
    end
    
    
    % Fin del ciclo con solo CHU!
%         for frame=1:length(CHU1Data_exists)
            
%             if CHU1Data_exists(frame)==1 && CHU2Data_exists(frame)==1
%                 indiim=frame;
%                 
%                 %Reconstruct Images cut:
%                 dim_row=chu_maximum_row(indiim)-chu_minimum_row(indiim)+1;
%                 dim_column=chu_maximum_column(indiim)-chu_minimum_column(indiim)+1;
%                        % GRÁFICO DE LOS CHU
%                 CHU1Data=CHU1_photon_flux;
%                 CHU2Data=CHU2_photon_flux;
%                 CHU1_pixel_longitude;
%                 CHU1_pixel_latitude;
%                 CHU2_pixel_latitude;
%                 CHU2_pixel_longitude;
%                 A_CHU1=reshape(CHU1Data(indiim,:),dim_column,dim_row);
%                 A_CHU2=reshape(CHU2Data(indiim,:),dim_column,dim_row);
%                 lat_chu1=reshape(CHU1_pixel_latitude(indiim,:),dim_column,dim_row);
%                 lon_chu1=reshape(CHU1_pixel_longitude(indiim,:),dim_column,dim_row);
%                 
%                 lat_chu2=reshape(CHU2_pixel_latitude(indiim,:),dim_column,dim_row);
%                 lon_chu2=reshape(CHU2_pixel_longitude(indiim,:),dim_column,dim_row);
%                 
%                 imagenFOV1=zeros(1026,1056);
%                 imagenFOV2=zeros(1026,1056);
%                 close all
%             end
%         end
        zz=zz+1;  %% Finaliza concatenación de los fiecheros MMIA
end

if exist('t_vectorL1')
    MMIA_all(:,1)=t_vectorL1;
    MMIA_all(:,2)=PHOT3Data_all_tmp';
    save('MMIA_data','MMIA_all');
    
    space_time = [min_lat, max_lat, min_lon, max_lon, min_t, max_t]
    save('space_time', 'MMIA_space_time');

end
%clearvars;
%close all;
%clear all;
