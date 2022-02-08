%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script sweeps across all directories with .cdf MMIA files and 
% extracts their data for photometers (time and signal vectors) as well as
% ISS position and initial and final corrected times.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SET cdf MATLAB patch
if isunix
    %addpath '/usr/local/MATLAB/matlab_cdf380_patch-64'
    addpath '/Users/jaimemorandominguez/Desktop/TFG/MMIA/cdf_patch/matlab_cdf380_patch-64'
else
    addpath 'C:\matlab_cdf370_patch'
end

% Get path to the directory with all the directories with .cdf files
dirs_path = fileread('mmia_dirs_path.txt');
dirs_path = dirs_path(1:end-1);  % Delete the new line character

% Get a list of all files and folders in this folder
files = dir(dirs_path);
% Get a logical vector that tells which is a directory
dirFlags = [files.isdir];
% Extract only those that are directories
directories = files(dirFlags); % A structure with extra info
% Get only the folder names into a cell array
dir_names = {directories(3:end).name}; % Start at 3 to skip . and ..

mats_path = fileread('mmia_mats_path.txt');
mats_path = mats_path(1:end-1);

for i = 1:length(dir_names)
    current_event = cell2mat(dir_names(i));
    str = append(dirs_path, current_event, '/');
    sentence = append(int2str(i),'. Extracting MMIA data from .cdf files for ',current_event,'...\n');
    fprintf(sentence)

% (still inside for i = 1:length(dir_names) loop)
tresh_frame=100; % UMBRAL DE TRESHOLD PARA CADA FOTOMETRO

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

    try
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

    frame=1;%find(MXGSTrigger);

    %Extraemos las curvas de los fotometros para 1 frame

    format long
    % CORRECCION DEL TIEMPO frame_time_phot segun el que se espcifica en
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

        % MODIFICACION DEL 2021-07-28
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

        trigger_length=length(PHOT1Data_all(:)); %%DEBE SER ESTE FOT�METRO YA QUE SE ELIMINA LOS FAKES FRAMES
        clear t_asim_corr % DEBO ELIMINAR EL VECTOR DE TIEMPO CREADO ANTERIORMENTE.

        t_asim_corr(1:trigger_length,1)=0;
        t_asim_corr(1)=t_ini_asim_corr(zz);

        for i=2:trigger_length
            t_asim_corr(i,1)= t_asim_corr(i-1,1)+sample_r;
        end

        t_vectorL1=[t_vector_tmp;t_asim_corr];
        t_vector_tmp=t_vectorL1;

        % Implementado para separa los tiempo de cada frame.
        [frame_size]=size(PHOT1_photon_flux);
        pos_frames=1:floor(length(t_asim_corr)/frame_size(1)):length(t_asim_corr);
        t_frame=t_asim_corr(pos_frames);
        t_frame_L1=[t_frame_L1;t_frame];

        CHU1Data=CHU1_photon_flux;
        CHU2Data=CHU2_photon_flux;
        CHU1_pixel_longitude;
        CHU1_pixel_latitude;
        CHU2_pixel_latitude;
        CHU2_pixel_longitude;
    end
    catch
        warning_message = append('Problem with the file ', dirandfile);
        warning(warning_message);
    end
        zz=zz+1;  %% Finaliza concatenacion de los fiecheros MMIA
end

if exist('t_vectorL1')
    MMIA_all(:,1)=t_vectorL1;
    MMIA_all(:,2)=PHOT1Data_all_tmp';
    MMIA_all(:,3)=PHOT2Data_all_tmp';
    MMIA_all(:,4)=PHOT3Data_all_tmp';
    mmia_data_save_name = append(mats_path,current_event, '_data');
    % File name // Variable name
    save(mmia_data_save_name,'MMIA_all');

    space_time = [min_lat, max_lat, min_lon, max_lon, min(t_vectorL1), max(t_vectorL1)];
    mmia_info_save_name = append(mats_path,current_event, '_info');
    save(mmia_info_save_name, 'space_time');
end
clearvars -except dirs_path dir_names mats_path
end

clearvars;
close all;
clear all;
