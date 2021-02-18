
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script imports the compressed Neurologger .dat files and converts
% sensor and ecg data to their appropriate units and scale.
% Sensor frames are aligned to each other and to the bird frame.
% 
% Code by Alexei Vyssotski and adaped by M. Conners for Conners et al 2021:
% "Hidden Markov models identify major movement modes in accelerometer and
% magnetometer data from four albatross species." Movement Ecology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set Environment

% Set folder directories
% fundir = set directory containing Functions_Toolboxes
% datadir = set directory containing compressed neurologger data

% Set path to required functions-------------------------------------------
% Add libraries ('d3matlab','tagmatlab','x3toolbox','tagtools') into path.
addpath(genpath([fundir,"/Functions_Toolboxes/"]))

% Set working directory to Neurologger Sensor Data--------------------------------
% Raw Compressed dat
CurrentPath =[datadir, "/L0_0_Raw_Compressed_dat/"] ;

cd(CurrentPath) 

fileList = dir('*.dat');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Routine needed to clear ghost files if working in Seagate Hard Drive
% Ignore except in very specific cases where ghost files are created and
% then read by the function dir.
% skipx=[];
% for k = 1:length(fileList)
% skipx(k)=~startsWith(fileList(k).name,'._','IgnoreCase',true);
% end
% 
% fileList = fileList(find(skipx));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%% Unpack Raw Neurologger .dat files

for i= 1:length(fileList)
    
FileName = fileList(i).name(1:length(fileList(i).name)-4);
FileOutputName = [pwd filesep 'OutputFolder' filesep FileName];
bird = fileList(i).name(1:8);
    
% set parameters-----------------------------------------------------------
Parameters.NChannels = 1; %number of electrophysiological channels
Parameters.sps = 600;     % Hz 
Parameters.Ku = 166.6;    %input range +/-6 mV
Parameters.Coef = 2000/Parameters.Ku/1024;    %2V, 10 bit probably +/-0.5mV range


% Import data--------------------------------------------------------------
[Calibration, EEG_t, Accel_t, Pressure_t, Temp_t, Gyro_t, Magn_t, Sync_t,HeaderPosition] = Datalogger4ChConverter_MARG_Altimeter([FileName '.dat'], Parameters);


%Write header positions in the text file-----------------------------------
fileID = fopen([FileOutputName '_HeaderPos.txt'],'w');
for I_t = 1:length(HeaderPosition)
  fprintf(fileID,'%d\n',HeaderPosition(I_t));  
end
fclose(fileID);

%Get tempterature and pressure---------------------------------------------
[ T, P ] = GetCalibTemperaturePressureBMP388( Temp_t, Pressure_t, Calibration );  %Celsius, Pascals

% Extract sensor data------------------------------------------------------

ECG_t = ECG_t(:,1);
ECG = single(cast(ECG_t,'uint16'))*Parameters.Coef;
ECG = ECG - repmat(median(ECG),size(ECG,1),1);

Accel = zeros(size(Accel_t),'int16');
for I_t = 1:3
  Accel(:,I_t) = typecast(Accel_t(:,I_t),'int16');  
end    
Accel = single(Accel)/256/128*16; % +/-16g

Gyro = zeros(size(Gyro_t),'int16');
for I_t = 1:3
  Gyro(:,I_t) = typecast(Gyro_t(:,I_t),'int16');  
end    
Gyro = single(Gyro)/256/128*2000*pi/180; %2000 degrees/s -> rad/s

Magn = zeros(size(Magn_t),'int16');
for I_t = 1:3
  Magn(:,I_t) = typecast(Magn_t(:,I_t),'int16');  
end    
Magn = single(Magn)/256/128*49.152; %49.152 Gauss

%% rotate acc+gyr axes to align with mag axes (mag frame is aligned with bird frame in 2019-2020 tag configuration)
% to see schematics of sensor orientations, see neurologger_sensor_frame_log.pdf 
% Rotations will be unique to tag configurations and orientation of tag on
% animal.
    ax_rotated = Accel(:,2).*-1;
    ay_rotated = Accel(:,1).*-1;
    az_rotated = Accel(:,3);
    
    gx_rotated = Gyro(:,2).*-1;
    gy_rotated = Gyro(:,1).*-1;
    gz_rotated = Gyro(:,3);
    
    Accel_r = [ax_rotated, ay_rotated, az_rotated];
    Gyro_r  = [gx_rotated, gy_rotated, gz_rotated];

%% save sensor files and EEG file separately

SensorMat = [Accel_r, Magn, Gyro_r, T, P]; % 75 Hz
ECG = ECG; % 600 Hz

sensorfilename = [bird,'_sensormatraw75hz.mat'];
ecgfilename = [bird,'_ecgmatraw600hz.mat'];
save(sensorfilename, 'SensorMat')
save(eegfilename,'ECG')

clearvars -except OutputFolder fileList i CurrentPath
end

%%

