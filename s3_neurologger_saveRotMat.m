%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script saves a prh (pitch roll heading) rotation matrix for each
% bird. The rotation matrix captures the tilt in the sensors and will be
% corrected in the script: s4_neurologger_prep4analysis.m
% 
% Code by M. Conners for Conners et al 2021:
% "Hidden Markov models identify major movement modes in accelerometer and
% magnetometer data from four albatross species." Movement Ecology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set folder directories
% fundir = set directory containing Functions_Toolboxes
% datadir = set directory containing compressed neurologger data

% Set path to required functions
% Add libraries into path.
addpath(genpath([fundir,"/Functions_Toolboxes/"]))

%% 1. Import raw sensor data (with calibrated Mag data)
cd([datadir,"L0_2_Raw_Decompressed_wCalibratedMag/"])

fileList = dir('*.txt');

df = table('Size',[length(fileList),2],'VariableTypes',{'string','double'});
df.bird = 'bird';
df.Q = [0,0,0;0,0,0;0,0,0]; %set empty Euler matrix

for i = 1:length(fileList)  
    
    m = readtable(fileList(i).name,'Delimiter',',','ReadVariableNames',true,'TreatAsEmpty',{'NA'});   
    bird=strsplit(fileList(i).name,"_");
    df(i).bird=char(strcat(bird{1},bird{2}));

    A = [Ax, Ay, Az];
    
    fs = 75; % sampling rate
    
    plott(A,fs) % will plot signal with hrs as the X axis
    grid on
    
    % identify focus region - find areas where bird is still on the water
    % or on nest at colony.
    g=ginput(2);
    g=g.*60*60*fs*24; % convert hrs to observation index
    Asub = A(g(1,1):g(2,1),:)    ;
    plot(Asub)  
    grid on
    
    % prh=[p0,r0,h0] are the Euler angles in radians of the second frame 
    % %      relative to the first
    Q = euler2rotmat(pi/180*[0 -65 0]); % change values here until you are satisfied with rotation
    V = rotate_vecs(Asub,Q);
    plott(Asub,25,V,25) % Check rotation after corrected with above rotation matrix
    grid on
    
     % if looks good, save Q
    df(i).Q = Q; 
     
    clearvars -except i fileList df
    close all
    
    
    
end
    
     
     save('Q_rotate_df.mat','df')
      
      
      
      