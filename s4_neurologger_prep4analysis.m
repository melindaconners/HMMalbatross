%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script does the final preparation of Neurologger sensor datasets for analysis. 
% Preparation is done in the following steps:
% 1. Add an index column 
% 2. Conform neurologger columns and units to AGM 
% 3. Untilt the tag-frame (tilt caused by sensors on an angle within the tag)
% 4. Decimate 75 Hz Neurologger to 25 Hz

%
% Code by M. Conners for Conners et al 2021:
% "Hidden Markov models identify major movement modes in accelerometer and
% magnetometer data from four albatross species." Movement Ecology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set folder directories---------------------------------------------------
% fundir = set directory containing Functions_Toolboxes
% datadir = set directory containing neurologger sensor data (with calibrated Mag data)
% metadir = set directory containing meta datatables required in script
% dropdir = set directory to store output

% Set path to required functions-------------------------------------------
% Add libraries into path.
addpath(genpath([fundir,"/Functions_Toolboxes/"]))


%% Loop through sensor data and prep for analysis
cd([datadir,'/L0_2_Raw_Decompressed_wCalibratedMag/'])

fileList = dir('*.txt');

for i = 1:length(fileList)  
    
    m = readtable(fileList(i).name,'Delimiter',',','ReadVariableNames',true,'TreatAsEmpty',{'NA'});   
    bird=strsplit(fileList(i).name,"_");
    bird=strcat(bird{1});

    
% 1. Add an index column (so that its stored with data index, instead of having to refer to an external metadatafile)
    % match bird
    k=find(strcmp(bird,meta.bird));

    % add index column to m
    index_v = meta.start(k):meta.xEnd(k);
    tmp = index_v';
    if length(tmp)>length(m.Ax)
        tmp=tmp(1:end-1);
    end
    m.index = tmp;


% 2. Conform neurologger columns and units to AGM 
 
    % convert P from Pascal to mbar (1 Pacscal = 0.01 mbar)
    P_mbar=m.PRPa.*0.01;
    
    newmat=array2table([m.index,ax_rotated, ay_rotated, az_rotated, mx_calibrated, my_calibrated, mz_calibrated, gx_rotated, gy_rotated, gz_rotated, m.TempC, P_mbar]); 
    newmat.Properties.VariableNames = {'index','Ax', 'Ay', 'Az', 'Mx', 'My', 'Mz', 'Gx', 'Gy', 'Gz', 'T_cel','P_mbar'};
    
    clear m P_mbar mx_calibrated my_calibrated mz_calibrated ax_rotated ay_rotated az_rotated gx_rotated gy_rotated gz_rotated tmp index_v 
    
% 3. Untilt the tag-frame (tilt caused by sensors on an angle within the tag)

    % Isolate A matrix and M matrix
    A = [newmat.Ax, newmat.Ay, newmat.Az];
    M = [newmat.Mx, newmat.My, newmat.Mz];
    G = [newmat.Gx, newmat.Gy, newmat.Gz];

    % Isolate Q for this bird
    % match bird
    birdca = {df(:).bird};
    k2=find(strcmp(bird,birdca));
    Q = df(k2).Q;

    Va = rotate_vecs(A,Q);
    Vm = rotate_vecs(M,Q);
    Vg = rotate_vecs(G,Q);
  
    newmat.Ax = Va(:,1);
    newmat.Ay = Va(:,2);
    newmat.Az = Va(:,3);
    newmat.Mx = Vm(:,1);
    newmat.My = Vm(:,2);
    newmat.Mz = Vm(:,3);
    newmat.Gx = Vg(:,1);
    newmat.Gy = Vg(:,2);
    newmat.Gz = Vg(:,3);

    clear Va Vm Vg

% 4. Decimate 75 Hz Neurologger to 25 Hz
    mm=table2array(newmat);
    m25=decdc(mm,3); %75/3 = 25
    clear mm
     
    m25(:,1)=[]; %remove index (not good after decimation)
    m25=array2table(m25);
    m25.Properties.VariableNames = {'Ax', 'Ay', 'Az', 'Mx', 'My', 'Mz', 'Gx', 'Gy', 'Gz', 'T_cel','P_mbar'};
    

% 5. Save files - these files are clean and ready for analysis.
    
    % full resolution 75Hz
    filename0=[bird,'_hrl_ready_75hz.txt'];
    writetable(newmat,[dropdir,'/L1_1_Analysis_Ready/75hz/',filename0])
    
    % 25Hz file (without index - save in filename)
    start_x = newmat.index(1);
    filenamei = [bird,'_',num2str(start_x),'_hrl_ready_25hz.txt'];    
    writetable(m25,[dropdir,'/L1_1_Analysis_Ready/25hz/',filenamei])
    

clearvars -except i fileList df meta

end