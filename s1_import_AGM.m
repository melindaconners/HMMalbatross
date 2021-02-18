
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script imports AGM .csv files and aligns sensor frames to each other and to the bird frame.
% Time is checked so that any irregularity in Hz results in expansion of rows as NaNs so that all observations
% are equivalent to 25 Hz. It also adds 3 (empty) columns that are footholders for Gyro data so that
% the files are compatible with s2_calibrate_magnetometer.m script. 
% Pressure and temperature data (originally sampled at 1Hz) are
% interpolated to 25 Hz to match resolution of acc and mag data. 
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
% dropdir = set directory to store output files.

fundir = '~/Dropbox/Academia/SUNY/Project_Components/1_NSF_Albatross_and_Wind/Manuscripts/ms1_suny_behavClass/ME2021_Connersetal_CodeUpload/' ;
datadir = '~/Dropbox/Academia/SUNY/Project_Components/1_NSF_Albatross_and_Wind/Manuscripts/ms1_suny_behavClass/ME2021_Connersetal_CodeUpload/';
dropdir = '~/Dropbox/Academia/SUNY/Project_Components/1_NSF_Albatross_and_Wind/Manuscripts/ms1_suny_behavClass/ME2021_Connersetal_CodeUpload/testout/';


% Set path to required functions-------------------------------------------
% Add libraries ('d3matlab','tagmatlab','x3toolbox','tagtools') into path.
addpath(genpath(strcat(fundir,"Functions_Toolboxes/")))

% Set working directory to Neurologger Sensor Data--------------------------------
% Raw Compressed dat
CurrentPath =strcat(datadir, "testfile/") ;

cd(CurrentPath) 

fileList = dir('*.csv');

for i= 1:length(fileList)

m = readtable(fileList(i).name,'Delimiter',',','ReadVariableNames',true,'TreatAsEmpty',{'NA'}); % read agm bird i's file
% Columns of these data: birdid,datetime_ms_gmt,Ax,Ay,Az,Mx,My,Mz,Pmbar_air,TempC,Activity

bird=strsplit(fileList(i).name,"_");
bird=strcat(bird{1},'_',bird{2});


% Rotate Acc and Mag sensor frames to align with each other and the bird frame
Acc=[m.Ax, m.Ay, m.Az];% Acc data in ENU (East-North-Up)
Mag=[m.My,m.Mx,-m.Mz]; % Mag data converted from NED (North East Down) to ENU (North East Up)

m.Ax = Acc(:,1);
m.Ay = Acc(:,2);
m.Az = Acc(:,3);

m.Mx = Mag(:,1);
m.My = Mag(:,2);
m.Mz = Mag(:,3);

clear Mag Acc


% TIME CHECK ----------------------------------------------------------
% Check to find instances where interval is not 25 Hz (e.g. in AGM data
% there are breaks, with some observations being a full second apart etc)
% Where breaks exist need to expand matrix by adding NaN rows to keep 
t = datetime(m.datetime_ms_gmt,'InputFormat','dd/MM/yyyy HH:mm:ss.SSS');
out = milliseconds(diff(t));
break_ix = find(out ~= 40, 1);

    %out should be 40 (1000/fs) so check if instances ~=40 (not equal to) and correct 
    if ~isempty(break_ix)
        df=table(break_ix,out(break_ix));
        df.Properties.VariableNames = {'index', 'value'};   
        format long
        writetable(df,strcat(dropdir,bird,'_timefailcheck.csv')) %record breaks and write table
        
         % where breaks exist, need to insert NaN rows --------------------
         
         % Loop through each break and 
            % 1. pre-emptively expand datatable to new size (fill with NaN and 
            % 2. fill out a dataframe with break durations that will be used for indexing
            
            break_df = cell2table(cell(0,4), 'VariableNames', {'breaknum', 'breaks_ix', 'durx', 'durcum'});
         
             for j=1:length(break_ix)

                nsec=out(break_ix(j))/1000; % number of seconds of break
                durx=round(nsec*24 + (nsec-1)); % number of NA rows to represent break
    
                break_df.breaknum(j) = j; % numeric number of break
                break_df.breaks_ix(j) = break_ix(j); % index of break from m
                break_df.durx(j) = durx; % duration of current break
                
                 % create NaN rows 
                 if j==1
                     expandlength=durx;
                 else
                     expandlength=expandlength+durx;
                 end

                break_df.durcum(j) = expandlength;   % cumulative duration of breaks up to and including current break
 
             end
         
            % create new datatable
            newm = NaN(height(m)+expandlength,width(m));
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % fill newm with old m info in chunks
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
         % ----------------------------------------------------------------
         if length(break_ix)==1 % simple indexing if only one break
             % how long (in seconds) is break?
             nsec=out(break_ix)/1000;
             durx=round(nsec*24 + (nsec-1));
             
             newm(1:break_ix,:)=m(1:break_ix,:);
             newm(break_ix+durx:height(newm),:) = m(break_ix+durx:height(newm),:);
         end
         
         clear durx nsec
         
         % ----------------------------------------------------------------
         if length(break_ix)>1 % more complicated indexing for multiple breaks
             
             for j=1:length(break_ix) % loop through each break
                 
                 if j==1 % for first break:
                     newm(1:break_ix,:)=m(1:break_ix,:);
                        
                 elseif j==length(break_ix) % for last break:
                     newm(break_df.break_ix(j-1)+1+break_df.durcum(j-1):height(newm),:) = m(break_df.break_ix(j-1)+1:height(m),:);
                 
                 else % for middle breaks
                     newm(break_df.break_ix(j-1)+1+break_df.durcum(j-1):break_df.break_ix(j)+breaks_df.durcum(j-1),:) = m(breaks_df.break_ix(j-1)+1:breaks_df.break_ix(j),:);
                     
                 end
   
             end   
         end
 
    end % ends loop for 'if there are breaks in timeseries regularity'
% END TIME CHECK ----------------------------------------------------------


% interpolate T and P------------------------------------------------------

m.Pmbar_air=fillmissing(m.Pmbar_air,'linear');
m.TempC=fillmissing(m.TempC,'linear');

% Create empty Gyro columns (so conform to match Neurologger data)
m.Gx = cell(height(m),1);
m.Gy = cell(height(m),1);
m.Gz = cell(height(m),1);

% rearrange matrix
m = m(:,[1:8 12:14 9:10]); % order: birdid, datetime, Ax, Ay, Az, Mx, My, Mz, Gx, Gy, Gz, T, P

% save file
filenamei = [bird,'_agm_25hz_imported.txt'];    
writetable(m,strcat(dropdir,filenamei))

clearvars -except i fileList fundir datadir dropdir CurrentPath

end