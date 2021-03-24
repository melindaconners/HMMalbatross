%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prep all Technosmart ACC tags for analysis:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% NOTE: This code ignores the mag data since we are not using it from
%%%% the 2019-2020 AGM tags. If need to prep AGM data for analysis, then
%%%% this requires another intermediate step of calibrating the
%%%% magnetometer data. Summary: MAG DATA HERE NOT CALIBRATED

%%%% Steps taken to prep:

% 1. Check for upside-down tag placements - for Technosmart tags, there
% should be none, but including this to have a conservative check. 

% 2. Rotate M frame to match A frame and bird frame (this is sort of a filler
% step since the M data is not calibrated and we are not using it, but don't want to get out of the
% habit of doing this. 

% 3. Check time intervals and expand with NAs if necessary

% 4. Interpolate P and T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set folder directories
datadir = '/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Bird_Island/2019_2020/Tag_Data/L1/2018-2019_ACC_Technosmart_Accel_Uniformat/'; % Bird Island 2019-2020
% datadir = '/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Midway/L1_Acc_Appended/'; % Midway 2018-2019

dropdir = '/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Analysis/HMM/data/s2_Acc_AnalysisReady/';

% Set path to required functions-------------------------------------------
% Add libraries ('d3matlab','tagmatlab','x3toolbox','tagtools') into path.
addpath(genpath("~/Dropbox/Academia/SUNY/Project_Components/BIRD_ISLAND/Analyses/Functions_Toolboxes/"))


% Set working directory to Neurologger Sensor Data--------------------------------
% Raw Compressed dat
CurrentPath =datadir;

cd(CurrentPath) 

fileList = dir('*.txt');

fs = 25; %sampling rate
fs_ms = milliseconds(seconds(1/fs));
fs_f = fs-1; % for indexing during mat expansino if breaks ixist in timeseries

% Pre-allocate loop vars to store table of z axis orientation issues --------------------------------
nfiles = length(fileList); % calculate length of this struct of files
T = table(cell(nfiles,1),cell(nfiles,1),zeros(nfiles,1),cell(nfiles,1),zeros(nfiles,1),zeros(nfiles,1),zeros(nfiles,1),'VariableNames', {'BirdID','TagType','Z-Mean-G','Z-Orientation','Dur-Days','Num_Breaks','Max_Break_Sec'});                     


for i= 1:length(fileList)

m = readtable(fileList(i).name,'Delimiter',',','ReadVariableNames',true,'Format','auto','TreatAsEmpty',{'NA'}); % read bird i's file
% Columns of these data: TagID,DateTime,Ax,Ay,Az,Mx,My,Mz,Pressure,Temperature,Activity

bird=strsplit(fileList(i).name,"_");
birdid=strcat(bird{1},'_',bird{2});
tagtype = bird{4};
T.("BirdID")(i,:) = {birdid};
T.("TagType")(i,:) = {tagtype};



% -------------------------------------------------------------------------
% 1. Pre-check to make sure tag was put on correctly 
% For all birds - do automated check for upside-down tag and automatically correct using euler rotation matrix if necessary:
% -------------------------------------------------------------------------

z_mean = mean(m.Az, 'omitnan')  ;
T.("Z-Mean-G")(i,:) = z_mean;

     if z_mean>0
         T.("Z-Orientation")(i,:) = {'Correct'};
         Q = euler2rotmat(pi/180*[0 0 0]); % no flip necessary
     else
         T.("Z-Orientation")(i,:) = {'Upside Down'};
         Q = euler2rotmat(pi/180*[0 180 180]); % Forward-Right-Up -> Flip upside down to right side up

     end
     
     % Use Q to correct tag frame 
     % Rotate Accelerometer Sensor Frame
     A = [m.Ax, m.Ay, m.Az];
     Av = rotate_vecs(A,Q); % Use Q to rotate A
     m.Ax = Av(:,1); % replace rotated Acc data
     m.Ay = Av(:,2); % replace rotated Acc data
     m.Az = Av(:,3); % replace rotated Acc data
     clear A Av
     
     % Rotate Magnetometer Sensor Frame
     M = [m.Mx, m.My, m.Mz];
     Mv = rotate_vecs(M,Q); % Use Q to rotate M
     m.Mx = Mv(:,1); % replace rotated Acc data
     m.My = Mv(:,2); % replace rotated Acc data
     m.Mz = Mv(:,3); % replace rotated Acc data
     clear M Mv
        
    
% -------------------------------------------------------------------------
% 2. SENSOR FRAME ROTATION
% Rotate Acc and Mag sensor frames to align with each other and the bird frame:
% -------------------------------------------------------------------------
    
     
    Acc=[m.Ax, m.Ay, m.Az];% Acc data in ENU (East-North-Up)
    Mag=[m.My,m.Mx,-m.Mz]; % Mag data converted from NED (North East Down) to ENU (North East Up)

    m.Ax = Acc(:,1);
    m.Ay = Acc(:,2);
    m.Az = Acc(:,3);

    m.Mx = Mag(:,1);
    m.My = Mag(:,2);
    m.Mz = Mag(:,3);

    clear Mag Acc 
     
% -------------------------------------------------------------------------
% 3. CHECK INTERVALS
% Check to find instances where interval is not 25 Hz (e.g. in AGM data
% there are breaks, with some observations being a full second apart etc)
% Where breaks exist need to expand matrix by adding NaN rows to keep
% intervals regular. 
% -------------------------------------------------------------------------

% CONFIRM DATEDATE FORMAT:
% For some reason a handful of files being read as YYY-dd-MM instead of
% YYY-MM-dd.
% To find, check if the number of unique months is greater than the
% number of unique days. if so, need to convert to YYYY-MM-dd.
if length(unique(day(m.DateTime))) < length(unique(month(m.DateTime))) 
   % First convert to string:
   dts=string(m.DateTime);
   % Then convert to datetime specifying format
   dtfix = datetime(dts,'InputFormat','yyyy-dd-MM HH:mm:ss.SSS'); % say what the format is when converting to datetime
   dtfix.Format = 'yyyy-MM-dd HH:mm:ss.SSS'; % write to make format same as other files
   m.DateTime = dtfix;
end


% Quick meta
daysdur = days(m.DateTime(end)-m.DateTime(1));
T.("Dur-Days")(i,:) = daysdur;

% First: Identify irregular intervals:

         out = milliseconds(diff(m.DateTime));
%          break_ix = find(out ~= 40, 1); % AGMs are consistently @40ms but
%          some AxyAirs are between 30-50 ms. This is "good enough" but
%          look for breaks that are > 80 ms and correct matrix by expanding
%          with NaNs.
         break_ix = find(out>80);
         

        % Loop through break indices-------------------------------
        if ~isempty(break_ix)
            
                 if height(m)-break_ix<5 %ignore bc break is at end of file.

                     
                        full_exp_mat=m(1:break_ix,:);
                        T.("Num_Breaks")(i,:) = length(break_ix);
                        T.("Max_Break_Sec")(i,:) =  999; %code for right at end

                 else

                        T.("Num_Breaks")(i,:) = length(break_ix);
                        T.("Max_Break_Sec")(i,:) =  seconds(milliseconds(max(out(break_ix),[],'omitnan')));

                        % -----------------------------------------------------
                        % If there is only one break:
                        % -----------------------------------------------------
                        if length(break_ix)==1 % if there is only one break:
                            % calculate number of rows to represent break
                            nsec=out(break_ix)/1000; % number of seconds of break

                            % RELATE BREAK in SEC to 25Hz for NUMBER OF ROWS - 
                            %For 1 second remember should be fs-1 for each second
                            if nsec==1 % special case for 1 second
                                durx = fs_f;
                            else
                                middles = fs_f*nsec;
                                sec_lines = nsec-1;
                                xx=mod(nsec,1); % how to check if even or fraction

                                if xx==0 %if even
                                    durx = middles + sec_lines;
                                else % if fraction
                                    fract_lines = milliseconds(seconds(xx))/fs_ms;
                                    durx = middles + sec_lines + fract_lines;
                                end
                            end


                            %append chunk with NAs
                            mchunk = m(1:break_ix(1),:);
                            nachunk = array2table(NaN(durx,width(mchunk)));
                            nachunk.Properties.VariableNames = mchunk.Properties.VariableNames;
                            nachunk.TagID=cellstr(repmat(birdid, height(nachunk), 1));
                            nachunk.Activity=cellstr(repmat("", height(nachunk), 1));

                            % create fake time vec
                            t1 = mchunk.DateTime(end) + milliseconds(40);
                            t2 = t1 + milliseconds(40)*(durx-1);
                            tfill = t1:milliseconds(40):t2;
                            ttfill = tfill';
                            nachunk.DateTime=ttfill;

                            % EXPAND MAT
                            exp_chunk = [mchunk;nachunk]; % Append NA chunk
                            full_exp_mat = [exp_chunk;m(break_ix(1)+1:end,:)]; %

                        else

                            % -----------------------------------------------------
                            % If there are multiple breaks, loop through each
                            % break.
                            % -----------------------------------------------------

                            for j = 1:length(break_ix) % start loop through break j

                                nsec=out(break_ix(j))/1000; % number of seconds of break

                                % isolate chunk of m before break:
                                if j==1 % for first break:
                                    mchunk = m(1:break_ix(j),:);
                                else    % for all other breaks
                                    mchunk = m(break_ix(j-1)+1:break_ix(j),:);
                                end


                                % Get durx - RELATE BREAK in SEC to 25Hz for NUMBER OF ROWS -

                                    if nsec==1 % special case for 1 second
                                        durx = fs_f;
                                    else
                                        middles = fs_f*nsec;
                                        sec_lines = nsec-1;
                                        xx=mod(nsec,1); % how to check if even or fraction

                                        if xx==0 %if even
                                            durx = middles + sec_lines;
                                        else % if fraction
                                            fract_lines = milliseconds(seconds(xx))/fs_ms;
                                            durx = middles + sec_lines + fract_lines;
                                        end
                                    end


                                % create a chunk of NAs to represent the break
                                nachunk = array2table(NaN(durx,width(mchunk)));
                                nachunk.Properties.VariableNames = mchunk.Properties.VariableNames;
                                nachunk.TagID=cellstr(repmat(birdid, height(nachunk), 1));
                                nachunk.Activity=cellstr(repmat("", height(nachunk), 1));

                                % create fake time vec
                                t1 = mchunk.DateTime(end) + milliseconds(40);
                                t2 = t1 + milliseconds(40)*(durx-1);
                                tfill = t1:milliseconds(40):t2;
                                ttfill = tfill';
                                nachunk.DateTime=ttfill;

                                % r bind mchunk_j with the NAs
                                exp_chunk = [mchunk;nachunk];

                                if j==1
                                    newmat = exp_chunk;
                                else
                                    newmat = [newmat;exp_chunk];
                                end


                            end % end loop through break j

                            % Add last chunk
                            last_chunk = m(break_ix(length(break_ix))+1:end,:);
                            full_exp_mat = [newmat;last_chunk];
                        end
                 end % end to check if loop is right at end
         

            
        else
            
            full_exp_mat = m;
            
        end
			


% -------------------------------------------------------------------------
% 4. BASIC LINEAR INTERPOLATION OF TEMP AND PRESSURE
% -------------------------------------------------------------------------

full_exp_mat.Pressure=fillmissing(full_exp_mat.Pressure,'linear');
full_exp_mat.Temperature=fillmissing(full_exp_mat.Temperature,'linear');

% -------------------------------------------------------------------------
% 5. WRITE ANALYSIS READY FILE
% -------------------------------------------------------------------------
filenamei = [dropdir,birdid,'_',tagtype,'_25hz_AnalysisReady.txt'];    
writetable(full_exp_mat,filenamei)

clearvars -except i fs fs_ms fs_f fileList CurrentPath dropdir T 
    
end
    


% A in G
% M in microTesla
% Pressure in mbar
% Temperature in Celsius



writetable(T,[dropdir,'Table_T_acc_s2_bi1920.csv'])







