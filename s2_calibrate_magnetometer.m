%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calibrates the magnetometer data using a data-driven method.
% Calibration diagnostics determine whether the files will be calibrated as
% a full file or piecewise in smaller segments.
% 
% Code by A. Vyssotski, E. Heywood, and M. Conners for Conners et al 2021:
% "Hidden Markov models identify major movement modes in accelerometer and
% magnetometer data from four albatross species." Movement Ecology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set folder directories
% fundir = set directory containing Functions_Toolboxes
% datadir = set directory containing txt files of sensor data

% Set path to required functions-------------------------------------------
% Add libraries ('d3matlab','tagmatlab','x3toolbox','tagtools') into path.
addpath(genpath([fundir,"/Functions_Toolboxes/"]))

% Set working directory to Neurologger Sensor Data-------------------------
CurrentPath =[datadir,'/L0_1_Raw_Decompressed_txt/1_SensorData/'];
cd(CurrentPath) 

fileList = dir('*.txt'); %

% Make table to store all indices and diagnostics of the calibration
% process

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import SensorMat from Bird i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except i CurrentPath fileList sps timesplit
i = i+1; % Calibrate files one at a time, rather than in a loop.
    
SensorData = readtable(fileList(i).name); % this will load 'SensorData' 
% Neurologger sensormats: 11 column dataframes: Ax  Ay  Az Mx  My  Mz  Gx  Gy  Gz T P
% AGM sensormats: 13 column dataframes: bird datetime Ax  Ay  Az Mx  My  Mz  Gx  Gy  Gz T P
% Adjust code accordingly

% Define the bird ID for later file naming and input into diagnostic tables
bird=strsplit(fileList(i).name,"_");
bird=strcat(bird{1},bird{2});

% Isolate triaxial mag data
Mag = SensorData(:,4:6); 
Mag = table2array(Mag); % Convert the mag data table to a matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Calibration Procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define sampling frequency and the timesplit for piecewise calibration
% sps = 25;   %25 Hz for AGMs
sps = 75;   %75 Hz for Neurologgers
timesplit = 21600; % 21600 seconds = 6 hour chunks; 14400=4 hr chunks

% ONLY EXECUTE THE FOLLOWING FOR VERY SPIKY AGM DATA
% Remove the huge spikes in AGM data
%ixy=find(Mag(:,1)>10 | Mag(:,1) <-10);
%ixx=find(Mag(:,2)>10 | Mag(:,2) <-10);
%ixz=find(Mag(:,3)>10 | Mag(:,3) <-10);
%Mag(ixy,1)=NaN;
%Mag(ixx,2)=NaN;
%Mag(ixz,3)=NaN;
%Mag(:,1)=fillmissing(Mag(:,1),'linear');
%Mag(:,2)=fillmissing(Mag(:,2),'linear');
%Mag(:,3)=fillmissing(Mag(:,3),'linear');

% 1. Cut off "bookends" to get rid of the messy start/end of the file (critical to save these indices)
plot(Mag(:,1:3))
[book_x, book_y] = ginput(2) % Select the main part of the track excluding the very beginning and end where there is clear interference

book_indices = book_x;

% Subset the sensor data table and the mag data by these bookend indices 
Mag_s = Mag(book_indices(1):book_indices(2),:); 
SensorData_s = table2array(SensorData);
SensorData_s = SensorData_s(book_indices(1):book_indices(2),:);

%% CALIBRATE PIECEWISE
chsize = timesplit * sps; % chunk size in sampling units
flen = length(Mag_s); % file length 

starts = 1:chsize:flen; % vector of start indices
ends = starts + chsize - 1; % vector of end indices

% make last end the ending index of the file
ends(end) = flen;

nch = length(starts); % number of chunks to iterate over

% CREATE TABLE TO STORE THIS BIRDS DIAGNOSTICS
%diags = table('VariableNames',{'BirdID','InitStartIndex', 'InitEndIndex', 'ChunkStartIndex','ChunkEndIndex','EllipsoidCenter','EllipsoidEvecs', 'AlgebraicForm', 'SDsphere', 'ResidualM2'}, 'VariableTypes', {'string', 'double', 'double', 'double','double','string', 'string', 'string', 'double', 'double'});
diags = table(cell(0,1),zeros(0,1),zeros(0,1), zeros(0,1), zeros(0,1), cell(0,1), cell(0,1), cell(0,1),cell(0,1), zeros(0,1), zeros(0,1), 'VariableNames',{'BirdID','InitStartIndex', 'InitEndIndex', 'ChunkStartIndex','ChunkEndIndex','EllipsoidCenter','EllipsoidRadii','EllipsoidEvecs', 'AlgebraicForm', 'SDsphere', 'ResidualM2'});

% FILL the table with appropriate BIRD ID, Master Indices (bookend) and
% chunk indices
diags.BirdID(1:nch,1) = {bird};
diags.InitStartIndex(1:nch) = round(book_indices(1));
diags.InitEndIndex(1:nch) = round(book_indices(2));
diags.ChunkStartIndex(:) = starts;
diags.ChunkEndIndex(:) = ends;


%% Loop through each chunk: calibrate and record diagnostics
    for q = 1:nch
        
        fnum = num2str(q);
        curM = Mag_s(starts(q):ends(q), :);  % grab q chunk 
        
        % Apply Median Filter to this chunk of mag data 'curM'
        clear Mf
        Mf(:,1)=median_filter(curM(:,1), 10); Mf(:,2)=median_filter(curM(:,2), 10); Mf(:,3)=median_filter(curM(:,3), 10);
        %plott(curM, 75, Mf,75) % check with plot - plot of unfiltered and filtered.
        
        % Isolate X, Y, and Z
        x = double(Mf(:,1));   
        y = double(Mf(:,2));   
        z = double(Mf(:,3));
        % -------------------------------------------------------------------------
        % Fit the Ellipsoid 
        % -------------------------------------------------------------------------
        [ center, radii, evecs, v, chi2, evals, stdsum ] = ellipsoid_fit([ x y z ], '' );  
        % CHECK SOME THINGS
        sqrt( stdsum / size( x, 1 ) )
        % Print some initial descriptors
        %fname=[bird,'_calibration_results_',fnum,'.txt'];
        %fileID = fopen(fname,'w');
        % -------------------------------------------------------------------------
        % Correct the Data with the Ellipsoid Fit Factors
        % -------------------------------------------------------------------------
        
        % Step 1: subtract center
        d = [ x - center(1), y - center(2), z - center(3) ]; % shift data to origin
        
        % Step 2: rotate to match axes of the ellipsoid to coordinate axes
        d = d * evecs; % rotate to cardinal axes of the conic;
        
        % Step 3: scale data to get radius 1 (one can scale later to local magnetic
        % field strength if needed)
        d = [ d(:,1) / radii(1), d(:,2) / radii(2), d(:,3) / radii(3) ]; % normalize to the conic radii
        
        % Step 4: rotate back to have zero rotation at the end - only scaling should be present!
        d = d * inv(evecs);

        xc = d(:,1); yc = d(:,2); zc = d(:,3);
        v = [1, 1, 1, 0, 0, 0, 0, 0, 0, -1];
        % %This is slightly different method to compute residuals used in Matlab Help to magcal()    
        N = length(xc); expMFS = 1;
        r = sum(d.^2,2) - expMFS.^2;
        E = sqrt(r.'*r./N)./(2*expMFS.^2);
        
        % IF RESIDUAL ERROR 'E' below 0.1 move on, if not, tweak until it is
        % Very Important!! Need to save indices noting where track ends were cut
        % off in order to sync times. The neurologgers do not come with a time
        % column, so we'll need the indices.
        
        % Save calibration details
        diags.EllipsoidCenter(q) = {center};
        diags.EllipsoidRadii(q) = {radii};
        c = [evecs(1), evecs(2), evecs(3), evecs(4), evecs(5), evecs(6), evecs(7), evecs(8), evecs(9)];
        diags.EllipsoidEvecs(q) = {c};
   
        diags.AlgebraicForm(q) = {v};
        diags.SDsphere(q) = sqrt( chi2 / size( x, 1 ) );
        
        % Visualize the uncorrected ellipsoid ------------------------------------
        St1 = 1;
        Step = round(1e3/St1);   %1e3; %every 1000th point will be plotted
        figure (1); %clf;
        subplot(1,2,1)
        plot3( x(1:Step:end), y(1:Step:end), z(1:Step:end), '.r' );
        hold on;

        mind = min( [ x y z ] );
        maxd = max( [ x y z ] );
        nsteps = 50;
        step = ( maxd - mind ) / nsteps;
        [ xt, yt, zt ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
        
        Ellipsoid = v(1) *xt.*xt +   v(2) * yt.*yt + v(3) * zt.*zt + ...
            2*v(4) *xt.*yt + 2*v(5)*xt.*zt + 2*v(6) * yt.*zt + ...
            2*v(7) *xt    + 2*v(8)*yt    + 2*v(9) * zt;
        p = patch( isosurface( xt, yt, zt, Ellipsoid, -v(10) ) );
        hold off;
        set( p, 'FaceColor', 'g', 'EdgeColor', 'none' );
        view( -70, 40 );
        axis vis3d equal;
        camlight;
        lighting phong;
        title('Original data');
        
        % PRINT Residual Error to diagnostic text file
        diags.ResidualM2(q) = E;
        
        %%Draw corrected data
        subplot(1,2,2)
        plot3( xc(1:Step:end), yc(1:Step:end), zc(1:Step:end), '.r' );
        hold on;
        
        %draw fit ----------------------------------------------------------------
        mind = min( [ xc yc zc ] );
        maxd = max( [ xc yc zc ] );
        nsteps = 50;
        step = ( maxd - mind ) / nsteps;
        [ xt, yt, zt ] = meshgrid( linspace( mind(1) - step(1), maxd(1) + step(1), nsteps ), linspace( mind(2) - step(2), maxd(2) + step(2), nsteps ), linspace( mind(3) - step(3), maxd(3) + step(3), nsteps ) );
        
        Ellipsoid = v(1) *xt.*xt +   v(2) * yt.*yt + v(3) * zt.*zt + ...
            2*v(4) *xt.*yt + 2*v(5)*xt.*zt + 2*v(6) * yt.*zt + ...
            2*v(7) *xt    + 2*v(8)*yt    + 2*v(9) * zt;
        p1 = patch( isosurface( xt, yt, zt, Ellipsoid, -v(10) ) );
        hold off;
        set( p1, 'FaceColor', 'g', 'EdgeColor', 'none' );
        view( -70, 40 );
        axis vis3d equal;
        camlight;
        lighting phong;
        title('Corrected data');
        
        set(gcf, 'PaperUnits', 'inches');
        x_width=6 ;y_width=4;
        set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
        saveas(gcf,[bird,'_cal_sphere_plots_', fnum,'.png'])
        
        % IFF ALL LOOKS OK, RECOMBINE ALL OF THE SENSOR DATA WITH THE INDICES USED
        % FOR CALIBRATING THE MAG DATA
       
        % Subsample for start and end indices for q chunk (remember SensorData_s is already a subset
        % of the bookend indices that were chopped off outside loop
        SensorMat = SensorData_s(starts(q):ends(q),:);
        
        % Isolate the sensor variables, accelerometer, magnetometer,
        Acc = SensorMat(:,1:3);
        Mag = SensorMat(:,4:6);
        G = SensorMat(:,7:9);
        T = SensorMat(:,10);
        P = SensorMat(:,11);
        
        % Calibrated Mag data
        Mc = [xc,yc,zc];
        
        % IF this is the first chunk, write out a txt file with the
        % headers, otherwise, simply append the new calibrated chunk to the
        % data using writematrix (this function preserves precision)
        if q == 1
            % RECOMBINE and give column names
            S = [Acc, Mag, Mc, G, T, P];
            ST = array2table(S,'VariableNames',{'Ax', 'Ay', 'Az', 'Mx', 'My', 'Mz', 'Mxcal','Mycal','Mzcal','Gx','Gy','Gz','T','P'});
            
            %Save Calibrated Mag Data (.txt)
            % In folder called 'Calibrated_Mag_Data' that is inside Current Directory
            writetable(ST,['Calibrated_Mag_Data/',bird,'_Calibrated_Mag_xyz.txt'])
        else
            S = [Acc, Mag, Mc, G, T, P];
            writematrix(S,['Calibrated_Mag_Data/',bird,'_Calibrated_Mag_xyz.txt'],'WriteMode','append');
            
        end

    end
    
    
    % Write out the diagnostics table
    writetable(diags,['Calibrated_Mag_Data/',bird,'_Calibration_Diagnostics.csv'])
    
   
    % DO A GUT CHECK 
        clearvars -except i CurrentPath fileList sps timesplit SensorData_s bird
        test = readtable(['Calibrated_Mag_Data/',bird,'/',bird,'_Calibrated_Mag_xyz.txt']);

        Mag = table2array(test(:,4:6));
        Mag_cal = table2array(test(:,7:9));
        Mag_og = SensorData_s(:,4:6);

        % gut check plot write out mag data and original mag data - these should
        % look identical
        plott(Mag, 75, Mag_og,75) % check with plot - plot of unfiltered and filtered.

        % plot calibrated stuff
        plott(Mag, 75, Mag_cal,75)





