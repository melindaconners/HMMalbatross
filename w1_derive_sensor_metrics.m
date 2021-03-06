%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script derives additional metrics to be summarized from the sensor
% data.
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
addpath(genpath(strcat(fundir,"/Functions_Toolboxes/")))

% Set directory Analysis_Ready data:
CurrentPath = strcat(datadir,'/L1_1_Analysis_Ready/25hz/');

cd(CurrentPath)

fileList = dir('*.txt'); 


tic()
for i = 1:length(fileList) 

       
    % Import bird i 
    m = readtable(fileList(i).name,'Delimiter',',','ReadVariableNames',true,'TreatAsEmpty',{'NA'});   
    bird=strsplit(fileList(i).name,"_");
    bird=strcat(bird{1});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART ONE: Create Additional Datastreams Derived from Sensordata
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Isolate A matrix and M matrix
    A = [m.Ax, m.Ay, m.Az];
    M = [m.Mx, m.My, m.Mz];
    fs = 25;
    
    %% Calculate Euler Angles - use low frequency M and A ----------------------------
    Mf=comp_filt(M,25,0.18); % filter to separate low and high frequency mag
    Mf=Mf{1,1}; % choose low frequency mag data
    clear M     % delete M because only will be using low frequency M moving forward
    Af=comp_filt(A,25,0.45); % filter to separate low and high frequency acc
    Af=Af{1,1}; % choose low frequency acc data for euler angle derivation
    
    [Pitch,Roll,~] = a2pr(Af);  % calculate pitch and roll
    [Heading,mm,incl] = m2h(Mf, Pitch, Roll);    % calculate heading
    h360=wrapTo360(rad2deg(Heading));   % convert heading to degrees
    pitch_deg=rad2deg(Pitch);           % convert pitch to degrees
    roll_deg = rad2deg(Roll);           % convert roll to degrees
    

    %% Classic Acceleration Metrics-----------------------------------------
    ns = 2; % number of seconds in smoothing window
    Astc = lowpassFilt(A, 'moveAvrg' , ns*fs); % static acceleration 
    Adyn = A - Astc;                           % dynamic acceleration 

    % ODBA and VeDBA
    odba = abs(Adyn(:,1)) + abs(Adyn(:,2)) + abs(Adyn(:,3)); % the wilson method 
    vedba = sqrt(Adyn(:,1).^2 + Adyn(:,2).^2 + Adyn(:,3).^2); % the vector norm

    
    %% Soaring/Mag Euler Body Rotation Metrics-----------------------------------------
    % AVEY, AVEP, AVER, AAV Gunner 2020
    
    % Calculate column of differentials with a stepping range of 1/2
    % second. Ideal sampling window to calculate AVeY should consider animal in question.
    % Window should be less than the time it takes the animal to rotate
    % through half a revolution. 
    
     ns = .5; % ns in window
     wavwin = round(ns*fs)-1; % half second time window
     
     for k = wavwin+1:length(h360)
       aver(k,1)= roll_deg(k)-roll_deg(k-wavwin);   % roll differential over time window
       avep(k,1)= pitch_deg(k)-pitch_deg(k-wavwin); % pitch differential over time window
       hdiff(k,1)= h360(k)-h360(k-wavwin);          % heading differential over time window
      end 
    
    % Correct Avey (hdiff) to reflect that it is far more likely that an animal
    % that caused the compass heading to change from 10 to 350 turned 20
    % degrees anticlockwise, rather than 330 degrees closckwise. 
    
    % This is slightly modified from Gunner et al because albatross do
    % contain some actual turns that are very large and quick. We need to
    % filter not just by looking for differentials (hdiff) > 180 but also for jumps
    % in the differential (hdiff2) that are > 180 that would reflect too
    % fast (a 1/25 sec) change of 180 degrees in heading. 
    
    hdiff2 = diff(hdiff);
    
    % Identify where hdiff > 180 or <-180 AND where there is an hdiff2 >
    % 180 within 1/2 second of i. If true, correct.
    
    % Add a column 'jumpcol' that identifies a large jump (>170degrees) in the
    % differential
    jumpcol=NaN(length(hdiff),1);
  
    for l=wavwin+1:length(hdiff2)-wavwin
        chunkl=hdiff2(l-wavwin:l+wavwin);
        if ~isempty(find(chunkl>170)) % if there is a jump within wavwin, jumpcol = 1
            jumpcol(l,1)=1;
        else
            jumpcol(l,1)=0; % else jumpcol = 0
        end    
    end
    
    clear l
    
    % Using both conditions now, correct hdiff. 
    for l=1:length(hdiff) 
        
        if hdiff(l)<-170  && jumpcol(l)==1 % If there is a large jump and differential is <-170, correct by +adding 360
            avey(l,1)=hdiff(l)+360;
        elseif hdiff(l) >170  && jumpcol(l)==1 % If there is a large jump and differential is >170, correct by -subtracting 360
            avey(l,1)=hdiff(l)-360;  
        else 
            avey(l,1)=hdiff(l); % If there is no large jump, avey does not need to be corrected
        end
    
    end

    
    % Absolute Angular Velocity (AAV):
    AAV = sqrt(avep.^2 + aver.^2 + avey.^2);
    
    % Save extended sendor dataframe with additional metrics: pitch_deg, roll_deg,
    % h360, avep, aver, avey, aav, odba, vedba, Adyn, Astc
    
    m.pitch_deg = pitch_deg; % pitch in degrees (-90-90)
    m.roll_deg = roll_deg;   % roll in degrees 
    m.h360 = h360;           % heading in degrees 0-360
    m.avep = avep;           % angular velocity about the pitch axis
    m.aver = aver;           % angular velocity about the roll axis
    m.avey = avey;           % angular velocity about the yaw axis
    m.odba = odba;           % odba
    m.vedba = vedba;         % vedba
    m.Adyn = Adyn;           % Dynamic Acceleration (_1,_2,_3)
    m.Astc = Astc;           % Static Acceleration (_1,_2,_3)
    
    % save as a text file
    writetable(m,strcat(datadir,bird,'_extsensormat.txt'),'delimiter',',')

   clearvars -except fileList i CurrentPath bird 
    
end
toc()
