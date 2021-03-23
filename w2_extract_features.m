%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Feature Extraction for Hidden Markov Models - end of June 2020

% This code is used to summarize sensordata within a fixed timewindow resulting in 
% data to be used by both the HMMs for behavioral classification 
% 
% Code by M. Conners for Conners et al 2021:
% "Hidden Markov models identify major movement modes in accelerometer and
% magnetometer data from four albatross species." Movement Ecology
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Interested in the following windows of time: 
    % 30, 15, 10
    
% Features to extract for the Hidden Markov Model analysis: 
        % 1. Dominant Frequency in Dynamic Heave Acceleration
        % 2. Highest Dominant Frequency in Dynamic Heave Acceleration 
        % 3. Mean Static Heave Acceleration
        % 4. Standard Deviation of Static Heave Acceleration
        % 5. Highest 10% Quantile in Static Heave Acceleration
        % 6. Circular Variance of the Heading Signal
        % 7. Circular Deviation of the Heading Signal
        % 8. IQR dynamic accel - sway
        % 9. IQR dynamic accel - heave
 
% Metrics to extract for HR ~ Metrics_EE analysis:
    % Traditional EE Metrics:
        % 10. Mean ODBA
        % 11. Mean VEDBA
    
    % Flapping Specific Metrics:
        % 12. Flapping Rate
        % 13. Mean Flapping Amplitude
    
    % Soaring Specific Metrics:
        % 14. Mean AAV
        % 15. Mean AvEY Degrees (deg/second)
        % 16. N Swoop Cycles
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set folder directories
% fundir = set directory containing Functions_Toolboxes
% datadir = set directory containing compressed neurologger data
% dropdir = set directory to store output files.

% Set path to required functions-------------------------------------------
% Add libraries ('d3matlab','tagmatlab','x3toolbox','tagtools') into path.
addpath(genpath(strcat(fundir,"Functions_Toolboxes/")))

% Set working directory to Neurologger Sensor Data--------------------------------
% Raw Compressed dat
CurrentPath = strcat(datadir,'/Extended_Sensor_Mats/');

cd(CurrentPath)


%% Feature Extraction 

fileList = dir('*.txt'); 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% routine needed to clear ghost files if working in Seagate Hard Drive 
skipx=[];
for k = 1:length(fileList)
skipx(k)=~startsWith(fileList(k).name,'._','IgnoreCase',true);
end

fileList = fileList(find(skipx));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   

fs=25; %sampling rate

% Loop through each bird
for i=1:length(fileList) 
tic()
        % Import bird i 
        m = readtable(fileList(i).name,'Delimiter',',','ReadVariableNames',true,'TreatAsEmpty',{'NA'});   
        bird=strsplit(fileList(i).name,"_");
        bird=strcat(bird{1});
                
        A = [m.Ax, m.Ay, m.Az];


        % Loop Through Each Time Window and Extract Features
       for n=[10,15,30] % 10, 15, and 30 sec windows
        
                    
                    % Set window params----------------------------------------------------
                    nsec = n; % number of seconds in window
                    w = nsec*fs;
                    winvec = 1:w:length(A);

                    % Set Metric Variables-------------------------------------------------
                    nrow= length(1:length(winvec)-1);   

                    % Features to extract for the Hidden Markov Model analysis:    
                    out1 = NaN(nrow,1);     % 1. Dominant Frequency in Total Heave Acceleration
                    out2 = NaN(nrow,1);     % 2. Highest Dominant Frequency in Total Heave Acceleration 
                    out3 = NaN(nrow,1);     % 3. Mean Static Heave Acceleration
                    out4 = NaN(nrow,1);     % 4. Standard Deviation of Static Heave Acceleration
                    out5 = NaN(nrow,1);     % 5. Highest 95% Quantile in Static Heave Acceleration
                    out6 = NaN(nrow,1);     % 6. Circular Variance of the Heading Signal
                    out7 = NaN(nrow,1);     % 7. Circular Deviation of the Heading Signal
                    out8 = NaN(nrow,1);     % 8. IQR dynamic accel - sway
                    out9 = NaN(nrow,1);     % 9. IQR dynamic accel - heave
                    % Traditional EE Metrics:
                    out10 = NaN(nrow,1);    % 10. Mean ODBA
                    out11 = NaN(nrow,1);    % 11. Mean VEDBA
                    % Flapping Specific Metrics:
                    out12 = NaN(nrow,1);    % 12. Flapping Rate
                    out13 = NaN(nrow,1);    % 13. Mean Flapping Amplitude
                    % Soaring Specific Metrics:
                    out14 = NaN(nrow,1);    % 14. Mean Swoop Degrees (mean of Avey Peak Values)
                    out15 = NaN(nrow,1);    % 15. Max Swoop Degrees (max of Avey Peak Values)
                    out16 = NaN(nrow,1);    % 16. N Swoop Cycles (number of peaks) >> 10s window is likely too small for soaring metrics.


                    % ---------------------------------------------------------------------  
                    % EXTRACT FEATURES:
                    % Set up FORLOOP to loop through moving time window
                    % Use parfor to maximize computer cores
                    % ---------------------------------------------------------------------  
                    tic()
                 
                   for j = 1:length(winvec)-1
                      ixs = winvec(j):winvec(j+1)-1;

                      if length(find(all(A(ixs,:) == 0,2))) > 2*fs 
                         % BBAL99 has two chunks with missing Acc data
                         % for any birds, like BBAL99, that somehow lost the acceleration signal for more than 2 seconds (50 observations), leave
                         % NaNs as output. 
                      else

                         %---------------------------------------------------------------------- 
                         % 1. Features to extract for the Hidden Markov Model analysis:    
                         %---------------------------------------------------------------------- 

                         % -----------------------------------------------------------
                         % Highest Dominant Frequency in Dynamic Heave Acceleration 

                         % Use a FAST FOURIER TRANSFORM to identify dominant frequencies in the total heave
                         % acceleration signal:

                        data = A(ixs,3); % Total Heave
                        data_dt = detrend(data,'constant');
                        data_n = numel(data);
                        PSD = abs(fft(data_dt, data_n))/data_n; % power spectral density
                        freq_scale = ((0:data_n/2-1)/data_n*fs)'; % frequency scale
                        plot(freq_scale,PSD(1:length(freq_scale)))

                        % Find dominant frequency in signal
                        [max_amp, max_freq] = max(PSD(1:length(freq_scale))); %maximum power spectral density
                        dom_freq = freq_scale(max_freq);

                        % Find multiple dominant frequencies identified in FFT

                         % smooth curve 
                            method = 'gaussian';
                            window = 30;
                            fft_s = smoothdata(PSD(1:length(freq_scale)),method,window);

                            % Find peaks in smoothed curve
                            [pks,locs,wdth,prom] = findpeaks(fft_s,freq_scale,'MinPeakDistance',1.5);

                            % Only keep prominent peaks
                            peak_ixs = find(pks>mean(fft_s)+(std(fft_s)/2));

                            if isempty(peak_ixs)
                                peak_ixs = find(fft_s == max(fft_s));

                            end

                         % Count Major Peaks 
                        npks = length(peak_ixs);
                        peaks = [locs(peak_ixs),pks(peak_ixs)];
                        total_amp = sum(peaks(:,2));
    
                        % What is Frequency at Dominant Peak?   
                        high_ix = find(peaks(:,1)==max(peaks(:,1)));
                        dom_ix  = find(peaks(:,2)==max(peaks(:,2)));
    
                        dom_freq        = peaks(dom_ix,1);
                        highest_freq    = peaks(high_ix,1);

                        % -----------------------------------------------------------
                        % FFT Frequency Domain Metrics       

                        out1(j,1) = dom_freq;                              % 1. dominant frequency component
                        out2(j,1) = highest_freq;                          % 2. frequency of peak at highest frequency

                        % -----------------------------------------------------------
                        % Mean and STD Static Heave Acceleration

                        out3(j,1) = nanmean(m.Astc_3(ixs));              % 3. Mean Static Heave Acceleration
                        out4(j,1) = nanstd(m.Astc_3(ixs));               % 4. Standard Deviation of Static Heave Acceleration  

                        % -----------------------------------------------------------
                        % Highest 95% Quantile in Static Heave Acceleration

                        qtl = quantile(m.Astc_3(ixs),[0.05 0.95]);
                        out5(j,1) = qtl(2);                               % 5. Top 5% Quantile in Static Heave Acceleration   

                        % -----------------------------------------------------------
                        % Circular Variance and Angular Deviation of the Heading Signal

                        out6(j,1) = rad2deg(circ_var(deg2rad(m.h360(ixs))));  % 6. Circular Variance of the Heading Signal
                        out7(j,1) = rad2deg(circ_std(deg2rad(m.h360(ixs))));  % 7. Circular Deviation of the Heading Signal 

                        % -----------------------------------------------------------
                        % IQR Dynamic Acceleration in the Sway and Heave Axes

                        out8(j,1) = iqr(m.Adyn_2(ixs));                  % 8. IQR dynamic accel - sway
                        out9(j,1) = iqr(m.Adyn_3(ixs));                  % 9. IQR dynamic accel - heave    

                        %---------------------------------------------------------------------- 
                        % 2. Traditional EE Metrics: Mean ODBA and VeDBA
                        %----------------------------------------------------------------------     

                        out10(j,1) = nanmean(m.odba(ixs));                % mean odba
                        out11(j,1) = nanmean(m.vedba(ixs));               % mean odba

                        %---------------------------------------------------------------------- 
                        % 3. Flapping-Specific EE Metrics:
                        %----------------------------------------------------------------------     
                        % if flapping:  
                          if highest_freq > 2.2 
                                              
                             heave = m.Adyn_3(ixs)*-1;% Flipping the timeseries because the downstroke is so much more pronounced.

                             % Find peaks in smoothed curve
                             MAXW = .30;% ignore peaks that are wider than MAXW (and are likely associated with soaring etc).
                             [flaps,locs2,wdth2,flap_amp] = findpeaks(heave,fs,'MinPeakDistance',.2,'MaxPeakWidth',MAXW,'MinPeakProminence',0.6);     

                             %flapping rate: flaps per minute:
                             flap_rate=length(flaps)*(60/nsec);

                             %mean flapping amplitude
                             flap_amp = mean(flap_amp);

                             out12(j,1) = flap_rate;                       % 12. Flapping Rate
                             out13(j,1) = flap_amp;                        % 13. Mean Flapping Amplitude

                          else

                             out12(j,1) = NaN;                              
                             out13(j,1) = NaN;                              
                          end           

                          %---------------------------------------------------------------------- 
                          % 4. Soaring-Specific EE Metrics:
                          %----------------------------------------------------------------------        

                          avey_seg=m.avey(ixs);
                          avey_seg_cen = detrend(avey_seg);

                          % local max separated by a min of 2 seconds
                          % local min separated by a min of 2 seconds

                          % local minima
                          TFmn = islocalmin(avey_seg_cen, 'MinSeparation',50,'FlatSelection','first');

                          % local maxima
                          TFmx = islocalmax(avey_seg_cen, 'MinSeparation',50,'FlatSelection','first');

                          out14(j,1) = mean([abs(avey_seg_cen(TFmn));abs(avey_seg_cen(TFmx))]);   % 14. Mean Swoop Degrees (mean of Avey Peak Values)
                          out15(j,1) = max([abs(avey_seg_cen(TFmn));abs(avey_seg_cen(TFmx))]);    % 15. Max Swoop Degrees (max of Avey Peak Values)           
                          out16(j,1) = length(avey_seg_cen(TFmn))+length(avey_seg_cen(TFmx));     % 16. N Swoop Cycles (number of peaks) >> 10s window is likely to small for soaring metrics. 
                      end
       
                   end
                      
                % create array of features extracted from bird i
                df=[out1,out2,out3,out4,out5,out6,out7,out8,out9,out10,out11,out12,out13,out14,out15,out16];

                % format array as datatable and provide column names
                T = array2table(df) ;
                T.Properties.VariableNames = {'df_heave','hf_heave', 'm_stheave', 'sd_stheave','q95_stheave','var_head', 'sd_head','iqr_sway', 'iqr_heave', 'm_odba', 'm_vedba', 'flap_rate', 'flap_amp', 'm_avey', 'mx_avey', 'n_swoops'}; 

                % save as a text file
                writetable(T,strcat(dropdir,num2str(n),'s/',bird,'_',num2str(n),'s_FeatExtr.txt'),'delimiter',',')
                
                % clear variables except global variables before looping into next bird.
                clearvars -except fileList i fs CurrentPath bird m A 

         end % Ends the 10,15,30 second loop
        
        clearvars -except fileList i fs CurrentPath 
toc()    
end % Ends the loop for birdi
    


  



