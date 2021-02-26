%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script "unpacks" 30-sec HMM-inferred states into full resolution
% sensor files. Objective is to plot raw sensor data with states, to
% inspect HMM-inferred behaviors as they relate to patterns in sensor data.

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
% hmmdir = set directory containing HMM output file
% metadir = set directory containing meta datatables required in script
% dropdir = set directory to store output

% Set path to required functions-------------------------------------------
% Add libraries into path.
addpath(genpath(strcat(fundir,"/Functions_Toolboxes/")))

% Directory for behavioral state HMM output summarized in fixed time
% windows:
% Model-30-2: 30s Acc+Mag+SppCovariate on Transition Probabilities only
hmm_data=readtable(strcat(hmmdir,'/model30_1.csv'));

% set current directory to sensor mats
CurrentPath= strcat(datadir,"/Extended_Sensor_Mats/");

cd(CurrentPath)
fileList=dir('*.txt'); 


% -------------------------------------------------------------------------------------
% set global parameters

fs=25;
nsec = 30;
w = nsec*fs;
str1 = '#F2AD00'; c1 = sscanf(str1(2:end),'%2x%2x%2x',[1 3])/255; % Colors used for plotting. Match those in R code
str2 = '#5BBCD6'; c2 = sscanf(str2(2:end),'%2x%2x%2x',[1 3])/255; % Colors used for plotting. Match those in R code.
str3 = '#00A08A'; c3 = sscanf(str3(2:end),'%2x%2x%2x',[1 3])/255; % Colors used for plotting. Match those in R code.

hmmbirds = unique(hmm_data.ID);
cmpcell={'BBAL','GHAL'};

cut = 7200/nsec; %two hours truncated off hmm dataset

% identify birds to plot
validationbirds = {'GHAL44','GHAL47','GHAL26','GHAL48','GHAL55','BBAL84','BBAL99','BFAL27','BFAL49','BBAL86','BBAL108','BBAL21','LAAL28','LAAL43'};


% -------------------------------------------------------------------------------------
% Loop through each sensor mat:


for i = 1:length(fileList)
    
   close all
   
        bird=strsplit(fileList(i).name,"_");
        bird=strcat(bird{1});
        spp = bird(1:4);
   
        % Does this bird have an hmm output? ---------------------------------------
        index = cell2mat(cellfun(@(a) strmatch(a,sensorlist{i}),hmmbirds,'uniform',false));
      
        if isempty(index)
            
        else
            
                % Import bird i 
                m = readtable(fileList(i).name,'Delimiter',',','ReadVariableNames',true,'TreatAsEmpty',{'NA'});   

                % set up indices ------------------------------------------------------
                winvec = 1:w:length(m.Ax);
                winvec(1:cut)=[]; % remove indices pertaining to first two hours of sensor file.
                
                % Create empty state column (make sure this is a matrix)
                statecol = NaN([length(m.Ax),1]);

                % import hmm results for file i ---------------------------------------
                
                hmm_i = hmm_data(strcmp(hmm_data.ID,char(sensorlist{i})),:);
                hmm_state = hmm_i.state;
                
                
                % Add state from HMM dataframe to sensormatrix
                for k = 1:length(hmm_i.ID) 
                    statecol(winvec(k):winvec(k+1)-1)=hmm_state(k);       
                end

                % Save statecolumn - save only this column bc if save full
                % matrix with added state column, you get datafiles 6GB+
                % large! 
           writematrix(statecol, strcat(dropdir,'/state_only/model30_1/',bird,'_statecolumn_fulllength.txt'),'delimiter',',')


               if ismember(bird,validationbirds)
               % Create figure 
                       if sum(strcmp(spp,cmpcell))==1
                            y=4; % Use a higher ranged y-axis if GHAL/BFAL bc sensor data had larger range
                        else
                            y=3;
                       end

                       Atot = [m.Ax, m.Ay, m.Az];
                       plot(Atot)
                       hold on

                       x=statecol;
                       Isame = @(x) [1; find(diff(x))+1; length(x)]; % index for same areas
                       Idx=Isame(x);


                       for k = 1:length(Idx)-1
                          state_k = statecol(Idx(k));

                          if isnan(state_k)

                          else
                                if state_k==1
                                    c_k=c1;
                                elseif state_k==2
                                    c_k=c2;
                                elseif state_k==3
                                    c_k=c3;
                                end

                                line([Idx(k),Idx(k+1)-1],[y,y],'LineWidth',15, 'Color',c_k)
                          end
                       end
                       
                       saveas(gcf,strcat(dropdir,'figs/2021_01_25/model30_1/',bird,'_hmm_statefig.fig'))
                       

               else
                   % skip figure
               end

        end
        

        clearvars -except sensordir fileList SensorList sensorlist hmmbirds hmm_data fs nsec w c1 c2 c3 i cmpcell cut validationbirds
end
   
   
   





