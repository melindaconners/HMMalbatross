% s1_AGM_import

% This code loops through each bird deployment folder and appends all the files
% associated with that deployment. Then writes a single file for each deployment.
% Add path
addpath(genpath("~/Dropbox/Academia/SUNY/Project_Components/BIRD_ISLAND/Analyses/Functions_Toolboxes/"))

% Set Working Directory to Bird Island 2019-2020 Technosmart Accelerometer
% Data
% Create a structure containing tag type and the respective directories:
acc_dir_list(1).tagtype = "AGM";
acc_dir_list(1).accdir = '/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Bird_Island/2019_2020/Tag_Data/L0_Raw_Data/AGM/L0_extracted_tag_data/';
acc_dir_list(2).tagtype = "AxyAir";
acc_dir_list(2).accdir = '/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Bird_Island/2019_2020/Tag_Data/L0_Raw_Data/AxyAir/L0_extracted_tag_data/';
acc_dir_list(3).tagtype = "AxyTrek";
acc_dir_list(3).accdir = '/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Bird_Island/2019_2020/Tag_Data/L0_Raw_Data/AxyTrek/L0_extracted_tag_data/';

% Set drop directory for writing files
dropdir = '/Volumes/GoogleDrive/My Drive/Thorne Lab/THORNE_LAB/Data!/Conners_Bird_Island/2019_2020/Tag_Data/L1/2018-2019_ACC_Technosmart_Accel_Uniformat/';
     
% Loop through each tag type and import:

for q=1:length(acc_dir_list)
    
   tagtype_q = acc_dir_list(q).tagtype;
   
   % Set directory to tag type q
   CurrentPath =  acc_dir_list(q).accdir;
   cd(CurrentPath) 
   
   % ----------------------------------------------------------------------
   % Deal with Occassional Cases of Multiple Files per Bird
   % ----------------------------------------------------------------------
   % some birds have multiple files that need to be appended so need to
   % create list of bird names:

   fileList = dir('*.csv');
   cellList = extractfield(fileList,'name');
   parsecells = cellfun(@(x) strsplit(x, '_'),cellList,'UniformOutput',false);
   for g = 1:length(parsecells)
    filesbyBird{g} = [parsecells{1, g}{1},'_',parsecells{1, g}{2}];
   end

   uniquebirds = unique(filesbyBird);
   
   % Loop through each unique bird ----------------------------------------
   for i = 1:length(uniquebirds)
      birdi = uniquebirds(i);
      match = find(strcmp(filesbyBird,birdi));
      
      % If match is greater than 1, there are multiple files for this bird
      % that need to be appended. 
      if length(match)>1
          
          for j = 1:length(match) % loop through each file for birdi
              mj = readtable(fileList(match(j)).name,'Delimiter',',','ReadVariableNames',true,'TreatAsEmpty',{'NA'}); % read agm bird i's file j
              % remove problematic column 'metadata':
              colnix = find(strcmp(mj.Properties.VariableNames,'Metadata'));
              mj(:,colnix)=[];
              
              if j==1
                  m = mj;
              else
                  m=[m;mj];
              end
          end
          
      else
        m = readtable(fileList(match).name,'Delimiter',',','ReadVariableNames',true,'TreatAsEmpty',{'NA'}); % read agm bird i's file j
        % remove problematic column 'metadata':
        colnix = find(strcmp(m.Properties.VariableNames,'Metadata'));
        m(:,colnix)=[];
      
      end 
      
   % ----------------------------------------------------------------------
   % Check Out Columns and Rearrange/Reorder/AddFiller If Needed
   % ----------------------------------------------------------------------
   
   % Desired output across all tagtypes:
            %    TagID
            %    DateTime
            %    Ax
            %    Ay
            %    Az
            %    Mx (NA if acc only)
            %    My (NA if acc only)
            %    Mz (NA if acc only)
            %    Pressure
            %    Temperature
            %    Activity (NA if axyair)
   
      
      % Reformat Dates and Times to a DateTime in ISO 8601
      
      DateTime = m.Date + m.Time;
      DateTime.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
      TagID = repmat(birdi,[height(m),1]);
      Pressure = m.Pressure;
      Temperature = m.Temp___C_;
          
      if tagtype_q == "AGM"
           %columns=TagID,Date,Time,accX,accY,accZ,compX,compY,compZ,Activity,Pressure,Temp. (?C),Battery Voltage (V)
           
           Ax = m.accX; 
           Ay = m.accY;
           Az = m.accZ;
           Mx = m.compX;
           My = m.compY;
           Mz = m.compZ;
           Activity = m.Activity;
           
           clear m

      elseif tagtype_q == "AxyAir"
          %columns=TagID,Date,Time,X,Y,Z,Pressure,Temp. (?C),Battery Voltage (V)
          
           Ax = m.X; 
           Ay = m.Y;
           Az = m.Z;
           Mx = NaN(height(m),1);
           My = NaN(height(m),1);
           Mz = NaN(height(m),1);
           Activity = NaN(height(m),1);
           
           clear m

      else % AxyTrek
          %columns=TagID,Date,Time,X,Y,Z,Activity,Pressure,Temp. (?C),location-lat,location-lon,height-msl,ground-speed,satellites,hdop,signal-strength,Sensor Raw,Battery (V)
      
          Ax = m.X;
          Ay = m.Y;
          Az = m.Z;
          Mx = NaN(height(m),1);
          My = NaN(height(m),1);
          Mz = NaN(height(m),1);
          Activity = m.Activity;
          
          clear m
      
      end
      
     
      
      
     T = table(TagID,DateTime,Ax,Ay,Az,Mx,My,Mz,Pressure,Temperature,Activity);

     % Save T 
     filename = strcat(dropdir,birdi,'_2018-2019_',tagtype_q,'_25hz-raw.txt');
     writetable(T,filename)

     clearvars -except fileList cellList parsecells filesbyBird uniquebirds i tagtype_q CurrentPath q dropdir acc_dir_list
   
   end % Loop to next bird
   
   clear filesbyBird
    
end % End loop through tagtype q


