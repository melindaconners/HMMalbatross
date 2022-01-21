# HMMalbatross
This repository stores the code and an example dataset used in Conners et al 2021 Movement Ecology: "Hidden Markov models identify major movement modes in accelerometer and magnetometer data from four albatross species."

Please address questions/comments to Melinda Conners: connersm@gmail.com

## Example dataset:
albatross-example-featurized30s-dataset.csv: This truncated dataset consists of featurized data from 4 individuals (2 black-browed albatross and 2 grey-headed albatross). Features include 'hf', 'p5',  and 'sh' (descriptions of their derivation in manuscript). Features were derived from 25 Hz data but were summarized to a 30-s timescale. Two additional columns for 'species' and 'ID' are also included. 

## Functions/Toolboxes:
For the following scripts to work, need to load functions and toolboxes into Matlab path. Files starting with 'ft' should be included in a Functions_Toolboxes folder.

## Scripts: 

#### 's1_import_AGM.m' (MATLAB): to import and pre-process sensor data from Technosmart AGM devices
This code loops through each bird deployment folder and appends all the files associated with that deployment. Then writes a single file for each deployment

#### 's1_import_neurologger.m' (MATLAB): to import and pre-process sensor data from the Evolocus Neurologger devices
  - data from raw .bin files are decompressed and converted to relevant units.
  - data is organized into separate 75 Hz sensor data and 600 Hz ECG data frames.
  - sensor frames are aligned to each other and to the bird frame.

#### 's2_calibrate_magnetometer.m' (MATLAB): uses a data-driven, piece-wise method to calibrate triaxial magnetometer data from both tag types
Typically, the calibration of magnetometer data is done using calibration roll data that is collected prior to each device's deployment that is then used to correct the magnetometer data of each deployment. However, this approach is potentially problematic for animals traveling over great distances as this method assumes the magnetic field insensity where the calibration rolls were collected will remain somewhat static along the full deployment. Since albatrosses fly vast distances across many degrees of latitude and longitude, the magnetic field can change substantially along the trip. To address this, we used a data-driven approach recommended in Mark Johnson's "Measuring the orientation and movement of marine animals using inertial and magnetic sensors - a tutorial" (https://synergy.st-andrews.ac.uk/soundtags/files/2013/01/animal_orientation_tutorial.pdf). This method requires triaxial data and animals that change orientation substantially throughout the deployment. In this script, triaxial magnetometer data from each bird deployment is divided into segments and each segment is calibrated separately using magnetometer data within that segment. Calibration diagnostics are run on each segment and if the calibrated data does not pass a diagnostic test, then the duration of the segments will be reduced to a shorter time interval (described in detail in Additional File 2). 

#### 's3_neurologger_saveRotMat.m' (MATLAB): 
creates a rotation matrix for each neurologger deployment to correct a slight tilt in the sensor frames relative to the bird frame inherent in these devices. 

#### 's4_AGM_prep4analysis.m' (MATLAB):  final step for creating analysis-ready AGM files
- Check for upside-down tag placements
- Rotate M frame to match A frame and bird frame
- Check time intervals and expand with NAs if necessary
- Interpolate P and T
  
#### 's4_neurologger_prep4analysis.m' (MATLAB): final step for creating analysis-ready neurologger files
- units are conformed to that of AGM
- tag frame is untilted using the rotation matrix from 's3_neurologger_saveRotMat.m'
- 75 Hz neurologger data is decimated to 25 Hz to match AGM sampling frequency

#### 'w1_derive_sensor_metrics.m' (MATLAB): Additional 25 Hz metrics are derived from accelerometer and magnetometer data 
- Euler angles: *pitch*, *roll*, *heading (yaw)*
- Angular velocities: *AvEP*, *AvER*, *AvEY*, *AAV* (see Gunner et al 2020)
- *ODBA*, *VeDBA* (traditional acceleration-derived energetic proxies)
- acceleration components: *dynamic acceleration*, *static acceleration*

#### 'w2_extract_features.m' (MATLAB): This script summarizes the 25 Hz sensor data into features using fixed time windows (10s, 15s, 30s). The following features are extracted: 
1.	Dominant frequency in dynamic heave acceleration
2.	Highest dominant frequency in dynamic heave acceleration
3.	Mean static heave acceleration
4.	Standard deviation of static heave acceleration
5.	95th quantile of static heave acceleration
6.	Circular variance of heading
7.	Circular deviation of heading
8.	IQR (inter-quartile range) of dynamic sway acceleration
9.	IQR of dynamic heave acceleration
10.	Mean ODBA
11.	Mean VeDBA
12.	Flapping Rate
13.	Mean Flapping Amplitude
14.	Mean AvEY Peak values (mean soaring arc amplitude (in degrees)
15.	Max AvEY peak values (max soaring arc amplitude)
16.	Number of AvEY peaks (number of soaring arcs)

#### 'x1_explore_feature_content.r' (R environment):
Explore features to determine final set of features. Plots histograms and creates correlation matrix

#### 'x2_run_HMM_optim.r' (R environment):
Construct HMM for Model-1: 3 state model on 3 features ('hf', 'p5', 'sh')
Par0 are optimized by identifying starting values resulting in best fit model from 25 iterations
Species is included as a fixed effect on transition probabilities.

#### 'z1_unpackHMMstate.m' (MATLAB): 
This script "unpacks" 30-sec HMM-inferred states into full resolution
sensor files. Objective is to plot raw sensor data with states, to
inspect HMM-inferred behaviors as they relate to patterns in sensor data.
