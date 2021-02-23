# HMMalbatross
This repository stores the code and an example dataset used in Conners et al 2021 Movement Ecology: "Hidden Markov models identify major movement modes in accelerometer and magnetometer data from four albatross species."

**This README and all files will be finalized and uploaded by Feb. 26 2021** \
Please address questions/comments to Melinda Conners: connersm@gmail.com

## Scripts: 

#### 's1_import_AGM.m' (MATLAB): to import and pre-process sensor data from Technosmart AGM devices
  - sensor frames are aligned to each other and the bird frame
  - times are checked for irregularities in sampling rate (25 Hz), and if irregularities exist, the data frame is expanded with NAs to result in a dataframe with regular observations.
  - temperature and pressure columns are interpolated from 1 Hz to 25 Hz

#### 's1_import_neurologger.m' (MATLAB): to import and pre-process sensor data from the Evolocus Neurologger devices
  - data from raw .bin files are decompressed and converted to relevant units.
  - data is organized into separate 75 Hz sensor data and 600 Hz ECG data frames.
  - sensor frames are aligned to each other and to the bird frame.
  
#### 's2_calibrate_magnetometer.m' (MATLAB): uses a data-driven, piece-wise method to calibrate triaxial magnetometer data from both tag types
Typically, the calibration of magnetometer data is done using calibration roll data that is collected prior to each device's deployment that is then used to correct the magnetometer data of each deployment. However, this approach is potentially problematic for animals traveling over great distances as this method assumes the magnetic field insensity where the calibration rolls were collected will remain somewhat static along the full deployment. Since albatrosses fly vast distances across many degrees of latitude and longitude, the magnetic field can change substantially along the trip. To address this, we used a data-driven approach recommended in Mark Johnson's "Measuring the orientation and movement of marine animals using inertial and magnetic sensors - a tutorial" (https://synergy.st-andrews.ac.uk/soundtags/files/2013/01/animal_orientation_tutorial.pdf). This method requires triaxial data and animals that change orientation substantially throughout the deployment. In this script, triaxial magnetometer data from each bird deployment is divided into segments and each segment is calibrated separately using magnetometer data within that segment. Calibration diagnostics are run on each segment and if the calibrated data does not pass a diagnostic test, then the duration of the segments will be reduced to a shorter time interval (described in detail in Additional File 2). 

#### 's3_neurologger_saveRotMat.m' (MATLAB): creates a rotation matrix for each neurologger deployment to correct a slight tilt in the sensor frames relative to the bird frame inherent in these devices. 

#### 's4_neurologger_prep4analysis.m' (MATLAB): final step for creating analysis-ready neurologger files
- units are conformed to that of AGM
- tag frame is untilted using the rotation matrix from 's3_neurologger_saveRotMat.m'
- 75 Hz neurologger data is decimated to 25 Hz to match AGM sampling frequency

#### 's4_AGM_prep4analysis.m' (MATLAB): final step for creating analysis-ready AGM files
All columns checked and edited to conform with that of neurologger files

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

#### 'x2_run_HMM_optim.r' (R environment):

#### 'z3_unpack_states_toFullSensorFiles.m' (MATLAB) and plot a subset of figures:
