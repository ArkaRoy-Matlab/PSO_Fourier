%Example code for getting depth profile from observed gravity anomalies and
%density contrast for any real case study inversion of sedimentary basin.

% ***************************************************************
% *** Matlab function for inversion of sedimentary basin from observed
%     gravity anomalies.
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Dr. Chandra Prakash Dubey (email:p.dubey48@gmail.com)
% ***       Mr. M. Prasad (email:prasadgoud333@gmail.com)
% ***       Crustal Processes Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ****************************************************************
%%Matlab fuction for inversion of sedimentary basin using PSO and adaptive

clear all
close all

%horizontal location of observation points for Godavari basin
%data files for observation points and gravity anomaly must be in single collumn array with multiple rows 

%horizontal location of observation points
x_obs_data='x_obs_Zhou2012_3a.dat';     

%corresponding gravity anomalies at observation points
gravity_anomaly_data='gravity_anomaly_Zhou2012_3a.dat';

%function for density contrast 
density=@(z) (-0.4692*exp(-4.078*10^-4*z)).*1000;

%Max depth 
max_depth=6000;
%calling the function for inversion
[depth_inverted,grav_inverted]=Fourier_PSO(x_obs_data,gravity_anomaly_data,density,max_depth);


            %%%%%%%%%%%%%%%%% End of Code %%%%%%%%%%%%%%%%%%

