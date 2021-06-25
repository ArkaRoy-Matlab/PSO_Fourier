%Example code for getting depth profile from observed gravity anomalies
%with incorrect density contrast for all synthetic models 

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
%%Matlab code for inversion of all synthetic basins with incorrect density
%%contrast
%%Effect of incorrect density on this algorithm for model 1 fixed density
clear all
close all
%location of true depth for model1 fixed density
%%function for Synthetic Example of depth of the basin
z_peaks=peaks(100);
y_peaks=((z_peaks(71,:))+0.5)*200;

%t and c are Legendre Gaussian quadrature points for numerical integration
    [t_leg,c_leg]=lgwt(10,0,1); 
    %synthetic depth and observation points
    x_obs=linspace(0,5000,100);
    depth=y_peaks;
    
    %function for incorrect density contrast 
    density1=@(z) -800+0.*z+(-800*0.1);
    %function for correct density contrast 
    density2=@(z) -800+0.*z;
    z_obs=0;   %height of observation point is in surface
    %Closed polygonic profile of the basin
    xx1=[x_obs x_obs(end) 0];
    yy1=[depth 0 0];
    %Gravity anomaly for fixed density model 
    gg=poly_gravityrho(x_obs,z_obs,xx1,yy1,density2,t_leg,c_leg);
    %saving temporarily all data
    x_obs=x_obs'; gg=gg';
    save x_obs_data.dat x_obs;
    save gravity_anomaly_data.dat gg;

 
    %horizontal location of observation points
    x_obs_data='x_obs_data.dat';     
    %corresponding gravity anomalies at observation points
    gravity_anomaly_data='gravity_anomaly_data.dat';   
    data=importdata(gravity_anomaly_data);
    
    %maximum depth bound
    max_depth=6000;
    %calling the function for inversion with incorrect density
    [depth_inverted,grav_inverted]=Fourier_PSO(x_obs_data,gravity_anomaly_data,density1,max_depth);
    %error norm
    err_depth=norm((depth-depth_inverted)/depth)*100;
    err_gravity=norm((data-grav_inverted')/data)*100;
    fprintf('For Model-1 fixed density the relative percentage error in depth profile is %f\n',err_depth)
    fprintf('For Model-1 fixed density the relative percentage error in gravity anomaly profile is %f\n',err_gravity)
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%Effect of incorrect density on this algorithm for model 1 variable density
clear all
close all
%location of true depth for model1 fixed density
%%function for Synthetic Example of depth of the basin
z_peaks=peaks(100);
y_peaks=((z_peaks(71,:))+0.5)*200;

%t and c are Legendre Gaussian quadrature points for numerical integration
    [t_leg,c_leg]=lgwt(10,0,1); 
    %synthetic depth and observation points
    x_obs=linspace(0,5000,100);
    depth=y_peaks;
    
    %function for incorrect density contrast 
    density1=@(z) (-0.38-0.42*exp(-0.5*z*10^-3))*1000 + ((-0.38-0.42*exp(-0.5*z*10^-3))*1000)*0.1; 
    %function for correct density contrast 
    density2=@(z) (-0.38-0.42*exp(-0.5*z*10^-3))*1000; 
    z_obs=0;   %height of observation point is in surface
    %Closed polygonic profile of the basin
    xx1=[x_obs x_obs(end) 0];
    yy1=[depth 0 0];
    %Gravity anomaly for fixed density model 
    gg=poly_gravityrho(x_obs,z_obs,xx1,yy1,density2,t_leg,c_leg);
    %saving temporarily all data
    x_obs=x_obs'; gg=gg';
    save x_obs_data.dat x_obs;
    save gravity_anomaly_data.dat gg;

 
    %horizontal location of observation points
    x_obs_data='x_obs_data.dat';     
    %corresponding gravity anomalies at observation points
    gravity_anomaly_data='gravity_anomaly_data.dat';   
    data=importdata(gravity_anomaly_data);
    
    %maximum depth bound
    max_depth=6000;
    %calling the function for inversion with incorrect density
    [depth_inverted,grav_inverted]=Fourier_PSO(x_obs_data,gravity_anomaly_data,density1,max_depth);
    %error norm
    err_depth=norm((depth-depth_inverted)/depth)*100;
    err_gravity=norm((data-grav_inverted')/data)*100;
    fprintf('For Model-1 variable density the relative percentage error in depth profile is %f\n',err_depth)
    fprintf('For Model-1 variable density the relative percentage error in gravity anomaly profile is %f\n',err_gravity)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%Effect of incorrect density on this algorithm for model 2 fixed density
clear all
close all
%location of true depth for model1 fixed density
%%function for Synthetic Example of depth of the basin
z_peaks=peaks(160);
yy=((z_peaks(80,:)))*400;
y_peaks=yy(61:160);
%t and c are Legendre Gaussian quadrature points for numerical integration
    [t_leg,c_leg]=lgwt(10,0,1); 
    %synthetic depth and observation points
    x_obs=linspace(0,5000,100);
    depth=y_peaks;
    
    %function for incorrect density contrast 
    density1=@(z) -800+0.*z+(-800*0.1);
    %function for correct density contrast 
    density2=@(z) -800+0.*z;
    z_obs=0;   %height of observation point is in surface
    %Closed polygonic profile of the basin
    xx1=[x_obs x_obs(end) 0];
    yy1=[depth 0 0];
    %Gravity anomaly for fixed density model 
    gg=poly_gravityrho(x_obs,z_obs,xx1,yy1,density2,t_leg,c_leg);
    %saving temporarily all data
    x_obs=x_obs'; gg=gg';
    save x_obs_data.dat x_obs;
    save gravity_anomaly_data.dat gg;

 
    %horizontal location of observation points
    x_obs_data='x_obs_data.dat';     
    %corresponding gravity anomalies at observation points
    gravity_anomaly_data='gravity_anomaly_data.dat';   
    data=importdata(gravity_anomaly_data);
    
    %maximum depth bound
    max_depth=6000;
    %calling the function for inversion with incorrect density
    [depth_inverted,grav_inverted]=Fourier_PSO(x_obs_data,gravity_anomaly_data,density1,max_depth);
    %error norm
    err_depth=norm((depth-depth_inverted)/depth)*100;
    err_gravity=norm((data-grav_inverted')/data)*100;
    fprintf('For Model-2 fixed density the relative percentage error in depth profile is %f\n',err_depth)
    fprintf('For Model-2 fixed density the relative percentage error in gravity anomaly profile is %f\n',err_gravity)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%Effect of incorrect density on this algorithm for model 2 variable density
clear all
close all
%location of true depth for model1 fixed density
%%function for Synthetic Example of depth of the basin
z_peaks=peaks(160);
yy=((z_peaks(80,:)))*400;
y_peaks=yy(61:160);
%t and c are Legendre Gaussian quadrature points for numerical integration
    [t_leg,c_leg]=lgwt(10,0,1); 
    %synthetic depth and observation points
    x_obs=linspace(0,5000,100);
    depth=y_peaks;
    
    %function for incorrect density contrast 
    density1=@(z) (-0.38-0.42*exp(-0.5*z*10^-3))*1000 + ((-0.38-0.42*exp(-0.5*z*10^-3))*1000)*0.1; 
    %function for correct density contrast 
    density2=@(z) (-0.38-0.42*exp(-0.5*z*10^-3))*1000; 
    z_obs=0;   %height of observation point is in surface
    %Closed polygonic profile of the basin
    xx1=[x_obs x_obs(end) 0];
    yy1=[depth 0 0];
    %Gravity anomaly for fixed density model 
    gg=poly_gravityrho(x_obs,z_obs,xx1,yy1,density2,t_leg,c_leg);
    %saving temporarily all data
    x_obs=x_obs'; gg=gg';
    save x_obs_data.dat x_obs;
    save gravity_anomaly_data.dat gg;

 
    %horizontal location of observation points
    x_obs_data='x_obs_data.dat';     
    %corresponding gravity anomalies at observation points
    gravity_anomaly_data='gravity_anomaly_data.dat';   
    data=importdata(gravity_anomaly_data);
    
    %maximum depth bound
    max_depth=6000;
    %calling the function for inversion with incorrect density
    [depth_inverted,grav_inverted]=Fourier_PSO(x_obs_data,gravity_anomaly_data,density1,max_depth);
    %error norm
    err_depth=norm((depth-depth_inverted)/depth)*100;
    err_gravity=norm((data-grav_inverted')/data)*100;
    fprintf('For Model-2 fixed density the relative percentage error in depth profile is %f\n',err_depth)
    fprintf('For Model-2 fixed density the relative percentage error in gravity anomaly profile is %f\n',err_gravity)
    
    delete  x_obs_data.dat
    delete gravity_anomaly_data.dat
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
