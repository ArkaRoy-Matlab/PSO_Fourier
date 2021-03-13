%Example code for getting depth profile from observed gravity anomalies
%with incorrect density contrast for Sayula basin

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
%%Matlab code for inversion of Sayula sedimentary basin with incorrect density
%%contrast
%%Effect of incorrect density on this algorithm
clear all
close all
%location of true depth for Sayula basin
    x_obs=(importdata('x_obs_sayula.dat'))';
    depth=(importdata('depth_sayula.dat'))';
    
    %corresponding gravity anomalies at observation points
    %horizontal location of observation points
    x_obs_data='x_obs_sayula.dat';     
    gravity_anomaly_data='gravity_anomaly_sayula.dat';
    
    data=importdata(gravity_anomaly_data);
    %function for incorrect density contrast 
    density1=@(z) (-0.4692*exp(-4.078*10^-4*z)).*1000;
    %function for correct density contrast 
    density2=@(z) (-0.8+ 0.7147.*(z*10^-3) -0.229.*(z*10^-3).^2)*1000;
    %calling the function for inversion with incorrect density
    [depth_inverted,grav_inverted]=Fourier_PSO(x_obs_data,gravity_anomaly_data,density1);
    %close all figures
    close all
    %patch plot for true model 
     figure(1)
     subplot(2,1,2)
     hold on
     xx_p=[x_obs(1) x_obs x_obs(end)];
     depth_p=[0 depth 0];
     %patch plot for true model 
     zd_p=density2(depth_p);
     patch(xx_p,depth_p,zd_p)
     colormap copper
     %optimized model 
     plot(x_obs,depth_inverted,'linewidth',2.25,'color','r')
     set(gca,'Ydir','reverse')
     c = colorbar('location','southoutside');
     c.Label.String = 'Density contrast (kg/m^3)';
     xlabel('Distance (m)')
     ylabel('Depth (m)')
     legend('True depth','Optimized depth','location','southeast')
     title('Depth profile of Sayula basin with incorrect density')
     
     subplot(2,1,1)
     hold on
     %Observed 
     plot(x_obs,data,'-o','linewidth',1.25)
     %Inverted 
     plot(x_obs,grav_inverted,'linewidth',1.25,'color','r')
     %Axis lebeling 
     xlabel('Observation points (m)')
     ylabel('Gravity anomaly (mGal)')
     title('Gravity anomaly of Sayula basin with incorrect density')
     legend('Observed Data','Inverted Data','location','southeast')
     box on
     