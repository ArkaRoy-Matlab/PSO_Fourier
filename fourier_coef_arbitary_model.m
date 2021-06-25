% ***************************************************************
% *** Matlab code for an arbritary model for finding correlation with depth and gravity profile.
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
%Synthetic model for gravity field and Fourier transform
clear all
close all

%data for an arbritary synthetic basin 
x_obs=(importdata('x_obs_Model3.dat'))'; %observation points
depth= (importdata('depth_Model3.dat'))';%Depth profile
%Finding Gravity field of the basin for exponential density in kg/m^3
density=@(z)(-0.38-0.42*exp(-0.5*z*10^-3))*1000; 
%t and c are Legendre Gaussian quadrature points for numerical integration
    [t_leg,c_leg]=lgwt(10,0,1); 
    %plotting of x_obs vs. depth
    figure(1)
    subplot(2,1,2)
     xx_p=x_obs;
     xx_p=[0 xx_p x_obs(end)];
     depth_p=[0 depth 0];
     %patch plot for true model 
     zd_p=density(depth_p);
     patch(xx_p,depth_p,zd_p)
     colormap copper
     hold on
    plot(x_obs,depth,'k','linewidth',2)
    set(gca,'Ydir','reverse')
    xlabel('Distance (m)')
    ylabel('Depth (m)')
    title('Depth profile of a complex synthetic basin')
    box on
    c = colorbar('location','southoutside');
    c.Label.String = 'Density contrast (kg/m^3)';
    z_obs=0;   %Vertical position of observation point is in surface
    %Closed polygonic profile of the basin
    xx1=[x_obs x_obs(end) 0];
    yy1=[depth 0 0];
    %Gravity anomaly for fixed density model 
    zz1=poly_gravityrho(x_obs,z_obs,xx1,yy1,density,t_leg,c_leg);
    zz2=zz1;
    %Plotting the Gravity field anomaly 
    figure(1)
    subplot(2,1,1)
    plot(x_obs,zz2)
    xlabel('Observation points (m)')
    ylabel('Gravity anomaly (mGal)')
    title('Gravity anomaly for complex synthetic basin')
    %Fourier coefficients for Gravity anomaly
    TT=1:2*length(zz1);
    ll=(TT(end)-TT(1))/2;
    %Gravity data for Fourier series 
    grav_data=[zz1 fliplr(zz1)];
    %Depth data for Fourier series 
    depth_data=[depth fliplr(depth)];
    %1st 100 Fourier Coefficients 
    for j=1:100
      %Fourier coefficients for Gravity field 
      ss1=grav_data.*cos(j*pi*TT/ll);  
      aa_grv(j)=(1/ll)*trapz(TT,ss1);
      ss2=grav_data.*sin(j*pi*TT/ll);  
      bb_grv(j)=(1/ll)*trapz(TT,ss2);
      vv_grv(j+1)=sqrt((aa_grv(j))^2+(bb_grv(j))^2);
      
      %Fourier coefficients for Depth Profile  
      dd1=depth_data.*cos(j*pi*TT/ll);  
      aa_dep(j)=(1/ll)*trapz(TT,dd1);
      dd2=depth_data.*sin(j*pi*TT/ll);  
      bb_dep(j)=(1/ll)*trapz(TT,dd2);
      vv_dep(j+1)=sqrt((aa_grv(j))^2+(bb_grv(j))^2);
      wn_num(j+1)=j*pi/ll;
    end
    %Normalized amplitude of Power spectrum for Gravity anomaly
    vv_grv(1)=abs((1/ll)*trapz(TT,grav_data));
    vv_grv=vv_grv./max(vv_grv);
    %Normalized amplitude of Power spectrum for Depth Profile 
    vv_dep(1)=abs((1/ll)*trapz(TT,depth_data));
    vv_dep=vv_dep./max(vv_dep);
    wn_num(1)=0;
    %Plotting the power spectrum for Gravity anomaly
    figure(2)
    subplot(2,1,1)
    semilogy(wn_num(1:51),vv_grv(1:51))
    hold on
    semilogy(wn_num(1:51),10^-2.*ones(51,1),'k--')
    title('Fourier power spectrum for gravity anomaly')
    xlabel('Wavenumbers')
    ylabel('Amplitude')
    
    %Plotting the power spectrum for Derpth Profile 
    subplot(2,1,2)
    semilogy(wn_num(1:51),vv_dep(1:51))
    
    title('Fourier power spectrum for depth profile')
    xlabel('Wavenumbers')
    ylabel('Amplitude')
    
%finding correlation coefficients between power spectrum for Gravity
%anomaly and Depth profile 
cc=corrcoef(vv_grv,vv_dep);
fprintf('Correlation coefficients of power spectrum for gravity and depth profile is %f\n',cc(1,2))

%selecting number of Fourier Coefficients less than 10^-2 order amplitude   
n1=find(vv_grv>=10^-2,1,'last');
%printing the number of parameters 
fprintf('Number of parameters for Fourier coefficients of Depth profile = %d\n',n1)
data_x=x_obs; %data for horizontal locations
data=zz1;     %data for gravity anomaly  
%correlation coefficients for Fourier spectrum of gravity anomaly and depth profile
cc=corrcoef(vv_grv,vv_dep);
fprintf('Correlation coefficients of power spectrum for gravity and depth profile is %f\n',cc(1,2))

%% creating Fourier Transformation matrix for multiplication 
T=1:2*length(data_x);
l=(T(end)-T(1))/2;
for i=1:length(T)
    jj1=0;jj2=0;
    for j=1:2*n1+1
        if j==1
            A_mat(i,j)=(1/2);
        elseif j>1 && j<=n1+1
            jj1=jj1+1;
            A_mat(i,j)=cos(jj1*pi*T(i)/l);
        else
            jj2=jj2+1;
            A_mat(i,j)=sin(jj2*pi*T(i)/l);
        end
    end
end
%coefficients for gravity anomaly reconstruction
x_var(1)=(1/ll)*trapz(TT,grav_data);
for i=1:n1
    x_var(i+1,1)=aa_grv(i);
end

for i=1:n1
    x_var(n1+i+1,1)=bb_grv(i);
end
dtaa_grav=A_mat*x_var;
figure(3)
subplot(2,1,1)
plot(x_obs,zz2,'linewidth',1.25)
hold on
plot(x_obs,dtaa_grav(1:length(x_obs)),'--','linewidth',1.25)
     xlabel('Observation points (m)')
     ylabel('Gravity anomaly (mGal)')
     legend('Synthetic Data','Reconstructed Data','location','southeast')
%coefficients for depth profile reconstruction
x_var(1)=(1/ll)*trapz(TT,depth_data);
for i=1:n1
    x_var(i+1,1)=aa_dep(i);
end

for i=1:n1
    x_var(n1+i+1,1)=bb_dep(i);
end
dtaa_depth=A_mat*x_var;

     subplot(2,1,2)
     %patch plot for true model 
     zd_p=density(depth_p);
     patch(xx_p,depth_p,zd_p)
     colormap copper
     hold on
    plot(x_obs,dtaa_depth(1:length(x_obs)),'r.-','linewidth',1)
    set(gca,'Ydir','reverse')
    xlabel('Distance (m)')
    ylabel('Depth (m)')
    legend('Synthetic Model','Reconstructed Model','location','southeast')
    box on
    c = colorbar('location','southoutside');
    c.Label.String = 'Density contrast (kg/m^3)';
    
    %selecting number of Fourier Coefficients less than 10^-5 order amplitude   
n1=100;
%printing the number of parameters 
fprintf('Number of parameters for Fourier coefficients of Depth profile = %d\n',n1)
data_x=x_obs; %data for horizontal locations
data=zz1;     %data for gravity anomaly  
%correlation coefficients for Fourier spectrum of gravity anomaly and depth profile
cc=corrcoef(vv_grv,vv_dep);
fprintf('Correlation coefficients of power spectrum for gravity and depth profile is %f\n',cc(1,2))

%% creating Fourier Transformation matrix for multiplication 
T=1:2*length(data_x);
l=(T(end)-T(1))/2;
for i=1:length(T)
    jj1=0;jj2=0;
    for j=1:2*n1+1
        if j==1
            A_mat(i,j)=(1/2);
        elseif j>1 && j<=n1+1
            jj1=jj1+1;
            A_mat(i,j)=cos(jj1*pi*T(i)/l);
        else
            jj2=jj2+1;
            A_mat(i,j)=sin(jj2*pi*T(i)/l);
        end
    end
end
%coefficients for gravity anomaly reconstruction
x_var(1)=(1/ll)*trapz(TT,grav_data);
for i=1:n1
    x_var(i+1,1)=aa_grv(i);
end

for i=1:n1
    x_var(n1+i+1,1)=bb_grv(i);
end
dtaa_grav=A_mat*x_var;
figure(4)
subplot(2,1,1)
plot(x_obs,zz2,'linewidth',1.25)
hold on
plot(x_obs,dtaa_grav(1:length(x_obs)),'--','linewidth',1.25)
     xlabel('Observation points (m)')
     ylabel('Gravity anomaly (mGal)')
     legend('Synthetic Data','Reconstructed Data','location','southeast')
%coefficients for depth profile reconstruction
x_var(1)=(1/ll)*trapz(TT,depth_data);
for i=1:n1
    x_var(i+1,1)=aa_dep(i);
end

for i=1:n1
    x_var(n1+i+1,1)=bb_dep(i);
end
dtaa_depth=A_mat*x_var;
     subplot(2,1,2)
     %patch plot for true model 
     zd_p=density(depth_p);
     patch(xx_p,depth_p,zd_p)
     colormap copper
     hold on
    plot(x_obs,dtaa_depth(1:length(x_obs)),'r.-','linewidth',1)
    set(gca,'Ydir','reverse')
    xlabel('Distance (m)')
    ylabel('Depth (m)')
    legend('Synthetic Model','Reconstructed Model','location','southeast')
    box on
    c = colorbar('location','southoutside');
    c.Label.String = 'Density contrast (kg/m^3)';
    