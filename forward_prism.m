% ***************************************************************
% *** Matlab function for prismatic model with fixed density 
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

%Code for Comparative study of Prismatic model with present scheme for Model-1 with
%depth varying density
clear all
close all
%%function for Synthetic Example of depth of the basin
z_peaks=peaks(100);
y_peaks=((z_peaks(71,:))+0.5)*200;

    %synthetic depth and observation points
    x_obs=linspace(0,5000,100);
    depth=y_peaks;
    
    %t and c are Legendre Gaussian quadrature points for numerical integration
    [t_leg,c_leg]=lgwt(10,0,1); 
    
    %plotting of x_obs vs. depth
    figure(1)
    plot(x_obs,depth,'r','linewidth',2)
    set(gca,'Ydir','reverse')
    xlabel('Distance (m)')
    ylabel('Depth (m)')
    %title('True Model and Prismatic Approximation of depth profile for synthetic basin (Model 1)')

    %Finding Gravity field of the basin for density rho(z)= (-0.55-2.5*10^-3.*z).*1000 kg/m^3
    density=-800; 
    z_obs=0; %height of observation point is in surface
    %polygonic profile of the basin
    xx1=[x_obs x_obs(end) 0];
    yy1=[depth 0 0];
    zz1=poly_gravity(x_obs,z_obs,xx1,yy1,density,t_leg,c_leg);
    

%prismatic model
%creating synthetic depth profile
    nn=50; %number of prism
    xx1=linspace(0,5000,nn);
    yy1=spline(x_obs,depth,xx1);
%loop for plotting all prisms along the depth profile
grav=0;
for ii=1:length(xx1)-1
    %vertices of each prism
    x1=[xx1(ii) xx1(ii) xx1(ii+1) xx1(ii+1)];
    y1=[0 yy1(ii+1) yy1(ii+1) 0];
    %plotting each prism
    figure(1)
    hold on
    pgon=polyshape(x1,y1);
    plot(pgon,'FaceColor','blue','FaceAlpha',0.1)
    %gravity anomaly for each prism at each observation points
    zz11=poly_gravityrho(x_obs,z_obs,x1,y1,@(z) (-0.55-2.5*10^-3.*z).*1000,t_leg,c_leg);
    %sum of gravity anomaly of each prism
    grav=grav+zz11;
end
%upper and lower limit for plotting window
xlim([0 5000])
ylim([0 1600])
legend('True Model','Stacked prisms','location','best')

%initialization of loop
nn=5; %number of prism
RMSE_g=10;
%loop for prismatic model
while RMSE_g>=0.50 %(desired stopping criterion same as best RMSE from SPoDEA)
    nn=nn+1; %Increament in number of prisms 
    xx1=linspace(0,5000,nn);
    yy1=spline(x_obs,depth,xx1);
    %gravity field due to nn prisms
    grav=0;
    for ii=1:length(xx1)-1
        x1=[xx1(ii) xx1(ii) xx1(ii+1) xx1(ii+1)];
        y1=[0 yy1(ii+1) yy1(ii+1) 0];
        zz11=poly_gravity(x_obs,z_obs,x1,y1,-800,t_leg,c_leg);
        grav=grav+zz11;
    end
    %RMSE error 
    N_g=length(grav);
    RMSE_g=(sqrt((sum((grav-zz1).^2))/N_g)/(max(grav(:))-min(grav(:))))*100;
    %printing RMSE for 10 prism
    if nn==8
        fprintf('For n=%d RMSE in gravity field=%f.\n',nn,RMSE_g)
    end
end
%RMSE of Prismatic model having desired accuracy 
fprintf('For n=%d RMSE in gravity field=%f.\n',nn,RMSE_g)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
%%function for Synthetic Example of depth of the basin
z_peaks=peaks(200);
yy=((z_peaks(51*2,:)))*450+abs(((z_peaks(9*2,:)))*450);
y_peaks=yy(74:200);


    %synthetic depth and observation points
    x_obs=linspace(0,5000,length(y_peaks));
    depth=y_peaks;
    
    %t and c are Legendre Gaussian quadrature points for numerical integration
    [t_leg,c_leg]=lgwt(10,0,1); 
    
    %plotting of x_obs vs. depth
    figure(2)
    plot(x_obs,depth,'r','linewidth',2)
    set(gca,'Ydir','reverse')
    xlabel('Distance (m)')
    ylabel('Depth (m)')
    %title('True Model and Prismatic Approximation of depth profile for synthetic basin (Model 1)')

    %Finding Gravity field of the basin for density rho(z)= (-0.55-2.5*10^-3.*z).*1000 kg/m^3
    density=-800; 
    z_obs=0; %height of observation point is in surface
    %polygonic profile of the basin
    xx1=[x_obs x_obs(end) 0];
    yy1=[depth 0 0];
    zz1=poly_gravity(x_obs,z_obs,xx1,yy1,density,t_leg,c_leg);
    

%prismatic model
%creating synthetic depth profile
    nn=50; %number of prism
    xx1=linspace(0,5000,nn);
    yy1=spline(x_obs,depth,xx1);
%loop for plotting all prisms along the depth profile
grav=0;
for ii=1:length(xx1)-1
    %vertices of each prism
    x1=[xx1(ii) xx1(ii) xx1(ii+1) xx1(ii+1)];
    y1=[0 yy1(ii+1) yy1(ii+1) 0];
    %plotting each prism
    figure(2)
    hold on
    pgon=polyshape(x1,y1);
    plot(pgon,'FaceColor','blue','FaceAlpha',0.1)
    %gravity anomaly for each prism at each observation points
    zz11=poly_gravityrho(x_obs,z_obs,x1,y1,@(z) (-0.55-2.5*10^-3.*z).*1000,t_leg,c_leg);
    %sum of gravity anomaly of each prism
    grav=grav+zz11;
end
%upper and lower limit for plotting window
xlim([0 5000])
ylim([0 1900])
legend('True Model','Stacked prisms','location','best')

%initialization of loop
nn=5; %number of prism
RMSE_g=10;
%loop for prismatic model
while RMSE_g>=0.50 %(desired stopping criterion same as best RMSE from SPoDEA)
    nn=nn+1; %Increament in number of prisms 
    xx1=linspace(0,5000,nn);
    yy1=spline(x_obs,depth,xx1);
    %gravity field due to nn prisms
    grav=0;
    for ii=1:length(xx1)-1
        x1=[xx1(ii) xx1(ii) xx1(ii+1) xx1(ii+1)];
        y1=[0 yy1(ii+1) yy1(ii+1) 0];
        zz11=poly_gravity(x_obs,z_obs,x1,y1,-800,t_leg,c_leg);
        grav=grav+zz11;
    end
    %RMSE error 
    N_g=length(grav);
    RMSE_g=(sqrt((sum((grav-zz1).^2))/N_g)/(max(grav(:))-min(grav(:))))*100;
    %printing RMSE for 10 prism
    if nn==8
        fprintf('For n=%d RMSE in gravity field=%f.\n',nn,RMSE_g)
    end
end
%RMSE of Prismatic model having desired accuracy 
fprintf('For n=%d RMSE in gravity field=%f.\n',nn,RMSE_g)

