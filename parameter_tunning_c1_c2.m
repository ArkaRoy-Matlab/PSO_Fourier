% ***************************************************************
% *** Matlab function for parameter tuning of cognitive (c1) and social (c2) component of pso
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
%parameter tuning of cognitive (c1) and social (c2) component of pso
clear all
close all

%function for Synthetic Example of depth of the basin
z_peaks=peaks(100);
y_peaks=((z_peaks(41,:))+1.5)*200;

%t and c are Legendre Gaussian quadrature points for numerical integration
    [t_leg,c_leg]=lgwt(10,0,1); 
    %synthetic depth and observation points
    x_obs=linspace(0,5000,100);
    depth=y_peaks;
    %Finding Gravity field of the basin for density -800 kg/m^3
    density=-800; 
    z_obs=0;   %height of observation point is in surface
    %Closed polygonic profile of the basin
    xx1=[x_obs x_obs(end) 0];
    yy1=[depth 0 0];
    %Gravity anomaly for fixed density model 
    zz1=poly_gravity(x_obs,z_obs,xx1,yy1,density,t_leg,c_leg);
    zz2=zz1;
    %Fourier coefficients for Gravity anomaly
    TT=1:2*length(zz1);
    ll=(TT(end)-TT(1))/2;
    %Gravity data for Fourier series 
    grav_data=[zz1 fliplr(zz1)];
    %Depth data for Fourier series 
    depth_data=[depth fliplr(depth)];
    %1st 100 Fourier Coefficients 
    for j=1:100
      %coefficients for Gravity field 
      ss1=grav_data.*cos(j*pi*TT/ll);  
      aa_grv(j)=(1/ll)*trapz(TT,ss1);
      ss2=grav_data.*sin(j*pi*TT/ll);  
      bb_grv(j)=(1/ll)*trapz(TT,ss2);
      vv_grv(j+1)=sqrt((aa_grv(j))^2+(bb_grv(j))^2);
      
      %coefficients for Depth Profile  
      dd1=depth_data.*cos(j*pi*TT/ll);  
      aa_dep(j)=(1/ll)*trapz(TT,dd1);
      dd2=grav_data.*sin(j*pi*TT/ll);  
      bb_dep(j)=(1/ll)*trapz(TT,dd2);
      vv_dep(j+1)=sqrt((aa_grv(j))^2+(bb_grv(j))^2);
      
    end
    %normalization of Fourier Coefficients of Gravity anomaly 
    vv_grv(1)=abs((1/ll)*trapz(TT,grav_data));
    vv_grv=vv_grv./max(vv_grv);
    %normalization of Fourier Coefficients of Depth Profile 
    vv_dep(1)=abs((1/ll)*trapz(TT,depth_data));
    vv_dep=vv_dep./max(vv_dep);
%selecting number of Fourier Coefficients for optimization   
n=find(vv_grv>=10^-2, 1, 'last' )-1;
data_x=x_obs;
data=zz1;

 %% creating Fourier Transformation matrix for multiplication 
T=1:2*length(data_x);
l=(T(end)-T(1))/2;
for i=1:length(T)
    jj1=0;jj2=0;
    for j=1:2*n+1
        if j==1
            A_mat(i,j)=(1/2);
        elseif j>1 && j<=n+1
            jj1=jj1+1;
            A_mat(i,j)=cos(jj1*pi*T(i)/l);
        else
            jj2=jj2+1;
            A_mat(i,j)=sin(jj2*pi*T(i)/l);
        end
    end
end
     
%% Problem Definition
%all parameters for c1 and c2 for tunning 
cc1=1:0.1:2; cc2=1:0.1:2; 
%loop for varyng c1 and c2
for ic1=1:length(cc1)
    for ic2=1:length(cc2)
            %Cost function with constraints for optimization
            CostFunction =@(x) myCostFunction(x,data,A_mat,t_leg,c_leg)+1000*(constrained1(x,A_mat)+constrained2(x,A_mat)); 
            nVar=(2*n+1);          %Number of Unknown Variables
            MaxIt = 300;           %Maximum number of iterations
            nPoP =  40;            %Population size or swarm size
            c1=cc1(ic1); c2=cc2(ic2);
            %Loop for 10 independent run
            for ii=1:10
                tic
                [best_var, best_cost,iter_count] = WIPSO(CostFunction,nVar,MaxIt,nPoP,c1,c2);
                time=toc;
                all_best(ii,1)=best_cost;
                all_best(ii,2)=time;
                all_best(ii,3)=iter_count;
            end
            %mean of all best cost, computation time and iteration count
            mean_beast=mean(all_best);
            costt(ic1,ic2)=mean_beast(1);
            timee(ic1,ic2)=mean_beast(2);
            iterr(ic1,ic2)=mean_beast(3);
            %printing results for all c1 and c2
            fprintf('For c1=%1.1f and c2=%1.1f, best cost =%f, best time==%2.2f and number of iteration=%d\n',c1,c2,mean_beast(1),mean_beast(2),mean_beast(3))
    end
end
%%    
%%%%    Objective functions and Constraints   %%%%
function val=myCostFunction(x,data,A_mat,t_leg,c_leg)
     %%inputs are 
     %x= parameters for PSO algorithm
     %data= observed gravity field 
     %A_mat= Matrix with Fourier coefficients 
     %t_leg and c_leg are  Legendre Gaussian quadrature points for numerical integration
     % subroutine for t_leg and c_leg evaluation is given in lgwt.m file 
    %observation points having 100 linearly spaced data points  
     xx=linspace(0,5*10^3,100);
    %40 linearly spaced data points in closed polygonal form
     xx1=linspace(min(xx),max(xx),40);
     yy1=(A_mat*x')';
     yy=yy1(1:length(yy1)/2);
     yy1=spline(xx,yy,xx1);
     %close polygonal form of depth profile 
     x1(1:length(xx1)+2)=[xx1 5000 0];
     y1(1:length(yy1)+2)=[yy1 0 0]; 
     %gravity field for depth profile 
     zz1=poly_gravity(xx,0,x1,y1,-800,t_leg,c_leg);
     %misfit functional for observed and inverted gravity anomaly
     val=norm(data-zz1);
end

%Constraints for upper bound 
function val=constrained1(x,A_mat)
    %%inputs are 
     %x= parameters for PSO algorithm
     %A_mat= Matrix with Fourier coefficients
     %depth profile from Fourier domain to time domain
     yy1=A_mat*x';
     yy=yy1(1:length(yy1)/2);
     %minimum of depth 
     m_yy=min(yy);
     %penalty barrier method for minimum bound 
     gg=(-m_yy+0);
     val=(max(0,gg))^2;
     
end

function val=constrained2(x,A_mat)
      %%inputs are 
     %x= parameters for DE algorithm
     %A_mat= Matrix with Fourier coefficients
     %depth profile from Fourier domain to time domain
     yy1=A_mat*x';
     yy=yy1(1:length(yy1)/2);
     %maximum of depth 
     m_yy=max(yy);
     %penalty barrier method for minimum bound 
     gg=(m_yy-6000);
     val=(max(0,gg))^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  