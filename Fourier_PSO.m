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
%%Matlab function for inversion of sedimentary basin using PSO and adaptive
%%Fourier coefficients
function [depth_inverted,grav_inverted]=Fourier_PSO(x_obs_data,gravity_anomaly_data,density)
    
    %input:
    %    1. x_obs           =  data file for horizontal location of
    %                          observation points (single collumn data file in ascii format in m.)
    %    2. gravity_anomaly =  corresponding gravity anomaly data file at the
    %                          observation points (single collumn data file in ascii format in mGal)
    %    3. density         =  density contrast as a function of z in
    %                          kg/m^3
    
    %output:
    %   1.depth_inverted    = inverted depth profile of sedimentary basin
    %                         at the locations of horizontal observation
    %                         points in m.
    %   1.grav_inverted    =  inverted gravity anomaly of sedimentary basin
    %                         at the locations of horizontal observation
    %                         points in mGal.

        %% Gravity anomaly and its inversion for any sedimentary basin
        
        %import data for observation points and gravity anomaly
        x_obs=(importdata(x_obs_data))';
        zz1=(importdata(gravity_anomaly_data))';
        
        sz1=size(x_obs);
        
        if sz1(1,1)~=1
            x_obs=x_obs';
        end
        
        sz2=size(zz1);
        if sz2(1,1)~=1
            zz1=zz1';
        end
        
        %t and c are Legendre Gaussian quadrature points for numerical integration
        [t_leg,c_leg]=lgwt(10,0,1); 
        z_obs=0;   %height of observation point is in surface

        %Fourier coefficients for Gravity anomaly
        TT=1:2*length(zz1);
        ll=(TT(end)-TT(1))/2;
        %Gravity data for Fourier series 
        grav_data=[zz1 fliplr(zz1)];

        %1st 100 Fourier Coefficients 
        for j=1:100
          %coefficients for Gravity field 
          ss1=grav_data.*cos(j*pi*TT/ll);  
          aa_grv(j)=(1/ll)*trapz(TT,ss1);
          ss2=grav_data.*sin(j*pi*TT/ll);  
          bb_grv(j)=(1/ll)*trapz(TT,ss2);
          vv_grv(j+1)=sqrt((aa_grv(j))^2+(bb_grv(j))^2);
          wn_num(j+1)=j*pi/ll;
        end
        %normalization of Fourier Coefficients of Gravity anomaly 
        vv_grv(1)=abs((1/ll)*trapz(TT,grav_data));
        vv_grv=vv_grv./max(vv_grv);

        wn_num(1)=0;
        figure(1)
        %Plotting the Fourier Spectrum of Gravity anomaly
        semilogy(wn_num(1:51),vv_grv(1:51))
        hold on
        semilogy(wn_num(1:51),10^-2.*ones(51,1),'k--')
        title('Fourier power spectrum for gravity anomaly')
        xlabel('Wavenumbers')
        ylabel('Amplitude')


    %selecting number of Fourier Coefficients    
    n=find(vv_grv>=10^-2, 1, 'last' )-1;
    %printing the number of parameters 
    fprintf('Number of parameters for Fourier coefficients of Depth profile = %d\n',n)
    data_x=x_obs; %data for horizontal locations
    data=zz1;     %data for gravity anomaly  

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
    %%%%    Optimization of Sayula using PSO     %%%%
    c1=1.4; c2=1.7; 
                %Cost function with constraints for optimization
                CostFunction =@(x) myCostFunction(x,x_obs,data,A_mat,t_leg,c_leg,density)+1000*(constrained1(x,A_mat)+constrained2(x,A_mat)); 
                nVar=(2*n+1);          %Number of Unknown Variables
                MaxIt = 300;           %Maximum number of iterations
                nPoP =  40;            %Population size or swarm size
                %Loop for 2 independent run
                for ii=1:2
                    %Calling the function for optimization using PSO
                    [bst_var, best_cost,iter_count,error_energy] = WIPSO(CostFunction,nVar,MaxIt,nPoP,c1,c2);
                    %Saving all best cost and parameters for each independent run 
                    all_best(ii)=best_cost;
                    all_var(:,ii)=bst_var;
                    err_enrgy(:,ii)=error_energy';
                end
                %finding best model 
                [val,id]=min(all_best);
                %Parameters for best Model 
                best_var=(squeeze(all_var(:,id)))';  
                %error energy for best model 
                best_error=(squeeze(err_enrgy(:,id)))';
                

    %% Plotting    
    %%%%    Plotting the results for Optimization model    %%%%
           x=best_var'; % Fourier coefficients for optimized model
           xx=x_obs;    % observation points
         %gravity field at 40 discrete data points  
         xx1=linspace(min(xx),max(xx),40);
         %depth profile from optimized Fourier coefficients
         yy_eval=(A_mat*x)';
         yy=yy_eval(1:length(yy_eval)/2);
         yy1=spline(xx,yy,xx1);
         %close polygonal form of depth profile 
         x1(1:length(xx1)+2)=[xx1 max(xx1) min(xx1)];
         y1(1:length(yy1)+2)=[yy1 0 0];
         %gravity field for optimized depth profile 
         zz_eval=poly_gravityrho(xx,0,x1,y1,density,t_leg,c_leg);
         depth_inverted=yy;
         grav_inverted=zz_eval;
         %true model 
         figure(2)
         subplot(2,1,2)
         hold on
         xx_p=[xx1(1) xx1 xx1(end)];
         depth_p=[0 yy1 0];
         %patch plot for true model 
         zd_p=density(depth_p);
         patch(xx_p,depth_p,zd_p)
         colormap copper
         %optimized model 
         plot(xx1,yy1,'linewidth',1.25,'color','r')
         set(gca,'Ydir','reverse')
         c = colorbar('location','southoutside');
         c.Label.String = 'Density contrast (kg/m^3)';
         xlabel('Distance (m)')
         ylabel('Depth (m)')
         title('Depth profile')
         box on
         subplot(2,1,1)
         hold on
         %Observed 
         plot(xx,data,'-o','linewidth',1.25)
         %Inverted 
         plot(xx,zz_eval,'linewidth',1.25,'color','r')
         %Axis lebeling 
         xlabel('Observation points(m)')
         ylabel('Gravity anomaly (mGal)')
         title('Gravity anomaly ')
         legend('Observed Data','Inverted Data','location','southeast')
         box on

         %RMSE for gravity  
         N_g=length(data); %Number of Observation points 
         RMSE_g=sqrt((sum((data-zz_eval).^2))/N_g)/(max(data(:))-min(data(:)));

     %Printing the RMSE error for depth and gravity profile 
         fprintf('RMSE in gravity field=%f\n',RMSE_g)



    %%%%    Objective functions and Constraints   %%%%
    function val=myCostFunction(x,x_obs,data,A_mat,t_leg,c_leg,density)
         %%inputs are 
         %x= parameters for PSO algorithm
         %data= observed gravity field 
         %A_mat= Matrix with Fourier coefficients 
         %t_leg and c_leg are  Legendre Gaussian quadrature points for numerical integration
         % subroutine for t_leg and c_leg evaluation is given in lgwt.m file 
        %observation points having 100 linearly spaced data points  
         xx=linspace(x_obs(1),x_obs(end),101);
        %40 linearly spaced data points in closed polygonal form
         xx1=linspace(min(xx),max(xx),40);
         yy1=(A_mat*x')';
         yy=yy1(1:length(yy1)/2);
         yy1=spline(xx,yy,xx1);
         %close polygonal form of depth profile 
         x1(1:length(xx1)+2)=[xx1 xx1(end) xx1(1)];
         y1(1:length(yy1)+2)=[yy1 0 0];
         %gravity field for given depth profile 
         zz1=poly_gravityrho(xx,0,x1,y1,density,t_leg,c_leg);
         %zz1=poly_gravityrho(xx,0,x1,y1,density,t_leg,c_leg);
         %misfit functional for observed and inverted gravity anomaly
         val=norm(data-zz1);
    end

    %Constraints for upper bound 
    function val=constrained1(x,A_mat)
        %%inputs are 
         %x= parameters for DE algorithm
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
    %Constraints for lower bound 
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

end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  