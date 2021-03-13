% ***************************************************************
% *** Matlab function for finding gravity anomaly of depth varying density distribution
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

function grav=poly_gravityrho(x_obs,z_obs,x,z,rho,t,c)

%poly_gravityrho function calculates z component of gravity field for any polygon 
%shape 2d body having depth varying density contrast. This program based on line
%integral in anticlockwise direction using Gauss Legendre quadrature
%integral formula. For more detail go through Zhou 2008. 

%Inputs 
%	x_obs = a vector containg observation points in x direction.
%	z_obs = level at which we are calculating gravity field, positive along
%			downward direction
%	x 	  = the x coordinates of polygon body in counterclockwise direction. 
%	z 	  = the z coordinates of polygon body in counterclockwise direction. 
%	roh   = the depth varying density contrast as a function of z.
%   t     = Legendre Gaussian quadrature integral points.
%   c     = Legendre Gaussian quadrature nodes.
%Outputs 
%	grav= gravity field for given inputs in mGal Unit 
% Always keep in mind x & z should always be taken in counter clockwise
% direction, otherwise sign convention will create problem while running. 

    n_poly=length(x); %length of the polygon
    x(length(x)+1)=x(1);% end point should be 1st point to close the integral
    z(length(z)+1)=z(1);% end point should be 1st point to close the integral
    G=6.67408*10^-11;% Gravitational constant in S.I
    for i=1:length(x_obs) % Loop for all observation points. 
        for j=1:n_poly % Loop for line integral over all sides of polygon 

            % Refer to Zhou 2008 paper for below steps, basically line
            % integral procedures. 
            
            ax1=(x(j).*(1-t)+x(j+1).*t-x_obs(i));
            ax2=(z(j).*(1-t)+z(j+1).*t-z_obs);
            rr=rho(ax2);
            ax=-2.*rr.*G.*(atan(ax1./ax2)).*(z(j+1)-z(j));
            value(j) = sum(c.*ax);
            
        end
        grav(i)=10^5*sum(value(:)); % Combined gravity field for all points. 
    end
end

