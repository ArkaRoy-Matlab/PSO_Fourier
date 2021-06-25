% ***************************************************************
% *** Matlab code for error energy plot   
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

%Error Energy plot 
%%Matlab code for Error energy plot
clear all
close all

%importing all data for error energy plot
%noise free data
y1=importdata('error_energy_model1_fixed_wo');
y2=importdata('error_energy_model1_varying_wo');
y3=importdata('error_energy_model2_fixed_wo');
y4=importdata('error_energy_model2_varying_wo');

%noisy data
z1=importdata('error_energy_model1_fixed_w');
z2=importdata('error_energy_model1_varying_w');
z3=importdata('error_energy_model2_fixed_w');
z4=importdata('error_energy_model2_varying_w');

%plotting the error energies for noise free data
figure(1)
semilogy(y1,'linewidth',2)
hold on
semilogy(y2,'linewidth',2)
semilogy(y3,'linewidth',2)
semilogy(y4,'linewidth',2)
ylabel('Error Energy (mGal^2)')
xlabel('Generation')
%title('Error energy plot for noise free case')
legend('Model1 fixed density','Model1 varying density','Model2 fixed density','Model2 varying density')
%plotting the error energies for noisy data

figure(2)
semilogy(z1,'linewidth',2)
hold on
semilogy(z2,'linewidth',2)
semilogy(z3,'linewidth',2)
semilogy(z4,'linewidth',2)
ylabel('Error Energy (mGal^2)')
xlabel('Generation')
%title('Error energy plot for noisy case')
legend('Model1 fixed density','Model1 varying density','Model2 fixed density','Model2 varying density')