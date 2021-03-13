% ***************************************************************
% *** Help file for running all codes
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
This is a help file for a description of all Data, Source Code, and Subroutine used for the implementation of our present paper 
'Gravity inversion of basement relief using Particle Swarm Optimization by automated parameter selection of Fourier coefficients.'  

(Copy all set of files including data in one folder and run)

	1. Data Files
		a. x_obs_Model3
		b. depth_Model3
		c. gravity_anomaly_sayula.dat
		d. gravity_anomaly_Zhou2013_3a.dat
		e. x_obs_sayula.dat
		f. x_obs_Zhou2013_3a.dat
		g. depth_sayula.dat
		h. depth_Zhou2013_3a.dat
		i. error_energy_mode1_fixed_w
		j. error_energy_mode1_fixed_wo
		k. error_energy_mode1_varying_w
		l. error_energy_mode1_varying_wo
		m. error_energy_mode2_fixed_w
		n. error_energy_mode2_fixed_wo
		o. error_energy_mode2_varying_w
		p. error_energy_mode2_varying_wo
		
	File (a) and (b) are the observation points and depth profile for the arbitary modle shown in figure 1 of the manuscript, file (c) and (d) are the data for observed gravity anomaly of Sayula Basin and Godavari Basin respectively.
	File (g) and (h) are the estimated depth profile by Garcia-Abdeslem, 2003, and Zhou, 2013 from the gravity anomaly of Sayula Basin and Godavari Basin. 
File (c) and (d) are used here for inversion using PSO and file (g) and (h) are used here for verification. File (e) and (f) are the location of observation points for both profiles. File (g), (h), (i) and (j) are the data file for error energy plot of Model 1 and Model 2 
with out noisy cases and (k), (l), (m) and (n) are the data file for error energy plot of noisy data case.

	2. Subroutines
		a. lgwt.m
		b. poly_gravity.m
		c. poly_gravityrho.m
		d. WIPSO.m
		e. Fourier_PSO.m
		f. pca_reduction
		
	a. lgwt.m - This script is for computing definite integrals using Legendre-Gauss 
 Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval [a,b] with truncation order N. Suppose you have a continuous function f(x) which is defined on [a,b]
which you can evaluate at any x in [a,b]. Simply evaluate it at all of the values contained in the x vector to obtain a vector f. Then compute the definite integral using sum(f.*w);

	This code is written by Greg von Winckel - 02/25/2004. Here we have used it for our calculation and cited in main manuscript. 

	b. poly_gravity.m - poly_gravity function calculates z component of gravity field for any polygon shape 2d body having finite density contrast. This program based on line
integral in anticlockwise direction using Gauss Legendre quadrature integral formula. For more detail go through Zhou 2008. Here we have used it for calculation of gravity field for forward model as well as for inversion. 
	
	c. poly_gravityrho.m - poly_gravityrho function calculates z component of gravity field for any polygon shape 2d body having depth varying density contrast. This program based on line integral in anticlockwise direction using Gauss Legendre quadrature
%integral formula. For more detail go through Zhou 2008. It is same as poly_gravity function but for depth varying density contrast. 

	d. WIPSO.m - WIPSO calculates the optimized parameters (best_var) for a given objective function (CostFunction) using Particle Swarm Optimization.

	e. Fourier_PSO.m - Matlab function for inversion of sedimentary basin using PSO 
and adaptive Fourier coefficients that can be applicable for any real sedimentary basin inversion for given input of gravity anomaly and density distributions. 
	
	f. pca_reduction.m -Matlab function for principal component analysis for reducing data dimension
	
	3. Source Codes
		a. fourier_coef_arbitary_model.m
		b. final_gravity_model1_fixed.m
		c. final_gravity_model1_varying.m
		d. final_gravity_model2_fixed.m
		e. final_gravity_model2_varying.m
		f. godavari_basin.m
		g. sayula_basin.m
		h. error_energy_plot.m
		i. parameter_tunning_c1_c2
		j. parameter_tunning_population.m
		k. incorrect_density_godavari
		l. incorrect_density_sayula
		m. example_arbritary_basin
	
	a. fourier_coef_arbitary_model.m - It calculates the Fourier power spectrum plot for a arbritary model shown in Figure 4. 
	
	b. final_gravity_model1_fixed.m- It calculates the inversion of gravity anomaly for a synthetic sedimentary basin having fixed density contrast with and without noise case (Model1) and the Fourier Spectrum. Output of the file shown in figure 5. 

	c. final_gravity_model1_varying.m- It calculates the inversion of gravity anomaly for a synthetic sedimentary basin having varying density contrast with and without noise case (Model1) and the Fourier Spectrum. Output of the file shown in figure 6. 

	d. final_gravity_model2_fixed.m- It calculates the inversion of gravity anomaly for a synthetic sedimentary basin having fixed density contrast with and without noise case (Model2) and the Fourier Spectrum. Output of the file shown in figure 7. 

	e. final_gravity_model2_varying.m- It calculates the inversion of gravity anomaly for a synthetic sedimentary basin having varying density contrast with and without noise case (Model2) and the Fourier Spectrum. Output of the file shown in figure 8. 

	f. godavari_basin.m - It calculates the depth profile of Godavari Basin using inversion of gravity anomaly and the Fourier power spectrum. Output is shown in figure 10.

	g. sayula_basin.m - It calculates the depth profile of Sayula Basin using inversion of gravity anomaly. Output is shown in figure 11.

	h. error_energy_plot.m - This code plots error energy as shown in figure 9. in the manuscript.  
	
	i. parameter_tunning_c1_c2 - This code for Parameter tuning of PSO for cognitive (c1) and social (c2) components. Output of the code is shown in Table 2. 
	
	j. parameter_tunning_population.m - This code for Parameter tuning of particle population. Output of the code is shown in Table 3. 
	
	k. incorrect_density_godavari.m - Produce plots for inversion of depth profile for 
Godavari basin with incorrect density distribution.

	l. incorrect_density_sayula.m - Produce plots for inversion of depth profile for 
Sayula basin with incorrect density distribution.

	m. example_arbritary_basin - Example code for inversion of any real sedimentary basin using the dedicated function of Fourier_PSO.m. 
