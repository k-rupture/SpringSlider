# Spring-Slider Modeling of Microseismic Events and Fault Slip during Hydraulic Fracturing

## Requirements
You might need g++ compiler. A sample of make file is uploaded as well.

## Examples
Here we explain the code through two examples

### Example 1
Rate and State Friction Law (RSF)

#### Step 1
In friction.cpp go to the following function

      setDataForRSF() ;
      
#### Step 3
In this function, set the parameters as you wish. These parameters have been already set for reproducing the results shown in our manuscript.

  	// Rate-and-state friction law
	fnot = 0.6 ; // Friction coefficient
	a = 0.015 ; // direct effect parameter
	b = 0.020 ; // state evolution parameter
	Vnot = pow(10.0 , -6.0); // m/s
	etha = 4.41*pow(10.,6.); // radation damping // MPa/(m/s)
  	dc = 14.0*10.*pow(10.,-6.); // state evolution distance // m
	
	// Hydrualic fracturing
	double pi = M_PI ; // Pi number
	alfa = 35.0*pi/180.0 ; // Fault orientation (\alpha) // degree
	pressure = 30.0*pow(10.0,6.0); // Fluid pressure inside the crack (p)  // Pa 
	sigma1i = 2.5*pressure; // Maximum principle (\sigma_1) // Pa
	Ki = 25.*pow(10. , 6.) ; // Stress intensity factor (K_I) // MPa*Sqrt(m) 
	U = 0.1 ; // Speed of hydrulic fracture tip // m/s
	mu = 30.0*pow(10, 9.0); // Shear modules (\mu) // Pa
	X_not = 13.5*pow(10.0 , 4.0) ; // Initial distance of hydraulic fracture from fault // m 

	double poisson_ratio = 0.25 ; // Poisson ratio
	double L_fault = 1.0 ; // Fault Lenght // m 
	K = mu/((1.0-poisson_ratio)*L_fault) ; // Fault stifness ;
	H = 1.3 ; // distance of fault center from Hydraulic fracturing path // m
	L_star = t_not_star = X_not/H ;

	// time integration
	final_time = 2.0*X_not/U ; // s
	step_time = pow(10.0, -5.0) ; // s (this is the initial time step; since we are using adaptive time stepping.)
	dt_max = 0.01*X_not/U; // s
	current_time = 0.00 ; //s
	tolerance = pow(10.0 , -6.0) ;

#### Step 4
In main.cpp call the following functions:

	   Friction C ;
	   C.setDataForRSF() ;
	   C.SpringSliderRSF() ;

### Example 2
Coulomb Failure Criterion (CFC)
