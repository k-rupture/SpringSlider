// === friction.cpp ===				                

#define _USE_MATH_DEFINES
#include "friction.h"

// -------------------- Set Parameters --------------------

void Friction :: setDataForRSF(){

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
	
	// memory allocation
	z1 = new double[2] ;
	z2 = new double[2] ;
	z3 = new double[2] ;

	return ;
}

// -------------------- Set Initial Condition --------------------

void Friction :: setInitialCondition(){

	z1[0] = 0.000 ; // slip
	z2[0] = 0.606 ; // state
	z3[0] = 1 ;

	z1[1] = 0.0 ;
	z2[1] = 0.0 ;
	z3[1] = 0.0 ;
}

// -------------------- Run Spring Slider --------------------

void Friction :: SpringSliderRSF(){

	// Equation: t_s(time) - K*delta - etha*V = t_n(time)*f(V,si)

	// Runge-Kutta parameters
	double k1 , k2 , k3 ;
	double l1 , l2 , l3 ;
	double Z1nPlus2nd , Z1nPlus3rd ;
	double Z2nPlus2nd , Z2nPlus3rd ;

	// set output files
	ofstream fout0("SlipVelocityVsTime.txt") ; // Slip Velocity vs. time
	ofstream fout1("SlipVsTime.txt") ; // Slip vs. time
	ofstream fout2("StateVsTime.txt") ; // State vs. time
	ofstream fout7("ShearLoadingVsNormalLoading.txt") ; // Shear loadding vs. normal loading 
	ofstream fout10("VelocityVsTimeStar.txt") ; // Slip Velocity vs. time* 
	ofstream fout11("TotalShearVsNormalLoading.txt") ; // Total shear vs. normal loading
	ofstream fout17("SlipVsTimeStar.txt") ; // Slip vs. time*
	ofstream fout18("ShearVsTimeStar.txt") ; // Shear vs. time*

	int counter = 0 ;
	setInitialCondition() ;

	// print initial condition
	velocity = BracketedNewtonRaphsonHF4(z1[0] , z2[0] , current_time) ;
	fout0 << current_time << "  " << velocity << endl ;
	fout1 << current_time << "  " << z1[0] << endl ;
	fout2 << current_time << "  " << z2[0] << endl ;

	double Error_rk = 0.0 ;

	// loop over time
	while(current_time < final_time){
		
			k1 = step_time*dZ1dtHF4(z1[0], z2[0] , 0. , current_time) ;
			l1 = step_time*dZ2dtHF(z1[0], z2[0] , 0. , current_time) ;
			
			k2 = step_time*dZ1dtHF4(z1[0]+ 0.5*k1 , z2[0] + 0.5*l1 , 0 , current_time + 0.5*step_time) ;
			l2 = step_time*dZ2dtHF(z1[0]+ 0.5*k1 , z2[0] + 0.5*l1 , 0 , current_time + 0.5*step_time) ;
			
			k3 = step_time*dZ1dtHF4(z1[0]+ 2.0*k2 - k1 , z2[0]+ 2.0*l2 - l1 , 0 , current_time + step_time) ;
			l3 = step_time*dZ2dtHF(z1[0]+ 2.0*k2 - k1 , z2[0]+ 2.0*l2 - l1 , 0 , current_time + step_time) ;
			
			Z1nPlus2nd = z1[0] + k2 ;
			Z1nPlus3rd = z1[0] + (1.0/6.0)*(k1 + 4.0*k2 + k3) ;
			
			Z2nPlus2nd = z2[0] + l2 ;
			Z2nPlus3rd = z2[0] + (1.0/6.0)*(l1 + 4.0*l2 + l3) ;
			
			Error_rk = max(absolute(Z1nPlus2nd - Z1nPlus3rd),absolute(Z2nPlus2nd - Z2nPlus3rd)) ;
			
			if(Error_rk > tolerance){
				
				cout << "Time Step Rejected" << endl ;
				step_time = 0.8*pow(tolerance/Error_rk , 1.0/3.0)*step_time ;
				step_time = min(step_time , dt_max) ;

				cout << "Update Time" << endl ;
				cout << step_time << " New Time Step" << endl ;
				continue ;

			} 
			
			cout << step_time << "  step" << endl ;
			current_time += step_time ;
			
			cout << current_time << " Current Time" << endl ;

			velocity = BracketedNewtonRaphsonHF4(Z1nPlus3rd , Z2nPlus3rd , current_time) ;
			
			double t_star = current_time/(H/U) ;
			double sigma_n = computeNormalStressV4(current_time) ; // from hydraulic fracturing
			double sigma_s = computeTangentStressV4(current_time) ; // from hydraulic fracturing
			if(sigma_n > 0){sigma_s = 0 ;}
			
			fout0 << current_time << "  " << velocity << endl ;
			fout1 << current_time << "  " << Z1nPlus3rd << endl ;
			fout2 << current_time << "  " << Z2nPlus3rd << endl ;
			fout7 << sigma_n/pressure << "  " << sigma_s/pressure << endl ;
			fout10 << t_star - t_not_star << "  " << velocity << endl ;
			fout11 << sigma_n/pressure << "  " << sigma_s/pressure - K*Z1nPlus3rd/pressure - etha*velocity/pressure << endl ; 
			fout17 << t_star - t_not_star << "  " << Z1nPlus3rd << endl ;
			fout18 << t_star - t_not_star << "  " << sigma_s/pressure - K*Z1nPlus3rd/pressure - etha*velocity/pressure << endl ;
				
			z1[1] = Z1nPlus3rd ;
			z2[1] = Z2nPlus3rd ;
			
			if(Error_rk != 0.0){step_time = 0.8*pow(tolerance/Error_rk , 1.0/3.0)*step_time ;}
			if(Error_rk == 0.0){step_time = 2.01*step_time ;}
			
			step_time = min(step_time , dt_max) ;

			counter ++ ;
			Error_rk = 0.0 ;
			updateMemory() ;
		}
	
	return ;
}

// ------------ Bracketed Newton Method + H.F. ---------------

double Friction :: BracketedNewtonRaphsonHF4(double z1 , double z2, double time){

	double sigma_n = computeNormalStressV4(time) ; // from H.F.
	double sigma_s = computeTangentStressV4(time) ; // from H.F.
	sigma_n = -sigma_n ;
	if(sigma_n < 0){sigma_s = 0 ;} /// Correct

	double aa = 0.0 ;
	double bb = (sigma_s - K*(z1))/etha ;
	
	if(sigma_s - K*(z1) < 0){aa = bb ; bb = 0.0 ;}

	double A = FunctionHF4(aa,z1,z2,time) ;
	double B = FunctionHF4(bb,z1,z2,time) ;

	double x, X ;
	double TOL = pow(10.0 , -6.0) ;
	double absoluteTOL = pow(10.0 , -22.0) ;

	while(true){
		
			x = bb - B/derivativeFunctionHF4(bb,z1,z2,time);
	
			if (x < min(aa,bb) || x > max(aa,bb)){x = bb + (aa-bb)/2.0 ;}
			if (absolute(x - bb) < absolute(x)*TOL + absoluteTOL){return x ;}

			X = FunctionHF4(x,z1,z2,time) ;
			
			if ((B<0 && X>0) || (B>0 && X<0)){aa = bb ; A = B ;}
			
			bb = x ; 
			B = X ;
		}
}

// -------------------- dZ1/dt + HF --------------------

double Friction :: dZ1dtHF4(double z1 , double z2 , double z3 , double time){

	double V = 0.0 ;
	V = BracketedNewtonRaphsonHF4(z1 , z2, time);
	velocity = V ;
	return  V ;
} 

// -------------------- dZ2/dt + HF --------------------

double Friction :: dZ2dtHF(double z1 , double z2 , double z3 , double time){
	
	double V = velocity ;
	if(V == 0){return 0 ;}
	if(V < 0){V = -V ;}
	double f = a*ArcSinh(V*exp(z2/a)/(2.0*Vnot)) ;
	return -V*(f-fnot+(b-a)*log(V/Vnot))/dc ;
} 

// -------------------- Function --------------------

double Friction :: FunctionHF4(double V, double z1 , double z2 , double time){

	double sigma_n = computeNormalStressV4(time) ; // from H.F.
	double sigma_s = computeTangentStressV4(time) ; // from H.F.
	sigma_n = -sigma_n ;
	
	return sigma_n*a*ArcSinh((V/(2*Vnot))*exp(z2/a)) + K*(z1) + etha*V - sigma_s ;

	// Note: if we would be in the tensile regime, it does not matter.
	// because then the code never uses this function.
}

// -------------------- Derivative Function --------------------

double Friction :: derivativeFunctionHF4(double V, double z1 , double z2 , double time){

	double sigma_n = computeNormalStressV4(time) ; // from H.F.
	double sigma_s = computeTangentStressV4(time) ; // from H.F.
	sigma_n = -sigma_n ;

	double part1 = a*exp(z2/a)*sigma_n ;
	double part2 = 4+exp(2*z2/a)*V*V/(Vnot*Vnot) ;
	double part3 = Vnot ;

	return etha + part1/(sqrt(part2)*part3) ;

	// Note: if we would be in the tensile regime, it does not matter.
	// because then the code never uses this function.
}

// -------------------- Arc Sinh --------------------

double Friction :: ArcSinh(double x){

	if (x >= 0){return log(x + sqrt(x*x+1.0)) ;}
	if (x < 0){return -log(-x + sqrt(x*x+1.0)) ;}
}

// -------------------- Update Memory --------------------

void Friction :: updateMemory(){

	z1[0] = z1[1] ;
	z2[0] = z2[1] ;
	z3[0] = z3[1] ;
}

// ------------- Compute Sigma Xp-Xp (With Dimension) -----------------

double Friction :: computeSigmaXpXpV4(double t){

	double pi = M_PI ;
	double thetaPrime = atan(1.0/((X_not - U*t)/H)) ;
	if(thetaPrime < 0){thetaPrime += pi ;}

	double radius = sqrt(H*H + (X_not-U*t)*(X_not-U*t)) ;

	return -sigma1i + (Ki/sqrt(2.0*pi*radius))*cos(thetaPrime/2.0)*(1.0 - sin(thetaPrime/2.0)*sin(3.0*thetaPrime/2.0)) ;
}

// ------------- Compute Sigma Yp-Yp (With Dimension) -----------------

double Friction :: computeSigmaYpYpV4(double t){

	double pi = M_PI ;
	double thetaPrime = atan(1.0/((X_not - U*t)/H)) ;
	if(thetaPrime < 0){thetaPrime += pi ;}

	double radius = sqrt(H*H + (X_not-U*t)*(X_not-U*t)) ;

	return -pressure + (Ki/sqrt(2.0*pi*radius))*cos(thetaPrime/2.0)*(1.0 + sin(thetaPrime/2.0)*sin(3.0*thetaPrime/2.0)) ;
}

// ------------- Compute Sigma Xp-Yp (With Dimension) -----------------

double Friction :: computeSigmaXpYpV4(double t){

	double pi = M_PI ;
	double thetaPrime = atan(1.0/((X_not - U*t)/H)) ;
	if(thetaPrime < 0){thetaPrime += pi ;}

	double radius = sqrt(H*H + (X_not-U*t)*(X_not-U*t)) ;

	return 0.0 + (Ki/sqrt(2.0*pi*radius))*cos(thetaPrime/2.0)*sin(thetaPrime/2.0)*cos(3.0*thetaPrime/2.0) ;
}

// -------------------- Compute Normal Stress --------------------

double Friction :: computeNormalStressV4(double t){

	double SigmaXpXp = computeSigmaXpXpV4(t);
	double SigmaYpYp = computeSigmaYpYpV4(t);
	double SigmaXpYp = computeSigmaXpYpV4(t);

	return (SigmaXpXp + SigmaYpYp)/2.0 - (SigmaXpXp - SigmaYpYp)*cos(2.0*alfa)/2.0 - SigmaXpYp*sin(2.0*alfa) ;
}

// -------------------- Compute Tangent Stress --------------------

double Friction :: computeTangentStressV4(double t){

	double SigmaXpXp = computeSigmaXpXpV4(t);
	double SigmaYpYp = computeSigmaYpYpV4(t);
	double SigmaXpYp = computeSigmaXpYpV4(t);

    return (SigmaYpYp - SigmaXpXp)*sin(2.0*alfa)/2.0 + SigmaXpYp*cos(2.0*alfa) ;
}

// -------------------- Compute Absolute Value --------------------

double Friction :: absolute(double value){

	if(value >= 0. ){return value ;}
	return -1.0*value ;
}

// -------------------- Set Data For Coulomb Failure Criteria --------------------

void Friction :: setDataForCFC(){

	// Coulomb Failure Criteria 
	fnot = 0.6 ; // Friction coefficient
	Vnot = pow(10.0 , -6.0); // Initial slip velocity
	etha = 4.41*pow(10.,6.); // Radiation damping // MPa/(m/s)

	// H.F. information
	double pi = M_PI ; // Pi number
	alfa = 35.0*pi/180.0 ; // Fault orientation // degree
	pressure = 30.0*pow(10.0,6.0); // Pressue inside the crack (p) // Pa 
	sigma1i = 2.5*pressure; // Maximum stress principle (sigma_1) // Pa 
	Ki = 25.*pow(10. , 6.) ;  // Stress intensity factor
	U = 0.1 ; // Fracture tip velocity // m/s
	mu = 30.0*pow(10, 9.0); // Shear modules // Pa
	X_not = 13.5*pow(10.0 , 4.0); // Initial distance // m
	
	double L_fault = 0.4 ; // Fault Length // m
	double poisson_ratio = 0.25 ;
	double mu_star = mu/(1.0 - poisson_ratio);
	K = mu/((1.0-poisson_ratio)*L_fault) ; // Fault Stiffness
	H = 1.3 ; // Fault distance from Hydraulic Fracture // m
	L_star = t_not_star = X_not/H ;
	
	// Time integration
	final_time = 2.0*X_not/U ; // time // s
	step_time = 100 ; // Starting Point // s
	dt_max = 0.01*X_not/U; // s
	current_time = 0.00 ; // s
	tolerance = pow(10.0 , -6.0) ; // Tolerance ;
	z1 = new double[2] ;
	z2 = new double[2] ;
	z3 = new double[2] ;

	return ;
}

// -------------------- Spring Slider for Coulomb Failure Criterion --------------------

void Friction :: SpringSliderCFC(){

	// Equation: t_s(time) - K*delta - etha*V = t_n(time)*fo
	
	double sss = 0 ; // slip

	ofstream fout7("ShearLoadingVsNormalLoading.txt") ; // Shear loadding vs. Normal loading 
	ofstream fout10("VelocityVsTimeStar.txt") ; // V vs. time* 
	ofstream fout11("TotalShearVsNormalLoading.txt") ;
	ofstream fout17("SlipVsTimeStar.txt") ; // Slip vs. time*
	ofstream fout18("ShearVsTimeStar.txt") ; // Shear vs. time*

	double shear_lock ;

	double sigma_n = computeNormalStressV4(current_time) ; // from H.F.
	double sigma_s = computeTangentStressV4(current_time) ; // from H.F.
			
	// Print Initial Condition
	shear_lock = sigma_s - K*0.0  ; 
	if(shear_lock > fnot*absolute(sigma_n)){velocity = (shear_lock - fnot*absolute(sigma_n))/etha ;}
	if(shear_lock <= fnot*absolute(sigma_n)){velocity = 0.0 ;}

	while(current_time < final_time){
				
			cout << step_time << "  step " << endl ;
			current_time += step_time ;
			
			cout << current_time << " Current Time" << endl ;

			double sigma_n = computeNormalStressV4(current_time) ; // from H.F.
			double sigma_s = computeTangentStressV4(current_time) ; // from H.F.
			
			if(sigma_n > 0){sigma_s = 0 ;}
			shear_lock = sigma_s - K*0.0  ; 
			
			if(shear_lock > fnot*absolute(sigma_n)){velocity = (shear_lock - fnot*absolute(sigma_n))/etha ;}
			if(shear_lock <= fnot*absolute(sigma_n)){velocity = 0.0 ;}

			sss = sss + step_time*velocity ;

			double t_star = current_time/(H/U) ;
			
			fout7 << sigma_n/pressure << "  " << sigma_s/pressure << endl ;
			fout10 << t_star - t_not_star << "  " << velocity << endl ;
			fout11 << sigma_n/pressure << "  " << sigma_s/pressure - K*sss/pressure - etha*velocity/pressure << endl ; 
			fout17 << t_star - t_not_star << "  " << sss << endl ;
			fout18 << t_star - t_not_star << "  " << sigma_s/pressure - K*sss/pressure - etha*velocity/pressure << endl ;

		}

	return ;
}

// -------------------- Destructor --------------------

Friction :: ~Friction(){}