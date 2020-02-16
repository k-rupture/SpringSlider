// === Friction Class ==============================				                
// === Description: The Rate and State Friction Law
// =================================================

#ifndef FRICTION_H
#define FRICTION_H

#include <time.h>
#include <math.h>
#include <complex> 
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using namespace std ;

class Friction{
  
 public:

	double *z1, *z2, *z3 ;
	double final_time , step_time , current_time ;
	double tolerance ;
	double sigma, a, Vnot, V_L, etha, K, fnot, b, dc ; 
	double velocity ;
	double ts , tn ;
	double k1_local , k2_local , k3_local ;
	double l1_local , l2_local , l3_local ;
	double Z1nPlus2nd_local , Z1nPlus3rd_local ; 
	double Z2nPlus2nd_local , Z2nPlus3rd_local ; 	
	double error ;
	double Delta , State ;
	double x , y ; // position of the node in 2D domain
	double KI, sigma1 ; // KI and sigma1 have no dimensions					  
	double Ki, sigma1i ; // Ki and sigma1i has dimension.
	double alfa ; // alfa is in "rad"
	double pressure, U, H, mu ; // p inside the crack
	double L_star, X_not, t_not_star ; // for H.F. method 2
	double dt_max , dt_min ;
	
 public:

	void updateMemory() ;
	void setInitialCondition() ;
	
	// Math Functions
	double ArcSinh(double x) ;
	double absolute(double value) ;

	// Functions for Runge-Kutta
	double dZ1dtHF4(double z1 , double z2 , double z3 , double time);
	double dZ2dtHF(double z1 , double z2 , double z3 , double time);
	
	// Functions for solving the non-linear equation
	double FunctionHF4(double V, double z1 , double z2 , double time);
	double derivativeFunctionHF4(double V, double z1 , double z2 , double time);
	double BracketedNewtonRaphsonHF4(double z1 , double z2, double time);

	// Normal and Shear Stress Along the Fault 
	double computeNormalStressV4(double t) ; 
	double computeTangentStressV4(double t) ; 
	double computeSigmaXpYpV4(double t) ; 
	double computeSigmaXpXpV4(double t) ; 
	double computeSigmaYpYpV4(double t) ;

	// Rate and State Friction Law
	void setDataForRSF() ;
	void SpringSliderRSF() ;

	// Coulomb Failure Criterion
	void setDataForCFC() ;
	void SpringSliderCFC() ;

	~Friction();
};

#endif