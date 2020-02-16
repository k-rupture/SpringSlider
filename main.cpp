// === main.cpp ===				                

#include "friction.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <sstream>

using namespace std ;

int main(){
	
	// Rate and State Friction Law (RSF)
	   Friction C ;
	   C.setDataForRSF() ;
	   C.SpringSliderRSF() ;

	// Coulomb Failure Criterion (CFC)
	// Friction C ;
	// C.setDataForCFC() ;
	// C.SpringSliderCFC() ;
	 
	return 0 ;
}
