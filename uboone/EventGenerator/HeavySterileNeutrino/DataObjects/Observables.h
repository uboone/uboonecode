#ifndef STERILE_OBSERVABLES_H_
#define STERILE_OBSERVABLES_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

namespace hsngen{

	// Define struct OBSERVABLES
	class Observables {
	public:
		// Constructors
		Observables();

		// Attributes
		// Heavy sterile neutrino attributes
		double E_sterile; 
		double Th_sterile; 
		// Collective attributes
		int chan_identifier;
		double E_sum; 	
		double Th_sum; 
		double AngSep;
		double Minvar;
		double FS_AngSep; //The foreshortened angular separation.
		// IP1 attributes
		double E1; 
		double Th1;
	  std::vector<double> P1;
		double mass1;
		int pdg1;
		// IP2 attributes
		double E2; 
		double Th2;
		std::vector<double> P2;
		double mass2;
		int pdg2;
		// Decay vertex coordinates in time/space
		double xPos;
		double yPos;
		double zPos;
		double time;

		// Functions
		void PrintHepEvt(int i);
	};

} // END namespace hsngen

#endif