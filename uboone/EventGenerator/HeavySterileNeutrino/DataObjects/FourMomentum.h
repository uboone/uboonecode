#ifndef STERILE_FOURMOMENTUM_H_
#define STERILE_FOURMOMENTUM_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>

namespace hsngen {

	// Class definining an initial sterile neutrino
	// This assumes on-shell 4-momenta (timelike or null)
	class FourMomentum {
	public: 
		// Constructors
		FourMomentum(double energy, std::vector<double> momentum);
		FourMomentum();

		// Attributes
		double mass;
		double E;
		double modp;
		std::vector<double> p;

		// Functions
		int Populate(double energy, std::vector<double> momentum);
		void Print(std::string name);
		std::vector<double> Direction();
		double Gamma();
		int RotBoostFromParent(FourMomentum * parentP);
	};

} // END namespace hsngen

#endif