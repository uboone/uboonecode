#ifndef STERILE_FLUX_H_
#define STERILE_FLUX_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <fstream>

#include "CLHEP/Random/RandomEngine.h"
#include "FourMomentum.h"

namespace hsngen {

	// Struct defining a sterile neutrino flux
	class FluxFile {
	public:
		// Constructor and destructor
		FluxFile();
		virtual ~FluxFile();
		FluxFile(std::string fileName, double sterileMass);

		// Attributes
		std::string fileName;
		std::vector<double> fluxValues;
		std::vector<double> energyValues;
		double maxFlux;
		double sterileMass;

		//Functions
		double GetFlux(double energy) const;
		double GetRandomEvent(CLHEP::HepRandomEngine& r) const;
	};

	double Interpolate(std::vector<double > xValues, std::vector<double> yValues, double eInterp);

} // END namespace hsngen

#endif