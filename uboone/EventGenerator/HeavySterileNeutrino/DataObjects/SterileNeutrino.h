#ifndef STERILE_NEUTRINO_H_
#define STERILE_NEUTRINO_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "FourMomentum.h"

namespace hsngen {

	// Class definining an initial sterile neutrino
	class SterileNeutrino {
	public:
		// Constructor
		SterileNeutrino(double M, double E, double CosTh, double Phi);

		// Attributes
		double mass;
		double energy;
		double cosTh;
		double phi;
		FourMomentum labFrameP;
	};

} // END namespace hsngen

#endif