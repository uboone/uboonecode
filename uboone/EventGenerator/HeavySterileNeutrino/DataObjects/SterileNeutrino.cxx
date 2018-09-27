#include "SterileNeutrino.h"

namespace hsngen {

	// SterileNeutrino constructor
	SterileNeutrino::SterileNeutrino(double M, double E, double CosTh, double Phi)
	{
		mass = M;
		energy = E;
		cosTh = CosTh; 
		phi = Phi;
		double totP = sqrt(energy*energy-M*M);	

		double temp[] = {totP*sqrt(1.0-cosTh*cosTh)*cos(phi), totP*sqrt(1.0-cosTh*cosTh)*sin(phi), totP*cosTh};

		std::vector<double> momentum(temp, temp + sizeof(temp)/sizeof(double));

		// Assign a fourmomentum to the labframe
		labFrameP.Populate(energy, momentum);
	} // END constructor SterileNeutrino

} // END namespace hsngen