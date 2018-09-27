#ifndef STERILE_HELPER_H_
#define STERILE_HELPER_H_

#include <iostream>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>


#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"

#include "../DataObjects/FourMomentum.h"
#include "../DataObjects/SterileNeutrino.h"
#include "../DataObjects/Flux.h"
#include "../DataObjects/Observables.h"
#include "../DataObjects/Channel.h"
#include "Settings.h"

#define MPION  0.13957
#define MPI0   0.13497
#define MKAON  0.49367
#define MMU    0.10566
#define	ME     0.00051

#define PDG_PI  	211
#define PDG_PI0   111
#define PDG_K		  321
#define PDG_MU    13
#define	PDG_E     11
#define	PDG_GAMMA 22
#define	PDG_NUE   12
#define	PDG_NUMU  14

namespace hsngen {

	void FillModel(CLHEP::HepRandomEngine& engine, twoIP_channel *&CHAN, std::vector<double> &model_params, const Settings &set);
	void GenerateObservables(CLHEP::HepRandomEngine& g, twoIP_channel * CHAN, const FluxFile &flux, const Settings &set, Observables &obs);

} // END namespace hsngen

#endif