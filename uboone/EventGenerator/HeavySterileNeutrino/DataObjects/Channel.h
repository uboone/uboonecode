#ifndef CHANNEL_H_
#define CHANNEL_H_

#include <iostream>
#include <cmath>
#include <vector>

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"
#include "FourMomentum.h"
#include "SterileNeutrino.h"
#include "Observables.h"

#define CHAN_ELECPOSI 0
#define CHAN_ELECPI 1
#define CHAN_MUONPI 2
#define CHAN_NUPI0 3
#define CHAN_GAMMA 4
#define CHAN_MUMU 5
#define CHAN_MUE 6

namespace hsngen {

	// Define mother class for all decay channels (to 2 ionizing particles)
	class twoIP_channel {
	public:
		// Constructor
		twoIP_channel(CLHEP::HepRandomEngine& g, std::vector<double> input);
		virtual ~twoIP_channel();

		// Attributes
		FourMomentum IP1;	//first outgoing particle 4 momentum.
		FourMomentum IP2;	//second outgoing particle 4 momentum.
		int chan_identifier;
		CLHEP::HepRandomEngine * r;
		std::vector<double> model_params;

		// Functions
		int observables(Observables * output, CLHEP::HepRandomEngine& g);
		virtual int decayfunction(SterileNeutrino nuS);	
		virtual int decayfunctionMassive(SterileNeutrino nuS, double m0, double m1, double m2);
	}; // END class twoIP_channel



	/* ###############################
	Below here we have a derived class for each channel
	############################### */



	class threebody : public twoIP_channel {
	public:
		threebody(CLHEP::HepRandomEngine& g, std::vector<double> input);
		int decayfunction(SterileNeutrino nuS);
		int decayfunctionMassive(SterileNeutrino nuS,double m0, double m1, double m2);

		struct PDF_CHOICE { 
			double Enu; 
			double cosThnu; 
			double Phinu; };
	private:
		int computeLabFrameVariablesMassive(SterileNeutrino nuS, double p0[4], double p1[4]);
		int computeLabFrameVariables(double mS, double Es, double costhS, double phiS, double restFrameParams[3]);
		double pdf_function(double x, double y, double mS, double mZprime, void * pointer);
		int rotor(double theta, double phi, double vector[3]);
		std::vector<double > generic_boost(double Ep, double px, double py, double pz, double Er, double rx, double ry, double rz);

		int drawRestFrameDist(CLHEP::HepRandomEngine& r, double mS, double mZprime, double output[3]);
		int drawRestFrameDistMassive(CLHEP::HepRandomEngine& r, double mS, double m0, double m1, double m2, double out0[4], double out1[4]);
	}; 


	/* ########################################################################
	This is the nu_s \to \nu Zprime \to \nu e+ e- channel (on-shell Zprime).
	######################################################################## */

	class Zprimeresonance : public twoIP_channel {
	public: 
		Zprimeresonance(CLHEP::HepRandomEngine& g, std::vector<double> input);
		int decayfunction(SterileNeutrino nuS);
	private:
		double fourvec_costheta(double FOURVEC[4]);
		double fourvec_cosphi(double FOURVEC[4]);
		double rot_boost(double costh, double phi, double gam, double FOURVEC[4]);	
	}; 

	/* ########################################################################
	Tradition 2body decay
	######################################################################## */

	class twobody : public twoIP_channel {
	public: 
		twobody(CLHEP::HepRandomEngine& g, std::vector<double> input);
		int decayfunction(SterileNeutrino nuS);
	}; 

} // END namespace hsngen

#endif
