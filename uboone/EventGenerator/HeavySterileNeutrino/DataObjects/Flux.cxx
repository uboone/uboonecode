#include "Flux.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandFlat.h"

namespace hsngen {

	// FluxFile constructors
	FluxFile::FluxFile()
	{}

	FluxFile::~FluxFile()
	{}

	FluxFile::FluxFile(std::string fileName, double inputSterileMass)
	{
		maxFlux = 0.0;
		sterileMass = inputSterileMass;

		if(fileName != "none"){
			// int k = 0;
			std::string strE;
			std::string strW;
			std::ifstream myfile (fileName);
			if (!myfile.is_open())
			{
				std::cout<<"#ERROR: flux::flux @flux.c, passed flux file does not exist"<<std::endl;
				exit(EXIT_FAILURE);
			}
			while(!myfile.eof()){
				myfile >> strE;
				myfile >> strW;
				double ien = atof(strE.c_str());
				double iflux = atof(strW.c_str());
				energyValues.push_back(ien);
				fluxValues.push_back(iflux);
				if(ien < sterileMass && iflux != 0.0  ){
					std::cout<<"#ERROR: flux file containts events with energy less than the sterile rest mass! mass "<<sterileMass<<" E "<<ien<<" flux "<<iflux<<std::endl;
					exit(EXIT_FAILURE);
				}
				if(iflux > maxFlux) maxFlux=iflux;
			}	

			energyValues.pop_back();
			fluxValues.pop_back();

			myfile.close();

			if(energyValues.size()!=fluxValues.size())
			{
				std::cout<<"WHAT? energyValues != fluxValues in size. error @interpolate() in sterileflux.cxx. Your flux file has uneven columns"<<std::endl;
				exit(EXIT_FAILURE);
			}
		}
	} // END constructor FluxFile


	// Get a flux value given a specific energy (by interpolating)
	double FluxFile::GetFlux(double energy) const
	{
		double ans = Interpolate(energyValues, fluxValues, energy);

		// If energy is lower than mass, there's obviously no flux. Linear interpolation gets confused when interpolating in that region, so enforce null flux
		if(energy < sterileMass) ans = 0.0;

		return ans;
	} // END function GetFlux


	// Get a random event from the flux distribution
	double FluxFile::GetRandomEvent(CLHEP::HepRandomEngine& r) const
	{
		// Initialize random flat generator
	  CLHEP::RandFlat flat(r);

	  // Parameters to start
		double yDraw = 0.0;
		double normFlux = 0.0;
		double energy = 0.0;
		// double current = 0;
		double eMax = energyValues.back();
		double eMin = energyValues.front();

		// Pick random points in flux spectrum space
		energy = (eMax-eMin)*flat()+eMin;
		normFlux = GetFlux(energy)/maxFlux;
		yDraw = flat();

		// When finally random draw is under flux curve, return energy value
		int num = 0;
		while(yDraw > normFlux)
		{
			energy = (eMax-eMin)*flat()+eMin;
			normFlux = GetFlux(energy)/maxFlux;
			yDraw = flat();
			num++;
		}

		return energy;
	} // END function GetRandomEvent


	// Interpolate between a set of x and y values
	double Interpolate(std::vector<double > xValues, std::vector<double> yValues, double eInterp)
	{
		// Deal with broken flux file
		if(xValues.size()!=yValues.size())
		{
			std::cout << "#ERROR: internal interpolate failure, vector of different size!" << std::endl;
			exit(EXIT_FAILURE);
		}

		// Rectangle in interpolation region
		double e0 = 0.0;
		double e1 = 0.0;
		double w0 = 0.0;
		double w1 = 0.0;
		
		// Deal with interpolation outside of range or at edges
		if(eInterp > xValues.back() || eInterp < xValues.front()) return 0.0;
		else if(xValues.front()==eInterp) return yValues.front();
		else if(xValues.back()==eInterp) return yValues.back();

		// Find rectangle around interpolation value
		for(unsigned int i=1; i<=xValues.size() ;i++)
		{
			if(xValues[i]==eInterp) return yValues[i];

			if(xValues[i] > eInterp)
			{
				e0 = xValues[i-1];
				w0 = yValues[i-1];

				e1 = xValues[i];
				w1 = yValues[i];
				break;
			}
		}

		// Calculate interpolated value
		double ans = w0 +(w1-w0)*(eInterp - e0)/(e1-e0);
			
		return ans; 
	} // END function Interpolate

} // END namespace hsngen