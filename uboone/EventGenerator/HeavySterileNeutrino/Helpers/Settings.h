#ifndef STERILE_SETTINGS_H_
#define STERILE_SETTINGS_H_

namespace hsngen {

	typedef struct Settings{
		bool printHepEvt;
		double sterileMass;
    int sterileType;
		int decayChannel;
		std::string fluxFile;
		double distance;
		double globalTimeOffset;
		double beamWindow;
		std::vector<double> boundariesX;
		std::vector<double> boundariesY;
		std::vector<double> boundariesZ;
    std::vector<double> generatedTimeWindow;
	} Settings;

} // END namespace hsngen

#endif