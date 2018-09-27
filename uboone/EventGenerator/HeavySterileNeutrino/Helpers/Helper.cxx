#include "Helper.h"

namespace hsngen {
	
	// Generate a single heavy sterile neutrino event
	void GenerateObservables(CLHEP::HepRandomEngine& engine, twoIP_channel * CHAN, const FluxFile &flux, const Settings &set, Observables &obs)
	{
		// Flat random provider
		CLHEP::RandFlat flat(engine);

		//The flux files don't provide phi or cos angles for the steriles, so we generate them here. Also set up positions and timings
		double	cosTh = 0.999999;
		double  phi = 2.0*M_PI*flat();
		double  energy  = flux.GetRandomEvent(engine);

		// Produce an initial heavy sterile neutrinos
		SterileNeutrino nus(set.sterileMass, energy, cosTh, phi);

		// Calling the appropriate functions from the channels.
		switch(CHAN->chan_identifier)
		{
			case CHAN_MUMU:
			CHAN->decayfunctionMassive(nus,MMU,MMU,0.0);
			break;
			case CHAN_MUE:
			CHAN->decayfunctionMassive(nus,MMU,ME,0.0);
			break;
			case CHAN_ELECPI:
			case CHAN_MUONPI:
			case CHAN_NUPI0:
			case CHAN_ELECPOSI:
			case CHAN_GAMMA:
			CHAN->decayfunction(nus);
			break;
		}

		// Fill the observables with the proper values from the decay
		CHAN->observables(&obs, engine);
		obs.E_sterile = nus.energy;
		obs.Th_sterile = nus.cosTh;
		// Calculate position and timing information
		double speedOfLight = 299792458.; // m/s
		double delay = set.distance/speedOfLight * (obs.E_sterile/(sqrt(pow(obs.E_sterile,2.) - pow(set.sterileMass,2.))) - 1.); // Delay with respect to active neutrinos

		double sterileDelay = delay*1e9; // Convert to ns
		double globalTimeOffset = set.globalTimeOffset;
		double randomTimeOffset = set.beamWindow*flat();

    // Generate decay vertex uniformly in detector box given by boundaries
    double widthX = set.boundariesX[1]-set.boundariesX[0];
    double widthY = set.boundariesY[1]-set.boundariesY[0];
    double widthZ = set.boundariesZ[1]-set.boundariesZ[0];
		obs.xPos = (widthX)*flat() + set.boundariesX[0];
		obs.yPos = (widthY)*flat() + set.boundariesY[0];
		obs.zPos = (widthZ)*flat() + set.boundariesZ[0];
		obs.time = globalTimeOffset + randomTimeOffset + sterileDelay;
		return;
	} // END Function GenerateObservables


	// Fill model with theoretical parameters depending on the chosen decay channel.
	void FillModel(CLHEP::HepRandomEngine& engine, twoIP_channel *&CHAN, std::vector<double> &model_params, const Settings &set)
	{
		//  Horrible way to deal with model parameters, but was inherited from previous code.
		// model_params is a vector (of doubles) containing always 5 elements.
		// For 2-body decays they are: (mass of particle 1, mass of particle 2, channel id, pdg1, pdg2)
		// For 3-body decays they are: (boson mediator mass, channel id, channel id... again [?], pdg1, pdg2)
		// Since it is a vector of double, channel ids and pdg codes are converted to double for convenience,
		// then reconverted to int. Decay channels to containing a neutrino in final states have a symbolic
		// muon neutrino pdg code in it (could be any flavour).
		

		switch(set.decayChannel)
		{
			case CHAN_ELECPOSI:
			// model_params.push_back(91.19); //mediator mass
			// model_params.push_back((double) CHAN_ELECPOSI);
			// model_params.push_back((double) CHAN_ELECPOSI);
			// model_params.push_back((double) PDG_E);
			// model_params.push_back(-1 * (double) PDG_E);
			// CHAN = new threebody(engine,model_params);
      printf("Error! The e-e channel has not been cross-checked since the InFlight/LArSoft porting and it would be unwise to generate events without making sure first that it's working the way it should. You can try commenting out this exception in the code in Helpers/Helper.cxx at your own risk. Correct execution is not guaranteed.");
      std::exit(1);
			break;
			case CHAN_ELECPI:
			model_params.push_back(ME); // the electron with an incorrect mass (to try and avoid spacelike vectors...).
			model_params.push_back(MPION); // the pion. 0.135
			model_params.push_back((double) CHAN_ELECPI);
			model_params.push_back((double) PDG_E);
			model_params.push_back((double) PDG_PI);
			CHAN = new twobody(engine,model_params);
			break;
			case CHAN_MUONPI:
			model_params.push_back(MMU); // the muon.
			model_params.push_back(MPION); // the pion.
			model_params.push_back((double) CHAN_MUONPI);
			model_params.push_back((double) PDG_MU);
			model_params.push_back((double) PDG_PI);
			CHAN = new twobody(engine,model_params);
			break;
			case CHAN_NUPI0:
			// model_params.push_back(0.00); // the neutrino
			// model_params.push_back(MPI0); // the pion0
			// model_params.push_back((double) CHAN_NUPI0);
			// model_params.push_back((double) PDG_NUMU);
			// model_params.push_back((double) PDG_PI0);
			// CHAN = new twobody(engine,model_params);
      printf("Error! The pi-nu channel has not been cross-checked since the InFlight/LArSoft porting and it would be unwise to generate events without making sure first that it's working the way it should. You can try commenting out this exception in the code in Helpers/Helper.cxx at your own risk. Correct execution is not guaranteed.");
      std::exit(1);
			break;
			case CHAN_GAMMA:
			// model_params.push_back(91.19); //mediator mass
			// model_params.push_back((double) CHAN_GAMMA); // the pion.
			// model_params.push_back((double) CHAN_GAMMA); // the pion.
			// model_params.push_back((double) PDG_GAMMA);
			// model_params.push_back((double) PDG_NUMU);
			// CHAN = new threebody(engine,model_params);
      printf("Error! The gamma-nu channel has not been cross-checked since the InFlight/LArSoft porting and it would be unwise to generate events without making sure first that it's working the way it should. You can try commenting out this exception in the code in Helpers/Helper.cxx at your own risk. Correct execution is not guaranteed.");
      std::exit(1);
			break;
			case CHAN_MUMU:
			// model_params.push_back(91.19); //mediator mass
			// model_params.push_back((double) CHAN_MUMU); // the pion.
			// model_params.push_back((double) CHAN_MUMU); // the pion.
			// model_params.push_back((double) PDG_MU);
			// model_params.push_back(-1 * (double) PDG_MU);
			// CHAN = new threebody(engine,model_params);
      printf("Error! The mu-mu-nu channel has not been cross-checked since the InFlight/LArSoft porting and it would be unwise to generate events without making sure first that it's working the way it should. You can try commenting out this exception in the code in Helpers/Helper.cxx at your own risk. Correct execution is not guaranteed.");
      std::exit(1);
			break;
			case CHAN_MUE:
			// model_params.push_back(91.19); //mediator mass
			// model_params.push_back((double) CHAN_MUE); // the pion.
			// model_params.push_back((double) CHAN_MUE); // the pion.
			// model_params.push_back((double) PDG_MU);
			// model_params.push_back(-1 * (double) PDG_E);
			// CHAN = new threebody(engine,model_params);
      printf("Error! The e-e-nu channel has not been cross-checked since the InFlight/LArSoft porting and it would be unwise to generate events without making sure first that it's working the way it should. You can try commenting out this exception in the code in Helpers/Helper.cxx at your own risk. Correct execution is not guaranteed.");
      std::exit(1);
			break;

			default:
			std::cout<<"ERROR: Bad channel specifier."<<std::endl;
			return;
		}

		return;
	} // END function FillModel
	
} // END namespace hsngen