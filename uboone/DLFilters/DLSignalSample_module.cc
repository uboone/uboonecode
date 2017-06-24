////////////////////////////////////////////////////////////////////////
// Class:       DLSignalSample
// Plugin Type: filter (art v2_06_03)
// File:        DLSignalSample_module.cc
//
// Generated at Thu Jun 22 18:02:02 2017 by Taritree Wongjirad using cetskelgen
// from cetlib version v2_03_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCTruth.h"

#include <memory>

class DLSignalSample;


class DLSignalSample : public art::EDFilter {
public:
  explicit DLSignalSample(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  DLSignalSample(DLSignalSample const &) = delete;
  DLSignalSample(DLSignalSample &&) = delete;
  DLSignalSample & operator = (DLSignalSample const &) = delete;
  DLSignalSample & operator = (DLSignalSample &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

private:

  // Declare member data here.
  float fdWall;
  float fDist2Wall;
  float fProtonMinKE;
  std::string fMCTruthProducer;
  std::vector<float> fEnuTrueRange;
  int verbosity;

  float dwall( const std::vector<float>& pos );

};


DLSignalSample::DLSignalSample(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  verbosity = p.get<int>("Verbosity",0);
  fdWall = p.get<float>("dWallcm",0.0);
  fDist2Wall = p.get<float>("Dist2Wall",10.0);
  fEnuTrueRange = p.get< std::vector<float> >("EnuRangeMeV");
  fMCTruthProducer = p.get< std::string >("MCTruthProducer","generator" );
  if ( fEnuTrueRange.size()!=2 || fEnuTrueRange[0]>fEnuTrueRange[1] ) {
    std::stringstream msg;
    msg << __FILE__ << ":" << __LINE__ << "EnuRangeMeV invalid." << std::endl;
    throw std::runtime_error(msg.str());
  }
  fProtonMinKE = p.get< float >( "ProtonMinKE_MeV", 30 );
}

bool DLSignalSample::filter(art::Event & e)
{
  // Implementation of required member function here.
  
  // Get the truth data out of the art event
  // Data products are typically stored as vectors of class instances
  art::Handle< std::vector<simb::MCTruth> > truthHandle;
  e.getByLabel(fMCTruthProducer,truthHandle);

  if ( !truthHandle.isValid() ) {
    std::stringstream msg;
    msg << __FILE__ << ":" << __LINE__ << " Unable to load MCTruth data" << std::endl;
    throw std::runtime_error( msg.str() );
  }


  // now we loop through the info in the MCTruth class. We look for the neutrino interaction and check the final state particles from the interaction.
  bool in_fv = false; // dwall cut
  bool contained = true; // dist2wall. only applies to muons
  bool minproton = false; // final state proton with largest KE has threshold KE
  bool enucut = false; // within Enu Range
  bool haslepton = false; // CC only for now

  for ( auto const& truthdata : (*truthHandle) ) {
    // check the origin. we want beam neutrinos.
    if ( truthdata.Origin()!=simb::Origin_t::kBeamNeutrino ) {
      continue;
    }

    // get the neutrino
    //auto const& neutrino = truthdata.GetNeutrino();

    // check the particles
    float max_proton_mom = 0.;
    for (int ipart=0; ipart<truthdata.NParticles(); ipart++) {
      auto const& particle = truthdata.GetParticle(ipart);

      // neutrino or other handled differently
      if ( particle.PdgCode()==14 || particle.PdgCode()==-14
	   || particle.PdgCode()==12 || particle.PdgCode()==-12
	   || particle.PdgCode()==16 || particle.PdgCode()==-16 ) {
	std::vector<float> pos(3,0);
	pos[0] = particle.Vx(0);
	pos[1] = particle.Vy(0);
	pos[2] = particle.Vz(0);
	float posdwall = dwall( pos );
	if ( posdwall>fdWall )
	  in_fv = true;

	if ( particle.E(0)*1000.0>fEnuTrueRange[0] && particle.E(0)*1000.0<fEnuTrueRange[1] )
	  enucut = true;

	if ( verbosity>0 ) {
	  std::cout << "Found Neutrino:" << std::endl;
	  std::cout << "  vertex (" << pos[0] << "," << pos[1] << "," << pos[2] << ")" << std::endl;
	  std::cout << "  dwall: " << posdwall << std::endl;
	  std::cout << "  Energy: " << particle.E(0)*1000.0 << " MeV" << std::endl;
	  std::cout << "  InFV: " << in_fv << std::endl;
	  std::cout << "  InEnuRange: " << enucut << std::endl;
	}
      }
      else if ( particle.PdgCode()==2212 ) {
	if ( particle.P(0)*1000.0>max_proton_mom )
	  max_proton_mom = particle.P(0)*1000.0;
	if ( verbosity>0 ) {
	  std::cout << "Found proton:" << std::endl;
	  std::cout << "  init momentum: " << particle.P(0)*1000.0 << std::endl;
	  std::cout << "  updated max proton mom: " << max_proton_mom << std::endl;
	}
      }
      else if ( particle.PdgCode()==13 || particle.PdgCode()==-13 ) {
	haslepton = true;
	std::vector<float> endpos(3,0);
	endpos[0] = particle.EndX();
	endpos[1] = particle.EndY();
	endpos[2] = particle.EndZ();
	float end_dwall = dwall( endpos );
	if ( end_dwall < fDist2Wall )
	  contained = false;
	if ( verbosity>0 ) {
	  std::cout << "Found muon:" << std::endl;
	  std::cout << "  end pos (" << endpos[0] << "," << endpos[1] << "," << endpos[2] << ")" << std::endl;
	  std::cout << "  dist2wall: " << end_dwall << std::endl;
	  std::cout << "  contained: " << contained << std::endl;
	}
      }
      else if ( particle.PdgCode()==11 || particle.PdgCode()==-11 ) {
	haslepton = true;
	if ( verbosity>0 ) {
	  std::cout << "Found electron." << std::endl;
	}
      }
    }//end of loop over MC particle data in truth class
    float max_proton_ke = sqrt(max_proton_mom*max_proton_mom + 938.20*938.20) - 938.20;
    if ( max_proton_ke>fProtonMinKE )
      minproton = true;
    if ( verbosity>0 ) {
      std::cout << "Max proton KE: " << max_proton_ke << std::endl;
      std::cout << "Min Proton KE cut passes: " << minproton << std::endl;
      std::cout << "Cut Summary: " << std::endl;
      std::cout << "  haslepton: " << haslepton << std::endl;
      std::cout << "  infv: " << in_fv << std::endl;      
      std::cout << "  inenurange: " << enucut << std::endl;      
      std::cout << "  contained: " << contained << std::endl;
      std::cout << "  minproton: " << minproton << std::endl;
    }

  }//end of loop over truth classes

  if ( verbosity>0 && (in_fv & contained & minproton & enucut & haslepton) ) {
    std::cout << "PASSES run subrun event: " << e.run() << " " << e.subRun() << " " << e.event() << std::endl;
  }
  
  return in_fv & contained & minproton & enucut & haslepton;
}

float DLSignalSample::dwall( const std::vector<float>& pos ) {
  
  float dx1 = pos[0]-0.0;
  float dx2 = 258.0-pos[0];
  float dy1 = 117.0-pos[1];
  float dy2 = pos[1] + 117.0; // - -117.0
  float dz1 = pos[2];
  float dz2 = 1036.0-pos[2];
  
  float fdwall = 1.0e9;
  
  if ( dy1<fdwall ) {
    fdwall = dy1;
  }
  if ( dy2<fdwall ) {
    fdwall = dy2;
  }
  if ( dz1<fdwall ) {
    fdwall = dz1;
  }
  if ( dz2<fdwall ) {
    fdwall = dz2;
  }
  if ( dx1<fdwall ) {
    fdwall = dx1;
  }
  if ( dx2<fdwall ) {
    fdwall = dx2;
  }

  return fdwall;

}

DEFINE_ART_MODULE(DLSignalSample)
