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
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/MCBase/MCTrack.h"

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
  int verbosity;
  std::string fMCTruthProducer;
  std::string fMCTrackProducer;

  bool  fMakeFVCut;
  float fdWall;

  bool fMakeContainmentCut;
  float fDist2Wall;

  bool fMakeEnergyCut;
  std::vector<float> fEnuTrueRange;

  bool fMakeProtonCut;
  float fProtonMinKE;

  bool fMakeLeptonCut;
  float fLeptonMinKE;

  bool fMakeNCCut;


  float dwall( const std::vector<float>& pos );

};


DLSignalSample::DLSignalSample(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  verbosity = p.get<int>("Verbosity",0);
  fMCTruthProducer = p.get< std::string >("MCTruthProducer","generator" );
  fMCTrackProducer = p.get< std::string >("MCTrackProducer","mcreco" );

  fMakeFVCut = p.get<bool>("MakeFVCut");
  fdWall = p.get<float>("dWallcm",0.0);

  fMakeContainmentCut = p.get<bool>("MakeContainmentCut");
  fDist2Wall = p.get<float>("Dist2Wall",10.0);

  fMakeEnergyCut = p.get<bool>("MakeEnergyCut");
  fEnuTrueRange = p.get< std::vector<float> >("EnuRangeMeV");
  if ( fEnuTrueRange.size()!=2 || fEnuTrueRange[0]>fEnuTrueRange[1] ) {
    std::stringstream msg;
    msg << __FILE__ << ":" << __LINE__ << "EnuRangeMeV invalid." << std::endl;
    throw std::runtime_error(msg.str());
  }

  fMakeProtonCut = p.get<bool>("MakeProtonCut");
  fProtonMinKE = p.get< float >( "ProtonMinKE_MeV", 30 );

  fMakeLeptonCut = p.get<bool>("MakeLeptonCut");
  fLeptonMinKE = p.get< float >( "LeptonMinKE_MeV", 30 );

  fMakeNCCut = p.get<bool>("MakeNCCut");
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
  bool ismuon = false; // if lepton produced is muon
  bool isnc   = false;
  std::vector<float> vertex( 3, 1e5 );

  for ( auto const& truthdata : (*truthHandle) ) {
    // check the origin. we want beam neutrinos.
    if ( truthdata.Origin()!=simb::Origin_t::kBeamNeutrino ) {
      continue;
    }

    // get the neutrino
    auto const& neutrino = truthdata.GetNeutrino();

    if ( fMakeNCCut ) {
      if ( neutrino.CCNC()==simb::kNC )
	isnc=true;
    }

    // check the particles
    float max_proton_mom = 0.;
    for (int ipart=0; ipart<truthdata.NParticles(); ipart++) {
      auto const& particle = truthdata.GetParticle(ipart);
      
      // neutrino or other handled differently
      if ( particle.PdgCode()==14 || particle.PdgCode()==-14
	   || particle.PdgCode()==12 || particle.PdgCode()==-12
	   || particle.PdgCode()==16 || particle.PdgCode()==-16 ) {
	vertex[0] = particle.Vx(0);
	vertex[1] = particle.Vy(0);
	vertex[2] = particle.Vz(0);
	float posdwall = dwall( vertex );
	if ( posdwall>fdWall )
	  in_fv = true;
	
	if ( particle.E(0)*1000.0>fEnuTrueRange[0] && particle.E(0)*1000.0<fEnuTrueRange[1] )
	  enucut = true;
	
	if ( verbosity>0 ) {
	  std::cout << "Found Neutrino:" << std::endl;
	  std::cout << "  vertex (" << vertex[0] << "," << vertex[1] << "," << vertex[2] << ")" << std::endl;
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
	ismuon = true;
	// std::vector<float> endpos(3,0);
	// endpos[0] = particle.EndX();
	// endpos[1] = particle.EndY();
	// endpos[2] = particle.EndZ();
	// float end_dwall = dwall( endpos );
	// if ( end_dwall < fDist2Wall )
	//   contained = false;
	if ( verbosity>0 ) {
	  std::cout << "Found muon:" << std::endl;
	  //   std::cout << "  end pos (" << endpos[0] << "," << endpos[1] << "," << endpos[2] << ")" << std::endl;
	  //   std::cout << "  dist2wall: " << end_dwall << std::endl;
	  //   std::cout << "  contained: " << contained << std::endl;
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
    }
    
  }//end of loop over truth classes
  
  
  if ( ismuon ) {
    // we need to search MCTrack information for end point to get containment
    art::Handle< std::vector<sim::MCTrack> > mctrackHandle;
    e.getByLabel(fMCTrackProducer,mctrackHandle);
    
    if ( !mctrackHandle.isValid() ) {
      std::stringstream msg;
      msg << __FILE__ << ":" << __LINE__ << " Unable to load MCTrack data with producer name " << fMCTrackProducer << std::endl;
      throw std::runtime_error( msg.str() );
    }

    // loop over the mctracks
    int closest_idx = -1;
    float closest_dist = -1;
    float endpt_dwall = 0;
    int n_within_1cm = 0;
    int idx = -1;
    for ( auto const& track : *mctrackHandle ) {
      // look for a muon whose start point is the same as the neutrino
      idx++;
      if ( track.PdgCode()==13 || track.PdgCode()==-13 ) {
	// a muon
	auto const&  start = track.Start();
	auto const&  end   = track.End();
	float dist = 0.;
	dist += (start.X()-vertex[0])*(start.X()-vertex[0]);
	dist += (start.Y()-vertex[1])*(start.Y()-vertex[1]);
	dist += (start.Z()-vertex[2])*(start.Z()-vertex[2]);
	dist = sqrt(dist);
	if ( dist<1.0 )
	  n_within_1cm++;
	if ( dist<5.0 && (closest_idx<0 || dist<closest_dist) ) {
	  closest_idx = idx;
	  closest_dist = dist;
	  std::vector<float> endpos(3,0);
	  endpos[0] = end.X();
	  endpos[1] = end.Y();
	  endpos[2] = end.Z();
	  endpt_dwall = dwall( endpos );
	  if ( verbosity>0 ) {
	    std::cout << "Found MCTRACK muon candidate:" << std::endl;
	    std::cout << "  dist2vertex: " << dist << std::endl;
	    std::cout << "  end pos (" << endpos[0] << "," << endpos[1] << "," << endpos[2] << ")" << std::endl;
	    std::cout << "  dist2wall: " << endpt_dwall << std::endl;
	  }
	}//if vertex match qualifies for check
      }//if vertex is muon
    }//end of track loop
    
    if ( n_within_1cm==0 || endpt_dwall<fDist2Wall )
      contained = false;
  }//end of ismuon block
  
  if ( verbosity>0 ) {
    std::cout << "Cut Summary: " << std::endl;
    std::cout << "  haslepton: " << haslepton << std::endl;
    std::cout << "  infv: " << in_fv << std::endl;      
    std::cout << "  inenurange: " << enucut << std::endl;      
    std::cout << "  contained: " << contained << std::endl;
    std::cout << "  minproton: " << minproton << std::endl;
    std::cout << "  isnc: " << isnc << std::endl;
  }
  

  // Event logic
  // ------------
  // Treat NC and CC differently
  
  if ( fMakeNCCut ) {
    if (!isnc)
      return false;

    // NC can make FV, Energy, proton cuts. 
    if ( fMakeFVCut && !in_fv )
      return false;

    if ( fMakeProtonCut && !minproton )
      return false;

    if ( fMakeEnergyCut && !enucut )
      return false;

    std::cout << "PASSES run subrun event: " << e.run() << " " << e.subRun() << " " << e.event() << std::endl;
    return true;
  }
  else {
    // No NC cut

    // When CC, can make FV, energy, lepton, proton cuts
    if ( fMakeFVCut && !in_fv )
      return false;

    if ( fMakeProtonCut && !minproton )
      return false;

    if ( fMakeEnergyCut && !enucut )
      return false;

    if ( fMakeLeptonCut && !haslepton )
      return false;

    if ( fMakeContainmentCut && !contained )
      return false;

    std::cout << "PASSES run subrun event: " << e.run() << " " << e.subRun() << " " << e.event() << std::endl;
    return true;
  }
  
  std::cout << "SHOULD NEVER GET HERE" << std::endl;
  throw std::runtime_error("DLSignalSample_module. should never get here.");

  return true;
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
