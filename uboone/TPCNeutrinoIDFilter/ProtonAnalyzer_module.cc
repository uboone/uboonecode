////////////////////////////////////////////////////////////////////////
// Class:       ProtonAnalyzer
// Module Type: analyzer
// File:        ProtonAnalyzer_module.cc
//
// Generated at Tue Jan 19 14:51:53 2016 by Katherine Woodruff using artmod
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Framework/Core/FindManyP.h"

// ROOT Includes
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "Simulation/SimChannel.h"
#include "Simulation/sim.h"
#include "CalibrationDBI/Interface/IChannelStatusService.h"
#include "CalibrationDBI/Interface/IChannelStatusProvider.h"

// LArSoft data definitions
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "AnalysisBase/CosmicTag.h"

#include <memory>

namespace ProtonAnalyzer
{

class ProtonAnalyzer : public art::EDAnalyzer {

public:
  explicit ProtonAnalyzer(fhicl::ParameterSet const& pset);
  virtual ~ProtonAnalyzer();

  // Required functions.
  void analyze(const art::Event& e);

  // Selected optional functions.
  void beginJob();
  void endJob();
  void beginRun(const art::Run&);
  void reconfigure(fhicl::ParameterSet const&);

private:

  // Declare member data here.
  std::string fTrackModuleLabel;
  std::string fClusterAssocLabel;
  std::string fCosmicTaggerAssocLabel;
  std::string fHitAssocLabel;

  TTree       *fProtonTracks;
  int         fevent;
  float       ftracklen;
  float       faverageintcharge;
  float       ftotalintcharge;
  float       fclosestfriend;
  float       ftrackphi;
};


ProtonAnalyzer::ProtonAnalyzer(fhicl::ParameterSet const & pset)
    : EDAnalyzer(pset)
{
  this->reconfigure(pset);
}

 ProtonAnalyzer::~ProtonAnalyzer()
{}

void ProtonAnalyzer::analyze(const art::Event& e)
{
  // Recover the handles to the track collection we want to analyze.
  art::Handle< std::vector<recob::Track> > trackVecHandle;
  e.getByLabel(fTrackModuleLabel, trackVecHandle);
 
  // Find dead wires to avoid
  const lariov::IChannelStatusProvider& ChannelStatusProvider 
     = art::ServiceHandle<lariov::IChannelStatusService>()->GetProvider();

  std::vector<unsigned int> deadwires;
  art::ServiceHandle<geo::Geometry> geo;
  const size_t N_CHANNELS = geo->Nchannels();
  for(unsigned int chan = 0; chan < N_CHANNELS; chan++) {
    if(ChannelStatusProvider.IsBad(chan) || !ChannelStatusProvider.IsPresent(chan))
    {
      deadwires.push_back(chan);
      std::cout << "deadwire: " << chan << std::endl;
    }
    std::cout << "total dead (or missing) wires found " << deadwires.size() << std::endl;
  }

  // Require valid handle, otherwise nothing to do
  if (trackVecHandle.isValid())
  {
    // Recover associations relating cosmic tags and track
    art::FindManyP<anab::CosmicTag> cosmicAssns(trackVecHandle, e, fCosmicTaggerAssocLabel);
    // Recover associations relating clusters and track
    //art::FindManyP<recob::Cluster> clusterAssns(trackVecHandle, e, fClusterAssocLabel);
    // Recover associations relating hits and track
    art::FindManyP<recob::Hit> hitAssns(trackVecHandle, e, fHitAssocLabel);
    // Make sure valid handles
    if (cosmicAssns.isValid() && hitAssns.isValid())
    {
      // Loop over input tracks
      for(size_t trackIdx = 0;trackIdx < trackVecHandle->size();trackIdx++)
      {
        art::Ptr<recob::Track> track(trackVecHandle,trackIdx);
        // Get cosmic tag score
        std::vector<art::Ptr<anab::CosmicTag> > cosmicVec = cosmicAssns.at(track.key());
        // Get assocciated clusters
        std::vector<art::Ptr<recob::Hit> > hitVec = hitAssns.at(track.key());
        // Skip if tagged cosmic
        if (!cosmicVec.empty())
        {
          art::Ptr<anab::CosmicTag>& cosmicTag(cosmicVec.front());
          if (cosmicTag->CosmicScore() > 0.4) continue;
        }
        // Get track length
        ftracklen = track->Length();
        // Get track phi
        ftrackphi = track->Phi();
        // Check average charge
        if (!hitVec.empty())
        {
          float totalcharge = 0.;
          for(auto const& hit: hitVec)
          {
            totalcharge += hit->Integral();
          }
          ftotalintcharge = totalcharge;
          faverageintcharge = totalcharge/float(hitVec.size());
        }
        else
        {
          ftotalintcharge = 0.;
          faverageintcharge = 0.;
        }
        // now find closest friend
        art::Handle< std::vector<recob::Track> > trackCompVecHandle;
        e.getByLabel(fTrackModuleLabel, trackCompVecHandle);
        // Require valid handle, otherwise nothing to do
        if (trackCompVecHandle.isValid())
        {
          Float_t mindist = 999.;
          std::cout << "trackCompVecHandle->size(): " << trackCompVecHandle->size() << std::endl;
          for(size_t compareTrackIdx = 0; compareTrackIdx < trackCompVecHandle->size(); compareTrackIdx++)
          {
            art::Ptr<recob::Track> ctrack(trackCompVecHandle,compareTrackIdx);
            std::cout << "track->ID(): " << track->ID() << "     ctrack->ID(): " << ctrack->ID() << std::endl;
            if (track->ID() != ctrack->ID())
            {
              // check begin/end of all tracks for some distance
              Float_t tmpdist;
              // Vertex -- Vertex
              tmpdist = TMath::Sqrt(TMath::Power(track->Vertex().X() - ctrack->Vertex().X(),2) + TMath::Power(track->Vertex().Y() - ctrack->Vertex().Y(),2) + TMath::Power(track->Vertex().Z() - ctrack->Vertex().Z(),2));
              mindist = TMath::Min(mindist,tmpdist);
              // Vertex -- End
              tmpdist = TMath::Sqrt(TMath::Power(track->Vertex().X() - ctrack->End().X(),2) + TMath::Power(track->Vertex().Y() - ctrack->End().Y(),2) + TMath::Power(track->Vertex().Z() - ctrack->End().Z(),2));
              mindist = TMath::Min(mindist,tmpdist);
              // End -- End
              tmpdist = TMath::Sqrt(TMath::Power(track->End().X() - ctrack->End().X(),2) + TMath::Power(track->End().Y() - ctrack->End().Y(),2) + TMath::Power(track->End().Z() - ctrack->End().Z(),2));
              mindist = TMath::Min(mindist,tmpdist);
              // End -- Vertex
              tmpdist = TMath::Sqrt(TMath::Power(track->End().X() - ctrack->Vertex().X(),2) + TMath::Power(track->End().Y() - ctrack->Vertex().Y(),2) + TMath::Power(track->End().Z() - ctrack->Vertex().Z(),2));
              std::cout << "tmpdist: " << tmpdist << std::endl;
              mindist = TMath::Min(mindist,tmpdist);
            }
          }
          fclosestfriend = mindist;
        }
        fevent = e.event();
        fProtonTracks->Fill();
      }
    }
  }
}

void ProtonAnalyzer::beginJob()
{
  // Access ART's TFileService, which will handle creating and writing
  // histograms and n-tuples for us. 
  art::ServiceHandle<art::TFileService> tfs;

  // Make tree and branches
  fProtonTracks = tfs->make<TTree>("ProtonTracks","ProtonTracks");
  fProtonTracks->Branch("event", &fevent, "event/I");
  fProtonTracks->Branch("tracklen", &ftracklen, "tracklen/F");
  fProtonTracks->Branch("trackphi", & ftrackphi, "trackphi/F");
  fProtonTracks->Branch("averageintcharge", &faverageintcharge, "averageintcharge/F");
  fProtonTracks->Branch("totalintcharge", &ftotalintcharge, "totalintcharge/F");
  fProtonTracks->Branch("closestfriend", &fclosestfriend, "closestfriend/F");
}

void ProtonAnalyzer::beginRun(art::Run const & r)
{
}

void ProtonAnalyzer::endJob()
{
}

void ProtonAnalyzer::reconfigure(fhicl::ParameterSet const & pset)
{
  fTrackModuleLabel       = pset.get<std::string> ("TrackModuleLabel", "trackkalmanhitc3d");
  fClusterAssocLabel      = pset.get<std::string> ("ClusterAssocLabel", "trackkalmanhitc3d");
  fCosmicTaggerAssocLabel = pset.get<std::string> ("CosmicTaggerAssocLabel", "trackkalmanhittagc3d");
  fHitAssocLabel          = pset.get<std::string> ("HitAssocLabel", "trackkalmanhitc3d");
}

DEFINE_ART_MODULE(ProtonAnalyzer)

}
