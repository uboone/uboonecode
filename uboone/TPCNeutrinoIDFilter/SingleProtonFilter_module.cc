////////////////////////////////////////////////////////////////////////
// Class:       SingleProtonFilter
// Module Type: filter
// File:        SingleProtonFilter_module.cc
//
// Generated at Mon Jan 18 12:15:00 2016 by Katherine Woodruff using artmod
// from cetpkgsupport v1_10_01.
////////////////////////////////////////////////////////////////////////

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
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

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "SimpleTypesAndConstants/geo_types.h"

// LArSoft data definitions
#include "RecoBase/Track.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "AnalysisBase/CosmicTag.h"

#include <memory>

namespace SingleProtonFilter
{

class SingleProtonFilter : public art::EDFilter {

public:

  explicit SingleProtonFilter(fhicl::ParameterSet const & pset);
  virtual ~SingleProtonFilter();

  // Required functions.
  bool filter(art::Event&);

  // Selected optional functions.
  void beginJob();
  void endJob();
  void beginRun(const art::Run&);
  void reconfigure(fhicl::ParameterSet const&);

private:

  bool testTrack(art::Ptr<recob::Track>&, std::vector<art::Ptr<anab::CosmicTag>>&, std::vector<art::Ptr<recob::Hit>>&, art::Event&);

  // Declare member data here.
  std::string fTrackModuleLabel;
  std::string fClusterAssocLabel;
  std::string fCosmicTaggerAssocLabel;
  std::string fHitAssocLabel;

  float       fMaxCosmicScore;
  float       fMaxLength;
  float       fMinCharge;
  float       fMinStartDist;
};

SingleProtonFilter::SingleProtonFilter(fhicl::ParameterSet const & pset)
{
  this->reconfigure(pset);
}

SingleProtonFilter::~SingleProtonFilter()
{
}

bool SingleProtonFilter::filter(art::Event & e)
{
  bool pass = false;

  // Recover the handles to the track collection we want to analyze.
  art::Handle< std::vector<recob::Track> > trackVecHandle;
  e.getByLabel(fTrackModuleLabel, trackVecHandle);
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
      size_t trackIdx = 0;
      while(trackIdx < trackVecHandle->size() && pass == false)
      {
        art::Ptr<recob::Track> track(trackVecHandle,trackIdx);
        // Get cosmic tag score
        std::vector<art::Ptr<anab::CosmicTag> > cosmicVec = cosmicAssns.at(track.key());
        // Get assocciated clusters
        std::vector<art::Ptr<recob::Hit> > hitVec = hitAssns.at(track.key());
        pass = testTrack(track, cosmicVec, hitVec, e);
        if(pass)
        {
          std::cout << "track vertex z: " << track->Vertex().Z() << std::endl;
          std::cout << "track vertex y: " << track->Vertex().Y() << std::endl;
        }
        trackIdx++;
      }
    }
  }
  std::cout << "pass: " << pass << std::endl;
  return pass;
}

bool SingleProtonFilter::testTrack(art::Ptr<recob::Track> &track, std::vector<art::Ptr<anab::CosmicTag>> &cosmicVec, std::vector<art::Ptr<recob::Hit>> &hitVec, art::Event & e)
{
  // Check cosmic tag
  if (!cosmicVec.empty())
  {
    art::Ptr<anab::CosmicTag>& cosmicTag(cosmicVec.front());
    if (cosmicTag->CosmicScore() > fMaxCosmicScore) return false;
  }
  // Skip long tracks
  if (track->Length() > fMaxLength) return false;
  // Check average charge
  if (!hitVec.empty())
  {
    float totalcharge = 0.;
    for(auto const& hit: hitVec)
    {
      totalcharge += hit->Integral();
    }
    if (totalcharge/float(hitVec.size()) < fMinCharge) return false;
  }
  else { return false; } // Skip if no clusters
  // Last, check for shared vertex
  art::Handle< std::vector<recob::Track> > trackCompVecHandle;
  e.getByLabel(fTrackModuleLabel, trackCompVecHandle);
  // Require valid handle, otherwise nothing to do
  if (trackCompVecHandle.isValid())
  {
    Float_t mindist = fMinStartDist;
    for(size_t compareTrackIdx = 0; compareTrackIdx < trackCompVecHandle->size(); compareTrackIdx++)
    {
      art::Ptr<recob::Track> ctrack(trackCompVecHandle,compareTrackIdx);
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
        mindist = TMath::Min(mindist,tmpdist);
      }
      if (mindist < fMinStartDist) continue;
    }
    if (mindist < fMinStartDist) return false;
  }

  return true;
}

void SingleProtonFilter::beginJob()
{
}

void SingleProtonFilter::beginRun(const art::Run & r)
{
}

void SingleProtonFilter::endJob()
{
}

void SingleProtonFilter::reconfigure(fhicl::ParameterSet const & pset)
{
  fTrackModuleLabel       = pset.get<std::string> ("TrackModuleLabel", "trackkalmanhitc3d");
  fClusterAssocLabel      = pset.get<std::string> ("ClusterAssocLabel", "trackkalmanhitc3d");
  fCosmicTaggerAssocLabel = pset.get<std::string> ("CosmicTaggerAssocLabel", "trackkalmanhittagc3d");
  fHitAssocLabel          = pset.get<std::string> ("HitAssocLabel", "trackkalmanhitc3d");
  fMaxCosmicScore         = pset.get<float>       ("MaxCosmicScore",  0.4);
  fMaxLength              = pset.get<float>       ("MaxLength",      0.25);
  fMinCharge              = pset.get<float>       ("MinCharge",       75.);
  fMinStartDist           = pset.get<float>       ("MinStartDist",    10.);
}

DEFINE_ART_MODULE(SingleProtonFilter)

}
