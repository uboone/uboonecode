
#include "uboone/AnalysisTree/MCTruth/IMCTruthMatching.h"

#include "uboone/AnalysisTree/MCTruth/MCTruthBase/MCTruthAssociations.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "nutools/ParticleNavigation/ParticleList.h"

#include <cmath>
#include <algorithm>

namespace truth
{
////////////////////////////////////////////////////////////////////////
//
// Class:       AssociationsTruth
// Module Type: art tool
// File:        AssociationsTruth.h
//
//              This provides MC truth information by using output
//              reco Hit <--> MCParticle associations
//
// Configuration parameters:
//
// TruncMeanFraction     - the fraction of waveform bins to discard when
//
// Created by Tracy Usher (usher@slac.stanford.edu) on November 21, 2017
//
////////////////////////////////////////////////////////////////////////

class AssociationsTruth : virtual public IMCTruthMatching
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  pset
     */
    explicit AssociationsTruth(fhicl::ParameterSet const & pset);
    
    /**
     *  @brief  Destructor
     */
    ~AssociationsTruth();
    
    // provide for initialization
    void reconfigure(fhicl::ParameterSet const & pset) override;
    
    /**
     *  @brief This rebuilds the internal maps
     */
    void Rebuild(const art::Event& evt) override;

    /**
     *  @brief Get a reference to the ParticleList
     */
    const sim::ParticleList& ParticleList() const override;
    
    // Set the EveIdCalculator for the owned ParticleList
//    void  SetEveIdCalculator(MCTruthEveIdCalculator *ec) { fParticleList.AdoptEveIdCalculator(ec); }
    
    // Return a pointer to the simb::MCParticle object corresponding to
    // the given TrackID
    const simb::MCParticle* TrackIDToParticle(int const& id)       const override {return fMCTruthAssociations.TrackIDToParticle(id);}
    const simb::MCParticle* TrackIDToMotherParticle(int const& id) const override {return fMCTruthAssociations.TrackIDToMotherParticle(id);}
    
    // Get art::Ptr<> to simb::MCTruth and related information
    const art::Ptr<simb::MCTruth>&                TrackIDToMCTruth(int const& id)                        const override;
    const art::Ptr<simb::MCTruth>&                ParticleToMCTruth(const simb::MCParticle* p)           const override;
    std::vector<const simb::MCParticle*>          MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const override;
    const std::vector< art::Ptr<simb::MCTruth> >& MCTruthVector()                                        const override;
    
    // this method will return the Geant4 track IDs of
    // the particles contributing ionization electrons to the identified hit
    std::vector<sim::TrackIDE> HitToTrackID(recob::Hit const& hit)           const override;
    std::vector<sim::TrackIDE> HitToTrackID(art::Ptr<recob::Hit> const& hit) const override;
    
    // method to return a subset of allhits that are matched to a list of TrackIDs
    const std::vector<std::vector<art::Ptr<recob::Hit>>> TrackIDsToHits(std::vector<art::Ptr<recob::Hit>> const& allhits,
                                                                                std::vector<int> const& tkIDs) const override;
    
    // method to return the EveIDs of particles contributing ionization
    // electrons to the identified hit
    std::vector<sim::TrackIDE> HitToEveID(art::Ptr<recob::Hit> const& hit) const override;
    
    // method to return the XYZ position of the weighted average energy deposition for a given hit
    std::vector<double>  HitToXYZ(art::Ptr<recob::Hit> const& hit) const override;
    
    // method to return the XYZ position of a space point (unweighted average XYZ of component hits).
    std::vector<double> SpacePointToXYZ(art::Ptr<recob::SpacePoint> const& spt,
                                        art::Event                  const& evt,
                                        std::string                 const& label) const override;
    
    // method to return the XYZ position of a space point (unweighted average XYZ of component hits).
    std::vector<double> SpacePointHitsToXYZ(art::PtrVector<recob::Hit> const& hits) const override;
    
    // method to return the fraction of hits in a collection that come from the specified Geant4 track ids
    double HitCollectionPurity(std::set<int>                              trackIDs,
                                       std::vector< art::Ptr<recob::Hit> > const& hits) const override;
    
    // method to return the fraction of all hits in an event from a specific set of Geant4 track IDs that are
    // represented in a collection of hits
    double HitCollectionEfficiency(std::set<int>                              trackIDs,
                                   std::vector< art::Ptr<recob::Hit> > const& hits,
                                   std::vector< art::Ptr<recob::Hit> > const& allhits,
                                   geo::View_t                         const& view) const override;
    
    // method to return the fraction of charge in a collection that come from the specified Geant4 track ids
    double HitChargeCollectionPurity(std::set<int>                              trackIDs,
                                             std::vector< art::Ptr<recob::Hit> > const& hits) const override;
    
    // method to return the fraction of all charge in an event from a specific set of Geant4 track IDs that are
    // represented in a collection of hits
    double HitChargeCollectionEfficiency(std::set<int>                              trackIDs,
                                         std::vector< art::Ptr<recob::Hit> > const& hits,
                                         std::vector< art::Ptr<recob::Hit> > const& allhits,
                                         geo::View_t                         const& view) const override;
    
    // method to return all EveIDs corresponding to the current sim::ParticleList
    std::set<int> GetSetOfEveIDs() const override;
    
    // method to return all TrackIDs corresponding to the current sim::ParticleList
    std::set<int> GetSetOfTrackIDs() const override;
    
    // method to return all EveIDs corresponding to the given list of hits
    std::set<int> GetSetOfEveIDs(std::vector< art::Ptr<recob::Hit> > const& hits) const override;
    
    // method to return all TrackIDs corresponding to the given list of hits
    std::set<int> GetSetOfTrackIDs(std::vector< art::Ptr<recob::Hit> > const& hits) const override;

private:
    
    // Fcl parameters.
    std::vector<art::InputTag>         fAssnsProducerLabels;  ///< tag for finding the tracks
    art::InputTag                      fG4ProducerLabel;      ///< Input tag for G4 producer (MCParticle/MCTruth)
    
    // The class that does all the work...
    MCTruthAssociations                fMCTruthAssociations;  ///< The class that does the work
    
    // Hopefully we can get rid of this soon
    sim::ParticleList                  fParticleList;
    
    // Useful services, keep copies for now (we can update during begin run periods)
    const geo::GeometryCore*           fGeometry;             ///< pointer to Geometry service
    const detinfo::DetectorProperties* fDetectorProperties;   ///< Detector properties service
};
    
//----------------------------------------------------------------------------
/// Constructor.
///
/// Arguments:
///
/// pset - Fcl parameters.
///
AssociationsTruth::AssociationsTruth(fhicl::ParameterSet const & pset) :
    fMCTruthAssociations(pset.get<fhicl::ParameterSet>("MCTruthAssociations"))
{
    fGeometry           = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    reconfigure(pset);
    
    // Report.
    mf::LogInfo("AssociationsTruth") << "AssociationsTruth configured\n";
}

//----------------------------------------------------------------------------
/// Destructor.
AssociationsTruth::~AssociationsTruth()
{}

//----------------------------------------------------------------------------
/// Reconfigure method.
///
/// Arguments:
///
/// pset - Fcl parameter set.
///
void AssociationsTruth::reconfigure(fhicl::ParameterSet const & pset)
{
    fAssnsProducerLabels = pset.get<std::vector<art::InputTag>>("AssnsProducerLabels");
    fG4ProducerLabel     = pset.get<art::InputTag>             ("G4ProducerLabel");
}

//----------------------------------------------------------------------------
/// Rebuild method -> rebuild the basic maps to get truth information
///
/// Arguments:
///
/// event - the art event used to extract all information
///
void AssociationsTruth::Rebuild(const art::Event& evt)
{
    // Create a container for testing
    HitParticleAssociationsVec partHitAssnsVec;

    // Get a handle for the associations...
    art::Handle<art::Assns<simb::MCParticle, recob::Hit, anab::BackTrackerHitMatchingData>> partHitAssnsHandle;
    
    for(const auto& assnsProducerLabel : fAssnsProducerLabels)
    {
        evt.getByLabel(assnsProducerLabel, partHitAssnsHandle);
    
        if (!partHitAssnsHandle.isValid())
        {
            throw cet::exception("AssociationsTruth") << "===>> NO MCParticle <--> Hit associations found for run/subrun/event: " << evt.run() << "/" << evt.subRun() << "/" << evt.id().event() << std::endl;
        }
    
        partHitAssnsVec.emplace_back(&*partHitAssnsHandle);
    }
    
    // Recover the associations between MCTruth and MCParticles
    art::Handle<MCTruthParticleAssociations> truthPartAssnsHandle;
    evt.getByLabel(fG4ProducerLabel, truthPartAssnsHandle);
    
    // Pass this to the truth associations code
    fMCTruthAssociations.setup(partHitAssnsVec, *truthPartAssnsHandle, *fGeometry, *fDetectorProperties);
    
    // Ugliness to follow! Basically, we need to build the "particle list" and the current implementation of
    // that code requires a copy...
    const MCTruthParticleList& locParticleList = fMCTruthAssociations.getParticleList();
    
    // Clear the current container
    fParticleList.clear();
    
    // Now we add particles back in one at a time...
    for(const auto& element : locParticleList)
    {
        fParticleList.Add(new simb::MCParticle(*(element.second)));
    }
}

//----------------------------------------------------------------------
const sim::ParticleList& AssociationsTruth::ParticleList() const
{
    // Unfortunately, this requires special handling at the moment...
    return fParticleList;
}
    
//----------------------------------------------------------------------
const art::Ptr<simb::MCTruth>& AssociationsTruth::TrackIDToMCTruth(int const& id) const
{
    return fMCTruthAssociations.TrackIDToMCTruth(id);
}
    
//----------------------------------------------------------------------
const art::Ptr<simb::MCTruth>& AssociationsTruth::ParticleToMCTruth(const simb::MCParticle* p) const
{
    return fMCTruthAssociations.ParticleToMCTruth(p);
}
    
//----------------------------------------------------------------------
std::vector<const simb::MCParticle*> AssociationsTruth::MCTruthToParticles(art::Ptr<simb::MCTruth> const& mct) const
{
    return fMCTruthAssociations.MCTruthToParticles(mct);
}
    
//----------------------------------------------------------------------
const std::vector< art::Ptr<simb::MCTruth> >&  AssociationsTruth::MCTruthVector() const
{
    return fMCTruthAssociations.MCTruthVector();
}

//----------------------------------------------------------------------
std::vector<sim::TrackIDE> AssociationsTruth::HitToTrackID(recob::Hit const& hit) const
{
    return fMCTruthAssociations.HitToTrackID(&hit);
}

//----------------------------------------------------------------------
std::vector<sim::TrackIDE> AssociationsTruth::HitToTrackID(art::Ptr<recob::Hit> const& hit) const
{
    return fMCTruthAssociations.HitToTrackID(hit);
}

//----------------------------------------------------------------------
const std::vector<std::vector<art::Ptr<recob::Hit>>> AssociationsTruth::TrackIDsToHits(std::vector<art::Ptr<recob::Hit>> const& allhits,
                                                                                       std::vector<int> const& tkIDs) const
{
    return fMCTruthAssociations.TrackIDsToHits(allhits,tkIDs);
}
    
//----------------------------------------------------------------------
// plist is assumed to have adopted the appropriate EveIdCalculator prior to
// having been passed to this method. It is likely that the EmEveIdCalculator is
// the one you always want to use
std::vector<sim::TrackIDE> AssociationsTruth::HitToEveID(art::Ptr<recob::Hit> const& hit) const
{
    return fMCTruthAssociations.HitToEveID(hit);
}
    
//----------------------------------------------------------------------
std::set<int> AssociationsTruth::GetSetOfEveIDs() const
{
    return fMCTruthAssociations.GetSetOfEveIDs();
}
    
//----------------------------------------------------------------------
std::set<int> AssociationsTruth::GetSetOfTrackIDs() const
{
    return fMCTruthAssociations.GetSetOfTrackIDs();
}
    
//----------------------------------------------------------------------
std::set<int> AssociationsTruth::GetSetOfEveIDs(std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    return fMCTruthAssociations.GetSetOfEveIDs(hits);
}
    
//----------------------------------------------------------------------
std::set<int> AssociationsTruth::GetSetOfTrackIDs(std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    return fMCTruthAssociations.GetSetOfTrackIDs(hits);
}
    
//----------------------------------------------------------------------
double AssociationsTruth::HitCollectionPurity(std::set<int> trackIDs, std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    return fMCTruthAssociations.HitCollectionPurity(trackIDs, hits);
}
    
//----------------------------------------------------------------------
double AssociationsTruth::HitChargeCollectionPurity(std::set<int> trackIDs, std::vector< art::Ptr<recob::Hit> > const& hits) const
{
    return fMCTruthAssociations.HitChargeCollectionPurity(trackIDs, hits);
}
    
    
//----------------------------------------------------------------------
double AssociationsTruth::HitCollectionEfficiency(std::set<int>                              trackIDs,
                                                  std::vector< art::Ptr<recob::Hit> > const& hits,
                                                  std::vector< art::Ptr<recob::Hit> > const& allhits,
                                                  geo::View_t const&                         view) const
{
    return fMCTruthAssociations.HitCollectionEfficiency(trackIDs, hits, allhits, view);
}
    
//----------------------------------------------------------------------
double AssociationsTruth::HitChargeCollectionEfficiency(std::set<int>                              trackIDs,
                                                        std::vector< art::Ptr<recob::Hit> > const& hits,
                                                        std::vector< art::Ptr<recob::Hit> > const& allhits,
                                                        geo::View_t                         const& view) const
{
    return fMCTruthAssociations.HitChargeCollectionEfficiency(trackIDs, hits, allhits, view);
}
    
//----------------------------------------------------------------------
std::vector<double> AssociationsTruth::HitToXYZ(art::Ptr<recob::Hit> const& hit) const
{
    return fMCTruthAssociations.HitToXYZ(hit);
}
    
//----------------------------------------------------------------------
std::vector<double> AssociationsTruth::SpacePointToXYZ(art::Ptr<recob::SpacePoint> const& spt,
                                                       art::Event                  const& evt,
                                                       std::string                 const& label) const
{
    // Get hits that make up this space point.
    art::PtrVector<recob::SpacePoint> spv;
    spv.push_back(spt);
    art::FindManyP<recob::Hit> fmh(spv, evt, label);
    std::vector< art::Ptr<recob::Hit> > hitv = fmh.at(0);
    
    // make a PtrVector
    art::PtrVector<recob::Hit> hits;
    for(size_t h = 0; h < hitv.size(); ++h) hits.push_back(hitv[h]);
    
    return this->SpacePointHitsToXYZ(hits);
}
    
//----------------------------------------------------------------------
std::vector<double> AssociationsTruth::SpacePointHitsToXYZ(art::PtrVector<recob::Hit> const& hits) const
{
    return fMCTruthAssociations.SpacePointHitsToXYZ(hits);
}

//----------------------------------------------------------------------------
    
DEFINE_ART_CLASS_TOOL(AssociationsTruth)
}
