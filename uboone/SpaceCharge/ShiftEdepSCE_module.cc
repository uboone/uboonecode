////////////////////////////////////////////////////////////////////////
// Class:       ShiftEdepSCE
// Plugin Type: producer (art v2_05_01)
// File:        ShiftEdepSCE_module.cc
//
// Generated at Thu Apr 19 00:41:18 2018 by Wesley Ketchum using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include <memory>

#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "TNtuple.h"

namespace spacecharge {
  class ShiftEdepSCE;
}


class spacecharge::ShiftEdepSCE : public art::EDProducer {
public:
  explicit ShiftEdepSCE(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ShiftEdepSCE(ShiftEdepSCE const &) = delete;
  ShiftEdepSCE(ShiftEdepSCE &&) = delete;
  ShiftEdepSCE & operator = (ShiftEdepSCE const &) = delete;
  ShiftEdepSCE & operator = (ShiftEdepSCE &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;
  void beginJob() override;

private:

  // Declare member data here.
  art::InputTag fEDepTag;
  bool          fMakeAnaTree;
  TNtuple*      fNtEdepAna;

};


spacecharge::ShiftEdepSCE::ShiftEdepSCE(fhicl::ParameterSet const & p)
  : fEDepTag(p.get<art::InputTag>("EDepTag")),
    fMakeAnaTree(p.get<bool>("MakeAnaTree",true))
{
  produces< std::vector<sim::SimEnergyDeposit> >();
}

void spacecharge::ShiftEdepSCE::beginJob()
{
  if(fMakeAnaTree){
    art::ServiceHandle<art::TFileService> tfs;
    fNtEdepAna = tfs->make<TNtuple>("nt_edep_ana","Edep PosDiff Ana Ntuple","orig_x:orig_y:orig_z:shift_x:shift_y:shift_z");
  }
}

void spacecharge::ShiftEdepSCE::produce(art::Event & e)
{
  auto sce = lar::providerFrom<spacecharge::SpaceChargeService>();

  std::unique_ptr< std::vector<sim::SimEnergyDeposit> > 
    outEdepVecPtr(new std::vector<sim::SimEnergyDeposit>() );
  auto & outEdepVec(*outEdepVecPtr);

  art::Handle< std::vector<sim::SimEnergyDeposit> > inEdepHandle;
  e.getByLabel(fEDepTag,inEdepHandle);
  auto const& inEdepVec(*inEdepHandle);
  
  outEdepVec.reserve(inEdepVec.size());
  
  std::vector<double> posOffsets{0.0,0.0,0.0};
  for(auto const& edep : inEdepVec){
    if(sce->EnableSimSpatialSCE())
      posOffsets = sce->GetPosOffsets(edep.X(),edep.Y(),edep.Z());
    outEdepVec.emplace_back(edep.NumPhotons(),
			    edep.NumElectrons(),
			    edep.Energy(),
			    sim::SimEnergyDeposit::Point_t{(float)(edep.StartX()+posOffsets[0]),
				(float)(edep.StartY()+posOffsets[1]),
				(float)(edep.StartZ()+posOffsets[2])},
			    sim::SimEnergyDeposit::Point_t{(float)(edep.EndX()+posOffsets[0]),
				(float)(edep.EndY()+posOffsets[1]),
				(float)(edep.EndZ()+posOffsets[2])},
			    edep.StartT(),
			    edep.EndT(),
			    edep.TrackID(),
			    edep.PdgCode());
    if(fMakeAnaTree)
      fNtEdepAna->Fill(edep.X(),edep.Y(),edep.Z(),
		       outEdepVec.back().X(),outEdepVec.back().Y(),outEdepVec.back().Z());
  }

  e.put(std::move(outEdepVecPtr));
}

DEFINE_ART_MODULE(spacecharge::ShiftEdepSCE)
