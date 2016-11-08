#include "uboone/CRTDetSim/CRTDetSim.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larcore/Geometry/AuxDetGeometry.h"

#include "ucorsika/Overlays/BernFEBFragment.hh"
#include "artdaq-core/Data/Fragment.hh"

#include <vector>
#include <memory>
#include <sstream>
#include <string>


namespace crt{

  std::string __log_name__ = "CRTDetSim";

  CRTDetSim::CRTDetSim(const fhicl::ParameterSet& pSet){
    produces< bernfebdaq::BernFEBFragmentMetadata >();
    produces< bernfebdaq::BernFEBEvent >();
    mf::LogInfo(__log_name__)<<"In construction: ";
  }

  CRTDetSim::~CRTDetSim(){
    mf::LogInfo(__log_name__)<<"In destruction: ";
  }

  void CRTDetSim::produce(art::Event& evt){

    std::unique_ptr< std::vector<bernfebdaq::BernFEBFragmentMetadata> > meta;
      //( new bernfebdaq::BernFEBFragmentMetadata(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0));
    std::unique_ptr< std::vector<bernfebdaq::BernFEBEvent> > febevent;
    //artdaq::Fragment artfrag;
    //create a new fragment instance
    //bernfebdaq::BernFEBFragment frag(artfrag);
    mf::LogInfo(__log_name__)<<"In produce ";
    //Get the geometry to figure out the channel mapping
    art::ServiceHandle<geo::AuxDetGeometry> geo;

    //houses all of the channels hit by the event
    art::Handle< std::vector<sim::AuxDetSimChannel> > channels;
    evt.getByLabel("",channels);
    if(!channels.isValid()) return;

    mf::LogInfo(__log_name__)<<" Number of Channels Hit: "<<channels->size();

    for(auto it = channels->begin(); it!= channels->end(); ++it){
      uint32_t id = it->AuxDetID();
      uint32_t sens_id = it->AuxDetSensitiveID();
      mf::LogInfo(__log_name__)<<"Found AuxDetData: "<<id<<" , "<<sens_id;
      std::vector< sim::AuxDetIDE > ides = it->AuxDetIDEs();

      for(auto ideIt = ides.begin(); ideIt!= ides.end(); ++it){
        int trackID = ideIt->trackID;
        float energyDep = ideIt->energyDeposited;

        float entryX = ideIt->entryX;
        float entryY = ideIt->entryY;
        float entryZ = ideIt->entryZ;
        float entryT = ideIt->entryT;
        //
        //
        //

        float exitX = ideIt->exitX;
        float exitY = ideIt->exitY;
        float exitZ = ideIt->exitZ;
        float exitT = ideIt->exitT;

        float exitMomX = ideIt->exitMomentumX;
        float exitMomY = ideIt->exitMomentumY;
        float exitMomZ = ideIt->exitMomentumZ;

      }
    }
    evt.put(std::move(meta));
    evt.put(std::move(febevent));
  }

  DEFINE_ART_MODULE(CRTDetSim)
}
