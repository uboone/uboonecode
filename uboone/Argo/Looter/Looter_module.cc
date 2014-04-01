////////////////////////////////////////////////////////////////////////
// Class:       Looter
// Module Type: analyzer
// File:        Looter_module.cc
//
// Generated at Wed Jan 30 14:54:31 2013 by Nathaniel Tagg using artmod
// from art v1_02_06.
////////////////////////////////////////////////////////////////////////

// LArSoft includes
#include "Simulation/SimChannel.h"
#include "Simulation/LArG4Parameters.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
// #include "Geometry/Geometry.h"
// #include "Utilities/DetectorProperties.h"
// #include "SimulationBase/MCParticle.h"
// #include "SimulationBase/MCTruth.h"
#include "SimpleTypesAndConstants/geo_types.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

#include "JsonElement.h"
#include <fstream>

namespace microboone {

class Looter : public art::EDAnalyzer {
public:
  explicit Looter(fhicl::ParameterSet const & p);
  virtual ~Looter();

  void analyze(art::Event const & e) override;

  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;

private:
  // art::ServiceHandle<geo::Geometry> fGeo;       // pointer to Geometry service
// >   art::ServiceHandle<util::DetectorProperties> fDetProp;      

  // Declare member data here.

};


Looter::Looter(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{
  JsonElement::SetPrettyPrint(true);
}

Looter::~Looter()
{
  // Clean up dynamic memory and other resources here.
}

using std::endl;
using std::cout;


void Looter::analyze(art::Event const & event)
{
  
  int thisEvent  = event.id().event(); 
  int thisRun    = event.run();
  int thisSubRun = event.subRun();

  // Create an output file to match.
  // ofstream of(Form("event_%06d_%06d_%06d.json",thisRun,thisSubRun,thisEvent));
  ofstream of(Form("event_%04d.json",thisEvent));

  JsonObject root;
  JsonObject header;
  header.add("run",thisRun);
  header.add("subrun",thisSubRun);
  header.add("event",thisEvent);
  root.add("header",header);

  // Get list of hits.
  typedef art::Handle< std::vector<recob::Hit> > hitHandle_t;
  std::vector<hitHandle_t> list_of_hitlists;
  try{
    event.getManyByType(list_of_hitlists);
  }catch(...) {
    cout << "problem getting types." << endl;
  }

  cout << "Got " << list_of_hitlists.size() << " types of hits." << endl;
  for(size_t i=0;i<list_of_hitlists.size(); i++) {
    cout << "--->" << list_of_hitlists[i].provenance()->moduleLabel() << endl;
  }


  hitHandle_t hitHandle;
  try {
    if(!event.getByLabel("rffhit", hitHandle)) return;
  }
  catch (...){
      return;
  }
  JsonArray hits;
  // For every Hit:
  // for ( auto const& hit : (*hitHandle) ) {
  const std::vector<recob::Hit>& hitvec(*hitHandle);
  for(size_t i=0;i<hitvec.size();i++){
    const recob::Hit& hit = hitvec[i];
    JsonObject jhit;
    jhit.add("t1",hit.StartTime()                );
    jhit.add("t2",hit.EndTime()                  );
    jhit.add("t",hit.PeakTime()                 );
    jhit.add("σt1",hit.SigmaStartTime()           );
    jhit.add("σt2",hit.SigmaEndTime()             );
    jhit.add("σt",hit.SigmaPeakTime()            );
    jhit.add("n",hit.Multiplicity()             );
    jhit.add("chan",hit.Channel()                  );
    jhit.add("q",hit.Charge()                   );
    jhit.add("σq",hit.SigmaCharge()              );
    jhit.add("gof",hit.GoodnessOfFit()            );
    jhit.add("sigtype",hit.SignalType());
    jhit.add("view",hit.View());

    // std::vector<geo::WireID> wireids = fGeo->ChannelToWire( hit.Channel() );
    //const geo::WireID& wireid = wireids[0];
    // JsonObject jwireid;
//     jwireid.add("Cryostat",wireid.Cryostat);
//     jwireid.add("TPC",wireid.TPC);
//     jwireid.add("Plane",wireid.Plane);
//     jwireid.add("Wire",wireid.Wire);
//     jhit.add("WireID",jwireid);
    hits.add(jhit);
  }
  root.add("hits",hits);

  of << root.str();
  of.close();
  
  
  // Loot the wire info.

    
  
  // Implementation of required member function here.
}

void Looter::beginJob()
{
  // Implementation of optional member function here.
}

void Looter::endJob()
{
  // Implementation of optional member function here.
}

void Looter::reconfigure(fhicl::ParameterSet const & p)
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(Looter)
}