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
#include "art/Framework/Core/FileBlock.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"

#include <TH1D.h>
#include <TProfile.h>
#include <TTree.h>

#include "JsonElement.h"
#include "RootToJson.h"
#include "LooterGlobals.h"
#include "Timer.h"
#include <fstream>


namespace microboone {

class TimeReporter
{
public:
  std::string fName;
  Timer t;
  TimeReporter(const std::string& name="") :fName(name), t() {};
  ~TimeReporter() { std::cout << "++TimeReporter " << fName << " " << t.Count() << " s" << std::endl;}
};


class ArgoLooter : public art::EDAnalyzer {
public:
  explicit ArgoLooter(fhicl::ParameterSet const & p);
  virtual ~ArgoLooter();

  void analyze(art::Event const & e) override;

  void beginJob() override;
  void endJob() override;
  void reconfigure(fhicl::ParameterSet const & p) override;
  void respondToOpenInputFile(art::FileBlock const&);
  
protected:
  void  composeHeader(art::Event const & event);

  void  composeHits(art::Event const & event);
  void  composeClusters(art::Event const & event);
  void  composeVertex2d(art::Event const & event);
  void  composeSpacepoints(art::Event const & event);
  void  composeTracks(art::Event const & event);

  // Optical
  void  composeOpFlashes(art::Event const & event);
  void  composeOpHits(art::Event const & event);
  void  composeOpPulses(art::Event const & event);

  // Wires
  void  composeCalAvailability(art::Event const & event); 
  void  composeRawAvailability(art::Event const & event); 
  void  composeCal(art::Event const & event);
  void  composeRaw(art::Event const & event);

  // Other
  void composeAuxDets(art::Event const & event);
  
  // Monte carlo
  void  composeMC(art::Event const & event);


private:
  // art::ServiceHandle<geo::Geometry> fGeo;       // pointer to Geometry service
// >   art::ServiceHandle<util::DetectorProperties> fDetProp;      

  // Declare member data here.

  // Configuration 
  std::string fOptions;
  
  // Output
  JsonObject fOutput;
  JsonObject fStats;

  int         source_fileEntries;
  std::string source_fileName;
};


ArgoLooter::ArgoLooter(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{
  
  JsonElement::SetPrettyPrint(false);
  reconfigure(p);
  source_fileEntries=-1;
  source_fileName = "(unknown)";
}

ArgoLooter::~ArgoLooter()
{
  // Clean up dynamic memory and other resources here.
}

using std::endl;
using std::cout;


void ArgoLooter::analyze(art::Event const & event)
{
  
  std::cout << "Running ArgoLooter! Options:" << fOptions << std::endl;

  fOutput = JsonObject();
  fStats = JsonObject();
  composeHeader(event);
  composeHits(event);


  int thisEvent  = event.id().event(); 
  // int thisRun    = event.run();
  // int thisSubRun = event.subRun();

  JsonObject source;
  source.add("file",source_fileName);
  source.add("numEntriesInFile",source_fileEntries);
  fOutput.add("source",source);


  // Create an output file to match.
  // ofstream of(Form("event_%06d_%06d_%06d.json",thisRun,thisSubRun,thisEvent));
  ofstream of(Form("event_%04d.json",thisEvent));
  of << fOutput.str();
  of.close();
  
  looterOutput() = fOutput.str();
  std::cout << "Output: " << looterOutput() << std::endl;
  // Loot the wire info.

    
  
  // Implementation of required member function here.
}



void ArgoLooter::composeHeader(art::Event const & event) 
{
  // Header data. Alas, this is the stuff that's nearly impossible to get!  
  JsonObject header;
  header.add("run",   event.run() );
  header.add("subrun",event.subRun() );
  header.add("event", event.id().event() );

  // Todo: build these into a proper javascript-style timestamp.
  uint32_t tlow  = event.time().timeLow();
  uint32_t thigh = event.time().timeHigh();

  header.add("time",tlow);
  header.add("timeHigh",thigh);

  header.add("isRealData",    event.isRealData());
  header.add("experimentType",event.experimentType());
  
  // Add my own things. 
  // FIXME: this should come from the event data, not be hardcoded, but this will have to do for the moment.
  header.add("TDCStart",0);
  header.add("TDCEnd",9600);
  
  fOutput.add("header",header);
}


void ArgoLooter::composeHits(art::Event const & event)
{
  JsonObject reco_list;
  JsonObject hist_list;

  TimeReporter hittimer("hits");
  // Get list of hits.
  typedef art::Handle< std::vector<recob::Hit> > hitHandle_t;
  std::vector<hitHandle_t> list_of_hitlists;
  try{
    event.getManyByType(list_of_hitlists);
  }catch(...) {
    cout << "problem getting types." << endl;
    return;
  }

  cout << "Got " << list_of_hitlists.size() << " types of hits." << endl;
  for(size_t i=0;i<list_of_hitlists.size(); i++) {
    std::string name = list_of_hitlists[i].provenance()->moduleLabel();
    cout << "--->" << list_of_hitlists[i].provenance()->moduleLabel() << endl;
    TimeReporter timer(name);
    

    hitHandle_t hitHandle(list_of_hitlists[i]);
    const std::vector<recob::Hit>& hitvec(*hitHandle);
    
    // Hit histograms.
    TH1D timeProfile("timeProfile","timeProfile",960,0,9600);
    std::vector<TH1*> planeProfile;
    planeProfile.push_back(new TH1D("planeProfile0","planeProfile0",218,0,2398));
    planeProfile.push_back(new TH1D("planeProfile1","planeProfile1",218,0,2398));
    planeProfile.push_back(new TH1D("planeProfile2","planeProfile2",432,0,3456));
    
    std::vector<JsonObject> v;
    for(size_t i=0;i<hitvec.size();i++){
      const recob::Hit& hit = hitvec[i];
      int wire   = hit.WireID().Wire;
      int plane  = hit.WireID().Plane;
      double q   = hit.Charge();
      double t   = hit.PeakTime();
      double t1  = hit.StartTime();
      double t2  = hit.EndTime();

      JsonObject h;
      
      if(plane==2)timeProfile.Fill(t,q);
      if(plane>=0 && plane<3) planeProfile[plane]->Fill(wire,q);
      h.add("wire",    wire  );
      h.add("plane",   plane );
      h.add("q",       JsonFixed(q,0)     );
      h.add("t",       JsonFixed(t,1)     );
      h.add("t1",      JsonFixed(t1,1)    );
      h.add("t2",      JsonFixed(t2,1)    );
      // h.add("view",    ftr.getJson(lhit_view   ,i) ); // View is redundant with plane.
      // h.add("m",       ftr.getJson(lhit_m      ,i) );  // unusued
      // h.add("\u03C3q", ftr.getJson(lhit_sigq   ,i) ); // unused
      // h.add("\u03C3t", ftr.getJson(lhit_sigt   ,i) ); //unused
      v.push_back(h);                              
    }
   
    // Attempt to find association data. Look for a match to the name above, take the last match.
    // SHOULD be good enough.
    //vector<string> assnames = findLeafOfType("art::Wrapper<art::Assns<recob::Cluster,recob::Hit");
    // ...
          
    JsonArray arr;
    for(size_t i=0;i<v.size();i++) arr.add(v[i]);
    reco_list.add(name,arr);
    
    JsonObject hists;
    hists.add("timeHist",TH1ToHistogram(&timeProfile));
    JsonArray jPlaneHists;
    jPlaneHists.add(TH1ToHistogram(planeProfile[0]));
    jPlaneHists.add(TH1ToHistogram(planeProfile[1]));
    jPlaneHists.add(TH1ToHistogram(planeProfile[2]));
    hists.add("planeHists",jPlaneHists);

    delete planeProfile[0];
    delete planeProfile[1];
    delete planeProfile[2];    
    hist_list.add(name,hists);
    fStats.add(name,timer.t.Count());
  }
  fOutput.add("hits",reco_list);
  fOutput.add("hit_hists",hist_list);
  fStats.add("hits",hittimer.t.Count());
  
}
  
void  composeClusters(art::Event const & event)
{}
  
void  composeVertex2d(art::Event const & event)
{}
  
void  composeSpacepoints(art::Event const & event)
{}
void  composeTracks(art::Event const & event)
{}

// Optical
void  composeOpFlashes(art::Event const & event)
{}
void  composeOpHits(art::Event const & event)
{}
void  composeOpPulses(art::Event const & event)
{}

// Wires
void  composeCalAvailability(art::Event const & event)
{}
void  composeRawAvailability(art::Event const & event)
{}
void  composeCal(art::Event const & event)
{}
void  composeRaw(art::Event const & event)
{}
  
// Other
void composeAuxDets(art::Event const & event)
{}

// Monte carlo
void  composeMC(art::Event const & event)
{}


void ArgoLooter::respondToOpenInputFile(art::FileBlock const& fb)
{
  if(fb.tree()) source_fileEntries=fb.tree()->GetEntriesFast();
  source_fileName = fb.fileName();
  
}


void ArgoLooter::beginJob()
{
  // Implementation of optional member function here.
}

void ArgoLooter::endJob()
{
  // Implementation of optional member function here.
}

void ArgoLooter::reconfigure(fhicl::ParameterSet const & p)
{
  fOptions = p.get< std::string >("options");
  // Implementation of optional member function here.
}



DEFINE_ART_MODULE(ArgoLooter)
}