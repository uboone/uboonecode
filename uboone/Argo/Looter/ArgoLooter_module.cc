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

#include "RecoBase/Wire.h"
#include "RawData/raw.h"
#include "RawData/RawDigit.h"
#include "Geometry/Geometry.h"

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
#include "ColorMap.h"
#include "MakePng.h"

#include "TSystem.h"

#include <fstream>
#include <unistd.h>
#include <string>

namespace microboone {

Timer gTimeStart;

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
  std::string fFileStoragePath;
  std::string fFileStorageUrl;
  
  // Output
  JsonObject fOutput;
  JsonObject fStats;

  int         source_fileEntries;
  std::string source_fileName;

  std::string source_selection;
  int         source_start;      
  int         source_end;      
  int         source_entry;    
  int         events_served;
  long        thisEventTimeStart;
  
  
  // Mapping from any adc value onto an 8-bit integer for crude color purposes.
  std::vector<unsigned char> fPalette;
  std::vector<unsigned char> fPaletteTrans;
  
  int inline static tanscale(float adc) 
  {
    return int(atan(adc/50.)/M_PI*256.) + 127;  
  }

  float inline static inv_tanscale(int y) 
  {
    return tan((y-127)*M_PI/256)*50.;
  }
  void wireOfChannel(int channel, int& plane, int& wire);
};


ArgoLooter::ArgoLooter(fhicl::ParameterSet const & p)
  : EDAnalyzer(p)
{
  
  JsonElement::SetPrettyPrint(false);
  reconfigure(p);
  source_fileEntries=-1;
  source_fileName = "(unknown)";
  fFileStoragePath = "../datacache";
  fFileStorageUrl = "datacache";
  
  fPalette.resize(256*3);
  fPaletteTrans.resize(256);

  unsigned char vv[] =   { 
    #include "palette.inc" 
  };
  fPalette.assign(&vv[0], &vv[0]+sizeof(vv));
  unsigned char vvt[] =   { 
    #include "palette_trans.inc" 
  };
  fPaletteTrans.assign(&vvt[0], &vvt[0]+sizeof(vvt));
  

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



  // Echo back the calling parameters, describe the source file.
  JsonObject source;
  // source.add("file",inRootFile);
  source.add("selection",source_selection);
  source.add("start",    source_start    );
  source.add("end",      source_end      );
  source.add("entry",    source_entry    );
  source.add("options",  fOptions);
  
  source.add("file",source_fileName);
  source.add("numEntriesInFile",source_fileEntries);
  fOutput.add("source",source);

  fOutput.add("converter","ComposeResult.cpp $Revision$ $Date$ ");

  // parse some options.
  int doCal = 1;
  int doRaw = 1;
  if( std::string::npos != fOptions.find("_NOCAL_")) doCal = 0;
  if( std::string::npos != fOptions.find("_NORAW_")) doRaw = 0;


  composeHeader(event);

  if(doRaw) { composeRaw(event); } else {composeRawAvailability(event);}
  if(doCal) { composeCal(event); } else {composeCalAvailability(event);}
  
  //reco
  composeHits(event);
  composeHits(event);
  composeClusters(event);
  composeVertex2d(event);
  composeSpacepoints(event);
  composeTracks(event);
  
  // Optical
  composeOpPulses(event);
  composeOpFlashes(event);
  composeOpHits(event);
  
  composeAuxDets(event);
  
  composeMC(event);
  

  JsonObject monitor;
  monitor.add("pid",getpid());
  monitor.add("ppid",getppid());
  monitor.add("events_served",events_served);
  SysInfo_t sysinfo;  gSystem->GetSysInfo(&sysinfo);
  // CpuInfo_t cpuinfo;  gSystem->GetCpuInfo(&cpuinfo);
  MemInfo_t meminfo;  gSystem->GetMemInfo(&meminfo);
  ProcInfo_t procinfo; gSystem->GetProcInfo(&procinfo);
  
  monitor.add("OS"           , sysinfo.fOS.Data());
  monitor.add("ComputerModel", sysinfo.fModel.Data());
  monitor.add("CpuType"      ,  sysinfo.fCpuType.Data());
  monitor.add("Cpus"         ,  sysinfo.fCpus);
  monitor.add("CpuSpeed"     ,  sysinfo.fCpuSpeed);
  monitor.add("PhysicalRam"  ,  sysinfo.fPhysRam);

  monitor.add("MemTotal"     ,  Form("%d MB",meminfo.fMemTotal));
  monitor.add("MemUsed"      ,  Form("%d MB",meminfo.fMemUsed));
  monitor.add("MemFree"      ,  Form("%d MB",meminfo.fMemFree));
  monitor.add("SwapTotal"    ,  Form("%d MB",meminfo.fSwapTotal));
  monitor.add("SwapUsed"     ,  Form("%d MB",meminfo.fSwapUsed));
  monitor.add("SwapFree"     ,  Form("%d MB",meminfo.fSwapFree));

  monitor.add("CpuTimeUser"     ,  procinfo.fCpuUser);
  monitor.add("CpuTimeSys"     ,   procinfo.fCpuSys);
  monitor.add("MemResident"    ,   Form("%f MB",procinfo.fMemResident/1000.));
  monitor.add("MemVirtual"     ,   Form("%f MB",procinfo.fMemVirtual/1000.));

  monitor.add("WallClockTime"  ,   gTimeStart.Count() );
  fOutput.add("backend_monitor",monitor);


  long ElapsedServerTime = ((long)(gSystem->Now()) - thisEventTimeStart);
  fOutput.add("ElapsedServerTime",ElapsedServerTime);
  std::cout << "ElapsedServerTime: " << ElapsedServerTime << std::endl;

  fOutput.add("stats",fStats);

  // Create an output file to match.
  // int thisEvent  = event.id().event();
  // int thisRun    = event.run();
  // int thisSubRun = event.subRun();
  // ofstream of(Form("event_%06d_%06d_%06d.json",thisRun,thisSubRun,thisEvent));
  // ofstream of(Form("event_%04d.json",thisEvent));
  // of << fOutput.str();
  // of.close();
  
  
  
  looterOutput() = fOutput.str();
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
  
void  ArgoLooter::composeClusters(art::Event const & event)
{}
  
void  ArgoLooter::composeVertex2d(art::Event const & event)
{}
  
void  ArgoLooter::composeSpacepoints(art::Event const & event)
{}
void  ArgoLooter::composeTracks(art::Event const & event)
{}

// Optical
void  ArgoLooter::composeOpFlashes(art::Event const & event)
{}
void  ArgoLooter::composeOpHits(art::Event const & event)
{}
void  ArgoLooter::composeOpPulses(art::Event const & event)
{}

// Wires
void  ArgoLooter::composeCalAvailability(art::Event const & event)
{
  JsonObject reco_list;

  // Get candidate list.
  typedef art::Handle< std::vector<recob::Wire> > wireHandle_t;
  std::vector<wireHandle_t> list;
  try{
    event.getManyByType(list);
  }catch(...) {
    cout << "problem getting types." << endl;
    return;
  }

  for(size_t i=0;i<list.size(); i++) {
    std::string name = list[i].provenance()->moduleLabel();
    reco_list.add(name,JsonElement());
  }
  fOutput.add("cal",reco_list);
  
}
  
  
void  ArgoLooter::composeRawAvailability(art::Event const & event)
{
  JsonObject reco_list;
  
  // Get candidate list.
  typedef art::Handle< std::vector<raw::RawDigit> > rawHandle_t;
  std::vector<rawHandle_t> list;
  try{
    event.getManyByType(list);
  }catch(...) {
    cout << "problem getting types." << endl;
    return;
  }

  for(size_t i=0;i<list.size(); i++) {
    std::string name = list[i].provenance()->moduleLabel();
    reco_list.add(name,JsonElement());

    // Remove from the list anything which looks like a prespill or a postspill window.
    if( std::string::npos != fOptions.find("_NoPreSpill_")) 
      if(std::string::npos != name.find("preSpill")) {
        reco_list.add(name,JsonElement()); // null it out.
        continue;
      }
    if( std::string::npos != fOptions.find("_NoPostSpill_"))
      if(std::string::npos != name.find("postSpill")) {
        reco_list.add(name,JsonElement()); // null it out.
        continue;
      }
    
    std::cout << "Looking at raw::RawDigit object " << name << endl;
    JsonObject r;
    ColorMap colormap;
      
      

  }
  fOutput.add("raw",reco_list);
    
}
  
void  ArgoLooter::composeCal(art::Event const & event)
{
  
}
  
  
  
void  ArgoLooter::composeRaw(art::Event const & event)
{
  JsonObject reco_list;
  
  // Get list of objects.
  typedef art::Handle< std::vector<raw::RawDigit> > rawHandle_t;
  std::vector<rawHandle_t> list;
  try{
    event.getManyByType(list);
  }catch(...) {
    cout << "problem getting types." << endl;
    return;
  }

  // art::ServiceHandle<geo::Geometry> geom;
  for(size_t i=0;i<list.size(); i++) {
    std::string name = list[i].provenance()->moduleLabel();
    TimeReporter timer(name);
    
    rawHandle_t handle(list[i]);
    const std::vector<raw::RawDigit>& vec(*handle);
    size_t ndig = vec.size();
    if(ndig==0) continue; // No digits.
    size_t width = vec[0].Samples();
    JsonObject r;
    ColorMap colormap;
    
    MakePng png(width,ndig, MakePng::palette_alpha,fPalette,fPaletteTrans);
    MakePng epng(width,ndig,MakePng::rgb);
    std::vector<unsigned char> imagedata(width);
    std::vector<unsigned char> encodeddata(width*3);
    TH1D timeProfile("timeProfile","timeProfile",width,0,width);
    std::vector<TH1*> planeProfile;
    std::vector<Double_t> timeProfileData(width+2,0);
    planeProfile.push_back(new TH1D("planeProfile0","planeProfile0",2398,0,2398));
    planeProfile.push_back(new TH1D("planeProfile1","planeProfile1",2398,0,2398));
    planeProfile.push_back(new TH1D("planeProfile2","planeProfile2",3456,0,3456));
    JsonArray jpedestals;

    // Get array by channel number
    std::vector<size_t> byWire(ndig);
    for(size_t i=0;i<ndig;i++) {
      unsigned int chan = vec[i].Channel();
      byWire[chan] = i;
    }
    
    std::vector<short> adcs(width);

    for(size_t c=0;c<byWire.size();c++) {
      size_t i = byWire[c];
      const raw::RawDigit& raw = vec[i];
      short pedestal = raw.GetPedestal();
      jpedestals.add(pedestal);
      
      raw::Compress_t compress = raw.Compression();
      if(compress!=raw::kNone)
        raw::Uncompress(raw.fADC, adcs, compress);

      double wiresum = 0;
      // if(adcs.size() < width) {
     //    std::cout << "Problem reading raw wires.. are they compressed?" << std::endl;
     //    JsonObject el;
     //    el.add("error","Error: adc size() less than width.");
     //    reco_list.add(name,el);
     //    goto continue_outer_loop;  // Yes, I used a goto. You wanna make something of it?
     //  }
      for(size_t k=0; k<width; k++) {
        short rawv;
        if(compress==raw::kNone) rawv =  raw.fADC[k];
        else                     rawv =  adcs[k]

        short pedcorr = rawv - pedestal;
        // colormap.get(&imagedata[k*3],float(raw)/4000.);
        imagedata[k] = tanscale(pedcorr);
      
        // Save bitpacked data as image map.
        int iadc = rawv + 0x8000;
        encodeddata[k*3]   = 0xFF&(iadc>>8);
        encodeddata[k*3+1] = iadc&0xFF;
        encodeddata[k*3+2] = 0;
        double val = fabs(pedcorr);
        wiresum += val;
        timeProfileData[k+1] += val;
      }
      png.AddRow(imagedata);
      epng.AddRow(encodeddata);

      // geo::WireID wid = geom->ChannelToWire(raw.Channel())[0];
      // int wire = wid.Wire;
      // int plane = wid.Plane;
      int wire, plane;
      wireOfChannel(i,plane,wire);
      planeProfile[plane]->Fill(wire,wiresum);      
    }
    timeProfile.SetContent(&timeProfileData[0]);
    png.Finish();
    epng.Finish();
    
    std::string wireimg = png.writeToUniqueFile(fFileStoragePath);
    std::string wireimg_thumb = wireimg+".thumb.png";
    BuildThumbnail(fFileStoragePath+wireimg,fFileStoragePath+wireimg_thumb);
    r.add("wireimg_url",fFileStorageUrl+wireimg);
    r.add("wireimg_url_thumb",fFileStorageUrl+wireimg_thumb);
    r.add("wireimg_encoded_url",fFileStorageUrl+
                              epng.writeToUniqueFile(fFileStoragePath)
                              );
    
    r.add("timeHist",TH1ToHistogram(&timeProfile));
    JsonArray jPlaneHists;
    jPlaneHists.add(TH1ToHistogram(planeProfile[0]));
    jPlaneHists.add(TH1ToHistogram(planeProfile[1]));
    jPlaneHists.add(TH1ToHistogram(planeProfile[2]));
    r.add("planeHists",jPlaneHists);
    r.add("pedestals",jpedestals);

    delete planeProfile[0];
    delete planeProfile[1];
    delete planeProfile[2];
    reco_list.add(name,r);
    fStats.add(name,timer.t.Count());

    // continue_outer_loop:
    // ; // no-op
  }
  
  fOutput.add("raw",reco_list);
  
  
}
  
// Other
void ArgoLooter::composeAuxDets(art::Event const & event)
{}

// Monte carlo
void  ArgoLooter::composeMC(art::Event const & event)
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
  fFileStoragePath  = p.get< std::string >("fileStoragePath");
  fFileStorageUrl   = p.get< std::string >("fileStorageUrl");
  fOptions         = p.get< std::string >("options");

  source_selection = p.get< std::string >("selection");
  source_start     = p.get< int >("entrystart");
  source_end       = p.get< int >("entryend");
  source_entry     = p.get< int >("entry");
  events_served    = p.get< int >("events_served");
  thisEventTimeStart = p.get< long >("thisEventTimeStart");
}

void ArgoLooter::wireOfChannel(int channel, int& plane, int& wire)
{
  if(channel < 2399) {
    plane = 0; wire= channel; return;
  }
  else if(channel <4798) {
    plane = 1; 
    wire = channel - 2399;
    return;
  }
  else{
    plane = 2;
    wire= channel-4798;
    return;
  }
}

DEFINE_ART_MODULE(ArgoLooter)
}