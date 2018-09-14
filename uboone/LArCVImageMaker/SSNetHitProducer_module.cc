////////////////////////////////////////////////////////////////////////
// Class:       SSNetHitProducer
// Plugin Type: producer (art v2_05_00)
// File:        SSNetHitProducer_module.cc
//
// Generated at Fri Nov 10 20:09:52 2017 by David Caratelli using cetskelgen
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

#include <memory>

#include "lardataobj/RecoBase/Hit.h"

#include "LArCVMetaMaker.h"
#include "LArCVSuperaDriver.h"
#include "Base/PSet.h"

#include "DataFormat/EventImage2D.h"
#include "DataFormat/IOManager.h"

#include "TFile.h"
#include "TTree.h"

class SSNetHitProducer;


class SSNetHitProducer : public art::EDProducer {
public:
  explicit SSNetHitProducer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SSNetHitProducer(SSNetHitProducer const &) = delete;
  SSNetHitProducer(SSNetHitProducer &&) = delete;
  SSNetHitProducer & operator = (SSNetHitProducer const &) = delete;
  SSNetHitProducer & operator = (SSNetHitProducer &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  larcv::IOManager *mgr;

  // location of LArCV file
  std::string fLArCVLocation;
  // pixel score threshold
  double fPxThresholdHigh, fPxThresholdLow;
  // in hit producer
  std::string fInHitProducer;
  // which hits to save? track (1) or shower (0)
  int fPxType;

};


SSNetHitProducer::SSNetHitProducer(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{


  fLArCVLocation   = p.get<std::string>("LArCVLocation"          );
  fPxType          = p.get<int>        ("PxType"       ,0        );
  fInHitProducer   = p.get<std::string>("InHitProducer","gaushit");
  fPxThresholdHigh = p.get<double>     ("PxThresholdHigh"        );
  fPxThresholdLow  = p.get<double>     ("PxThresholdLow"         );
  
  produces< std::vector< recob::Hit > >(Form("shrhit%zu",(size_t)(fPxThresholdHigh*100.)));
  produces< std::vector< recob::Hit > >(Form("shrhit%zu",(size_t)(fPxThresholdLow*100.) ));
  
}


void SSNetHitProducer::produce(art::Event & e)
{

  // load already existing hits
  art::Handle<std::vector<recob::Hit> > hit_h;
  e.getByLabel(fInHitProducer,hit_h);

  int currentevt = (int)e.event();

  // make sure flash look good
  if(!hit_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate input hit!"<<std::endl;
    throw std::exception();
  }

  // vector of various hit producers to store
  std::unique_ptr< std::vector<recob::Hit> >  Hit_High_v(new std::vector<recob::Hit>);
  std::unique_ptr< std::vector<recob::Hit> >  Hit_Low_v(new std::vector<recob::Hit>);

  auto nentries = mgr->get_n_entries();

  bool FOUND = false;

  // loop through all entries in LArCV file and find 
  for (size_t i=0; i < nentries; i++) {

    if (FOUND == true) break;
      
      mgr->read_entry(i);
      
      // get LArCV event information
      auto evtimage = (larcv::EventImage2D*)mgr->get_data(larcv::kProductImage2D,"wire");
      if (!evtimage) { std::cout << "NO EventImage2D !!!" << std::endl; break; }
      int evtnum = (int)evtimage->event();

      std::cout << "CURRENT evt : " << currentevt << "\t LArCV evt : " << evtnum << std::endl;
      
      // if event matches across files
      if ( currentevt == evtnum) {

	FOUND = true;

	// loop through planes
	for(size_t plane=0; plane<3; ++plane) {

	  auto ev_image = (larcv::EventImage2D*)mgr->get_data(larcv::kProductImage2D,Form("uburn_plane%zu",plane));
	  if (!ev_image) { std::cout << "NO EventImage2D !!!" << std::endl; break; }
	  
	  auto image2D = ev_image->Image2DArray();
	  
	  // grab enry 0 for tracks, entry 1 for showers
	  const auto& seg_img = image2D.at(fPxType);
	  const auto& meta    = image2D.at(fPxType).meta();
	  
	  int rows = (int)meta.rows();
	  int cols = (int)meta.cols();
	  
	  for (auto const& hit : *hit_h){
	    
	    if (hit.WireID().Plane != plane) continue;
	    
	    auto time = hit.PeakTime() + 2400;
	    auto wire = hit.WireID().Wire;
	    
	    auto xpixel = (wire - meta.min_x()) / meta.pixel_width();
	    auto ypixel = (meta.max_y() - time) / meta.pixel_height();
	    
	    int xx = (int)(xpixel+0.5);
	    int yy = (int)(ypixel+0.5);
	    
	    if ( (xx < 0) || (yy < 0) ) continue;
	    if ( (xx >= cols) || (yy >= rows) ) continue;
	    
	    float pixel_type  = seg_img.pixel(yy,xx);
	    
	    if ( pixel_type > fPxThresholdHigh ) 
	      Hit_High_v->emplace_back(hit);

	    if ( pixel_type > fPxThresholdLow ) 
	      Hit_Low_v->emplace_back(hit);
	    
	  }// for all hits
	} // end this plane 
      }// if event from LArCV file and artroot files match
  }// for all entries in LArCV file

  e.put(std::move(Hit_High_v),Form("shrhit%zu",(size_t)(fPxThresholdHigh*100.)));  
  e.put(std::move(Hit_Low_v) ,Form("shrhit%zu",(size_t)(fPxThresholdLow*100.) ));  

}

void SSNetHitProducer::beginJob()
{

  mgr = new larcv::IOManager();
  mgr->add_in_file(fLArCVLocation);
  mgr->initialize();

}

void SSNetHitProducer::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(SSNetHitProducer)
