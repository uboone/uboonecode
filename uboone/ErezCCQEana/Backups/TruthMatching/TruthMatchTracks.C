/*************************************************************
 * 
 * This is a quick demonstrator macro for performing track/particle
 * matching using the new hit/particle associations.
 *
 * root [0] .L TruthMatchTracks.C++
 * root [1] TruthMatchTracks("input_file.root","track_tag","hit_tag","hitparticleassns_tag")
 *
 * Wesley Ketchum (wketchum@fnal.gov), Nov17, 2017
 * 
 *************************************************************/


//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TStopwatch.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

//"larsoft" object includes
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

//I like doing this to not get fooled by underflow/overflow
void ShowUnderOverFlow(TH1* h1){
  h1->SetBinContent(1, h1->GetBinContent(0)+h1->GetBinContent(1));
  h1->SetBinContent(0,0);

  int nbins = h1->GetNbinsX();
  h1->SetBinContent(nbins, h1->GetBinContent(nbins)+h1->GetBinContent(nbins+1));
  h1->SetBinContent(nbins+1,0);
}

#include <boost/algorithm/string/predicate.hpp>
std::vector<std::string> GetFileList(std::string input_file)
{
  std::vector<std::string> filenames;
  if(boost::algorithm::ends_with(input_file,".root")){
    filenames.emplace_back(input_file);
    std::cout << "Added file " << filenames.back() << " to be processed." << std::endl;
  }
  else{
    std::ifstream input(input_file);
    for(std::string line; std::getline(input,line);){
      filenames.emplace_back(line);
      std::cout << "Added file " << filenames.back() << " to be processed." << std::endl;
    }
  }
  return filenames;
}
void TruthMatchTracks(std::string input_file, art::InputTag tk_tag, art::InputTag hit_tag, art::InputTag hitparticleassns_tag) {

  gStyle->SetOptStat(0);

  //format our file list
  std::vector<std::string> filenames = GetFileList(input_file);

  art::InputTag trk_tag { "pandoraNu" };
  
  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {
    
    //to get run and event info, you use this "eventAuxillary()" object.
    cout << "Processing "
	 << "Run " << ev.eventAuxiliary().run() << ", "
	 << "Event " << ev.eventAuxiliary().event() << endl;
        
    auto const& hit_handle = ev.getValidHandle<std::vector<recob::Hit>>(hit_tag);
    auto const& hit_vec(*hit_handle);

    auto const& trk_handle = ev.getValidHandle<std::vector<recob::Track>>(trk_tag);
    auto const& trk_vec(*trk_handle);

    std::cout << "\tThere are " << trk_vec.size() << " tracks in this event." << std::endl;
    art::FindManyP<recob::Hit> hits_per_track(trk_handle, ev, trk_tag);

    for(size_t i_t=0; i_t<trk_vec.size(); ++i_t){

      std::vector< art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(i_t);
      std::cout << "\tThere are " << trk_hits_ptrs.size() << " associated hits." << std::endl;

      std::unordered_map<int,double> trkide;
      double maxe=-1, tote=0;
      simb::MCParticle const* maxp_me; //pointer for the particle match we will calculate

      art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hit_handle,ev,hitparticleassns_tag);
      std::vector<simb::MCParticle const*> particle_vec;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

      //loop only over our hits
      for(size_t i_h=0; i_h<trk_hits_ptrs.size(); ++i_h){

	particle_vec.clear(); match_vec.clear();
	particles_per_hit.get(trk_hits_ptrs[i_h].key(),particle_vec,match_vec);
	//the .key() gives us the index in the original collection

	//std::cout << "\t\tThere are " << particle_vec.size() << " particles matched to hit " << i_h << std::endl;

	//loop over particles
	for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
	  trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy; //store energy per track id
	  tote += match_vec[i_p]->energy; //calculate total energy deposited
	  if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){ //keep track of maximum
	    maxe = trkide[ particle_vec[i_p]->TrackId() ];
	    maxp_me = particle_vec[i_p];
	  }
	}//end loop over particles per hit

      }


      std::cout << "Final Match (from my loop) is " << maxp_me->TrackId() << " with energy " << maxe << " over " << tote << " (" << maxe/tote << ")"
		<< " \n\tpdg=" << maxp_me->PdgCode()
		<< " trkid=" << maxp_me->TrackId()
		<< " ke=" << maxp_me->E()-maxp_me->Mass()
		<< "\n\tstart (x,y,z)=(" << maxp_me->Vx()
		<< "," << maxp_me->Vy()
		<< "," << maxp_me->Vz()
		<< ")\tend (x,y,z)=(" << maxp_me->EndX()
		<< "," << maxp_me->EndY()
		<< "," << maxp_me->EndZ() << ")" << std::endl;
      
    }
      
  }
    
}
