#ifndef __SUPERASIMCH_CXX__
#define __SUPERASIMCH_CXX__

#include "SuperaSimCh.h"
#include "LAr2Image.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/DataFormatUtil.h"
namespace larcv {

  static SuperaSimChProcessFactory __global_SuperaSimChProcessFactory__;

  SuperaSimCh::SuperaSimCh(const std::string name)
    : SuperaBase(name)
  {}
    
  void SuperaSimCh::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsImage2D::configure(cfg);
    supera::ImageMetaMaker::configure(cfg);
    _origin = cfg.get<unsigned short>("Origin",0);
  }

  void SuperaSimCh::initialize()
  {}

  void SuperaSimCh::PdgCode2ROIType(int & this_pdg, std::vector<int> & pdg_counter){
    
    if (this_pdg < 0 ) this_pdg=this_pdg*(-1);
    
    std::cout<<"this pdg is"<<this_pdg<<std::endl;

    if (!(this_pdg == 11   ||
	  this_pdg == 13   ||
	  this_pdg == 22   ||
	  this_pdg == 111  ||//Pi0
	  this_pdg == 211  ||
	  this_pdg == 321  ||//Kaon
	  this_pdg == 2212 )) this_pdg=2213; //other set to 2213 (1st prime number after 2212)
    
    if ((unsigned)this_pdg > pdg_counter.size())
      pdg_counter.resize(this_pdg + 1, 0);
    else 
      pdg_counter[this_pdg]++;
  }
  
  bool SuperaSimCh::process(IOManager& mgr)
  {
    SuperaBase::process(mgr);

    if(supera::PulledPork3DSlicer::Is(supera::ImageMetaMaker::MetaMakerPtr())) {
      auto ptr = (supera::PulledPork3DSlicer*)(supera::ImageMetaMaker::MetaMakerPtr());
      ptr->ClearEventData();
      ptr->AddConstraint(LArData<supera::LArMCTruth_t>());
      ptr->GenerateMeta(LArData<supera::LArSimCh_t>(),TimeOffset());
    }
    
    auto const& meta_v = Meta();
    
    if(meta_v.empty()) {
      LARCV_CRITICAL() << "Meta not created!" << std::endl;
      throw larbys();
    }
    auto ev_image = (EventImage2D*)(mgr.get_data(kProductImage2D,OutImageLabel()));
    if(!ev_image) {
      LARCV_CRITICAL() << "Output image could not be created!" << std::endl;
      throw larbys();
    }
    if(!(ev_image->Image2DArray().empty())) {
      LARCV_CRITICAL() << "Output image array not empty!" << std::endl;
      throw larbys();
    }

    //std::vector<larcv::ROIType_t> track2type_v;
    std::vector<int> track2type_v;
    std::vector<int> pdg_counter;

    std::cout<<"starting with track2type_v size of "<<track2type_v.size()<<std::endl;

    //Loop over track particles
    for(auto const& mctrack : LArData<supera::LArMCTrack_t>()) {

      if(_origin && ((unsigned short)(mctrack.Origin())) != _origin) continue;
      
      std::cout<<"mctrack.TrackID()          "<<mctrack.TrackID()<<std::endl;
      std::cout<<"mctrack.PdgCode()          "<<mctrack.PdgCode()<<std::endl;
      std::cout<<"mctrack.MotherTrackID()    "<<mctrack.MotherTrackID()<<std::endl;
      std::cout<<"mctrack.MotherPdgCode()    "<<mctrack.MotherPdgCode()<<std::endl;
      std::cout<<"mctrack.AncestorTrackID()  "<<mctrack.AncestorTrackID()<<std::endl;
      std::cout<<"mctrack.AncestorPdgCode()  "<<mctrack.AncestorPdgCode()<<std::endl;

      if(mctrack.TrackID() >= track2type_v.size())
	track2type_v.resize(mctrack.TrackID()+1, 0);
	//track2type_v.resize(mctrack.TrackID()+1,larcv::ROIType_t::kROIUnknown);
      if (mctrack.PdgCode()==11 or mctrack.PdgCode()==22)
	track2type_v[mctrack.TrackID()] = larcv::PdgCode2ROIType(mctrack.PdgCode());
      else {
	int mctrack_pdgcode = mctrack.PdgCode();
	PdgCode2ROIType(mctrack_pdgcode, pdg_counter);
	track2type_v[mctrack.TrackID()]=abs(mctrack_pdgcode)*100 + pdg_counter[abs(mctrack_pdgcode)];
	std::cout<<"track2type_v[mctrack.TrackID()] is "<<track2type_v[mctrack.TrackID()]<<std::endl;
	std::cout<<"mctrack_pdgcode "<<mctrack_pdgcode<<std::endl;
      }
      
      //track2type_v[mctrack.TrackID()] = larcv::PdgCode2ROIType(mctrack_pdgcode);
      //std::cout<<"mctrack Pdgcode is "<<larcv::PdgCode2ROIType(mctrack_pdgcode)<<std::endl;

      if(mctrack.MotherTrackID() >= track2type_v.size())
	track2type_v.resize(mctrack.MotherTrackID()+1,0);
      //track2type_v.resize(mctrack.MotherTrackID()+1,larcv::ROIType_t::kROIUnknown);
      if (mctrack.MotherPdgCode()==11 or mctrack.MotherPdgCode()==22)
	track2type_v[mctrack.MotherTrackID()] = larcv::PdgCode2ROIType(mctrack.MotherPdgCode());
      else {
	int mctrack_motherpdgcode = mctrack.MotherPdgCode();
	PdgCode2ROIType(mctrack_motherpdgcode, pdg_counter);
	track2type_v[mctrack.MotherTrackID()]=abs(mctrack_motherpdgcode)*100 + pdg_counter[abs(mctrack_motherpdgcode)];
	std::cout<<"track2type_v[mctrack.MotherTrackID()] is "<<track2type_v[mctrack.MotherTrackID()]<<std::endl;
	std::cout<<"mctrack_motherpdgcode "<<mctrack_motherpdgcode<<std::endl;
      }
      //track2type_v[mctrack.MotherTrackID()] = larcv::PdgCode2ROIType(mctrack_motherpdgcode);
      //std::cout<<"mctrack MotherPdgcode is "<<larcv::PdgCode2ROIType(mctrack_motherpdgcode)<<std::endl;
      

      if(mctrack.AncestorTrackID() >= track2type_v.size())
	track2type_v.resize(mctrack.AncestorTrackID()+1,0);
      //track2type_v.resize(mctrack.AncestorTrackID()+1,larcv::ROIType_t::kROIUnknown);
      if (mctrack.AncestorPdgCode()==11 or mctrack.AncestorPdgCode()==22)
	track2type_v[mctrack.AncestorTrackID()] = larcv::PdgCode2ROIType(mctrack.AncestorPdgCode());
      else {
	int mctrack_ancestorpdgcode = mctrack.AncestorPdgCode();
	PdgCode2ROIType(mctrack_ancestorpdgcode, pdg_counter);
	track2type_v[mctrack.AncestorTrackID()]=abs(mctrack_ancestorpdgcode)*100 + pdg_counter[abs(mctrack_ancestorpdgcode)];
	std::cout<<"track2type_v[mctrack.AncestorTrackID()] is "<<track2type_v[mctrack.AncestorTrackID()]<<std::endl;
	std::cout<<"mctrack_ancestorpdgcode "<<mctrack_ancestorpdgcode<<std::endl;
      }
      //track2type_v[mctrack.AncestorTrackID()] = larcv::PdgCode2ROIType(mctrack_ancestorpdgcode);
      //std::cout<<"mctrack AncestorPdgcode is "<<larcv::PdgCode2ROIType(mctrack_ancestorpdgcode)<<std::endl;
    }
    
    //Loop over shower particles
    for(auto const& mcshower : LArData<supera::LArMCShower_t>()) {

      if(_origin && ((unsigned short)(mcshower.Origin())) != _origin) continue;
      
      std::cout<<"mcshower.TrackID()          "<<mcshower.TrackID()<<std::endl;
      std::cout<<"mcshower.PdgCode()          "<<mcshower.PdgCode()<<std::endl;
      std::cout<<"mcshower.MotherTrackID()    "<<mcshower.MotherTrackID()<<std::endl;
      std::cout<<"mcshower.MotherPdgCode()    "<<mcshower.MotherPdgCode()<<std::endl;
      std::cout<<"mcshower.AncestorTrackID()  "<<mcshower.AncestorTrackID()<<std::endl;
      std::cout<<"mcshower.AncestorPdgCode()  "<<mcshower.AncestorPdgCode()<<std::endl;

      if(mcshower.TrackID() >= track2type_v.size())
	track2type_v.resize(mcshower.TrackID()+1, 0 );
        //track2type_v.resize(mcshower.TrackID()+1,larcv::ROIType_t::kROIUnknown);
      if (mcshower.PdgCode()==11 or mcshower.PdgCode()==22)
	track2type_v[mcshower.TrackID()] = larcv::PdgCode2ROIType(mcshower.PdgCode());
      else{
	int mcshower_pdgcode = mcshower.PdgCode();
	PdgCode2ROIType(mcshower_pdgcode, pdg_counter);
	track2type_v[mcshower.TrackID()]=abs(mcshower_pdgcode)*100 + pdg_counter[abs(mcshower_pdgcode)];
      }
      
      //std::cout<<"track2type_v[mcshower.TrackID()] is "<<track2type_v[mcshower.TrackID()]<<std::endl;
      
      if(mcshower.MotherTrackID() >= track2type_v.size())
	track2type_v.resize(mcshower.MotherTrackID()+1,0);
        //track2type_v.resize(mcshower.MotherTrackID()+1,larcv::ROIType_t::kROIUnknown);
      if (mcshower.MotherPdgCode()==11 or mcshower.MotherPdgCode()==22)
	track2type_v[mcshower.MotherTrackID()] = larcv::PdgCode2ROIType(mcshower.MotherPdgCode());
      else{
	int mcshower_motherpdgcode = mcshower.MotherPdgCode();
	PdgCode2ROIType(mcshower_motherpdgcode, pdg_counter);
	track2type_v[mcshower.MotherTrackID()]=abs(mcshower_motherpdgcode)*100 + pdg_counter[abs(mcshower_motherpdgcode)];
      }	
      std::cout<<"track2type_v[mcshower.MotherTrackID()] is "<<track2type_v[mcshower.MotherTrackID()]<<std::endl;

      if(mcshower.AncestorTrackID() >= track2type_v.size())
	track2type_v.resize(mcshower.AncestorTrackID()+1,0);
        //track2type_v.resize(mcshower.AncestorTrackID()+1,larcv::ROIType_t::kROIUnknown);
      if (mcshower.AncestorPdgCode()==11 or mcshower.AncestorPdgCode()==22)
	track2type_v[mcshower.AncestorTrackID()] = larcv::PdgCode2ROIType(mcshower.AncestorPdgCode());
      else{
      	int mcshower_ancestorpdgcode = mcshower.AncestorPdgCode();
	PdgCode2ROIType(mcshower_ancestorpdgcode, pdg_counter);
	track2type_v[mcshower.AncestorTrackID()]=abs(mcshower_ancestorpdgcode)*100 + pdg_counter[abs(mcshower_ancestorpdgcode)];
      }
      std::cout<<"track2type_v[mcshower.AncestorTrackID()] is "<<track2type_v[mcshower.AncestorTrackID()]<<std::endl;

      for(auto const& daughter_track_id : mcshower.DaughterTrackID()) {
	if(daughter_track_id == larcv::kINVALID_UINT)
	  continue;
	if(daughter_track_id >= track2type_v.size())
	  track2type_v.resize(daughter_track_id+1,larcv::ROIType_t::kROIUnknown);
	//std::cout<<"in daughter, track2type_v[mcshower.TrackID()]"<<track2type_v[mcshower.TrackID()]<<std::endl;
	track2type_v[daughter_track_id] = track2type_v[mcshower.TrackID()];
      }
    }
    
    //std::cout<<"size of track2type_v is "<<track2type_v.size()<<std::endl;
    //for (auto shit:track2type_v) std::cout<<"fuck, roi"<<shit<<std::endl;

    auto image_v = supera::SimCh2Image2D(meta_v, track2type_v, LArData<supera::LArSimCh_t>(), TimeOffset());

    for(size_t plane=0; plane<image_v.size(); ++plane) {
      auto& image = image_v[plane];
      image.compress(image.meta().rows() / RowCompressionFactor().at(plane),
		     image.meta().cols() / ColCompressionFactor().at(plane),
		     larcv::Image2D::kMaxPool);
    }
    
    ev_image->Emplace(std::move(image_v));
    
    return true;
  }

  void SuperaSimCh::finalize()
  {}

}
#endif
