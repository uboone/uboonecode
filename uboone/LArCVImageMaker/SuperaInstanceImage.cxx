#ifndef __SUPERALARFLOW_CXX__
#define __SUPERALARFLOW_CXX__

#include "SuperaInstanceImage.h"
#include "Instance2Image.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/DataFormatUtil.h"
namespace larcv {

  static SuperaInstanceImageProcessFactory __global_SuperaInstanceImageProcessFactory__;
  
  SuperaInstanceImage::SuperaInstanceImage(const std::string name)
    : SuperaBase(name)
  {}
  
  void SuperaInstanceImage::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsImage2D::configure(cfg);
    supera::ImageMetaMaker::configure(cfg);
    _origin = cfg.get<unsigned short>("Origin",0);
    m_ancestor_label = cfg.get<std::string>("AncestorImageLabel");
  }

  void SuperaInstanceImage::initialize()
  {}

  bool SuperaInstanceImage::process(IOManager& mgr)
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
    
    auto ev_ancestor = (EventImage2D*)(mgr.get_data(kProductImage2D,m_ancestor_label));
    if(!ev_ancestor) {
      LARCV_CRITICAL() << "Ancestor image could not be created!" << std::endl;
      throw larbys();
    }
    if(!(ev_ancestor->Image2DArray().empty())) {
      LARCV_CRITICAL() << "Ancestor image array not empty!" << std::endl;
      throw larbys();
    }


    // the map I need to make
    // [trackid] -> [ancenstor id]

    std::cout << "==============================================" << std::endl;
    std::cout << "MC Track Scraping" << std::endl;

    std::map<int,int> trackid2ancestorid;

    for(auto const& mctrack : LArData<supera::LArMCTrack_t>()) {
      // std::cout << "mctrack: "
      // 		<< " id=" << mctrack.TrackID() 
      // 		<< " ancestorid=" << mctrack.AncestorTrackID() 
      // 		<< " motherid=" << mctrack.MotherTrackID() 
      // 		<< " pdg=" << mctrack.PdgCode() 
      // 		<< " origin=" << mctrack.Origin() 
      // 		<< std::endl;

      if(_origin && ((unsigned short)(mctrack.Origin())) != _origin) continue;
      
      trackid2ancestorid[mctrack.TrackID()] = mctrack.AncestorTrackID();
    }

    std::cout << "==============================================" << std::endl;
    std::cout << "MC Shower Scraping" << std::endl;
    for(auto const& mcshower : LArData<supera::LArMCShower_t>()) {

      // std::cout << "mcshower: "
      // 		<< " id=" << mcshower.TrackID() 
      // 		<< " ancestorid=" << mcshower.AncestorTrackID() 
      // 		<< " motherid=" << mcshower.MotherTrackID() 
      // 		<< " pdg=" << mcshower.PdgCode() 
      // 		<< " origin=" << mcshower.Origin() 
      // 		<< std::endl;

      if(_origin && ((unsigned short)(mcshower.Origin())) != _origin) continue;
      
      trackid2ancestorid[mcshower.TrackID()] = mcshower.AncestorTrackID();
    }


    std::vector<float> row_compression_factor;
    std::vector<float> col_compression_factor;
    for ( auto const& meta : meta_v ) {
      row_compression_factor.push_back( RowCompressionFactor().at(meta.plane()) );
      col_compression_factor.push_back( ColCompressionFactor().at(meta.plane()) );
    }

    std::vector<larcv::Image2D> idimg_v;
    std::vector<larcv::Image2D> ancestor_v;
    supera::Instance2Image(meta_v, trackid2ancestorid, LArData<supera::LArSimCh_t>(),
			   row_compression_factor, col_compression_factor, TimeOffset(),
			   idimg_v, ancestor_v );

          
    ev_image->Emplace(std::move(idimg_v));
    ev_ancestor->Emplace(std::move(ancestor_v));
    
    return true;
  }

  void SuperaInstanceImage::finalize()
  {}

}
#endif
