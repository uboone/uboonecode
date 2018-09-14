#ifndef __SUPERALARFLOW_CXX__
#define __SUPERALARFLOW_CXX__

#include "SuperaLArFlow.h"
#include "LArFlowLAr2Image.h"
#include "ImageMetaMakerFactory.h"
#include "PulledPork3DSlicer.h"
#include "DataFormat/EventImage2D.h"
#include "DataFormat/DataFormatUtil.h"
namespace larcv {

  static SuperaLArFlowProcessFactory __global_SuperaLArFlowProcessFactory__;

  SuperaLArFlow::SuperaLArFlow(const std::string name)
    : SuperaBase(name)
  {}
    
  void SuperaLArFlow::configure(const PSet& cfg)
  {
    SuperaBase::configure(cfg);
    supera::ParamsImage2D::configure(cfg);
    supera::ImageMetaMaker::configure(cfg);
    _origin = cfg.get<unsigned short>("Origin",0);
    _chstatus_producer = cfg.get<std::string>("ChStatusProducer");
  }

  void SuperaLArFlow::initialize()
  {}

  bool SuperaLArFlow::process(IOManager& mgr)
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

    std::vector<larcv::ROIType_t> track2type_v;
    for(auto const& mctrack : LArData<supera::LArMCTrack_t>()) {
      //std::cout << "mctrack id=" << mctrack.TrackID() << " pdg=" << mctrack.PdgCode() << " origin=" << mctrack.Origin() << std::endl;

      if(_origin && ((unsigned short)(mctrack.Origin())) != _origin) continue;
      
      if(mctrack.TrackID() >= track2type_v.size())
	track2type_v.resize(mctrack.TrackID()+1,larcv::ROIType_t::kROIUnknown);
      track2type_v[mctrack.TrackID()] = larcv::PdgCode2ROIType(mctrack.PdgCode());

      if(mctrack.MotherTrackID() >= track2type_v.size())
	track2type_v.resize(mctrack.MotherTrackID()+1,larcv::ROIType_t::kROIUnknown);
      track2type_v[mctrack.MotherTrackID()] = larcv::PdgCode2ROIType(mctrack.MotherPdgCode());

      if(mctrack.AncestorTrackID() >= track2type_v.size())
	track2type_v.resize(mctrack.AncestorTrackID()+1,larcv::ROIType_t::kROIUnknown);
      track2type_v[mctrack.AncestorTrackID()] = larcv::PdgCode2ROIType(mctrack.AncestorPdgCode());
    }
    for(auto const& mcshower : LArData<supera::LArMCShower_t>()) {

      if(_origin && ((unsigned short)(mcshower.Origin())) != _origin) continue;
      
      if(mcshower.TrackID() >= track2type_v.size())
	track2type_v.resize(mcshower.TrackID()+1,larcv::ROIType_t::kROIUnknown);
      track2type_v[mcshower.TrackID()] = larcv::PdgCode2ROIType(mcshower.PdgCode());

      if(mcshower.MotherTrackID() >= track2type_v.size())
	track2type_v.resize(mcshower.MotherTrackID()+1,larcv::ROIType_t::kROIUnknown);
      track2type_v[mcshower.MotherTrackID()] = larcv::PdgCode2ROIType(mcshower.MotherPdgCode());
      
      if(mcshower.AncestorTrackID() >= track2type_v.size())
	track2type_v.resize(mcshower.AncestorTrackID()+1,larcv::ROIType_t::kROIUnknown);
      track2type_v[mcshower.AncestorTrackID()] = larcv::PdgCode2ROIType(mcshower.AncestorPdgCode());

      for(auto const& daughter_track_id : mcshower.DaughterTrackID()) {
	if(daughter_track_id == larcv::kINVALID_UINT)
	  continue;
	if(daughter_track_id >= track2type_v.size())
	  track2type_v.resize(daughter_track_id+1,larcv::ROIType_t::kROIUnknown);
	track2type_v[daughter_track_id] = track2type_v[mcshower.TrackID()];
      }
    }


    std::vector<float> row_compression_factor;
    std::vector<float> col_compression_factor;
    for ( auto const& meta : meta_v ) {
      row_compression_factor.push_back( RowCompressionFactor().at(meta.plane()) );
      col_compression_factor.push_back( ColCompressionFactor().at(meta.plane()) );
    }

    // get ChStatus Data
    auto ev_chstatus = (EventChStatus*)(mgr.get_data(kProductChStatus,_chstatus_producer));

    auto image_v = supera::SimCh2LArFlowImages(meta_v, track2type_v, LArData<supera::LArSimCh_t>(), *ev_chstatus,
					       row_compression_factor, col_compression_factor, TimeOffset() );
          
    ev_image->Emplace(std::move(image_v));
    
    return true;
  }

  void SuperaLArFlow::finalize()
  {}

}
#endif
