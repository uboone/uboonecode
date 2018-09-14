#ifndef __SUPERA_LARFLOW_LAR2IMAGE_H__
#define __SUPERA_LARFLOW_LAR2IMAGE_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "FMWKInterface.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/EventChStatus.h"

namespace supera {

  //
  // SimChannel => PixelFlowMaps
  // 
  std::vector<larcv::Image2D>
    SimCh2LArFlowImages( const std::vector<larcv::ImageMeta>& meta_v,
			 const std::vector<larcv::ROIType_t>& track2type_v,			 
			 const std::vector<supera::LArSimCh_t>& sch_v,
			 const larcv::EventChStatus& ev_chstatus,
			 const std::vector<float>& row_compression_factor,
			 const std::vector<float>& col_compression_factor,
			 const int time_offset );
  
}
#endif
//#endif
//#endif
