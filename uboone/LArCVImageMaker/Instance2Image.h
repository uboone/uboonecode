#ifndef __SUPERA_INSTANCE_LAR2IMAGE_H__
#define __SUPERA_INSTANCE_LAR2IMAGE_H__
//#ifndef __CINT__
//#ifndef __CLING__
#include "FMWKInterface.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/EventChStatus.h"
#include <map>

namespace supera {

  //
  // SimChannel => PixelFlowMaps
  // 
  void Instance2Image( const std::vector<larcv::ImageMeta>& meta_v,
		       const std::map<int,int>& trackid2ancestorid,
		       const std::vector<supera::LArSimCh_t>& sch_v,
		       const std::vector<float>& row_compression_factor,
		       const std::vector<float>& col_compression_factor,
		       const int time_offset,
		       std::vector<larcv::Image2D>& img_out_v,
		       std::vector<larcv::Image2D>& ancestor_out_v );

  
}
#endif
//#endif
//#endif
