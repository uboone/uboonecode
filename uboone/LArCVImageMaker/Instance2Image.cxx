#ifndef __SUPERA_INSTANCE_LAR2IMAGE_CXX__
#define __SUPERA_INSTANCE_LAR2IMAGE_CXX__

#include "Instance2Image.h"
#include "Base/larcv_logger.h"
#include "LArUtil/Geometry.h"


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
		       std::vector<larcv::Image2D>& ancestor_out_v ) {
    
    LARCV_SINFO() << "Instance ID Image ..." << std::endl;

    // we pack truth info about the ancestor particle type
    // we label energy deposition by ancestor track id
    // this groups secondaries with their primary parent
    
    // create images we are going to fill
    std::vector<larcv::Image2D> img_v;      // Filling track ID per pixel
    std::vector<larcv::Image2D> ancestor_v; // Filling ancestor ID per pixel
    std::vector<larcv::Image2D> energy_v;   // largest energy deposition
    for ( auto const& meta : meta_v ) {
      LARCV_SINFO() << meta.dump();      
      larcv::Image2D img1(meta);
      img1.paint(-1.0);
      img_v.emplace_back( std::move(img1) );
      
      larcv::Image2D img2(meta);
      img2.paint(-1.0);
      ancestor_v.emplace_back( std::move(img2) );

      larcv::Image2D Eimg(meta);
      Eimg.paint(-1.0);
      energy_v.emplace_back( std::move(Eimg) );
    }
    
    // loop over sim channel information
    //int num_no_ancestor = 0;
    //std::set<int> noancestorids;
    for (auto const& sch : sch_v) {
      auto ch = sch.Channel();
      auto const& wid   = ::supera::ChannelToWireID(ch);
      auto const& plane = wid.Plane;

      
      auto& imgmap1     = img_v.at(plane);
      auto& ancestormap = ancestor_v.at(plane);
      auto& Eimg        = energy_v.at(plane);
      
      auto const& meta = imgmap1.meta();

      // is the channel inside the meta window
      size_t col = wid.Wire;
      if (col < meta.min_x()) continue;
      if (meta.max_x() <= col) continue;
      if (plane != meta.plane()) continue;
      // remove offset to get position inside image
      col -= (size_t)(meta.min_x());
      
      // loop over energy deposition
      for (auto const tick_ides : sch.TDCIDEMap()) {
	int tick = supera::TPCTDC2Tick((double)(tick_ides.first)) + time_offset; // true deposition tick
	if (tick <= meta.min_y()) continue;
	if (tick >= meta.max_y()) continue;
	// Where is this tick in the image
	int row   = (int)meta.row(tick); // compressed position
	if ( row<0 || row>=(int)meta.rows() ) continue;
	
	// now we loop over the energy depositions in this tick
	double energy = (double)Eimg.pixel(row,col); // use to keep track the most energetic energy deposition at this point
	int ancestorid   = -1;
	
	// energy deposition at tick
	for (auto const& edep : tick_ides.second) {
	  
	  // for the Y-plane (ipass==0), we store the deposition with the highest energy
	  if (edep.energy < energy ) continue; // lower deposition than before
	  // if ( trackid2ancestorid.find( edep.trackID )==trackid2ancestorid.end() ) {
	  //   num_no_ancestor++;
	  //   noancestorids.insert( edep.trackID );
	  //   continue;
	  // }
	  
	  energy = edep.energy;
	  auto it_map = trackid2ancestorid.find( edep.trackID );

	  if ( it_map!=trackid2ancestorid.end() ) {
	    // found ID in ID -> ancestor map	  
	    ancestorid = it_map->second;
	    ancestormap.set_pixel( row, col, ancestorid );
	  }
	  
	  // regardless if we find ancestor, we store trackID
	  imgmap1.set_pixel( row, col, edep.trackID ); // raw ID
	  // we have non-zero energy deposition, valid edep
	  Eimg.set_pixel( row, col, energy );
	  
	}
      }
      
    }//end of pass loop
    // std::cout << "Number of no ancestor pixels: " << num_no_ancestor << std::endl;
    // std::cout << "Number of no ancestor IDs: " << noancestorids.size() << std::endl;
    // std::cout << "   ";
    // for ( auto const& id : noancestorids ) {
    //   std::cout << id << ", ";
    // }
    // std::cout << std::endl;
    
    // make output, compressed images
    //std::vector<larcv::Image2D> img_out_v;
    //std::vector<larcv::Image2D> ancestor_out_v;
    img_out_v.clear();
    ancestor_out_v.clear();
    for ( auto const& img : img_v ) {
      const larcv::ImageMeta& meta = img.meta();
      larcv::ImageMeta meta_out(meta.width(), meta.height(), 
				int( meta.rows()/row_compression_factor.at(meta.plane()) ), int( meta.cols()/col_compression_factor.at(meta.plane()) ),
				meta.min_x(), meta.max_y(), 
				meta.plane() );
      larcv::Image2D img_out( meta_out );
      img_out.paint(-1.0);
      img_out_v.emplace_back( std::move(img_out) );
      
      larcv::Image2D ancestor_out( meta_out );
      ancestor_out.paint(-1.0);
      ancestor_out_v.emplace_back( std::move(ancestor_out) );
    }
    
    int compressed_pixels_filled = 0;
    for (size_t iidx=0; iidx<img_out_v.size(); iidx++) {
      const larcv::Image2D& img       = img_v[iidx];
      const larcv::Image2D& ancestor  = ancestor_v[iidx];
      larcv::Image2D& imgout          = img_out_v[iidx];
      larcv::Image2D& ancestorout     = ancestor_out_v[iidx];
      size_t plane = img.meta().plane();
      const larcv::Image2D& energyimg = energy_v.at( plane );
      
      for (int rout=0; rout<(int)imgout.meta().rows(); rout++) {
	for (int clout=0; clout<(int)imgout.meta().cols(); clout++) {
	  // find max pixel to transfer
	  int rmax = 0;
	  int cmax = 0;
	  float enmax = 0;
	  for (int dr=0; dr<(int)row_compression_factor.at(plane); dr++) {
	    for (int dc=0; dc<(int)col_compression_factor.at(plane); dc++) {
	      
	      int r = rout *int(row_compression_factor.at(plane)) + dr;
	      int c = clout*int(col_compression_factor.at(plane)) + dc;
	      
	      // no label, we continue
	      if ( img.pixel(r,c)<0 )
		continue;
	      
	      float pixenergy = energyimg.pixel( r, c );
	      if ( pixenergy>enmax ) {
		enmax = pixenergy;
		rmax = r;
		cmax = c;
	      }
	      
	    }
	  }
	  
	  // set the output
	  if (enmax>0 ) {
	    imgout.set_pixel( rout, clout, img.pixel( rmax, cmax ) );
	    ancestorout.set_pixel( rout, clout, ancestor.pixel( rmax, cmax ) );
	    compressed_pixels_filled++;
	  }
	}
      }
    }//end of loop over index
    
    
    //std::cout << "compressed pixels filled: " << compressed_pixels_filled << std::endl;
    
  }
 


}
#endif
