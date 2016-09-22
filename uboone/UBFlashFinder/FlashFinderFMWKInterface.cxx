#ifndef __FLASHFINDERFMWKINTERFACE_CXX__
#define __FLASHFINDERFMWKINTERFACE_CXX__

#include "FlashFinderFMWKInterface.h"
namespace pmtana {

  size_t NOpDets() {
    ::art::ServiceHandle<geo::Geometry> geo;
    return geo->NOpDets();
  }

  size_t OpDetFromOpChannel(size_t opch) {
    ::art::ServiceHandle<geo::Geometry> geo;
    return geo->OpDetFromOpChannel(opch);
  }

}
#endif
