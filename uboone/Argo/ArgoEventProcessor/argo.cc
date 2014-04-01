#include "art/Framework/Art/artapp.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>

int main( int argc, char* argv[] ) {
  int result = artapp(argc,argv);
  // Make sure the message logger is shut down so our message really is
  // the last one put out.
  mf::MessageFacilityService::instance().MFPresence.reset();
  // End message.
  std::cout
    << "Art has completed and will exit with status "
    << result
    << "."
    << std::endl;
  return result;
}
