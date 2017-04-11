#ifndef COMPRESS_ALGORITHMFACTORY_CXX
#define COMPRESS_ALGORITHMFACTORY_CXX

#include "AlgorithmFactory.h"

// Algorithms
#include "MicrobooneFirmware.h"

namespace compress {

  std::unique_ptr< CompressionAlgoBase > AlgorithmFactory::MakeCompressionAlgo(fhicl::ParameterSet const& p) {

    std::cout << "creating algorithm" << std::endl;

    std::string algname = p.get<std::string>("CompressionAlgoName");

    std::unique_ptr< CompressionAlgoBase > ptr;

    if (algname.compare("MicrobooneFirmware") == 0) {
      std::unique_ptr< CompressionAlgoBase> new_ptr(new MicrobooneFirmware(p) );
      ptr.swap(new_ptr);
    }
    
    else {
      std::cout << "Algorithm name provided is " << algname << std::endl;
      throw std::runtime_error("ERROR in AlgorithmFactory: no registered algorithm by that name.");
    }

    return std::move(ptr);
  }

}

#endif
