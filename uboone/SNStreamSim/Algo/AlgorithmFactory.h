#ifndef COMPRESS_ALGORITHMFACTORY_H
#define COMPRESS_ALGORITHMFACTORY_H

#include <string>
#include <cstdlib>

// Abstract algorithm class include
#include "uboone/SNStreamSim/Fmwk/CompressionAlgoBase.h"

// Framework Includes
#include "fhiclcpp/ParameterSet.h"

namespace compress {


  class AlgorithmFactory {

  public:
    
    AlgorithmFactory() {}
    ~AlgorithmFactory() {}

    std::unique_ptr< CompressionAlgoBase > MakeCompressionAlgo(fhicl::ParameterSet const& p);
  
  };
  
}

#endif
