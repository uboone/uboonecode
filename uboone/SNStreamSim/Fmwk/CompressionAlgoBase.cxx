#ifndef COMPRESSIONALGOBASE_CXX
#define COMPRESSIONALGOBASE_CXX

#include "CompressionAlgoBase.h"

namespace compress {

  CompressionAlgoBase::CompressionAlgoBase(){

    _verbose = false; 
  }

  CompressionAlgoBase::CompressionAlgoBase(fhicl::ParameterSet  const &pset){
    CompressionAlgoBase();
  }
  
  void CompressionAlgoBase::ApplyCompression(const std::vector<short> &waveform, int mode, const unsigned int ch){

    if (_verbose) { std::cout << "Exploring plane " << mode << std::endl; }

    /// This algorithm simply returns the original waveform as the compressed one! Not very useful!
    std::cout << "Defaut Base Compression Algo called..." << std::endl;
    //Set the Output waveform to be identical to the input waveform
    //_OutWF.push_back(waveform);
    // Set the start tick for the output waveform to be 0
    // if this were something actually serious/useful
    // we should set the start-tick to be the time-tick
    // in the old waveform at which the new waveform starts
    //_OutWFStartTick.push_back(0);

    // make a pair that contains the entire vector
    _timeRange.push_back(std::make_pair(waveform.begin(),waveform.end()));
    _begin = waveform.begin();
    _end   = waveform.end();

    return;
  }
  
}

#endif
