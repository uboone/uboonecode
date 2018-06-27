/**
 * \file CompressionAlgoBase.h
 *
 * @author David Caratelli
 */

/** \addtogroup ClusterRecoUtil

    @{*/
#ifndef COMPRESSIONALGOBASE_H
#define COMPRESSIONALGOBASE_H

#include <iostream>
#include <vector>
#include <utility> // for pair

// Art Framework
#include "fhiclcpp/ParameterSet.h"

namespace compress {

  typedef std::vector<short>::const_iterator tick;

  /**
     \class CMAlgoBase
     ...
  */
  class CompressionAlgoBase {
    
  public:
    
    /// Default constructor
    CompressionAlgoBase();

    /// Constructor with input
    CompressionAlgoBase(fhicl::ParameterSet const&pset);
    
    /// Default destructor
    virtual ~CompressionAlgoBase(){}

    /// Function to reset the algorithm instance 
    virtual void Reset() { _InWF.clear(); _timeRange.clear(); }

    /// Function where compression is performed
    virtual void ApplyCompression(const std::vector<short> &waveform, const int mode, const unsigned int ch);

    /// Get range of output waveforms
    virtual const std::vector<std::pair< compress::tick, compress::tick> > GetOutputRanges() { return _timeRange; }
    
    /// Setter function for verbosity
    virtual void SetVerbose(bool doit=true) { _verbose = doit; }

    /// Setter function for debug mode
    virtual void SetDebug(bool doit=true) { _debug = doit; }

    /// Get Baselines vector
    virtual const std::vector<double> GetBaselines() { return _baselines; }

    /// Get Variances vector
    virtual const std::vector<double> GetVariances() { return _variances; }

    /// Get iterator to beginning of input WF
    virtual const tick GetInputBegin() { return _begin; }
    /// Get iterator to end of input WF
    virtual const tick GetInputEnd() { return _end; }


  protected:

    /// Boolean to choose verbose mode.
    bool _verbose;

    /// Boolean to choose debug mode.
    bool _debug;

    /// Waveform read in by algorithm...to be scanned and compressed
    std::vector<short> _InWF;

    /// vector of pairs of iterators marking the beginning/end of the segment to be saved
    std::vector<std::pair< compress::tick, compress::tick> > _timeRange;
    /// iterator pointing to current tick
    tick _thisTick;
    /// iterator pointing to beginning of input wf
    tick _begin;
    /// iterator pointing to end of input wf
    tick _end;

    /// Vector where to hold the various baselines for the various blocks
    std::vector<double> _baselines;
    /// Vector where to hold the variance measured per block
    std::vector<double> _variances;

  };

}

#endif
/** @} */ // end of doxygen group 

