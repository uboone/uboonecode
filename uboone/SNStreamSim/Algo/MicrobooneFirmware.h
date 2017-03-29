/**
 * \file compressionalgosncompress.h
 *
 * \ingroup SNCompression
 * 
 * @author David Caratelli [dcaratelli@houston.nevis.columbia.edu]
 * @author Connor Callahan
 * @author Anya Fadeeva
 */

/** \addtogroup SNCompression

    @{*/
#ifndef MICROBOONEFIRMWARE_H
#define MICROBOONEFIRMWARE_H

#include "uboone/SNStreamSim/Fmwk/CompressionAlgoBase.h"
#include <math.h>

namespace compress {

  /**
     \class compressionalgosncompress
     A Class where to write a compressiona algorithm for TPC wire signals.
  */
  class MicrobooneFirmware : public CompressionAlgoBase {
    
  public:

    MicrobooneFirmware(fhicl::ParameterSet const& pset);
    
    /// Function where compression is performed
    void ApplyCompression(const std::vector<short> &waveform, const int mode, const unsigned int ch);



    // Decide if to fill tree with info or not
    void SetFillTree(bool on) { _fillTree = on; }

  protected:

    /// Function that determines if we passed the threshold. Per plane
    bool PassThreshold(double thisADC, double base);

    // setter function for algo specifications
    void SetCompressThresh(int tU, int tV, int tY) { _thresh[0] = tU; _thresh[1] = tV; _thresh[2] = tY; }
    void SetPolarity(int polu, int polv, int poly) {_pol[0] = polu; _pol[1] = polv; _pol[2] = poly; }
    void SetUVYplaneBuffer(int upre, int upost, int vpre, int vpost, int ypre, int ypost);

    int _block;

    // Value of Baseline Threshold for initial stage of compression
    double _deltaB; // max difference for Baseline tolerared
    double _deltaV; // max difference for Variance tolerated
    std::vector<int> _thresh;
    std::vector<int> _pol;
    


    // Buffer Tick number for various planes
    // structure: [ [Upre, Upost], [Vpre, Vpost], [Ypre, Ypost] ]
    std::vector<std::vector<int> > _buffer;

    // what is the dynamic range of the ADC? This is the MAX ADC tick value.
    int _maxADC;

    // Keep track of the per-channel baseline in a map
    // The baseline will be re-set per run.
    // If the baseline is not found in the first 3 blocks
    // (because noisy) then no output can be saved.
    // this will go on until a quiet region in the run is found
    std::map<int, int> _baselineMap;

    // boolean to decide if to fill tree
    bool _fillTree;
    int _pl;
    int _v1, _v2, _v3;
    int _b1, _b2, _b3;
    double _max;
    int _interesting;
    int _save;
  };

}

#endif
/** @} */ // end of doxygen group 

