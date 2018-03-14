#ifndef __uboone_TreeInterface__
#define __uboone_TreeInterface__

/**
 * Base interface for importing flux files.
 *
 * \author Z. Pavlovic
 */

#include "nusimdata/SimulationBase/MCFlux.h"
#include "TLorentzVector.h"

namespace uboone {

/**
 * \class TreeInterface
 * \brief Wrapper to access flux trees
 */
class TreeInterface {
public:    
  virtual bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux) = 0;
  virtual const float GetPOT() = 0;
  virtual const Long64_t GetEntries() = 0;
  virtual const int GetRun() = 0;

  virtual const TLorentzVector GetNuPosition() = 0;
  virtual const TLorentzVector GetNuMomentum() = 0;
};

}  // namespace uboone

#endif  // __uboone_TreeInterface__

