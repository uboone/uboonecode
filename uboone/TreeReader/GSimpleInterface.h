#ifndef __uboone_GSimpleInterface__
#define __uboone_GSimpleInterface__

/**
 * Interface for importing GSimple files.
 *
 * \author Z. Pavlovic
 */

#include "TreeInterface.h"
#include "FluxDrivers/GNuMIFlux.h"
#include "FluxDrivers/GSimpleNtpFlux.h"

class TTree;
class TFile;

namespace uboone {

/**
 * \class GSimpleInterface
 * \brief Wrapper to access GSimple flux trees
 */
class GSimpleInterface : public TreeInterface {
public:
  GSimpleInterface();
  ~GSimpleInterface();
  
  const Long64_t       GetEntries()              { return fNEntries; }
  const int            GetRun()                  { return fRun; }
  const float          GetPOT()                  { return fPOT; }
  const TLorentzVector GetNuPosition()           { return fNuPos; }
  const TLorentzVector GetNuMomentum()           { return fNuMom; }

  void SetRootFile(TFile* rootFileName);
  bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux);

private:
  TTree*                        fFluxTree;
  TTree*                        fMetaTree;
  genie::flux::GSimpleNtpEntry* fGSimpleEntry;
  genie::flux::GSimpleNtpNuMI*  fGSimpleNuMI;
  genie::flux::GSimpleNtpAux*   fGSimpleAux;
  genie::flux::GSimpleNtpMeta*  fGSimpleMeta;
  Long64_t                      fNEntries;
  int                           fRun;
  float                         fPOT;
  TLorentzVector                fNuPos;
  TLorentzVector                fNuMom;
};

}  // namespace uboone

#endif  // __uboone_GSimpleInterface__

