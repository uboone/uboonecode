#ifndef __uboone_NTupleInterface__
#define __uboone_NTupleInterface__

/**
 * Base interface for importing flux files.
 *
 * \author Z. Pavlovic
 */

#include "TreeInterface.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "FluxDrivers/GNuMIFlux.h"
#include "FluxDrivers/GSimpleNtpFlux.h"

class TTree;
class TFile;

namespace uboone {

class NTupleInterface : public TreeInterface {
  public:
    NTupleInterface();
    ~NTupleInterface();

    const Long64_t GetEntries()                    { return fNEntries; }
    const int      GetRun()                        { return fRun; }
    const float    GetPOT()                        { return fPOT; }
    const TLorentzVector GetNuPosition()           { return fNuPos; }
    const TLorentzVector GetNuMomentum()           { return fNuMom; }

    void SetRootFile(TFile* rootFileName);

    void SetBranches();

    bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux);

  private:
    TTree*                        fFluxTree;
    TTree*                        fMetaTree;
    double                        fvtxx, fvtxy, fvtxz;
    double                        fpx, fpy, fpz, fE;
    int                           fpdg;
    double                        fwgt;
    double                        fdist;
    int                           fptype;
    int                           frun;
    double                        fnenergyn;
    Long64_t                      fNEntries;
    int                           fRun;
    float                         fPOT;
    TLorentzVector                fNuPos;
    TLorentzVector                fNuMom;
};

}  // namespace uboone

#endif // __uboone_NTupleInterface__

