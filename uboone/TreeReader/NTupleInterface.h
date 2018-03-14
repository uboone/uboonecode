#ifndef __uboone_NTupleInterface__
#define __uboone_NTupleInterface__

/**
 * Base interface for importing flux files.
 *
 * \author Z. Pavlovic
 */

#include "TreeInterface.h"

#include "fhiclcpp/ParameterSet.h"
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

    void SetRootFile(TFile* rootFileName, TString treeName, fhicl::ParameterSet& branchDef);

    bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux);

  private:
    TTree*                        fTree;
    TTree*                        MetaTree;
    double                        vtxx, vtxy, vtxz;
    double                        px, py, pz, E;
    int                           pdg;
    double                        wgt;
    double                        dist;
    int                           ptype;
    int                           run;
    double                        nenergyn;
    Long64_t                      fNEntries;
    int                           fRun;
    float                         fPOT;
    TLorentzVector                fNuPos;
    TLorentzVector                fNuMom;
};

}  // namespace uboone

#endif // __uboone_NTupleInterface__

