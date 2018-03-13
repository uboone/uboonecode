#ifndef _NTUPLEINTERFACE_H_
#define _NTUPLEINTERFACE_H_

#include "FluxInterface.h"

#include "FluxDrivers/GNuMIFlux.h"
#include "FluxDrivers/GSimpleNtpFlux.h"

class TTree;
class TFile;

namespace fluxr {
  class NTupleInterface : public FluxInterface
  {
    public:
      NTupleInterface();
      ~NTupleInterface();

      const Long64_t GetEntries()                    {return fNEntries;};
      const int      GetRun()                        {return fRun;};
      const float    GetPOT()                        {return fPOT;};
      const TLorentzVector GetNuPosition()           {return fNuPos;};
      const TLorentzVector GetNuMomentum()           {return fNuMom;};

      void SetRootFile(TFile* rootFileName);
      bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux);

    private:
      TTree*                        fFluxTree;
      TTree*                        fMetaTree;
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

}

#endif // _NTUPLEINTERFACE_H_
