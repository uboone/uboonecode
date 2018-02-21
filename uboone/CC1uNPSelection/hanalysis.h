//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Feb 18 21:44:14 2018 by ROOT version 6.06/08
// from TTree fMC_TrunMean/Data Holder
// found on file: /pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/MCC8.6_dev/out/cc1unp_output.root
//////////////////////////////////////////////////////////

#ifndef hanalysis_h
#define hanalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>



class hanalysis : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> fRun = {fReader, "fRun"};
   TTreeReaderValue<Int_t> fSubRun = {fReader, "fSubRun"};
   TTreeReaderValue<Int_t> fEvent = {fReader, "fEvent"};
   TTreeReaderValue<Int_t> _fTrueccnc = {fReader, "_fTrueccnc"};
   TTreeReaderValue<Int_t> _fTruemode = {fReader, "_fTruemode"};
   TTreeReaderValue<Int_t> _fTrueinttype = {fReader, "_fTrueinttype"};
   TTreeReaderValue<Int_t> _fTruenupdg = {fReader, "_fTruenupdg"};
   TTreeReaderValue<Float_t> _fTrueenu = {fReader, "_fTrueenu"};
   TTreeReaderValue<Float_t> _fTrueq2truth = {fReader, "_fTrueq2truth"};
   TTreeReaderValue<Float_t> _fTrueWtruth = {fReader, "_fTrueWtruth"};
   TTreeReaderValue<Float_t> _fTrueXtruth = {fReader, "_fTrueXtruth"};
   TTreeReaderValue<Float_t> _fTrueYtruth = {fReader, "_fTrueYtruth"};
   TTreeReaderValue<Double_t> trueMuonTrueMomentum = {fReader, "trueMuonTrueMomentum"};
   TTreeReaderValue<Double_t> trueMuonTrueTheta = {fReader, "trueMuonTrueTheta"};
   TTreeReaderValue<Double_t> trueMuonTruePhi = {fReader, "trueMuonTruePhi"};
   TTreeReaderArray<double> trueProtonsTrueMomentum = {fReader, "trueProtonsTrueMomentum"};
   TTreeReaderArray<double> trueProtonsTrueTheta = {fReader, "trueProtonsTrueTheta"};
   TTreeReaderArray<double> trueProtonsTruePhi = {fReader, "trueProtonsTruePhi"};
   TTreeReaderValue<Float_t> fvtxx = {fReader, "fvtxx"};
   TTreeReaderValue<Float_t> fvtxy = {fReader, "fvtxy"};
   TTreeReaderValue<Float_t> fvtxz = {fReader, "fvtxz"};
   TTreeReaderValue<Float_t> fLlep = {fReader, "fLlep"};
   TTreeReaderValue<Float_t> fLhad = {fReader, "fLhad"};
   TTreeReaderValue<Float_t> fPlep = {fReader, "fPlep"};
   TTreeReaderValue<Float_t> fPhad = {fReader, "fPhad"};
   TTreeReaderValue<Float_t> fThetaLep = {fReader, "fThetaLep"};
   TTreeReaderValue<Float_t> fThetaHad = {fReader, "fThetaHad"};
   TTreeReaderValue<Float_t> fPhiLep = {fReader, "fPhiLep"};
   TTreeReaderValue<Float_t> fPhiHad = {fReader, "fPhiHad"};
   TTreeReaderValue<Float_t> fCosThetaLep = {fReader, "fCosThetaLep"};
   TTreeReaderValue<Float_t> fCosThetaHad = {fReader, "fCosThetaHad"};
   TTreeReaderValue<Float_t> trackendxcandidate = {fReader, "trackendxcandidate"};
   TTreeReaderValue<Float_t> trackendycandidate = {fReader, "trackendycandidate"};
   TTreeReaderValue<Float_t> trackendzcandidate = {fReader, "trackendzcandidate"};
   TTreeReaderValue<Float_t> trackstartxcandidate = {fReader, "trackstartxcandidate"};
   TTreeReaderValue<Float_t> trackstartycandidate = {fReader, "trackstartycandidate"};
   TTreeReaderValue<Float_t> trackstartzcandidate = {fReader, "trackstartzcandidate"};
   TTreeReaderValue<Int_t> fNRecoTrks = {fReader, "fNRecoTrks"};
   TTreeReaderValue<Float_t> TrunMean_cand = {fReader, "TrunMean_cand"};
   TTreeReaderValue<Float_t> TrunMean_pcand = {fReader, "TrunMean_pcand"};
   TTreeReaderValue<Int_t> fNRecoPTrks = {fReader, "fNRecoPTrks"};
   TTreeReaderValue<Int_t> fNTruePTrks = {fReader, "fNTruePTrks"};
   TTreeReaderValue<Float_t> ftrklenmuoncand = {fReader, "ftrklenmuoncand"};
   TTreeReaderValue<Float_t> ftrklenprotoncand = {fReader, "ftrklenprotoncand"};
   TTreeReaderValue<Float_t> trackmomcandidate = {fReader, "trackmomcandidate"};
   TTreeReaderValue<Float_t> trackmomcandidate_mcs = {fReader, "trackmomcandidate_mcs"};
   TTreeReaderValue<Float_t> trackmomprotoncandidate = {fReader, "trackmomprotoncandidate"};
   TTreeReaderValue<Float_t> fopflashtime = {fReader, "fopflashtime"};
   TTreeReaderValue<Float_t> fopflashmax = {fReader, "fopflashmax"};
   TTreeReaderValue<Float_t> flstrkdist = {fReader, "flstrkdist"};
   TTreeReaderValue<Int_t> Nhits_muoncand = {fReader, "Nhits_muoncand"};
   TTreeReaderValue<Int_t> Nhits_protoncand = {fReader, "Nhits_protoncand"};
   TTreeReaderValue<Int_t> truthtop = {fReader, "truthtop"};
   TTreeReaderValue<Int_t> truthtop_200thresh = {fReader, "truthtop_200thresh"};
   TTreeReaderValue<Int_t> truthtop_300thresh = {fReader, "truthtop_300thresh"};
   TTreeReaderValue<Int_t> truthtop_400thresh = {fReader, "truthtop_400thresh"};
   TTreeReaderValue<Int_t> trackcand_origin = {fReader, "trackcand_origin"};
   TTreeReaderValue<Int_t> trackcand_nuset = {fReader, "trackcand_nuset"};
   TTreeReaderValue<Int_t> trackcand_parPDG = {fReader, "trackcand_parPDG"};
   TTreeReaderValue<Int_t> trackcand_parStatusCode = {fReader, "trackcand_parStatusCode"};
   TTreeReaderValue<Float_t> trackcand_parTheta = {fReader, "trackcand_parTheta"};
   TTreeReaderValue<Float_t> trackcand_parCosTheta = {fReader, "trackcand_parCosTheta"};
   TTreeReaderValue<Float_t> trackcand_parSinTheta = {fReader, "trackcand_parSinTheta"};
   TTreeReaderValue<Float_t> trackcand_parE = {fReader, "trackcand_parE"};
   TTreeReaderValue<Float_t> trackcand_parMass = {fReader, "trackcand_parMass"};
   TTreeReaderValue<Float_t> trackcand_parKE = {fReader, "trackcand_parKE"};
   TTreeReaderValue<Float_t> trackcand_parEndE = {fReader, "trackcand_parEndE"};
   TTreeReaderValue<Float_t> trackcand_parPx = {fReader, "trackcand_parPx"};
   TTreeReaderValue<Float_t> trackcand_parPy = {fReader, "trackcand_parPy"};
   TTreeReaderValue<Float_t> trackcand_parPz = {fReader, "trackcand_parPz"};
   TTreeReaderValue<Float_t> trackcand_parPhi = {fReader, "trackcand_parPhi"};
   TTreeReaderValue<Float_t> trackcand_parCosPhi = {fReader, "trackcand_parCosPhi"};
   TTreeReaderValue<Float_t> trackcand_parSinPhi = {fReader, "trackcand_parSinPhi"};
   TTreeReaderValue<Int_t> trackpcand_origin = {fReader, "trackpcand_origin"};
   TTreeReaderValue<Int_t> trackpcand_nuset = {fReader, "trackpcand_nuset"};
   TTreeReaderValue<Int_t> trackpcand_parPDG = {fReader, "trackpcand_parPDG"};
   TTreeReaderValue<Int_t> trackpcand_parStatusCode = {fReader, "trackpcand_parStatusCode"};
   TTreeReaderValue<Float_t> trackpcand_parTheta = {fReader, "trackpcand_parTheta"};
   TTreeReaderValue<Float_t> trackpcand_parCosTheta = {fReader, "trackpcand_parCosTheta"};
   TTreeReaderValue<Float_t> trackpcand_parSinTheta = {fReader, "trackpcand_parSinTheta"};
   TTreeReaderValue<Float_t> trackpcand_parE = {fReader, "trackpcand_parE"};
   TTreeReaderValue<Float_t> trackpcand_parMass = {fReader, "trackpcand_parMass"};
   TTreeReaderValue<Float_t> trackpcand_parKE = {fReader, "trackpcand_parKE"};
   TTreeReaderValue<Float_t> trackpcand_parEndE = {fReader, "trackpcand_parEndE"};
   TTreeReaderValue<Float_t> trackpcand_parPx = {fReader, "trackpcand_parPx"};
   TTreeReaderValue<Float_t> trackpcand_parPy = {fReader, "trackpcand_parPy"};
   TTreeReaderValue<Float_t> trackpcand_parPz = {fReader, "trackpcand_parPz"};
   TTreeReaderValue<Float_t> trackpcand_parPhi = {fReader, "trackpcand_parPhi"};
   TTreeReaderValue<Float_t> trackpcand_parCosPhi = {fReader, "trackpcand_parCosPhi"};
   TTreeReaderValue<Float_t> trackpcand_parSinPhi = {fReader, "trackpcand_parSinPhi"};
   TTreeReaderValue<Float_t> Evis = {fReader, "Evis"};
   TTreeReaderValue<Float_t> Q2cal = {fReader, "Q2cal"};
   TTreeReaderValue<Float_t> Wcal = {fReader, "Wcal"};
   TTreeReaderArray<int> protoncandidate_id = {fReader, "protoncandidate_id"};
   TTreeReaderArray<float> protoncandidate_startx = {fReader, "protoncandidate_startx"};
   TTreeReaderArray<float> protoncandidate_starty = {fReader, "protoncandidate_starty"};
   TTreeReaderArray<float> protoncandidate_startz = {fReader, "protoncandidate_startz"};
   TTreeReaderArray<float> protoncandidate_endx = {fReader, "protoncandidate_endx"};
   TTreeReaderArray<float> protoncandidate_endy = {fReader, "protoncandidate_endy"};
   TTreeReaderArray<float> protoncandidate_endz = {fReader, "protoncandidate_endz"};
   TTreeReaderArray<double> protoncandidate_momentum = {fReader, "protoncandidate_momentum"};
   TTreeReaderArray<double> protoncandidate_length = {fReader, "protoncandidate_length"};
   TTreeReaderArray<double> protoncandidate_theta = {fReader, "protoncandidate_theta"};
   TTreeReaderArray<double> protoncandidate_phi = {fReader, "protoncandidate_phi"};
   TTreeReaderArray<double> protoncandidate_trunmeandqdx = {fReader, "protoncandidate_trunmeandqdx"};
   TTreeReaderArray<double> protoncandidate_pida = {fReader, "protoncandidate_pida"};


   hanalysis(TTree * /*tree*/ =0) { }
   virtual ~hanalysis() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(hanalysis,0);

};

#endif

#ifdef hanalysis_cxx
void hanalysis::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t hanalysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef hanalysis_cxx
