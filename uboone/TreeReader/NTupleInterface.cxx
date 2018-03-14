#include <cassert>
#include "NTupleInterface.h"
#include "TFile.h"
#include "TTree.h"

namespace uboone {

NTupleInterface::NTupleInterface() {}

NTupleInterface::~NTupleInterface() {}

void NTupleInterface::SetRootFile(TFile* inputFile) {
  fFluxTree=dynamic_cast<TTree*>(inputFile->Get("testmod/flux"));
  
  fFluxTree->SetBranchAddress("vtxx" ,&vtxx);
  fFluxTree->SetBranchAddress("vtxy" ,&vtxy);
  fFluxTree->SetBranchAddress("vtxz" ,&vtxz);
  fFluxTree->SetBranchAddress("px"   ,&px);
  fFluxTree->SetBranchAddress("py"   ,&py);
  fFluxTree->SetBranchAddress("pz"   ,&pz);
  fFluxTree->SetBranchAddress("E"    ,&E);
  fFluxTree->SetBranchAddress("pdg"  ,&pdg);
  fFluxTree->SetBranchAddress("ptype",&ptype);
  fFluxTree->SetBranchAddress("wgt"  ,&wgt);
  fFluxTree->SetBranchAddress("dist" ,&dist);
  fFluxTree->SetBranchAddress("evtno"  ,&run);
  fFluxTree->SetBranchAddress("nenergyn", &nenergyn);

  fNEntries = fFluxTree->GetEntries();
  assert(fNEntries > 0);
  fFluxTree->GetEntry(0);
  fRun = run;
}

bool NTupleInterface::FillMCFlux(Long64_t ientry, simb::MCFlux& flux) {
  if (!fFluxTree->GetEntry(ientry))
    return false;

  fNuPos=TLorentzVector(vtxx, vtxy, vtxz, 0);
  fNuMom=TLorentzVector(px  , py  , pz  , E);

  flux.fptype = ptype;
  flux.fntype  = pdg;
  flux.fnimpwt = wgt;
  flux.fdk2gen = dist;
  flux.fnenergyn = nenergyn;

  std::cout << ">> vtxx: " << vtxx << std::endl;
  std::cout << ">> vtxy: " << vtxy << std::endl;
  std::cout << ">> vtxz: " << vtxz << std::endl;
  std::cout << ">> px  : " << px << std::endl;
  std::cout << ">> py  : " << py << std::endl;
  std::cout << ">> pz  : " << pz << std::endl;
  std::cout << ">> E   : " << E << std::endl;
  std::cout << ">> ptype:" << ptype << std::endl;
  std::cout << ">> ntype : " << pdg << std::endl;
  std::cout << ">> wgt : " << wgt << std::endl;
  std::cout << ">> dist: " << dist << std::endl;
  std::cout << ">> nenergyn: " << nenergyn << std::endl;

  return true;
}

}  // namespace uboone

