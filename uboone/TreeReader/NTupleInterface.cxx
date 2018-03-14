#include <cassert>
#include "NTupleInterface.h"
#include "TFile.h"
#include "TTree.h"

namespace uboone {

NTupleInterface::NTupleInterface() {}

NTupleInterface::~NTupleInterface() {}

void NTupleInterface::SetRootFile(TFile* inputFile, TString /*treeName*/, fhicl::ParameterSet& /*branchDef*/) {
  fTree=dynamic_cast<TTree*>(inputFile->Get("testmod/flux"));
  
  fTree->SetBranchAddress("vtxx" ,&vtxx);
  fTree->SetBranchAddress("vtxy" ,&vtxy);
  fTree->SetBranchAddress("vtxz" ,&vtxz);
  fTree->SetBranchAddress("px"   ,&px);
  fTree->SetBranchAddress("py"   ,&py);
  fTree->SetBranchAddress("pz"   ,&pz);
  fTree->SetBranchAddress("E"    ,&E);
  fTree->SetBranchAddress("pdg"  ,&pdg);
  fTree->SetBranchAddress("ptype",&ptype);
  fTree->SetBranchAddress("wgt"  ,&wgt);
  fTree->SetBranchAddress("dist" ,&dist);
  fTree->SetBranchAddress("evtno"  ,&run);
  fTree->SetBranchAddress("nenergyn", &nenergyn);

  fNEntries = fTree->GetEntries();
  assert(fNEntries > 0);
  fTree->GetEntry(0);
  fRun = run;
}

bool NTupleInterface::FillMCFlux(Long64_t ientry, simb::MCFlux& flux) {
  if (!fTree->GetEntry(ientry))
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

