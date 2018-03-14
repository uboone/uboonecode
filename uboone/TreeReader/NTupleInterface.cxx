#include <cassert>
#include "NTupleInterface.h"
#include "TFile.h"
#include "TTree.h"

namespace uboone {

NTupleInterface::NTupleInterface() {}

NTupleInterface::~NTupleInterface() {}

void NTupleInterface::SetRootFile(TFile* inputFile, TString treeName, fhicl::ParameterSet& branchDef) {
  std::cout << "TREE NAME: " << treeName << std::endl;
  fTree=dynamic_cast<TTree*>(inputFile->Get(treeName));
  
  fTree->SetBranchAddress(branchDef.get<std::string>("vtxx").c_str()    , &vtxx);
  fTree->SetBranchAddress(branchDef.get<std::string>("vtxy").c_str()    , &vtxy);
  fTree->SetBranchAddress(branchDef.get<std::string>("vtxz").c_str()    , &vtxz);
  fTree->SetBranchAddress(branchDef.get<std::string>("px").c_str()      , &px);
  fTree->SetBranchAddress(branchDef.get<std::string>("py").c_str()      , &py);
  fTree->SetBranchAddress(branchDef.get<std::string>("pz").c_str()      , &pz);
  fTree->SetBranchAddress(branchDef.get<std::string>("E").c_str()       , &E);
  fTree->SetBranchAddress(branchDef.get<std::string>("pdg").c_str()     , &pdg);
  fTree->SetBranchAddress(branchDef.get<std::string>("ptype").c_str()   , &ptype);
  fTree->SetBranchAddress(branchDef.get<std::string>("wgt").c_str()     , &wgt);
  fTree->SetBranchAddress(branchDef.get<std::string>("dist").c_str()    , &dist);
  fTree->SetBranchAddress(branchDef.get<std::string>("evtno").c_str()   , &run);
  fTree->SetBranchAddress(branchDef.get<std::string>("nenergyn").c_str(), &nenergyn);

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

