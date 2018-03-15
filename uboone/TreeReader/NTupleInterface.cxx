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

  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_NParticles"              , &MCTruth_NParticles);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_TrackId"       , &MCTruth_particles_TrackId);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_PdgCode"       , &MCTruth_particles_PdgCode);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Mother"        , &MCTruth_particles_Mother);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_StatusCode"    , &MCTruth_particles_StatusCode);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_NDaughters"    , &MCTruth_particles_NDaughters);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Daughters"     , &MCTruth_particles_Daughters);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Gvx"           , &MCTruth_particles_Gvx);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Gvy"           , &MCTruth_particles_Gvy
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Gvz"           , &MCTruth_particles_Gvz);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Gvt"           , &MCTruth_particles_Gvt);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_px0"           , &MCTruth_particles_px0);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_py0"           , &MCTruth_particles_py0);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_pz0"           , &MCTruth_particles_pz0);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_e0"            , &MCTruth_particles_e0);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Rescatter"     , &MCTruth_particles_Rescatter);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_polx"          , &MCTruth_particles_polx);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_poly"          , &MCTruth_particles_poly);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_polz"          , &MCTruth_particles_polz);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_CCNC"           , &MCTruth_neutrino_CCNC);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_mode"           , &MCTruth_neutrino_mode);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_interactionType", &MCTruth_neutrino_interactionType);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_target"         , &MCTruth_neutrino_target);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_nucleon"        , &MCTruth_neutrino_nucleon);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_quark"          , &MCTruth_neutrino_quark);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_W"              , &MCTruth_neutrino_W);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_X"              , &MCTruth_neutrino_X);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_Y"              , &MCTruth_neutrino_Y);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_Q2"             , &MCTruth_neutrino_Q2);

  fNEntries = fTree->GetEntries();
  assert(fNEntries > 0);
  fTree->GetEntry(0);
  fRun = run;
}


bool NTupleInterface::FillMCFlux(Long64_t ientry, simb::MCFlux& flux) {
  if (!fTree->GetEntry(ientry)) {
    return false;
  }

  fNuPos = TLorentzVector(vtxx, vtxy, vtxz, 0);
  fNuMom = TLorentzVector(px, py, pz, E);

  flux.fptype = ptype;
  flux.fntype  = pdg;
  flux.fnimpwt = wgt;
  flux.fdk2gen = dist;
  flux.fnenergyn = nenergyn;

  return true;
}


bool NTupleInterface::FillMCTruth(Long64_t ientry, simb::MCruth& mctruth) {
  if (!fTree->GetEntry(ientry)) {
    return false;
  }

  // Add MCParticles list with neutrino first
  for (int i=0; i<MCTruth_NParticles; i++) {
    // New MCParticle
    simb::MCParticle part(
      MCTruth_particles_TrackId[i];
      MCTruth_particles_PdgCode[i];
      "",
      MCTruth_particles_Mother[i],
      simb::MCParticle::s_uninitialized,
      MCTruth_particles_StatusCode[i]);

    // Add daughter track IDs for this particle
    for (int j=0; j<MCTruth_particles_NDaughters[i]; j++) {
      part.AddDaugher(MCTruth_particles_Daughters[i][j]);
    }

    // Vertex
    part.SetGvtx(MCTruth_particles_Gvx[i],
                 MCTruth_particles_Gvy[i],
                 MCTruth_particles_Gvz[i],
                 MCTruth_particles_Gvt[i]);

    // First trajectory point, the only point used in reweighting
    part.AddTrajectoryPoint(
      TLorentzVector(),
      TLorentzVector(MCTruth_particles_px0[i],
                     MCTruth_particles_py0[i],
                     MCTruth_particles_pz0[i],
                     MCTruth_particles_e0[i])
    );

    part.SetRescatter(MCTruth_particles_Rescatter[i]);

    part.SetPolarization(TVector3(MCTruth_particles_polx[i],
                                  MCTruth_particles_poly[i],
                                  MCTruth_particles_polz[i]));

    mctruth.Add(part);
  }

  // Neutrino (automatically sets lepton from the MCParticles list)
  mctruth.SetNeutrino(MCTruth_neutrino_CCNC,
                      MCTruth_neutrino_mode,
                      MCTruth_neutrino_interactionType,
                      MCTruth_neutrino_target,
                      MCTruth_neutrino_nucleon,
                      MCTruth_neutrino_quark,
                      MCTruth_neutrino_W,
                      MCTruth_neutrino_X,
                      MCTruth_neutrino_Y,
                      MCTruth_neutrino_Q2);

  return true;
}  


bool NTupleInterface::FillMCTruth(Long64_t ientry, simb::MCruth& mctruth) {
  if (!fTree->GetEntry(ientry)) {
    return false;
  }

  //gtruth.fProbePDG
  //gtruth.fIsSeaQuark
  //gtruth.ftgtPDG
  //gtruth.gweight
  //gtruth.fprobability
  //gtruth.fXsec
  //gtruth.fDissXsec
  //gtruth.fVertex (4vec)
  //gtruth.fGscatter
  //gtruth.fGint
  //gtruth.fResNum
  //gtruth.fNumPiPlus
  //gtruth.fNumPi0
  //gtruth.fNumPiMinus
  //gtruth.fNumProton
  //gtruth.fNumNeutron
  //gtruth.fIsCharm
  //gtruth.fgX
  //gtruth.fgY
  //gtruth.fgZ
  //gtruth.fgT
  //gtruth.fgW
  //gtruth.fgQ2
  //gtruth.fgq2
  //gtruth.fShadSystP4.Px
  //gtruth.fShadSystP4.Py
  //gtruth.fShadSystP4.Pz
  //gtruth.fShadSystP4.E

  return true;
}

}  // namespace uboone

