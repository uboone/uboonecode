#include <cassert>
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "NTupleInterface.h"
#include "TFile.h"
#include "TTree.h"

namespace uboone {

NTupleInterface::NTupleInterface() {}


NTupleInterface::~NTupleInterface() {}


void NTupleInterface::SetRootFile(TFile* inputFile, TString treeName, fhicl::ParameterSet& branchDef) {
  std::cout << "[NTUPLEINTERFACE] TREE NAME: " << treeName << std::endl;

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
  fTree->SetBranchAddress(branchDef.get<std::string>("tpx").c_str()     , &tpx);
  fTree->SetBranchAddress(branchDef.get<std::string>("tpy").c_str()     , &tpy);
  fTree->SetBranchAddress(branchDef.get<std::string>("tpz").c_str()     , &tpz);
  fTree->SetBranchAddress(branchDef.get<std::string>("tptype").c_str()  , &tptype);
  fTree->SetBranchAddress(branchDef.get<std::string>("vx").c_str()      , &vx);
  fTree->SetBranchAddress(branchDef.get<std::string>("vy").c_str()      , &vy);
  fTree->SetBranchAddress(branchDef.get<std::string>("vz").c_str()      , &vz);

  std::cout << "[NTUPLEINTERFACE] >> Tree configured." << std::endl;

  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_NParticles").c_str()              , &MCTruth_NParticles);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_TrackId").c_str()       , &MCTruth_particles_TrackId);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_PdgCode").c_str()       , &MCTruth_particles_PdgCode);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Mother").c_str()        , &MCTruth_particles_Mother);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_StatusCode").c_str()    , &MCTruth_particles_StatusCode);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_NDaughters").c_str()    , &MCTruth_particles_NDaughters);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Daughters").c_str()     , &MCTruth_particles_Daughters);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Gvx").c_str()           , &MCTruth_particles_Gvx);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Gvy").c_str()           , &MCTruth_particles_Gvy);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Gvz").c_str()           , &MCTruth_particles_Gvz);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Gvt").c_str()           , &MCTruth_particles_Gvt);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_px0").c_str()           , &MCTruth_particles_px0);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_py0").c_str()           , &MCTruth_particles_py0);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_pz0").c_str()           , &MCTruth_particles_pz0);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_e0").c_str()            , &MCTruth_particles_e0);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Rescatter").c_str()     , &MCTruth_particles_Rescatter);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_polx").c_str()          , &MCTruth_particles_polx);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_poly").c_str()          , &MCTruth_particles_poly);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_polz").c_str()          , &MCTruth_particles_polz);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_CCNC").c_str()           , &MCTruth_neutrino_CCNC);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_mode").c_str()           , &MCTruth_neutrino_mode);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_interactionType").c_str(), &MCTruth_neutrino_interactionType);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_target").c_str()         , &MCTruth_neutrino_target);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_nucleon").c_str()        , &MCTruth_neutrino_nucleon);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_quark").c_str()          , &MCTruth_neutrino_quark);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_W").c_str()              , &MCTruth_neutrino_W);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_X").c_str()              , &MCTruth_neutrino_X);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_Y").c_str()              , &MCTruth_neutrino_Y);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_Q2").c_str()             , &MCTruth_neutrino_Q2);

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

  flux.fptype     = ptype;
  flux.fntype     = pdg;
  flux.fnimpwt    = wgt;
  flux.fdk2gen    = dist;
  flux.fnenergyn  = nenergyn;
  flux.ftpx       = tpx;
  flux.ftpy       = tpy;
  flux.ftpz       = tpz;
  flux.fvx        = vx;
  flux.fvy        = vy;
  flux.fvz        = vz;
  flux.ftptype    = tptype;

  return true;
}


bool NTupleInterface::FillMCTruth(Long64_t ientry, simb::MCTruth& mctruth) {
  if (!fTree->GetEntry(ientry)) {
    return false;
  }

  // Add MCParticles list with neutrino first
  for (int i=0; i<MCTruth_NParticles; i++) {
    // New MCParticle
    simb::MCParticle part(
      MCTruth_particles_TrackId[i],
      MCTruth_particles_PdgCode[i],
      "",
      MCTruth_particles_Mother[i],
      simb::MCParticle::s_uninitialized,
      MCTruth_particles_StatusCode[i]);

    // Add daughter track IDs for this particle
    for (int j=0; j<MCTruth_particles_NDaughters[i]; j++) {
      part.AddDaughter(MCTruth_particles_Daughters[i][j]);
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


bool NTupleInterface::FillGTruth(Long64_t ientry, simb::GTruth& gtruth) {
  if (!fTree->GetEntry(ientry)) {
    return false;
  }

  gtruth.fProbePDG = GTruth_ProbePDG;
  gtruth.fIsSeaQuark = GTruth_IsSeaQuark;
  gtruth.ftgtPDG = GTruth_tgtPDG;
  gtruth.gweight = GTruth_weight;
  gtruth.fprobability = GTruth_probability;
  gtruth.fXsec = GTruth_Xsec;
  gtruth.fDiffXsec = GTruth_fDiffXsec;
  gtruth.fVertex = TLorentzVector(
    GTruth_vertexX,
    GTruth_vertexY,
    GTruth_vertexZ,
    GTruth_vertexT);
  gtruth.fGscatter = GTruth_Gscatter;
  gtruth.fGint = GTruth_Gint;
  gtruth.fResNum = GTruth_ResNum;
  gtruth.fNumPiPlus = GTruth_NumPiPlus;
  gtruth.fNumPi0 = GTruth_NumPi0;
  gtruth.fNumPiMinus = GTruth_NumPiMinus;
  gtruth.fNumProton = GTruth_NumProton;
  gtruth.fNumNeutron = GTruth_NumNeutron;
  gtruth.fIsCharm = GTruth_IsCharm;
  gtruth.fgX = GTruth_gX;
  gtruth.fgY = GTruth_gY;
  gtruth.fgZ = GTruth_gZ;
  gtruth.fgT = GTruth_gT;
  gtruth.fgW = GTruth_gW;
  gtruth.fgQ2 = GTruth_gQ2;
  gtruth.fgq2 = GTruth_gq2;
  gtruth.fProbePDG = GTruth_ProbePDG;
  gtruth.fProbeP4 = TLorentzVector(
    GTruth_ProbeP4x,
    GTruth_ProbeP4y,
    GTruth_ProbeP4z,
    GTruth_ProbeP4E);
  gtruth.fHitNucP4 = TLorentzVector(
    GTruth_HitNucP4x,
    GTruth_HitNucP4y,
    GTruth_HitNucP4z,
    GTruth_HitNucP4E);
  gtruth.fFShadSystP4 = TLorentzVector(
    GTruth_FShadSystP4x,
    GTruth_FShadSystP4y,
    GTruth_FShadSystP4z,
    GTruth_FShadSystP4E);

  return true;
}

}  // namespace uboone

