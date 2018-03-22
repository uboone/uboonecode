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

  std::cout << "[NTUPLEINTERFACE] Tree retrieved. Setting branch addresses." << std::endl;

  // Metadata
  fTree->SetBranchAddress(branchDef.get<std::string>("run").c_str()                             , &run);
  fTree->SetBranchAddress(branchDef.get<std::string>("subrun").c_str()                          , &subrun);
  fTree->SetBranchAddress(branchDef.get<std::string>("event").c_str()                           , &event);

  // MCFlux
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_evtno").c_str()                    , &MCFlux_evtno);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_NuPosX").c_str()                   , &MCFlux_NuPosX);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_NuPosY").c_str()                   , &MCFlux_NuPosY);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_NuPosZ").c_str()                   , &MCFlux_NuPosZ);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_NuMomX").c_str()                   , &MCFlux_NuMomX);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_NuMomY").c_str()                   , &MCFlux_NuMomY);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_NuMomZ").c_str()                   , &MCFlux_NuMomZ);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_NuMomE").c_str()                   , &MCFlux_NuMomE);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_ntype").c_str()                    , &MCFlux_ntype);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_ptype").c_str()                    , &MCFlux_ptype);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_nimpwt").c_str()                   , &MCFlux_nimpwt);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_dk2gen").c_str()                   , &MCFlux_dk2gen);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_nenergyn").c_str()                 , &MCFlux_nenergyn);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_tpx").c_str()                      , &MCFlux_tpx);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_tpy").c_str()                      , &MCFlux_tpy);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_tpz").c_str()                      , &MCFlux_tpz);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_tptype").c_str()                   , &MCFlux_tptype);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_vx").c_str()                       , &MCFlux_vx);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_vy").c_str()                       , &MCFlux_vy);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCFlux_vz").c_str()                       , &MCFlux_vz);

  // MCTruth
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_NParticles").c_str()              , &MCTruth_NParticles);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_TrackId").c_str()       , &MCTruth_particles_TrackId);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_PdgCode").c_str()       , &MCTruth_particles_PdgCode);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_Mother").c_str()        , &MCTruth_particles_Mother);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_StatusCode").c_str()    , &MCTruth_particles_StatusCode);
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_particles_NumberDaughters").c_str()    , &MCTruth_particles_NumberDaughters);
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
  fTree->SetBranchAddress(branchDef.get<std::string>("MCTruth_neutrino_QSqr").c_str()             , &MCTruth_neutrino_QSqr);

  // GTruth
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_ProbePDG").c_str()                 , &GTruth_ProbePDG);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_IsSeaQuark").c_str()               , &GTruth_IsSeaQuark);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_tgtPDG").c_str()                   , &GTruth_tgtPDG);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_weight").c_str()                   , &GTruth_weight);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_probability").c_str()              , &GTruth_probability);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_Xsec").c_str()                     , &GTruth_Xsec);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_DiffXsec").c_str()                , &GTruth_DiffXsec);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_vertexX").c_str()                  , &GTruth_vertexX);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_vertexY").c_str()                  , &GTruth_vertexY);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_vertexZ").c_str()                  , &GTruth_vertexZ);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_vertexT").c_str()                  , &GTruth_vertexT);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_Gscatter").c_str()                 , &GTruth_Gscatter);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_Gint").c_str()                     , &GTruth_Gint);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_ResNum").c_str()                   , &GTruth_ResNum);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_NumPiPlus").c_str()                , &GTruth_NumPiPlus);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_NumPi0").c_str()                   , &GTruth_NumPi0);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_NumPiMinus").c_str()               , &GTruth_NumPiMinus);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_NumProton").c_str()                , &GTruth_NumProton);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_NumNeutron").c_str()               , &GTruth_NumNeutron);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_IsCharm").c_str()                  , &GTruth_IsCharm);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_gX").c_str()                       , &GTruth_gX);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_gY").c_str()                       , &GTruth_gY);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_gT").c_str()                       , &GTruth_gT);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_gW").c_str()                       , &GTruth_gW);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_gQ2").c_str()                      , &GTruth_gQ2);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_gq2").c_str()                      , &GTruth_gq2);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_ProbeP4x").c_str()                 , &GTruth_ProbeP4x);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_ProbeP4y").c_str()                 , &GTruth_ProbeP4y);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_ProbeP4z").c_str()                 , &GTruth_ProbeP4z);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_ProbeP4E").c_str()                 , &GTruth_ProbeP4E);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_HitNucP4x").c_str()                , &GTruth_HitNucP4x);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_HitNucP4y").c_str()                , &GTruth_HitNucP4y);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_HitNucP4z").c_str()                , &GTruth_HitNucP4z);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_HitNucP4E").c_str()                , &GTruth_HitNucP4E);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_FShadSystP4x").c_str()             , &GTruth_FShadSystP4x);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_FShadSystP4y").c_str()             , &GTruth_FShadSystP4y);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_FShadSystP4z").c_str()             , &GTruth_FShadSystP4z);
  fTree->SetBranchAddress(branchDef.get<std::string>("GTruth_FShadSystP4E").c_str()             , &GTruth_FShadSystP4E);

  fNEntries = fTree->GetEntries();
  assert(fNEntries > 0);
  fTree->GetEntry(0);
  fRun = run;

}


bool NTupleInterface::FillMCFlux(Long64_t ientry, simb::MCFlux& flux) {
  if (!fTree->GetEntry(ientry)) {
    std::cout << "[TEST] DIDN'T PASS" << std::endl;
    return false;
  }
  std::cout << "[TEST] DID PASS" << std::endl;

  fNuPos = TLorentzVector(
    MCFlux_NuPosX,
    MCFlux_NuPosY,
    MCFlux_NuPosZ,
    0);

  fNuMom = TLorentzVector(
    MCFlux_NuMomX,
    MCFlux_NuMomY,
    MCFlux_NuMomZ,
    MCFlux_NuMomE);

  flux.fptype     = MCFlux_ptype;
  flux.fntype     = MCFlux_ntype;
  flux.fnimpwt    = MCFlux_nimpwt;
  flux.fdk2gen    = MCFlux_dk2gen;
  flux.fnenergyn  = MCFlux_nenergyn;
  flux.ftpx       = MCFlux_tpx;
  flux.ftpy       = MCFlux_tpy;
  flux.ftpz       = MCFlux_tpz;
  flux.fvx        = MCFlux_vx;
  flux.fvy        = MCFlux_vy;
  flux.fvz        = MCFlux_vz;
  flux.ftptype    = MCFlux_tptype;

  return true;
}


bool NTupleInterface::FillMCTruth(Long64_t ientry, simb::MCTruth& mctruth) {
  std::cout << "[NTUPLEINTERFACE] Filling MCTruth" << std::endl; 
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
    for (int j=0; j<MCTruth_particles_NumberDaughters[i]; j++) {
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
                      MCTruth_neutrino_QSqr);

  return true;
}  


bool NTupleInterface::FillGTruth(Long64_t ientry, simb::GTruth& gtruth) {
  std::cout << "[NTUPLEINTERFACE] Filling GTruth" << std::endl; 
  if (!fTree->GetEntry(ientry)) {
    return false;
  }

  gtruth.fProbePDG = GTruth_ProbePDG;
  gtruth.fIsSeaQuark = GTruth_IsSeaQuark;
  gtruth.ftgtPDG = GTruth_tgtPDG;
  gtruth.fweight = GTruth_weight;
  gtruth.fprobability = GTruth_probability;
  gtruth.fXsec = GTruth_Xsec;
  gtruth.fDiffXsec = GTruth_DiffXsec;
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
  gtruth.fgT = GTruth_gT;
  gtruth.fgW = GTruth_gW;
  gtruth.fgQ2 = GTruth_gQ2;
  gtruth.fgq2 = GTruth_gq2;
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

