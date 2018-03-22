#include "EventWeightTreeUtility.h"

namespace uboone {
  int WriteTree(const art::Event e, 
      const art::Ptr<simb::MCFlux> mcFlux, 
      const art::Ptr<simb::MCTruth> mcTruth, 
      const art::Ptr<simb::GTruth> gTruth){

    // create TFile
    TFile* file = new TFile("eventWeightFile.root", "UPDATE");
    file->cd();

    // create tree
    if (!file->Get("tree")){
      tree = new TTree("tree", "tree");
      InitialiseTree(tree);
    }
    else{
      tree = (TTree*)file->Get("tree");
      ConfigureTree(tree);
    }

    const simb::MCParticle& nu = mcTruth->GetNeutrino().Nu();

    run = e.run();
    subrun = e.subRun();
    event = e.event();

    // Fill MCFlux information
    MCFlux_evtno     = mcFlux->fevtno;
    MCFlux_NuPosX    = nu.Vx();
    MCFlux_NuPosY    = nu.Vy();
    MCFlux_NuPosZ    = nu.Vz();
    MCFlux_NuMomX    = nu.Px(); 
    MCFlux_NuMomY    = nu.Py(); 
    MCFlux_NuMomZ    = nu.Pz(); 
    MCFlux_NuMomE    = nu.E();
    MCFlux_genx      = mcFlux->fgenx;
    MCFlux_geny      = mcFlux->fgeny;
    MCFlux_genz      = mcFlux->fgenz;
    MCFlux_ntype     = mcFlux->fntype;
    MCFlux_ptype     = mcFlux->fptype;
    MCFlux_nimpwt    = mcFlux->fnimpwt;
    MCFlux_dk2gen    = mcFlux->fdk2gen;
    MCFlux_nenergyn  = mcFlux->fnenergyn;
    MCFlux_tpx       = mcFlux->ftpx;
    MCFlux_tpy       = mcFlux->ftpy;
    MCFlux_tpz       = mcFlux->ftpz;
    MCFlux_tptype    = mcFlux->ftptype;
    MCFlux_vx        = mcFlux->fvx;
    MCFlux_vy        = mcFlux->fvy;
    MCFlux_vz        = mcFlux->fvz;

    // loop MCParticle info for MCTruth object

    MCTruth_NParticles = mcTruth->NParticles();

    for (int i = 0; i < MCTruth_NParticles; i++){

      const simb::MCParticle& mcParticle = mcTruth->GetParticle(i);

      MCTruth_particles_TrackId[i] = mcParticle.TrackId();
      MCTruth_particles_PdgCode[i] = mcParticle.PdgCode();
      MCTruth_particles_Mother[i]  = mcParticle.Mother();
      MCTruth_particles_StatusCode[i] = mcParticle.StatusCode();
      MCTruth_particles_NumberDaughters[i] = mcParticle.NumberDaughters();

      for (int j = 0; j < MCTruth_particles_NumberDaughters[i]; j++){

        const simb::MCParticle& daughterMcParticle = mcTruth->GetParticle(j);
        MCTruth_particles_Daughters[i][j] = daughterMcParticle.TrackId();

      }

      MCTruth_particles_Gvx[i] = mcParticle.Gvx();
      MCTruth_particles_Gvy[i] = mcParticle.Gvy();
      MCTruth_particles_Gvz[i] = mcParticle.Gvz();
      MCTruth_particles_Gvt[i] = mcParticle.Gvt();
      MCTruth_particles_px0[i] = mcParticle.Px(0);
      MCTruth_particles_py0[i] = mcParticle.Py(0);
      MCTruth_particles_pz0[i] = mcParticle.Pz(0);
      MCTruth_particles_e0[i] = mcParticle.E(0);
      MCTruth_particles_Rescatter[i] = mcParticle.Rescatter();
      MCTruth_particles_polx[i] = mcParticle.Polarization().X();
      MCTruth_particles_poly[i] = mcParticle.Polarization().Y();
      MCTruth_particles_polz[i] = mcParticle.Polarization().Z();
    }

    const simb::MCNeutrino& mcNeutrino = mcTruth->GetNeutrino();

    MCTruth_neutrino_CCNC = mcNeutrino.CCNC();
    MCTruth_neutrino_mode = mcNeutrino.Mode();
    MCTruth_neutrino_interactionType = mcNeutrino.InteractionType();
    MCTruth_neutrino_target = mcNeutrino.Target();
    MCTruth_neutrino_nucleon = mcNeutrino.HitNuc();
    MCTruth_neutrino_quark = mcNeutrino.HitQuark();
    MCTruth_neutrino_W = mcNeutrino.W();
    MCTruth_neutrino_X = mcNeutrino.X();
    MCTruth_neutrino_Y = mcNeutrino.Y();
    MCTruth_neutrino_QSqr = mcNeutrino.QSqr();

    GTruth_IsSeaQuark = gTruth->fIsSeaQuark;
    GTruth_tgtPDG = gTruth->ftgtPDG;
    GTruth_weight = gTruth->fweight;
    GTruth_probability = gTruth->fprobability;
    GTruth_Xsec = gTruth->fXsec;
    GTruth_DiffXsec = gTruth->fDiffXsec;
    GTruth_vertexX = gTruth->fVertex.X();
    GTruth_vertexY = gTruth->fVertex.Y();
    GTruth_vertexZ = gTruth->fVertex.Z();
    GTruth_vertexT = gTruth->fVertex.T();
    GTruth_Gscatter = gTruth->fGscatter;
    GTruth_Gint = gTruth->fGint;
    GTruth_ResNum = gTruth->fResNum;
    GTruth_NumPiPlus = gTruth->fNumPiPlus;
    GTruth_NumPi0 = gTruth->fNumPi0;
    GTruth_NumPiMinus = gTruth->fNumPiMinus;
    GTruth_NumProton = gTruth->fNumProton;
    GTruth_NumNeutron = gTruth->fNumNeutron;
    GTruth_IsCharm = gTruth->fIsCharm;
    GTruth_gX = gTruth->fgX;
    GTruth_gY = gTruth->fgY;
    GTruth_gT = gTruth->fgT;
    GTruth_gW = gTruth->fgW;
    GTruth_gQ2 = gTruth->fgQ2;
    GTruth_gq2 = gTruth->fgq2;
    GTruth_ProbePDG = gTruth->fProbePDG;
    tree->Fill();
    tree->Write("tree", TObject::kOverwrite);
    file->Close();

    return 0;
  }

  int WriteTree(const art::Event e, 
      const art::Ptr<simb::MCFlux> mcFlux, 
      const art::Ptr<simb::MCTruth> mcTruth, 
      const art::Ptr<simb::GTruth> gTruth
      TTree* t){

    tree = t;
    WriteTree(e, mcFlux, mcTruth, gTruth);

  }

  void InitialiseTree(TTree* t){

    t->Branch("run"     , &run);
    t->Branch("subrun"  , &subrun);
    t->Branch("event"   , &event);

    // MCFlux
    t->Branch("MCFlux_evtno"   , &MCFlux_evtno);
    t->Branch("MCFlux_NuPosX"  , &MCFlux_NuPosX);
    t->Branch("MCFlux_NuPosY"  , &MCFlux_NuPosY);
    t->Branch("MCFlux_NuPosZ"  , &MCFlux_NuPosZ);
    t->Branch("MCFlux_NuMomX"  , &MCFlux_NuMomX);
    t->Branch("MCFlux_NuMomY"  , &MCFlux_NuMomY);
    t->Branch("MCFlux_NuMomZ"  , &MCFlux_NuMomZ);
    t->Branch("MCFlux_NuMomE"  , &MCFlux_NuMomE);
    t->Branch("MCFlux_genx"    , &MCFlux_genx);
    t->Branch("MCFlux_geny"    , &MCFlux_geny);
    t->Branch("MCFlux_genz"    , &MCFlux_genz);
    t->Branch("MCFlux_ntype"   , &MCFlux_ntype);
    t->Branch("MCFlux_ptype"   , &MCFlux_ptype);
    t->Branch("MCFlux_nimpwt"  , &MCFlux_nimpwt);
    t->Branch("MCFlux_dk2gen"  , &MCFlux_dk2gen);
    t->Branch("MCFlux_nenergyn", &MCFlux_nenergyn);
    t->Branch("MCFlux_tpx"     , &MCFlux_tpx);
    t->Branch("MCFlux_tpy"     , &MCFlux_tpy);
    t->Branch("MCFlux_tpz"     , &MCFlux_tpz);
    t->Branch("MCFlux_tptype"  , &MCFlux_tptype);
    t->Branch("MCFlux_vx"      , &MCFlux_vx);
    t->Branch("MCFlux_vy"      , &MCFlux_vy);
    t->Branch("MCFlux_vz"      , &MCFlux_vz);

    // MCTruth
    t->Branch("MCTruth_NParticles", &MCTruth_NParticles);
    t->Branch("MCTruth_particles_TrackId", MCTruth_particles_TrackId, "MCTruth_particles_TrackId[MCTruth_NParticles]/I");
    t->Branch("MCTruth_particles_PdgCode", &MCTruth_particles_PdgCode, "MCTruth_particles_PdgCode[MCTruth_NParticles]/I");
    t->Branch("MCTruth_particles_Mother", &MCTruth_particles_Mother, "MCTruth_particles_Mother[MCTruth_NParticles]/I"); 
    t->Branch("MCTruth_particles_StatusCode", MCTruth_particles_StatusCode, "MCTruth_particles_StatusCode[MCTruth_NParticles]/I"); 
    t->Branch("MCTruth_particles_NumberDaughters", MCTruth_particles_NumberDaughters, "MCTruth_particles_NumberDaughters[MCTruth_NParticles]/I");
    t->Branch("MCTruth_particles_Daughters", MCTruth_particles_Daughters, "MCTruth_particles_Daughters[MCTruth_NParticles][100]/I");
    t->Branch("MCTruth_particles_Gvx", MCTruth_particles_Gvx, "MCTruth_particles_Gvx[MCTruth_NParticles]/D");
    t->Branch("MCTruth_particles_Gvy", MCTruth_particles_Gvy, "MCTruth_particles_Gvy[MCTruth_NParticles]/D");
    t->Branch("MCTruth_particles_Gvz", MCTruth_particles_Gvz, "MCTruth_particles_Gvz[MCTruth_NParticles]/D");
    t->Branch("MCTruth_particles_Gvt", MCTruth_particles_Gvt, "MCTruth_particles_Gvt[MCTruth_NParticles]/D");
    t->Branch("MCTruth_particles_px0", MCTruth_particles_px0, "MCTruth_particles_px0[MCTruth_NParticles]/D");
    t->Branch("MCTruth_particles_py0", MCTruth_particles_py0, "MCTruth_particles_py0[MCTruth_NParticles]/D");
    t->Branch("MCTruth_particles_pz0", MCTruth_particles_pz0, "MCTruth_particles_pz0[MCTruth_NParticles]/D");
    t->Branch("MCTruth_particles_e0", MCTruth_particles_e0, "MCTruth_particles_e0[MCTruth_NParticles]/D");
    t->Branch("MCTruth_particles_Rescatter", MCTruth_particles_Rescatter, "MCTruth_particles_Rescatter[MCTruth_NParticles]/I");
    t->Branch("MCTruth_particles_polx", MCTruth_particles_polx, "MCTruth_particles_polx[MCTruth_NParticles]/D");
    t->Branch("MCTruth_particles_poly", MCTruth_particles_poly, "MCTruth_particles_poly[MCTruth_NParticles]/D");
    t->Branch("MCTruth_particles_polz", MCTruth_particles_polz, "MCTruth_particles_polz[MCTruth_NParticles]/D");
    t->Branch("MCTruth_neutrino_CCNC", &MCTruth_neutrino_CCNC);
    t->Branch("MCTruth_neutrino_mode", &MCTruth_neutrino_mode);
    t->Branch("MCTruth_neutrino_interactionType", &MCTruth_neutrino_interactionType);
    t->Branch("MCTruth_neutrino_target", &MCTruth_neutrino_target);
    t->Branch("MCTruth_neutrino_nucleon", &MCTruth_neutrino_nucleon);
    t->Branch("MCTruth_neutrino_quark", &MCTruth_neutrino_quark);
    t->Branch("MCTruth_neutrino_W", &MCTruth_neutrino_W);
    t->Branch("MCTruth_neutrino_X", &MCTruth_neutrino_X);
    t->Branch("MCTruth_neutrino_Y", &MCTruth_neutrino_Y);
    t->Branch("MCTruth_neutrino_QSqr",&MCTruth_neutrino_QSqr);

    //GTruth
    t->Branch("GTruth_IsSeaQuark"   , &GTruth_IsSeaQuark);
    t->Branch("GTruth_tgtPDG"       , &GTruth_tgtPDG);
    t->Branch("GTruth_weight"       , &GTruth_weight);
    t->Branch("GTruth_probability"  , &GTruth_probability);
    t->Branch("GTruth_Xsec"         , &GTruth_Xsec);
    t->Branch("GTruth_DiffXsec"     , &GTruth_DiffXsec);
    t->Branch("GTruth_vertexX"      , &GTruth_vertexX);
    t->Branch("GTruth_vertexY"      , &GTruth_vertexY);
    t->Branch("GTruth_vertexZ"      , &GTruth_vertexZ);
    t->Branch("GTruth_vertexT"      , &GTruth_vertexT);
    t->Branch("GTruth_Gscatter"     , &GTruth_Gscatter);
    t->Branch("GTruth_Gint"         , &GTruth_Gint);
    t->Branch("GTruth_ResNum"       , &GTruth_ResNum);
    t->Branch("GTruth_NumPiPlus"    , &GTruth_NumPiPlus);
    t->Branch("GTruth_NumPi0"       , &GTruth_NumPiMinus);
    t->Branch("GTruth_NumPiMinus"   , &GTruth_NumPiMinus);
    t->Branch("GTruth_NumProton"    , &GTruth_NumProton);
    t->Branch("GTruth_NumNeutron"   , &GTruth_NumNeutron);
    t->Branch("GTruth_IsCharm"      , &GTruth_IsCharm);
    t->Branch("GTruth_gX"           , &GTruth_gX);
    t->Branch("GTruth_gY"           , &GTruth_gY);
    t->Branch("GTruth_gT"           , &GTruth_gT);
    t->Branch("GTruth_gW"           , &GTruth_gW);
    t->Branch("GTruth_gQ2"          , &GTruth_gQ2);
    t->Branch("GTruth_gq2"          , &GTruth_gq2);
    t->Branch("GTruth_ProbePDG"     , &GTruth_ProbePDG);
    t->Branch("GTruth_ProbeP4x"     , &GTruth_ProbeP4x);
    t->Branch("GTruth_ProbeP4y"     , &GTruth_ProbeP4y);
    t->Branch("GTruth_ProbeP4z"     , &GTruth_ProbeP4z);
    t->Branch("GTruth_ProbeP4E"     , &GTruth_ProbeP4E);
    t->Branch("GTruth_HitNucP4x"    , &GTruth_HitNucP4x);
    t->Branch("GTruth_HitNucP4y"    , &GTruth_HitNucP4y);
    t->Branch("GTruth_HitNucP4z"    , &GTruth_HitNucP4z);
    t->Branch("GTruth_HitNucP4E"    , &GTruth_HitNucP4E);
    t->Branch("GTruth_FShadSystP4x" , &GTruth_FShadSystP4x);
    t->Branch("GTruth_FShadSystP4y" , &GTruth_FShadSystP4y);
    t->Branch("GTruth_FShadSystP4z" , &GTruth_FShadSystP4z);
    t->Branch("GTruth_FShadSystP4E" , &GTruth_FShadSystP4E);


  }

  void ConfigureTree(TTree* t){
    t->SetBranchAddress("run"     , &run);
    t->SetBranchAddress("subrun"  , &subrun);
    t->SetBranchAddress("event"   , &event);

    // MCFlux
    t->SetBranchAddress("MCFlux_evtno"   , &MCFlux_evtno);
    t->SetBranchAddress("MCFlux_NuPosX"  , &MCFlux_NuPosX);
    t->SetBranchAddress("MCFlux_NuPosY"  , &MCFlux_NuPosY);
    t->SetBranchAddress("MCFlux_NuPosZ"  , &MCFlux_NuPosZ);
    t->SetBranchAddress("MCFlux_NuMomX"  , &MCFlux_NuMomX);
    t->SetBranchAddress("MCFlux_NuMomY"  , &MCFlux_NuMomY);
    t->SetBranchAddress("MCFlux_NuMomZ"  , &MCFlux_NuMomZ);
    t->SetBranchAddress("MCFlux_NuMomE"  , &MCFlux_NuMomE);
    t->SetBranchAddress("MCFlux_genx"    , &MCFlux_genx);
    t->SetBranchAddress("MCFlux_geny"    , &MCFlux_geny);
    t->SetBranchAddress("MCFlux_genz"    , &MCFlux_genz);
    t->SetBranchAddress("MCFlux_ntype"   , &MCFlux_ntype);
    t->SetBranchAddress("MCFlux_ptype"   , &MCFlux_ptype);
    t->SetBranchAddress("MCFlux_nimpwt"  , &MCFlux_nimpwt);
    t->SetBranchAddress("MCFlux_dk2gen"  , &MCFlux_dk2gen);
    t->SetBranchAddress("MCFlux_nenergyn", &MCFlux_nenergyn);
    t->SetBranchAddress("MCFlux_tpx"     , &MCFlux_tpx);
    t->SetBranchAddress("MCFlux_tpy"     , &MCFlux_tpy);
    t->SetBranchAddress("MCFlux_tpz"     , &MCFlux_tpz);
    t->SetBranchAddress("MCFlux_tptype"  , &MCFlux_tptype);
    t->SetBranchAddress("MCFlux_vx"      , &MCFlux_vx);
    t->SetBranchAddress("MCFlux_vy"      , &MCFlux_vy);
    t->SetBranchAddress("MCFlux_vz"      , &MCFlux_vz);

    // MCTruth
    t->SetBranchAddress("MCTruth_NParticles", &MCTruth_NParticles);
    t->SetBranchAddress("MCTruth_particles_TrackId", MCTruth_particles_TrackId, "MCTruth_particles_TrackId[MCTruth_NParticles]/I");
    t->SetBranchAddress("MCTruth_particles_PdgCode", &MCTruth_particles_PdgCode, "MCTruth_particles_PdgCode[MCTruth_NParticles]/I");
    t->SetBranchAddress("MCTruth_particles_Mother", &MCTruth_particles_Mother, "MCTruth_particles_Mother[MCTruth_NParticles]/I"); 
    t->SetBranchAddress("MCTruth_particles_StatusCode", MCTruth_particles_StatusCode, "MCTruth_particles_StatusCode[MCTruth_NParticles]/I"); 
    t->SetBranchAddress("MCTruth_particles_NumberDaughters", MCTruth_particles_NumberDaughters, "MCTruth_particles_NumberDaughters[MCTruth_NParticles]/I");
    t->SetBranchAddress("MCTruth_particles_Daughters", MCTruth_particles_Daughters, "MCTruth_particles_Daughters[MCTruth_NParticles][100]/I");
    t->SetBranchAddress("MCTruth_particles_Gvx", MCTruth_particles_Gvx, "MCTruth_particles_Gvx[MCTruth_NParticles]/D");
    t->SetBranchAddress("MCTruth_particles_Gvy", MCTruth_particles_Gvy, "MCTruth_particles_Gvy[MCTruth_NParticles]/D");
    t->SetBranchAddress("MCTruth_particles_Gvz", MCTruth_particles_Gvz, "MCTruth_particles_Gvz[MCTruth_NParticles]/D");
    t->SetBranchAddress("MCTruth_particles_Gvt", MCTruth_particles_Gvt, "MCTruth_particles_Gvt[MCTruth_NParticles]/D");
    t->SetBranchAddress("MCTruth_particles_px0", MCTruth_particles_px0, "MCTruth_particles_px0[MCTruth_NParticles]/D");
    t->SetBranchAddress("MCTruth_particles_py0", MCTruth_particles_py0, "MCTruth_particles_py0[MCTruth_NParticles]/D");
    t->SetBranchAddress("MCTruth_particles_pz0", MCTruth_particles_pz0, "MCTruth_particles_pz0[MCTruth_NParticles]/D");
    t->SetBranchAddress("MCTruth_particles_e0", MCTruth_particles_e0, "MCTruth_particles_e0[MCTruth_NParticles]/D");
    t->SetBranchAddress("MCTruth_particles_Rescatter", MCTruth_particles_Rescatter, "MCTruth_particles_Rescatter[MCTruth_NParticles]/I");
    t->SetBranchAddress("MCTruth_particles_polx", MCTruth_particles_polx, "MCTruth_particles_polx[MCTruth_NParticles]/D");
    t->SetBranchAddress("MCTruth_particles_poly", MCTruth_particles_poly, "MCTruth_particles_poly[MCTruth_NParticles]/D");
    t->SetBranchAddress("MCTruth_particles_polz", MCTruth_particles_polz, "MCTruth_particles_polz[MCTruth_NParticles]/D");
    t->SetBranchAddress("MCTruth_neutrino_CCNC", &MCTruth_neutrino_CCNC);
    t->SetBranchAddress("MCTruth_neutrino_mode", &MCTruth_neutrino_mode);
    t->SetBranchAddress("MCTruth_neutrino_interactionType", &MCTruth_neutrino_interactionType);
    t->SetBranchAddress("MCTruth_neutrino_target", &MCTruth_neutrino_target);
    t->SetBranchAddress("MCTruth_neutrino_nucleon", &MCTruth_neutrino_nucleon);
    t->SetBranchAddress("MCTruth_neutrino_quark", &MCTruth_neutrino_quark);
    t->SetBranchAddress("MCTruth_neutrino_W", &MCTruth_neutrino_W);
    t->SetBranchAddress("MCTruth_neutrino_X", &MCTruth_neutrino_X);
    t->SetBranchAddress("MCTruth_neutrino_Y", &MCTruth_neutrino_Y);
    t->SetBranchAddress("MCTruth_neutrino_QSqr",&MCTruth_neutrino_QSqr);

    //GTruth
    //t->SetBranchAddress("GTruth_ProbePDG"   , &GTruth_ProbePDG);
    t->SetBranchAddress("GTruth_IsSeaQuark"   , &GTruth_IsSeaQuark);
    t->SetBranchAddress("GTruth_tgtPDG"       , &GTruth_tgtPDG);
    t->SetBranchAddress("GTruth_weight"       , &GTruth_weight);
    t->SetBranchAddress("GTruth_probability"  , &GTruth_probability);
    t->SetBranchAddress("GTruth_Xsec"         , &GTruth_Xsec);
    t->SetBranchAddress("GTruth_DiffXsec"     , &GTruth_DiffXsec);
    t->SetBranchAddress("GTruth_vertexX"      , &GTruth_vertexX);
    t->SetBranchAddress("GTruth_vertexY"      , &GTruth_vertexY);
    t->SetBranchAddress("GTruth_vertexZ"      , &GTruth_vertexZ);
    t->SetBranchAddress("GTruth_vertexT"      , &GTruth_vertexT);
    t->SetBranchAddress("GTruth_Gscatter"     , &GTruth_Gscatter);
    t->SetBranchAddress("GTruth_Gint"         , &GTruth_Gint);
    t->SetBranchAddress("GTruth_ResNum"       , &GTruth_ResNum);
    t->SetBranchAddress("GTruth_NumPiPlus"    , &GTruth_NumPiPlus);
    t->SetBranchAddress("GTruth_NumPi0"       , &GTruth_NumPiMinus);
    t->SetBranchAddress("GTruth_NumPiMinus"   , &GTruth_NumPiMinus);
    t->SetBranchAddress("GTruth_NumProton"    , &GTruth_NumProton);
    t->SetBranchAddress("GTruth_NumNeutron"   , &GTruth_NumNeutron);
    t->SetBranchAddress("GTruth_IsCharm"      , &GTruth_IsCharm);
    t->SetBranchAddress("GTruth_gX"           , &GTruth_gX);
    t->SetBranchAddress("GTruth_gY"           , &GTruth_gY);
    t->SetBranchAddress("GTruth_gT"           , &GTruth_gT);
    t->SetBranchAddress("GTruth_gW"           , &GTruth_gW);
    t->SetBranchAddress("GTruth_gQ2"          , &GTruth_gQ2);
    t->SetBranchAddress("GTruth_gq2"          , &GTruth_gq2);
    t->SetBranchAddress("GTruth_ProbePDG"     , &GTruth_ProbePDG);
    t->SetBranchAddress("GTruth_ProbeP4x"     , &GTruth_ProbeP4x);
    t->SetBranchAddress("GTruth_ProbeP4y"     , &GTruth_ProbeP4y);
    t->SetBranchAddress("GTruth_ProbeP4z"     , &GTruth_ProbeP4z);
    t->SetBranchAddress("GTruth_ProbeP4E"     , &GTruth_ProbeP4E);
    t->SetBranchAddress("GTruth_HitNucP4x"    , &GTruth_HitNucP4x);
    t->SetBranchAddress("GTruth_HitNucP4y"    , &GTruth_HitNucP4y);
    t->SetBranchAddress("GTruth_HitNucP4z"    , &GTruth_HitNucP4z);
    t->SetBranchAddress("GTruth_HitNucP4E"    , &GTruth_HitNucP4E);
    t->SetBranchAddress("GTruth_FShadSystP4x" , &GTruth_FShadSystP4x);
    t->SetBranchAddress("GTruth_FShadSystP4y" , &GTruth_FShadSystP4y);
    t->SetBranchAddress("GTruth_FShadSystP4z" , &GTruth_FShadSystP4z);
    t->SetBranchAddress("GTruth_FShadSystP4E" , &GTruth_FShadSystP4E);

  }

} // namespace uboone
