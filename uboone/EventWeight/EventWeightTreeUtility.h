#ifndef EVENTWEIGHTTREEUTILITY_H
#define EVENTWEIGHTTREEUTILITY_H

/** 
 * Utilities to be used to EventWeighting purposes.
 *  
 * \author A. Lister 
 */

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"

#include "TFile.h"
#include "TTree.h"
#include "TObject.h"

namespace uboone {

  /**
   *  \class EWTreeUtil
   *  \brief Automatically builds TTree with variables for reweighting
   */
  class EWTreeUtil {
    public:

      /** Constuctor */
      EWTreeUtil();

      /** Destructor */
      ~EWTreeUtil();

      /**
       * Take variables from mcFlux, mcTruth, gTruth
       * and dump them to new TTree for later reweighting
       */
      int WriteTree(const art::Event & e, 
        const art::Ptr<simb::MCFlux> mcFlux, 
        const art::Ptr<simb::MCTruth> mcTruth, 
        const art::Ptr<simb::GTruth> gTruth);

      /**
       * Take variables from mcFlux, mcTruth, gTruth
       * and dump them to user specified TTree for later
       * reweighting
       */
      int WriteTree(const art::Event & e, 
          const art::Ptr<simb::MCFlux> mcFlux, 
          const art::Ptr<simb::MCTruth> mcTruth, 
          const art::Ptr<simb::GTruth> gTruth,
          TTree* t);

      /**
       * If TTree doesn't already exist then create tree 
       */
      void InitialiseTree(TTree* t);
      
      /**
       * If TTree already exists, then set branch addresses
       */
      void ConfigureTree(TTree* t);

    private:
      TTree* tree;
      bool isDebug;

      static const unsigned kMaxMCParticles = 50;

      // Metadata
      int run;
      int subrun;
      int event;

      // MCFlux
      int MCFlux_evtno;
      double MCFlux_NuPosX, MCFlux_NuPosY, MCFlux_NuPosZ;
      double MCFlux_NuMomX, MCFlux_NuMomY, MCFlux_NuMomZ, MCFlux_NuMomE;
      double MCFlux_genx, MCFlux_geny, MCFlux_genz;
      int MCFlux_ntype;
      int MCFlux_ptype;
      double MCFlux_nimpwt;
      double MCFlux_dk2gen;
      double MCFlux_nenergyn;
      double MCFlux_tpx, MCFlux_tpy, MCFlux_tpz;
      double MCFlux_vx, MCFlux_vy, MCFlux_vz;
      int MCFlux_tptype;

      // MCTruth
      int MCTruth_NParticles;
      int MCTruth_particles_TrackId[kMaxMCParticles];
      int MCTruth_particles_PdgCode[kMaxMCParticles];
      int MCTruth_particles_Mother[kMaxMCParticles];
      int MCTruth_particles_StatusCode[kMaxMCParticles];
      int MCTruth_particles_NumberDaughters[kMaxMCParticles];
      // only one multi-dimensional array allowed... thanks root.
      int MCTruth_particles_Daughters[kMaxMCParticles][100]; 
      double MCTruth_particles_Gvx[kMaxMCParticles];
      double MCTruth_particles_Gvy[kMaxMCParticles];
      double MCTruth_particles_Gvz[kMaxMCParticles];
      double MCTruth_particles_Gvt[kMaxMCParticles];
      double MCTruth_particles_px0[kMaxMCParticles];
      double MCTruth_particles_py0[kMaxMCParticles];
      double MCTruth_particles_pz0[kMaxMCParticles];
      double MCTruth_particles_e0[kMaxMCParticles];
      int MCTruth_particles_Rescatter[kMaxMCParticles];
      double MCTruth_particles_polx[kMaxMCParticles];
      double MCTruth_particles_poly[kMaxMCParticles];
      double MCTruth_particles_polz[kMaxMCParticles];
      int MCTruth_neutrino_CCNC;
      int MCTruth_neutrino_mode;
      int MCTruth_neutrino_interactionType;
      int MCTruth_neutrino_target;
      int MCTruth_neutrino_nucleon;
      int MCTruth_neutrino_quark;
      double MCTruth_neutrino_QSqr;
      double MCTruth_neutrino_W;
      double MCTruth_neutrino_X;
      double MCTruth_neutrino_Y;

      // GTruth
      bool GTruth_IsSeaQuark;
      int GTruth_tgtPDG;
      double GTruth_weight;
      double GTruth_probability;
      double GTruth_Xsec;
      double GTruth_DiffXsec;
      double GTruth_vertexX;
      double GTruth_vertexY;
      double GTruth_vertexZ;
      double GTruth_vertexT;
      int GTruth_Gscatter;
      int GTruth_Gint;
      int GTruth_ResNum;
      int GTruth_NumPiPlus;
      int GTruth_NumPi0;
      int GTruth_NumPiMinus;
      int GTruth_NumProton;
      int GTruth_NumNeutron;
      bool GTruth_IsCharm;
      double GTruth_gX;
      double GTruth_gY;
      //double GTruth_gZ;
      double GTruth_gT;
      double GTruth_gW;
      double GTruth_gQ2;
      double GTruth_gq2;
      int GTruth_ProbePDG;
      double GTruth_ProbeP4x;
      double GTruth_ProbeP4y;
      double GTruth_ProbeP4z;
      double GTruth_ProbeP4E;
      double GTruth_HitNucP4x;
      double GTruth_HitNucP4y;
      double GTruth_HitNucP4z;
      double GTruth_HitNucP4E;
      double GTruth_FShadSystP4x;
      double GTruth_FShadSystP4y;
      double GTruth_FShadSystP4z;
      double GTruth_FShadSystP4E;
  };

} // namespace uboone

#endif // EVENTWEIGHTTREEUTILITY


