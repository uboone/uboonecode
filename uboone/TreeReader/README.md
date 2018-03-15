Tree Reader
===========
This module reads in data from a (semi-)arbitrary ROOT tree, to populate
an `art::Event`.

The primary goal is to load in truth information from post-selection analysis
trees and calculate weights with EventWeight, without need to go back to the
the original art ROOT files.

The initial implementation of this module reads in MCTruth, MCFlux, and GTruth
objects, which are required for flux and cross section reweighting. The output
MCEventWeights can then be attached to the analysis trees, for propagation.

NTuple Input
------------
In Ntuple input mode, the TreeReader source accepts a ROOT file, where the
user specifies the name of the ROOT tree and a list of branch names where the
required data is stored.

An example configuration, where the names of the branches in the input TTree
happen to match the key names, which are fixed.

    source: {
      fileNames: []
      module_type: TreeReader
      skipEvents: 0
      maxEvents: -1
      inputType: "ntuple"
      treeName: "NAME OF TREE"
      branches: {
        vtxx: "vtxx"           # neutrino position at the window (fgenx in MCFlux)
        vtxy: "vtxy" 
        vtxz: "vtxz"
        px: "px"               # neutrino momentum (taken from Neutrino MCParticle)
        py: "py"
        pz: "pz"
        E: "E"                 # neutrino energy
        pdg: "pdg"             # neutrino pdg
        ptype: "ptype"         # parent pdg
        wgt: "wgt"             # weight (always set to 1 for files generated with gSimple)
        dist: "dist"
        evtno: "evtno"         # POT
        nenergyn: "nenergyn"   # neutrino energy for a neutrino forced at the center of the near detector
        tpx: "tpx"             # parent momentum exiting the target
        tpy: "tpy"
        tpz: "tpz"
        vx: "vx"               # exit point at the target
        vy: "vy"
        vz: "vz"
        tptype: "tptype"       # parent particle ID exiting the target
        }
    }

Note that the tree can be inside a `TDirectory`, e.g. `dir/treename`.


    int MCTruth_neutrino_CCNC;
    int MCTruth_neutrino_mode;
    int MCTruth_neutrino_interactionType;
    int MCTruth_neutrino_target;
    int MCTruth_neutrino_nucleon;
    int MCTruth_neutrino_quark;
    double MCTruth_neutrino_W;
    double MCTruth_neutrino_X;
    double MCTruth_neutrino_Y;
    double MCTruth_neutrino_Q2;

    int MCTruth_NParticles;
    int MCTruth_particles_TrackId[kMaxParticles];
    int MCTruth_particles_PdgCode[kMaxParticles];
    int MCTruth_particles_Mother[kMaxParticles];
    int MCTruth_particles_StatusCode[kMaxParticles];
    int MCTruth_particles_NDaughters[kMaxParticles];
    int MCTruth_particles_Daughters[kMaxParticles][kMaxParticles];
    double MCTruth_particles_Gvx[kMaxParticles];
    double MCTruth_particles_Gvy[kMaxParticles];
    double MCTruth_particles_Gvz[kMaxParticles];
    double MCTruth_particles_Gvt[kMaxParticles];
    double MCTruth_particles_px0[kMaxParticles];
    double MCTruth_particles_py0[kMaxParticles];
    double MCTruth_particles_pz0[kMaxParticles];
    double MCTruth_particles_e0[kMaxParticles];
    int MCTruth_particles_Rescatter[kMaxParticles];
    double MCTruth_particles_polx[kMaxParticles];
    double MCTruth_particles_poly[kMaxParticles];
    double MCTruth_particles_polz[kMaxParticles];

Authors
-------
This module is based on the FluxReader by Zarko Pavlovic. It was generalized
to handle other analysis tree inputs by Adam Lister and Andy Mastbaum.

