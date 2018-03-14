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
        vtxx: "vtxx"
        vtxy: "vtxy"
        vtxz: "vtxz"
        px: "px"
        py: "py"
        pz: "pz"
        E: "E"
        pdg: "pdg"
        ptype: "ptype"
        wgt: "wgt"
        dist: "dist"
        evtno: "evtno"
        nenergyn: "nenergyn"
      }
    }

Note that the tree can be inside a `TDirectory`, e.g. `dir/treename`.

Authors
-------
This module is based on the FluxReader by Zarko Pavlovic. It was generalized
to handle other analysis tree inputs by Adam Lister and Andy Mastbaum.

