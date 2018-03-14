#ifndef __uboone_NTupleInterface__
#define __uboone_NTupleInterface__

/**
 * Interface for importing Ntuples.
 *
 * \author A. Lister, A. Mastbaum
 */

#include "fhiclcpp/ParameterSet.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "TreeInterface.h"

class TTree;

namespace uboone {

/**
 * \class NTupleInterface
 * \brief Wrapper to access user TTrees
 */
class NTupleInterface : public TreeInterface {
public:
  /** Constructor. */
  NTupleInterface();

  /** Destructor. */
  ~NTupleInterface();

  /** Get the number of entries in the input tree. */
  const Long64_t GetEntries() { return fNEntries; }

  /** Get the current run number. */
  const int GetRun() { return fRun; }

  /** Get the POT count. */
  const float GetPOT() { return fPOT; }

  /** Get the neutrino position for the current entry. */
  const TLorentzVector GetNuPosition() { return fNuPos; }

  /** Get the neutrino momentum for the current entry. */
  const TLorentzVector GetNuMomentum() { return fNuMom; }

  /**
   * Set up the input ROOT file.
   *
   * The branch names are FHICL-configurable, with the mapping specified
   * in branchDef. These names are used to set the branch pointers for the
   * data required to populate the truth-level objects.
   *
   * \param rootFile The input ROOT file
   * \param treeName Path to the input TTree object
   * \param branchDef A FHICL config specifying the branch name mapping
   */
  void SetRootFile(TFile* rootFile,
                   TString treeName,
                   fhicl::ParameterSet& branchDef);

  /** 
   * Fill in the MCFlux object from the tree data.
   *
   * \param ientry Tree entry index
   * \param mcflux MCFlux object to fill in (by reference)
   * \returns True if successful
   */
  bool FillMCFlux(Long64_t ientry, simb::MCFlux& mcflux);

private:
  TTree* fTree;  //!< Input TTree
  double vtxx, vtxy, vtxz;  //!< Neutrino vertex
  double px, py, pz, E;  //!< Neutrino momentum
  int pdg;  //!< Neutrino PDG code
  double wgt;  //!< Neutrino weight (i.e. from upstream flux simulation)
  double dist;  //!< ???
  int ptype;  //!< Parent PDG?
  int run;  //!< Run ID
  double nenergyn;  //!< ???
  Long64_t fNEntries;  //!< Number of entries in the tree
  int fRun;  //!< Run ID
  TLorentzVector fNuPos;  //!< Neutrino vertex vector
  TLorentzVector fNuMom;  //!< Neutrino momentum vector
};

}  // namespace uboone

#endif // __uboone_NTupleInterface__

