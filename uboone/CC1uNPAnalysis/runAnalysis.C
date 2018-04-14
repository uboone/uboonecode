#include "TROOT.h"
#include "TChain.h"
#include "TSystem.h"

void runAnalysis () {
  
  TChain chain("cc1unpselana/fMC_TrunMean");
  //chain.Add("/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/MCC8.7_tune1/out/cc1unp_MC_tune1.root");  //MC Tune1
  //chain.Add("/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/BNB_mcc8.8/out/cc1unp_onbeam_mcc8.8.root"); //On-Beam Data
  //chain.Add("/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/EXTBNB_mcc8.8/out/cc1unp_offbeam_mcc8.8.root"); // Off-beam data
  chain.Add("/pnfs/uboone/persistent/users/jiangl/CC1uNPSelection/MCC8.7_tune3/out/cc1unp_MC_tune3.root"); // Tune3

  std::cout<<"start processing the analysis"<<std::endl;


  chain.Process("hanalysis.C");
}
  
