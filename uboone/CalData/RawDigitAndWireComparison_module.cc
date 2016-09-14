////////////////////////////////////////////////////////////////////////
// Class:       RawDigitAndWireComparison
// Module Type: analyzer
// File:        RawDigitAndWireComparison_module.cc
//
// Generated at Sat Jan  3 21:14:09 2015 by Wesley Ketchum using artmod
// from cetpkgsupport v1_07_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include <numeric>
#include "RawDigitAndWireComparisonAlg.h"

namespace caldata{
  class RawDigitAndWireComparison;
}

namespace caldata{
  
  class RawDigitAndWireComparison : public art::EDAnalyzer {
  public:
    explicit RawDigitAndWireComparison(fhicl::ParameterSet const & p);
    // The destructor generated by the compiler is fine for classes
    // without bare pointers or other resource use.
    
    // Plugins should not be copied or assigned.
    RawDigitAndWireComparison(RawDigitAndWireComparison const &) = delete;
    RawDigitAndWireComparison(RawDigitAndWireComparison &&) = delete;
    RawDigitAndWireComparison & operator = (RawDigitAndWireComparison const &) = delete;
    RawDigitAndWireComparison & operator = (RawDigitAndWireComparison &&) = delete;
    
    void analyze(art::Event const & e) override;
    void beginJob();
    
  private:

    std::string                   fRawDigitModuleLabel;
    std::string                   fWireModuleLabel;
    RawDigitAndWireComparisonAlg  fAlg;

    TTree *fROICompareTree;
  };
  
  
  RawDigitAndWireComparison::RawDigitAndWireComparison(fhicl::ParameterSet const & p)
    : EDAnalyzer(p)  ,
      fRawDigitModuleLabel(p.get<std::string>("RawDigitModuleLabel")),
      fWireModuleLabel(p.get<std::string>("WireModuleLabel")),
      fAlg(p.get<fhicl::ParameterSet>("RawDigitAndWireComparisonAlg"))
  {}
  
  void RawDigitAndWireComparison::beginJob()
  {
    art::ServiceHandle<art::TFileService> tfs;
    fROICompareTree = tfs->make<TTree>("ROICompareTree","ROICompareTree");
    fAlg.SetROIOutputTree(fROICompareTree);
  }

  void RawDigitAndWireComparison::analyze(art::Event const & e)
  {

    //get event and run numbers
    unsigned int eventNumber = e.id().event();
    unsigned int runNumber   = e.run();    //get the raw digits

    //get the raw data
    art::Handle< std::vector<raw::RawDigit> > digitHandle;
    e.getByLabel(fRawDigitModuleLabel,digitHandle);
    std::vector<raw::RawDigit> const& digitVector(*digitHandle);
    
    //get the wire data
    art::Handle< std::vector<recob::Wire> > wireHandle;
    e.getByLabel(fWireModuleLabel,wireHandle);
    std::vector<recob::Wire> const& wireVector(*wireHandle);

    //we need to create a wire-->digit association vector
    std::vector<unsigned int> assocVector(wireVector.size());
    for(size_t i_wire=0; i_wire<wireVector.size(); i_wire++)
      assocVector[i_wire] = wireVector[i_wire].Channel();
    
    
    fAlg.RunROICompare(wireVector,digitVector,assocVector,runNumber,eventNumber);
    
  }
  
}

DEFINE_ART_MODULE(caldata::RawDigitAndWireComparison)
  
