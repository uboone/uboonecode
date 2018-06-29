#ifndef CRT_MERGER_HH
#define CRT_MERGER_HH

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include <artdaq-core/Data/Fragment.hh>
#include "fhiclcpp/ParameterSet.h"
#include "canvas/Utilities/InputTag.h"
//#include "gallery/Event.h"
//#include "gallery/ValidHandle.h"
#include <string>
#include <istream>
#include <set>
//#include "ifdh.h"
#include "IFDH_service.h"

using namespace art;
using namespace std;

namespace crt
{
  class CRTMerger : public art::EDProducer
  {
    std::vector<std::string> fFileNames;
    std::vector<std::string> fTestFiles;
    std::string previouscrtrootfile;
		 
    // Producer tag of the CRT events
    art::InputTag fTag;
    std::string fUBversion_CRTHits;
    std::string Merged_Object;
		
    // Time window
    unsigned long fTimeWindow;
    unsigned int fMaxCount;
    std::vector<unsigned long> fTimeOffSet;
    bool _debug;
		
  public:
    CRTMerger(const fhicl::ParameterSet&);
    ~CRTMerger();
    //explicit CRTMerger(fhicl::ParameterSet const &p);
		
    void produce( art::Event &evt ) override;
		
    void reconfigure(fhicl::ParameterSet const & p) override;
  private:
		
    std::vector< std::vector< artdaq::Fragment > > w;
    ifdh_ns::ifdh* tIFDH=0;
    ifdh_ns::ifdh* fIFDH=0;
		
    std::string data_label_DAQHeader_;
    std::string cTag;

    std::set<std::string> fCRTFiles;
		
  };
}
#endif // CRT_MERGER_HH
