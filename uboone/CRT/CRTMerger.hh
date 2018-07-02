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
#include <map>
//#include "ifdh.h"
#include "IFDH_service.h"
#include "boost/date_time/posix_time/posix_time.hpp"

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

    // Nested struct with information about CRT files.

    struct CRTFileInfo
    {
      std::string fFileName;                 // Name of CRT binary file.
      boost::posix_time::ptime fStartTime;   // Start time of CRT binary file.
      boost::posix_time::ptime fEndTime;     // End time of CRT binary file.
      std::vector<std::string> fSwizzled;    // Names of CRT swizzled files.
    };
    CRTMerger(const fhicl::ParameterSet&);
    ~CRTMerger();
    //explicit CRTMerger(fhicl::ParameterSet const &p);
		
    void produce( art::Event &evt ) override;
		
    void reconfigure(fhicl::ParameterSet const & p) override;

    // Find CRT swizzled files that matches a particular event time.

    std::vector<std::string> findMatchingCRTFiles(boost::posix_time::ptime event_time);

  private:
		
    std::vector< std::vector< artdaq::Fragment > > w;
    ifdh_ns::ifdh* tIFDH=0;
    ifdh_ns::ifdh* fIFDH=0;
		
    std::string data_label_DAQHeader_;
    std::string cTag;

    // Files that we know about.

    std::map<std::string, CRTFileInfo> fCRTFiles;
    std::set<std::string> fCRTSwizzledFiles;		
  };
}
#endif // CRT_MERGER_HH
