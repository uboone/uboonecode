#ifndef CRT_MERGER_EXTRA_HH
#define CRT_MERGER_EXTRA_HH

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
#include "gallery/Event.h"
//#include "gallery/ValidHandle.h"
#include <string>
#include <istream>
#include <set>
#include <map>
//#include "ifdh.h"
#include "IFDH_service.h"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "uboone/CRT/CRTProducts/CRTHit.hh"

using namespace art;
using namespace std;

namespace crt
{
  class CRTMergerExtra : public art::EDProducer
  {
    std::vector<std::string> fFileNames;
    std::vector<std::string> fTestFiles;
    std::string previouscrtrootfile;
		 
    // Producer tag of the CRT events
    std::string fUBversion_CRTHits;
    std::string Merged_Object;
		
    // Time window
    unsigned long fTimeWindow;
    unsigned int fMaxCount;
    std::vector<unsigned long> fTimeOffSet;
    bool _debug;
		
  public:

    // Nested struct with information about CRT binary files.

    struct CRTFileInfo
    {
      std::string fFileName;                 // Name of CRT binary file.
      boost::posix_time::ptime fStartTime;   // Start time of CRT binary file.
      boost::posix_time::ptime fEndTime;     // End time of CRT binary file.
      std::vector<std::string> fSwizzled;    // Names of CRT swizzled files.
    };

    CRTMergerExtra(const fhicl::ParameterSet&);
    ~CRTMergerExtra();
    //explicit CRTMerger(fhicl::ParameterSet const &p);
		
    void produce( art::Event &evt ) override;
		
    void reconfigure(fhicl::ParameterSet const & p) override;

    // Find CRT swizzled files that matches a particular event time.

    std::vector<std::string> findMatchingCRTFiles(boost::posix_time::ptime event_time);

    // Reposiiton gallery event to tpc event time.

    bool reposition(gallery::Event& event, unsigned long evt_time_sec, int delta_t);

    // Filter CRT hit collection.

    void filter_crt_hits(const std::vector<crt::CRTHit>& input_hits, 
			 std::vector<crt::CRTHit>& output_hits) const;

  private:
		
    std::vector< std::vector< artdaq::Fragment > > w;
    ifdh_ns::ifdh* tIFDH=0;
    ifdh_ns::ifdh* fIFDH=0;
		
    std::string data_label_DAQHeader_;
    std::string cTag;

    // Files that we know about.

    std::map<std::string, CRTFileInfo> fCRTFiles;
    std::set<std::string> fCRTSwizzledFiles;		

    // Open CRT swizzled gallery files.
    // After a CRT swizzled file is open, we never close it.
    // In normal circumstances, we should never have more than 18 
    // CRT swizzled files open (one seven-hour TPC run can overlap
    // with up to three four-hour CRT runs, and there are six CRT
    // swizzled streams).

    std::map<std::string, std::shared_ptr<gallery::Event> > fCRTEvents;
  };
}
#endif // CRT_MERGER_EXTRA_HH
