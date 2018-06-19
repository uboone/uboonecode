#include "uboone/CRT/CRTMerger.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
////#include "uboone/CRT/CRTProducts/MSetCRTFrag.hh"
#include <artdaq-core/Data/Fragment.hh>
#include "gallery/Event.h"
#include "CRTBernFEBDAQCore/Overlays/BernZMQFragment.hh"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "art/Framework/Principal/Handle.h"
#include "IFDH_service.h"
#include "uboone/RawData/utils/DAQHeaderTimeUBooNE.h"
#include "uboone/CRT/CRTProducts/CRTHit.hh"

#include <memory>
#include <string>
#include <vector>
#include <exception>
#include <sstream>
#include <unistd.h>
#include <ctime>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/local_time_adjustor.hpp"
#include <boost/date_time/c_local_time_adjustor.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

namespace {

  // Local function to calculate absolute difference between two unsigned longs.

  unsigned int absdiff(unsigned long a, unsigned long b) {
    return (a>b ? a-b : b-a);
  }
}

using namespace boost::posix_time;
crt::CRTMerger::CRTMerger(const fhicl::ParameterSet& pset): data_label_DAQHeader_(pset.get<std::string>("data_label_DAQHeader_"))
{
	std::cout<<"crt::CRTMerger::CRTMerger"<<std::endl;
	
	setenv("TZ", "CST+6CDT", 1);
	tzset();
	
	//setenv("IFDH_DATA_DIR","/uboone/data/users/kolahalb/MicroBooNE/",1);
	//std::cout<<"ifdh_data_dir "<<getenv("IFDH_DATA_DIR")<<std::endl;
	



	this->reconfigure(pset);
	produces< std::vector<crt::CRTHit> >();
	fTag		= pset.get<art::InputTag> ("InputTagName");
	_debug		= pset.get<bool>		 ("debug");
	fTimeOffSet	= pset.get<std::vector< unsigned long > > ("test_t_offset");

	previouscrtrootfile = "";
	if ( ! tIFDH ) tIFDH = new ifdh_ns::ifdh;
}

crt::CRTMerger::~CRTMerger()
{
	//fIFDH->cleanup();
	std::cout<<"crt::CRTMerger::~CRTmerger"<<std::endl;
	
	if ( ! tIFDH )
	  delete tIFDH;
}

void crt::CRTMerger::produce(art::Event& event)
{
	if (_debug)
	{
		std::cout <<"crt::CRTMerger::produce, NEW EVENT" << std::endl;
	}

	const int T0(5); // Number of files/seconds of CRT wrt TPC evt time less than which we won't check for merge candidates..
	
	//For this event
	if(_debug)
	  std::cout<<fMaxCount<<std::endl;
	
	//get DAQ Header
	art::Handle< raw::DAQHeaderTimeUBooNE > rawHandle_DAQHeader;
	event.getByLabel(data_label_DAQHeader_, rawHandle_DAQHeader);
	
	if(!rawHandle_DAQHeader.isValid())
	{
		std::cout << "Run " << event.run() << ", subrun " << event.subRun()<< ", event " << event.event() << " has zero"<< " DAQHeaderTimeUBooNE  " << " in with label " << data_label_DAQHeader_ << std::endl;
    		return;
	}
	
	raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
	art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
	art::Timestamp evtTimeNTP = my_DAQHeader.ntp_time();
	
	std::cout<<"evt_timeGPS_sec "<<evtTimeGPS.timeHigh()<<",  evt_timeGPS_nsec "<<evtTimeGPS.timeLow()<<",  evt_timeNTP_sec "<<evtTimeNTP.timeHigh()<<",  evt_timeNTP_nsec "<<evtTimeNTP.timeLow()<<std::endl;
	
	// First find the art event time stamp
	//art::Timestamp evtTime = event.time();
	art::Timestamp evtTime = evtTimeGPS;
	
	unsigned long evt_time_sec = evtTime.timeHigh();//+fTimeOffSet[2];	
	unsigned long evt_time_nsec = evtTime.timeLow();
	
	// Use the information for configuring SAM query about the coincident crt-binrary-raw files
	unsigned long time_tpc = evt_time_sec*1000000000 + evt_time_nsec;
	unsigned long time_tpc1= time_tpc/1000;
	
	const char* tz = getenv("TZ");
	std::string tzs(tz);
	if (_debug)
	  std::cout<<"time-zone "<<tzs<<std::endl;
	
	if(tz != 0 && *tz != 0)
	{
		std::string tzs(tz);
		if(tzs != std::string("CST+6CDT"))
		{
			// Timezone is wrong, throw exception.
			throw cet::exception("CRTMerger") << "Wrong timezone: " << tzs;
		}
	}
	else
	{
		// Timezone is not set.  Throw exception.
		throw cet::exception("CRTMerger") << "Timezone not set.";
	}
	
	boost::posix_time::ptime const time_epoch(boost::gregorian::date(1970, 1, 1));
	boost::posix_time::ptime this_event_time = time_epoch + boost::posix_time::microseconds(time_tpc1);
	typedef boost::date_time::c_local_adjustor<ptime> local_adj;
	ptime this_event_localtime = local_adj::utc_to_local(this_event_time);
	std::cout<<"Local time of event: "<< boost::posix_time::to_iso_extended_string(this_event_localtime)<<std::endl;
	

	std::cout << boost::posix_time::to_iso_extended_string(this_event_time)<<std::endl;
	
	std::string stringTime = boost::posix_time::to_iso_extended_string(this_event_localtime);
	stringTime = "'"+stringTime+"'";
	
	struct tm tm;
	memset(&tm, 0, sizeof(tm));
	tm.tm_isdst = -1;
	
	if (_debug)
	  std::cout<<"TPC event time (s) "<<evt_time_sec<<"   event time (ns) "<<evt_time_nsec<<std::endl;
	
	// At this stage, find out the crtdaq files that may be time co-incident
	std::ostringstream dim;
	dim<<"file_format "<<"crt-binaryraw"<<" and file_type "<<"data"<<" and start_time < "<<stringTime<<" and end_time > "<<stringTime;
	
	if (_debug)
	  std::cout<<"dim = "<<dim.str()<<std::endl;
	
	// List those crtdaq files:
	std::vector< std::string > crtfiles = tIFDH->translateConstraints(dim.str());
	
	// Find the corresponding child artroot file
	std::vector< std::string > crtrootfile;
	std::vector< std::string > tmprootfile;


	if ( fUBversion_CRTHits == NULL ) {
	  std::cout << "Did not retrieve value for UBOONECODE_VERSION nor is there a specified CRTHits version to use. Will not find any proper CRT artroot daughters to merge." << std::endl;
	}

	std::string ubversion(fUBversion_CRTHits);

	for(unsigned k =0; k<crtfiles.size(); k++)
	{
		std::ostringstream dim1;
		// add constraint that current CRTMerge job release must match that in which the CRTHits were created/swizzled.
		dim1<<"file_format "<<"artroot"<<" and ub_project.version " << ubversion << " and ischildof: (file_name "<<crtfiles[k]<<" with availability physical )"<<std::endl;
		
		if (_debug)
		  std::cout<<"dim1 = "<<dim1.str()<<std::endl;
		
		tmprootfile = tIFDH->translateConstraints(dim1.str());
		
		if (tmprootfile.size()>0) {
		  std::cout << "We get " << tmprootfile.size()  << " daughters of " << crtfiles[k] << ". Pushing them/it onto vector of artroot files "  << std::endl;
		  for (const auto& artrootchild : tmprootfile)
		    crtrootfile.push_back(artrootchild);
		}
	}
	std::cout<<"total: "<<crtrootfile.size()<<std::endl;
	if (!crtrootfile.size())
	  std::cout << "\n\t CRTMerger_module: No child CRT files found that conform to constraints: " << "file_format "<<"artroot"<<" and ub_project.version " << ubversion  << std::endl;

	// Throw exception if there are fewer than six CRT files.

	if(crtrootfile.size() < 6) {
	  throw cet::exception("CRTMerger") << "Too few matching CRT files: " 
					    << crtrootfile.size() << "\n";
	}
	  
	
	std::unique_ptr<std::vector<crt::CRTHit> > CRTHitEventsSet(new std::vector<crt::CRTHit>); //collection of CRTHits for this event
	
	for(unsigned crf_index = 0; crf_index < crtrootfile.size(); crf_index++)
	{
	  if (_debug)
	    std::cout<<"The child artroot file is "<<crtrootfile[crf_index]<<std::endl;

	  // Add this file to set of seen CRT files.

	  if(fCRTFiles.count(crtrootfile[crf_index]) == 0) {
	    std::cout << "Adding CRT parent file: " << crtrootfile[crf_index] << std::endl;
	    art::ServiceHandle<art::FileCatalogMetadata> md;
	    std::ostringstream ostr;
	    ostr << "mixparent" << fCRTFiles.size();
	    md->addMetadataString(ostr.str(), crtrootfile[crf_index]);
	    fCRTFiles.insert(crtrootfile[crf_index]);
	  }
	
	  // Read off the root file by streaming over internet. Use xrootd URL
	  // This is alternative to use gsiftp URL. In that approach we use the URL to ifdh::fetchInput the crt file to local directory and keep it open as long as we make use of it
	  //xrootd URL
	  std::string schema = "root";
	  std::vector< std::string > crtrootFile_xrootd_url;
	  try{
	    std::vector<std::string> tmp(tIFDH->locateFile(crtrootfile[crf_index], schema));
	    crtrootFile_xrootd_url.swap(tmp);
	  }
	  catch(...)
	    {
	      std::cout << "This Root File does not exist. No Merger CRTHit candidates put on event for this TPC evt for this chunk of the CRT." << std::endl;
	      throw cet::exception("CRTMerger") << "Could not locate CRT file: " 
						<< crtrootfile[crf_index] << "\n";
	    }
	  std::cout<<"xrootd URL: "<<crtrootFile_xrootd_url[0]<<std::endl;
	
	  // gallery, when fed the list of the xrootd URL, internally calls TFile::Open() to open the file-list
	  // In interactive mode, you have to get your proxy authenticated by issuing:
	  // voms-proxy-init -noregen -rfc -voms fermilab:/fermilab/uboone/Role=Analysis
	  // when you would like to launch a 'lar -c ... ... ...'
	  // In batch mode, this step is automatically done
	  
	  gallery::Event fCRTEvent(crtrootFile_xrootd_url);
	  std::cout<<"Opened the CRT root file from xrootd URL"<<std::endl;
	
	  if (_debug)
	    std::cout<<"Start merging attempts"<<std::endl;
	  
	  unsigned int merging = 0;
	  //int triscuit = 0;
	  if (_debug)
	    std::cout << "New [collection of] CRTEvent[s]. Its size is  " << fCRTEvent.numberOfEventsInFile() << "." << std::endl;

	  unsigned long jump = 1; 
	  bool firstE(true);
	  unsigned long cnt(0);
	  int cnt2(0);
	  for(fCRTEvent.toBegin(); !fCRTEvent.atEnd(); ++fCRTEvent)
	    {
	      // would rather be allowed to advance the iterator by more than 1, but gallery doesn't seem to allow it. So must do the below.

	      if (jump!=1 && cnt!=jump)
		{
		  cnt++;
		  continue;
		}
	      if (cnt==jump && _debug) 
		std::cout << "Okay we've jumped fwd " << cnt << " CRT sub events. " << std::endl;
	      //art::Handle< std::vector< crt::CRTHit> > rawHandle;
	      //fCRTEvent.getByLabel(cTag, rawHandle);

		
	      std::vector< crt::CRTHit >  const& CRTHitCollection = *(fCRTEvent.getValidHandle< std::vector<crt::CRTHit> >(cTag));
	      jump = 1;
	      // If our TPC evt is more than T0 sec beyond the end of the first CRT event's last Hit, then ffwd
	      // I am assuming these CRT events are precisely 1 second long each. 
	      if ( ( evt_time_sec > (CRTHitCollection[CRTHitCollection.size()-1].ts0_s + T0) ) && firstE )
		{
		  firstE = false;

		  // Let's stop T0 sec short of where we'd like to land.
		  jump = evt_time_sec - CRTHitCollection[CRTHitCollection.size()-1].ts0_s - T0; 
		  // if jump is a large number, don't zip off end of list. Stop short by 2 files.
		  if (jump > (unsigned long)(fCRTEvent.numberOfEventsInFile()-2) && ((fCRTEvent.numberOfEventsInFile()-2)>0) ) jump = fCRTEvent.numberOfEventsInFile()-2;
		  // Sometimes the Top panel has tiny (0,1) values, so make sure seconds are not crazy since the time since Epoch.
		  if (CRTHitCollection[0].ts0_s<1300000000 && CRTHitCollection[CRTHitCollection.size()-1].ts0_s<1300000000) 
		    {
		      // don't ffwd any sub-events, but don't rule out doing it later (meaning, keep firstE set to true)
		      jump = 1;
		      firstE = true;
		      cnt2++;
		      if( _debug)
			std::cout << "Count of hits w bad times in this file is " << cnt2 << std::endl;
		      continue;
		    }
		  if (_debug)
		    std::cout << "Within this CRT sub event the first and last CRT Hit times in seconds are " << CRTHitCollection[0].ts0_s << " and " << CRTHitCollection[CRTHitCollection.size()-1].ts0_s << " whereas TPC event sec is " << evt_time_sec << " ... and therefore we are going to jump fwd "<< jump << " sub events." <<  std::endl;
		  continue;

		}
	      if ( evt_time_sec < (CRTHitCollection[0].ts0_s-T0))
		{
		  std::cout << "Within this CRT sub event the last CRT Hit times in seconds is "  << CRTHitCollection[CRTHitCollection.size()-1].ts0_s << " whereas TPC event sec is " << evt_time_sec << " ... and therefore we are moving to next full collection of CRT Events (new file)." <<  std::endl;

		  break;
		}



	      if (_debug) std::cout << "CRT Event second falls inside CRT sub-event Window ..." << std::endl;

	      

	      bool exitCollection(false);
	      for(std::vector<int>::size_type hit_index=0; hit_index!=CRTHitCollection.size(); hit_index++)
		{
		  crt::CRTHit CRTHitevent = CRTHitCollection[hit_index];
		  unsigned long CRTtime_s	= CRTHitevent.ts0_s;
		  unsigned long CRTtime_ns= CRTHitevent.ts0_ns;
		  
		  //unsigned long TPCtime_ns= evt_time_sec*1000000000+evt_time_nsec;
		  unsigned long TPCtime_s	= evt_time_sec;
		  unsigned long TPCtime_ns= evt_time_nsec;
		  
		  unsigned long MergingWindow_start = TPCtime_ns - 2000000;
		  unsigned long MergingWindow_end	  = TPCtime_ns + 4000000;
		  

		  if (_debug && 0) // TMI, even in debug mode
		    {
		      std::cout<<"TPC_ns: "<<TPCtime_s<<", CRT_ns: "<<CRTtime_s<<std::endl;
		    }

		  if(absdiff(CRTtime_s, TPCtime_s)<1 ) // was <2. Change at  DL's request. EC, 12-Apr-2018.
		    {
		      if ((CRTtime_ns > MergingWindow_start) && (CRTtime_ns < MergingWindow_end))
			{
			  if (_debug)
			    std::cout<<"found match"<<std::endl;
			  CRTHitEventsSet->emplace_back(CRTHitevent);
			  merging += 1;

			}
		    }
		} // end loop on Hit Collection within CRT evt
	      if (exitCollection) break;
	      
	    } // end loop on CRT evts
	  std::cout<<"# merging in the stream: "<<merging<<std::endl;
	} // end loop on ROOT File
	std::cout<<"# of Merged CRTHits in CRTHitEventsSet being written to event:: "<<CRTHitEventsSet->size()<<std::endl;
	event.put(std::move(CRTHitEventsSet));
	
	if (_debug)
	  std::cout<<"---X---"<<std::endl;
}


void crt::CRTMerger::reconfigure(fhicl::ParameterSet const & pset)
{
	std::cout<<"crt::CRTMerger::reconfigure"<<std::endl;
	char const* ubchar = std::getenv( "UBOONECODE_VERSION" );
	
	cTag = {pset.get<std::string>("data_label_CRTHit_")};
	fTag = {pset.get<std::string>("InputTagName","crthit")};
	fTimeWindow = pset.get<unsigned>("TimeWindow",5000000);
	fUBversion_CRTHits   = pset.get<std::string>   ("ubversion_CRTHits",ubchar);
}
DEFINE_ART_MODULE(crt::CRTMerger)
