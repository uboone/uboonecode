#include <algorithm>
#include "uboone/CRT/CRTMergerExtra.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"
////#include "uboone/CRT/CRTProducts/MSetCRTFrag.hh"
#include <artdaq-core/Data/Fragment.hh>
#include "CRTBernFEBDAQCore/Overlays/BernZMQFragment.hh"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "art/Framework/Principal/Handle.h"
#include "IFDH_service.h"
#include "uboone/RawData/utils/DAQHeaderTimeUBooNE.h"

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

#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/local_time_adjustor.hpp"
#include <boost/date_time/c_local_time_adjustor.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

namespace {

  // Local function to compare CRT hits according to time.

  bool CRTHitComp(const crt::CRTHit& crthit1, const crt::CRTHit& crthit2)
  {
    bool result = (crthit1.ts0_s < crthit2.ts0_s ||
		   (crthit1.ts0_s == crthit2.ts0_s && crthit1.ts0_ns < crthit2.ts0_ns));
    return result;
  }
}

crt::CRTMergerExtra::CRTMergerExtra(const fhicl::ParameterSet& pset) :
  data_label_DAQHeader_(pset.get<std::string>("data_label_DAQHeader_"))
{
  std::cout<<"crt::CRTMergerExtra::CRTMergerExtra"<<std::endl;

  setenv("TZ", "CST+6CDT", 1);  // Fermilab time zone.
  tzset();

  this->reconfigure(pset);
  produces< std::vector<crt::CRTHit> >();
  _debug = pset.get<bool>("debug");
  fTimeOffSet = pset.get<std::vector< unsigned long > > ("test_t_offset");

  if ( ! tIFDH ) tIFDH = new ifdh_ns::ifdh;
}

crt::CRTMergerExtra::~CRTMergerExtra()
{
  //fIFDH->cleanup();
  std::cout<<"crt::CRTMergerExtra::~CRTMergerExtra"<<std::endl;
  
  if ( ! tIFDH )
    delete tIFDH;
}

void crt::CRTMergerExtra::produce(art::Event& event)
{
  if (_debug) {
    std::cout <<"crt::CRTMergerExtra::produce, NEW EVENT" << std::endl;
  }

  const int T0(5); // Number of files/seconds of CRT wrt TPC evt time
                   //less than which we won't check for merge candidates..
  
  //For this event
  if(_debug)
    std::cout<<fMaxCount<<std::endl;
  
  //get DAQ Header
  art::Handle< raw::DAQHeaderTimeUBooNE > rawHandle_DAQHeader;
  event.getByLabel(data_label_DAQHeader_, rawHandle_DAQHeader);
  
  if(!rawHandle_DAQHeader.isValid()) {
    std::cout << "Run " << event.run() 
	      << ", subrun " << event.subRun()
	      << ", event " << event.event()
	      << " has zero DAQHeaderTimeUBooNE with label "
	      << data_label_DAQHeader_ << std::endl;
    return;
  }
  
  raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
  art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
  art::Timestamp evtTimeNTP = my_DAQHeader.ntp_time();
  
  std::cout<<"evt_timeGPS_sec " << evtTimeGPS.timeHigh()
	   <<",  evt_timeGPS_nsec " << evtTimeGPS.timeLow()
	   <<",  evt_timeNTP_sec " << evtTimeNTP.timeHigh()
	   <<",  evt_timeNTP_nsec "<< evtTimeNTP.timeLow() << std::endl;

  // Quit if GPS time is zero.
  // This is somehow normal for the first event of a run...
  if(evtTimeGPS.timeHigh() == 0 && evtTimeGPS.timeLow() == 0) {

    // We promised to put a collection of CRTHits into the event.

    std::unique_ptr<std::vector<crt::CRTHit> > p(new std::vector<crt::CRTHit>);
    event.put(std::move(p));
    return;
  }
  
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
  
  if(tz != 0 && *tz != 0) {
    std::string tzs(tz);
    if(tzs != std::string("CST+6CDT")) {
      // Timezone is wrong, throw exception.
      throw cet::exception("CRTMergerExtra") << "Wrong timezone: " << tzs;
    }
  }
  else {
    // Timezone is not set.  Throw exception.
    throw cet::exception("CRTMergerExtra") << "Timezone not set.";
  }
  
  boost::posix_time::ptime const time_epoch(boost::gregorian::date(1970, 1, 1));
  boost::posix_time::ptime this_event_time = 
    time_epoch + boost::posix_time::microseconds(time_tpc1);
  typedef boost::date_time::c_local_adjustor<boost::posix_time::ptime> local_adj;
  boost::posix_time::ptime this_event_localtime = local_adj::utc_to_local(this_event_time);
  std::cout << "Local time of event: " 
	    << boost::posix_time::to_iso_extended_string(this_event_localtime)<<std::endl;

  std::cout << "UTC time of event:   "
	    << boost::posix_time::to_iso_extended_string(this_event_time)<<std::endl;

  // Get matching swizzled CRT files.

  std::vector<std::string> crtrootfiles = findMatchingCRTFiles(this_event_localtime);
  
  // collection of CRTHits for this event

  std::unique_ptr<std::vector<crt::CRTHit> > CRTHitEventsSet(new std::vector<crt::CRTHit>);

  // Loop over swizzled CRT files.
  
  for(auto const& crtrootfile : crtrootfiles) {
    std::cout << "\nThe child artroot file is " << crtrootfile << std::endl;

    // Add this file to set of seen CRT files.

    if(fCRTSwizzledFiles.count(crtrootfile) == 0) {
      std::cout << "Adding CRT parent file: " << crtrootfile << std::endl;
      art::ServiceHandle<art::FileCatalogMetadata> md;
      std::ostringstream ostr;
      ostr << "mixparent" << fCRTSwizzledFiles.size();
      md->addMetadataString(ostr.str(), crtrootfile);
      fCRTSwizzledFiles.insert(crtrootfile);
    }
  
    // Read off the root file by streaming over internet. Use xrootd URL
    // This is alternative to use gsiftp URL. In that approach we use the URL to
    // ifdh::fetchInput the crt file to local directory and keep it open as long
    // as we make use of it
    //xrootd URL

    std::shared_ptr<gallery::Event> crt_event;

    // Do we have an open gallery event for this file?

    if(fCRTEvents.count(crtrootfile) > 0) {

      // Yes, reuse the open gallery event.

      std::cout << "Reusing open gallery event for file " << crtrootfile << std::endl;
      crt_event = fCRTEvents[crtrootfile];
    }

    else {

      // No, open a new gallery event.

      std::string schema = "root";
      std::vector< std::string > crtrootfiles_xrootd_url;
      try {
	std::vector<std::string> tmp(tIFDH->locateFile(crtrootfile, schema));
	crtrootfiles_xrootd_url.swap(tmp);
      }
      catch(...) {
	std::cout << "This Root File does not exist. No Merger CRTHit candidates put on event for this TPC evt for this chunk of the CRT." << std::endl;
	throw cet::exception("CRTMergerExtra") << "Could not locate CRT file: " 
					  << crtrootfile << "\n";
      }

      // Ignore all locations except the first.

      if(crtrootfiles_xrootd_url.size() > 1)
	crtrootfiles_xrootd_url.erase(crtrootfiles_xrootd_url.begin()+1,
				      crtrootfiles_xrootd_url.end());
      std::cout<<"xrootd URL: "<<crtrootfiles_xrootd_url[0]<<std::endl;
  
      // gallery, when fed the list of the xrootd URL, internally calls TFile::Open()
      // to open the file-list
      // In interactive mode, you have to get your proxy authenticated by issuing:
      // voms-proxy-init -noregen -rfc -voms fermilab:/fermilab/uboone/Role=Analysis
      // when you would like to launch a 'lar -c ... ... ...'
      // In batch mode, this step is automatically done
    
      // Make several attempts to open file.

      int mtry = 5;
      int wait = 0;
      bool open_ok = false;

      while(mtry-- > 0 && !open_ok) {
	if(wait != 0) {
	  std::cout << "Waiting " << wait << " seconds." << std::endl;
	  sleep(wait);
	  wait *= 2;
	}
	else
	  wait = 1;

	try {
	  std::cout << "Open file." << std::endl;
	  crt_event = std::shared_ptr<gallery::Event>(new gallery::Event(crtrootfiles_xrootd_url));
	  open_ok = true;
	}
	catch(...) {
	  open_ok = false;
	}

	if(open_ok)
	  std::cout << "Open succeeded." << std::endl;
	else {
	  std::cout << "Open failed." << std::endl;
	}
      }
      std::cout<<"Opened the CRT root file from xrootd URL"<<std::endl;
      fCRTEvents[crtrootfile] = crt_event;
  
      std::cout << "New [collection of] CRTEvent[s]. Its size is  "
		<< crt_event->numberOfEventsInFile() << "." << std::endl;
    }

    // Attempt to reposition the CRT event.
    // In case of initial failure, rewind the file and retry once.

    std::cout << "CRT event entry before reposition = " << crt_event->eventEntry() << std::endl;
    bool ok = reposition(*crt_event, evt_time_sec, T0);

    // If reposition failed, skip this file.a

    if(!ok) {
      std::cout << "Failed to reposition file." << std::endl;
      continue;
    }

    std::cout << "CRT event entry after reposition = " << crt_event->eventEntry() << std::endl;

    int merging = 0;
    long double TPCtime = evt_time_sec + 1.e-9L * evt_time_nsec;
    long double MergingWindow_start = TPCtime - 0.002L;
    long double MergingWindow_end   = TPCtime + 0.004L;

    // Loop over CRT events.

    bool done = false;
    for(; !crt_event->atEnd(); ++*crt_event) {

      // Loop over events in the current CRT hit collection.

      gallery::ValidHandle<std::vector<crt::CRTHit> > h =
	crt_event->getValidHandle< std::vector<crt::CRTHit> >(cTag);
      std::vector<crt::CRTHit> CRTHitCollection;
      filter_crt_hits(*h, CRTHitCollection);
      if(_debug) {
	uint32_t first_sec = CRTHitCollection.front().ts0_s;
	uint32_t last_sec = CRTHitCollection.back().ts0_s;
	std::cout << "\nTPC event time = " << evt_time_sec << std::endl;
	std::cout << "First event time = " << first_sec << std::endl;
	std::cout << "Last event time = " << last_sec << std::endl;
      }
      for(auto const& CRTHitevent : CRTHitCollection) {
	long double CRTtime =  CRTHitevent.ts0_s - 1 + 1.e-9L * CRTHitevent.ts0_ns;//one second substracted for adding 1 extra second
	if ( CRTtime  >= MergingWindow_start && CRTtime <= MergingWindow_end) {
	  if (_debug)
	    std::cout<<"found match"<<std::endl;
	  CRTHitEventsSet->emplace_back(CRTHitevent);
	  merging += 1;
        }
	else if(CRTtime > MergingWindow_end) {
	  done = true;
	  break;
	}
      } // end loop on Hit Collection within CRT evt
      if(done)
	break;
    } // end loop on CRT evts
    std::cout<<"# merging in the stream: "<<merging<<std::endl;
  } // end loop on ROOT File
  std::cout << "# of Merged CRTHits in CRTHitEventsSet being written to event:: "
	    <<CRTHitEventsSet->size()<<std::endl;
  event.put(std::move(CRTHitEventsSet));
  
  if (_debug)
    std::cout<<"---X---"<<std::endl;
}


void crt::CRTMergerExtra::reconfigure(fhicl::ParameterSet const & pset)
{
  std::cout<<"crt::CRTMergerExtra::reconfigure"<<std::endl;
  cTag = {pset.get<std::string>("data_label_CRTHit_")};
  fTimeWindow = pset.get<unsigned>("TimeWindow",5000000);
  fUBversion_CRTHits   = pset.get<std::string>   ("ubversion_CRTHits");
}

// Find matching CRT files for a TPC event time.
//
// Notes:
//
// 1.  The argument should be the event time in the local (Fermilab) time zone.
// 2.  Start and end times in sam metadata are in the local time zone.

std::vector<std::string> crt::CRTMergerExtra::findMatchingCRTFiles(boost::posix_time::ptime event_time)
{

  // Result vector.
  std::vector< std::string > crtrootfiles;

  // Check cached files.

  boost::posix_time::time_duration one_minute(0, 1, 0, 0);
  for(const auto& crtfileinfo : fCRTFiles) {
    const std::string& crtfile = crtfileinfo.first;
    const CRTFileInfo& fileinfo = crtfileinfo.second;
    if(event_time >= fileinfo.fStartTime + one_minute && 
       event_time <= fileinfo.fEndTime - one_minute) {
      std::cout << "Found cached file " << crtfile << std::endl;
      for(auto const& crt_swizzled : fileinfo.fSwizzled)
	crtrootfiles.push_back(crt_swizzled);
    }
  }
  if(crtrootfiles.size() < 6) {

    // Didn't match enough cached files.

    crtrootfiles.erase(crtrootfiles.begin(), crtrootfiles.end());

    // Query CRT binery files.

    std::string stringTime = boost::posix_time::to_iso_extended_string(event_time);
    stringTime = "'"+stringTime+"'";
    std::ostringstream dim;
    dim << "file_format " << "crt-binaryraw"
	<<" and file_type " << "data"
	<<" and start_time < " << stringTime 
	<< " and end_time > " << stringTime;
  
    if (_debug)
      std::cout<<"dim = "<<dim.str()<<std::endl;
  
    // List those crtdaq files:
    std::vector< std::string > crtfiles = tIFDH->translateConstraints(dim.str());
    if(_debug)
      std::cout << "Found " << crtfiles.size() << " CRT binary files." << std::endl;

    // Get metadata of binary CRT files.

    for(auto const & crtfile : crtfiles) {
      if(_debug)
	std::cout << "\n" << crtfile << std::endl;
      std::string md = tIFDH->getMetadata(crtfile);

      // Extract CRT binary file start and end time from metadata as strings.

      std::string start_time;
      size_t n = md.find("Start Time: ");
      size_t m = 0;
      if(n < std::string::npos) {
	n += 12;
	m = md.find("+", n);
	if(m > n && m < std::string::npos) {
	  start_time = md.substr(n, m-n);
	  if(_debug)
	    std::cout << "Start time = " << start_time << std::endl;
	}
      }
      std::string end_time;
      n = md.find("End Time: ");
      if(n < std::string::npos) {
	n += 10;
	m = md.find("+", n);
	if(m > n && m < std::string::npos) {
	  end_time = md.substr(n, m-n);
	  if(_debug)
	    std::cout << "End time = " << end_time << std::endl;
	}
      }

      // Convert start and end times to boost ptime.

      boost::posix_time::ptime start_time_ptime =
	boost::date_time::parse_delimited_time<boost::posix_time::ptime>(start_time, 'T');
      boost::posix_time::ptime end_time_ptime =
	boost::date_time::parse_delimited_time<boost::posix_time::ptime>(end_time, 'T');

      // Make sure we closed the loop.

      std::string start_time2 = boost::posix_time::to_iso_extended_string(start_time_ptime);
      std::string end_time2 = boost::posix_time::to_iso_extended_string(end_time_ptime);
      if(_debug) {
	std::cout << "Start ptime " << start_time2 << std::endl;
	std::cout << "End ptime " << end_time2 << std::endl;
      }
      if(start_time2 != start_time || end_time2 != end_time) {
	std::cout << "Start ptime " << start_time2 << std::endl;
	std::cout << "End ptime " << end_time2 << std::endl;
	throw cet::exception("CRTMergerExtra") << "Problem converting start and end time.";
      }
  
      // Query CRT swizzled files that are children of this CRT binary file.

      std::ostringstream dim1;
      dim1 << "file_format " << "artroot"
	   <<" and ub_project.version " << fUBversion_CRTHits
	   << " and ischildof: (file_name " << crtfile
	   <<" with availability physical )";

      if (_debug)
	std::cout << "dim1 = " << dim1.str() << std::endl;
    
      std::vector< std::string > tmprootfiles = tIFDH->translateConstraints(dim1.str());
    
      std::cout << "Found " << tmprootfiles.size()
		<< " daughters of " << crtfile << std::endl;
      if (tmprootfiles.size()>0) {
	for (const auto& artrootchild : tmprootfiles) {
	  crtrootfiles.push_back(artrootchild);
	  if(_debug)
	    std::cout << artrootchild << std::endl;
	}
      }

      // Construct a CRTFileInfo struct and stash it in our cache, if not already there.

      if(fCRTFiles.count(crtfile) == 0) {
	std::cout << "Adding " << crtfile << " to file cache." << std::endl;
	CRTFileInfo fileinfo;
	fileinfo.fFileName = crtfile;
	fileinfo.fStartTime = start_time_ptime;
	fileinfo.fEndTime = end_time_ptime;
	fileinfo.fSwizzled = tmprootfiles;
	fCRTFiles[crtfile] = fileinfo;
	std::cout << "File cache contains " << fCRTFiles.size() << " files." << std::endl;
      }
    }
  }
  std::cout<<"\nNumber of matching CRT swizzled files: "<<crtrootfiles.size()<<std::endl;
  if(_debug) {
    for(const auto& crt_swizzled : crtrootfiles)
      std::cout << crt_swizzled << std::endl;
  }
  if (!crtrootfiles.size())
    std::cout << "\n\t CRTMergerExtra_module: No child CRT files found that conform to constraints: "
	      << "file_format " << "artroot" << " and ub_project.version "
	      << fUBversion_CRTHits << std::endl;

  // Throw exception if there are fewer than six CRT files.

  if(crtrootfiles.size() < 6) {
    throw cet::exception("CRTMergerExtra") << "Too few matching CRT files: " 
					   << crtrootfiles.size() << "\n";
  }

  // Done.

  return crtrootfiles;
}

// Reposiiton gallery event to approximate time in seconds.
// Return true if success, false if fail.

bool crt::CRTMergerExtra::reposition(gallery::Event& event, unsigned long evt_time_sec, int delta_t)
{
  // Result.

  bool ok = false;

  // Number of retries.

  int ntry = 2;
  while(ntry-- > 0 && !ok) {

    try {

      // Loop over events until we get a non-empty CRT his collection or an exception.

      for(;;event.next()) {

	// Look for a non-empty CRT hit collection.

	gallery::ValidHandle<std::vector<crt::CRTHit> > h =
	  event.getValidHandle< std::vector<crt::CRTHit> >(cTag);
	std::vector< crt::CRTHit > CRTHitCollection;
	filter_crt_hits(*h, CRTHitCollection);
	if(_debug)
	  std::cout << "Number of CRT hits in event = " << CRTHitCollection.size() << std::endl;
	if(CRTHitCollection.size() > 0) {

	  // Got a non-empty CRT hit collection.

	  uint32_t first_sec = CRTHitCollection.front().ts0_s;
	  uint32_t last_sec = CRTHitCollection.back().ts0_s;
	  if(_debug) {
	    std::cout << "TPC event time = " << evt_time_sec << std::endl;
	    std::cout << "First event time = " << first_sec << std::endl;
	    std::cout << "Last event time = " << last_sec << std::endl;
	  }

	  // Compare TPC and CRT times.

	  if(evt_time_sec > last_sec + delta_t) {

	    // TPC event is in the future relative to CRT event.
	    // Skip events assuming each CRT event is exactly one second.
	    // Return success.

	    int jump = evt_time_sec - last_sec - delta_t;
	    if(_debug)
	      std::cout << "Skipping ahead " << jump << " CRT events." << std::endl;
	    for(; jump>0; --jump)
	      event.next();
	    ok = true;
	    break;  // Break out of loop over events.
	  }
	  else if(evt_time_sec < first_sec) {

	    // TPC event is in the past relative to CRT event.
	    // Return failure, which may signal the calling program to rewind
	    // the file and try again.

	    ok = false;
	    break;  // Break out of loop over events.
	  }
	  else {

	    // TPC event is close to CRT event.
	    // Stay at the current position and return success.

	    ok = true;
	    break;  // Break out of loop over events.
	  }
	}
      }
    }
    catch(...) {
      ok = false;
    }

    if(!ok) {

      // Make several attempts to rewind/reopen file.

      int mtry = 5;
      int wait = 0;
      bool rewind_ok = false;

      while(mtry-- > 0 && !rewind_ok) {
	if(wait != 0) {
	  std::cout << "Waiting " << wait << " seconds." << std::endl;
	  sleep(wait);
	  wait *= 2;
	}
	else
	  wait = 1;

	try {
	  std::cout << "Rewinding file." << std::endl;
	  event.toBegin();
	  rewind_ok = true;
	}
	catch(...) {
	  rewind_ok = false;
	}

	if(rewind_ok)
	  std::cout << "Rewind succeeded." << std::endl;
	else
	  std::cout << "Rewind failed." << std::endl;
      }

      // If rewind failed, break out of overall retry loop.  Reposition has failed.

      if(!rewind_ok)
	break;
    }
  }
  return ok;
}

// Filter CRT hit collection.
// Remove bad CRT hits.
// Sort CRT hits into increasing time order.

void crt::CRTMergerExtra::filter_crt_hits(const std::vector<crt::CRTHit>& input_hits, 
				     std::vector<crt::CRTHit>& output_hits) const
{
  // Copy hits from input collectio to output collection, removing bad hits.

  output_hits.reserve(input_hits.size());

  for(const auto& crthit : input_hits) {

    // Filter out hits with bad times.

    if(crthit.ts0_s > 1300000000)
      output_hits.push_back(crthit);
  }

  // Sort hits into increaseing time order.

  std::sort(output_hits.begin(), output_hits.end(), CRTHitComp);
}


DEFINE_ART_MODULE(crt::CRTMergerExtra)
