////////////////////////////////////////////////////////////////////////
// Class:       CRTFileManager
// Plugin Type: service (art v2_05_01)
// File:        CRTFileManager_service.cc
//
// Generated at Thu Jan  3 14:09:49 2019 by Herbert Greenlee using cetskelgen
// from cetlib version v1_21_00.
//
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <ctime>
#include "CRTFileManager.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "boost/date_time/c_local_time_adjustor.hpp"
#include "IFDH_service.h"


namespace {

  // Local function to compare CRT hits according to time.

  bool CRTHitComp(const crt::CRTHit& crthit1, const crt::CRTHit& crthit2)
  {
    bool result = (crthit1.ts0_s < crthit2.ts0_s ||
		   (crthit1.ts0_s == crthit2.ts0_s && crthit1.ts0_ns < crthit2.ts0_ns));
    return result;
  }
}

// Constructor.

crt::CRTFileManager::CRTFileManager(fhicl::ParameterSet const & p, art::ActivityRegistry & areg) :
  fDebug(p.get<bool>("debug")),
  fMaxFiles(p.get<unsigned int>("maxFiles")),
  fCRTHitLabel(p.get<std::string>("CRTHitLabel")),
  fCRTVersion(p.get<std::string>("ubversion_CRTHits"))
{
  setenv("TZ", "CST6CDT", 1);  // Fermilab time zone.
  tzset();
}

// Use sam to find swizzled CRT files that match the event time stamp.

std::vector<std::string> crt::CRTFileManager::findMatchingCRTFiles(art::Timestamp event_time) const
{
  // Result vector.

  std::vector<std::string> crtrootfiles;

  // Get ifdh service.

  art::ServiceHandle<ifdh_ns::IFDH> ifdh;

  // Convert the art time stamp into a string format that is understandable to sam.

  unsigned long evt_time_sec = event_time.timeHigh();
  unsigned long evt_time_nsec = event_time.timeLow();
  
  unsigned long time_tpc = evt_time_sec*1000000000 + evt_time_nsec;   // Nanoseconds
  unsigned long time_tpc1= time_tpc/1000;                             // Microseconds.
  
  const char* tz = getenv("TZ");
  std::string tzs(tz);
  if (fDebug)
    std::cout<<"time-zone "<<tzs<<std::endl;
  
  if(tz != 0 && *tz != 0) {
    std::string tzs(tz);
    if(tzs != std::string("CST6CDT")) {
      // Timezone is wrong, throw exception.
      throw cet::exception("CRTFileManager") << "Wrong timezone: " << tzs;
    }
  }
  else {
    // Timezone is not set.  Throw exception.
    throw cet::exception("CRTFileManager") << "Timezone not set.";
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

  // Check cached files.

  boost::posix_time::time_duration one_minute(0, 1, 0, 0);
  for(const auto& crtfileinfo : fCRTFiles) {
    const std::string& crtfile = crtfileinfo.first;    // crt swizzled file.
    const CRTFileInfo& fileinfo = crtfileinfo.second;  // crt binary file.
    if(this_event_localtime >= fileinfo.fStartTime + one_minute && 
       this_event_localtime <= fileinfo.fEndTime - one_minute) {
      std::cout << "Found cached file " << crtfile << std::endl;
      for(auto const& crt_swizzled : fileinfo.fSwizzled)
	crtrootfiles.push_back(crt_swizzled);
    }
  }
  if(crtrootfiles.size() < 6) {

    // Didn't match enough cached files.

    crtrootfiles.erase(crtrootfiles.begin(), crtrootfiles.end());

    // Query CRT binery files.

    std::string stringTime = boost::posix_time::to_iso_extended_string(this_event_localtime);
    stringTime = "'"+stringTime+"'";
    std::ostringstream dim;
    dim << "file_format " << "crt-binaryraw"
	<<" and file_type " << "data"
	<<" and start_time < " << stringTime 
	<< " and end_time > " << stringTime;
  
    if (fDebug)
      std::cout<<"dim = "<<dim.str()<<std::endl;
  
    // List those crtdaq files:
    std::vector< std::string > crtfiles = ifdh->translateConstraints(dim.str());
    if(fDebug)
      std::cout << "Found " << crtfiles.size() << " CRT binary files." << std::endl;

    // Get metadata of binary CRT files.

    for(auto const & crtfile : crtfiles) {
      if(fDebug)
	std::cout << "\n" << crtfile << std::endl;
      std::string md = ifdh->getMetadata(crtfile);

      // Extract CRT binary file start and end time from metadata as strings.

      std::string start_time;
      size_t n = md.find("Start Time: ");
      size_t m = 0;
      if(n < std::string::npos) {
	n += 12;
	m = md.find("+", n);
	if(m > n && m < std::string::npos) {
	  start_time = md.substr(n, m-n);
	  if(fDebug)
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
	  if(fDebug)
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
      if(fDebug) {
	std::cout << "Start ptime " << start_time2 << std::endl;
	std::cout << "End ptime " << end_time2 << std::endl;
      }
      if(start_time2 != start_time || end_time2 != end_time) {
	std::cout << "Start ptime " << start_time2 << std::endl;
	std::cout << "End ptime " << end_time2 << std::endl;
	throw cet::exception("CRTFileManager") << "Problem converting start and end time.";
      }
  
      // Query CRT swizzled files that are children of this CRT binary file.

      std::ostringstream dim1;
      dim1 << "file_format " << "artroot"
	   <<" and ub_project.version " << fCRTVersion
	   << " and ischildof: (file_name " << crtfile
	   <<" with availability physical )";

      if (fDebug)
	std::cout << "dim1 = " << dim1.str() << std::endl;
    
      std::vector< std::string > tmprootfiles = ifdh->translateConstraints(dim1.str());
    
      std::cout << "Found " << tmprootfiles.size()
		<< " daughters of " << crtfile << std::endl;
      if (tmprootfiles.size()>0) {
	for (const auto& artrootchild : tmprootfiles) {
	  crtrootfiles.push_back(artrootchild);
	  if(fDebug)
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
	std::cout << "File metadata cache contains " << fCRTFiles.size() << " files." << std::endl;
      }
    }
  }
  std::cout<<"\nNumber of matching CRT swizzled files: "<<crtrootfiles.size()<<std::endl;
  if(fDebug) {
    for(const auto& crt_swizzled : crtrootfiles)
      std::cout << crt_swizzled << std::endl;
  }
  if (!crtrootfiles.size())
    std::cout << "\n\t CRTFileManager_module: No child CRT files found that conform to constraints: "
	      << "file_format " << "artroot" << " and ub_project.version "
	      << fCRTVersion << std::endl;

  // Throw exception if there are fewer than six CRT files.

  if(crtrootfiles.size() < 6) {
    throw cet::exception("CRTFileManager") << "Too few matching CRT files: " 
					   << crtrootfiles.size() << "\n";
  }

  // Done.

  return crtrootfiles;
}

// Return an open gallery file for the specified file.
// If we already have an open gallery file in the file cache, return that.
// Otherwise open a new file and add it to the cache.

gallery::Event& crt::CRTFileManager::openFile(std::string file_name)
{
  // Do we have an open gallery event for this file?
  // If yes, return it.

  for(auto& crtfile : fCRTEvents) {
    if(crtfile.first == file_name) {

      // Yes, reuse the open gallery event.

      std::cout << "Reusing open gallery event for file " << file_name << std::endl;
      return *(crtfile.second);
    }
  }

  // If we get here, we don't have a cached open file, so we will need to open one.

  // Make sure that the open file cache won't exceed the maximum size.

  if(fCRTEvents.size() >= fMaxFiles) {
    std::cout << "Closing file " << fCRTEvents.front().first << std::endl;
    fCRTEvents.pop_front();
  }

  // Get xrootd url of file.

  art::ServiceHandle<ifdh_ns::IFDH> ifdh;
  std::string schema = "root";
  std::vector< std::string > xrootd_urls;
  bool locate_ok = false;
  try {
    xrootd_urls = ifdh->locateFile(file_name, schema);
    locate_ok = true;
  }
  catch(...) {
    locate_ok = false;
  }
  if(!locate_ok || xrootd_urls.size() == 0) {
    std::cout << "No locations found for root file " << file_name << std::endl;
    throw cet::exception("CRTFileManager") << "Could not locate CRT file: " 
					   << file_name << "\n";
  }

  // Use first location only.

  if(xrootd_urls.size() > 1)
    xrootd_urls.erase(xrootd_urls.begin()+1, xrootd_urls.end());
  std::cout<<"xrootd URL: " << xrootd_urls[0] << std::endl;
  
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
  std::unique_ptr<gallery::Event> crt_event;

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
      crt_event = std::move(std::unique_ptr<gallery::Event>(new gallery::Event(xrootd_urls)));
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

  if(open_ok) {
    std::cout<< "Opened the CRT root file from xrootd URL" << std::endl;
    fCRTEvents.emplace_back(file_name, std::move(crt_event));
    std::cout << "New [collection of] CRTEvent[s]. Its size is  "
	      << fCRTEvents.back().second->numberOfEventsInFile() << "." << std::endl;
  }
  else {
    std::cout << "Failed to open CRT root file xrootd URL." << std::endl;
    throw cet::exception("CRTFileManager") << "Could not open CRT file: " 
					   << file_name << "\n";
  }

  // Done.

  return *(fCRTEvents.back().second);
}

// Rewind gallery file to first event.
// Return true if success, false if fail.

bool crt::CRTFileManager::rewind(gallery::Event& event) const
{
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

  // Done.

  return rewind_ok;
}

// Get next nonempty filtered and sorted CRT hit collection starting from 
// current gallery event.
// Return empty collection if end of file is reached.

void crt::CRTFileManager::get_crt_hits(gallery::Event& event,
				       std::vector<crt::CRTHit>& output_hits) const
{
  // Make sure output hit container is initially empty.

  output_hits.erase(output_hits.begin(), output_hits.end());

  // Event loop.

  for(; !event.atEnd() && output_hits.size() == 0; event.next()) {

    // Get crt hit handle for the current event..

    gallery::ValidHandle<std::vector<crt::CRTHit> > h =
      event.getValidHandle< std::vector<crt::CRTHit> >(fCRTHitLabel);

    // Copy and filter hits into output hits collection.

    output_hits.reserve(h->size());
    for(const auto& crthit : *h) {

      // Filter out hits with bad times.

      if(crthit.ts0_s > 1300000000)
	output_hits.push_back(crthit);
    }
    if(output_hits.size() != 0)
      break;
  }

  // Sort hits into increaseing time order.

  std::sort(output_hits.begin(), output_hits.end(), CRTHitComp);
}


DEFINE_ART_SERVICE(crt::CRTFileManager)
