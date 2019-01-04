#include <vector>
#include <sstream>
#include <iostream>

#include "uboone/CRT/CRTMergerExtra.hh"
#include "uboone/CRT/CRTFileManager.h"
#include "uboone/CRT/CRTProducts/CRTHit.hh"
#include "uboone/RawData/utils/DAQHeaderTimeUBooNE.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "art/Framework/Principal/Handle.h"

crt::CRTMergerExtra::CRTMergerExtra(const fhicl::ParameterSet& pset) :
  fDebug(pset.get<bool>("debug")),
  fDAQHeaderTimeUBooNELabel(pset.get<std::string>("DAQHeaderTimeUBooNELabel")),
  fCRTHitLabel(pset.get<std::string>("CRTHitLabel"))
{
  std::cout << "CRTMergerExtra module constructed." << std::endl;
  produces< std::vector<crt::CRTHit> >();
}

void crt::CRTMergerExtra::produce(art::Event& event)
{
  if (fDebug) {
    std::cout <<"crt::CRTMergerExtra::produce, NEW EVENT" << std::endl;
  }

  // Get CRT file manager service.

  art::ServiceHandle<crt::CRTFileManager> crtman;

  const int T0(5); // Number of event/seconds of CRT wrt TPC evt time
                   // more than which we won't check for merge candidates..
  
  //get DAQ Header GPS time.
  art::Handle< raw::DAQHeaderTimeUBooNE > rawHandle_DAQHeader;
  event.getByLabel(fDAQHeaderTimeUBooNELabel, rawHandle_DAQHeader);
  
  if(!rawHandle_DAQHeader.isValid()) {
    std::cout << "Run " << event.run() 
	      << ", subrun " << event.subRun()
	      << ", event " << event.event()
	      << " has zero DAQHeaderTimeUBooNE with label "
	      << fDAQHeaderTimeUBooNELabel << std::endl;
    return;
  }
  
  raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
  art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
  art::Timestamp evtTimeNTP = my_DAQHeader.ntp_time();
  
  std::cout<<"evt_timeGPS_sec " << evtTimeGPS.timeHigh()
	   <<",  evt_timeGPS_nsec " << evtTimeGPS.timeLow()
	   <<",  evt_timeNTP_sec " << evtTimeNTP.timeHigh()
	   <<",  evt_timeNTP_nsec "<< evtTimeNTP.timeLow() << std::endl;

  unsigned long evt_time_sec = evtTimeGPS.timeHigh();
  unsigned long evt_time_nsec = evtTimeGPS.timeLow();

  // Quit if GPS time is zero.
  // This is somehow normal for the first event of a run...
  if(evt_time_sec == 0 && evt_time_nsec == 0) {

    // We promised to put a collection of CRTHits into the event.

    std::unique_ptr<std::vector<crt::CRTHit> > p(new std::vector<crt::CRTHit>);
    event.put(std::move(p));
    return;
  }
  
  // Get matching swizzled CRT files.

  std::vector<std::string> crtrootfiles = crtman->findMatchingCRTFiles(evtTimeGPS);
  
  // Create collection for matching CRTHits for this event that we will fill.

  std::unique_ptr<std::vector<crt::CRTHit> > CRTHitEventsSet(new std::vector<crt::CRTHit>);

  // Loop over swizzled CRT files.
  
  for(auto const& crtrootfile : crtrootfiles) {
    std::cout << "\nThe child artroot file is " << crtrootfile << std::endl;

    // Add this file to set of seen CRT files for sam metadata.

    if(fCRTSwizzledFiles.count(crtrootfile) == 0) {
      std::cout << "Adding CRT parent file: " << crtrootfile << std::endl;
      art::ServiceHandle<art::FileCatalogMetadata> md;
      std::ostringstream ostr;
      ostr << "mixparent" << fCRTSwizzledFiles.size();
      md->addMetadataString(ostr.str(), crtrootfile);
      fCRTSwizzledFiles.insert(crtrootfile);
    }

    // Get open gallery file for this file.

    gallery::Event& crt_event = crtman->openFile(crtrootfile);
  
    // Attempt to reposition the CRT event.
    // In case of initial failure, rewind the file and retry once.

    std::cout << "CRT event entry before reposition = " << crt_event.eventEntry() << std::endl;
    bool ok = crtman->reposition(crt_event, evt_time_sec, T0);

    // If reposition failed, skip this file.

    if(!ok) {
      std::cout << "Failed to reposition file." << std::endl;
      continue;
    }

    std::cout << "CRT event entry after reposition = " << crt_event.eventEntry() << std::endl;

    int merging = 0;
    long double TPCtime = evt_time_sec + 1.e-9L * evt_time_nsec;
    long double MergingWindow_start = TPCtime - 0.002L;
    long double MergingWindow_end   = TPCtime + 0.004L;

    // Loop over CRT events.

    bool done = false;
    for(; !crt_event.atEnd(); ++crt_event) {

      // Loop over events in the current CRT hit collection.

      gallery::ValidHandle<std::vector<crt::CRTHit> > h =
	crt_event.getValidHandle< std::vector<crt::CRTHit> >(fCRTHitLabel);
      std::vector<crt::CRTHit> CRTHitCollection;
      crtman->filter_crt_hits(*h, CRTHitCollection);
      if(fDebug) {
	uint32_t first_sec = CRTHitCollection.front().ts0_s;
	uint32_t last_sec = CRTHitCollection.back().ts0_s;
	std::cout << "\nTPC event time = " << evt_time_sec << std::endl;
	std::cout << "First event time = " << first_sec << std::endl;
	std::cout << "Last event time = " << last_sec << std::endl;
      }
      for(auto const& CRTHitevent : CRTHitCollection) {
	long double CRTtime = CRTHitevent.ts0_s - 1 + 1.e-9L * CRTHitevent.ts0_ns;//one second substracted for adding 1 extra second
	if (CRTtime >= MergingWindow_start && CRTtime <= MergingWindow_end) {
	  if (fDebug)
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
  
  if (fDebug)
    std::cout<<"---X---"<<std::endl;
}


DEFINE_ART_MODULE(crt::CRTMergerExtra)
