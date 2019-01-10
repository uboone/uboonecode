#include <vector>
#include <sstream>
#include <iostream>

#include "uboone/CRT/CRTMerger.hh"
#include "uboone/CRT/CRTFileManager.h"
#include "uboone/RawData/utils/DAQHeaderTimeUBooNE.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "art/Framework/Principal/Handle.h"

crt::CRTMerger::CRTMerger(const fhicl::ParameterSet& pset) :
  fDebug(pset.get<bool>("debug")),
  fDAQHeaderTimeUBooNELabel(pset.get<std::string>("DAQHeaderTimeUBooNELabel")),
  fCRTHitLabel(pset.get<std::string>("CRTHitLabel")),
  fTimeStart(pset.get<long double>("TimeStart")),
  fTimeEnd(pset.get<long double>("TimeEnd")),
  fTimeOffset(pset.get<long double>("TimeOffset"))
{
  std::cout << "CRTMerger module constructed." << std::endl;
  produces< std::vector<crt::CRTHit> >();
}

void crt::CRTMerger::produce(art::Event& event)
{
  if (fDebug) {
    std::cout <<"crt::CRTMerger::produce, NEW EVENT" << std::endl;
  }

  // Get CRT file manager service.

  art::ServiceHandle<crt::CRTFileManager> crtman;

  // Get DAQ Header GPS time.

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
  
  // Calculate event time and merge window as floating point seconds.

  long double TPCtime = evt_time_sec + 1.e-9L * evt_time_nsec;
  long double MergingWindow_start = TPCtime + fTimeOffset + fTimeStart;
  long double MergingWindow_end   = TPCtime + fTimeOffset + fTimeEnd;

  // Get matching swizzled CRT files.

  std::vector<std::string> crtrootfiles = crtman->findMatchingCRTFiles(evtTimeGPS);

  // Purge CRT hit cache of files not included in above list.

  std::set<std::string> cache_files;
  for(auto const& item : fCRTHitCache)
    cache_files.insert(item.first);
  for(auto const& file_name : crtrootfiles)
    cache_files.erase(file_name);
  for(auto const& file_name : cache_files) {
    std::cout << "Erasing file " << file_name << " from hit cache." << std::endl;
    fCRTHitCache.erase(file_name);
  }
  
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

    // Get CRT hit cache for this file.

    std::map<long long, std::vector<crt::CRTHit> >& file_cache = fCRTHitCache[crtrootfile];

    // Purge CRT hit cache of nonrelevant entries.

    for(auto const& item : file_cache) {
      long long entry = item.first;
      const crt::CRTHit& first_hit = item.second.front();
      const crt::CRTHit& last_hit = item.second.back();
      long double first_time = first_hit.ts0_s + 1.e-9L * first_hit.ts0_ns;
      long double last_time = last_hit.ts0_s + 1.e-9L * last_hit.ts0_ns;
      if(first_time > MergingWindow_end || last_time < MergingWindow_start) {
	std::cout << "Erasing entry " << entry << " from hit cache." << std::endl;
	file_cache.erase(entry);
      }
    }

    // Summarize file cache contents.

    std::cout << "Hit cache contains " << file_cache.size() << " entries." << std::endl;
    for(auto const& item : file_cache)
      std::cout << "Entry " << item.first << std::endl;

    if(file_cache.size() != 0) {

      // If at this point the cache is not empty, extend the cache by one entry 
      // into the future, if necessary.
      // Check last hit of last entry.

      long long last_entry = file_cache.rbegin()->first;
      const crt::CRTHit& last_hit = file_cache.rbegin()->second.back();
      long double last_time = last_hit.ts0_s + 1.e-9L * last_hit.ts0_ns;
      if(last_time < MergingWindow_end) {

	// add one entry to cache.

	long long new_entry = last_entry + 1;
	bool rewind_ok = true;
	if(new_entry < crt_event.eventEntry())
	  rewind_ok = crtman->rewind(crt_event);
	if(rewind_ok) {
	  while(!crt_event.atEnd() && crt_event.eventEntry() < new_entry)
	    crt_event.next();
	  std::vector<crt::CRTHit> new_hits;
	  if(!crt_event.atEnd())
	    crtman->get_crt_hits(crt_event, new_hits);
	  if(new_hits.size() != 0) {
	    long long entry = crt_event.eventEntry();
	    std::cout << "Adding entry " << entry << " to hit cache." << std::endl;
	    file_cache.emplace(entry, std::move(new_hits));
	  }
	}
      }
    }
    else {

      // The file hit cache is empty.
      // Do an ab initio positioning of the file and populate the hit cache.

      // Starting from the current file position, get the next nonempty collection of crt hits.

      std::cout << "CRT event entry before reposition = " << crt_event.eventEntry() << std::endl;
      std::vector<crt::CRTHit> CRTHitCollection;
      crtman->get_crt_hits(crt_event, CRTHitCollection);

      // If we didn't find anything, rewind file and try again.

      if(CRTHitCollection.size() == 0) {
	bool rewind_ok = crtman->rewind(crt_event);
	if(rewind_ok)
	  crtman->get_crt_hits(crt_event, CRTHitCollection);
      }

      // If we still got nothing, skip this file.  This file has no hits.

      if(CRTHitCollection.size() == 0) {
	std::cout << "No CRT hits in file." << std::endl;
	continue;
      }

      // Find the time of the earliest crt hit.
    
      long double CRTtime = CRTHitCollection[0].ts0_s + 1.e-9L * CRTHitCollection[0].ts0_ns;

      // If the earliest hit time falls after the merging window, rewind the file and try again.

      if(CRTtime > MergingWindow_end) {
	CRTHitCollection.erase(CRTHitCollection.begin(), CRTHitCollection.end());
	bool rewind_ok = crtman->rewind(crt_event);
	if(rewind_ok)
	  crtman->get_crt_hits(crt_event, CRTHitCollection);

	// We should always find some hits, since we found some before.

	if(CRTHitCollection.size() == 0) {
	  std::cout << "No CRT hits in file on second attempt." << std::endl;
	  throw cet::exception("CRTFileMerger") << "No CRT hits";
	}

	// New earliest time.

	CRTtime = CRTHitCollection[0].ts0_s + 1.e-9L * CRTHitCollection[0].ts0_ns;
      }

      // If the earliest crt hit time is still after the merging window, then there
      // are no mergable hits in this file.  Skip file.

      if(CRTtime > MergingWindow_end) {
	std::cout << "No mergable hits in file." << std::endl;
	continue;
      }

      // Now we know for sure that there are some CRT hits before the end of the merging window.
      // If the tpc event is far in the future, skip over some events
      // without reading the hits.  Keep doing this until we get close to the merging
      // window.

      while(MergingWindow_start - CRTtime > 2.) {
	long double nskip = 0.8 * (MergingWindow_start - CRTtime - 2.);
	if(nskip > 1.e6) {

	  // Something wrong...

	  std::cout << "Ridiculously large time difference between CRT hit time " << CRTtime
		    << " and start of merging window " << MergingWindow_start << std::endl;
	  throw cet::exception("CRTFileMerger") << "Suspicious times";
	}
	long long new_entry = crt_event.eventEntry() + (long long)nskip;
	long long skip_entries = new_entry - crt_event.eventEntry();
	if(skip_entries == 0)
	  break;
	//std::cout << "Skipping " << skip_entries << " events." << std::endl;
	while(!crt_event.atEnd() && crt_event.eventEntry() < new_entry)
	  crt_event.next();
	crtman->get_crt_hits(crt_event, CRTHitCollection);

	// If we got no hits here, there are no more mergable hits in this file.
	// Skip further processing on this file.

	if(CRTHitCollection.size() == 0) {
	  std::cout << "No mergable hits in file." << std::endl;
	  break;
	}

	// New earliest time.

	CRTtime = CRTHitCollection[0].ts0_s + 1.e-9L * CRTHitCollection[0].ts0_ns;
      }
      std::cout << "CRT event entry after reposition = " << crt_event.eventEntry() << std::endl;
      CRTHitCollection.erase(CRTHitCollection.begin(), CRTHitCollection.end());

      // Now we are ready for a detailed look at CRT hits starting from the current collection.
      // Loop over entries starting from the current position and populate the hit cache.
      // Keep going until the merging window is completely covered or we reach the end of 
      // the file.

      for(; !crt_event.atEnd(); crt_event.next()) {
	std::vector<crt::CRTHit> crt_hits;
	crtman->get_crt_hits(crt_event, crt_hits);

	// An empty collection signals end of file.

	if(crt_hits.size() == 0)
	  break;

	// Get the time of the first and last hit in this entry.

	long double first_time = crt_hits.front().ts0_s + 1.e-9L * crt_hits.front().ts0_ns;	
	long double last_time = crt_hits.back().ts0_s + 1.e-9L * crt_hits.back().ts0_ns;	

	// If this entry's time range overlaps with the merging window, add it to the cache.

	long long entry = crt_event.eventEntry();
	if(first_time <= MergingWindow_end && last_time >= MergingWindow_start) {
	  std::cout << "Adding entry " << entry << " to hit cache." << std::endl;
	  file_cache.emplace(entry, std::move(crt_hits));
	}
	else
	  std::cout << "Skipping entry " << entry << std::endl;

	// If the last hit time is after the end of the merging window,
	// the merging window is covered.  We can quit.

	if(last_time > MergingWindow_end)
	  break;
      }
    }

    // At this point, the CRT hit cache should contain all hits needed for merging this event.
    // Loop over all hits in the cache and merge them.

    bool done = false;
    int merging = 0;
    for(const auto& item : file_cache) {

      const std::vector<crt::CRTHit>& CRTHitCollection = item.second;

      if(fDebug) {
	uint32_t first_sec = CRTHitCollection.front().ts0_s;
	uint32_t last_sec = CRTHitCollection.back().ts0_s;
	std::cout << "\nTPC event time = " << evt_time_sec << std::endl;
	std::cout << "First event time = " << first_sec << std::endl;
	std::cout << "Last event time = " << last_sec << std::endl;
      }

      // Check whether any hits in this entry are relevant for merging.
      // If not skip this entry.

      const crt::CRTHit& first_hit = CRTHitCollection.front();
      const crt::CRTHit& last_hit = CRTHitCollection.back();
      long double first_time = first_hit.ts0_s + 1.e-9L * first_hit.ts0_ns;
      long double last_time = last_hit.ts0_s + 1.e-9L * last_hit.ts0_ns;
      if(first_time > MergingWindow_end || last_time < MergingWindow_start)
	continue;

      // Loop over CRT hits in this entry.

      for(auto const& CRTHitevent : CRTHitCollection) {
	long double CRTtime = CRTHitevent.ts0_s + 1.e-9L * CRTHitevent.ts0_ns;
	//std::cout << CRTtime-MergingWindow_start << std::endl;
	if (CRTtime >= MergingWindow_start && CRTtime <= MergingWindow_end) {
	  if (fDebug)
	    std::cout<<"found match"<<std::endl;
	  CRTHitEventsSet->push_back(CRTHitevent);
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
    std::cout << "CRT event entry after processing = " << crt_event.eventEntry() << std::endl;
    std::cout<<"# merging in the stream: "<<merging<<std::endl;
  } // end loop on ROOT File
  std::cout << "# of Merged CRTHits in CRTHitEventsSet being written to event:: "
	    <<CRTHitEventsSet->size()<<std::endl;
  event.put(std::move(CRTHitEventsSet));
  
  if (fDebug)
    std::cout<<"---X---"<<std::endl;
}


DEFINE_ART_MODULE(crt::CRTMerger)
