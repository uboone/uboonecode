#include "bernfebdaq-core/BernZMQBinarySource/BernZMQBinaryInputStreamReader.hh"

#include <iostream>
using namespace std;
bernfebdaq::BernZMQBinaryInputStreamReader::BernZMQBinaryInputStreamReader(fhicl::ParameterSet const& ps)
{
	fNFEBBuffers = ps.get<unsigned int>("NFEBBuffers",2);
	fVerbosity = ps.get<int>("Verbosity",-1);
	fMaxTSDiff = ps.get<unsigned int>("MaxTSDiff",20);
	fMaxStuckEvents = ps.get<unsigned int>("MaxStuckEvents",20);
	fPullTimeLatency = ps.get<uint32_t>("PullTimeLatency",500);
	fFragsPerSecond = ps.get<uint32_t>("FragsPerSecond",200);
	fLastPullTime = 0;
	fEventsPerPullCounter = 0;
	fNanosecondsPerFragment = 1000000000 / fFragsPerSecond;
}

art::Timestamp bernfebdaq::BernZMQBinaryInputStreamReader::ConvertMyTimebType( struct mytimeb_t const& tb)
{
	return ( ( (uint64_t)(tb.time) << 32 ) +
	   ( (uint64_t)(tb.millitm)*1000000 ) );
}

art::Timestamp bernfebdaq::BernZMQBinaryInputStreamReader::GetNextPullTime()
{
  bernfebdaq::BernZMQEvent zmq_ev;
  zmq_ev.mac5 = 0x1234;

  uint64_t current_position = fInputStream->tellg();

  while(true){
  
    fInputStream->read((char*)(&zmq_ev),sizeof(bernfebdaq::BernZMQEvent));
    
    if(!fInputStream)
      return art::Timestamp(0);
    
    if(zmq_ev.MAC5()==0xffff) break;
  }

  fInputStream->seekg(current_position, fInputStream->beg);

  mytimeb_t time_pull_started  = *((mytimeb_t*)((char*)(zmq_ev.adc)+sizeof(int)));
  return ConvertMyTimebType(time_pull_started);
}

art::Timestamp bernfebdaq::BernZMQBinaryInputStreamReader::SpecialEventTimestamp(BernZMQEvent const& zmq_ev)
{
	if (zmq_ev.MAC5()!=0xffff) return 0;
	
	mytimeb_t time_pull_started  = *((mytimeb_t*)((char*)(zmq_ev.adc)+sizeof(int)));
	auto ts = ConvertMyTimebType(time_pull_started);
	
	//mytimeb_t time_poll_finished = *((mytimeb*)((char*)(zmq_ev.adc)+sizeof(int)+sizeof(struct mytimeb)));
	if(fVerbosity>0)
	std::cout << "--Pulled " << fEventsPerPullCounter
		<< " events from "
		<< fLastPullTime.timeHigh() << " " << fLastPullTime.timeLow()
		<< " to "
		<< ts.timeHigh() << " " << ts.timeLow() << std::endl;
	
	fEventsPerPullCounter=0;
	fLastPullTime = ts;
	
	for(auto const& buffer : fFEBBuffers)
	{
		auto mac = buffer.first;
		auto cbuf = fFEBCurrentBuffer[mac];
		if(buffer.second[cbuf].size()>0)
		fPullLastEvent[mac] = buffer.second[cbuf].back().Time_TS0();
		else
		fPullLastEvent[mac] = 1000000000;
	}
	return ts;
}

void bernfebdaq::BernZMQBinaryInputStreamReader::InitializeBuffer(uint16_t mac)
{
	fFEBCurrentBuffer[mac]=0;
	fFEBBuffers[mac].resize(fNFEBBuffers);
}

bool bernfebdaq::BernZMQBinaryInputStreamReader::KeepCurrentEvent(std::vector<BernZMQEvent>& myFEBBuffer,
								  BernZMQEvent const& zmq_ev,
								  bool & last_pull_found)
{
	if(myFEBBuffer.size()<1) return true;
	
	last_pull_found = false;
	auto const& last_ts1_ref = fLastTS1Ref[zmq_ev.MAC5()];
	auto const& last_pull_event = fPullLastEvent[zmq_ev.MAC5()];
	
	uint32_t prev1_ts0 = myFEBBuffer.back().Time_TS0();
	uint32_t prev1_ts1 = myFEBBuffer.back().Time_TS1();
	bool prev_gps =  myFEBBuffer.back().IsReference_TS0();
	
	//things aren't out of order, or we saw the gps on the previous go
	if(zmq_ev.Time_TS0()>prev1_ts0 || prev_gps) return true;
	
	//if we've previous recognized this event ...
	if(fFEBStuckEventMaps[zmq_ev.MAC5()][zmq_ev.Time_TS0()]==zmq_ev.Time_TS1())
	return false;
	
	uint32_t diff_0_1_ts0 = (1000000000 + zmq_ev.Time_TS0() - prev1_ts0);
	uint32_t diff_0_1_ts1 = (0x40000000 + zmq_ev.Time_TS1() - prev1_ts1);
	if(zmq_ev.Time_TS1() > prev1_ts1) diff_0_1_ts1 = zmq_ev.Time_TS1() - prev1_ts1;
	if(zmq_ev.Time_TS0() > prev1_ts0) diff_0_1_ts0 = zmq_ev.Time_TS0() - prev1_ts0;
	int diff_diff_0_1 = (long int)(diff_0_1_ts0) - (long int)(diff_0_1_ts1);
	
	if(fVerbosity>0)
	std::cout << "\tROLLOVER CHECK: (buf_size=" << myFEBBuffer.size() << ")"
		<< "\n\t" << prev1_ts0 << " " << prev1_ts1
		<< "\n\t" << zmq_ev.Time_TS0() << " " << zmq_ev.Time_TS1()
		<< "\n\t"
		<< " " << diff_0_1_ts0 << " " << diff_0_1_ts1 << " " << diff_diff_0_1
		<< std::endl;
	
	if(std::abs(diff_diff_0_1)>(int)fMaxTSDiff)
	{
		size_t max_multi_stuck_events = fMaxStuckEvents;
		if(myFEBBuffer.size()<fMaxStuckEvents)
		max_multi_stuck_events = myFEBBuffer.size();
		
		bool current_is_bad=true;
		uint32_t prev_ts0,prev_ts1,diff_ts0,diff_ts1;
		uint32_t pts0,pts1;
		
		int diff_diff;
		for(size_t i_e=1; i_e<max_multi_stuck_events; ++i_e)
		{
			prev_ts0 = myFEBBuffer[myFEBBuffer.size()-i_e].Time_TS0();
			prev_ts1 = myFEBBuffer[myFEBBuffer.size()-i_e].Time_TS1();
			
			if(prev_ts0==last_pull_event)
			{
				last_pull_found = true;
				if(fVerbosity>0) std::cout << "\tFOUND LAST PULL EVENT." << std::endl;
			}
			
			diff_ts0 = (1000000000   + zmq_ev.Time_TS0() - prev_ts0);
			diff_ts1 = (last_ts1_ref + zmq_ev.Time_TS1() - prev_ts1);
			
			if(zmq_ev.Time_TS0() > prev_ts0) diff_ts0 = zmq_ev.Time_TS0() - prev_ts0;
			if(zmq_ev.Time_TS1() > prev_ts1) diff_ts1 = zmq_ev.Time_TS1() - prev_ts1;
			
			diff_diff = (long int)(diff_ts0) - (long int)(diff_ts1);
			
			//OK, so, if we find a matching one ...
			if(std::abs(diff_diff)<=(int)fMaxStuckEvents)
			{
				current_is_bad=false;
				for(size_t i_r=1; i_r<i_e; ++i_r)
				{
					pts0 = myFEBBuffer.back().Time_TS0();
					pts1 = myFEBBuffer.back().Time_TS1();
					fFEBStuckEventMaps[zmq_ev.MAC5()][pts0] = pts1;
					
					if(fVerbosity>0)
					std::cout << "\t\tInserting -" << i_r << " (" << pts0 << "," << pts1
						<< ") as stuck event in " << zmq_ev.MAC5() << std::endl;
					myFEBBuffer.pop_back();
				}
				return true;
				//we want to post this event to the buffer
			}
		}
		
		if(current_is_bad)
		{
			fFEBStuckEventMaps[zmq_ev.MAC5()][zmq_ev.Time_TS0()] = zmq_ev.Time_TS1();
			if(fVerbosity>0)
			std::cout << "\t\tInserting 0 (" << zmq_ev.Time_TS0() << "," << zmq_ev.Time_TS1()
			 << ") as stuck event in " << zmq_ev.MAC5() << std::endl;      
			return false;
		}
	}
	else
	{
		size_t max_multi_stuck_events = fMaxStuckEvents;
		if(myFEBBuffer.size()<fMaxStuckEvents)
		max_multi_stuck_events = myFEBBuffer.size();
		
		uint32_t prev_ts0;
		for(size_t i_e=1; i_e<max_multi_stuck_events; ++i_e)
		{
			prev_ts0 = myFEBBuffer[myFEBBuffer.size()-i_e].Time_TS0();
			if(prev_ts0==last_pull_event)
			{
				last_pull_found = true;
				if(fVerbosity>0) std::cout << "\tFOUND LAST PULL EVENT." << std::endl;
			}
		}  
	}
	return true;
} 

void bernfebdaq::BernZMQBinaryInputStreamReader::ClearNextBuffer(BernZMQEvent const& zmq_ev)
{
  size_t i_next_buffer = fFEBCurrentBuffer[zmq_ev.MAC5()]+1;
  if(i_next_buffer==fNFEBBuffers) i_next_buffer=0;
  fFEBBuffers[zmq_ev.MAC5()][i_next_buffer].clear();
  fFEBCurrentBuffer[zmq_ev.MAC5()] = i_next_buffer;
}

void bernfebdaq::BernZMQBinaryInputStreamReader::ClearNextBufferAndInsert(BernZMQEvent const& zmq_ev)
{
  size_t i_next_buffer = fFEBCurrentBuffer[zmq_ev.MAC5()]+1;
  if(i_next_buffer==fNFEBBuffers) i_next_buffer=0;
  fFEBBuffers[zmq_ev.MAC5()][i_next_buffer].clear();
  fFEBBuffers[zmq_ev.MAC5()][i_next_buffer].push_back(zmq_ev);
  fFEBCurrentBuffer[zmq_ev.MAC5()] = i_next_buffer;
}

void bernfebdaq::BernZMQBinaryInputStreamReader::FillFragments(std::vector<BernZMQEvent> const& myFEBBuffer,
							       FragMap_t & fragMap)
{
  //get the right time stamp here...
  art::Timestamp pulltime = GetNextPullTime();
  uint64_t pulltime_ns = ((uint64_t)pulltime.timeHigh())*1000000000 + pulltime.timeLow() - fPullTimeLatency*1000000;
  art::Timestamp mytime = ( (pulltime_ns/1000000000) << 32 ) + ( (pulltime_ns%1000000000) & 0xffffffff );

  uint32_t sec_strt = mytime.timeHigh()-1;
  uint32_t sec_end  = mytime.timeHigh();
  
  auto const& last_ev = myFEBBuffer.back();

  if(fVerbosity>0)
    std:: cout << last_ev.MAC5() << " : "
	       << pulltime.timeHigh() << " " << pulltime.timeLow() << " "
	       << mytime.timeHigh() << " " << mytime.timeLow()
	       << " ( " << last_ev.Time_TS0() << ", " << last_ev.Time_TS1() << " )"
	       << std::endl;

  int time_correction_diff = 0;
  if(last_ev.IsReference_TS0())
    time_correction_diff = (int)(last_ev.Time_TS0()) - 1000000000;

  //now let's go through and make all the fragments
  uint32_t f_nsec_strt=0, f_nsec_end=fNanosecondsPerFragment;
  uint32_t f_sec_end=sec_strt;

  size_t i_ev=0;
  for(size_t i_f=0; i_f<fFragsPerSecond; ++i_f){

    //f_nsec_end += fNanosecondsPerFragment;
    if( (1000000000 - f_nsec_end)<fNanosecondsPerFragment ){
      f_sec_end = sec_end;
      f_nsec_end = 0;
    }

    BernZMQFragmentMetadata metadata(sec_strt,f_nsec_strt,
				     f_sec_end,f_nsec_end,
				     time_correction_diff, 0,
				     0,0,
				     last_ev.MAC5(),0,
				     32,12);
    uint64_t frag_time =  ( ((uint64_t)metadata.time_start_seconds() << 32 ) + metadata.time_start_nanosec() );

    size_t i_frag_start=i_ev;
    uint32_t end_time_ns = f_nsec_end;
    if(end_time_ns==0) end_time_ns=1000000000;
    while(i_ev < myFEBBuffer.size()){
      if( GetCorrectedTime(myFEBBuffer[i_ev].Time_TS0(),metadata) > f_nsec_end )
	break;

      metadata.inc_Events();
      ++i_ev;
    }

    std::unique_ptr<artdaq::Fragment> myfragptr =
      artdaq::Fragment::FragmentBytes(metadata.n_events()*sizeof(BernZMQEvent),
				      metadata.time_start_seconds(),metadata.feb_id(),
				      bernfebdaq::detail::FragmentType::BernZMQ, metadata,
				      metadata.time_start_seconds());
    artdaq::Fragment myfrag(*myfragptr);

    if(fVerbosity>0)
      std::cout << " There are " << fragMap.count(frag_time) << " fragment collections at " << frag_time << std::endl;

    if(fragMap.count(frag_time)==0){
      if(fVerbosity>0)
	std::cout << "\tCreating new frag list ..." << std::endl;
      fragMap[frag_time].reset( new std::vector<artdaq::Fragment>);
    }
    fragMap[frag_time]->emplace_back( myfrag );
    std::copy((char*)(&myFEBBuffer[i_frag_start]),
	      (char*)(&myFEBBuffer[i_frag_start])+metadata.n_events()*sizeof(BernZMQEvent),
	      (char*)(fragMap[frag_time]->back().dataBegin()));
    if(fVerbosity>0)
      std::cout << " Now there are " << fragMap[frag_time]->size() << " fragments in that collection." << std::endl;

    f_nsec_end += fNanosecondsPerFragment;
    f_nsec_strt += fNanosecondsPerFragment;
  }
}

art::Timestamp bernfebdaq::BernZMQBinaryInputStreamReader::ReadUntilSpecialEvent(FragMap_t & fragMap)
{
	bernfebdaq::BernZMQEvent zmq_ev;
	
	if(fVerbosity>0)
	std::cout << "Reading file..." << std::endl;
	
	while(*fInputStream)
	{
		fInputStream->read((char*)(&zmq_ev),sizeof(bernfebdaq::BernZMQEvent));
		
		auto ts = SpecialEventTimestamp(zmq_ev);
		if(ts!=0) return ts;
		
		if(fFEBBuffers.count(zmq_ev.MAC5())==0)
		{
			if(fVerbosity>0)
			std::cout << "Initializing buffer for mac=" << zmq_ev.MAC5() << std::endl;
			InitializeBuffer(zmq_ev.MAC5());
		}
		
		size_t i_buffer = fFEBCurrentBuffer[zmq_ev.MAC5()];
		auto & myFEBBuffers = fFEBBuffers[zmq_ev.MAC5()];
		auto & myFEBBuffer = myFEBBuffers[i_buffer];
		
		bool last_pull_found;
		if (!KeepCurrentEvent(myFEBBuffer,zmq_ev,last_pull_found))
		continue;
		
		fEventsPerPullCounter++;
		
		if(fVerbosity>1)
		std::cout << zmq_ev.hdr_str() << std::endl;
		
		if(zmq_ev.IsReference_TS1())
		fLastTS1Ref[zmq_ev.MAC5()] = zmq_ev.Time_TS1();
		else if(zmq_ev.IsOverflow_TS1())
		fLastTS1Ref[zmq_ev.MAC5()] = 0x40000000;
		
		uint32_t prev1_ts0 = 0;
		bool     prev_gps =  false;
		if(myFEBBuffer.size()>0)
		{
			myFEBBuffer.back().IsReference_TS0();
			myFEBBuffer.back().Time_TS0();
		}
		
		if(zmq_ev.IsReference_TS0())
		{
			myFEBBuffer.push_back(zmq_ev);
			if(fVerbosity>0)
			std::cout << "FOUND GPS PULSE FOR FEB " << std::endl;
			FillFragments(myFEBBuffer,fragMap);
			ClearNextBuffer(zmq_ev);
		}
		else if(zmq_ev.Time_TS0()<prev1_ts0 && !prev_gps && !last_pull_found)
		{
			if(fVerbosity>0)
			std::cout << "FOUND GPS ROLLOVER FOR FEB " << std::endl;
			FillFragments(myFEBBuffer,fragMap);
			ClearNextBufferAndInsert(zmq_ev);
		}
		else
		myFEBBuffer.push_back(zmq_ev);
	}
	return 0;
}
