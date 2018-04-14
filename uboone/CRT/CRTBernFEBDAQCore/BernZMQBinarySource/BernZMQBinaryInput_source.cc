#include "art/Framework/Core/InputSourceMacros.h" 
#include "art/Framework/IO/Sources/Source.h" 
#include "art/Framework/IO/Sources/SourceTraits.h" 

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/IO/Sources/SourceHelper.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "canvas/Persistency/Provenance/SubRunID.h"
#include "art/Framework/IO/Sources/put_product_in_principal.h"

#include "../Overlays/BernZMQFragment.hh"
#include "../Overlays/FragmentType.hh"
#include "artdaq-core/Data/Fragment.hh"

#include "BernZMQBinaryInputStreamReader.hh"

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>

namespace bernfebdaq {
  // Forward declaration of detail class.
  class BernZMQBinaryInputDetail;
}
////////////////////////////////////////////////////////////////////////
// Any declarations of specializations of Source traits (see
// art/Framework/IO/Sources/SourceTraits.h) should go here.
//
// ...
//
////////////////////////////////////////////////////////////////////////

// Definition of detail class.
class bernfebdaq::BernZMQBinaryInputDetail {

public: 
  BernZMQBinaryInputDetail(fhicl::ParameterSet const &pset,
			   art::ProductRegistryHelper &helper,
			   art::SourceHelper const &src_hlpr);
  
  void readFile(std::string const & filename, art::FileBlock*& fb);

  bool readNext(art::RunPrincipal const* const inR,
		art::SubRunPrincipal const* const inSR,
		art::RunPrincipal*& outR,
		art::SubRunPrincipal*& outSR,
		art::EventPrincipal*& outE);
  void closeCurrentFile();



private:
  
  art::SourceHelper  fSourceHelper;

  std::string        fModuleLabel;
  std::string        fInstanceLabel;

  std::ifstream      fInputStreamList;
  
  size_t             fEventNumber;
  size_t             fPrevRunNumber;
  size_t             fPrevSubRunNumber;

  size_t             fFragsPerEvent;

  fhicl::ParameterSet fConfigPSet;
  
  std::vector<BernZMQBinaryInputStreamReader> fStreamReaders;
  std::vector<std::ifstream>                  fInputStreams;
  std::vector<art::Timestamp>                 fInputStreamLastPullTime;
  
  FragMap_t          fFragMap;

  void writeEvent(FragMap_t::value_type & e,
		  art::RunPrincipal const* const,
		  art::SubRunPrincipal const* const,
		  art::RunPrincipal*& outR,
		  art::SubRunPrincipal*& outSR,
		  art::EventPrincipal*& outE);
  size_t getNextFileIndex();
  
};

bernfebdaq::BernZMQBinaryInputDetail::BernZMQBinaryInputDetail(fhicl::ParameterSet const & ps,
							       art::ProductRegistryHelper & helper,
							       art::SourceHelper const & src_hlpr)
  : fSourceHelper(src_hlpr),
    fModuleLabel(ps.get<std::string>("ModuleLabel")),
    fInstanceLabel(ps.get<std::string>("InstanceLabel")),
    fEventNumber(0),
    fPrevRunNumber(0),
    fPrevSubRunNumber(0),
    fFragsPerEvent(ps.get<size_t>("FragsPerEvent")),
    fConfigPSet(ps.get<fhicl::ParameterSet>("BernZMQBinaryInputStreamReaderConfig"))
{ 
	helper.reconstitutes< std::vector<artdaq::Fragment>, art::InEvent >(fModuleLabel,fInstanceLabel);
}


void bernfebdaq::BernZMQBinaryInputDetail::readFile(std::string const & filename, art::FileBlock*& fb)
{
	fb = new art::FileBlock(art::FileFormatVersion(1, "RawEvent2011"), filename);
	
	fInputStreamList.open(filename.c_str(),std::ios_base::in);
	char file_name[512];
	
	while(fInputStreamList.getline(file_name,512))
	{
		std::cout << "Got file " << file_name << std::endl;
		fInputStreams.emplace_back(file_name,std::ios_base::in | std::ios_base::binary);
	}
	for(size_t i_f=0; i_f<fInputStreams.size(); ++i_f)
	{
		fStreamReaders.emplace_back(fConfigPSet,fInputStreams[i_f]);
		fInputStreamLastPullTime.emplace_back(0);    
	}
	std::cout << "Opened input streams: " << fInputStreams.size() << std::endl;
	std::cout << "Opened StreamReaders: " << fStreamReaders.size()<< std::endl;
}

void bernfebdaq::BernZMQBinaryInputDetail::closeCurrentFile()
{
	for(auto & stream : fInputStreams)
	stream.close();
	fInputStreamList.close();
}

size_t bernfebdaq::BernZMQBinaryInputDetail::getNextFileIndex()
{
	return (size_t)(std::distance(fInputStreamLastPullTime.begin(),std::min_element(fInputStreamLastPullTime.begin(),fInputStreamLastPullTime.end())));
}

void bernfebdaq::BernZMQBinaryInputDetail::writeEvent(FragMap_t::value_type & e,
						      art::RunPrincipal const* const,
						      art::SubRunPrincipal const* const,
						      art::RunPrincipal*& outR,
						      art::SubRunPrincipal*& outSR,
						      art::EventPrincipal*& outE)
{
	uint64_t const& frag_time = e.first;
	
	size_t run_num = 1;
	size_t subrun_num = 1;
	
	if(run_num!=fPrevRunNumber)
	{
		outR = fSourceHelper.makeRunPrincipal(run_num,frag_time);
		outSR = fSourceHelper.makeSubRunPrincipal(run_num,subrun_num,frag_time);
		fPrevRunNumber = run_num;
		fPrevSubRunNumber =subrun_num;
	}
	else if(subrun_num!=fPrevSubRunNumber)
	{
		outSR = fSourceHelper.makeSubRunPrincipal(run_num,subrun_num,frag_time);
		fPrevSubRunNumber =subrun_num;
	}
	
	outE = fSourceHelper.makeEventPrincipal(run_num,
					  subrun_num,
					  fEventNumber++,
					  frag_time);
  
	art::put_product_in_principal(std::move(e.second),*outE, fModuleLabel, fInstanceLabel);
}

bool bernfebdaq::BernZMQBinaryInputDetail::readNext(art::RunPrincipal const* const inR,
						    art::SubRunPrincipal const* const inSR,
						    art::RunPrincipal*& outR,
						    art::SubRunPrincipal*& outSR,
						    art::EventPrincipal*& outE)
{
	size_t i_file = getNextFileIndex();
	
	if(!fInputStreams[i_file])
	{
		if(fFragMap.size()==0)
		return false;
		auto & e = *(fFragMap.begin());
		writeEvent(e,inR,inSR,outR,outSR,outE);      
		fFragMap.erase(fFragMap.begin());
		return true;
	}
	
	while(fInputStreams[i_file])
	{
		//std::cout << "Total of " << fFragMap.size() << " fragment events registered for file # " << i_file << std::endl;
		for(FragMap_t::iterator i_fm=fFragMap.begin(); i_fm!=fFragMap.end(); i_fm++)
		{
			if(i_fm->second->size()<fFragsPerEvent) continue;
			auto & e = *(fFragMap.begin());
			writeEvent(e,inR,inSR,outR,outSR,outE);      
			fFragMap.erase(fFragMap.begin());
			return true;
		}
		auto ts = fStreamReaders[i_file].ReadUntilSpecialEvent(fFragMap);
		
		long test = fInputStreams[i_file].tellg();
		std::cout<<"test: "<<test<<std::endl;
		
		std::cout << "File " << i_file << ":\treached pull event ... "
			<< ts.timeHigh() << ", " << ts.timeLow()
			<< "\t(" << fFragMap.size() << " fragment events registered)"
			<< std::endl;
		
		if(ts!=0)
		fInputStreamLastPullTime[i_file] = ts;
		else
		{
			std::cout<< "\tFile is closed? " << fInputStreams[i_file].eof() << " "
				<< fInputStreams[i_file].rdstate() << std::endl;
			if(fInputStreams[i_file].eof())
			fInputStreamLastPullTime[i_file] = 0xffffffffffffffff;
			else if(fInputStreams[i_file].fail())
			fInputStreams[i_file].clear();
		}
		i_file = getNextFileIndex();
		std::cout<<"i_file: "<<i_file<<std::endl;
		std::cout<<"................"<<std::endl;
	}
	
	if(fFragMap.size()==0)
	return false;
	auto & e = *(fFragMap.begin());
	writeEvent(e,inR,inSR,outR,outSR,outE);      
	fFragMap.erase(fFragMap.begin());
	return true;
}
 
// Optional typedef.
namespace bernfebdaq
{
  using BernZMQBinaryInputSource = art::Source<BernZMQBinaryInputDetail>;
}

// Define the input source to the art system.
DEFINE_ART_INPUT_SOURCE(bernfebdaq::BernZMQBinaryInputSource)
