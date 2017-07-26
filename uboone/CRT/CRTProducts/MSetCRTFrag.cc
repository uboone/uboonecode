#include "uboone/CRT/CRTProducts/MSetCRTFrag.hh"

namespace crt
{
	MSetCRTFrag::MSetCRTFrag()
	{
	}
	/*
	MSetCRTFrag::MSetCRTFrag(bernfebdaq::BernZMQFragment frag, uint32_t GPS_s, uint32_t TPC_ns, uint32_t crt_frag_beg_ns, uint32_t crt_frag_end_ns): fSec(GPS_s), fTPC_ns(TPC_ns), fCRT_frag_beg_ns(crt_frag_beg_ns), fCRT_frag_end_ns(crt_frag_end_ns) 
	{	
	}
	*/
	MSetCRTFrag::MSetCRTFrag(artdaq::Fragment const & f)
	{
		bernfebdaq::BernZMQFragment testfrag(f);	
	}
	MSetCRTFrag::~MSetCRTFrag()
        {
        }	
}


