#ifndef MSetCRTFrag_hh_
#define MSetCRTFrag_hh_
#include <cstdint>
#include "bernfebdaq-core/Overlays/BernZMQFragment.hh"

namespace crt
{
	class MSetCRTFrag
	{
		uint32_t fSec;
                uint32_t fTPC_ns;
                uint32_t fCRT_frag_beg_ns;
                uint32_t fCRT_frag_end_ns;
	     public:
		MSetCRTFrag();
                //MSetCRTFrag(bernfebdaq::BernZMQFragment f, unsigned tpc_s, unsigned tpc_ns, unsigned crt_frag_beg_ns, unsigned crt_frag_end_ns);
		MSetCRTFrag(artdaq::Fragment const& i_crt_frag);
                virtual ~MSetCRTFrag();
		/*
	     private:
		bernfebdaq::BernZMQFragment fFrag;
		*/
	};
}

#endif
