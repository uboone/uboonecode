/**
  \defgroup crt All things Cosmic Ray Tagger related
**/

#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"
#include "uboone/CRT/CRTData.hh"
#include "uboone/CRT/CRTProducts/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include <vector>

template class std::vector<crt::CRTData>;
template class art::Wrapper< std::vector<crt::CRTData> >;
template class art::Wrapper< crt::CRTData >;
template class art::Assns<crt::CRTHit,recob::OpFlash,void>;
template class art::Assns<recob::OpFlash,crt::CRTHit,void>;
template class art::Wrapper<art::Assns<crt::CRTHit,recob::OpFlash,void> >;
template class art::Wrapper<art::Assns<recob::OpFlash,crt::CRTHit,void> >;
