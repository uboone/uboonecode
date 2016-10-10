
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Common/Assns.h"

#include "lardata/RecoBase/OpFlash.h"
#include "lardata/RecoBase/OpHit.h"

template class std::pair< art::Ptr<recob::OpHit>,        art::Ptr<recob::OpFlash>    >;
template class std::pair< art::Ptr<recob::OpFlash>,      art::Ptr<recob::OpHit>      >;

template class art::Assns<recob::OpHit,      recob::OpFlash,    void>;
template class art::Assns<recob::OpFlash,    recob::OpHit,      void>;

template class art::Wrapper< art::Assns<recob::OpHit,      recob::OpFlash,    void> >;
template class art::Wrapper< art::Assns<recob::OpFlash,    recob::OpHit,      void> >;
