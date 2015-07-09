
#include "art/Persistency/Common/Wrapper.h"

#include "MuCSData.h"

template class std::vector<MuCS::MuCSData>;

template class art::Wrapper< MuCS::MuCSData >; 

template class art::Wrapper< std::vector<MuCS::MuCSData>  >;
