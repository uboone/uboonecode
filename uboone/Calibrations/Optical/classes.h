#include "art/Persistency/Common/Wrapper.h"

#include "SWTriggerBase/Result.h"

#include <vector>

template class std::vector<trigger::Result>;
template class art::Wrapper<std::vector<trigger::Result> >;
