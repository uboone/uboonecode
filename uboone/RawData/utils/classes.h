#include "art/Persistency/Common/Wrapper.h"

#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"
#include <string>
#include <utility>
#include <vector>

template class art::Wrapper< raw::ubdaqSoftwareTriggerData >;
//template class std::pair<std::basic_string<char>,bool>;
template class std::pair<std::string,bool>;
template class std::vector< std::pair<std::string,bool> >;

