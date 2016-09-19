#include <map>
#include <string>
#include "art/Persistency/Common/Wrapper.h"

#include "MCEventWeight.h"

template class art::Wrapper<evwgh::MCEventWeight>;
template class art::Wrapper<std::vector<evwgh::MCEventWeight> >;
template class std::map<std::string, double>;
template class std::map<std::string, std::vector<double> >;
