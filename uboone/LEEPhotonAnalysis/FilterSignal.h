

#ifndef FILTERSIGNAL_H
#define FILTERSIGNAL_H


#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"



class FilterSignal {
 public:
  bool RunOld(art::Event const & e, size_t & delta_mct_index);
  bool RunSingle(art::Event const & e, size_t const & mct_index, size_t & exiting_photon_index);
  bool Run(art::Event const & e, size_t & delta_mct_index);
};


#endif
