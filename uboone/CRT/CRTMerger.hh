#ifndef CRT_MERGER_HH
#define CRT_MERGER_HH

#include <string>
#include <set>
#include <map>
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "gallery/Event.h"
#include "uboone/CRT/CRTProducts/CRTHit.hh"

namespace crt
{
  class CRTMerger : public art::EDProducer
  {
  public:

    // Constructor.

    CRTMerger(const fhicl::ParameterSet&);

    // Overrides.

    void produce( art::Event &evt ) override;

  private:

    // FCL parameters.

    bool fDebug;
    std::string fDAQHeaderTimeUBooNELabel;
    std::string fCRTHitLabel;
    long double fTimeStart;    // Start of crt matching interval (seconds).
    long double fTimeEnd;      // End of crt matching interval (seconds).
    long double fTimeOffset;   // Offset of crt matching interval (second).

    // CRT files that we have seen (for sam metadata).

    std::set<std::string> fCRTSwizzledFiles;		

    // CRT hit cache.
    // _crthit_cache[file_name][entry][hit_number].

    std::map<std::string, std::map<long long, std::vector<crt::CRTHit> > > fCRTHitCache;
  };
}
#endif // CRT_MERGER_HH
