#ifndef CRT_MERGER_EXTRA_HH
#define CRT_MERGER_EXTRA_HH

#include <string>
#include <set>
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "gallery/Event.h"

namespace crt
{
  class CRTMergerExtra : public art::EDProducer
  {
  public:

    // Constructor, destructor.

    CRTMergerExtra(const fhicl::ParameterSet&);

    // Overrides.

    void produce( art::Event &evt ) override;
		
  private:
		
    // FCL parameters.

    bool fDebug;
    std::string fDAQHeaderTimeUBooNELabel;
    std::string fCRTHitLabel;

    // CRT files that we have seen (for sam metadata).

    std::set<std::string> fCRTSwizzledFiles;		
  };
}
#endif // CRT_MERGER_EXTRA_HH
