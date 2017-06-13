////////////////////////////////////////////////////////////////////////
// Class:       QuietEventFilter
// Plugin Type: filter (art v2_05_00)
// File:        QuietEventFilter_module.cc
//
// Generated at Thu Feb 16 02:50:08 2017 by Taritree Wongjirad using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "lardataobj/RecoBase/Wire.h"

class QuietEventFilter;


class QuietEventFilter : public art::EDFilter {
public:
  explicit QuietEventFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  QuietEventFilter(QuietEventFilter const &) = delete;
  QuietEventFilter(QuietEventFilter &&) = delete;
  QuietEventFilter & operator = (QuietEventFilter const &) = delete;
  QuietEventFilter & operator = (QuietEventFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

private:

  // Declare member data here.
  std::string fWireModuleLabel;
  float fThreshold;
  float fMaxADCSum;

};


QuietEventFilter::QuietEventFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  fWireModuleLabel = p.get<std::string>("WireModuleLabel");
  fThreshold       = p.get<float>("WireThreshold");
  fMaxADCSum       = p.get<float>("MaxADCSum");
}

bool QuietEventFilter::filter(art::Event & e)
{
  // Implementation of required member function here.

  art::Handle< std::vector<recob::Wire> > wireHandle;
  e.getByLabel(fWireModuleLabel,wireHandle);

  std::vector<recob::Wire> const& wireVector(*wireHandle);

  float ADCSum = 0.;
  for ( auto const& wire : wireVector ) {
    if ( wire.View()!=geo::kZ ) continue;

    std::vector<float> signal = wire.Signal();

    for ( auto &val : signal ) {
      if ( val>fThreshold ) {
	ADCSum += val;
      }
    }

  }

  if ( ADCSum<fMaxADCSum )
    return true;

  return false;

}

DEFINE_ART_MODULE(QuietEventFilter)
