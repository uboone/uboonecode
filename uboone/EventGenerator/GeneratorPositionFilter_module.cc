////////////////////////////////////////////////////////////////////////
// Class:       GeneratorPositionFilter
// Plugin Type: filter (art v2_05_00)
// File:        GeneratorPositionFilter_module.cc
//
// Generated at Fri Nov 17 11:28:05 2017 by Giuseppe Cerati using cetskelgen
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
#include "nusimdata/SimulationBase/GTruth.h"

#include <memory>

class GeneratorPositionFilter;


class GeneratorPositionFilter : public art::EDFilter {
public:
  explicit GeneratorPositionFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  GeneratorPositionFilter(GeneratorPositionFilter const &) = delete;
  GeneratorPositionFilter(GeneratorPositionFilter &&) = delete;
  GeneratorPositionFilter & operator = (GeneratorPositionFilter const &) = delete;
  GeneratorPositionFilter & operator = (GeneratorPositionFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;

private:
  art::InputTag genInputTag;
  double minx;
  double maxx;
  double miny;
  double maxy;
  double minz;
  double maxz;
};


GeneratorPositionFilter::GeneratorPositionFilter(fhicl::ParameterSet const & p)
  :genInputTag(p.get<art::InputTag>("genInputTag")),
   minx(p.get<double>("minX")),
   maxx(p.get<double>("maxX")),
   miny(p.get<double>("minY")),
   maxy(p.get<double>("maxY")),
   minz(p.get<double>("minZ")),
   maxz(p.get<double>("maxZ"))
{
}

bool GeneratorPositionFilter::filter(art::Event & e)
{
  auto gt = e.getValidHandle<std::vector<simb::GTruth> >(genInputTag);
  for (auto v : (*gt) ) {
    auto x = v.fVertex.X();
    auto y = v.fVertex.Y();
    auto z = v.fVertex.Z();
    if (x>minx && x<maxx && y>miny && y<maxy && z>minz && z<maxz) return true;
  }
  return false;
}

DEFINE_ART_MODULE(GeneratorPositionFilter)
