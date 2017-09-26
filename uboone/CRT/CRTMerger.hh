#ifndef CRT_MERGER_HH
#define CRT_MERGER_HH


#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "fhiclcpp/ParameterSet.h"

#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"

#include <string>
#include <istream>

namespace crt{


  class CRTMerger : public art::EDProducer
  {
    std::vector<std::string> fFileNames;


    //The gallery events to check for the  crt data
    gallery::Event fCRTEvent;

    // Producer tag of the CRT events
    art::InputTag fTag;

    // Time window
    unsigned fTimeWindow;

    //See which T0 to use
    bool fUseT0;

  public:

    CRTMerger(const fhicl::ParameterSet&);

    ~CRTMerger();

    virtual void beginRun(art::Run &){}
 
    virtual void endRun(art::Run &){}

    virtual void produce (art::Event&);

    void reconfigure(fhicl::ParameterSet const & p) override;

  };

}


#endif // CRT_MERGER_HH
