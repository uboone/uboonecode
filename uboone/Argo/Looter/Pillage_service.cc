///

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/InputSourceDescription.h"

#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Persistency/Provenance/MasterProductRegistry.h"

#include "cetlib/exception.h"
#include "art/Framework/Core/InputSourceFactory.h"
#include "art/Framework/IO/Root/RootInput.h"


#include "art/Framework/Principal/Handle.h"
#include "RecoBase/Hit.h"

namespace microboone
{
  class Pillage {
  public:
    Pillage(fhicl::ParameterSet const& pset,
  			     art::ActivityRegistry& reg);
    virtual ~Pillage() {};
    
    void postBeginJob();
    art::ActivityRegistry& areg;
  };

  Pillage::Pillage(fhicl::ParameterSet const& pset,
			     art::ActivityRegistry& reg):
           areg(reg)
  {
    //   evdb::DisplayWindow::Register("Test1","Test display #1",600,900,mk_canvas1);
    //   evdb::DisplayWindow::OpenWindow(0);

    // this->reconfigure(pset);
    reg.sPostBeginJob.watch       (this, &Pillage::postBeginJob);
    std::cout << "PILLAGE!" << std::endl;
    

  }

  //......................................................................
  void Pillage::postBeginJob() 
  {    
    std::cout << "postBeginJob" << std::endl;
    
    std::string filename = "standard_reco_uboone.root";
    
    fhicl::ParameterSet main_input;
    main_input.put("module_type", "RootInput");
    main_input.put("module_label","source");
    main_input.put("max_events",-1);
    std::vector<std::string> source_list;
    source_list.push_back(filename); // Add the file.
    main_input.put("fileNames", source_list);
  
    art::MasterProductRegistry preg;
  
    art::ModuleDescription md(main_input.id(),
                              "RootInput",
                              "source",
                               art::ProcessConfiguration("blah",main_input.id(),"release","passid"));
    art::InputSourceDescription isd(md, preg, areg);
    std::unique_ptr<art::InputSource> source(art::InputSourceFactory::make(main_input,isd).release());

    art::RootInput* ri = dynamic_cast<art::RootInput*>(source.get());
    if(ri) std::cout << "Got a file, fuckers!" << std::endl;
  }
  


} 

DECLARE_ART_SERVICE(microboone::Pillage, LEGACY)

DEFINE_ART_SERVICE(microboone::Pillage)

