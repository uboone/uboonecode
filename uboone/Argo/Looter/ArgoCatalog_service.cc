#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/FileServiceInterfaces/CatalogInterface.h"
#include <vector>
#include "art/Framework/Services/FileServiceInterfaces/FileDeliveryStatus.h"
#include "art/Framework/Services/FileServiceInterfaces/CatalogInterface.h"
#include "art/Framework/Services/FileServiceInterfaces/FileDisposition.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/detail/ServiceHelper.h"
#include "art/Framework/IO/Root/RootInput.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Persistency/Provenance/BranchIDListHelper.h"

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Event.h"
#include "RecoBase/Hit.h"

namespace microboone
{
  class ArgoCatalog : public art::CatalogInterface {
  public:
    ArgoCatalog(fhicl::ParameterSet const& pset,
  			       art::ActivityRegistry& reg);
    virtual ~ArgoCatalog() {};

  private:
    virtual void doConfigure(std::vector<std::string> const & item);
    virtual int  doGetNextFileURI(std::string & uri, double & waitTime);
    virtual void doUpdateStatus(std::string const & uri, art::FileDisposition status);
    virtual void doOutputFileOpened(std::string const & module_label);
    virtual void doOutputModuleInitiated(std::string const & module_label,
                                         fhicl::ParameterSet const & pset);
    virtual void doOutputFileClosed(std::string const & module_label,
                                    std::string const & fileFQname);
    virtual void doEventSelected(std::string const & module_label,
                                 art::EventID const & event_id,
                                 art::HLTGlobalStatus const & acceptance_info);
    virtual bool doIsSearchable();
    virtual void doRewind();    
    
    void preSourceRun();
    void preSourceEvent();
    void postBeginJob();
    void postBeginJobWorkers(art::InputSource* inputs,
                             std::vector<art::Worker*> const& workers);
    void preProcessEvent(art::Event const&);
    void postProcessEvent(art::Event const&);
    art::InputSource* inputSource_; ///< Input source of events
    std::string nextFile_;
    
    size_t iter_;
    size_t eventsToSkip_;
    bool   skipRestOfFile_;
  };



  ArgoCatalog::ArgoCatalog(fhicl::ParameterSet const& pset,
			     art::ActivityRegistry& reg)
             : inputSource_(0)
             , nextFile_("file1.root")
             , iter_(0)
             , eventsToSkip_(0)
             , skipRestOfFile_(false)
  {
    std::cout << "ArgoCatalog!" << std::endl;
    reg.sPreSourceRun.watch       (this, &ArgoCatalog::preSourceRun);
    reg.sPreSource.watch          (this, &ArgoCatalog::preSourceEvent);
    reg.sPostBeginJob.watch       (this, &ArgoCatalog::postBeginJob);
    reg.sPostBeginJobWorkers.watch(this, &ArgoCatalog::postBeginJobWorkers);
    reg.sPreProcessEvent.watch    (this, &ArgoCatalog::preProcessEvent);
    reg.sPostProcessEvent.watch   (this, &ArgoCatalog::postProcessEvent);
    

  }
  
  void ArgoCatalog::doConfigure(std::vector<std::string> const & item)
  {
    std::cout << "AgrgoCatalog doConfigure " << item.size() << std::endl;
    
  }


  int  ArgoCatalog::doGetNextFileURI(std::string & uri, double & waitTime){
    
    if(iter_==0) uri = "file1.root";
    if(iter_==1) uri = "file2.root";
    if(iter_==2) uri = "file3.root";
    if(iter_>2)  uri = "file4.root";
    std::cout << "ArgoCatalog doGetNextFileURI returning " << uri << std::endl;
    iter_++;
    std::cout << "  and attempting to skip to event " << iter_ << std::endl;
    eventsToSkip_= iter_;
    skipRestOfFile_ = false;
    // art::BranchIDListHelper::clearRegistries();
   
    return art::FileDeliveryStatus::SUCCESS;
  }
  
  void ArgoCatalog::doUpdateStatus(std::string const & uri, art::FileDisposition status){
    std::cout << "AgrgoCatalog doUpdateStatus " << uri  
              <<  " status:" << translateFileDisposition(status) << std::endl;
  }

  void ArgoCatalog::doOutputFileOpened(std::string const & module_label)
  {
    std::cout << "AgrgoCatalog doOutputFileOpened " << module_label << std::endl;
    
  }
  void ArgoCatalog::doOutputModuleInitiated(std::string const & module_label,
                               fhicl::ParameterSet const & pset)
   {
     std::cout << "AgrgoCatalog doOutputModuleInitiated " << module_label << std::endl;
     
   }
                                 
                                 
  void ArgoCatalog::doOutputFileClosed(std::string const & module_label,
                          std::string const & fileFQname)
  {
    std::cout << "AgrgoCatalog doOutputFileClosed " << module_label << " " <<fileFQname << std::endl;
    
  }
  void ArgoCatalog::doEventSelected(std::string const & module_label,
                       art::EventID const & event_id,
                       art::HLTGlobalStatus const & acceptance_info) {
                         
                         std::cout << "AgrgoCatalog doEventSelected " << module_label << " "  << std::endl;
                         
                         
                       }
  bool ArgoCatalog::doIsSearchable() { return false; }  // looks important - segfaults if turned to true.
  void ArgoCatalog::doRewind() {}   
  
  void ArgoCatalog::preSourceRun()
  {
    std::cout << "ArgoCatalog preSourceRun()" << std::endl;
    art::RootInput* rootInput = dynamic_cast<art::RootInput*>(inputSource_);
    if(rootInput) std::cout << "Have the InputSource." << std::endl;
    if(rootInput){ // fails if first event through.
     // std::cout << "Resetting next event..." << iter_ << std::endl;
      // art::EventID nextEvent(1,0,iter_);
      // rootInput->seekToEvent(nextEvent);

    // rootInput->skipEvents(iter_); skip
    }

  }
  
  void ArgoCatalog::preSourceEvent()
  {
    std::cout << "ArgoCatalog preSourceEvent()" << std::endl;
      // art::RootInput* rootInput = dynamic_cast<art::RootInput*>(inputSource_);
  }
  
  void ArgoCatalog::postBeginJobWorkers(art::InputSource* input_source,
                                        std::vector<art::Worker*> const&)
  {
    // Called beore the input source is called for the first time.
    std::cout << "ArgoCatalog postBeginJobWorkers" << std::endl;
    inputSource_ = input_source;
    // art::RootInput* rootInput = dynamic_cast<art::RootInput*>(inputSource_);
    // assert(rootInput);
    // art::EventID nextEvent(999999,999999,999999);
    // rootInput->seekToEvent(nextEvent);
    
  }


  void ArgoCatalog::postBeginJob() {
    std::cout << "ArgoCatalog postBeginJob" << std::endl;
  }

  void ArgoCatalog::preProcessEvent(art::Event const&)
  {
    std::cout << "ArgoCatalog preProcessEvent" << std::endl;
  }
  
  void ArgoCatalog::postProcessEvent(art::Event const& event)
  {
    // We may not be at the correct event. Let's see how long it takes to skip some.
    if(eventsToSkip_) {
      std::cout << eventsToSkip_ << " events left to skip" << std::endl;
      eventsToSkip_--;
      return;
    } 
    if(skipRestOfFile_) {
      std::cout << "skipping to end of file" << std::endl;
      return;
    }

    
    
    // This is where I can make the json.
    std::cout << "ArgoCatalog postProcessEvent" << std::endl;
    std::cout << event.run() 
              <<  "|" << event.subRun()  
              <<  "|" << event.id().event()    << std::endl;   
                                       
    typedef art::Handle< std::vector<recob::Hit> > hitHandle_t;
    std::vector<hitHandle_t> list_of_hitlists;
    try{
      event.getManyByType(list_of_hitlists);
    }catch(...) {
      std::cout << "problem getting types." << std::endl;
    }

    std::cout << "Got " << list_of_hitlists.size() << " types of hits." << std::endl;
    for(size_t i=0;i<list_of_hitlists.size(); i++) {
      std::cout << "--->" << list_of_hitlists[i].provenance()->moduleLabel() << std::endl;
    }
    
    
    skipRestOfFile_ = true;
    
    // Force the next file to load.
    // art::RootInput* rootInput = dynamic_cast<art::RootInput*>(inputSource_);
    // assert(rootInput);
    // art::EventID nextEvent(999999,999999,999999);
    // rootInput->seekToEvent(nextEvent);
  }

  
  
} 
DECLARE_ART_SERVICE_INTERFACE_IMPL(microboone::ArgoCatalog, art::CatalogInterface, LEGACY)
DEFINE_ART_SERVICE_INTERFACE_IMPL(microboone::ArgoCatalog, art::CatalogInterface)
  

