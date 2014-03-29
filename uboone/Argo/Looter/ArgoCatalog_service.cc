#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/FileServiceInterfaces/CatalogInterface.h"
#include <vector>
#include "art/Framework/Services/FileServiceInterfaces/FileDeliveryStatus.h"
#include "art/Framework/Services/FileServiceInterfaces/CatalogInterface.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/detail/ServiceHelper.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Principal/Handle.h"
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
  };



  ArgoCatalog::ArgoCatalog(fhicl::ParameterSet const& pset,
			     art::ActivityRegistry& reg)
  {
    std::cout << "ArgoCatalog!" << std::endl;
    

  }
  
  void ArgoCatalog::doConfigure(std::vector<std::string> const & item)
  {
    std::cout << "AgrgoCatalog doConfigure " << item.size() << std::endl;
    
  }

  static int iter=-1;

  int  ArgoCatalog::doGetNextFileURI(std::string & uri, double & waitTime){
    std::vector<std::string> list;
    list.push_back("file1.root");
    list.push_back("file2.root");
    list.push_back("file3.root");
    list.push_back("file3.root");
    iter++;
    iter = iter%list.size();
    uri = list[iter];
    std::cout << "AgrgoCatalog doGetNextFileURI returning " << uri << std::endl;
    return art::FileDeliveryStatus::SUCCESS;
  }
  
  void ArgoCatalog::doUpdateStatus(std::string const & uri, art::FileDisposition status){
    std::cout << "AgrgoCatalog doUpdateStatus " << uri << std::endl;
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
  bool ArgoCatalog::doIsSearchable() { return true; } 
  void ArgoCatalog::doRewind() {}   
  
  
  
  
  
  
  
} 
DECLARE_ART_SERVICE_INTERFACE_IMPL(microboone::ArgoCatalog, art::CatalogInterface, LEGACY)
DEFINE_ART_SERVICE_INTERFACE_IMPL(microboone::ArgoCatalog, art::CatalogInterface)
  

