#include "art/Framework/Art/artapp.h"

#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Art/BasicOutputOptionsHandler.h"
#include "art/Framework/Art/BasicSourceOptionsHandler.h"
#include "art/Framework/Art/DebugOptionsHandler.h"
#include "art/Framework/Art/FileCatalogOptionsHandler.h"
#include "art/Framework/Art/OptionsHandlers.h"
#include "art/Framework/Art/run_art.h"
#include "art/Utilities/FirstAbsoluteOrLookupWithDotPolicy.h"

#include "art/Framework/Art/BasicOptionsHandler.h"
#include "art/Framework/Art/BasicPostProcessor.h"
#include "art/Framework/Art/InitRootHandlers.h"
// #include "art/Framework/EventProcessor/EventProcessor.h"
#include "EventProcessor.h"
#include "art/Framework/Core/RootDictionaryManager.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceRegistry.h"
#include "art/Framework/Services/Registry/ServiceToken.h"
#include "art/Utilities/ExceptionMessages.h"
#include "art/Utilities/RootHandlers.h"
#include "art/Utilities/UnixSignalHandlers.h"
#include "cetlib/container_algorithms.h"
#include "cetlib/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/parse.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "boost/program_options.hpp"
#include "boost/regex.hpp"
#include "TError.h"


#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>

namespace bpo = boost::program_options;

namespace art {
  int run_phoenix(int argc, char* argv[]);

  int onePhoenix(fhicl::ParameterSet& main_pset);
}



int main( int argc, char* argv[] ) {
  return art::run_phoenix(argc,argv);
}

// -----------------------------------------------
namespace {
  struct RootErrorHandlerSentry {
    RootErrorHandlerSentry(bool reset) {
      art::setRootErrorHandler(reset);
    }
    ~RootErrorHandlerSentry() {
      SetErrorHandler(DefaultErrorHandler);
    }
  };

  class EventProcessorWithSentry {
  public:
    explicit EventProcessorWithSentry() : ep_(), callEndJob_(false) { }
    explicit EventProcessorWithSentry(std::unique_ptr<art::EventProcessor> && ep) :
    ep_(std::move(ep)),
    callEndJob_(false) { }
    EventProcessorWithSentry(EventProcessorWithSentry &&) = default;
    EventProcessorWithSentry &
      operator =(EventProcessorWithSentry &&) & = default;

    ~EventProcessorWithSentry() {
      if (callEndJob_ && ep_.get()) {
        try {
          ep_->endJob();
        }
        catch (cet::exception & e) {
          //art::printArtException(e, kProgramName);
        }
        catch (std::bad_alloc & e) {
          //art::printBadAllocException(kProgramName);
        }
        catch (std::exception & e) {
          //art::printStdException(e, kProgramName);
        }
        catch (...) {
          //art::printUnknownException(kProgramName);
        }
      }
    }
    void on() {
      callEndJob_ = true;
    }
    void off() {
      callEndJob_ = false;
    }

    art::EventProcessor * operator->() {
      return ep_.get();
    }
  private:
    std::unique_ptr<art::EventProcessor> ep_;
    bool callEndJob_;
  }; // EventProcessorWithSentry

} // namespace



namespace art{

  int run_phoenix(int argc, char* argv[])
  {
    //
    // from   // int result = artapp(argc,argv);
    //

    // Configuration file lookup policy.
    char const * fhicl_env = getenv("FHICL_FILE_PATH");
    std::string search_path;
    if (fhicl_env == nullptr) {
      std::cerr
        << "Expected environment variable FHICL_FILE_PATH is "
          << "missing or empty: using \".\"\n";
      search_path = ".:";
    }
    else {
      search_path = std::string(fhicl_env) + ":";
    }
    art::FirstAbsoluteOrLookupWithDotPolicy lookupPolicy(search_path);
    // Empty options_description.
    bpo::options_description in_desc;
    // Create and store options handlers.
    art::OptionsHandlers handlers;
    handlers.reserve(4); // -ish.
    // Add new handlers here. Do *not* add a BasicOptionsHandler: it will
    // be done for you.
    handlers.emplace_back(new art::BasicSourceOptionsHandler(in_desc));
    handlers.emplace_back(new art::BasicOutputOptionsHandler(in_desc));
    handlers.emplace_back(new art::DebugOptionsHandler(in_desc, true));
    handlers.emplace_back(new art::FileCatalogOptionsHandler(in_desc));
  
    // 
    // from 
    //   return art::run_art(argc, argv, all_desc, lookupPolicy, std::move(handlers));
    //
    std::ostringstream descstr;
    descstr << "Usage: "
      << argv[0]
        << " <-c <config-file>> <other-options> [<source-file>]+\n\n"
          << "Allowed options";
    bpo::options_description all_desc(descstr.str());
    all_desc.add(in_desc);
    // BasicOptionsHandler should always be first in the list!
    handlers.emplace(handlers.begin(),
    new BasicOptionsHandler(all_desc, lookupPolicy));
    handlers.emplace_back(new BasicPostProcessor);
    // This must be added separately: how to deal with any non-option arguments.
    bpo::positional_options_description pd;
    // A single non-option argument will be taken to be the source data file.
    pd.add("source", -1);
    // Parse the command line.
    bpo::variables_map vm;
    try {
      bpo::store(bpo::command_line_parser(argc, argv).options(all_desc).positional(pd).run(),
      vm);
      bpo::notify(vm);
    }
    catch (bpo::error const & e) {
      std::cerr << "Exception from command line processing in " << argv[0]
        << ": " << e.what() << "\n";
      return 7000;
    }
    // Preliminary argument checking.
    for (auto & handler : handlers) {
      auto result = handler->checkOptions(vm);
      if (result != 0) {
        return result;
      }
    }
    // Processing of arguments and post-processing of config.
    fhicl::intermediate_table raw_config;
    for (auto & handler : handlers) {
      auto result = handler->processOptions(vm, raw_config);
      if (result != 0) {
        return result;
      }
    }
    //
    // Make the parameter set from the intermediate table:
    //
    // fhicl::ParameterSet main_pset;
    // try {
    //   make_ParameterSet(raw_config, main_pset);
    // }
    // catch (cet::exception & e) {
    //   std::cerr << "ERROR: Failed to create a parameter set from parsed configuration with exception "
    //     << e.what()
    //       << ".\n";
    //   std::cerr << "       Intermediate configuration state follows:\n"
    //     << "------------------------------------"
    //       << "------------------------------------"
    //         << "\n";
    //   for (auto const & item : raw_config) {
    //     std::cerr << item.first << ": " << item.second.to_string() << "\n";
    //   }
    //   std::cerr
    //     << "------------------------------------"
    //       << "------------------------------------"
    //         << "\n";
    //   return 7003;
    // }
    // // Main parameter set must be placed in registry manually.
    // try {
    //   fhicl::ParameterSetRegistry::put(main_pset);
    // }
    // catch (...) {
    //   std::cerr << "Uncaught exception while inserting main parameter set into registry.\n";
    //   throw;
    // }
    //   
    // 
    // from 
    //  return run_art_common_(main_pset);
    //
  
  
    // fhicl::ParameterSet
    //   services_pset(main_pset.get<fhicl::ParameterSet>("services",
    // fhicl::ParameterSet()));
    // fhicl::ParameterSet
    //   scheduler_pset(services_pset.get<fhicl::ParameterSet>("scheduler",
    // fhicl::ParameterSet()));
    //
    // Start the messagefacility
    //
    mf::MessageDrop::instance()->jobMode = std::string("analysis");
    mf::MessageDrop::instance()->runEvent = std::string("JobSetup");
    // mf::StartMessageFacility(mf::MessageFacilityService::MultiThread,
    // services_pset.get<fhicl::ParameterSet>("message",
    // fhicl::ParameterSet()));
    mf::LogInfo("MF_INIT_OK") << "Messagelogger initialization complete.";
    //
    // Initialize:
    //   unix signal facility
    // art::setupSignals(scheduler_pset.get<bool>("enableSigInt", true));
    // //   init root handlers facility
    // if (scheduler_pset.get<bool>("unloadRootSigHandler", true)) {
    //   art::unloadRootSigHandler();
    // }
    // RootErrorHandlerSentry re_sentry(scheduler_pset.get<bool>("resetRootErrHandler", true));
    // Load all dictionaries.
    art::RootDictionaryManager rdm;
    art::completeRootHandlers();
    art::ServiceToken dummyToken;
    //
    // Now create the EventProcessor
    //
    int rc;
    while(1) {
            
      std::string filename;
      int event;
      std::cout << "Filename: ";
      std::cin >> filename;
      std::cout << std::endl << "Event: ";
      std::cin >> event;
      std::cout << std::endl << "Trying filename " << filename << std::endl;
      std::vector<std::string> source_list;
      source_list.push_back(filename);
      raw_config.put("source.module_type","RootInput");
      raw_config.put("source.module_label","source");
      raw_config.put("source.fileNames",source_list);
      raw_config.put("source.firstEvent",event);
      raw_config.put("source.maxEvents",1);
      
            
      fhicl::ParameterSet main_pset;
      make_ParameterSet(raw_config, main_pset);
      fhicl::ParameterSetRegistry::put(main_pset);
      
      std::cout << main_pset.to_indented_string() << "\n";

      rc = onePhoenix(main_pset);
    }
    
    
    // Make sure the message logger is shut down so our message really is
    // the last one put out.
    mf::MessageFacilityService::instance().MFPresence.reset();
    // End message.
    std::cout
      << "Art has completed and will exit with status "
        << rc
          << "."
            << std::endl;
    return rc;
  }




  int onePhoenix(fhicl::ParameterSet& the_pset)
  {
    EventProcessorWithSentry proc;
    int rc = -1;
    try {
      std::unique_ptr<art::EventProcessor>
        procP(new
          art::EventProcessor(the_pset));
      EventProcessorWithSentry procTmp(std::move(procP));
      proc = std::move(procTmp);
      proc->beginJob();
      proc.on();
      if (proc->runToCompletion() == EventProcessor::epSignal) {
        std::cerr << "Art caught and handled signal "
          << art::shutdown_flag
            << ".\n";
      }
      proc.off();
      proc->endJob();
      rc = 0;
    }
    catch (art::Exception & e) {
      rc = e.returnCode();
      art::printArtException(e, "art"); // , "Thing1", rc);
    }
    catch (cet::exception & e) {
      rc = 8001;
      art::printArtException(e, "art"); // , "Thing2", rc);
    }
    catch (std::bad_alloc & bda) {
      rc = 8004;
      art::printBadAllocException("art"); // , "Thing3", rc);
    }
    catch (std::exception & e) {
      rc = 8002;
      art::printStdException(e, "art"); // , "Thing4", rc);
    }
    catch (...) {
      rc = 8003;
      art::printUnknownException("art"); // , "Thing5", rc);
    }
    return rc;
  }



} // namespace

