///
/// Notes:
/// This is the main() function for the new argo backend.
/// It requires 

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
#include "art/Framework/EventProcessor/EventProcessor.h"
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
#include "TTimeStamp.h"
#include "TSystem.h"


#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <signal.h>


#include "SocketServer.h"

#include "../Looter/LooterGlobals.h"

int run_larsoft(std::string filename,int event,std::string options);

namespace bpo = boost::program_options;

class MySocketServer : public SocketServer{
public:
  MySocketServer( int port ) : SocketServer(port) {};
  virtual ~MySocketServer() {};
  
  size_t SufficientData( const unsigned char* inData, size_t inSize, size_t inMaxSize) 
  { 
    static char options[100];
    static char filename[990];
    static char selection[990];
    Long64_t entrystart;
    Long64_t entryend;
    int bytes;
    int r =  sscanf((char*)inData,"%99[^,\n],%900[^,\n],%900[^,\n],%lld,%lld\n%n",options,filename,selection,&entrystart,&entryend,&bytes);
    if(r==5) return bytes;
    return 0;
  };
};

MySocketServer* ss = 0;

int gEventsServed = 0;
TTime gTimeStart;

// static std::string gLooterOutput = "";
// extern std::string gLooterOutput;




int main( int argc, char* argv[] ) {
  signal(SIGCHLD, SIG_IGN); // Tell system to kill zombie processes immediately, no need to wait()  
    int tcpPortNumber = 9092;
    if(argc>1) {
      sscanf(argv[1],"%d",&tcpPortNumber);
    }
    std::cout << argv[0] << " starting up at " <<  TTimeStamp().AsString() << " on port " << tcpPortNumber << std::endl;
  
    ss = new MySocketServer(tcpPortNumber);
    if(ss->Setup()) exit(1);  // Quit if socket won't bind.

    // write a PID file.
    {
      std::string pidfilename(argv[0]);
      pidfilename+=".pid";
      ofstream pidfile(pidfilename.c_str());
      pidfile << gSystem->GetPid();
      pidfile.close();
    }

    gTimeStart = gSystem->Now();
    
    // Do a manual load of libraries before we fork. This will help forking speed, I think.
    std::cout << "Loading libraries." << std::endl;
    art::RootDictionaryManager rdm;
    art::completeRootHandlers();
    std::cout << "...done" << std::endl;
    
    // run_larsoft("/Users/tagg/Argo/server/prod_piminus_0.1-2.0GeV_isotropic_3window_uboone_15314288_17_gen_15314703_17_g4_15342663_17_detsim_tpc_15368710_17_reco2D.root",0,"_NoPreSpill_NoPostSpill_");
    
    
    
    
    while (1) {
      //cout << ".";
      //cout.flush();
      ss->CheckForBadClients();
    
      unsigned char* dataRecvd;
      int   client;
      bool  newClient;
      ss->Listen(100., // seconds timeout
                1000, // bytes max
                dataRecvd, // data recieved in message
                client,    // fd number of client.
                newClient  // true if a new client connected.
                );
      if(newClient) {
        std::cout << "New client " << client << std::endl;
      }
    
      if(dataRecvd) {
        // Try to parse format of FILENAME,GATE\n
        char options[100];
        char filename[990];
        char selection[990];
        Long64_t entrystart;
        Long64_t entryend;
        int r =  sscanf((char*)dataRecvd,"%99[^,\n],%900[^,\n],%900[^,\n],%lld,%lld\n",options,filename,selection,&entrystart,&entryend);
        if(r==5) {
          //Successful conversion. Give it a try.
          std::cout << "Got a valid request at " << TTimeStamp().AsString() << std::endl;
          std::cout << "    Filename: --" << filename << "--" << std::endl;
          std::cout << "    Selection:--" << selection << "--" << std::endl;
          std::cout << "    From:     --" << entrystart << " to " << entryend << std::endl;
          std::cout << "    Options:  --" << options << std::endl;
          
          // fork a process to cope.
          pid_t pid = fork();
          if(pid ==0) {
            // Child process
            std::cout << "Child process: " << getpid() << std::endl;
            long t1 = gSystem->Now();
            // Now do your stuff.

            // FIXME: Correctly use the event file to convert the entrystart and entryend and selection variables into a unique event id.

            // Run Larsoft once.
            run_larsoft(filename,entrystart,options);            
            std::cout << "Sending: " << looterOutput() << std::endl;

            looterOutput().append("\n");
            long t2 = gSystem->Now();
            // Send it out.
            ss->SendTo(client, (unsigned char*)looterOutput().c_str(),  looterOutput().length() );
            std::cout << "Request served." << std::endl;
            long t3 = gSystem->Now();
          
            ss->Close(client);
            long t4 = gSystem->Now();

            std::cout << "Time to compose: " << t2-t1 << "  Time to Serve: " << t3-t2 << " Total: " << t4-t1 << std::endl;
            _exit(0);
          }
          
          gEventsServed++;
        }

      }
    }
    delete ss;
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




int run_larsoft(std::string filename, int event, std::string options)
{
  using namespace art;
    //
    // from   // int result = artapp(argc,argv);
    //
    
    
    // Do our own custom FCL file, hardcoded here.
    std::string myfcl = "\
      outputs: {}                                                                          \
      physics: { ana: [ looter ]                                                           \
                 analyzers: { looter: { module_label: looter                               \
                                        module_type: ArgoLooter                                \
                                      }                                                    \
                            }                                                              \
                 end_paths: [ ana ]                                                        \
                 filters: {}                                                               \
                 producers: {}                                                             \
               }                                                                           \
      process_name: ArgoLooter                                                             \
      services: {                                                                          \
        message: { destinations: { STDOUT: { categories: { ArtReport: { limit: 100 }       \
                                                           default:   { limit: -1 }        \
                                                         }                                 \
                                             threshold: INFO                               \
                                             type: cout                                    \
                                           }                                               \
                                }                                                          \
                 }                                                                         \
       scheduler: { defaultExceptions: false }                                             \
       user: { CatalogInterface: { service_provider: TrivialFileDelivery }                 \
               FileTransfer: { service_provider: TrivialFileTransfer }                     \
             }                                                                             \
      }                                                                                    \
      ";
    
    std::string search_path = ".:";
    art::FirstAbsoluteOrLookupWithDotPolicy lookupPolicy(search_path);
    
    fhicl::intermediate_table raw_config;
    fhicl::parse_document(myfcl, raw_config);
    
    // Configure the input source and the event to read.
    std::vector<std::string> source_list;
    source_list.push_back(filename);
    raw_config.put("source.module_type","RootInput");
    raw_config.put("source.module_label","source");
    raw_config.put("source.fileNames",source_list);
    // raw_config.put("source.firstEvent",event);
    raw_config.put("source.skipEvents",event);
    raw_config.put("source.maxEvents",1);
    
    raw_config.put("physics.analyzers.looter.options",options);
    
    
    
    // Start the messagefacility
    //
    mf::MessageDrop::instance()->jobMode = std::string("analysis");
    mf::MessageDrop::instance()->runEvent = std::string("JobSetup");
    // mf::StartMessageFacility(mf::MessageFacilityService::MultiThread,
    // services_pset.get<fhicl::ParameterSet>("message",
    // fhicl::ParameterSet()));
    mf::LogInfo("MF_INIT_OK") << "Messagelogger initialization complete.";

    // Now create the EventProcessor
    //
    int rc = 1;
            
          
    fhicl::ParameterSet main_pset;
    make_ParameterSet(raw_config, main_pset);
    fhicl::ParameterSetRegistry::put(main_pset);
    
    std::cout << main_pset.to_indented_string() << "\n";

    EventProcessorWithSentry proc;
    try {
      std::unique_ptr<art::EventProcessor>
        procP(new
          art::EventProcessor(main_pset));
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


