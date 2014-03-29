#include <stdlib.h>
#include <stdio.h>

#include <TROOT.h>
#include <TRint.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

#include "art/Persistency/Provenance/Timestamp.h"
#include "RecoBase/Hit.h"


#include "art/Framework/Art/run_art.h"

#include "art/Framework/Art/BasicOptionsHandler.h"
#include "art/Framework/Art/BasicPostProcessor.h"
#include "art/Framework/Art/InitRootHandlers.h"
// #include "art/Framework/EventProcessor/EventProcessor.h"
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

// #include "art/Framework/IO/Root/RootInputFile.h"
// #include "art/Framework/Core/GroupSelectorRules.h"
// #include "art/Persistency/Provenance/ProductMetaData.h"
// #include "art/Persistency/Provenance/ProcessConfiguration.h"
// #include "art/Persistency/Provenance/ModuleDescription.h"
// #include "art/Persistency/Provenance/MasterProductRegistry.h"
// 
// #include "art/Framework/Principal/EventPrincipal.h"
// #include "art/Framework/Principal/Event.h"
// #include "art/Framework/Principal/Handle.h"

#include "art/Framework/Core/InputSourceDescription.h"
#include "art/Framework/Core/InputSourceFactory.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceRegistry.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "art/Persistency/Provenance/MasterProductRegistry.h"
#include "art/Persistency/Provenance/ModuleDescription.h"
#include "art/Persistency/Provenance/ModuleDescriptionID.h"
#include "art/Persistency/Provenance/ProcessConfiguration.h"
// #include "art/Version/GetReleaseVersion.h"
#include "art/Utilities/GetPassID.h"
#include "art/Framework/IO/Root/RootInput.h"

#include "RecoBase/Hit.h"
#include <memory>

int Error; 

extern void InitGui(); // Initializer for GUI
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };

// Initialize the ROOT system
// TROOT root("Rint", "The ROOT Interactive Interface", initfuncs);

// CHEAT
// class MPRGlobalTestFixture
// {
// public:
//   MPRGlobalTestFixture() {
//     art::MasterProductRegistry masterProductRegistry;
//     art::ProductMetaData::create_instance(masterProductRegistry);
//     
//   }
// };


int main(int argc, char **argv)
{
 std::string filename = "standard_reco_uboone.root";
 
  art::RootDictionaryManager dictLoader;
  
   art::RootDictionaryManager rdm;
   art::completeRootHandlers();
   art::ServiceToken dummyToken;

   // Configuration:
  //  cet::filepath_lookup policy("FHICL_FILE_PATH");
  //  fhicl::intermediate_table tbl;
  //  fhicl::parse_document(configName() , policy, tbl);
  // 
  //  fhicl::ParameterSet params;
  // fhicl::make_ParameterSet(tbl, params);
  // 
  //   std::vector<fhicl::ParameterSet> groupParams = params.get<std::vector<fhicl::ParameterSet> >(groupName());
   // fhicl::ParameterSet pset;
   

   // MPRGlobalTestFixture cheat;

   // art::ProcessConfiguration pc;
   // art::EventID origEventID;
   // art::FastCloningInfoProvider fcip;
   // art::GroupSelectorRules groupSelectorRules(pset,"blah_GSR","blah_GSR2");
   // std::shared_ptr<art::DuplicateChecker> dupChecker;
   // art::ModuleDescription moduleDesc;
   // // Attempt to open a file by hand!.
   //
   // 
   // std::shared_ptr<TFile> filePtr(new TFile(filename.c_str(),"READ"));
   // 
   // art::RootInputFile rif(
   //   filename,                    // fileName
   //   "",                          // catalogName
   //   pc,                          // process configuration
   //   filename,                    // logical file name
   //   filePtr,                     // file pointer
   //   origEventID ,                // origEventID
   //   0,                           // eventsToSkip
   //   std::vector<art::SubRunID>(), // which SubRunsToSkip
   //   fcip,                         // FastCloningInfoProvider
   //   0,                           // treeCacheSize
   //   -1,                          // treeMaxVirtualSize
   //   art::InputSource::ProcessingMode::RunsSubRunsAndEvents,  // InputSource::ProcessingMode
   //   0,                           // forcedRunOffset
   //   std::vector<art::EventID>(),  // which events to process
   //   false,                         // noEventSort
   //   groupSelectorRules,          // GroupSelectorRules
   //   false,                       // dropMergable
   //   dupChecker,
   //   false                        // dropDescendants     
   //   );
   //   
   //   rif.skipEvents(2);
   //   std::unique_ptr<art::EventPrincipal> result = rif.readEvent();
   //   art::Event event(*result, moduleDesc);
  
   fhicl::ParameterSet pset;
   std::string processName=argv[0];

   // art::EvProcInitHelper helper_;

   art::MasterProductRegistry preg_;
   art::ServiceToken serviceToken_;
   art::ActivityRegistry actReg_;
   // art::ServiceDirector serviceDirector_(pset,actReg_,serviceToken_);
     
   // Setup the fake fcl file.
   pset.put("processName", processName);

   art::ServiceRegistry::Operate operate(serviceToken_); // Make usable.
   // addSystemServices_(pset); // maybe don't need? .. yess, maybe I do, since serviceToken.forceCreation fails.
   
   // ParameterSet const fpc_pset = helper_.servicesPS().get<ParameterSet>("floating_point_control", ParameterSet());
   
   // serviceDirector_.addSystemService(std::unique_ptr<CurrentModule>(new CurrentModule(actReg_)));
   // // special construction
   // serviceDirector_.addSystemService(std::unique_ptr<TriggerNamesService>
   //                   (new TriggerNamesService(pset, pathManager_.triggerPathNames())));
   // serviceDirector_.addSystemService(std::unique_ptr<FloatingPointControl>(new FloatingPointControl(fpc_pset, actReg_)));
   // serviceDirector_.addSystemService(std::unique_ptr<ScheduleContext>(new ScheduleContext));
   // ParameterSet pathSelection;
   // if (helper_.servicesPS().get_if_present("PathSelection", pathSelection)) {
   //   serviceDirector_.addSystemService(std::unique_ptr<PathSelection>(new PathSelection(*this)));
   // }
   // 
   // 
   serviceToken_.forceCreation();
   // System service FileCatalogMetadata needs to know about the process name.
   art::ServiceHandle<art::FileCatalogMetadata>()->addMetadata("process_name", processName);
   
   fhicl::ParameterSet main_input;
   main_input.put("module_type", "RootInput");
   main_input.put("module_label","source");
   main_input.put("max_events",-1);
   std::vector<std::string> source_list;
   source_list.push_back(filename); // Add the file.
   main_input.put("fileNames", source_list);
  
   art::ModuleDescription md(main_input.id(),
                             "RootInput",
                             "source",
                              art::ProcessConfiguration(processName,
                                                       pset.id(),
                                                       "Whothefuckcares",//art::getReleaseVersion(),
                                                       art::getPassID()));
   art::InputSourceDescription isd(md, preg_, actReg_);

   std::unique_ptr<art::InputSource> source(art::InputSourceFactory::make(main_input,isd).release());

   art::RootInput* ri = dynamic_cast<art::RootInput*>(source.get());
   
   if(ri) std::cout << "Got a RootInput" << std::endl;
   
     // // Get list of hits.
  //    typedef art::Handle< std::vector<recob::Hit> > hitHandle_t;
  //    std::vector<hitHandle_t> list_of_hitlists;
  //    try{
  //      event.getManyByType(list_of_hitlists);
  //    }catch(...) {
  //      std::cout << "problem getting types." << std::endl;
  //    }
  // 
  //    std::cout << "Got " << list_of_hitlists.size() << " types of hits." << std::endl;
  //    for(size_t i=0;i<list_of_hitlists.size(); i++) {
  //      std::cout << "--->" << list_of_hitlists[i].provenance()->moduleLabel() << std::endl;
  //    }
     
  
  // Try.
  // TFile f("single_gen.root");
  // f.ls();
  
  // gSystem->Load("libart_Persistency_Provenance_map.dylib");
  // gSystem->Load("libRecoBase.dylib");
  // gSystem->Load("libRecoBase_dict.dylib");
  // gSystem->Load("libRecoBase_map.dylib");
   
  // art::Timestamp *ts = new art::Timestamp;
  // recob::Hit* hit = new recob::Hit();
  // 
  // TFile f("temp.root","RECREATE");
  // TTree* mytree = new TTree("mytree","mytree");
  // mytree->Branch("hits",&hit,3200,2);
  // mytree->Fill();
  // mytree->Fill();
  // mytree->Fill();
  // mytree->Write();
  // f.Write();
  // f.WriteObject(hit,"recob::Hit","blah");
  // f.WriteObject(ts,"art::Timestamp","blah_ts");
  // f.Write();
  // 
  // f.Close();
  // 
  // 
  // {
  //   TFile infile("/Users/tagg/Argo/server/newpionreco_s1_uboone.root","READ");
  //   if(infile.IsZombie()) {
  //     std::cout << "Couldn't open file." << std::endl;
  //   }
  //   TTree* tree =  (TTree*) infile.Get("Events");
  //   tree->Show(0);
  // }
  // 
  // TRint *theApp = new TRint("ROOT example", &argc, argv, NULL, 0, 0);

//  gStyle->SetPalette(1,0);
  
  // Run interactive interface
  // theApp->Run();
  return(0);
}


/*
#include "art/Framework/Art/artapp.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <iostream>

int main( int argc, char* argv[] ) {
  int result = artapp(argc,argv);
  // Make sure the message logger is shut down so our message really is
  // the last one put out.
  mf::MessageFacilityService::instance().MFPresence.reset();
  // End message.
  std::cout
    << "Art has completed and will exit with status "
    << result
    << "."
    << std::endl;
  return result;
}

// Local Variables:
// mode: c++
// End:
*/