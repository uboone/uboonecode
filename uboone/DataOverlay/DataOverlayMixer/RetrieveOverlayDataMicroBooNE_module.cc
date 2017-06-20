////////////////////////////////////////////////////////////////////////
// Class:       RetrieveOverlayDataMicroBooNE
// Module Type: producer
// File:        RetrieveOverlayDataMicroBooNE_module.cc
//
// This borrows a lot from the Mu2e mixing module:
//      EventMixing/src/MixMCEvents_module.cc
//
// The purpose of this module is simply to grab raw data products from the 
// secondary stream for storage.  Mixing will occur downstream once 
// the relevant primary stream raw data products are available
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Core/PtrRemapper.h"
#include "art/Persistency/Common/CollectionUtilities.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "IFDH_service.h"

#include <memory>
#include <string>
#include <vector>
#include <exception>
#include <sstream>
#include <unistd.h>

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"

#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCShower.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/TriggerData.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "lardataobj/RawData/OpDetWaveform.h"

#include "DataOverlayProducts/EventMixingSummary.h"

namespace mix {
  class RetrieveOverlayDataMicroBooNE;
}

class mix::RetrieveOverlayDataMicroBooNE : public boost::noncopyable {
public:

  RetrieveOverlayDataMicroBooNE(fhicl::ParameterSet const& p,
				 art::MixHelper &helper);
  ~RetrieveOverlayDataMicroBooNE();

  void startEvent(const art::Event&);  //called at the start of every event
  void finalizeEvent(art::Event &);    //called at the end of every event
  
  size_t nSecondaries() { return fEventsToMix; } 

  void processEventAuxiliaries(art::EventAuxiliarySequence const& seq); //bookkepping for event IDs

  // Mixing Functions

  // For now, allow exactly one  input. Assume MC inputs have been merged
  // previously and one detsim output created if needed. This could be changed
  // but would require mixing functions for MC here.

  //a lot of MC collections are just simple copies of the collections...
  template<typename T>
  bool MixSimpleCopy( std::vector< std::vector<T> const*> const& inputs,
		      std::vector< T > & output,
		      art::PtrRemapper const &);

  // Choose mix file.
  std::string getMixFile();
		   
  
private:

  // Declare member data here.
  
  fhicl::ParameterSet  fpset;
  bool                 fInputFileIsData;

  std::string          fRawDigitDataModuleLabel;
  std::string          fOpDetDataModuleLabel;
  std::string          fTriggerDataModuleLabel;
  std::string          fRawDigitMCModuleLabel;
  std::string          fOpDetMCModuleLabel;
  std::string          fTriggerMCModuleLabel;

  std::string          fRawDigitMixerSourceModuleLabel;
  std::string          fOpDetMixerSourceModuleLabel;
  std::string          fTriggerMixerSourceModuleLabel;

  std::string          fG4InputModuleLabel;
  std::string          fGeneratorInputModuleLabel;

  bool                 fDoMCReco;
  std::string          fMCRecoInputModuleLabel;

  size_t               fEventsToMix;

  std::string          fSamDefname;
  std::string          fSamProject;
  std::string          fSamStation;
  std::string          fSamAppFamily;
  std::string          fSamAppName;
  std::string          fSamAppVersion;
  std::string          fSamUser;
  std::string          fSamDescription;
  int                  fSamFileLimit;
  std::string          fSamSchema;

  std::string          fSamProjectURI;
  std::string          fSamProcessID;
  std::string          fSamCurrentFileURI;
  std::string          fSamCurrentFileName;

  std::unique_ptr< std::vector<mix::EventMixingSummary> > fEventMixingSummary;
  
};


mix::RetrieveOverlayDataMicroBooNE::RetrieveOverlayDataMicroBooNE(fhicl::ParameterSet const& p,
								    art::MixHelper &helper)
  :
  fpset(p.get<fhicl::ParameterSet>("detail")),
  fInputFileIsData(fpset.get<bool>("InputFileIsData")),
  fRawDigitDataModuleLabel(fpset.get<std::string>("RawDigitDataModuleLabel")),
  fOpDetDataModuleLabel(fpset.get<std::string>("OpDetDataModuleLabel")),
  fTriggerDataModuleLabel(fpset.get<std::string>("TriggerDataModuleLabel")),
  fRawDigitMCModuleLabel(fpset.get<std::string>("RawDigitMCModuleLabel")),
  fOpDetMCModuleLabel(fpset.get<std::string>("OpDetMCModuleLabel")),
  fTriggerMCModuleLabel(fpset.get<std::string>("TriggerMCModuleLabel")),
  fEventsToMix(fpset.get<size_t>("EventsToMix",1)),


  // Get sam related parameters.
  // These parameters should normally be set by the work flow.
  // Usually, the only ones that should need to be set are "SamDefname" and "SamProject."

  fSamDefname(fpset.get<std::string>("SamDefname", "")),
  fSamProject(fpset.get<std::string>("SamProject", "")),
  fSamStation(fpset.get<std::string>("SamStation", "")),
  fSamAppFamily(fpset.get<std::string>("SamAppFamily", "art")),
  fSamAppName(fpset.get<std::string>("SamAppName", "retrieve")),
  fSamAppVersion(fpset.get<std::string>("SamAppVersion", "1")),
  fSamUser(fpset.get<std::string>("SamUser", "")),
  fSamDescription(fpset.get<std::string>("SamDescription", "")),
  fSamFileLimit(fpset.get<int>("SamFileLimit", 100)),
  fSamSchema(fpset.get<std::string>("SamSchema", "root")),    // xrootd by default.

  fEventMixingSummary(nullptr)
{
  
  if(fEventsToMix!=1){
    std::stringstream err_str;
    err_str << "ERROR! Really sorry, but we can only do mixing for one collection right now! ";
    err_str << "\nYep. We're gonna throw an exception now. You should change your fcl to set 'EventsToMix' to 1";
    throw cet::exception("OverlayRawDataMicroBooNE") << err_str.str() << std::endl;;
  }

  if(fInputFileIsData){
    fRawDigitMixerSourceModuleLabel = fRawDigitMCModuleLabel;
    fOpDetMixerSourceModuleLabel    = fOpDetMCModuleLabel;
    fTriggerMixerSourceModuleLabel  = fTriggerMCModuleLabel;
  }
  else if(!fInputFileIsData){
    fRawDigitMixerSourceModuleLabel = fRawDigitDataModuleLabel;
    fOpDetMixerSourceModuleLabel    = fOpDetDataModuleLabel;
    fTriggerMixerSourceModuleLabel  = fTriggerDataModuleLabel;
  }
  
  if(fInputFileIsData){
    fDoMCReco = fpset.get_if_present<std::string>("MCRecoInputModuleLabel",fMCRecoInputModuleLabel);
    fG4InputModuleLabel = fpset.get<std::string>("G4InputModuleLabel");
    fGeneratorInputModuleLabel = fpset.get<std::string>("GeneratorInputModuleLabel");
    
    //MC generator info is a simple copy
    helper.declareMixOp( art::InputTag(fGeneratorInputModuleLabel),
			 &RetrieveOverlayDataMicroBooNE::MixSimpleCopy<simb::MCTruth>,
			 *this );
    
    //Simple copies of G4 SimPhotons, MCParticles, SimChannels, and SimAuxDetChannel
    helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
			 &RetrieveOverlayDataMicroBooNE::MixSimpleCopy<simb::MCParticle>,
			 *this );
    helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
			 &RetrieveOverlayDataMicroBooNE::MixSimpleCopy<sim::SimPhotons>,
			 *this );
    helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
			 &RetrieveOverlayDataMicroBooNE::MixSimpleCopy<sim::SimChannel>,
			 *this );
    helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
			 &RetrieveOverlayDataMicroBooNE::MixSimpleCopy<sim::AuxDetSimChannel>,
			 *this );
    /*
    //Associations of MCParticles to MCTruth...hopefully a simple copy is enough
    helper.declareMixOp( art::InputTag(fG4InputModuleLabel),
    &RetrieveOverlayDataMicroBooNE::MixSimpleCopy
		       < art::Assns<simb::MCTruth,simb::MCParticle,void> >,
		       *this );
		       */
    
    //Copies of MCShower and MCTrack
    if(fDoMCReco){
      helper.declareMixOp( art::InputTag(fMCRecoInputModuleLabel),
			 &RetrieveOverlayDataMicroBooNE::MixSimpleCopy<sim::MCShower>,
			   *this );
      helper.declareMixOp( art::InputTag(fMCRecoInputModuleLabel),
			   &RetrieveOverlayDataMicroBooNE::MixSimpleCopy<sim::MCTrack>,
			   *this );
    }
  }//end if file is input data
  else if(!fInputFileIsData){
    helper.declareMixOp( art::InputTag(fOpDetMixerSourceModuleLabel,"OpdetCosmicHighGain"),
			 &RetrieveOverlayDataMicroBooNE::MixSimpleCopy<raw::OpDetWaveform>,
			 *this );
  }//end if file is input mc
  
  //Simply copy the RawDigits and OpDetWaveforms from the Mixer Source
  helper.declareMixOp( art::InputTag(fRawDigitMixerSourceModuleLabel),
		       &RetrieveOverlayDataMicroBooNE::MixSimpleCopy<raw::RawDigit>,
		       *this );
  
  helper.declareMixOp( art::InputTag(fOpDetMixerSourceModuleLabel,"OpdetBeamHighGain"),
		       &RetrieveOverlayDataMicroBooNE::MixSimpleCopy<raw::OpDetWaveform>,
		       *this );
  helper.declareMixOp( art::InputTag(fOpDetMixerSourceModuleLabel,"OpdetBeamLowGain"),
		       &RetrieveOverlayDataMicroBooNE::MixSimpleCopy<raw::OpDetWaveform>,
		       *this );

  helper.declareMixOp( art::InputTag(fTriggerMixerSourceModuleLabel),
		       &RetrieveOverlayDataMicroBooNE::MixSimpleCopy<raw::Trigger>,
		       *this );

  //If it produces something on its own, declare it here
  helper.produces< std::vector<mix::EventMixingSummary> >();

  // Following block handles case of mix input from sam.

  if(!fSamDefname.empty()) {

    // Register getMixFile method with MixHelper.

    helper.registerSecondaryFileNameProvider(std::bind(&mix::RetrieveOverlayDataMicroBooNE::getMixFile, this));

    // Get IFDH art service.

    art::ServiceHandle<ifdh_ns::IFDH> ifdh;

    // Get sam station.
    // If the station was not specified by a fcl parameter, use environment variable
    // $SAM_STATION, or else use a default value of "uboone."

    if(fSamStation.empty()) {
      const char* c = getenv("SAM_STATION");
      if(c == 0 || *c == 0)
	c = "uboone";
      fSamStation = c;
      //std::cout << "Mix SAM: Station = " << fSamStation << std::endl;
    }

    // Find project uri.

    fSamProjectURI = ifdh->findProject(fSamProject, fSamStation);
    //std::cout << "Mix SAM: project uri = " << fSamProjectURI << std::endl;
    if(fSamProjectURI.empty())
      throw cet::exception("OverlayRawDataMicroBooNE") << "Failed to find project uri.";

    // Get hostname.

    char hostname[256];
    gethostname(hostname, sizeof hostname);

    // Get user.
    // If the user was not specified by a fcl parameter, use environment variable
    // $SAM_USER (this should work on grid), or else use environment variable $LOGNAME.

    if(fSamUser.empty()) {
      const char* c = getenv("SAM_USER");
      if(c == 0 || *c == 0)
	c = getenv("LOGNAME");
      if(c != 0 && *c != 0)
	fSamUser = c;
      //std::cout << "Mix SAM: User = " << fSamUser << std::endl;
    }

    // Join project.

    fSamProcessID = ifdh->establishProcess(fSamProjectURI,
					   fSamAppName,
					   fSamAppVersion,
					   hostname,
					   fSamUser,
					   fSamAppFamily,
					   fSamDescription,
					   fSamFileLimit,
					   fSamSchema);
    mf::LogInfo("OverlayRawDigitMicroBooNE") << "Overlay sam definition: " << fSamDefname << "\n"
					     << "Overlay sam project: " << fSamProject;

    //std::cout << "Mix SAM: process id = " << fSamProcessID << std::endl;
    if(fSamProcessID.empty())
      throw cet::exception("OverlayRawDataMicroBooNE") << "Failed to start sam process.";
  }
}

// Destructor.
mix::RetrieveOverlayDataMicroBooNE::~RetrieveOverlayDataMicroBooNE()
{
  if(!fSamProcessID.empty()) {

    // Get IFDH art service.

    art::ServiceHandle<ifdh_ns::IFDH> ifdh;

    // Mark current file as consumed.

    if(!fSamCurrentFileName.empty()) {
      ifdh->updateFileStatus(fSamProjectURI,
			     fSamProcessID,
			     fSamCurrentFileName,
			     "consumed");
      //std::cout << "Mix SAM: File " << fSamCurrentFileName << " status changed to consumed." << std::endl;
    }

    // Stop process.

    ifdh->endProcess(fSamProjectURI, fSamProcessID);
    //std::cout << "Mix SAM: End process." << std::endl;
  }
}

//Initialize for each event
void mix::RetrieveOverlayDataMicroBooNE::startEvent(const art::Event& event) {

  if(!( (event.isRealData() && fInputFileIsData) || (!event.isRealData() && !fInputFileIsData)))
    throw cet::exception("OverlayRawDataMicroBooNE") << "Input file claimed to be data/not data, but it's not." << std::endl;;

  fEventMixingSummary.reset(new std::vector<mix::EventMixingSummary>);
}

//For each of the mixed in events...bookkepping for event IDs
void mix::RetrieveOverlayDataMicroBooNE::processEventAuxiliaries(art::EventAuxiliarySequence const & seq){
  for (auto const& ev : seq)
    fEventMixingSummary->emplace_back(ev.id().event(),ev.id().subRun(),ev.id().run(),ev.time());
}

//End each event
void mix::RetrieveOverlayDataMicroBooNE::finalizeEvent(art::Event& event) {
  event.put(std::move(fEventMixingSummary));
}


template<typename T>
bool mix::RetrieveOverlayDataMicroBooNE::MixSimpleCopy( std::vector< std::vector<T> const*> const& inputs,
							 std::vector< T > & output,
							 art::PtrRemapper const &){
  art::flattenCollections(inputs,output);
  return true;
}


// Return next file to mix.
std::string mix::RetrieveOverlayDataMicroBooNE::getMixFile()
{
  std::string result;

  if(!fSamProcessID.empty()) {

    // Get IFDH art service.

    art::ServiceHandle<ifdh_ns::IFDH> ifdh;

    // Update status of current file, if any, to "consumed."

    if(!fSamCurrentFileName.empty()) {
      ifdh->updateFileStatus(fSamProjectURI,
			     fSamProcessID,
			     fSamCurrentFileName,
			     "consumed");

      //std::cout << "Mix SAM: File " << fSamCurrentFileName << " status changed to consumed." << std::endl;
      fSamCurrentFileURI = std::string();
      fSamCurrentFileName = std::string();
    }

    // Get next file uri.

    fSamCurrentFileURI = fSamCurrentFileURI = ifdh->getNextFile(fSamProjectURI,	fSamProcessID);
    unsigned int n = fSamCurrentFileURI.find_last_of('/') + 1;
    fSamCurrentFileName = fSamCurrentFileURI.substr(n);
    //std::cout << "Mix SAM: Next file uri = " << fSamCurrentFileURI << std::endl;
    //std::cout << "Mix SAM: Next file name = " << fSamCurrentFileName << std::endl;
    mf::LogInfo("OverlayRawDigitMicroBooNE") << "Next mix file uri: " << fSamCurrentFileURI << "\n"
					     << "Next mix file name: " << fSamCurrentFileName;

    // Throw an exception if we didn't get a next file.

    if(fSamCurrentFileURI.empty() || fSamCurrentFileName.empty())
      throw cet::exception("OverlayRawDataMicroBooNE") << "Failed to get next mix file.";

    // Here is where we would copy the file to the local node, if that were necessary.
    // Since we are using schema "root" (i.e. xrootd) to stream files, copying the 
    // file is not necessary.
    // Note further that we should not update the file status to "transferred" for
    // streaming files, since that can in principle allow the file to be deleted from
    // disk cache.

    // Update metadata.

    art::ServiceHandle<art::FileCatalogMetadata> md;
    md->addMetadataString("mixparent", fSamCurrentFileName);

    // Done.

    result = fSamCurrentFileURI;
  }

  return result;
}


DEFINE_ART_MODULE(art::MixFilter<mix::RetrieveOverlayDataMicroBooNE>)
