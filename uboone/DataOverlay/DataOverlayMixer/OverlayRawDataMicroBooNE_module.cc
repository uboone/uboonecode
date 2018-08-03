////////////////////////////////////////////////////////////////////////
// Class:       OverlayRawDataDetailMicroBooNE
// Module Type: producer
// File:        OverlayRawDataDetailMicroBooNE_module.cc
//
// This borrows a lot from the Mu2e mixing module:
//      EventMixing/src/MixMCEvents_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Persistency/Common/CollectionUtilities.h"
#include "art/Framework/Services/System/FileCatalogMetadata.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <exception>
#include <sstream>
#include <unistd.h>

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"
#include "larevt/CalibrationDBI/Interface/ChannelStatusProvider.h"

#include "DataOverlay/RawDigitMixer.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/TriggerData.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "DataOverlay/OpDetWaveformMixer.h"
#include "lardataobj/RawData/OpDetWaveform.h"

namespace mix {

class OverlayRawDataMicroBooNE : public art::EDProducer {
  public:

    OverlayRawDataMicroBooNE(fhicl::ParameterSet const& p);
    ~OverlayRawDataMicroBooNE() {}

    void produce (art::Event& evt);

  private:

    // Declare member data here.
    RawDigitMixer              fRDMixer;
    OpDetWaveformMixer         fODMixer;

    short                fDefaultRawDigitSatPoint;
    short                fDefaultOpDetSatPoint;
    size_t               fOpDetMinSampleSize;

    std::string          fRawDigitDataModuleLabel;
    std::string          fOpDetDataModuleLabel;
    std::string          fTriggerDataModuleLabel;
    std::string          fRawDigitMCModuleLabel;
    std::string          fOpDetMCModuleLabel;
    std::string          fTriggerMCModuleLabel;

    float                fDefaultMCRawDigitScale;
    float                fDefaultMCOpDetScale;

    void GenerateMCRawDigitScaleMap(std::vector<raw::RawDigit> const&);
    std::unordered_map<raw::ChannelID_t,float> fMCRawDigitScaleMap;

    void GenerateMCOpDetHighGainScaleMap(std::vector<raw::OpDetWaveform> const&);
    std::unordered_map<raw::Channel_t,float> fMCOpDetHighGainScaleMap;

    void GenerateMCOpDetLowGainScaleMap(std::vector<raw::OpDetWaveform> const&);
    std::unordered_map<raw::Channel_t,float> fMCOpDetLowGainScaleMap;

    bool MixRawDigits( const art::Event& evt, std::vector<raw::RawDigit> & output);
    
    bool MixTriggerData( const art::Event& evt, std::vector<raw::Trigger> & output);

    bool MixOpDetWaveforms_HighGain( const art::Event& evt, std::vector<raw::OpDetWaveform> & output);
    bool MixOpDetWaveforms_LowGain( const art::Event& evt, std::vector<raw::OpDetWaveform> & output);
  };
}//end namespace mix


mix::OverlayRawDataMicroBooNE::OverlayRawDataMicroBooNE(fhicl::ParameterSet const& p)
  :
  fRDMixer(false), //print warnings turned off
  fODMixer(false), //print warnings turned off

  fDefaultRawDigitSatPoint(p.get<short>("DefaultRawDigitSaturationPoint",4096)),
  fDefaultOpDetSatPoint(p.get<short>("DefaultOpDetSaturationPoint",4096)),
  fOpDetMinSampleSize(p.get<size_t>("OpDetMinSampleSize",100)),
  fRawDigitDataModuleLabel(p.get<std::string>("RawDigitDataModuleLabel")),
  fOpDetDataModuleLabel(p.get<std::string>("OpDetDataModuleLabel")),
  fTriggerDataModuleLabel(p.get<std::string>("TriggerDataModuleLabel")),
  fRawDigitMCModuleLabel(p.get<std::string>("RawDigitMCModuleLabel")),
  fOpDetMCModuleLabel(p.get<std::string>("OpDetMCModuleLabel")),
  fTriggerMCModuleLabel(p.get<std::string>("TriggerMCModuleLabel")),
  fDefaultMCRawDigitScale(p.get<float>("DefaultMCRawDigitScale",1)),
  fDefaultMCOpDetScale(p.get<float>("DefaultMCOpDetScale",1))
{
  fRDMixer.SetSaturationPoint(fDefaultRawDigitSatPoint);
  fODMixer.SetSaturationPoint(fDefaultOpDetSatPoint);
  fODMixer.SetMinSampleSize(fOpDetMinSampleSize);
  
  produces< std::vector<raw::RawDigit> >();
  produces< std::vector<raw::OpDetWaveform> >("OpdetBeamHighGain");
  produces< std::vector<raw::OpDetWaveform> >("OpdetBeamLowGain");
  produces< std::vector<raw::Trigger> >();
}

void mix::OverlayRawDataMicroBooNE::produce(art::Event& evt) {

  //make output containers 
  std::unique_ptr<std::vector<raw::RawDigit> >     rawdigits(new std::vector<raw::RawDigit>);
  std::unique_ptr<std::vector<raw::OpDetWaveform> > opdet_hg(new std::vector<raw::OpDetWaveform>);
  std::unique_ptr<std::vector<raw::OpDetWaveform> > opdet_lg(new std::vector<raw::OpDetWaveform>);
  std::unique_ptr<std::vector<raw::Trigger> > triggerdata(new std::vector<raw::Trigger>);
  
  //get output digits
  MixRawDigits(evt, *rawdigits);
  MixOpDetWaveforms_HighGain(evt, *opdet_hg);
  MixOpDetWaveforms_LowGain(evt, *opdet_lg);
  MixTriggerData(evt, *triggerdata);

  
  //put output digits
  evt.put(std::move(rawdigits));
  evt.put(std::move(opdet_hg),"OpdetBeamHighGain");
  evt.put(std::move(opdet_lg),"OpdetBeamLowGain");
  evt.put(std::move(triggerdata));
}


void mix::OverlayRawDataMicroBooNE::GenerateMCRawDigitScaleMap(std::vector<raw::RawDigit> const& dataDigitVector){
  //right now, assume the number of channels is the number in the collection
  //and, loop through the channels one by one to get the right channel number
  //note: we will put here access to the channel database to determine dead channels
  fMCRawDigitScaleMap.clear();

  const lariov::ChannelStatusProvider& chanStatus = art::ServiceHandle<lariov::ChannelStatusService>()->GetProvider();
  
  for(auto const& d : dataDigitVector){
    if(chanStatus.IsBad(d.Channel()))
      fMCRawDigitScaleMap[d.Channel()] = 0.0;
    else
      fMCRawDigitScaleMap[d.Channel()] = fDefaultMCRawDigitScale;
  }
}

bool mix::OverlayRawDataMicroBooNE::MixRawDigits( const art::Event& event, std::vector<raw::RawDigit> & output) {
  
  output.clear();
  
  art::Handle< std::vector<raw::RawDigit> > mcDigitHandle;
  event.getByLabel( fRawDigitMCModuleLabel,mcDigitHandle);
  
  art::Handle< std::vector<raw::RawDigit> > dataDigitHandle;
  event.getByLabel( fRawDigitDataModuleLabel,dataDigitHandle);

  GenerateMCRawDigitScaleMap(*dataDigitHandle);
  fRDMixer.DeclareData(*dataDigitHandle);
  fRDMixer.Mix(*mcDigitHandle,fMCRawDigitScaleMap);
  fRDMixer.FillRawDigitOutput(output);
  
  return true;
}

void mix::OverlayRawDataMicroBooNE::GenerateMCOpDetHighGainScaleMap(std::vector<raw::OpDetWaveform> const& dataVector){
  //right now, assume the number of channels is the number in the collection
  //and, loop through the channels one by one to get the right channel number
  //note: we will put here access to the channel database to determine dead channels
  fMCOpDetHighGainScaleMap.clear();
  for(auto const& d : dataVector)
    fMCOpDetHighGainScaleMap[d.ChannelNumber()] = fDefaultMCOpDetScale;
}

void mix::OverlayRawDataMicroBooNE::GenerateMCOpDetLowGainScaleMap(std::vector<raw::OpDetWaveform> const& dataVector){
  //right now, assume the number of channels is the number in the collection
  //and, loop through the channels one by one to get the right channel number
  //note: we will put here access to the channel database to determine dead channels
  fMCOpDetLowGainScaleMap.clear();
  for(auto const& d : dataVector)
    fMCOpDetLowGainScaleMap[d.ChannelNumber()] = fDefaultMCOpDetScale;
}

bool mix::OverlayRawDataMicroBooNE::MixOpDetWaveforms_HighGain( const art::Event& event, std::vector<raw::OpDetWaveform> & output) {
  
  output.clear();
  
  art::Handle< std::vector<raw::OpDetWaveform> > mcOpDetHandle_HighGain;
  event.getByLabel(fOpDetMCModuleLabel,"OpdetBeamHighGain",mcOpDetHandle_HighGain);
  
  art::Handle< std::vector<raw::OpDetWaveform> > dataOpDetHandle_HighGain;
  event.getByLabel(fOpDetDataModuleLabel,"OpdetBeamHighGain",dataOpDetHandle_HighGain);  

  GenerateMCOpDetHighGainScaleMap(*dataOpDetHandle_HighGain); 
  fODMixer.DeclareData(*dataOpDetHandle_HighGain,output);
  fODMixer.Mix(*mcOpDetHandle_HighGain, fMCOpDetHighGainScaleMap, output);
  
  return true;
}

bool mix::OverlayRawDataMicroBooNE::MixOpDetWaveforms_LowGain( const art::Event& event, std::vector<raw::OpDetWaveform> & output) {

  output.clear();

  art::Handle< std::vector<raw::OpDetWaveform> > mcOpDetHandle_LowGain;
  event.getByLabel(fOpDetMCModuleLabel,"OpdetBeamLowGain",mcOpDetHandle_LowGain);
  
  art::Handle< std::vector<raw::OpDetWaveform> > dataOpDetHandle_LowGain;
  event.getByLabel(fOpDetDataModuleLabel,"OpdetBeamLowGain",dataOpDetHandle_LowGain);  

  GenerateMCOpDetLowGainScaleMap(*dataOpDetHandle_LowGain); 
  fODMixer.DeclareData(*dataOpDetHandle_LowGain,output);
  fODMixer.Mix(*mcOpDetHandle_LowGain, fMCOpDetLowGainScaleMap, output);
  
  return true;
}

bool mix::OverlayRawDataMicroBooNE::MixTriggerData( const art::Event& event, std::vector<raw::Trigger> & output) {
  
  output.clear();
  
  art::Handle< std::vector<raw::Trigger> > mcTriggerHandle;
  event.getByLabel(fTriggerMCModuleLabel, mcTriggerHandle);
  
  art::Handle< std::vector<raw::Trigger> > dataTriggerHandle;
  event.getByLabel(fTriggerDataModuleLabel, dataTriggerHandle);
  
  unsigned int trig_num; double trig_time,gate_time; unsigned int trig_bits; 
  trig_num  = dataTriggerHandle->at(0).TriggerNumber();
  trig_time = dataTriggerHandle->at(0).TriggerTime();
  gate_time = dataTriggerHandle->at(0).BeamGateTime();
  trig_bits = dataTriggerHandle->at(0).TriggerBits() | mcTriggerHandle->at(0).TriggerBits();
  output.emplace_back(trig_num,trig_time,gate_time,trig_bits);
  
  return true;
}


DEFINE_ART_MODULE(mix::OverlayRawDataMicroBooNE)
