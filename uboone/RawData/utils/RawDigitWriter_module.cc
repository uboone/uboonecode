/**
 * @file   RawDigitWriter_module.cc
 * @brief  Saves RawDigits in a flat ROOT tree. For quick-access to digits for diagnostics.
 * @author Taritree Wongjirad (taritree@mit.edu)
 * @date   June 13th, 2015
 */

// C//C++ standard libraries
#include <string>
#include <vector>
#include <algorithm> // std::min(), std::copy_n()
#include <ios> // std::fixed
#include <iomanip> // std::setprecision(), std::setw()
#include <memory> // std::unique_ptr<>
#include <cstdlib> // FILE
#include <stdio.h> // remove

// support libraries
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// art libraries
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LArSoft includes
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::View_t
#include "larevt/Filters/ChannelFilter.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
// RawDigits
#include "lardataobj/RawData/raw.h" // raw::Uncompress()
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RawData/TriggerData.h"
// Optical Channel Maps
#include "uboone/Geometry/UBOpChannelTypes.h"
#include "uboone/Geometry/UBOpReadoutMap.h"
// TPC Channel Map
#include "lardata/Utilities/DatabaseUtil.h"

// ROOT
#include "TTree.h"

namespace zmqds {

  /**
   * @brief Prints the content of all the raw digits on screen
   *
   * This analyser prints the content of all the raw digits into the
   * LogInfo/LogVerbatim stream.
   * 
   * <b>Configuration parameters</b>
   * 
   * - <b>DetSimModuleLabel</b> (string, default: "daq"): label of the
   *   producer used to create the raw::RawDigits collection
   * - <b>OutputCategory</b> (string, default: "DumpDigits"): the category used
   *   for the output (useful for filtering)
   * - <b>DigitsPerLine</b> (integer, default: 20): the dump of digits and ticks
   *   will put this many of them for each line
   * - <b>IgnoreFilters</b> (boolean, default: false): if true, channel filters
   *   will be ignored; by default, only raw digits on channels that are not bad
   *   are printed out
   * - <b>Pedestal</b> (integer, default: 0): digit values are written relative
   *   to this number
   *
   */

  class RawDigitWriter: public art::EDAnalyzer {
  public:
    
    /// Default constructor
    explicit RawDigitWriter(fhicl::ParameterSet const& pset); 
    ~RawDigitWriter();
    
    /// Does the printing
    void analyze (const art::Event& evt);
    
  private:

    // Channel Map
    const util::UBChannelReverseMap_t& fChannelReverseMap;

    // FICHL Parameters

    // Trigger Data
    double fTrigTimestamp;
    double fBeamGateTime;
    void GetTrigger( const art::Event& evt);    
    
    // TPC Raw Digits
    int fRun;
    int fSubRun;
    int fEvent;
    int fCrate;
    int fSlot;
    int fFemCH;
    int fWireID;
    std::vector< short > rawdigits;
    TBranch* fb_rawdigits;
    TTree* fTRawDigits;
    void GetRawDigits(const art::Event& evt);
    
    // Optical Raw Digits
    int fCategory;
    int fType;
    int fFrame;
    int fSample;
    double fTimeStamp;
    double fTrigTimeStamp;
    double fBeamTimeStamp;
    int fOpCrate;
    int fOpSlot;
    int fOpFemCH;
    int fOpReadoutCH;
    TTree* fTOpDetWaveforms;
    std::vector< short > opdetwaveforms;
    void GetOpticalRawDigits(const art::Event& evt);

    bool fWriteTPCdata;
    bool fWritePMTdata;

    // Index tree variables
    TTree* fTRawDigitIndex;
    TTree* fTOpDigitIndex;
    int idx_raw_start;
    int idx_raw_nentries;
    int idx_op_start;
    int idx_op_nentries;
    int idx_raw;
    int idx_op;


  }; // class RawDigitWriter
  
} // namespace zmqds


namespace zmqds {

  //-------------------------------------------------
  RawDigitWriter::RawDigitWriter(fhicl::ParameterSet const& pset) 
    : EDAnalyzer(pset), fChannelReverseMap( art::ServiceHandle<util::DatabaseUtil>()->GetUBChannelReverseMap() )
    {

      art::ServiceHandle<art::TFileService> tfs;
      //art::TFileDirectory tfbeamdir = tfs->mkdir( "RawData" );
      fTRawDigits = tfs->make<TTree>("RawDigits","TPC Waveform data");
      fTRawDigits->Branch( "run", &fRun, "run/I" );
      fTRawDigits->Branch( "subrun", &fSubRun, "subrun/I" );
      fTRawDigits->Branch( "event", &fEvent, "event/I" );
      fTRawDigits->Branch( "crate", &fCrate, "crate/I" );
      fTRawDigits->Branch( "slot", &fSlot, "slot/I" );
      fTRawDigits->Branch( "channel", &fFemCH, "channel/I" );
      fTRawDigits->Branch( "wireid", &fWireID, "wireid/I" );
      fTRawDigits->Branch( "adcs", &rawdigits );

      fTOpDetWaveforms = tfs->make<TTree>("OpDetWaveforms","PMT Readout Waveforms");
      fTOpDetWaveforms->Branch( "run", &fRun, "run/I" );
      fTOpDetWaveforms->Branch( "subrun", &fSubRun, "subrun/I" );
      fTOpDetWaveforms->Branch( "event", &fEvent, "event/I" );
      fTOpDetWaveforms->Branch( "opcrate", &fOpCrate, "opcrate/I" );
      fTOpDetWaveforms->Branch( "opslot", &fOpSlot, "opslot/I" );
      fTOpDetWaveforms->Branch( "opfemch", &fOpFemCH, "opfemch/I" );
      fTOpDetWaveforms->Branch( "frame", &fFrame, "frame/I" );
      fTOpDetWaveforms->Branch( "sample", &fSample, "sample/I" );
      fTOpDetWaveforms->Branch( "readoutch", &fOpReadoutCH, "readoutch/I" );
      fTOpDetWaveforms->Branch( "category", &fCategory, "category/I" );
      fTOpDetWaveforms->Branch( "gaintype", &fType, "gaintype/I" );
      fTOpDetWaveforms->Branch( "timestamp", &fTimeStamp, "timestamp/D" );
      fTOpDetWaveforms->Branch( "trig_timestamp", &fTrigTimeStamp, "trig_timestamp/D" );
      fTOpDetWaveforms->Branch( "beam_timestamp", &fBeamTimeStamp, "beam_timestamp/D" );
      fTOpDetWaveforms->Branch( "adcs", &opdetwaveforms );

      fTRawDigitIndex = tfs->make<TTree>("IndexRawDigits","Event to entry indices for RawDigits");
      fTRawDigitIndex->Branch( "idx_run", &fRun, "idx_run/I");
      fTRawDigitIndex->Branch( "idx_subrun", &fSubRun, "idx_subrun/I");
      fTRawDigitIndex->Branch( "idx_event", &fEvent, "idx_event/I");
      fTRawDigitIndex->Branch( "entrystart", &idx_raw_start, "   idx_start/I");
      fTRawDigitIndex->Branch( "nentries",   &idx_raw_nentries, "idx_nentries/I");
      fTOpDigitIndex  = tfs->make<TTree>("IndexOpDetWfms","Event to entry indices for RawDigits");
      fTOpDigitIndex->Branch( "idx_run", &fRun, "idx_run/I");
      fTOpDigitIndex->Branch( "idx_subrun", &fSubRun, "idx_subrun/I");
      fTOpDigitIndex->Branch( "idx_event", &fEvent, "idx_event/I");
      fTOpDigitIndex->Branch( "entrystart", &idx_op_start, "   idx_start/I");
      fTOpDigitIndex->Branch( "nentries",   &idx_op_nentries, "idx_nentries/I");
      idx_raw = 0; // set these to zero, counter of entries filled into rawdigits tree
      idx_op  = 0; // set these to zero, counter of entries filled into opdet tree      

      mf::LogInfo("")<<"Fetched channel map from DB";

      fWriteTPCdata = pset.get< bool >( "WriteTPCdata", true );
      fWritePMTdata = pset.get< bool >( "WritePMTdata", true );
    }

  RawDigitWriter::~RawDigitWriter() {
    
  }


  //-------------------------------------------------
  void RawDigitWriter::analyze(const art::Event& evt) {

    fRun = (int)evt.run();
    fSubRun = (int)evt.subRun();
    fEvent = (int)evt.event();

    GetTrigger(evt);

    if ( fWriteTPCdata )
      GetRawDigits(evt);

    if ( fWritePMTdata )
      GetOpticalRawDigits( evt );

    // if ( fSendOpticalRawDigits )
    //   SendOpticalRawDigits(evt);
    // else
    //   std::cout << "Skipping OpticalRawDigits..." << std::endl;

    // if ( m_data_files.size()>=(unsigned int)fEventsPerMessage ) {
    //   TransferData();
    // }
    
  } // RawDigitWriter::analyze()

  // =============================================================== 
  // RawDigits

  void RawDigitWriter::GetTrigger( const art::Event& evt) {
   
    //auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    art::Handle<std::vector<raw::Trigger> > trig_handle;

    //evt.getByLabel(fTrigModuleName, trig_handle);
    evt.getByLabel("daq", trig_handle);

    if(!trig_handle.isValid()) {
      // Trigger simulation has not run yet!
      std::cout << "No Trigger Data!!!" << std::endl;
      return;
    }

    if ( trig_handle->size()>0 ) {
      art::Ptr<raw::Trigger> trigVec0( trig_handle, 0 );
      fTrigTimestamp = trigVec0->TriggerTime();
      fBeamGateTime = trigVec0->BeamGateTime();
      //timeService->SetTriggerTime( trigVec0->TriggerTime(), trigVec0->BeamGateTime() );
      //std::cout << "Retrieved Trigger: Trigger Timestamp=" << fTrigTime << " Beam Window Timestamp=" << fBeamGateTime << std::endl;
      //timeService->debugReport();
    }

    if(trig_handle->size()>1)
      throw cet::exception("TimeService::preProcessEvent")
	<< "Found " << trig_handle->size() << " triggers (only 1 trigger/event supported)\n";
  }

  // =============================================================== 
  // RawDigits
  void RawDigitWriter::GetRawDigits(const art::Event& evt) {

    // fetch the data
    art::Handle< std::vector<raw::RawDigit> > digitVecHandle;
    evt.getByLabel("daq", digitVecHandle);

    if ( !digitVecHandle.isValid() ) {
      std::cout << "Missing daq info. skipping." << std::endl;
      return;
    }

    // Use the handle to get a particular (0th) element of collection.
    art::Ptr<raw::RawDigit> digitVec0(digitVecHandle, 0);
    
    std::cout << "Processing RawDigits..." << std::endl;

    idx_raw_start = idx_raw;

    // data size
    unsigned int dataSize = digitVec0->Samples(); //size of raw data vectors
    rawdigits.reserve( dataSize );
    
    // loop over all wires
    for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){
      art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);
      fWireID = digitVec->Channel();
      util::UBChannelReverseMap_t::const_iterator it_mapentry = fChannelReverseMap.find( (util::UBLArSoftCh_t)fWireID );
      if ( it_mapentry!=fChannelReverseMap.end() ) {
	fCrate = it_mapentry->second.crate;
	fSlot  = it_mapentry->second.card;
	fFemCH = it_mapentry->second.channel;
      }

      // uncompress the data, copy
      //raw::Uncompress(digitVec->ADCs(), rawdigits, digitVec->Compression());
      rawdigits = digitVec->ADCs();

      fTRawDigits->Fill();
      idx_raw++;

    }// end of wire loop

    idx_raw_nentries = idx_raw-idx_raw_start;
    fTRawDigitIndex->Fill();

  }
  
  // =============================================================== 
  // OpticalRawDigits
  void RawDigitWriter::GetOpticalRawDigits(const art::Event& evt) {
    // we have two type of arrays to packup
    // (1) frames of waveforms
    // (2) table with information on the frames
    // They all go into the same file

    art::ServiceHandle<geo::UBOpReadoutMap> ub_pmt_channel_map;
    art::Handle< std::vector< raw::OpDetWaveform > > wfHandle;
    auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    fTrigTimeStamp = ts->TriggerTime();
    fBeamTimeStamp = ts->BeamGateTime();

    std::cout << "OpticalDRAM: Trigger time=" << ts->TriggerTime() << " Beam gate time=" << ts->BeamGateTime() << std::endl;
    //ts->debugReport();
    idx_op_start = idx_op;
    
    for ( unsigned int cat=0; cat<(unsigned int)opdet::NumUBOpticalChannelCategories; cat++ ) {
      //std::stringstream ss;
      //ss << "pmtreadout" << opdet::UBOpChannelEnumName( (opdet::UBOpticalChannelCategory_t)cat );

      evt.getByLabel( "pmtreadout", opdet::UBOpChannelEnumName( (opdet::UBOpticalChannelCategory_t)cat ), wfHandle);

      if ( !wfHandle.isValid() ) {
	std::cout << "Missing pmtreadout/" << opdet::UBOpChannelEnumName( (opdet::UBOpticalChannelCategory_t)cat ) << " info. skipping." << std::endl;
	return;
      }

      std::vector<raw::OpDetWaveform> const& opwfms(*wfHandle);

      if ( opwfms.size()==0 ) {
	std::cout << opdet::UBOpChannelEnumName( (opdet::UBOpticalChannelCategory_t)cat ) << " zero waveforms." << std::endl;
	continue;
      }

      for(auto &wfm : opwfms )  {
	fOpReadoutCH = wfm.ChannelNumber();
	fCategory = (int)cat;
	fType = (int)ub_pmt_channel_map->GetChannelType( (unsigned int)fOpReadoutCH );
	unsigned int c, s, f;
	ub_pmt_channel_map->GetCrateSlotFEMChFromReadoutChannel( fOpReadoutCH, c, s, f );
	fOpCrate = (int)c;
	fOpSlot = (int)s;
	fOpFemCH = (int)f;
	fTimeStamp = wfm.TimeStamp();
	fFrame = ts->OpticalClock().Frame( fTimeStamp );
	fSample = ts->OpticalClock().Sample( fTimeStamp );

	opdetwaveforms.clear();
	opdetwaveforms.reserve( wfm.size() );
	for ( auto &adc : wfm )
	  opdetwaveforms.push_back( (short)adc );
	fTOpDetWaveforms->Fill();
	idx_op++;
      }
    }// end of category loop
    idx_op_nentries = idx_op-idx_op_start;
    fTOpDigitIndex->Fill();
  }

  DEFINE_ART_MODULE(RawDigitWriter)

} // namespace zmqds
 
