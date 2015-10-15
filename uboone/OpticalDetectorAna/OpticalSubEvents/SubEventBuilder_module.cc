////////////////////////////////////////////////////////////////////////
// Class:       SubEventBuilder
// Module Type: SubEventBuilder
// File:        SubEventBuilder_module.cc
//
// Generated at Thu Oct  8 13:18:22 2015 by Taritree Wongjirad using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <set>
#include <map>
#include <vector>
#include <iostream>

// LArSoft
#include "RawData/OpDetWaveform.h"
#include "RecoBase/OpFlash.h"
#include "RecoBase/OpHit.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/TimeService.h"
#include "Geometry/Geometry.h"
#include "Geometry/OpDetGeo.h"

#include "uboone/OpticalDetectorAna/OpticalSubEvents/subevent_algo/FlashList.hh"
#include "uboone/OpticalDetectorAna/OpticalSubEvents/subevent_algo/SubEvent.hh"
#include "uboone/OpticalDetectorAna/OpticalSubEvents/subevent_algo/SubEventList.hh"
#include "uboone/OpticalDetectorAna/OpticalSubEvents/subevent_algo/SubEventModule.hh"
#include "uboone/OpticalDetectorAna/OpticalSubEvents/subevent_algo/SubEventModConfig.hh"
#include "uboone/OpticalDetectorAna/OpticalSubEvents/subevent_algo/WaveformData.hh"
#include "uboone/OpticalDetectorAna/OpticalSubEvents/subevent_algo/pedestal.hh"

class SubEventBuilder;

class SubEventBuilder : public art::EDProducer {
public:
  explicit SubEventBuilder(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  SubEventBuilder(SubEventBuilder const &) = delete;
  SubEventBuilder(SubEventBuilder &&) = delete;
  SubEventBuilder & operator = (SubEventBuilder const &) = delete;
  SubEventBuilder & operator = (SubEventBuilder &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;


private:

  // Declare member data here.
  std::string fOpDetInputModule;
  subevent::SubEventModConfig fConfig;
  bool fMakeOpFlash;
  double fTrigCoinc;

  bool prepWaveforms( art::Event& event, subevent::WaveformData& wfms );
  bool gatherWaveforms( art::Event& event, subevent::WaveformData& wfms );
  void baselineCorrectWaveforms( subevent::WaveformData& wfms );
  void makeOpFlashes( art::Event& e, subevent::SubEventList& subevents, subevent::WaveformData& wfms );
  void GetHitGeometryInfo(subevent::Flash const& flash,
			  geo::Geometry const& geom,
			  std::vector<double> & sumw,
			  std::vector<double> & sumw2,
			  double & sumy, double & sumy2,
			  double & sumz, double & sumz2);
  
  
};


SubEventBuilder::SubEventBuilder(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  fOpDetInputModule           = p.get<std::string>( "inputModule");
  fMakeOpFlash                = p.get<bool>( "makeOpFlash", true );
  fTrigCoinc                  = p.get<double>       ("TrigCoinc");

  fConfig.cfdconfig.threshold = p.get<int>( "threshold" );
  fConfig.cfdconfig.deadtime  = p.get<int>( "deadtime" );
  fConfig.cfdconfig.delay     = p.get<int>( "delay" );
  fConfig.cfdconfig.width     = p.get<int>( "width" );
  fConfig.cfdconfig.gate      = p.get<int>( "gate" );
  fConfig.cfdconfig_pass2.threshold = p.get<int>( "pass2_threshold" );
  fConfig.cfdconfig_pass2.deadtime  = p.get<int>( "pass2_deadtime" );
  fConfig.cfdconfig_pass2.delay     = p.get<int>( "pass2_delay" );
  fConfig.cfdconfig_pass2.width     = p.get<int>( "pass2_width" );
  fConfig.cfdconfig_pass2.gate      = p.get<int>( "pass2_gate" );
  fConfig.fastfraction        = p.get<double>( "fastfraction", 0.8 );
  fConfig.slowfraction        = p.get<double>( "slowfraction", 0.3 );
  fConfig.fastconst_ns        = p.get<double>( "fastconst_ns", 6.0 );
  fConfig.slowconst_ns        = p.get<double>( "slowconst_ns", 1600.0 );
  fConfig.noslowthreshold     = p.get<double>( "noslowthreshold", 40.0 );
  fConfig.pedsamples          = p.get<int>( "pedsamples", 100 );
  fConfig.npresamples         = p.get<int>( "npresamples", 5 );
  fConfig.maxchflashes        = p.get<int>( "maxchflashes", 30 );
  fConfig.pedmaxvar           = p.get<double>( "pedmaxvar", 1.0 );
  fConfig.spe_sigma           = p.get< double >( "spe_sigma", 15.625*4.0 );
  fConfig.hgslot              = p.get<int>( "hgslot", 5 );
  fConfig.lgslot              = p.get<int>( "lgslot", 6 );
  fConfig.flashgate           = p.get<int>( "flashgate", 10 );
  fConfig.maxsubeventloops    = p.get<int>( "maxsubeventloops", 50 );
  fConfig.ampthresh           = p.get<double>( "ampthresh" );
  fConfig.hitthresh           = p.get<int>( "hitthresh" );
  fConfig.nspersample         = p.get<double>( "nspersample", 15.625 );

  // Call appropriate produces<>() functions here.
  produces< std::vector< subevent::SubEvent > >();
  produces< subevent::FlashList >( "unclaimedFlashes" );

  if ( fMakeOpFlash ) {
    produces<std::vector< recob::OpFlash> >();
    produces<std::vector< recob::OpHit> >();
    produces<art::Assns<recob::OpFlash, recob::OpHit> >();
  }
}

void SubEventBuilder::produce(art::Event & e)
{
  // Implementation of required member function here.  

  // get waveforms and prep them
  subevent::WaveformData wfms;
  bool ok = prepWaveforms( e, wfms );
  if (!ok) {
    std::cout << "trouble loading waveforms!" << std::endl;
    return;
  }

  // declare containers for output
  subevent::SubEventList subevents;
  std::unique_ptr< subevent::FlashList > unclaimed_flashes( new subevent::FlashList() );

  // hack: make spe calibration table
  std::map< int, double > pmtspemap;
  for (int i=0; i<36; i++)
    pmtspemap[i] = 20.0;
  
  // Find subevents
  formSubEvents( wfms, fConfig, pmtspemap, subevents, *unclaimed_flashes );
  std::cout << "subevents formed. now store in event." << std::endl;
  
  // make opflash
  if ( fMakeOpFlash )
    makeOpFlashes( e, subevents, wfms );

  // transfer subevents to the event
  std::unique_ptr< std::vector< subevent::SubEvent > > subeventvec( new std::vector< subevent::SubEvent > );
  for ( subevent::SubEventListIter it=subevents.begin(); it!=subevents.end(); it++ ) {
    subeventvec->emplace_back( *it );
  }
  e.put( std::move( subeventvec ) );
  e.put( std::move( unclaimed_flashes ), "unclaimedFlashes" );
  
}

bool SubEventBuilder::prepWaveforms( art::Event& event, subevent::WaveformData& wfms ) {
  bool ok = gatherWaveforms( event, wfms );
  //baselineCorrectWaveforms( wfms );
  return ok;
}

bool SubEventBuilder::gatherWaveforms( art::Event& event, subevent::WaveformData& wfms ) {

  art::ServiceHandle<util::TimeService> ts;

  art::Handle< std::vector< raw::OpDetWaveform > > hgwfmHandle;
  bool loadedhg = event.getByLabel( fOpDetInputModule, "OpdetBeamHighGain", hgwfmHandle );
  
  art::Handle< std::vector< raw::OpDetWaveform > > lgwfmHandle;
  bool loadedlg = event.getByLabel( fOpDetInputModule, "OpdetBeamLowGain", lgwfmHandle );

  if ( !loadedhg || !loadedlg ) {
    std::cout << "Could not load optdet waveforms!" << std::endl;
    return false;
  }

  std::vector< double > wfmstore;
  std::set< int > use_lowgain_wfm;

  // High-Gain Waveforms
  unsigned int biggestcosmiclen = 0;
  unsigned int beamwinlen = 0;
  double shortestdt = 1e10;
  for ( auto const& opdetData: (*hgwfmHandle) ) {
    if ( opdetData.size()<600 ) { // arbitrary!
      if ( opdetData.size()>biggestcosmiclen )
	biggestcosmiclen = opdetData.size();
      //std::cout << "skipping cosmic window: " << opdetData.size() << std::endl;
      double dtcosmic = opdetData.TimeStamp() - ts->BeamGateTime();
      if ( fabs( dtcosmic )<shortestdt )
	shortestdt = fabs( dtcosmic );
      continue; // skip cosmic windows
    }
    raw::Channel_t channel = opdetData.ChannelNumber();
    wfmstore.clear();

    if ( wfmstore.size() != opdetData.size() ) {
      wfmstore.reserve( opdetData.size() );
    }
    if ( beamwinlen<opdetData.size() )
      beamwinlen = opdetData.size();

    bool marked_for_lg = false;
    double adcmax = 0.0;
    //std::cout << " ch " << channel << ": ";
    for ( auto adc: opdetData ) {
      if ( !marked_for_lg && adc>4090 ) {
	marked_for_lg = true;
	use_lowgain_wfm.insert( (int)channel );
      }
      wfmstore.push_back( (double)adc );
      //std::cout << adc << " ";
      if ( adcmax < (double)adc )
	adcmax = (double)adc;
    }
    //std::cout << std::endl;
    //double ped =  subevent::removePedestal( wfmstore, 20, 2.0, 2047.0 );
    subevent::removePedestal( wfmstore, 20, 2.0, 2047.0 );
    //std::cout << "  adc max: " << adcmax-ped << " pedestal=" << ped << std::endl;
    wfms.set( (int)channel, wfmstore, false );
    // store time stamp for later
    wfms.storeTimeInfo( (int)channel,  ts->OpticalClock().Frame( opdetData.TimeStamp() ), opdetData.TimeStamp() );
  }

  // Low-Gain Waveforms
  if ( use_lowgain_wfm.size() > 0 ) {
    for ( auto const& opdetData: (*lgwfmHandle) ) {
      if ( opdetData.size()<600 )
	continue; // skip cosmic windows
      raw::Channel_t channel = opdetData.ChannelNumber();
      if ( use_lowgain_wfm.find( (int)channel )!=use_lowgain_wfm.end() ) {
	// marked to replace hg
	wfmstore.clear();
	double adcmax = 0.0;
	std::cout << " LG ch " << channel << ": ";
	for ( auto adc: opdetData ) {
	  wfmstore.push_back( (double)adc );
	  if ( adcmax < (double)adc )
	    adcmax = (double)adc;
	}
	double ped = subevent::removePedestal( wfmstore, 20, 2.0, 2047.0 );
	std::cout << "  adc max: " << 10.0*(adcmax-ped) << " pedestal=" << ped << std::endl;
	
	// scale up
	for (int i=0; i<(int)wfmstore.size(); i++ )
	  wfmstore.at(i) *= 10.0;
	
	wfms.set( (int)channel, wfmstore, true ); // low gain waveform
	

      }// if marked for lg replacement
    }// loop over lowgain waveforms
  } // if replacements necessary

  std::cout << "gathered beam waveforms: beamwin length=" << beamwinlen << " cosmicwin length=" << biggestcosmiclen << " shorteddt=" << shortestdt*1000 << " ns" << std::endl;
  
  return true;
}

void SubEventBuilder::baselineCorrectWaveforms( subevent::WaveformData& wfms ) {
  double RC = 50000.0; // ns
  double f = exp( -15.625/RC );
  std::vector< double > qcap;
  
  for ( subevent::ChannelSetIter itch=wfms.chbegin(); itch!=wfms.chend(); itch++ ) {

    // calc charge on capacitor
    qcap.clear();
    if ( qcap.size()!=wfms.get( *itch ).size() )
      qcap.reserve( wfms.get( *itch ).size() );
    qcap.push_back( 0.0 );
    for ( unsigned int tdc=1; tdc<wfms.get( *itch ).size(); tdc++ ) {
      double adc = wfms.get( *itch ).at( tdc );
      double q = 50.0*adc/RC + f*qcap.at( tdc-1 );
      qcap.push_back ( q );
    }
    
    // subtract correction
    for ( unsigned int tdc=0; tdc<wfms.get( *itch ).size(); tdc++ ) {
      wfms.get( *itch ).at( tdc ) -= qcap.at( tdc );
    }
  }
}

void SubEventBuilder::GetHitGeometryInfo(subevent::Flash const& flash,
					 geo::Geometry const& geom,
					 std::vector<double> & sumw,
					 std::vector<double> & sumw2,
					 double & sumy, double & sumy2,
					 double & sumz, double & sumz2)
{
  double xyz[3];
  geom.OpDetGeoFromOpChannel( (unsigned int)flash.ch ).GetCenter(xyz);

  double PEThisHit = flash.area;
  for(size_t p=0; p!=geom.Nplanes(); p++){
    unsigned int w = geom.NearestWire(xyz,p);
    sumw.at(p)  += w*PEThisHit;
    sumw2.at(p) += w*w*PEThisHit;
  }
  
  sumy+=xyz[1]*PEThisHit; sumy2+=xyz[1]*xyz[1]*PEThisHit;
  sumz+=xyz[2]*PEThisHit; sumz2+=xyz[2]*xyz[2]*PEThisHit;
}

void SubEventBuilder::makeOpFlashes( art::Event& e, subevent::SubEventList& subevents, subevent::WaveformData& wfms ) {
  
  art::ServiceHandle<util::TimeService> ts;
  art::ServiceHandle<geo::Geometry> geom;
  geo::Geometry const& Geometry(*geom);

  std::unique_ptr< std::vector< recob::OpFlash > >  opflashes( new std::vector< recob::OpFlash > );
  std::unique_ptr< std::vector< recob::OpHit > > ophits( new std::vector< recob::OpHit > );
  std::unique_ptr< art::Assns<recob::OpFlash, recob::OpHit> >  AssnPtr( new art::Assns<recob::OpFlash, recob::OpHit> );

  for ( subevent::SubEventListIter isubevent=subevents.begin(); isubevent!=subevents.end(); isubevent++ ) {
    subevent::SubEvent& asubevent = (*isubevent);
    // gather info to make an opflash for this subevent:
    //
    // OpFlash(double time, double timewidth, double abstime, unsigned int frame,
    // 	    std::vector< double > PEperOpDet,
    // 	    bool InBeamFrame=0, int OnBeamTime=0, double FastToTotal=1,
    // 	    double yCenter=0, double yWidth=0,
    // 	    double zCenter=0, double zWidth=0,
    // 	    std::vector<double> WireCenters = std::vector<double>(0),
    // 	    std::vector<double> WireWidths  = std::vector<double>(0));
    double timestamp = wfms.getTimestamp(asubevent.flashes.get(0).ch); // us
    double abstime = timestamp + asubevent.tmax_sample*ts->OpticalClock().TickPeriod(); // us
    double reltime = abstime - ts->BeamGateTime(); // us
    double width = ( asubevent.tend_sample - asubevent.tstart_sample )*ts->OpticalClock().TickPeriod(); // us
    std::cout << "opflash: timestamp=" << timestamp << " us"
	      << " abstime=" << abstime << " us"
	      << " reltime=" << reltime << " us"  
	      << " width=" << width << " us" 
	      << " beamgatetime=" << ts->BeamGateTime() << " us"
	      << " tickperiod=" << ts->OpticalClock().TickPeriod() << " us" << std::endl;
    // Emprical corrections to get the Frame right
    // // Eventual solution - remove frames
    // taken from OpFlashAlg.cxx
    unsigned int frame = ts->OpticalClock().Frame( abstime );
    unsigned int trigframe = ts->OpticalClock().Frame( ts->BeamGateTime() );
    bool InBeamFrame = (frame==trigframe);
    int OnBeamTime =0;
    if( std::abs(reltime) < fTrigCoinc ) OnBeamTime=1; // this can't be right, can it?
    double FastToTotal = 0.;

    // geo stuff I copied
    //std::vector<double> PEs(geom->MaxOpChannel()+1,0.0);
    unsigned int Nplanes = geom->Nplanes();
    std::vector<double> sumw(Nplanes,0), sumw2(Nplanes,0);
    double sumy=0, sumz=0, sumy2=0, sumz2=0;
    
    // sum pe for each opdet and do geo stuff
    std::vector< double > PEperOpDet( geom->NOpDets(), 0.0 );
    for ( subevent::FlashListIter iflash=asubevent.flashes.begin(); iflash!=asubevent.flashes.end(); iflash++ ) {
      int iopdet = geom->OpDetFromOpChannel( (unsigned int)(*iflash).ch );
      PEperOpDet.at( iopdet ) += (*iflash).area; // need calibration constants here
      GetHitGeometryInfo( (*iflash), Geometry, sumw, sumw2, sumy, sumy2, sumz, sumz2 );
    }
    double meany = sumy/asubevent.totpe;
    double meanz = sumz/asubevent.totpe;
    double widthy  = sumy2*asubevent.totpe - sumy*sumy;
    if ( widthy>0.0 )
      widthy = sqrt( widthy );
    else
      widthy = 0.0;
    double widthz = sumz2*asubevent.totpe - sumz*sumz;
    if ( widthz>0.0 )
      widthz = sqrt( widthz );
    widthz = 0.0;

    std::vector<double> WireCenters(Nplanes,0);
    std::vector<double> WireWidths(Nplanes,0);
    for(size_t p=0; p!=Nplanes; ++p){
      WireCenters.at(p) = sumw.at(p)/asubevent.totpe;
      WireWidths.at(p)  = sumw2.at(p)*asubevent.totpe - sumw.at(p)*sumw.at(p);
      if ( WireWidths.at(p)> 0 ) WireWidths.at(p) = sqrt( WireWidths.at(p) );
      else WireWidths.at(p) = 0.;
    }

    // finally, make the opflash
    recob::OpFlash aopflash( reltime, width, abstime, frame, PEperOpDet, 
			     InBeamFrame, OnBeamTime, FastToTotal,
			     meany, widthy, meanz, widthz,
			     WireCenters, WireWidths );

    // make ophits and associate it with the flashes
    for ( subevent::FlashListIter iflash=asubevent.flashes.begin(); iflash!=asubevent.flashes.end(); iflash++ ) {
      double flash_abstime = (*iflash).tmax*ts->OpticalClock().TickPeriod() + wfms.getTimestamp( (*iflash).ch );
      double flash_reltime = flash_abstime -  ts->BeamGateTime();
      unsigned int flash_frame = ts->OpticalClock().Frame( flash_abstime );
      double flash_width = ( (*iflash).tend-(*iflash).tstart )*ts->OpticalClock().TickPeriod();
      ophits->emplace_back( (int)(*iflash).ch, flash_reltime, flash_abstime, flash_frame, flash_width, (*iflash).area, (*iflash).maxamp, (*iflash).area/100.0, 0.0 ); // wants microseconds
    }
    opflashes->emplace_back( std::move( aopflash ) );
  }//end of subevent loop
  
  
  e.put( std::move( opflashes ) );
  e.put( std::move( ophits ) );
  e.put( std::move( AssnPtr ) );
  
}

DEFINE_ART_MODULE(SubEventBuilder)
