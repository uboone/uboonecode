#ifndef UBOPTDETCHANNELMAP_CXX
#define UBOPTDETCHANNELMAP_CXX

#include "UBOptDetChannelMap.h"
#include "UBOpticalException.h"
#include "LArG4/OpDetLookup.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/OpDetGeo.h"

namespace opdet {

  //---------------------------------------------------------------------------------
  UBOptDetChannelMap::UBOptDetChannelMap(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg)
  //---------------------------------------------------------------------------------
  {
    kDumpG4OptDetPos = false;
    this->reconfigure(pset);
  }
  
  //-----------------------------------------------------------
  void UBOptDetChannelMap::reconfigure(fhicl::ParameterSet const& pset)
  //-----------------------------------------------------------
  {
    
    kNDetectors      = pset.get< int >("numberOfDetectors");
    kDumpG4OptDetPos = pset.get< bool >("dumpG4OptDetPos");

    // sanity check
    art::ServiceHandle<geo::Geometry> geom;
    if((unsigned int)kNDetectors != geom->NOpChannels())
      throw UBOpticalException(Form("Number of detectors mapped in FICHL file (%d) is not the number of OpDets loaded into th G4 geometry (%u)!", kNDetectors, geom->NOpChannels() ));

    // read in positions that will be used to map physical detector to geant 4 representation
    kDetPositions.clear();
    for (int iop=0; iop<kNDetectors; iop++) {
      char entryname[50];
      
      sprintf(entryname,"OpDet%d",iop+1);
      std::vector< double > inputpos = pset.get< std::vector<double> >( entryname );
      kDetPositions.push_back( G4ThreeVector( inputpos.at(0), inputpos.at(1), inputpos.at(2) ) );
      
      sprintf(entryname,"OpDet%d_name",iop+1);
      std::string inputname = pset.get< std::string >( entryname );
      pmtname2index[ inputname ] = iop;
      pmtindex2name[ iop ] = inputname;
      
    }

    // read in low gain, high gain ranges
    std::vector< int > lorange_input = pset.get< std::vector<int> >( "ReadoutLowChannelRange" );
    std::vector< int > hirange_input = pset.get< std::vector<int> >( "ReadoutHighChannelRange" );
    for (int i=0; i++; i+=2 ) {
      int lo[2] = { lorange_input.at(i), lorange_input.at(i+1) };
      int hi[2] = { hirange_input.at(i), hirange_input.at(i+1) };
      std::vector< int > vlo( lo, lo+2 );
      std::vector< int > vhi( hi, hi+2 );
      lowgain_channel_ranges.push_back( vlo );
      highgain_channel_ranges.push_back( vhi );
    }

    // now make channel map
    for (int iop=0; iop<kNDetectors; iop++) {
      char entryname[50];
      sprintf(entryname,"OpDet%d_channels",iop+1);
      std::vector< int > chinput =  pset.get< std::vector<int> >( entryname );
      pmt2channel[ iop ] = chinput;
      for (std::vector<int>::iterator it_ch=chinput.begin(); it_ch!=chinput.end(); it_ch++) {
	channel2pmt[ *it_ch ] = iop;
      }
    }

    // we get the geo service to get access to the G4 OpDets.
    // we want their positions
    int icryo = 0;
    std::cout << "Matching positions" << std::endl;
    larg4::OpDetLookup* odlookup = larg4::OpDetLookup::Instance();
    for (int iop=0; iop<(int)geom->NOpChannels(); iop++) {
      //if ( kDumpG4OptDetPos ) {
      double xyz[3];
      geom->Cryostat(icryo).OpDet(iop).GetCenter(xyz);
      char out[1000];
      sprintf( out, "G4OptDet ID%d: center=(%.2f, %.2f, %.2f)",iop,xyz[0], xyz[1], xyz[2]);
      //std::cout << Form("G4OptDet ID%d: center=(%.2f, %.2f, %.2f)",iop,xyz[0], xyz[1], xyz[2]) << std::endl;
      std::cout << out << std::endl;
      //}
      // use listed position 
      G4ThreeVector pos = kDetPositions.at(iop);
      double dist = 0.;
      int closest_cryo = 0;
      int g4opid = odlookup->FindClosestOpDet( pos, dist, closest_cryo);
      g4opdet2pmt[g4opid] = iop;
      pmt2g4opdet[iop] = g4opid;
    }

  }

  // ===============================================================
  // Public Functions

  // From G4 Opt Det
  int UBOptDetChannelMap::PMTIDfromG4OptDet( int g4optdet ) {
    return g4opdet2pmt[ g4optdet ];
  }

  std::vector<int> UBOptDetChannelMap::ReadoutChannelsFromG4OptDet(int g4optdet ) {
    return ReadoutChannelsFromPMTID( PMTIDfromG4OptDet(g4optdet) );
  }
  
  std::vector<int> UBOptDetChannelMap::HighReadoutChannelsFromG4OptDet(int g4optdet ) {
    return HighReadoutChannelsFromPMTID( PMTIDfromG4OptDet(g4optdet) );
  }

  std::vector<int> UBOptDetChannelMap::LowReadoutChannelsFromG4OptDet(int g4optdet ) {
    return LowReadoutChannelsFromPMTID( PMTIDfromG4OptDet(g4optdet) );
  }

  G4ThreeVector UBOptDetChannelMap::GetPositionFromG4OptDet( int g4optdet ) {
    art::ServiceHandle<geo::Geometry> geom;
    int icryo = 0;
    double pos[3];
    geom->Cryostat(icryo).OpDet(g4optdet).GetCenter(pos);
    return G4ThreeVector( pos[0]*cm, pos[1]*cm, pos[2]*cm );
  }

  // From PMT ID
  int UBOptDetChannelMap::G4OptDetIDfromPMTID( int pmtid ) {
    return pmt2g4opdet[ pmtid ];
  }

  std::vector<int> UBOptDetChannelMap::ReadoutChannelsFromPMTID( int pmtid ) {
    return pmt2channel[ pmtid ];
  }

  std::vector<int> UBOptDetChannelMap::HighReadoutChannelsFromPMTID(int pmtid ) {
    std::vector< int > highs;
    std::vector< int > channels = ReadoutChannelsFromPMTID( pmtid );
    for ( std::vector< int >::iterator it=channels.begin(); it!=channels.end(); it++ )
      if ( isHighGainChannel( *it ) )
	highs.push_back( *it );
    return highs;
  }

  std::vector<int> UBOptDetChannelMap::LowReadoutChannelsFromPMTID(int pmtid ) {
    std::vector< int > lows;
    std::vector< int > channels = ReadoutChannelsFromPMTID( pmtid );
    for ( std::vector< int >::iterator it=channels.begin(); it!=channels.end(); it++ )
      if ( isLowGainChannel( *it ) )
	lows.push_back( *it );
    return lows;
  }

  G4ThreeVector UBOptDetChannelMap::GetPositionFromPMTID( int pmtid ) {
    return GetPositionFromG4OptDet( G4OptDetIDfromPMTID( pmtid ) );
  }

  // From Readout Channel
  int UBOptDetChannelMap::PMTIDfromReadoutChannel( int roch ) {
    return channel2pmt[ roch ];
  }

  int UBOptDetChannelMap::G4OptDetIDfromReadoutChannel( int roch ) {
    return G4OptDetIDfromPMTID( PMTIDfromReadoutChannel( roch ) );
  }

  bool UBOptDetChannelMap::isLowGainChannel( int ich ) {
    for ( std::vector< std::vector<int> >::iterator it=lowgain_channel_ranges.begin(); it!=lowgain_channel_ranges.end(); it++ ) {
      if ( (*it).at(0)>=ich && ich<=(*it).at(1) )
	return true;
    }
    return false;
  }

  bool UBOptDetChannelMap::isHighGainChannel( int ich ) {
    for ( std::vector< std::vector<int> >::iterator it=highgain_channel_ranges.begin(); it!=highgain_channel_ranges.end(); it++ ) {
      if ( (*it).at(0)>=ich && ich<=(*it).at(1) )
	return true;
    }
    return false;
  }
  
  G4ThreeVector UBOptDetChannelMap::GetPositionFromReadoutChannel( int channel ) {
    return GetPositionFromG4OptDet( G4OptDetIDfromReadoutChannel( channel ) );
  }
  // ===============================================================

  DEFINE_ART_SERVICE(UBOptDetChannelMap)
    
}
#endif
