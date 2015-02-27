/**
 * \file UBOptDetChannelMap.h
 *
 * \ingroup OpticalDetectorSim
 * 
 * \brief Service providing map between OpDet and PMT Readout Channels
 *
 * @author taritree
 */

#ifndef UBOPTDETCHANNELMAP_H
#define UBOPTDETCHANNELMAP_H

#include <vector>
#include <string>
#include <map>
#include "Geometry/Geometry.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "Geant4/G4ThreeVector.hh"

namespace opdet {
  /**
     \class UBOptDetChannelMap
     FICHL Modfiable Map between LArG4 OpDet ID, UB PMT ID,  Readout Chanel ID.
     * LArG4 Opdet ID is index of Geant4 sens. det. volume. Numbered in the order it is
        loaded into memory -- treated as random. Needed to get access to hits.
     * UB PMT ID: ID number assigned to physical PMT. We match this PMT and the sens. volume
        repesentation using geometry.  User lists location of physical PMT and we match ths position
	to the center of the closest PMT in the simulation.
	Don't know the PMT coordinates? Set DumpG4OptDetLocation to True! Then label.
     * Readout Channel: Each PMT is split into 2 high and 2 low gain signals.
        For one copy, the signals from the 32 PMTs and 4 paddles go into 3 Shaper cards.
        There are also several trigger channels!  Assignment of channels can be set in the FICHL file.
  */
  class UBOptDetChannelMap {
    
  public:
    
    /// Default constructor
    UBOptDetChannelMap(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg);
    
    /// Default destructor
    virtual ~UBOptDetChannelMap(){};

    void reconfigure(fhicl::ParameterSet const& pset);

    // From G4 Opt Det
    int PMTIDfromG4OptDet( int g4optdet );
    std::vector<int> ReadoutChannelsFromG4OptDet(int g4optdet );
    std::vector<int> HighReadoutChannelsFromG4OptDet(int g4optdet );
    std::vector<int> LowReadoutChannelsFromG4OptDet(int g4optdet );
    G4ThreeVector GetPositionFromG4OptDet( int g4optdet );

    // From PMT ID
    int G4OptDetIDfromPMTID( int pmtid );
    std::vector<int> ReadoutChannelsFromPMTID( int pmtid );
    std::vector<int> HighReadoutChannelsFromPMTID( int pmtid );
    std::vector<int> LowReadoutChannelsFromPMTID( int pmtid );
    G4ThreeVector GetPositionFromPMTID( int pmtid );

    // From Readout Channel
    int PMTIDfromReadoutChannel( int roch );
    int G4OptDetIDfromReadoutChannel( int roch );
    bool isLowGainChannel( int roch );
    bool isHighGainChannel( int roch );
    G4ThreeVector GetPositionFromReadoutChannel( int channel );
    
  private:
    
    int kNDetectors;
    bool kDumpG4OptDetPos;
    bool kWriteMapToTFile;
    std::vector< G4ThreeVector > kDetPositions;
    std::map< int, int > g4opdet2pmt; //< [ g4, pmtid ]
    std::map< int, int > pmt2g4opdet; //< [ pmtid, g4 ]
    std::map< std::string, int > pmtname2index; //< [ name, pmtid ]
    std::map< int, std::string > pmtindex2name; //< [ pmtid, name ]

    std::map< int, std::vector< int > > pmt2channel; //< [ pmtid, channel vector ]
    std::map< int, int > channel2pmt;                //< [ channel id, pmtid ]
    std::vector< std::vector<int> > lowgain_channel_ranges;  //< [ [lower, upper] ], ...
    std::vector< std::vector<int> > highgain_channel_ranges; //< [ [lower, upper] ], ...

  };

}


DECLARE_ART_SERVICE(opdet::UBOptDetChannelMap, LEGACY)

#endif
/** @} */ // end of doxygen group 

