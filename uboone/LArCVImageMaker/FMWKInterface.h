#ifndef __FMWKINTERFACE_H__
#define __FMWKINTERFACE_H__

#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larcore/Geometry/Geometry.h"
#include "fhiclcpp/ParameterSet.h"
#include "Base/PSet.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/Simulation/SimChannel.h"

namespace supera {
  typedef larcv::PSet        Config_t;
  typedef recob::Wire        LArWire_t;
  typedef raw::OpDetWaveform LArOpDigit_t;
  typedef recob::Hit         LArHit_t;
  typedef simb::MCTruth      LArMCTruth_t;
  typedef sim::MCTrack       LArMCTrack_t;
  typedef sim::MCShower      LArMCShower_t;
  typedef sim::SimChannel    LArSimCh_t;
  typedef sim::MCStep        LArMCStep_t;
}
//
// Utility functions (geometry, lar properties, etc.)
//
namespace supera {
  
  //typedef ::fhicl::ParameterSet Config_t;

  /// Channel number to wire ID
  ::geo::WireID ChannelToWireID(unsigned int ch);
  
  /// DriftVelocity in cm/us
  double DriftVelocity();
  
  /// Number of time ticks
  unsigned int NumberTimeSamples();
  
  /// Number of channels
  unsigned int Nchannels();
  
  /// Number of planes
  unsigned int Nplanes();
  
  /// Number of wires
  unsigned int Nwires(unsigned int plane);
  
  /// Nearest wire
  unsigned int NearestWire(const TVector3& xyz, unsigned int plane);
  
  /// G4 time to TPC tick
  int TPCG4Time2Tick(double ns);
  
  /// TPC TDC to Tick
  double TPCTDC2Tick(double tdc);

  /// per-plane tick offset
  double PlaneTickOffset(size_t plane0, size_t plane1);

  void ApplySCE(geo::Point_t& pt);

  void ApplySCE(double x, double y, double z);

  void ApplySCE(double* xyz);

  void ApplySCE(TVector3& xyz);
  
}

#endif
