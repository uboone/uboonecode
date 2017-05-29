#ifndef __FMWKINTERFACE_CXX__
#define __FMWKINTERFACE_CXX__

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "Base/larbys.h"
#include "FMWKInterface.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesServiceStandard.h"
#include "lardata/DetectorInfoServices/LArPropertiesServiceStandard.h"
#include "lardata/DetectorInfoServices/DetectorClocksServiceStandard.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

namespace supera {

  ::geo::WireID ChannelToWireID(unsigned int ch)
  { 
      auto const* geom = ::lar::providerFrom<geo::Geometry>();
      return geom->ChannelToWire(ch).front();
  }
  
  double DriftVelocity()
  { 
    auto const* detp = ::lar::providerFrom<detinfo::DetectorPropertiesService>();
    return detp->DriftVelocity(); 
  }
  
  unsigned int NumberTimeSamples()
  {
    throw ::larcv::larbys("NumberTimeSamples function not available!");
    auto const* detp = ::lar::providerFrom<detinfo::DetectorPropertiesService>();
    return detp->NumberTimeSamples(); 
  }
  
  unsigned int Nchannels()
  {
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    return geom->Nchannels();
  }
  
  unsigned int Nplanes()
  { 
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    return geom->Nplanes();
  }
  
  unsigned int Nwires(unsigned int plane)
  { 
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    return geom->Nwires(plane); 
  }
  
  unsigned int NearestWire(const TVector3& xyz, unsigned int plane)
  {
    double min_wire=0;
    double max_wire=Nwires(plane)-1;
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    
    double wire = geom->WireCoordinate(xyz[1],xyz[2],plane,0,0) + 0.5;
    if(wire<min_wire) wire = min_wire;
    if(wire>max_wire) wire = max_wire;
    
    return (unsigned int)wire;
  }

  unsigned int NearestWire(const double* xyz, unsigned int plane)
  {
    double min_wire=0;
    double max_wire=Nwires(plane)-1;
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    
    double wire = geom->WireCoordinate(xyz[1],xyz[2],plane,0,0) + 0.5;
    if(wire<min_wire) wire = min_wire;
    if(wire>max_wire) wire = max_wire;
    
    return (unsigned int)wire;
  }

  double WireAngleToVertical(unsigned int plane)
  {
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    return geom->WireAngleToVertical(geo::View_t(plane));
  }

  double WirePitch(size_t plane)
  {
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    return geom->WirePitch(0,1,plane);
  }

  double DetHalfWidth() 
  {
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    return geom->DetHalfWidth();
  }

  double DetHalfHeight() 
  {
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    return geom->DetHalfHeight();
  }

  double DetLength() 
  {
    auto const* geom = ::lar::providerFrom<geo::Geometry>();
    return geom->DetLength();
  }
  
  int TPCG4Time2Tick(double ns)
  { 
    auto const* ts = ::lar::providerFrom<detinfo::DetectorClocksService>();      
    return ts->TPCG4Time2Tick(ns); 
  }

  int TPCG4Time2TDC(double ns)
  {
    auto const* ts = ::lar::providerFrom<detinfo::DetectorClocksService>();
    return ts->TPCG4Time2TDC(ns);
  }
  
  double TPCTDC2Tick(double tdc)
  { 
    auto const* ts = ::lar::providerFrom<detinfo::DetectorClocksService>();
    return ts->TPCTDC2Tick(tdc); 
  }

  double TPCTickPeriod()
  {
    auto const* ts = ::lar::providerFrom<detinfo::DetectorClocksService>();
    return ts->TPCClock().TickPeriod();
  }

  double TriggerOffsetTPC()
  {
    auto const* ts = ::lar::providerFrom<detinfo::DetectorClocksService>();
    return ts->TriggerOffsetTPC();
  }
  
  double PlaneTickOffset(size_t plane0, size_t plane1)
  {
    static double pitch = ::lar::providerFrom<geo::Geometry>()->PlanePitch();
    static double tick_period = ::lar::providerFrom<detinfo::DetectorClocksService>()->TPCClock().TickPeriod();
    return (plane1 - plane0) * pitch / DriftVelocity() / tick_period;
  }

  void ApplySCE(geo::Point_t& pt)
  {
    auto pos = ::lar::providerFrom<spacecharge::SpaceChargeService>()->GetPosOffsets(pt);
    pt.SetX(pt.X() - pos.X() + 0.7);
    pt.SetY(pt.Y() + pos.Y());
    pt.SetZ(pt.Z() + pos.Z());
  }

  void ApplySCE(double& x, double& y, double& z)
  {
    static geo::Point_t pt;
    pt.SetX(x);
    pt.SetY(y);
    pt.SetZ(z);
    ApplySCE(pt);
    x = pt.X();
    y = pt.Y();
    z = pt.Z();
  }

  void ApplySCE(double* xyz)
  {
    static geo::Point_t pt;
    pt.SetX(xyz[0]);
    pt.SetY(xyz[1]);
    pt.SetZ(xyz[2]);
    ApplySCE(pt);
    xyz[0] = pt.X();
    xyz[1] = pt.Y();
    xyz[2] = pt.Z();
  }

  void ApplySCE(TVector3& xyz)
  {
    static geo::Point_t pt;
    pt.SetX(xyz[0]);
    pt.SetY(xyz[1]);
    pt.SetZ(xyz[2]);
    ApplySCE(pt);
    xyz[0] = pt.X();
    xyz[1] = pt.Y();
    xyz[2] = pt.Z();
  }

}

#endif
