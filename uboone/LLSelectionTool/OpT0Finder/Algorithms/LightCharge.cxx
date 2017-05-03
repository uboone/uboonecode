#ifndef LIGHTCHARGE_CXX
#define LIGHTCHARGE_CXX

#include "LightCharge.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"

namespace flashana {

  static LightChargeFactory __global_LightChargeFactory__;

  LightCharge::LightCharge(const std::string name)
    : BaseAlgorithm(kCustomAlgo, name)
    , _plane_no        ( 2      )
    , _charge_to_light ( 40000. )
  {}

  void LightCharge::_Configure_(const Config_t &pset)
  {
    _plane_no         = pset.get< int >    ( "PlaneNumber"   );
    _charge_to_light  = pset.get< double > ( "ChargeToLight" );
  }


  QCluster_t LightCharge::FlashHypothesis(const std::vector<flashana::Hit3D_t> hit3d_v) const {

    QCluster_t result;
    result.clear();

    QPoint_t q_pt;

    for (auto hit3d : hit3d_v) {
      if (hit3d.plane != _plane_no) continue;
      q_pt.x = hit3d.x;
      q_pt.y = hit3d.y;
      q_pt.z = hit3d.z;
      q_pt.q = _charge_to_light * hit3d.q;
      result.emplace_back(q_pt); 
    }


    return result;
  }

  // More general Flash Hypothesis for tracks and showers, needs to be given a conversion factor to photons
  QCluster_t LightCharge::FlashHypothesisCharge(const std::vector<flashana::Hit3D_t> hit3d_v, double charge_to_light) {
    SetChargeToLight(charge_to_light);
    return LightCharge::FlashHypothesis(hit3d_v);
  }

}


#endif
