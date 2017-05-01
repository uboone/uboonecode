/**
 * \file LightCharge.h
 *
 * \ingroup Algorithms
 * 
 * \brief Creates flash hypothesis from charge
 *
 * @author Marco Del Tutto
 */

/** \addtogroup Algorithms

    @{*/
#ifndef LIGHTCHARGE_H
#define LIGHTCHARGE_H

#include <iostream>
#include <numeric>
#include <functional>
#include <algorithm>
#include "uboone/LLSelectionTool/OpT0Finder/Base/BaseAlgorithm.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/CustomAlgoFactory.h"

namespace flashana{
/**
   \class LightCharge
   This class constructs a flash hypothesis from a collection of 3D hits. 
   3D hits need to have space location (x,y,z), deposited charge q
   and the plane from which the charge is taken.
 */

  class LightCharge : public flashana::BaseAlgorithm {
    
  public:
    
    /// Default constructor
    LightCharge(const std::string name="LightCharge");
    
    /// Default destructor
    ~LightCharge(){}

    // Setter function
    double SetChargeToLight(double x) {_charge_to_light = x; return _charge_to_light;}
      
    // Flash Hypothesis for Trajectory (Track)
    flashana::QCluster_t FlashHypothesis(const std::vector<flashana::Hit3D_t>) const;

    // Getter for light yield configured paramater
    double GetChargeToLight() const { return _charge_to_light; }


  protected:

    void _Configure_(const Config_t &pset);
 
    int    _plane_no;   
    double _charge_to_light;

  };
  
  /**
     \class flashana::LightChargeFactory
  */
  class LightChargeFactory : public CustomAlgoFactoryBase {
  public:
    /// ctor
    LightChargeFactory() { CustomAlgoFactory::get().add_factory("LightCharge",this); }
    /// dtor
    ~LightChargeFactory() {}
    /// creation method
    BaseAlgorithm* create(const std::string instance_name) { return new LightCharge(instance_name); }
  };
} 

#endif
/** @} */ // end of doxygen group 

