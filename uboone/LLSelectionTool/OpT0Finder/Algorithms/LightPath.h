/**
 * \file LightPath.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class LightPath
 *
 * @author Rui
 */

/** \addtogroup Algorithms

    @{*/
#ifndef LIGHTPATH_H
#define LIGHTPATH_H

#include <iostream>
#include <numeric>
//#include "GeoAlgo/GeoAlgo.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoAlgo.h"
#include <functional>
#include <algorithm>
//#include "OpT0Finder/Base/BaseAlgorithm.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/BaseAlgorithm.h"
namespace flashana{
/**
   \class LightPath
   User defined class LightPath ... these comments are used to generate
   doxygen documentation!
 */

  class LightPath : public flashana::BaseAlgorithm {
    
  public:
    
    /// Default constructor
    LightPath(const std::string name="LightPath");
    
    /// Default destructor
    ~LightPath(){}

    //void Configure(const ::fcllite::PSet &pset);
    void Configure(const ::fhicl::ParameterSet &pset);
    
    // Setter function
    double Set_Gap      ( double x) { _gap   =x;      return _gap;}
      
    // Flash Hypothesis for Trajectory (Track)
    flashana::QCluster_t FlashHypothesis(const ::geoalgo::Trajectory& trj) const;

    void QCluster(const ::geoalgo::Vector& pt_1,
                  const ::geoalgo::Vector& pt_2,
                  flashana::QCluster_t& Q_cluster) const;

    // Getter for light yield configured paramater
    double GetLightYield() const { return _light_yield; }


  protected:
    double _gap;
    double _light_yield;
    double _dEdxMIP;
    ::geoalgo::GeoAlgo _geoAlgo;
  };
} 

#endif
/** @} */ // end of doxygen group 

