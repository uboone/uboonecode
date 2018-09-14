/**
 * \file SuperaLArFlow.h
 *
 * \ingroup Package_Name
 * 
 * \brief Class def header for a class SuperaLArFlow
 *
 * @author taritree
 */

/** \addtogroup Package_Name

    @{*/
#ifndef __SUPERALARFLOW_H__
#define __SUPERALARFLOW_H__

#include <string>

#include "SuperaBase.h"
#include "FMWKInterface.h"
#include "ImageMetaMaker.h"
#include "ParamsImage2D.h"
#include "DataFormat/Image2D.h"

namespace larcv {

  /**
     \class ProcessBase
     User defined class SuperaLArFlow ... these comments are used to generate
     doxygen documentation!
  */
  class SuperaLArFlow : public SuperaBase,
                        public supera::ParamsImage2D,
                        public supera::ImageMetaMaker {

  public:
    
    /// Default constructor
    SuperaLArFlow(const std::string name="SuperaLArFlow");
    
    /// Default destructor
    ~SuperaLArFlow(){}

    void configure(const PSet&);

    void initialize();

    bool process(IOManager& mgr);

    void finalize();

  private:

    unsigned short _origin;
    std::string    _chstatus_producer;

  };

  /**
     \class larcv::SuperaLArFlowFactory
     \brief A concrete factory class for larcv::SuperaLArFlow
  */
  class SuperaLArFlowProcessFactory : public ProcessFactoryBase {
  public:
    /// ctor
    SuperaLArFlowProcessFactory() { ProcessFactory::get().add_factory("SuperaLArFlow",this); }
    /// dtor
    ~SuperaLArFlowProcessFactory() {}
    /// creation method
    ProcessBase* create(const std::string instance_name) { return new SuperaLArFlow(instance_name); }
  };

}
#endif
/** @} */ // end of doxygen group 

