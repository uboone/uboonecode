/**
 * \file Box.h
 *
 * \ingroup GBDTprotonsPackage
 *
 * \brief Class def header for a class Box
 *
 * @author erezcohen
 */

/** \addtogroup GBDTprotonsPackage
 
 @{*/
#ifndef box_H
#define box_H

#include "Rtypes.h"
#include <iostream>



class box{
    
public:
    
    /// Default constructor
    box () = default;
    box (Int_t, Int_t, Int_t, Int_t);

    
    void                     Print () const;
    
    /// SETters
    void              SetStartWire (Int_t _wire)      {start_wire = _wire;};
    void                SetEndWire (Int_t _wire)      {end_wire = _wire;};
    void              SetStartTime (Int_t _time)      {start_time = _time;};
    void                SetEndTime (Int_t _time)      {end_time = _time;};
    
    
    // GETters
    Int_t              GetStartWire () const  {return start_wire;};
    Int_t                GetEndWire () const  {return end_wire;};
    Int_t              GetStartTime () const  {return start_time;};
    Int_t                GetEndTime () const  {return end_time;};

    
    
    
private:
    
    // Int_t features
    Int_t       start_wire=0 , start_time=0 , end_wire=0 , end_time=0;

    
    
};
#endif
/** @} */ // end of doxygen group

