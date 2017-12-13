/**
 * \file flash.h
 *
 * \ingroup GBDTprotonsPackage
 *
 * \brief Class def header for a class flash
 *
 * @author erezcohen
 */

/** \addtogroup GBDTprotonsPackage
 
 @{*/
#ifndef flash_H
#define flash_H

#include "Rtypes.h"
#include <iostream>



class flash{
    
public:
    
    /// Default constructor
    flash () = default;
    flash (Float_t ftime, Float_t ftimewidth,           
           Float_t fZcenter , Float_t fZwidth,
           Float_t fYcenter , Float_t fYwidth,
           Float_t ftotalPE);

    
    void                     Print () const;
    
    /// SETters
    void                   SetTime (Float_t _t)     {time = _t;};
    void              SetTimeWidth (Float_t _tw)    {timewidth = _tw;};
    void                SetZcenter (Float_t _z)     {Zcenter = _z;};
    void                 SetZwidth (Float_t _zw)    {Zwidth = _zw;};
    void                SetYcenter (Float_t _y)     {Ycenter = _y;};
    void                 SetYwidth (Float_t _yw)    {Ywidth = _yw;};
    void                SetTotalPE (Float_t _totpe) {totalPE = _totpe;};
    
    
    // GETters
    Float_t              GetTime () const  {return time;};
    Float_t         GetTimeWidth () const  {return timewidth;};
    Float_t           GetZcenter () const  {return Zcenter;};
    Float_t            GetZwidth () const  {return Zwidth;};
    Float_t           GetYcenter () const  {return Ycenter;};
    Float_t            GetYwidth () const  {return Ywidth;};
    Float_t           GetTotalPE () const  {return totalPE;};

    
    
    
private:
    
    // Float_t features
    Float_t time=-1, timewidth=-1;
    Float_t Zcenter=-1, Zwidth=-1, Ycenter=-1, Ywidth=-1, totalPE=-1;

    
    
};
#endif
/** @} */ // end of doxygen group

