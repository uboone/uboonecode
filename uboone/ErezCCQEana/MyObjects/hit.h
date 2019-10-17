/**
 * \file hit.h
 *
 * \ingroup GBDTprotonsPackage
 *
 * \brief Class def header for a class Hit
 *
 * @author erezcohen
 */

/** \addtogroup GBDTprotonsPackage
 
 @{*/
#ifndef hit_H
#define hit_H

#include "Rtypes.h"
#include "box.h"




class hit{
    
public:
    
    /// Default constructor
    hit () = default;
    hit (Int_t fplane, Int_t fwire, Int_t fid, float fpeakT, float fcharge);
    
    
    void    Print () const;
    
    /// SETters
    void               SetTrackKey (Int_t _key)      {hit_trkKey = _key;};
    
    // GETters
    Int_t              GetTrackKey ()   const {return hit_trkKey;};
    Int_t                 GetPlane ()   const {return hit_plane;};
    Int_t                  GetWire ()   const {return hit_wire;};
    Int_t                    GetID ()   const {return hit_id;};
    float                GetCharge ()   const {return hit_charge;};
    float                 GetPeakT ()   const {return hit_peakT;};

    bool                   InPlane (int fplane) const {return (hit_plane==fplane);};

    // check if hit is inside a box
    bool                     InBox (box fbox)   const;

//    // check if hit is inside a sphere of radius r [cm]
//    // localized around (cwire,ctime)
//    bool                  InSphere (float cwire, float ctime,
//                                    float r, // radius in [cm]
//                                    float f_timeticks2cm,
//                                    float f_wire2cm ) const ;

    
private:
    
    // Int_t features
    Int_t       hit_plane=-1,   hit_wire=-1, hit_id=-1,  hit_trkKey=-1;
    
    
    // Float_t features
    float       hit_peakT=-1,    hit_charge=-1;
    
};
#endif
/** @} */ // end of doxygen group

