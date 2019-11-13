/**
 * \file tripleVertex.h
 *
 * \ingroup CCQEPackage
 * 
 * \brief Class def header for a class tripleVertex
 *
 * @author erezcohen
 */

/** \addtogroup CCQEPackage

    @{*/
#ifndef tripleVertex_H
#define tripleVertex_H

#include <iostream>
#include "PandoraNuTrack.h"
#include "hit.h"
#include "box.h"
#include "flash.h"
#include "GENIEinteraction.h"




/**
   \class tripleVertex
   User defined class tripleVertex ... these comments are used to generate
   doxygen documentation!
 */
class tripleVertex {

public:

    
    tripleVertex() = default;
    tripleVertex (Int_t frun,Int_t fsubrun,Int_t fevent,Int_t fID);
   
    
    
    
    // functionallity
    vector<size_t>                                sort_l (const std::vector<PandoraNuTrack> &v);
    bool                              SortTracksByLength ();
    vector<size_t>                             sort_pida (const std::vector<PandoraNuTrack> &v);
    
    void                                        AddTrack (PandoraNuTrack ftrack);
    void                                      AddTrackID (Int_t ftrack_id)               {track_id.push_back(ftrack_id);};
    
    void                                          Print (bool DoPrintTracks=false, bool DoPrintFull=true) const;
    bool                                  IncludesTrack (Int_t ftrack_id);
    bool                                RemoveFarTracks (float max_mu_p_distance );
    std::vector<PandoraNuTrack>   RemoveTrackFromVector ( std::vector<PandoraNuTrack> AllTracks
                                                         , PandoraNuTrack TrackToBeRemoved );
    std::vector<PandoraNuTrack>  RemoveTracksFromVector ( std::vector<PandoraNuTrack> AllTracks
                                                         , std::vector<PandoraNuTrack> TracksToBeRemoved );
    
    // @brief check if the recostructed vertex is inside the TPC
    bool                                  IsVertexInTPC (float max_FV_y = 116.5,
                                                         float min_FV_z = 0, float max_FV_z = 1037,
                                                         float min_FV_x = 0, float max_FV_x = 257)  const{
        if( ( position.x() < min_FV_x )    | ( position.x() > max_FV_x ) )    return false;
        if( ( position.y() < -max_FV_y )   | ( position.y() > max_FV_y ) )    return false;
        if( ( position.z() < min_FV_z )    | ( position.z() > max_FV_z ) )    return false;
        return true;
    };
    void                             CheckIsVertexInTPC () { IsRecoVertexInTPC = IsVertexInTPC(); };
    
    
    
    
    
    
    // SETters
    void                SetVertexID (Int_t fvertex_id)                              {vertex_id = fvertex_id;};
    void                     SetRSE (Int_t r, Int_t s, Int_t e)                     {run=r; subrun=s; event=e;};
    void                SetPosition (TVector3 fposition)                            {position = fposition;};
    void         SetTracksRelations ();

    void         SetPlaneProjection (int plane , float _wire , float _time )        {vertex_wire[plane]=_wire; vertex_time[plane]=_time;};
    void       SetIsVertexContained (bool fIsCon)                                   {IsVertexContained = fIsCon;};


    
    
    
    // GETters
    TString         GetTruthTopologyString () const {return TruthTopologyString;}
    
    bool              GetIsVertexContained () const {return IsVertexContained;};
    bool          GetIsVertexReconstructed () const {return IsVertexReconstructed;};
    bool              GetIsRecoVertexInTPC () const {return IsRecoVertexInTPC;};
    
    Int_t                           GetRun () const {return run;};
    Int_t                        GetSubrun () const {return subrun;};
    Int_t                         GetEvent () const {return event;};
    Int_t                      GetVertexID () const {return vertex_id;};
    Int_t                       GetNtracks () const {return (Int_t)tracks.size();};
 
    Int_t                    GetVertexWire (int plane) const {return vertex_wire[plane];};
    Int_t                    GetVertexTime (int plane) const {return vertex_time[plane];};
    

    std::vector<float>    Get_delta_phi_ij () const {return delta_phi_ij;};
    std::vector<float>    Get_distances_ij () const {return distances_ij;};
    std::vector<float>  Get_delta_theta_ij () const {return delta_theta_ij;};

    
    TVector3                   GetPosition () const {return position;};
    
    PandoraNuTrack        GetShortestTrack () const {return ShortestTrack;};
    PandoraNuTrack         GetLongestTrack () const {return LongestTrack;};
    
    std::vector<PandoraNuTrack>  GetTracks () const {return tracks;};
    
    // operators
    inline bool operator==(const tripleVertex & v) {
        return  (vertex_id == v.GetVertexID());
    }
    inline bool operator!=(const tripleVertex & v) {
        return  (vertex_id != v.GetVertexID());
    }
    
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    // debug
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
    Int_t debug=0;
    void Debug (Int_t verobosity_level, std::string text){
        if ( debug > verobosity_level ) cout << text << endl;
    }
    // ---- - - -- -- - -- -- -- -- --- - - - - -- --- - - - --- -- - -
   
private:
    
    
    // -------------------------------------------------------
    
    // variables
    TString             TruthTopologyString="unknown truth topology";
    
    bool                IsVertexContained=false;
    bool                IsVertexReconstructed=false;
    bool                IsRecoVertexInTPC=false;

    Int_t               run=-1 , subrun=-1 , event=-1, vertex_id=-1;
    Int_t               Ntracks=-1;
    
    // location in each plane
    float               vertex_wire[3]={0,0,0} , vertex_time[3]={0,0,0}; // these are floating point numbers since they are projections
    
    TVector3            position=TVector3();
    
    PandoraNuTrack      ShortestTrack,  LongestTrack;

    std::vector<Int_t>  track_id;

    std::vector <std::vector<float> >   tracks_distances;
    std::vector <std::vector<float> >   tracks_delta_phi;
    std::vector <std::vector<float> >   tracks_delta_theta;
    std::vector<float>                  tracks_dis_from_vertex, delta_phi_ij,    distances_ij , delta_theta_ij;
    std::vector<PandoraNuTrack>         tracks, tracks_lengthsorted ;
    
};

#endif
/** @} */ // end of doxygen group 

