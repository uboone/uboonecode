/**
 * \file pairVertex.h
 *
 * \ingroup CCQEPackage
 * 
 * \brief Class def header for a class pairVertex
 *
 * @author erezcohen
 */

/** \addtogroup CCQEPackage

    @{*/
#ifndef PAIRVERTEX_H
#define PAIRVERTEX_H

#include <iostream>
#include "PandoraNuTrack.h"
#include "hit.h"
#include "box.h"
#include "flash.h"
#include "GENIEinteraction.h"




/**
   \class pairVertex
   User defined class pairVertex ... these comments are used to generate
   doxygen documentation!
 */
class pairVertex {

public:

    
    pairVertex() = default;
    pairVertex (Int_t frun,Int_t fsubrun,Int_t fevent,Int_t fID);

    
   
    
    
    
    // functionallity
    bool                           MatchGenieInteraction ( std::vector<GENIEinteraction> , PandoraNuTrack );
    vector<size_t>                                sort_l (const std::vector<PandoraNuTrack> &v);
    bool                              SortTracksByLength ();
    vector<size_t>                             sort_pida (const std::vector<PandoraNuTrack> &v);
    bool                                SortTracksByPIDA ();
    
    void                             FixTracksDirections ();
    void                                        AddTrack (PandoraNuTrack ftrack);
    void                                      AddTrackID (Int_t ftrack_id)               {track_id.push_back(ftrack_id);};
    
    void                                          Print (bool DoPrintTracks=false) const;
    bool                                  IncludesTrack (Int_t ftrack_id);
    bool                                RemoveFarTracks (float max_mu_p_distance );
    
    std::vector<PandoraNuTrack>  CloseSemiContainedTracks ( std::vector<PandoraNuTrack> AllTracksInTheEvent , float fmax_distance=100 );
    std::vector<PandoraNuTrack>     RemoveTrackFromVector ( std::vector<PandoraNuTrack> AllTracks , PandoraNuTrack TrackToBeRemoved );
    std::vector<PandoraNuTrack>    RemoveTracksFromVector ( std::vector<PandoraNuTrack> AllTracks , std::vector<PandoraNuTrack> TracksToBeRemoved );
    
    void                                ReconstructBeam ();
    void                          ReconstructKinematics ();
    
    void                          AssociateHitsToTracks (std::vector<hit> hits);
    
    
    // SETters
    void                SetVertexID (Int_t fvertex_id)                              {vertex_id = fvertex_id;};
    void                     SetRSE (Int_t r, Int_t s, Int_t e)                     {run=r; subrun=s; event=e;};
    void                SetPosition (TVector3 fposition)                            {position = fposition;};
    void                 SetAs1mu1p ();
    void               SetAsCC1p0pi ();
    void              SetAsNon1mu1p ();
    void                SetAsCosmic ();
    void    SetTrueMuonProtonTracks (PandoraNuTrack , PandoraNuTrack );
    void          SetTrueMuonProton (PandoraNuTrack , PandoraNuTrack );
    void         SetTracksRelations ();
    bool         SetIsReconstructed (float fmax_mu_p_distance = 11 );
    void               SetGENIEinfo (GENIEinteraction fgenie_interaction)           {genie_interaction = fgenie_interaction; };
    void            SetClosestGENIE (GENIEinteraction fgenie_interaction)           {closest_genie_interaction = fgenie_interaction; };
    void       SetReconstructedInfo ();
    void            AssignMuonTrack (PandoraNuTrack ftrack)                         {AssignedMuonTrack = ftrack; };
    void          AssignProtonTrack (PandoraNuTrack ftrack)                         {AssignedProtonTrack = ftrack; };
    void    SetReconstructedMomenta (float PmuFromRange = 0, float PpFromRange = 0 );
    void   SetReconstructedFeatures (float PmuFromRange = 0, float PpFromRange = 0 );
    void         SetPlaneProjection (int plane , float _wire , float _time )        {vertex_wire[plane]=_wire; vertex_time[plane]=_time;};
    void            SetClosestFlash (flash _flash)                                  {ClosestFlash = _flash;};

    // GETters
    TString         GetTruthTopologyString () const {return TruthTopologyString;}
    
    bool                        GetIs1mu1p () const {return Is1mu1p;};
    bool       GetIsGENIECC_1p_200MeVc_0pi () const {return IsGENIECC_1p_200MeVc_0pi;};
    bool                     GetIsNon1mu1p () const {return IsNon1mu1p;};
    bool                       GetIsCosmic () const {return IsCosmic;};
    bool              GetIsVertexContained () const {return IsVertexContained;};
    bool       GetIs_mu_TrackReconstructed () const {return Is_mu_TrackReconstructed;};
    bool        GetIs_p_TrackReconstructed () const {return Is_p_TrackReconstructed;};
    bool          GetIsVertexReconstructed () const {return IsVertexReconstructed;};

    
    Int_t                           GetRun () const {return run;};
    Int_t                        GetSubrun () const {return subrun;};
    Int_t                         GetEvent () const {return event;};
    Int_t                      GetVertexID () const {return vertex_id;};
    Int_t                       GetNtracks () const {return (Int_t)tracks.size();};
 
    Int_t                    GetVertexWire (int plane) const {return vertex_wire[plane];};
    Int_t                    GetVertexTime (int plane) const {return vertex_time[plane];};
    
    float           GetAngleBetween2tracks () const; // return the angle between the muon and proton candidates, in degrees (!)
    float                        GetRecoEv () const {return reco_Pnu.E();};
    float                        GetRecoQ2 () const {return reco_Q2;};
    float                        GetRecoXb () const {return reco_Xb;};
    float                         GetRecoY () const {return reco_y;};
    float                        GetRecoW2 () const {return reco_W2;};
    float                        GetRecoPt () const {return (reco_Pmu + reco_Pp).Pt();};
    float                 GetReco_theta_pq () const {return reco_theta_pq;};
    float                 GetTruthDeltaPhi () const;
    float               GetDistanceToGENIE () const {return (genie_interaction.GetVertexPosition()-position).Mag();};
    float        GetDistanceToClosestGENIE () const {return (closest_genie_interaction.GetVertexPosition()-position).Mag();};
    
    // get the ratio of tracks-charge deposited to total-charge deposited
    // in a box of N(wires) x N(time-ticks) around the vertex in plane i=0,1,2
    // input: plane, N(wires) & N(time-ticks) for the box, hits in event
    // see docdb-10958
    float               GetRdQaroundVertex (int plane, int Nwires, int Nticks, std::vector<hit> hits) const ;
    float                   GetChargeInBox (int plane, std::vector<hit> hits, box VertexBox) const;
    
    
    
    Float_t                   GetDis2Flash (flash) const;
    Float_t            GetDis2ClosestFlash () const ;

    std::vector<float>    Get_delta_phi_ij () const {return delta_phi_ij;};
    std::vector<float>    Get_distances_ij () const {return distances_ij;};
    std::vector<float>  Get_delta_theta_ij () const {return delta_theta_ij;};

    
    
    TVector3                   GetPosition () const {return position;};
    
    TLorentzVector              GetRecoPnu () const {return reco_Pnu;};
    TLorentzVector              GetRecoPmu () const {return reco_Pmu;};
    TLorentzVector               GetRecoPp () const {return reco_Pp;};
    TLorentzVector            GetRecoPmiss () const {return (reco_Pp - reco_q);};
    
    
    PandoraNuTrack        GetShortestTrack () const {return ShortestTrack;};
    PandoraNuTrack         GetLongestTrack () const {return LongestTrack;};
    PandoraNuTrack       GetSmallPIDaTrack () const {return SmallPIDaTrack;};
    PandoraNuTrack       GetLargePIDaTrack () const {return LargePIDaTrack;};
    PandoraNuTrack    GetAssignedMuonTrack () const {return AssignedMuonTrack;};
    PandoraNuTrack  GetAssignedProtonTrack () const {return AssignedProtonTrack;};
    
    std::vector<PandoraNuTrack>  GetTracks () const {return tracks;};
    
    GENIEinteraction          GetGENIEinfo () const {return genie_interaction;};
    GENIEinteraction       GetClosestGENIE () const {return closest_genie_interaction;};
    flash                  GetClosestFlash () const {return ClosestFlash;};

    std::vector<hit>           GetMuonHits (int plane) const {return hits_muon[plane];};
    std::vector<hit>         GetProtonHits (int plane) const {return hits_proton[plane];};
    
    
    // operators
    inline bool operator==(const pairVertex & v) {
        return  (vertex_id == v.GetVertexID());
    }
    inline bool operator!=(const pairVertex & v) {
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
    
    bool                Is1mu1p=false,    IsGENIECC_1p_200MeVc_0pi=false,   IsNon1mu1p=false,   IsCosmic=false;
    bool                IsVertexContained=false, Is_mu_TrackReconstructed=false, Is_p_TrackReconstructed=false, IsVertexReconstructed=false;

    Int_t               run=-1 , subrun=-1 , event=-1, vertex_id=-1;
    Int_t               Ntracks=-1;
    
    // location in each plane
    float               vertex_wire[3]={0,0,0} , vertex_time[3]={0,0,0}; // these are floating point numbers since they are projections
    
    // reconstructed features
    // calorimentric reconstruction:
    // Ev = Eµ + Tp + Sn + T(A-1)
    float               reco_mu_p_distance=-1;
    float               reco_BeamE=-1,   reco_theta_pq=-1, reco_Pp_3momentum=-1, reco_Pmu_3momentum=-1;
    float               reco_p_over_q=-1, reco_Q2=-1;
    float               reco_omega=-1;
    float               reco_Xb=-1, reco_y=-1, reco_W2=-1, reco_s=-1;
    float               reco_alpha_p=-1 , reco_alpha_q=-1 , reco_alpha_mu=-1, reco_alpha_miss=-1;
    // --- - - --- -- - -- -- -- -- --- -- - --- - -- - -- -- -- --- - -- - --- - - -- - -- -
    
    float               truth_alpha_q, truth_alpha_p, truth_alpha_mu, truth_alpha_miss;
    
    TVector3            position=TVector3();
    TVector3            reco_Pp_3vect=TVector3(), reco_Pmu_3vect=TVector3();
    
    // Tp + Eµ
    TLorentzVector      reco_Pnu=TLorentzVector(-1,-1,-1,-1),  reco_Pp=TLorentzVector(-1,-1,-1,-1);
    TLorentzVector      reco_Pmu=TLorentzVector(-1,-1,-1,-1),  reco_q=TLorentzVector(-1,-1,-1,-1);
    TLorentzVector      reco_n_miss=TLorentzVector(-1,-1,-1,-1);
    
    PandoraNuTrack      muonTrueTrack,  protonTrueTrack;
    PandoraNuTrack      ShortestTrack,  LongestTrack;
    PandoraNuTrack      LargePIDaTrack, SmallPIDaTrack;
    PandoraNuTrack      AssignedMuonTrack, AssignedProtonTrack;
    
    GENIEinteraction    genie_interaction=GENIEinteraction();
    GENIEinteraction    closest_genie_interaction=GENIEinteraction();

    std::vector<Int_t>  track_id;

    std::vector <std::vector<float> >   tracks_distances;
    std::vector <std::vector<float> >   tracks_delta_phi;
    std::vector <std::vector<float> >   tracks_delta_theta;
    std::vector<float>                  tracks_dis_from_vertex, delta_phi_ij,    distances_ij , delta_theta_ij;
    std::vector<PandoraNuTrack>         tracks, tracks_lengthsorted,  tracks_pidasorted ;
    
    std::vector<hit>    hits_muon[3], hits_proton[3]; // in 3 wire planes
    flash               ClosestFlash=flash();
    
};

#endif
/** @} */ // end of doxygen group 

