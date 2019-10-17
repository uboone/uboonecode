/**
 * \file GENIEinteraction.h
 *
 * \ingroup GBDTprotonsPackage
 *
 * \brief Class def header for a class GENIEinteraction
 *
 * @author erezcohen
 */

/** \addtogroup GBDTprotonsPackage
 
 @{*/
#ifndef GENIEinteraction_H
#define GENIEinteraction_H

#include "PandoraNuTrack.h"
#include "TLorentzVector.h"
#include "TVector3.h"

/**
 \class GENIEinteraction
 User defined class GENIEinteraction ... these comments are used to generate
 doxygen documentation!
 */
using namespace std;






class GENIEinteraction{
    
public:
    
    /// Default constructor
    GENIEinteraction () = default;
    GENIEinteraction (Int_t frun, Int_t fsubrun, Int_t fevent, Int_t fmcevent_id);
    
    
    
    
    
    
    // running

    vector<size_t> sort_by_momentum_magnitude(const vector<TVector3> &v);
    bool            SortNucleons ();
    bool       ComputeKinematics ();
    bool        ComputePmissPrec ();
    bool CheckContainement (float max_y = 116.5,
                            float min_z = 0, float max_z = 1037,
                            float min_x = 0, float max_x = 256){
        // start with contained = true
        IsVertexContained=true;
        // and then check, to see if contained = false
        if( ( vertex_position.x() < min_x )    | ( vertex_position.x() > max_x ) )    IsVertexContained=false;
        if( ( vertex_position.y() < -max_y )   | ( vertex_position.y() > max_y ) )    IsVertexContained=false;
        if( ( vertex_position.z() < min_z )    | ( vertex_position.z() > max_z ) )    IsVertexContained=false;
        return  IsVertexContained;
    }
    void                   Print (bool DoPrintTracks=false, bool DoPrintAssignedMuonProtonTracks=false) const;
    void                AddTrack (PandoraNuTrack ftrack);
    bool           IncludesTrack (Int_t ftrackID=0) const;    // ask if the track (ftrackID) is in the list of tracks belonging to this genie interaction

   
    
    
    
    // SETters
    void         SetTruthTopology (); // set TRUTH topology DEFINITIONs
    void SetReconstructedTopology (); // set RECONSTRUCTED topology DEFINITIONs
    void                   SetRSE (int frun , int fsubrun , int fevent)   {run=frun; subrun=fsubrun;event=fevent;};
    void                  SetCCNC (int fccnc)                             {ccnc = fccnc;};
    void                  SetMode (int fmode)                             {mode = fmode;};
    void       SetVertexContained (bool fcontained)                       {IsVertexContained = fcontained;};
    void            SetNprimaries (Int_t n)                               {Nprimaries = n;};
    void                 SetNuPDG (Int_t pdg)                             {nuPDG = pdg;};
    void            SetNuMomentum (TLorentzVector fnu)                    {nu = fnu;};
    void        SetLeptonMomentum (TLorentzVector fmuon)                  {muon = fmuon;}; // if the neutrino is v(e) this will be the e!
    void      SetMomentumTransfer (); // should only be called after SetNuMomentum() and SetLeptonMomentum()
    void        SetVertexPosition (TVector3 fpos);
    void            SetKinematics (Float_t fQ2, Float_t fXb, Float_t fW, Float_t fy, Int_t fccnc, Int_t fmode);
    bool               AddPrimary ( Int_t fpdg                  // pdg code
                                   ,TLorentzVector fmomentum    // 4-momentum
                                   ,Int_t fstatus_code          // status code
                                   ,Int_t fmother               // mother
                                   ,std::string fprocess        // process
    );
    
    void      SetIsInActiveVolume (bool fIsInActiveVolume=false)        {IsInActiveVolume = fIsInActiveVolume;};
    void                SetIsCCQE (bool fIsCCQE=false)                  {IsCCQE = fIsCCQE;};
    void   SetIs_mu_Reconstructed (bool fIs_mu_TrackReconstructed=false){Is_mu_TrackReconstructed = fIs_mu_TrackReconstructed;};
    void    SetIs_p_Reconstructed (bool fIs_p_TrackReconstructed=false) {Is_p_TrackReconstructed = fIs_p_TrackReconstructed;};
    void SetIsVertexReconstructed (bool fIsVertexReconstructed=false)   {IsVertexReconstructed = fIsVertexReconstructed;};
    void               SetIs1mu1p (bool fIs1mu1p=false)                 {Is1mu1p = fIs1mu1p;};
    
    void    SetReco_mu_p_distance ();
    void           SetNuIntInBeam (TVector3 fpos)                       {posNuIntInBeam = fpos;};
    
    
    
    // GETters
    bool                        GetVertexContained () const {return IsVertexContained;};
    bool                        GetIsCC_Np_200MeVc () const {return IsCC_Np_200MeVc;};
    bool                        GetIsCC_1p_200MeVc () const {return IsCC_1p_200MeVc;};
    bool                    GetIsCC_1p_200MeVc_0pi () const {return IsCC_1p_200MeVc_0pi;};
    
    bool               GetIs_mu_TrackReconstructed () const {return Is_mu_TrackReconstructed;};
    bool                        GetIs_mu_TrackInFV () const {return Is_mu_TrackInFV;};
    bool                GetIs_p_TrackReconstructed () const {return Is_p_TrackReconstructed;};
    bool                         GetIs_p_TrackInFV () const {return Is_p_TrackInFV;};
    bool                       GetIsInActiveVolume () const {return IsInActiveVolume;};
    bool                  GetIsVertexReconstructed () const {return IsVertexReconstructed;};
    bool                           GetIsVertexInFV () const {return IsVertexInFV;};
    bool                                GetIs1mu1p () const {return Is1mu1p;};
    bool                                 GetIsCCQE () const {return IsCCQE;};
    
    
    Int_t                                   GetRun () const {return run;};
    Int_t                                GetSubrun () const {return subrun;};
    Int_t                                 GetEvent () const {return event;};
    
    Int_t                             GetMCeventID () const {return mcevent_id;};
    Int_t                                  GetCCNC () const {return ccnc;};
    Int_t                                  GetMode () const {return mode;};
    Int_t                            GetNprimaries () const {return Nprimaries;};
    Int_t                             GetPrimaries () const {return Nprimaries;};
    Int_t                              GetNprotons () const {return protons.size();};
    Int_t                             GetNneutrons () const {return neutrons.size();};
    Int_t                                GetNpions () const {return Npi;};
    Int_t                                  GetNpi0 () const {return Npi_0;};

    Float_t                                  GetPt () const {return (protons.size()) ? (muon+protons.at(0)).Pt() : -1;};
    Float_t                           Get_theta_pq () const {return theta_pq;};
    Float_t                                  GetEv () const {return nu.E();};
    Float_t                                  Get_q () const {return q.P();};
    Float_t                              Get_omega () const {return q.E();};
    Float_t                                  GetQ2 () const {return Q2;};
    Float_t                                  GetXb () const {return Xb;};
    Float_t                                   GetY () const {return y;};
    Float_t                                  GetW2 () const {return W*W;};
    Float_t                  GetReco_mu_p_distance () const {return reco_mu_p_distance;};

    TVector3                     GetVertexPosition () const {return vertex_position;};
    TLorentzVector               GetLeptonMomentum () const {return muon;}; // if the neutrino is v(e) this will be the e!
    TLorentzVector             GetMomentumTransfer () const {return q;}; // if the neutrino is v(e) this will be the e!
    TLorentzVector                          GetPmu () const {return GetLeptonMomentum();};
    TLorentzVector                           GetPp () const {return ((protons.size()>0) ? protons.at(0) : TLorentzVector());};
    TLorentzVector                           GetPv () const {return nu;};
    
    std::vector<TLorentzVector>  GetProtonsMomenta () const {return protons;}; // when taking this, check if the vector is not empty...
    
    
    
    
    PandoraNuTrack                    Get_mu_track () const {return muonTrack;};
    PandoraNuTrack                     Get_p_track () const {return protonTrack;};
    std::vector<PandoraNuTrack>          GetTracks () const {return tracks;}; // when taking this, check if the vector is not empty...
    
    TVector3                        GetNuIntInBeam () const {return posNuIntInBeam;};

    
    
    
    

private:
    
    // booleans on the genie interaction
    // TRUTH topology DEFINITIONs
    bool                    IsCC_Np_200MeVc=false; // an interaction with at least 1 muon and N protons > 200 MeV/c
    bool                    IsCC_1p_200MeVc=false; // an interaction with at least 1 muon and exactly 1 proton > 200 MeV/c and no charged pions with momentum greater than 70 MeV/c, and no neutral pions in the final state
    bool                    IsCC_1p_200MeVc_0pi=false; // an interaction with at least 1 muon and 1 proton > 200 MeV/c and no pions, and no photons or electrons outside the nucleus
    bool                    IsCCQE=false; // is QE or not - genie's "mode" flag: QE=0
    // RECONSTRUCTED topology DEFINITIONs
    bool                    IsVertexContained=false, IsInActiveVolume=false; // should be practically the same
    bool                    Is_mu_TrackReconstructed=false, Is_p_TrackReconstructed=false;
    bool                    IsVertexReconstructed=false; // vertex reconstructed = muon reconstructed && proton reconsutructed
    bool                    Is1mu1p=false; // 1 muon + 1 proton reconstructed, and nothing else
    bool                    Is_mu_TrackInFV=false, Is_p_TrackInFV=false;
    bool                    IsVertexInFV=false;

    // Int_t
    Int_t                   nuPDG=-9999;
    Int_t                   ccnc=-9999; // CC=0/NC=1
    Int_t                   mode=-9999; // neutrino nucleus 0=Quasi-elastic or Elastic, 1=Resonant (RES), 2=DIS, 3=Coherent production, 10=MEC
    Int_t                   run=-9999, subrun=-9999, event=-9999, mcevent_id=-9999;
    Int_t                   Nprimaries=0;
    Int_t                   Nnu=0, Nnu_e=0, Nnu_mu=0;
    Int_t                   Np=0, Nn=0, Npi=0;
    Int_t                   Nmu=0, Nel=0, Ngamma=0;
    Int_t                   Nmu_minus=0, Nmu_plus=0;
    Int_t                   Npi_minus=0, Npi_plus=0, Npi_0=0;
    Int_t                   Ne_plus=0, Ne_minus=0;
    Int_t                   pos_wire_u=-9999, pos_wire_v=-9999, pos_wire_y=-9999;
    Int_t                   pos_time_u=-9999, pos_time_v=-9999, pos_time_y=-9999;
    
    // Flaot_t
    Float_t                 Xb=-9999 , Q2=-9999, W=-9999, y=-9999;
    Float_t                 theta_pq=-9999, p_over_q=-9999, Mmiss=-9999;
    Float_t                 reco_mu_p_distance=-9999;
    Float_t                 fastest_proton_momentum = -1;
    
    // TVector3
    TVector3                vertex_position=TVector3();
    TVector3                posNuIntInBeam=TVector3();
    
    // TLorentzVector
    TLorentzVector          momentum=TLorentzVector(), Plead=TLorentzVector() ;
    TLorentzVector          nu=TLorentzVector(-1,-1,-1,-1), muon=TLorentzVector()  , q=TLorentzVector() ;
    TLorentzVector          n_miss=TLorentzVector()  , Pcm=TLorentzVector()   , Prec=TLorentzVector();
    
    // std::vector-s
    std::vector<Int_t>      pdg;           //particle type (pdg) of the GENIE particle
    std::vector<Int_t>      trackID;       //trackID of the GENIE particle (different from the GEANT-assigned track ID)
    std::vector<Int_t>      ND;            //number of daughters of the GENIE particle
    std::vector<Int_t>      mother;        //mother trackID of the GENIE particle
    std::vector<Int_t>      status_code;   //particle status code of the GENIE particle
    
    std::vector<Float_t>    Eng;           //Energy of the GENIE particle in GeV
    std::vector<Float_t>    Px;            //Px of the GENIE particle in GeV
    std::vector<Float_t>    Py;            //Py of the GENIE particle in GeV
    std::vector<Float_t>    Pz;            //Pz of the GENIE particle in GeV
    std::vector<Float_t>    P;             //Magnitude of the momentum vector of the GENIE particle in GeV
    std::vector<Float_t>    mass;          //mass of the GENIE particle in GeV

    std::vector<TVector3>       p3vect  , n3vect, pion3vect;
    std::vector<TLorentzVector> protons , neutrons;

    
    std::vector<std::string>    process;   
    
    // PandoraNuTrack
    PandoraNuTrack              muonTrack=PandoraNuTrack(), protonTrack=PandoraNuTrack();
    std::vector<PandoraNuTrack> tracks; // pandoraNu tracks that are associated with the genie interacion

};
#endif
/** @} */ // end of doxygen group

