/**
 * \file PandoraNuTrack.h
 *
 * \ingroup GBDTprotonsPackage
 *
 * \brief Class def header for a class PandoraNuTrack
 *
 * @author erezcohen
 */

/** \addtogroup GBDTprotonsPackage
 
 @{*/
#ifndef PANDORANUTRACK_H
#define PANDORANUTRACK_H


#include "Rtypes.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "box.h"
#include "flash.h"
#include <iostream>
#include <iomanip>

using namespace std;

// prints....
#define EndEventBlock() std::cout << "\033[32m"<< "....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......" << "\033[0m"<< endl
#define PrintLine() std::cout << "--------------------------------------------------------------" << "\033[0m" << std::endl
#define PrintXLine() std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << "\033[0m" << std::endl
#define PrintXLine() std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << "\033[0m" << std::endl
#define PrintHeader(header) std::cout << "\033[33m" << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n" << (header) << "\n" << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" << std::endl



#define SHOW(a) std::cout << setprecision(2) << fixed << #a << ": " << (a) << std::endl;
#define SHOW2(a,b) std::cout <<"\033[34m"<<#a<<": "<<(a)<<"," << #b <<": "<<(b)<< "\033[0m"<< std::endl;
#define SHOW3(a,b,c) std::cout <<"\033[36m"<<#a<<": "<<(a)<<"," << #b <<": "<<(b)<<","<<#c<<": "<<(c)<< "\033[0m"<< std::endl;
#define SHOW4(a,b,c,d) std::cout <<"\033[31m"<<#a<<": "<<(a)<<"," << #b <<": "<<(b)<<","<<#c<<": "<<(c)<<","<<#d<<": "<<(d)<< "\033[0m"<< std::endl;
#define SHOWstdVector(v){ if (v.size()<1) {std::cout << #v << " is empty" << std::endl;} else {std::cout << #v << "( " << v.size() << " entries):\t"; for (auto it:v) std::cout << it << "\t"; std::cout << endl;}}
#define SHOWTVector3(v){ std::cout << #v << ": (" << v.X() << "," << v.Y() << "," << v.Z() << ")" << endl;}
#define SHOWTLorentzVector(v) std::cout << #v << ": " << "\t(" << setprecision(2) << fixed << v.Px() << ","  << v.Py() << "," << v.Pz()  << "," << v.E() << ")" << ", P = " << v.P() << ", M = " << v.M() << std::endl

#define PrintPhys(a,units) std::cout  << setprecision(2) << fixed << #a << ": " << (a) <<  " " << (units) << std::endl


#define r2d TMath::RadToDeg()
#define d2r TMath::DegToRad()
#define PI TMath::Pi()



/**
 \class PandoraNuTrack
 User defined class PandoraNuTrack ... these comments are used to generate
 doxygen documentation!
 */





class PandoraNuTrack{
    
public:
    
    /// Default constructor
    PandoraNuTrack () = default;
    PandoraNuTrack (Int_t frun, Int_t fsubrun, Int_t fevent
                    ,Int_t ftrack_id
                    ,Float_t flength
                    ,TVector3 fstart_pos, TVector3 fend_pos );
    
    
    
    // -- - - -- -- -- -- ----- - - -- -- -- -- ----- - - -- -- -- -- ---
    /// SETters
    // -- - - -- -- -- -- ----- - - -- -- -- -- ----- - - -- -- -- -- ---
    void      SetTruthProcess (std::string _process)        {truth_process = _process;};
    void            SetOrigin (std::string _origin)         {truth_origin = _origin;};

    void SetCC_1p_200MeVc_0pi ()                            {IsGENIECC_1p_200MeVc_0pi=true;};
    
    void               SetRun (Int_t _run)                  {run = _run;};
    void            SetSubrun (Int_t _subrun)               {subrun = _subrun;};
    void             SetEvent (Int_t _event)                {event = _event;};
    void         SetMCeventID (Int_t fmcevent_id)           {mcevent_id = fmcevent_id;};
    void              SetCCNC (Int_t fccnc)                 {truth_ccnc = fccnc;};
    void              SetMode (Int_t fmode)                 {truth_mode = fmode;};
    void         SetMCpdgCode (Int_t _mcpdg)                {MCpdgCode = _mcpdg;};
    void         SetBestPlane (Int_t best)                  {BestPlane = best;};
    void          SetMaxNHits (Int_t maxnhits)              {MaxNHits = maxnhits;};
    void       SetTruthMother (Int_t _mother)               {truth_mother = _mother;};
    void   SetTruthMCParticle (Int_t _mcparticle)           {truth_mcparticle = _mcparticle;};
    void    SetCaloKEPerPlane (Int_t plane , Float_t ke )   {CaloKEPerPlane[plane] = ke;};
    void      SetPIDaPerPlane (Int_t plane , Float_t pida ) {PIDaPerPlane[plane] = pida;};
    
    
    void            SetLength (Float_t l)                   {length = l;};
    void      SetPandoraTheta (Float_t ftheta)              {pandora_theta = ftheta;};
    void        SetPandoraPhi (Float_t fphi)                {pandora_phi = fphi;};
    void   SetMomentumMCS_fwd (Float_t fmomentum)           {momentum_MSCMu_fwd = fmomentum; momentum_MSCMu=fmomentum;};
    void   SetMomentumMCS_bwd (Float_t fmomentum)           {momentum_MSCMu_bwd = fmomentum;};

    
    
    void       SetTruthLength ()                            {truth_length = (truth_start_pos-truth_end_pos).Mag();};
    void       SetTruthPurity (Float_t fpur)                {truth_purity = fpur;};

    void          SetStartPos (TVector3 pos)                {start_pos = pos;};
    void            SetEndPos (TVector3 pos)                {end_pos = pos;};
    void      SetRecDirection ()                            {rec_dir = end_pos-start_pos;};
    void     SetTruthStartPos (TVector3 pos)                {truth_start_pos = pos;};
    void       SetTruthEndPos (TVector3 pos)                {truth_end_pos = pos;};
    void    SetTruthDirection ()                            {truth_dir = truth_end_pos-truth_start_pos;};
    void              SetPIDa ()                            {PIDa = PIDaPerPlane[BestPlane]; };
    void     SetTruthMomentum (TLorentzVector fmomentum)    {truth_momentum=fmomentum; };
    
    void     SetStartEndPlane (Int_t plane ,
                               Int_t start_wire, Int_t start_time ,
                               Int_t end_wire, Int_t end_time );
    void      SetClosestFlash (flash _flash)                {ClosestFlash = _flash;};
    
    // calibrated PIDa
    void      SetMaxNHitsCali (Int_t maxnhits)              {MaxNHitsCali = maxnhits;};
    void     SetBestPlaneCali (Int_t best)                  {BestPlaneCali = best;};
    void  SetPIDaCaliPerPlane (Int_t plane , Float_t pida ) {PIDaCaliPerPlane[plane] = pida;};
    void          SetPIDaCali ()                            {PIDaCali = PIDaCaliPerPlane[BestPlaneCali]; };


    
    
    // completeness of the track MC-truth matching
    void    SetMaxdQinTruthMatchedHits (float fmax)         {max_dQinTruthMatchedHits=fmax;};
    void        SetdQinAllHits (float fdQ)                  {dQinAllHits = fdQ;};
    // dE/dx
    void               SetdEdx (int plane, std::vector<float> fdEdx)    {dEdx[plane]=fdEdx;};
    void           FilldEdxHit (int plane, float fdEdx_hit)             {dEdx[plane].push_back(fdEdx_hit);};
    void               SetdQdx (int plane, std::vector<float> fdQdx)    {dQdx[plane]=fdQdx;};
    void           FilldQdxHit (int plane, float fdQdx_hit)             {dQdx[plane].push_back(fdQdx_hit);};
    void           SetResRange (int plane, std::vector<float> fResRange){ResRange[plane]=fResRange;};
    void          FillResRange (int plane, float fResRange_hit)         {ResRange[plane].push_back(fResRange_hit);};

    // truncated dE/dx
    void          SetdEdxTrunc (int plane, std::vector<float> fdEdx)    {dEdxTrunc[plane]=fdEdx;};
    void          SetdQdxTrunc (int plane, std::vector<float> fdQdx)    {dQdxTrunc[plane]=fdQdx;};

    // PIDa from PandoraNu
    void       SetPandoraNuPID (int plane,
                                int fPdg , float fMinChi2, float fChi2Proton, float fChi2Kaon, float fChi2Pion, float fChi2Muon, float fPIDA )
    {
        PandoraNuPID_Pdg[plane]         = fPdg;
        PandoraNuPID_MinChi2[plane]     = fMinChi2;
        PandoraNuPID_Chi2Proton[plane]  = fChi2Proton;
        PandoraNuPID_Chi2Kaon[plane]    = fChi2Kaon;
        PandoraNuPID_Chi2Pion[plane]    = fChi2Pion;
        PandoraNuPID_Chi2Muon[plane]    = fChi2Muon;
        PandoraNuPID_PIDA[plane]        = fPIDA;
    };
    // PIDaCali from PandoraNu
    void   SetPandoraNuCaliPID (int plane,
                                int fPdg , float fMinChi2, float fChi2Proton, float fChi2Kaon, float fChi2Pion, float fChi2Muon, float fPIDA )
    {
        std::cout<<"plane "<<plane<<" fPdg "<<fPdg<<" fChi2Proton "<<fChi2Proton<<" fChi2Muon "<<fChi2Muon<<std::endl;
        PandoraNuCaliPID_Pdg[plane]         = fPdg;
        PandoraNuCaliPID_MinChi2[plane]     = fMinChi2;
        PandoraNuCaliPID_Chi2Proton[plane]  = fChi2Proton;
        PandoraNuCaliPID_Chi2Kaon[plane]    = fChi2Kaon;
        PandoraNuCaliPID_Chi2Pion[plane]    = fChi2Pion;
        PandoraNuCaliPID_Chi2Muon[plane]    = fChi2Muon;
        PandoraNuCaliPID_PIDA[plane]        = fPIDA;
    };
    
    // momentum from momentum calculator
    void SetMomCalc (
                     Float_t fPmuHypothesisCalc=0, // momentum for muon hypothesis
                     Float_t fPpHypothesisCalc=0 // momentum for proton hypothesis
    )
    {
        PmuHypothesisCalc = fPmuHypothesisCalc;
        PpHypothesisCalc = fPpHypothesisCalc;
    };


    
    
    
    // -- - - -- -- -- -- ----- - - -- -- -- -- ----- - - -- -- -- -- ---
    // GETters
    // -- - - -- -- -- -- ----- - - -- -- -- -- ----- - - -- -- -- -- ---

    std::string     GetTruthProcess () const {return truth_process;};
    std::string           GetOrigin () const {return truth_origin;};    // "unknown origin" / "beam neutrino" / "cosmic ray"

    
    Int_t                    GetRun () const {return run;};
    Int_t                 GetSubrun () const {return subrun;};
    Int_t                  GetEvent () const {return event;};
    Int_t                GetTrackID () const {return track_id;};
    Int_t              GetMCeventID () const {return mcevent_id;};
    Int_t            GetTruthMother () const {return truth_mother;};
    Int_t        GetTruthMCParticle () const {return truth_mcparticle;};
    
    Int_t              GetMCpdgCode () const {return MCpdgCode;};
    Int_t                   GetCCNC () const {return truth_ccnc;};
    Int_t              GetBestPlane () const {return BestPlane;};
    Int_t               GetMaxNHits () const {return MaxNHits;};

    Int_t              GetStartWire (int plane) const;
    Int_t              GetStartTime (int plane) const;
    Int_t                GetEndWire (int plane) const;
    Int_t                GetEndTime (int plane) const;

    
    Float_t               GetLength () const {return length;};
    Float_t                GetTheta () const {return rec_dir.Theta();};
    Float_t                  GetPhi () const {return rec_dir.Phi();};
    Float_t         GetPandoraTheta () const {return pandora_theta;};
    Float_t           GetPandoraPhi () const {return pandora_phi;};
    Float_t      GetMomentumMCS_fwd () const {return momentum_MSCMu_fwd;};
    Float_t      GetMomentumMCS_bwd () const {return momentum_MSCMu_bwd;};
    Float_t          GetMomentumMCS () const {return momentum_MSCMu;};

    Float_t          GetTruthLength () const {return truth_length;};
    Float_t           GetTruthTheta () const {return truth_dir.Theta();};
    Float_t             GetTruthPhi () const {return truth_dir.Phi();};
    Float_t          GetTruthPurity () const {return truth_purity;};

    Float_t       GetCaloKEPerPlane ( Int_t plane ) const { return CaloKEPerPlane[plane];};
    Float_t         GetPIDaPerPlane ( Int_t plane ) const { return PIDaPerPlane[plane];};
    Float_t                 GetPIDa ()              const { return PIDa; };
    Float_t     GetPIDaCaliPerPlane ( Int_t plane ) const { return PIDaCaliPerPlane[plane];};
    Float_t             GetPIDaCali ()              const { return PIDaCali; };
    
    // pandoraNu features
    Float_t                   GetPID_PIDA ( Int_t plane )   const { return PandoraNuPID_PIDA[plane];};
    Float_t                   GetPID_PIDA ()                const { return PandoraNuPID_PIDA[BestPlane];}; // best plane
    Float_t               GetCaliPID_PIDA ( Int_t plane )   const { return PandoraNuCaliPID_PIDA[plane];};
    Float_t               GetCaliPID_PIDA ()                const { return PandoraNuCaliPID_PIDA[BestPlane];}; // best plane
    Float_t         GetCaliPID_Chi2Proton ( Int_t plane )   const { return PandoraNuCaliPID_Chi2Proton[plane];};
    Float_t           GetCaliPID_Chi2Muon ( Int_t plane )   const { return PandoraNuCaliPID_Chi2Muon[plane];};

    
    Float_t                   GetDis2Flash (flash)  const;
    Float_t            GetDis2ClosestFlash ()       const ;
    Float_t     GetMaxdQinTruthMatchedHits ()       const {return max_dQinTruthMatchedHits;};
    Float_t                 GetdQinAllHits ()       const {return dQinAllHits;};
    Float_t   GetRatiodQinTruthMatchedHits ()       const
    {
        if (fabs(dQinAllHits)>0.01) {
        return (max_dQinTruthMatchedHits/dQinAllHits);
        } else {
            return -9999.;
        }
    }
    
    
    TVector3            GetStartPos () const {return start_pos;};
    TVector3              GetEndPos () const {return end_pos;};
    TVector3        GetRecDirection () const {return rec_dir;};

    TVector3       GetTruthStartPos () const {return truth_start_pos;};
    TVector3         GetTruthEndPos () const {return truth_end_pos;};
    TVector3      GetTruthDirection () const {return truth_dir;};
    
    TLorentzVector GetTruthMomentum () const {return truth_momentum; };
    
    flash           GetClosestFlash () const {return ClosestFlash;};
    
    // dE/dx
    std::vector<float>      GetdEdx (int plane) const {return dEdx[plane];};
    std::vector<float>      GetdQdx (int plane) const {return dQdx[plane];};
    std::vector<float>  GetResRange (int plane) const {return ResRange[plane];};
    // truncated dE/dx
    std::vector<float> GetdEdxTrunc (int plane) const {return dEdxTrunc[plane];};
    std::vector<float> GetdQdxTrunc (int plane) const {return dQdxTrunc[plane];};
    
    // momentum for muon hypothesis from calculator
    Float_t     GetPmuHypothesisCalc () const { return PmuHypothesisCalc;};
    Float_t      GetPpHypothesisCalc () const { return PpHypothesisCalc;};

    
    
    
    // functionallity
    void            FlipTrack ();
    void           CreateROIs ();
    void                Print (bool DoPrintPandoraNuFeatures = true , bool DoPrintAllPIDa = true ) const;

    
    // operators
    bool              AskIfCC1p0pi () const {return IsGENIECC_1p_200MeVc_0pi;};
    bool IsTrackContainedSoft (float max_FV_y = 115,
                              float min_FV_z = 5, float max_FV_z = 1045,
                              float min_FV_x = 3, float max_FV_x = 257)  const{
        if( ( start_pos.x() < min_FV_x )    | ( start_pos.x() > max_FV_x ) )    return false;
        if( ( start_pos.y() < -max_FV_y )   | ( start_pos.y() > max_FV_y ) )    return false;
        if( ( start_pos.z() < min_FV_z )    | ( start_pos.z() > max_FV_z ) )    return false;
        if( ( end_pos.x() < min_FV_x )      | ( end_pos.x() > max_FV_x ) )      return false;
        if( ( end_pos.y() < -max_FV_y )     | ( end_pos.y() > max_FV_y ) )      return false;
        if( ( end_pos.z() < min_FV_z )      | ( end_pos.z() > max_FV_z ) )      return false;
        return true;
    };
    bool IsTrackInFV (float max_FV_y = 110,
                      float min_FV_z = 5, float max_FV_z = 1037,
                      float min_FV_x = 3, float max_FV_x = 250)  const{
        if( ( start_pos.x() < min_FV_x )    | ( start_pos.x() > max_FV_x ) )    return false;
        if( ( start_pos.y() < -max_FV_y )   | ( start_pos.y() > max_FV_y ) )    return false;
        if( ( start_pos.z() < min_FV_z )    | ( start_pos.z() > max_FV_z ) )    return false;
        if( ( end_pos.x() < min_FV_x )      | ( end_pos.x() > max_FV_x ) )      return false;
        if( ( end_pos.y() < -max_FV_y )     | ( end_pos.y() > max_FV_y ) )      return false;
        if( ( end_pos.z() < min_FV_z )      | ( end_pos.z() > max_FV_z ) )      return false;
        return true;
    };
    
    
    
    
    Float_t           DistanceFromPoint ( TVector3 position , std::string * StartOrEnd=nullptr );
    Float_t ClosestDistanceToOtherTrack ( PandoraNuTrack other_track , std::string * StartOrEnd=nullptr );
    Float_t    ClosestLine2LineDistance ( PandoraNuTrack other_track
                                                      , TVector3 * QuasiIntersectionPoint=nullptr );
    Float_t   ClosestLine2PointDistance ( TVector3 x0 );
    
    inline bool operator==(const PandoraNuTrack & t) {
        return  (run==t.GetRun() && subrun==t.GetSubrun() && event==t.GetEvent() && track_id==t.GetTrackID());
    }
    inline bool operator!=(const PandoraNuTrack & t) {
        return  (run!=t.GetRun() || subrun!=t.GetSubrun() || event!=t.GetEvent() || track_id!=t.GetTrackID());
    }

    
private:

    // std::string
    std::string truth_process="unknown process";
    std::string truth_origin="unknown origin";
    
    // bool
    bool        IsGENIECC_1p_200MeVc_0pi=false;
    bool        IsFlipped=false;
    // Int_t
    Int_t       run=0, subrun=0, event=0, track_id=-9999;
    Int_t       truth_ccnc=0, truth_mode=0;
    Int_t       MCpdgCode=-9999;
    Int_t       start_wire_u=0, start_wire_v=0, start_wire_y=0;
    Int_t       start_time_u=0, start_time_v=0, start_time_y=0;
    Int_t       end_wire_u=0,   end_wire_v=0,   end_wire_y=0;
    Int_t       end_time_u=0,   end_time_v=0,   end_time_y=0;
    Int_t       BestPlane=2,    BestPlaneCali=2;
    Int_t       MaxNHits=0,     MaxNHitsCali=0;
    Int_t       mcevent_id=-1;
    Int_t       truth_mother=-1;
    Int_t       truth_mcparticle=-1;
    
    // Float_t
    Float_t     length=0, pandora_theta=0, pandora_phi=0;
    Float_t     CaloKEPerPlane[3]={0,0,0};
    Float_t     PIDaPerPlane[3]={0,0,0};
    Float_t     PIDa=-1;
    Float_t     PIDaCaliPerPlane[3]={0,0,0};
    Float_t     PIDaCali=-1;
    Float_t     momentum_MSCMu_fwd=0, momentum_MSCMu_bwd=0, momentum_MSCMu=0;

    // pandoraNu objects
    Float_t     PandoraNuPID_Pdg[3]={-9999,-9999,-9999};
    Float_t     PandoraNuPID_MinChi2[3]={-9999,-9999,-9999};
    Float_t     PandoraNuPID_Chi2Proton[3]={-9999,-9999,-9999};
    Float_t     PandoraNuPID_Chi2Kaon[3]={-9999,-9999,-9999};
    Float_t     PandoraNuPID_Chi2Pion[3]={-9999,-9999,-9999};
    Float_t     PandoraNuPID_Chi2Muon[3]={-9999,-9999,-9999};
    Float_t     PandoraNuPID_PIDA[3]={-9999,-9999,-9999};
    Float_t     PandoraNuCaliPID_Pdg[3]={-9999,-9999,-9999};
    Float_t     PandoraNuCaliPID_MinChi2[3]={-9999,-9999,-9999};
    Float_t     PandoraNuCaliPID_Chi2Proton[3]={-9999,-9999,-9999};
    Float_t     PandoraNuCaliPID_Chi2Kaon[3]={-9999,-9999,-9999};
    Float_t     PandoraNuCaliPID_Chi2Pion[3]={-9999,-9999,-9999};
    Float_t     PandoraNuCaliPID_Chi2Muon[3]={-9999,-9999,-9999};
    Float_t     PandoraNuCaliPID_PIDA[3]={-9999,-9999,-9999};
    
    // momentum calculator
    Float_t     PmuHypothesisCalc=0;
    Float_t     PpHypothesisCalc=0;

    
    // completeness of the track MC-truth matching
    Float_t     max_dQinTruthMatchedHits=-1, dQinAllHits=-1;

    // truth information - only valid for MC data
    Float_t     truth_length = -9999;
    Float_t     truth_purity = -9999;// maxe/tote

    // dE/dx in three planes
    std::vector<Float_t> ResRange[3];
    std::vector<Float_t> dEdx[3];
    std::vector<Float_t> dQdx[3];
    // truncated dE/dx
    std::vector<Float_t> dEdxTrunc[3];
    std::vector<Float_t> dQdxTrunc[3];
    

    
    
    // TVector3
    TVector3    start_pos=TVector3(), end_pos=TVector3(), rec_dir=TVector3();
    TVector3    truth_start_pos=TVector3(), truth_end_pos=TVector3(), truth_dir=TVector3();
    
    TLorentzVector truth_momentum=TLorentzVector();
    
    box         roi[3]={box(),box(),box()};
    
    flash       ClosestFlash=flash();
};
#endif
/** @} */ // end of doxygen group

