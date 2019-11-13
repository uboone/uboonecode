#ifndef PANDORANUTRACK_CXX
#define PANDORANUTRACK_CXX

#include "PandoraNuTrack.h"
#include "TFile.h"
#include "TProfile.h"
#include <stdexcept> // std::range_error
#include <vector>
#include <string>
#include <algorithm> // std::fill()
#include <functional>
#include <random>
#include <chrono>

#include <string>
#include <vector>
#include <memory> // std::unique_ptr<>
#include <cassert>

// // Framework includes
#include "cetlib/exception.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PandoraNuTrack::PandoraNuTrack( Int_t frun, Int_t fsubrun, Int_t fevent
                               ,Int_t ftrack_id
                               ,Float_t flength
                               ,TVector3 fstart_pos, TVector3 fend_pos):
run(frun),
subrun(fsubrun),
event(fevent),
track_id(ftrack_id),
length(flength),
start_pos(fstart_pos),
end_pos(fend_pos)
{
    SetRecDirection();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PandoraNuTrack::SetStartEndPlane(Int_t plane ,
                                      Int_t start_wire, Int_t start_time ,
                                      Int_t end_wire, Int_t end_time ){
    
    roi[plane] = box( start_wire , start_time , end_wire , end_time );
    switch (plane) {
        case 0:
            start_wire_u = start_wire;
            start_time_u = start_time;
            end_wire_u = end_wire;
            end_time_u = end_time;
            break;
        case 1:
            start_wire_v = start_wire;
            start_time_v = start_time;
            end_wire_v = end_wire;
            end_time_v = end_time;
            break;
        case 2:
            start_wire_y = start_wire;
            start_time_y = start_time;
            end_wire_y = end_wire;
            end_time_y = end_time;
            break;
            
        default:
            break;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Int_t PandoraNuTrack::GetStartWire (int plane) const{
    switch (plane) {
        case 0:
            return start_wire_u;
            break;
        case 1:
            return start_wire_v;
            break;
        case 2:
        default:
            return start_wire_y;
            break;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Int_t PandoraNuTrack::GetStartTime (int plane) const{
    switch (plane) {
        case 0:
            return start_time_u;
            break;
        case 1:
            return start_time_v;
            break;
        case 2:
        default:
            return start_time_y;
            break;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Int_t PandoraNuTrack::GetEndWire (int plane) const{
    switch (plane) {
        case 0:
            return end_wire_u;
            break;
        case 1:
            return end_wire_v;
            break;
        case 2:
        default:
            return end_wire_y;
            break;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Int_t PandoraNuTrack::GetEndTime (int plane) const{
    switch (plane) {
        case 0:
            return end_time_u;
            break;
        case 1:
            return end_time_v;
            break;
        case 2:
        default:
            return end_time_y;
            break;
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::DistanceFromPoint( TVector3 position, std::string * fStartOrEnd  ){
    Float_t DistanceStart, DistanceEnd , distance = 1000;
    std::string StartOrEnd = "None";
    
    DistanceStart = ( start_pos - position).Mag();
    DistanceEnd = ( end_pos - position).Mag();
    if ( DistanceStart < DistanceEnd ){
        StartOrEnd = "Start";
        distance = DistanceStart;
    }
    else{
        StartOrEnd = "End";
        distance = DistanceEnd;
    }
    
    return distance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::ClosestDistanceToOtherTrack( PandoraNuTrack other_track, std::string * fStartOrEnd ){
    Float_t MinDistanceToOtherTrack = 10000;
    std::string StartOrEnd = "None";
    Float_t DistanceStartStart = (start_pos - other_track.start_pos).Mag();
    if (MinDistanceToOtherTrack>DistanceStartStart)     {MinDistanceToOtherTrack = DistanceStartStart; StartOrEnd = "Start";}
    
    Float_t DistanceStartEnd = (start_pos - other_track.end_pos).Mag();
    if (MinDistanceToOtherTrack>DistanceStartEnd)       {MinDistanceToOtherTrack = DistanceStartEnd; StartOrEnd = "Start";}
    
    Float_t DistanceEndStart = (end_pos - other_track.start_pos).Mag();
    if (MinDistanceToOtherTrack>DistanceEndStart)       {MinDistanceToOtherTrack = DistanceEndStart; StartOrEnd = "End";}
    
    Float_t DistanceEndEnd = (end_pos - other_track.end_pos).Mag();
    if (MinDistanceToOtherTrack>DistanceEndEnd)         {MinDistanceToOtherTrack = DistanceEndEnd; StartOrEnd = "End";}
    
    
    if (fStartOrEnd!=nullptr) *fStartOrEnd = StartOrEnd;
    
    return MinDistanceToOtherTrack;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PandoraNuTrack::FlipTrack(){
    
    // flip start and end positions
    TVector3 tmp_pos = start_pos;
    start_pos   = end_pos;
    end_pos     = tmp_pos;
    
    // change angles
    rec_dir = -rec_dir;
    pandora_theta   = PI - pandora_theta;
    pandora_phi     = (pandora_phi > 0 ? pandora_phi - PI : pandora_phi + PI );
    
    IsFlipped  = !IsFlipped;
    momentum_MSCMu = momentum_MSCMu_bwd;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::GetDis2Flash (flash Flash) const {
    float TrackZcenter = (start_pos.z() + end_pos.z())/2.;
    float FlashZcenter = Flash.GetZcenter();
    float Zdis = TrackZcenter - FlashZcenter;
    float TrackYcenter = (start_pos.y() + end_pos.y())/2.;
    float FlashYcenter = Flash.GetYcenter();
    float Ydis = TrackYcenter - FlashYcenter;
    float YZdistance = sqrt( Zdis*Zdis + Ydis*Ydis );
    return YZdistance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::GetDis2ClosestFlash () const {
    return GetDis2Flash(ClosestFlash);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::ClosestLine2LineDistance ( PandoraNuTrack other_track
                                                  , TVector3 * QuasiIntersectionPoint){

    // @brief find closest line-line distance, and the intersection
    // based on [http://mathworld.wolfram.com/Line-LineDistance.html]
    // for the 3D objects:
    // track-1: start -> end, track-2: start -> end
    // define
    // a = track-1(end) - track-1(start)
    // b = track-2(end) - track-2(start)
    // c = track-2(end) - track-1(start)
    // and then
    // D = |c * (a X b)| / |a X b|
    
    TVector3 a = this->GetEndPos() - this->GetStartPos();
    TVector3 b = other_track.GetEndPos() - other_track.GetStartPos();
    TVector3 c = other_track.GetEndPos() - this->GetStartPos();
    
    return ( fabs( c.Dot( a.Cross(b) ))
            / (a.Cross(b)).Mag() );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::ClosestLine2PointDistance ( TVector3 x0 ){
    
    // @brief find closest line-point distance
    // based on [http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html]
    // for the 3D objects:
    // track-1: start -> end, point: x0
    // the closest distance is
    // D = |(x0-start) X (x0-end)| / |start-end|
    
    return (( (x0-this->GetStartPos()).Cross(x0-this->GetEndPos())  ).Mag()
            / ( this->GetStartPos()-this->GetEndPos()).Mag() );
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PandoraNuTrack::Print( bool DoPrintPandoraNuFeatures , bool DoPrintAllPIDa ) const{
    
    cout << "\033[31m" << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl
    << "track " << track_id << endl << "-------------------"    << "\033[0m" << endl;
    
    SHOWTVector3(start_pos);
    SHOWTVector3(end_pos);
    SHOWTVector3(rec_dir);
    
    if (DoPrintPandoraNuFeatures){
        SHOW3( PandoraNuCaliPID_PIDA[0] , PandoraNuCaliPID_PIDA[1]  , PandoraNuCaliPID_PIDA[2]  );
        PrintPhys(length,"cm");
        Printf("rec theta=%.1f, phi=%.1f deg.", r2d*rec_dir.Theta(), r2d*rec_dir.Phi() );
        Printf("pandora theta=%.1f, pandora phi=%.1f deg.", r2d*pandora_theta, r2d*pandora_phi );
        SHOW(IsFlipped);
    }
    SHOW3( momentum_MSCMu_fwd, momentum_MSCMu_bwd , momentum_MSCMu );
    if (MCpdgCode!=-9999){
        cout << "........................" << endl << "track MC information " << endl ;
        SHOW(MCpdgCode);
        SHOW(mcevent_id);
        SHOWTVector3(truth_start_pos);
        SHOWTVector3(truth_end_pos);
        SHOWTVector3(truth_dir);
        SHOWTLorentzVector(truth_momentum);
        SHOW3(truth_mother, truth_process, truth_origin);
        SHOW(truth_mcparticle);
        SHOW(truth_purity);
        Printf("truth theta=%.1f, phi=%.1f deg.", r2d*truth_dir.Theta(), r2d*truth_dir.Phi() );
        if ( r2d*fabs(truth_dir.Theta()-rec_dir.Theta()) > 90 ) {
            Printf("theta - truth-theta = %.1f(!)",r2d*fabs(truth_dir.Theta()-rec_dir.Theta()));
        }
        cout << "........................" << endl;
    }
    // cout << "closest-flash:" << endl;
    // ClosestFlash.Print();
    cout << "\033[31m" << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv" << "\033[0m" << endl;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Float_t PandoraNuTrack::Dedx(float dqdx){/// Box
  double Rho = 1.383;
  double betap = 0.183592;//from Minuit fit
  double alpha = 0.921969;///from Minuit fit
  double Wion = 23.6e-6;
  double Efield = 0.273;
  double ConversionFactor = 0.00411911;

  return (exp((dqdx/ConversionFactor)*(betap/(Rho*Efield))*Wion)-alpha)/(betap/(Rho*Efield));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//void PandoraNuTrack::SetPandoraNuNewCaliPID(int plane, TFile *dEdxrr) 
void PandoraNuTrack::SetPandoraNuNewCaliPID(int plane, TProfile *dedx_range_pro, TProfile *dedx_range_ka, TProfile *dedx_range_pi, TProfile *dedx_range_mu) 
{

  //TProfile* dedx_range_pro = (TProfile*) dEdxrr->Get("dedx_range_pro");
  //TProfile* dedx_range_ka = (TProfile*) dEdxrr->Get("dedx_range_ka");
  //TProfile* dedx_range_pi = (TProfile*) dEdxrr->Get("dedx_range_pi");
  //TProfile* dedx_range_mu = (TProfile*) dEdxrr->Get("dedx_range_mu");

  int npt = 0;
  double chi2pro = 0;
  double chi2ka = 0;
  double chi2pi = 0;
  double chi2mu = 0;

  double avgdedx = 0;
  double PIDA = 0;
  std::vector<double> vpida;
  std::vector<double> trkdedx;
  std::vector<double> trkres;

  vpida.clear();
  trkdedx.clear();
  trkres.clear();
  for(size_t itrk=0; itrk<dQdx[plane].size();itrk++){
    trkdedx.push_back(Dedx(dQdx[plane][itrk]));
    if(Dedx(dQdx[plane][itrk])<0){
      std::cout<<".........................................."<<std::endl;
      std::cout<<"dqdx= "<<dQdx[plane][itrk]<<" dEdx= "<<Dedx(dQdx[plane][itrk])<<std::endl;
      std::cout<<".........................................."<<std::endl;
    }
 }
 for(size_t itrk=0; itrk<ResRange[plane].size();itrk++){
    trkres.push_back(ResRange[plane][itrk]);
 }
 int used_trkres = 0;
 for (unsigned int i = 0; i<trkdedx.size(); i++){//hits
    if (i==0 || i==trkdedx.size()-1) continue;
    avgdedx += trkdedx[i];
    if(trkres[i] < 30) {
      PIDA += trkdedx[i]*pow(trkres[i],0.42);
      vpida.push_back(trkdedx[i]*pow(trkres[i],0.42));
      used_trkres++;
    }
    if (trkdedx[i]>1000) continue; //protect against large pulse height
    int bin = dedx_range_pro->FindBin(trkres[i]);
    if (bin>=1&&bin<=dedx_range_pro->GetNbinsX()){
      double bincpro = dedx_range_pro->GetBinContent(bin);
      if (bincpro<1e-6){//for 0 bin content, using neighboring bins
        bincpro = (dedx_range_pro->GetBinContent(bin-1)+dedx_range_pro->GetBinContent(bin+1))/2.;
      }
      double bincka = dedx_range_ka->GetBinContent(bin);
      if (bincka<1e-6){
        bincka = (dedx_range_ka->GetBinContent(bin-1)+dedx_range_ka->GetBinContent(bin+1))/2.;
      }
      double bincpi = dedx_range_pi->GetBinContent(bin);
      if (bincpi<1e-6){
        bincpi = (dedx_range_pi->GetBinContent(bin-1)+dedx_range_pi->GetBinContent(bin+1))/2.;
      }
      double bincmu = dedx_range_mu->GetBinContent(bin);
      if (bincmu<1e-6){
        bincmu = (dedx_range_mu->GetBinContent(bin-1)+dedx_range_mu->GetBinContent(bin+1))/2.;
      }
      double binepro = dedx_range_pro->GetBinError(bin);
      if (binepro<1e-6){
        binepro = (dedx_range_pro->GetBinError(bin-1)+dedx_range_pro->GetBinError(bin+1))/2.;
      }
      double bineka = dedx_range_ka->GetBinError(bin);
      if (bineka<1e-6){
        bineka = (dedx_range_ka->GetBinError(bin-1)+dedx_range_ka->GetBinError(bin+1))/2.;
      }
      double binepi = dedx_range_pi->GetBinError(bin);
      if (binepi<1e-6){
        binepi = (dedx_range_pi->GetBinError(bin-1)+dedx_range_pi->GetBinError(bin+1))/2.;
      }
      double binemu = dedx_range_mu->GetBinError(bin);
      if (binemu<1e-6){
        binemu = (dedx_range_mu->GetBinError(bin-1)+dedx_range_mu->GetBinError(bin+1))/2.;
      }
      double errdedx = 0.04231+0.0001783*trkdedx[i]*trkdedx[i]; //resolution on dE/dx
      errdedx *= trkdedx[i];
      chi2pro += pow((trkdedx[i]-bincpro)/std::sqrt(pow(binepro,2.)+pow(errdedx,2.)),2.);
      chi2ka += pow((trkdedx[i]-bincka)/std::sqrt(pow(bineka,2.)+pow(errdedx,2.)),2.);
      chi2pi += pow((trkdedx[i]-bincpi)/std::sqrt(pow(binepi,2.)+pow(errdedx,2.)),2.);
      chi2mu += pow((trkdedx[i]-bincmu)/std::sqrt(pow(binemu,2.)+pow(errdedx,2.)),2.);
      //if (plane ==2) std::cout<<"chi2pro "<<chi2pro<<" trkdedx[i] "<<trkdedx[i]<<" bincpro "<<bincpro<<" binepro "<<binepro<<" errdedx "<<errdedx<<std::endl;
      ++npt;
    }
  }//hits
  //std::cout<<"chi2pro "<<chi2pro<<" npt "<<npt<<" chi2 normalized "<<chi2pro/npt<<" trkdedx.size(): "<<trkdedx.size()<<" trkres.size(): "<<trkres.size()<<std::endl;
  //if (plane ==2) std::cout<<"trkdedx.size() "<<trkdedx.size()<<" ResRange[plane].size() "<<ResRange[plane].size()<<" npt "<<npt<<" chi2pro "<<chi2pro<<std::endl;
  if (npt>0) PandoraNuCaliPID_NewChi2Proton[plane] = chi2pro/npt;
  if (npt>0) PandoraNuCaliPID_NewChi2Muon[plane] = chi2mu/npt;
  else PandoraNuCaliPID_NewChi2Proton[plane] = -9998;
} 

#endif

























