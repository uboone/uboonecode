#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"

#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>

using namespace std;

namespace microboone{

class Datamuonvalidation : public art::EDAnalyzer {
public:

    explicit Datamuonvalidation(fhicl::ParameterSet const& pset);
    virtual ~Datamuonvalidation();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void reset();
    
private:
    TTree* fEventTree;
    Int_t    run;                  
    Int_t    subrun;               
    Int_t    event;
    Int_t    ismu_p0;
    Int_t    ismu_p1;
    Int_t    ismu_p2;
    Float_t  mu_trklen_p0;
    Float_t  mu_trklen_p1;
    Float_t  mu_trklen_p2;
    Float_t  mu_thetaxz_p0;
    Float_t  mu_thetaxz_p1;
    Float_t  mu_thetaxz_p2;
    Float_t  mu_thetayz_p0;
    Float_t  mu_thetayz_p1;
    Float_t  mu_thetayz_p2;
    Int_t    mu_TrkID_p0;
    Int_t    mu_TrkID_p1;
    Int_t    mu_TrkID_p2;
    Float_t  mu_start_info_p0[3];
    Float_t  mu_start_info_p1[3];
    Float_t  mu_start_info_p2[3];
    Float_t  mu_end_info_p0[3];
    Float_t  mu_end_info_p1[3];
    Float_t  mu_end_info_p2[3];
    Int_t    mu_hits_p0;
    Int_t    mu_hits_p1;
    Int_t    mu_hits_p2;
    Float_t  mu_trkdqdx_p0[5000];
    Float_t  mu_trkdqdx_p1[5000];
    Float_t  mu_trkdqdx_p2[5000];
    Float_t  mu_resrange_p0[5000];
    Float_t  mu_resrange_p1[5000];
    Float_t  mu_resrange_p2[5000];
    Float_t  mu_trkdedx_p0[5000];
    Float_t  mu_trkdedx_p1[5000];
    Float_t  mu_trkdedx_p2[5000];
    
    
    TH2F *dqdx_res_range_p0_hist;
    TH2F *dqdx_res_range_p1_hist;
    TH2F *dqdx_res_range_p2_hist;
    TH2F *dedx_res_range_p0_hist;
    TH2F *dedx_res_range_p1_hist;
    TH2F *dedx_res_range_p2_hist;
    
    std::string fDigitModuleLabel;
    std::string fHitsModuleLabel;
    std::string fLArG4ModuleLabel;
    std::string fGenieGenModuleLabel;
    std::string fTrackModuleLabel;
    std::string fPOTModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fParticleIDModuleLabel;
    std::string fClusterModuleLabel;
    std::string fVertexModuleLabel;
    bool  fSaveCaloInfo;
    bool  fSaveTrackInfo;
    bool  fSaveClusterInfo;
    bool  fSaveClusterHitInfo;
    bool  fSaveGenieInfo;
    bool  fSaveVertexInfo; 
    float fG4minE;
}; 

//========================================================================
Datamuonvalidation::Datamuonvalidation(fhicl::ParameterSet const& pset) :
EDAnalyzer(pset),
  fDigitModuleLabel         (pset.get< std::string >("DigitModuleLabel","")        ),
  fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","")         ),
  fLArG4ModuleLabel         (pset.get< std::string >("LArGeantModuleLabel","")     ),
  fGenieGenModuleLabel      (pset.get< std::string >("GenieGenModuleLabel","")     ),  
  fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
  fPOTModuleLabel           (pset.get< std::string >("POTModuleLabel","")          ),
  fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
  fParticleIDModuleLabel    (pset.get< std::string >("ParticleIDModuleLabel","")   ),
  fClusterModuleLabel       (pset.get< std::string >("ClusterModuleLabel","")      ),
  fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel","")       ),
  fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false)),
  fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false)),  
  fSaveClusterInfo          (pset.get< bool>("SaveClusterInfo",false)),
  fSaveClusterHitInfo       (pset.get< bool>("SaveClusterHitInfo",false)),
  fSaveGenieInfo            (pset.get< bool>("SaveGenieInfo",false)),
  fSaveVertexInfo           (pset.get< bool>("SaveVertexInfo",false)),
  fG4minE                   (pset.get< float>("G4minE",0.01))  
{
  if (fSaveTrackInfo == false) fSaveCaloInfo = false;
}
 
//========================================================================
Datamuonvalidation::~Datamuonvalidation(){
}
//========================================================================

//========================================================================
void Datamuonvalidation::beginJob(){
  std::cout<<"job begin..."<<std::endl;
  art::ServiceHandle<art::TFileService> tfs;
  fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
  fEventTree->Branch("event", &event,"event/I");
  fEventTree->Branch("run", &run,"run/I");
  fEventTree->Branch("subrun", &subrun,"surbrun/I");
  fEventTree->Branch("mu_trklen_p0",&mu_trklen_p0,"mu_trklen_p0/I");
  fEventTree->Branch("mu_trklen_p1",&mu_trklen_p1,"mu_trklen_p1/I");
  fEventTree->Branch("mu_trklen_p2",&mu_trklen_p2,"mu_trklen_p2/I");
  fEventTree->Branch("mu_thetaxz_p0",&mu_thetaxz_p0,"mu_thetaxz_p0/I");
  fEventTree->Branch("mu_thetaxz_p1",&mu_thetaxz_p1,"mu_thetaxz_p1/I");
  fEventTree->Branch("mu_thetaxz_p2",&mu_thetaxz_p2,"mu_thetaxz_p2/I");
  fEventTree->Branch("mu_thetayz_p0",&mu_thetayz_p0,"mu_thetayz_p0/I");
  fEventTree->Branch("mu_thetayz_p1",&mu_thetayz_p1,"mu_thetayz_p1/I");
  fEventTree->Branch("mu_thetayz_p2",&mu_thetayz_p2,"mu_thetayz_p2/I");
  fEventTree->Branch("mu_TrkID_p0",&mu_TrkID_p0,"mu_TrkID_p0/I");
  fEventTree->Branch("mu_TrkID_p1",&mu_TrkID_p1,"mu_TrkID_p1/I");
  fEventTree->Branch("mu_TrkID_p2",&mu_TrkID_p2,"mu_TrkID_p2/I");
  fEventTree->Branch("mu_start_info_p0",mu_start_info_p0,"mu_start_info_p0[3]/F");
  fEventTree->Branch("mu_start_info_p1",mu_start_info_p1,"mu_start_info_p1[3]/F");
  fEventTree->Branch("mu_start_info_p2",mu_start_info_p2,"mu_start_info_p2[3]/F");
  fEventTree->Branch("mu_end_info_p0",mu_end_info_p0,"mu_end_info_p0[3]/F");
  fEventTree->Branch("mu_end_info_p1",mu_end_info_p1,"mu_end_info_p1[3]/F");
  fEventTree->Branch("mu_end_info_p2",mu_end_info_p2,"mu_end_info_p2[3]/F");
  fEventTree->Branch("mu_hits_p0",&mu_hits_p0,"mu_hits_p0/I");
  fEventTree->Branch("mu_hits_p1",&mu_hits_p1,"mu_hits_p1/I");
  fEventTree->Branch("mu_hits_p2",&mu_hits_p2,"mu_hits_p2/I");
  fEventTree->Branch("ismu_p0",&ismu_p0,"ismu_p0/I");
  fEventTree->Branch("ismu_p1",&ismu_p1,"ismu_p1/I");
  fEventTree->Branch("ismu_p2",&ismu_p2,"ismu_p2/I");
  fEventTree->Branch("mu_trkdqdx_p0",mu_trkdqdx_p0,"mu_trkdqdx_p0[5000]/F");
  fEventTree->Branch("mu_trkdqdx_p1",mu_trkdqdx_p1,"mu_trkdqdx_p1[5000]/F");
  fEventTree->Branch("mu_trkdqdx_p2",mu_trkdqdx_p2,"mu_trkdqdx_p2[5000]/F");
  fEventTree->Branch("mu_resrange_p0",mu_resrange_p0,"mu_resrange_p0[5000]/F");
  fEventTree->Branch("mu_resrange_p1",mu_resrange_p1,"mu_resrange_p1[5000]/F");
  fEventTree->Branch("mu_resrange_p2",mu_resrange_p2,"mu_resrange_p2[5000]/F");
  fEventTree->Branch("mu_trkdedx_p0",mu_trkdedx_p0,"mu_trkdedx_p0[5000]/F");
  fEventTree->Branch("mu_trkdedx_p1",mu_trkdedx_p1,"mu_trkdedx_p1[5000]/F");
  fEventTree->Branch("mu_trkdedx_p2",mu_trkdedx_p2,"mu_trkdedx_p2[5000]/F");
  
  
  dqdx_res_range_p0_hist=tfs->make<TH2F>("dqdx_res_range_p0_hist","Plane 0",300,0,100,3000,0,1000);
  dqdx_res_range_p0_hist->GetXaxis()->SetTitle("Residual range (cm)");
  dqdx_res_range_p0_hist->GetYaxis()->SetTitle("dQ/dx (ADC/cm)");
  
  dqdx_res_range_p1_hist=tfs->make<TH2F>("dqdx_res_range_p1_hist","Plane 1",300,0,100,3000,0,1000);
  dqdx_res_range_p1_hist->GetXaxis()->SetTitle("Residual range (cm)");
  dqdx_res_range_p1_hist->GetYaxis()->SetTitle("dQ/dx (ADC/cm)");
  
  dqdx_res_range_p2_hist=tfs->make<TH2F>("dqdx_res_range_p2_hist","Plane 2",300,0,100,3000,0,1000);
  dqdx_res_range_p2_hist->GetXaxis()->SetTitle("Residual range (cm)");
  dqdx_res_range_p2_hist->GetYaxis()->SetTitle("dQ/dx (ADC/cm)");
  
  dedx_res_range_p0_hist=tfs->make<TH2F>("dedx_res_range_p0_hist","Plane 0",300,0,100,150,0,50);
  dedx_res_range_p0_hist->GetXaxis()->SetTitle("Residual range (cm)");
  dedx_res_range_p0_hist->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
  
  dedx_res_range_p1_hist=tfs->make<TH2F>("dedx_res_range_p1_hist","Plane 1",300,0,100,150,0,50);
  dedx_res_range_p1_hist->GetXaxis()->SetTitle("Residual range (cm)");
  dedx_res_range_p1_hist->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
  
  dedx_res_range_p2_hist=tfs->make<TH2F>("dedx_res_range_p2_hist","Plane 2",300,0,100,150,0,50);
  dedx_res_range_p2_hist->GetXaxis()->SetTitle("Residual range (cm)");
  dedx_res_range_p2_hist->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
  
}

//========================================================================
void Datamuonvalidation::endJob(){     

}

//========================================================================
void Datamuonvalidation::beginRun(const art::Run&){
  mf::LogInfo("Datamuonvalidation")<<"begin run..."<<std::endl;
}
//========================================================================

//========================================================================

//========================================================================

void Datamuonvalidation::analyze( const art::Event& evt){
     reset();  
     
     art::Handle< std::vector<recob::Track> > trackListHandle;
     std::vector<art::Ptr<recob::Track> > tracklist;
  
     if(evt.getByLabel(fTrackModuleLabel,trackListHandle)) art::fill_ptr_vector(tracklist, trackListHandle);
  
     art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
     auto const &assoc_handle=evt.getValidHandle<art::Assns<recob::Vertex, recob::Track>>("NuMuCCSelectionII");
     
     run = evt.run();
     subrun = evt.subRun();
     event = evt.id().event();
     
     int cc_mu_index=-1;
     
     std::vector<recob::Vertex> Sel2_vtx_list;
     std::vector<recob::Track> Sel2_trk_list;
     
     for(auto &ass : *assoc_handle){
         art::Ptr<recob::Vertex> Sel2_vtx=ass.first;
         Sel2_vtx_list.emplace_back(*Sel2_vtx);

         art::Ptr<recob::Track> Sel2_trk = ass.second;
         Sel2_trk_list.emplace_back(*Sel2_trk);
      
         cc_mu_index=Sel2_trk.key();
    }
     
    if(cc_mu_index!=-1){
        float st_x;
        float st_y;
        float st_z;
        float en_x;
        float en_y;
        float en_z;
	art::Ptr<recob::Track> ptrack(trackListHandle,cc_mu_index);
        const recob::Track& track = *ptrack;
	TVector3 pos,end,dir_start,dir_end;
	pos = track.Vertex();
	end = track.End();
	dir_start=track.VertexDirection();
     	dir_end=track.EndDirection();
	st_x=pos.X();st_y=pos.Y();st_z=pos.Z();
	en_x=end.X();en_y=end.Y();en_z=end.Z();
	double theta_xz=std::atan2(dir_start.X(), dir_start.Z());
        double theta_yz=std::atan2(dir_start.Y(), dir_start.Z());
	float mu_tlen=track.Length();
	if((st_x>10 && st_x<250) && (en_x>10 && en_x<250) && (st_y>-110 && st_y<110) && (en_y>-110 && en_y<110) && (st_z>10 && st_z<1030) && (en_z>10 && en_z<1030)){
	    if(mu_tlen>150){
	       if(!(TMath::Abs(theta_xz)>1.31 && TMath::Abs(theta_xz)<1.83)){
	            std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(cc_mu_index);
	            for(size_t ical = 0; ical<calos.size(); ++ical){
	                if(!calos[ical]) continue;
	                if(!calos[ical]->PlaneID().isValid) continue;
	                int planenum = calos[ical]->PlaneID().Plane;
	                if(planenum<0||planenum>2) continue;
	                const size_t NHits = calos[ical] -> dEdx().size();
	    
	                if(planenum==0){
	                   if((!(theta_yz>-2.967 && theta_yz<-2.269)) && (!(theta_yz>0.1745 && theta_yz<0.873))){
			         std::vector<float> dqdx_large_res,dqdx_small_res;
				 for(size_t iHit = 0; iHit < NHits; ++iHit){
	                             if((calos[ical]->ResidualRange())[iHit]<=5) dqdx_small_res.push_back((calos[ical] -> dQdx())[iHit]);
				     if((calos[ical]->ResidualRange())[iHit]<=150 && (calos[ical]->ResidualRange())[iHit]>145) dqdx_large_res.push_back((calos[ical] -> dQdx())[iHit]);
			         }
			         float ratio=0;
			         float denominator=0;
			         float numerator=0;
			         if(dqdx_small_res.size()) numerator=TMath::Median(dqdx_small_res.size(),&dqdx_small_res[0]);
			         if(dqdx_large_res.size()) denominator=TMath::Median(dqdx_large_res.size(),&dqdx_large_res[0]);
			         if(numerator!=0 && denominator!=0) ratio=float(numerator)/denominator;
				 if(ratio>=1.5){
				    mu_trklen_p0=mu_tlen;
				    mu_thetaxz_p0=theta_xz;
				    mu_thetayz_p0=theta_yz;
				    mu_TrkID_p0=track.ID();
				    mu_start_info_p0[0]=st_x;mu_start_info_p0[1]=st_y;mu_start_info_p0[2]=st_z;
				    mu_end_info_p0[0]=en_x;mu_end_info_p0[1]=en_y;mu_end_info_p0[2]=en_z;
				    mu_hits_p0=int(NHits);
				    ismu_p0=1;
				    for(size_t iHit = 0; iHit < NHits; ++iHit){
				        if((calos[ical] -> TrkPitchVec())[iHit]>=0.3 && (calos[ical] -> TrkPitchVec())[iHit]<=0.4){
					    dqdx_res_range_p0_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical] -> dQdx())[iHit]);
		                            dedx_res_range_p0_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical] -> dEdx())[iHit]);
					    mu_trkdqdx_p0[iHit]=(calos[ical] -> dQdx())[iHit];
					    mu_resrange_p0[iHit]=(calos[ical]->ResidualRange())[iHit];
					    mu_trkdedx_p0[iHit]=(calos[ical] -> dEdx())[iHit];
					}
				    }
				}
			    }
		        }
	    
	                if(planenum==1){
	                   if((!(theta_yz<2.967 && theta_yz>2.269)) && (!(theta_yz<-0.1745 && theta_yz>-0.873))){
			         std::vector<float> dqdx_large_res,dqdx_small_res;
				 for(size_t iHit = 0; iHit < NHits; ++iHit){
				     if((calos[ical]->ResidualRange())[iHit]<=5) dqdx_small_res.push_back((calos[ical] -> dQdx())[iHit]);
				     if((calos[ical]->ResidualRange())[iHit]<=150 && (calos[ical]->ResidualRange())[iHit]>145) dqdx_large_res.push_back((calos[ical] -> dQdx())[iHit]);
			         }
				 float ratio=0;
			         float denominator=0;
			         float numerator=0;
			         if(dqdx_small_res.size()) numerator=TMath::Median(dqdx_small_res.size(),&dqdx_small_res[0]);
			         if(dqdx_large_res.size()) denominator=TMath::Median(dqdx_large_res.size(),&dqdx_large_res[0]);
			         if(numerator!=0 && denominator!=0) ratio=float(numerator)/denominator;
				 if(ratio>=1.5){
				    mu_trklen_p1=mu_tlen;
				    mu_thetaxz_p1=theta_xz;
				    mu_thetayz_p1=theta_yz;
				    mu_TrkID_p1=track.ID();
				    mu_start_info_p1[0]=st_x;mu_start_info_p1[1]=st_y;mu_start_info_p1[2]=st_z;
				    mu_end_info_p1[0]=en_x;mu_end_info_p1[1]=en_y;mu_end_info_p1[2]=en_z;
				    mu_hits_p1=int(NHits);
				    ismu_p1=1;
				    for(size_t iHit = 0; iHit < NHits; ++iHit){
				        if((calos[ical] -> TrkPitchVec())[iHit]>=0.3 && (calos[ical] -> TrkPitchVec())[iHit]<=0.4){
					    dqdx_res_range_p1_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical] -> dQdx())[iHit]);
		                            dedx_res_range_p1_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical] -> dEdx())[iHit]);
					    mu_trkdqdx_p1[iHit]=(calos[ical] -> dQdx())[iHit];
					    mu_resrange_p1[iHit]=(calos[ical]->ResidualRange())[iHit];
					    mu_trkdedx_p1[iHit]=(calos[ical] -> dEdx())[iHit];
					}
				    }
				}   
			    } 
	                }
	    
	                if(planenum==2){
	                   if(!(TMath::Abs(theta_yz)>1.396 && TMath::Abs(theta_yz)<1.745)){
			        std::vector<float> dqdx_large_res,dqdx_small_res;
			        for(size_t iHit = 0; iHit < NHits; ++iHit){
	                            if((calos[ical]->ResidualRange())[iHit]<=5) dqdx_small_res.push_back((calos[ical] -> dQdx())[iHit]);
				    if((calos[ical]->ResidualRange())[iHit]<=150 && (calos[ical]->ResidualRange())[iHit]>145) dqdx_large_res.push_back((calos[ical] -> dQdx())[iHit]);
			        }
				float ratio=0;
			        float denominator=0;
			        float numerator=0;
			        if(dqdx_small_res.size()) numerator=TMath::Median(dqdx_small_res.size(),&dqdx_small_res[0]);
			        if(dqdx_large_res.size()) denominator=TMath::Median(dqdx_large_res.size(),&dqdx_large_res[0]);
			        if(numerator!=0 && denominator!=0) ratio=float(numerator)/denominator;
				if(ratio>=1.5){
				   mu_trklen_p2=mu_tlen;
				   mu_thetaxz_p2=theta_xz;
				   mu_thetayz_p2=theta_yz;
				   mu_TrkID_p2=track.ID();
				   mu_start_info_p2[0]=st_x;mu_start_info_p2[1]=st_y;mu_start_info_p2[2]=st_z;
				   mu_end_info_p2[0]=en_x;mu_end_info_p2[1]=en_y;mu_end_info_p2[2]=en_z;
				   mu_hits_p2=int(NHits);
				   ismu_p2=1;
				   for(size_t iHit = 0; iHit < NHits; ++iHit){
				       if((calos[ical] -> TrkPitchVec())[iHit]>=0.3 && (calos[ical] -> TrkPitchVec())[iHit]<=0.4){
					    dqdx_res_range_p2_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical] -> dQdx())[iHit]);
		                            dedx_res_range_p2_hist->Fill((calos[ical]->ResidualRange())[iHit],(calos[ical] -> dEdx())[iHit]);
					    mu_trkdqdx_p2[iHit]=(calos[ical] -> dQdx())[iHit];
					    mu_resrange_p2[iHit]=(calos[ical]->ResidualRange())[iHit];
					    mu_trkdedx_p2[iHit]=(calos[ical] -> dEdx())[iHit];
					}
				    }
				}   
			    }
	               }
	    
	           } // loop over cal objects
                } // theta XZ cut
	     } // tracklen cut
         } // muon is contained
    } // cc_mu_index !=-1
     
    fEventTree->Fill();
} // end of analyze function
	   
 /////////////////// Defintion of reset function ///////////
void Datamuonvalidation::reset(){
     run = -9999;
     subrun = -9999;
     event = -9999;
     mu_trklen_p0=-9999;
     mu_trklen_p1=-9999;
     mu_trklen_p2=-9999;
     mu_thetaxz_p0=-9999;
     mu_thetaxz_p1=-9999;
     mu_thetaxz_p2=-9999;
     mu_thetayz_p0=-9999;
     mu_thetayz_p1=-9999;
     mu_thetayz_p2=-9999;
     mu_TrkID_p0=-9999;
     mu_TrkID_p1=-9999;
     mu_TrkID_p2=-9999;
     mu_hits_p0=-9999;
     mu_hits_p1=-9999;
     mu_hits_p2=-9999;
     ismu_p0=-9999;
     ismu_p1=-9999;
     ismu_p2=-9999;
     for(int i=0; i<3; i++){
         mu_start_info_p0[i]=-9999;
	 mu_start_info_p1[i]=-9999;
	 mu_start_info_p2[i]=-9999;
	 mu_end_info_p0[i]=-9999;
	 mu_end_info_p1[i]=-9999;
	 mu_end_info_p2[i]=-9999;
     }
     
     for(int i=0; i<5000; i++){
         mu_trkdqdx_p0[i]=-9999;
	 mu_trkdqdx_p1[i]=-9999;
	 mu_trkdqdx_p2[i]=-9999;
	 mu_resrange_p0[i]=-9999;
	 mu_resrange_p1[i]=-9999;
	 mu_resrange_p2[i]=-9999;
	 mu_trkdedx_p0[i]=-9999;
	 mu_trkdedx_p1[i]=-9999;
	 mu_trkdedx_p2[i]=-9999;
     }
     
}
 //////////////////////// End of definition ///////////////	
	  
DEFINE_ART_MODULE(Datamuonvalidation)
}


