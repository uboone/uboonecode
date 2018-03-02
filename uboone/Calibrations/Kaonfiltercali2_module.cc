#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
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
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
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

//////////////////////////// Local function definitions ///////////////////////////

float my_minimum_distance(vector<float>& par_1,vector<float>& par_2){
      float par1_st_par2_st=TMath::Sqrt((par_1[0]-par_2[0])*(par_1[0]-par_2[0]) + (par_1[1]-par_2[1])*(par_1[1]-par_2[1]) + (par_1[2]-par_2[2])*(par_1[2]-par_2[2]));
      float par1_st_par2_en=TMath::Sqrt((par_1[0]-par_2[3])*(par_1[0]-par_2[3]) + (par_1[1]-par_2[4])*(par_1[1]-par_2[4]) + (par_1[2]-par_2[5])*(par_1[2]-par_2[5]));
      float par1_en_par2_st=TMath::Sqrt((par_1[3]-par_2[0])*(par_1[3]-par_2[0]) + (par_1[4]-par_2[1])*(par_1[4]-par_2[1]) + (par_1[5]-par_2[2])*(par_1[5]-par_2[2]));
      float par1_en_par2_en=TMath::Sqrt((par_1[3]-par_2[3])*(par_1[3]-par_2[3]) + (par_1[4]-par_2[4])*(par_1[4]-par_2[4]) + (par_1[5]-par_2[5])*(par_1[5]-par_2[5]));
      std::vector<float> dis_vector;
      dis_vector.push_back(par1_st_par2_st);
      dis_vector.push_back(par1_st_par2_en);
      dis_vector.push_back(par1_en_par2_st);
      dis_vector.push_back(par1_en_par2_en);
      float min_dis=5e20;
      for(unsigned int i=0; i<dis_vector.size(); i++){
          if(dis_vector[i]<min_dis) min_dis=dis_vector[i];
      }
      return min_dis;
}

int my_close_distance_index(vector<float>& par_1,vector<float>& par_2){
    float par1_st_par2_st=TMath::Sqrt((par_1[0]-par_2[0])*(par_1[0]-par_2[0]) + (par_1[1]-par_2[1])*(par_1[1]-par_2[1]) + (par_1[2]-par_2[2])*(par_1[2]-par_2[2]));
    float par1_st_par2_en=TMath::Sqrt((par_1[0]-par_2[3])*(par_1[0]-par_2[3]) + (par_1[1]-par_2[4])*(par_1[1]-par_2[4]) + (par_1[2]-par_2[5])*(par_1[2]-par_2[5]));
    float par1_en_par2_st=TMath::Sqrt((par_1[3]-par_2[0])*(par_1[3]-par_2[0]) + (par_1[4]-par_2[1])*(par_1[4]-par_2[1]) + (par_1[5]-par_2[2])*(par_1[5]-par_2[2]));
    float par1_en_par2_en=TMath::Sqrt((par_1[3]-par_2[3])*(par_1[3]-par_2[3]) + (par_1[4]-par_2[4])*(par_1[4]-par_2[4]) + (par_1[5]-par_2[5])*(par_1[5]-par_2[5]));
    std::vector<float> dis_vector;
    dis_vector.push_back(par1_st_par2_st);
    dis_vector.push_back(par1_st_par2_en);
    dis_vector.push_back(par1_en_par2_st);
    dis_vector.push_back(par1_en_par2_en);
    float min_dis=5e20;
    int indicator=-1;
    for(unsigned int i=0; i<dis_vector.size(); i++){
        if(dis_vector[i]<min_dis){ 
	    min_dis=dis_vector[i];
	    indicator=i;
	}
    }
    return indicator;
}

////////////////////////// End of local function definitions /////////////////////

class Kaonfiltercali2 : public art::EDAnalyzer {
public:

    explicit Kaonfiltercali2(fhicl::ParameterSet const& pset);
    virtual ~Kaonfiltercali2();

    void beginJob();
    void endSubRun(const art::SubRun& sr);
    void analyze(const art::Event& evt);
    void reset();
    
private:
    TTree* fPOT;
    Double_t potbnb;
    TTree* fEventTree;
    Int_t  run;                  
    Int_t  subrun;               
    Int_t  event;
    Int_t  rec_k_mu;
    Int_t  rec_pri_mu;
    Int_t  n_cc_reco;
    Int_t  n_k_true_rec_match;
    Int_t  n_primu_true_rec_match;
    Int_t  n_decmu_true_rec_match;
    Int_t  n_true_events;
    
    TH1F *flash_z_distance_hist;
    TH1F *large_flash_hist;
    
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
    std::string fOpFlashModuleLabel;
    bool  fSaveCaloInfo;
    bool  fSaveTrackInfo;
    bool  fSaveClusterInfo;
    bool  fSaveClusterHitInfo;
    bool  fSaveGenieInfo;
    bool  fSaveVertexInfo; 
    bool  fSaveFlashInfo;
   float fG4minE;
}; 

//========================================================================
Kaonfiltercali2::Kaonfiltercali2(fhicl::ParameterSet const& pset) :
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
  fOpFlashModuleLabel       (pset.get< std::string >("OpFlashModuleLabel","")      ),
  fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false)),
  fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false)),  
  fSaveClusterInfo          (pset.get< bool>("SaveClusterInfo",false)),
  fSaveClusterHitInfo       (pset.get< bool>("SaveClusterHitInfo",false)),
  fSaveGenieInfo            (pset.get< bool>("SaveGenieInfo",false)),
  fSaveVertexInfo           (pset.get< bool>("SaveVertexInfo",false)),
  fSaveFlashInfo            (pset.get< bool>("SaveFlashInfo",false)),
  fG4minE                   (pset.get< float>("G4minE",0.01))  

{
  if (fSaveTrackInfo == false) fSaveCaloInfo = false;
}
 
//========================================================================
Kaonfiltercali2::~Kaonfiltercali2(){
  //destructor
}
//========================================================================

//========================================================================

void Kaonfiltercali2::beginJob(){
  std::cout<<"job begin..."<<std::endl;

  art::ServiceHandle<art::TFileService> tfs;
  
  fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
  fPOT = tfs->make<TTree>("POT", "Pot Tree from Reco");
  fPOT->Branch("potbnb",&potbnb,"potbnb/D");
  fEventTree->Branch("event", &event);
  fEventTree->Branch("run", &run);
  fEventTree->Branch("subrun", &subrun);
  fEventTree->Branch("rec_k_mu",&rec_k_mu,"rec_k_mu/I");
  fEventTree->Branch("rec_pri_mu",&rec_pri_mu,"rec_pri_mu/I");
  fEventTree->Branch("n_cc_reco",&n_cc_reco,"n_cc_reco/I");
  fEventTree->Branch("n_k_true_rec_match",&n_k_true_rec_match,"n_k_true_rec_match/I");
  fEventTree->Branch("n_primu_true_rec_match",&n_primu_true_rec_match,"n_primu_true_rec_match/I");
  fEventTree->Branch("n_decmu_true_rec_match",&n_decmu_true_rec_match,"n_decmu_true_rec_match/I");
  fEventTree->Branch("n_true_events",&n_true_events,"n_true_events/I");
  
  flash_z_distance_hist=tfs->make<TH1F>("flash_z_distance_hist","",100,0,1000);
  flash_z_distance_hist->GetXaxis()->SetTitle("Flash Z distance");
  
  large_flash_hist=tfs->make<TH1F>("large_flash_hist","Flash PE of largest flash in the event",200,0,5000);
  large_flash_hist->GetXaxis()->SetTitle("Flash PE");
  
}

void Kaonfiltercali2::analyze( const art::Event& evt){

     reset(); 
     bool isMC = !evt.isRealData();
  
     art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
     std::vector<art::Ptr<simb::MCTruth> > mclist;
     if(isMC){
        if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
        art::fill_ptr_vector(mclist, mctruthListHandle);
     } 
  
     art::Handle< std::vector<recob::Track> > trackListHandle;
     std::vector<art::Ptr<recob::Track> > tracklist;
     if(fSaveTrackInfo){
     if(evt.getByLabel(fTrackModuleLabel,trackListHandle))
        art::fill_ptr_vector(tracklist, trackListHandle);
     }
  
     art::Handle< std::vector<recob::Hit> > hitListHandle;
     std::vector<art::Ptr<recob::Hit> > hitlist;
     if(evt.getByLabel(fHitsModuleLabel,hitListHandle))
     art::fill_ptr_vector(hitlist, hitListHandle);

     art::ServiceHandle<cheat::BackTracker> bt;
  
     art::Handle< std::vector<recob::Vertex> > vertexListHandle;
     std::vector<art::Ptr<recob::Vertex> > vertexlist;
     if(fSaveVertexInfo){
         if(evt.getByLabel(fVertexModuleLabel,vertexListHandle))
            art::fill_ptr_vector(vertexlist,vertexListHandle);
     }
     
     art::Handle<std::vector<recob::OpFlash> > flashListHandle;
     std::vector<art::Ptr<recob::OpFlash> > flashlist;
     if(fSaveFlashInfo){
        if(evt.getByLabel(fOpFlashModuleLabel,flashListHandle))
           art::fill_ptr_vector(flashlist, flashListHandle);
     }	 
  
     art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
     art::FindManyP<recob::Hit> fmht(trackListHandle,evt,fTrackModuleLabel);
     art::FindManyP<anab::ParticleID> fmpid(trackListHandle, evt, fParticleIDModuleLabel);
     
     run = evt.run();
     subrun = evt.subRun();
     event = evt.id().event();
     
     rec_k_mu=0;
     rec_pri_mu=0;
     n_cc_reco=0;
     n_k_true_rec_match=0;
     n_primu_true_rec_match=0;
     n_decmu_true_rec_match=0;
     
     if(isMC){
       if(fSaveFlashInfo){
	   int NFlashes=flashlist.size();
	   int large_flash=0;
	   float large_flash_z=-1;
	   int dec_mu_index=-1;
	   int pri_k_index=-1;
	   int pri_mu_index=-1;
	   bool k_match=false;
	   bool dec_mu_match=false;
	   bool pri_mu_match=false;
	   int global_k_mu_min_dis_index=-1;
	   
	   ////////////////////////////////// loop over flashes ///////////////////////
	   
	   for(int i=0; i<NFlashes;++i){
	      float flash_time=flashlist[i]->Time();
	      if(flash_time<=4.8 && flash_time>=3.2){
	         if(flashlist[i]->TotalPE()>large_flash){ 
		    large_flash=flashlist[i]->TotalPE();
		    large_flash_z=flashlist[i]->ZCenter();
		 }
	      }
	   }
	   
	   large_flash_hist->Fill(large_flash);
	   flash_z_distance_hist->Fill(large_flash_z);
	   
	   /////////////////////////////// End of looping over flashes //////////////// 
	   
	   if(fSaveTrackInfo){
	      int NTracks=tracklist.size();
	      if(NTracks>=3){ // (NTracks>=3 && large_flash>=75 && large_flash_z!=-1)
	         
		 std::vector<int> k_id_vec,dec_mu_id_vec;
		 std::vector<float> k_pid_vec;
		 std::vector<float> k_mu_min_dis_vec;
		 std::vector<int> k_mu_min_dis_index_vec;
		 
		 for(int i=0; i<NTracks;++i){
		     art::Ptr<recob::Track> ptrack(trackListHandle, i);
		     const recob::Track& track = *ptrack;
	             TVector3 pos,end;
	             pos = track.Vertex();
		     end = track.End();
		     if((pos.X()>5 && pos.X()<251.3) && (end.X()>5 && end.X()<251.3) && (pos.Y()>-111.5 && pos.Y()<111.5) && (end.Y()>-111.5 && end.Y()<111.5) && (pos.Z()>5 && pos.Z()<1031.8) && (end.Z()>5 && end.Z()<1031.8)){
		         std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(i);
			 float k_tlen=track.Length();
			 int k_hits_p2=0;
			 int k_hits_p1=0;
			 int k_hits_p0=0;
			 
			 for(unsigned int ical=0; ical<calos.size(); ++ical){
			     if(!calos[ical]) continue;
			     if(!calos[ical]->PlaneID().isValid) continue;
			     int planenum = calos[ical]->PlaneID().Plane;
			     if(planenum<0||planenum>2) continue; 
			     if(planenum==2) k_hits_p2=calos[ical]->dEdx().size();
			     if(planenum==1) k_hits_p1=calos[ical]->dEdx().size();
			     if(planenum==0) k_hits_p0=calos[ical]->dEdx().size();
			 }
			 
			 if(k_tlen>=5 && (k_hits_p2>=10 || k_hits_p1>=10 || k_hits_p0>=10)){
			 
			    float k_pida=-1;
			    std::vector<art::Ptr<anab::ParticleID>> pids=fmpid.at(i);
			    for(unsigned int ipid=0; ipid<pids.size(); ++ipid){
			        if(!pids[ipid]->PlaneID().isValid) continue;
				int planenum=pids[ipid]->PlaneID().Plane;
				if(planenum<0||planenum>2) continue;
				if(planenum==2) k_pida=float(pids[ipid]->PIDA());
			    }
			    
			    if(k_pida>=10 && k_pida<=20){
			    
			       ////////////////////////////////////// Getting decaying muon ////////////////////////////////////
			       
			       for(int j=0; j<NTracks;++j){
			           if(j==i) continue;
				   art::Ptr<recob::Track> ptrack(trackListHandle, j);
		                   const recob::Track& track = *ptrack;
	                           TVector3 pos,end;
	                           pos = track.Vertex();
		                   end = track.End();
				   if((pos.X()>5 && pos.X()<251.3) && (end.X()>5 && end.X()<251.3) && (pos.Y()>-111.5 && pos.Y()<111.5) && (end.Y()>-111.5 && end.Y()<111.5) && (pos.Z()>5 && pos.Z()<1031.8) && (end.Z()>5 && end.Z()<1031.8)){
				       std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(j);
			               float dec_mu_tlen=track.Length();
			               int dec_mu_hits_p2=0;
			               int dec_mu_hits_p1=0;
			               int dec_mu_hits_p0=0;
				       
				       for(unsigned int ical=0; ical<calos.size(); ++ical){
			                   if(!calos[ical]) continue;
			                   if(!calos[ical]->PlaneID().isValid) continue;
			                   int planenum = calos[ical]->PlaneID().Plane;
			                   if(planenum<0||planenum>2) continue; 
			                   if(planenum==2) dec_mu_hits_p2=calos[ical]->dEdx().size();
			                   if(planenum==1) dec_mu_hits_p1=calos[ical]->dEdx().size();
			                   if(planenum==0) dec_mu_hits_p0=calos[ical]->dEdx().size();
			               }
				       
				       if(dec_mu_tlen>=45 && dec_mu_tlen<=60 && (dec_mu_hits_p2>=10 || dec_mu_hits_p1>=10 || dec_mu_hits_p0>=10)){
				       
				          float dec_mu_pida=-1;
			                  std::vector<art::Ptr<anab::ParticleID>> pids=fmpid.at(j);
			                  for(unsigned int ipid=0; ipid<pids.size(); ++ipid){
			                      if(!pids[ipid]->PlaneID().isValid) continue;
				              int planenum=pids[ipid]->PlaneID().Plane;
				              if(planenum<0||planenum>2) continue;
				              if(planenum==2) dec_mu_pida=float(pids[ipid]->PIDA());
			                  }
					  
					  if(dec_mu_pida<10 && dec_mu_pida!=-1){
					     vector<float> dec_muon_vector;
					     float dec_mu_st_x=pos.X();float dec_mu_st_y=pos.Y();float dec_mu_st_z=pos.Z();
					     float dec_mu_en_x=end.X();float dec_mu_en_y=end.Y();float dec_mu_en_z=end.Z();
					     dec_muon_vector.push_back(dec_mu_st_x);dec_muon_vector.push_back(dec_mu_st_y);dec_muon_vector.push_back(dec_mu_st_z);
					     dec_muon_vector.push_back(dec_mu_en_x);dec_muon_vector.push_back(dec_mu_en_y);dec_muon_vector.push_back(dec_mu_en_z);
					     art::Ptr<recob::Track> ptrack(trackListHandle,i);
			                     const recob::Track& track = *ptrack;
			                     TVector3 k_pos,k_end;
				             k_pos = track.Vertex();
		                             k_end = track.End();
					     float k_st_x=k_pos.X();float k_st_y=k_pos.Y();float k_st_z=k_pos.Z();float k_en_x=k_end.X();float k_en_y=k_end.Y();float k_en_z=k_end.Z();
					     vector<float> kaon_vector;
					     kaon_vector.push_back(k_st_x);kaon_vector.push_back(k_st_y);kaon_vector.push_back(k_st_z);
					     kaon_vector.push_back(k_en_x);kaon_vector.push_back(k_en_y);kaon_vector.push_back(k_en_z);
					     float k_dec_mu_min_dis=my_minimum_distance(kaon_vector,dec_muon_vector);
					     float k_dec_mu_index=my_close_distance_index(kaon_vector,dec_muon_vector);
						 
				             float k_st_dec_mu_st=TMath::Sqrt((k_st_x-dec_mu_st_x)*(k_st_x-dec_mu_st_x) + (k_st_y-dec_mu_st_y)*(k_st_y-dec_mu_st_y) + (k_st_z-dec_mu_st_z)*(k_st_z-dec_mu_st_z));
		                             float k_st_dec_mu_en=TMath::Sqrt((k_st_x-dec_mu_en_x)*(k_st_x-dec_mu_en_x) + (k_st_y-dec_mu_en_y)*(k_st_y-dec_mu_en_y) + (k_st_z-dec_mu_en_z)*(k_st_z-dec_mu_en_z));
		                             float k_en_dec_mu_st=TMath::Sqrt((k_en_x-dec_mu_st_x)*(k_en_x-dec_mu_st_x) + (k_en_y-dec_mu_st_y)*(k_en_y-dec_mu_st_y) + (k_en_z-dec_mu_st_z)*(k_en_z-dec_mu_st_z));
		                             float k_en_dec_mu_en=TMath::Sqrt((k_en_x-dec_mu_en_x)*(k_en_x-dec_mu_en_x) + (k_en_y-dec_mu_en_y)*(k_en_y-dec_mu_en_y) + (k_en_z-dec_mu_en_z)*(k_en_z-dec_mu_en_z));
					         
				             float av_common_open_angle_deg=-1;
					     
					     if(k_dec_mu_min_dis<=10){
					        if(k_dec_mu_index==0 || k_dec_mu_index==1){
				                   if(k_st_dec_mu_st<k_st_dec_mu_en){
			                              float average_common_point_x = float(k_st_x+dec_mu_st_x)/2;
	                                              float average_common_point_y = float(k_st_y+dec_mu_st_y)/2;
			                              float average_common_point_z = float(k_st_z+dec_mu_st_z)/2;
			                              float k_av_dis=TMath::Sqrt((k_en_x-average_common_point_x)*(k_en_x-average_common_point_x) + (k_en_y-average_common_point_y)*(k_en_y-average_common_point_y) + (k_en_z-average_common_point_z)*(k_en_z-average_common_point_z));
		                                      float mu_av_dis=TMath::Sqrt((dec_mu_en_x-average_common_point_x)*(dec_mu_en_x-average_common_point_x) + (dec_mu_en_y-average_common_point_y)*(dec_mu_en_y-average_common_point_y) + (dec_mu_en_z-average_common_point_z)*(dec_mu_en_z-average_common_point_z));
			                              float av_common_k_cos_alpha  = float(average_common_point_x-k_en_x)/k_av_dis; 
                                                      float av_common_k_cos_beta   = float(average_common_point_y-k_en_y)/k_av_dis; 
			                              float av_common_k_cos_gamma  = float(average_common_point_z-k_en_z)/k_av_dis; 
		                                      float av_common_mu_cos_alpha = float(dec_mu_en_x-average_common_point_x)/mu_av_dis; 
			                              float av_common_mu_cos_beta  = float(dec_mu_en_y-average_common_point_y)/mu_av_dis; 
	                                              float av_common_mu_cos_gamma = float(dec_mu_en_z-average_common_point_z)/mu_av_dis; 
			                              av_common_open_angle_deg = (TMath::ACos(av_common_k_cos_alpha*av_common_mu_cos_alpha + av_common_k_cos_beta*av_common_mu_cos_beta + av_common_k_cos_gamma*av_common_mu_cos_gamma))*57.3;
			                           } 
			          
				                   if(k_st_dec_mu_st>k_st_dec_mu_en){
			                              float average_common_point_x = float(k_st_x+dec_mu_en_x)/2;
		                                      float average_common_point_y = float(k_st_y+dec_mu_en_y)/2;
			                              float average_common_point_z = float(k_st_z+dec_mu_en_z)/2;
			                              float k_av_dis=TMath::Sqrt((k_en_x-average_common_point_x)*(k_en_x-average_common_point_x) + (k_en_y-average_common_point_y)*(k_en_y-average_common_point_y) + (k_en_z-average_common_point_z)*(k_en_z-average_common_point_z));
		                                      float mu_av_dis=TMath::Sqrt((dec_mu_st_x-average_common_point_x)*(dec_mu_st_x-average_common_point_x) + (dec_mu_st_y-average_common_point_y)*(dec_mu_st_y-average_common_point_y) + (dec_mu_st_z-average_common_point_z)*(dec_mu_st_z-average_common_point_z));
			                              float av_common_k_cos_alpha  = float(average_common_point_x-k_en_x)/k_av_dis; 
                                                      float av_common_k_cos_beta   = float(average_common_point_y-k_en_y)/k_av_dis; 
			                              float av_common_k_cos_gamma  = float(average_common_point_z-k_en_z)/k_av_dis; 
			                              float av_common_mu_cos_alpha = float(dec_mu_st_x-average_common_point_x)/mu_av_dis; 
			                              float av_common_mu_cos_beta  = float(dec_mu_st_y-average_common_point_y)/mu_av_dis; 
	                                              float av_common_mu_cos_gamma = float(dec_mu_st_z-average_common_point_z)/mu_av_dis; 
			                              av_common_open_angle_deg = (TMath::ACos(av_common_k_cos_alpha*av_common_mu_cos_alpha + av_common_k_cos_beta*av_common_mu_cos_beta + av_common_k_cos_gamma*av_common_mu_cos_gamma))*57.3;
			                          } 
				               }
				
				               if(k_dec_mu_index==2 || k_dec_mu_index==3){
			                          if(k_en_dec_mu_st<k_en_dec_mu_en){
			                             float average_common_point_x = float(k_en_x+dec_mu_st_x)/2;
	                                             float average_common_point_y = float(k_en_y+dec_mu_st_y)/2;
			                             float average_common_point_z = float(k_en_z+dec_mu_st_z)/2;
			                             float k_av_dis=TMath::Sqrt((k_st_x-average_common_point_x)*(k_st_x-average_common_point_x) + (k_st_y-average_common_point_y)*(k_st_y-average_common_point_y) + (k_st_z-average_common_point_z)*(k_st_z-average_common_point_z));
		                                     float mu_av_dis=TMath::Sqrt((dec_mu_en_x-average_common_point_x)*(dec_mu_en_x-average_common_point_x) + (dec_mu_en_y-average_common_point_y)*(dec_mu_en_y-average_common_point_y) + (dec_mu_en_z-average_common_point_z)*(dec_mu_en_z-average_common_point_z));
			                             float av_common_k_cos_alpha  = float(average_common_point_x-k_st_x)/k_av_dis; 
                                                     float av_common_k_cos_beta   = float(average_common_point_y-k_st_y)/k_av_dis; 
			                             float av_common_k_cos_gamma  = float(average_common_point_z-k_st_z)/k_av_dis; 
		                                     float av_common_mu_cos_alpha = float(dec_mu_en_x-average_common_point_x)/mu_av_dis; 
			                             float av_common_mu_cos_beta  = float(dec_mu_en_y-average_common_point_y)/mu_av_dis; 
	                                             float av_common_mu_cos_gamma = float(dec_mu_en_z-average_common_point_z)/mu_av_dis; 
			                             av_common_open_angle_deg = (TMath::ACos(av_common_k_cos_alpha*av_common_mu_cos_alpha + av_common_k_cos_beta*av_common_mu_cos_beta + av_common_k_cos_gamma*av_common_mu_cos_gamma))*57.3; 
			                          } 
			           
				                  if(k_en_dec_mu_st>k_en_dec_mu_en){
			                             float average_common_point_x = float(k_en_x+dec_mu_en_x)/2;
	                                             float average_common_point_y = float(k_en_y+dec_mu_en_y)/2;
			                             float average_common_point_z = float(k_en_z+dec_mu_en_z)/2;
			                             float k_av_dis=TMath::Sqrt((k_st_x-average_common_point_x)*(k_st_x-average_common_point_x) + (k_st_y-average_common_point_y)*(k_st_y-average_common_point_y) + (k_st_z-average_common_point_z)*(k_st_z-average_common_point_z));
		                                     float mu_av_dis=TMath::Sqrt((dec_mu_st_x-average_common_point_x)*(dec_mu_st_x-average_common_point_x) + (dec_mu_st_y-average_common_point_y)*(dec_mu_st_y-average_common_point_y) + (dec_mu_st_z-average_common_point_z)*(dec_mu_st_z-average_common_point_z));
			                             float av_common_k_cos_alpha  = float(average_common_point_x-k_st_x)/k_av_dis; 
                                                     float av_common_k_cos_beta   = float(average_common_point_y-k_st_y)/k_av_dis; 
			                             float av_common_k_cos_gamma  = float(average_common_point_z-k_st_z)/k_av_dis; 
		                                     float av_common_mu_cos_alpha = float(dec_mu_st_x-average_common_point_x)/mu_av_dis; 
			                             float av_common_mu_cos_beta  = float(dec_mu_st_y-average_common_point_y)/mu_av_dis; 
	                                             float av_common_mu_cos_gamma = float(dec_mu_st_z-average_common_point_z)/mu_av_dis;
			                             av_common_open_angle_deg = (TMath::ACos(av_common_k_cos_alpha*av_common_mu_cos_alpha + av_common_k_cos_beta*av_common_mu_cos_beta + av_common_k_cos_gamma*av_common_mu_cos_gamma))*57.3; 
			                         } 
			                      } 
						 
					      if(av_common_open_angle_deg<=150 && av_common_open_angle_deg!=-1){
					      
					         k_id_vec.push_back(i);
						 dec_mu_id_vec.push_back(j);
						 k_pid_vec.push_back(k_pida);
						 k_mu_min_dis_vec.push_back(k_dec_mu_min_dis);
						 k_mu_min_dis_index_vec.push_back(k_dec_mu_index);
					      
					      } // k/mu angle... 
					    } // k/mu minimum distance	 
					 } // dec mu pida cut
				       } // dec mu track len and track hit cut...
				    } // decay muon contained
			        } // loop over tracks to get decay muon
			      
			      //////////////////////////////////// End of decaying muon //////////////////////////////////////
			    
			    } // Kaon pida cut
			 } // k track len cut and track hit cut
		     } // kaon track contained
		 } // track loop to get Kaons
		 
		 /////////////////////////////////////////// looking into K/Mu information /////////////////////////////////
		 
		 if(k_mu_min_dis_vec.size()){
		    
		    float low_k_mu_dis=5e20;
		    int k_mu_low_index=-1;
		    int less_5_count=0;
		    
		    for(unsigned int k=0; k<k_mu_min_dis_vec.size(); k++){
		        if(k_mu_min_dis_vec[k]<low_k_mu_dis){
			   low_k_mu_dis=k_mu_min_dis_vec[k];
			   k_mu_low_index=k;
			}
			if(k_mu_min_dis_vec[k]<=5) less_5_count++;
		    }
		    
		    if(less_5_count>1){
		       std::vector<int> index_vec;
		       for(unsigned int k=0; k<k_mu_min_dis_vec.size(); k++){
		           if(k_mu_min_dis_vec[k]<=5) index_vec.push_back(k);
		       }
		    
		       float k_pida_diff=5e20;
		       for(int k=0; index_vec.size(); k++){
		           if(TMath::Abs(14-k_pid_vec[index_vec[k]])<k_pida_diff){
			      k_pida_diff=TMath::Abs(14-k_pid_vec[index_vec[k]]);
			      k_mu_low_index=index_vec[k];
			   }
		       }
		    
		    }
		    
		    pri_k_index=k_id_vec[k_mu_low_index];
		    dec_mu_index=dec_mu_id_vec[k_mu_low_index];
		    global_k_mu_min_dis_index=k_mu_min_dis_index_vec[k_mu_low_index];
		 
		 } // min dis info there
		 
		 //////////////////////////////////////// end of K/Mu information ////////////////////////////////////////
		 
		 
		 if(pri_k_index!=-1 && dec_mu_index!=-1){
		    rec_k_mu++;
		    std::vector<int> pri_mu_index_vector;
		    std::vector<float> pri_mu_trklen;
		    std::vector<float> k_pri_mu_min_dis;
		    for(int k=0; k<NTracks;++k){
		        if(k==pri_k_index || k==dec_mu_index) continue;
			std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(k);
			int pri_mu_hits_p2=0;
			int pri_mu_hits_p1=0;
			int pri_mu_hits_p0=0;
			
			for(unsigned int ical=0; ical<calos.size(); ++ical){
			    if(!calos[ical]) continue;
			    if(!calos[ical]->PlaneID().isValid) continue;
			    int planenum = calos[ical]->PlaneID().Plane;
			    if(planenum<0||planenum>2) continue; 
			    if(planenum==2) pri_mu_hits_p2=calos[ical]->dEdx().size();
			    if(planenum==1) pri_mu_hits_p1=calos[ical]->dEdx().size();
			    if(planenum==0) pri_mu_hits_p0=calos[ical]->dEdx().size();
			}
			
			
			if(pri_mu_hits_p2>=10 || pri_mu_hits_p1>=10 || pri_mu_hits_p0>=10){
			   
			   float pri_mu_pida=-1;
			   std::vector<art::Ptr<anab::ParticleID>> pids=fmpid.at(k);
			   for(unsigned int ipid=0; ipid<pids.size(); ++ipid){
			       if(!pids[ipid]->PlaneID().isValid) continue;
			       int planenum=pids[ipid]->PlaneID().Plane;
			       if(planenum<0||planenum>2) continue;
			       if(planenum==2) pri_mu_pida=float(pids[ipid]->PIDA());
			   }
			   
			   if(pri_mu_pida<10 && pri_mu_pida!=-1){
			      art::Ptr<recob::Track> ptrack(trackListHandle,k);
		              const recob::Track& track = *ptrack;
	                      TVector3 pos,end;
	                      pos = track.Vertex();
		              end = track.End();
			      float pri_mu_tlen=track.Length();
			      float mu_st_x=pos.X();float mu_st_y=pos.Y();float mu_st_z=pos.Z();float mu_en_x=end.X();float mu_en_y=end.Y();float mu_en_z=end.Z();
			      art::Ptr<recob::Track> ptrack_1(trackListHandle, pri_k_index);
			      const recob::Track& track_1 = *ptrack_1;
			      TVector3 pos_1,end_1;
			      pos_1 = track_1.Vertex();
		              end_1 = track_1.End();
			      float k_st_x=pos_1.X();float k_st_y=pos_1.Y();float k_st_z=pos_1.Z();float k_en_x=end_1.X();float k_en_y=end_1.Y();float k_en_z=end_1.Z();
			      float k_st_mu_st=TMath::Sqrt((k_st_x-mu_st_x)*(k_st_x-mu_st_x) + (k_st_y-mu_st_y)*(k_st_y-mu_st_y) + (k_st_z-mu_st_z)*(k_st_z-mu_st_z));
		              float k_st_mu_en=TMath::Sqrt((k_st_x-mu_en_x)*(k_st_x-mu_en_x) + (k_st_y-mu_en_y)*(k_st_y-mu_en_y) + (k_st_z-mu_en_z)*(k_st_z-mu_en_z));
		              float k_en_mu_st=TMath::Sqrt((k_en_x-mu_st_x)*(k_en_x-mu_st_x) + (k_en_y-mu_st_y)*(k_en_y-mu_st_y) + (k_en_z-mu_st_z)*(k_en_z-mu_st_z));
		              float k_en_mu_en=TMath::Sqrt((k_en_x-mu_en_x)*(k_en_x-mu_en_x) + (k_en_y-mu_en_y)*(k_en_y-mu_en_y) + (k_en_z-mu_en_z)*(k_en_z-mu_en_z));
			      
			      if(global_k_mu_min_dis_index==0 || global_k_mu_min_dis_index==1){
			         if(k_en_mu_st<=10 || k_en_mu_en<=10){
				    pri_mu_index_vector.push_back(k);
				    pri_mu_trklen.push_back(pri_mu_tlen);
				    if(k_en_mu_st<k_en_mu_en) k_pri_mu_min_dis.push_back(k_en_mu_st);
				    if(k_en_mu_en<k_en_mu_st) k_pri_mu_min_dis.push_back(k_en_mu_en);
				 }
			      }
			      
			      if(global_k_mu_min_dis_index==2 || global_k_mu_min_dis_index==3){
			         if(k_st_mu_st<10 || k_st_mu_en<10){
			            pri_mu_index_vector.push_back(k);
				    pri_mu_trklen.push_back(pri_mu_tlen);
				    if(k_st_mu_st<k_st_mu_en) k_pri_mu_min_dis.push_back(k_st_mu_st);
			            if(k_st_mu_en<k_st_mu_st) k_pri_mu_min_dis.push_back(k_st_mu_en);
			         }
			      }
			   } // pida cut on priamry muon
			} // trk hit cut on primary muon
		     } // loop over tracks to get primary mu
		     
		     ////////////////////////////////////// K/pri mu ///////////////////////////////////////////////////
		     
		     if(pri_mu_trklen.size()){
		     
		        float small_length=5e20;
		        int long_mu_index=-1;
		        int more_25_count=0;
			
			for(unsigned int k=0; k<pri_mu_trklen.size(); k++){
			    if(pri_mu_trklen[k]>small_length){
			       small_length=pri_mu_trklen[k];
			       long_mu_index=k;
			    }
			    if(pri_mu_trklen[k]>25) more_25_count++;
			}
			
			if(more_25_count>1){
			   std::vector<int> pri_mu_index_vec;
			   for(unsigned int k=0; k<pri_mu_trklen.size(); k++){
			       if(pri_mu_trklen[k]>25) pri_mu_index_vec.push_back(k);
			   }
			   
			   float large_dis=5e20;
			   
			   for(unsigned int k=0; k<pri_mu_index_vec.size(); k++){
			       if(k_pri_mu_min_dis[pri_mu_index_vec[k]]<large_dis){
			          large_dis=k_pri_mu_min_dis[pri_mu_index_vec[k]];
				  long_mu_index=pri_mu_index_vec[k];
			       } 
			   }
			   
			}
			
			pri_mu_index=pri_mu_index_vector[long_mu_index];
			if(pri_mu_index!=-1) rec_pri_mu++;
		     
		     } // pri mu trklen vec non empty...
		     
		     //////////////////////////////////// end of K/pri mu info ///////////////////////////////////////
		     
		  } // found match for pri k and dec mu
	       } // more than 3 trks and relible flash info exits
	   } // track info is saved
	   
	   if(pri_mu_index!=-1 && pri_k_index!=-1 && dec_mu_index!=-1) n_cc_reco++;
	   
	   //////////////////////////////////////////////////// Koan truth matching ////////////////////////////////
	       
	       int k_id=-1;
	       
	       if(pri_k_index!=-1){
	          std::vector< art::Ptr<recob::Hit> > allKHits = fmht.at(pri_k_index);
	          std::map<int,double> trk_k_ide;
	          for(size_t h = 0; h < allKHits.size(); ++h){
		      art::Ptr<recob::Hit> hit = allKHits[h];
		      std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(hit);
		      for(size_t e = 0; e < TrackIDs.size(); ++e){
		          trk_k_ide[TrackIDs[e].trackID] += TrackIDs[e].energy;
		      }
	           }
	      
	           double maxke = -1;
	           double totke = 0;
	           int Track_k_id = 0;
	      
	           for(std::map<int,double>::iterator ii=trk_k_ide.begin();ii!=trk_k_ide.end(); ++ii){
		       totke += ii->second;
		       if((ii->second)>maxke){
		           maxke = ii->second;
		           Track_k_id=ii->first;
		       }
	            }
	      
	            int k_origin=bt->TrackIDToMCTruth(Track_k_id)->Origin();
	            const simb::MCParticle* particle=bt->TrackIDToParticle(Track_k_id);
	            std::string pri("primary");
	            bool k_isPrimary=0;
	            k_isPrimary=particle->Process()==pri;
	            int k_pdg=particle->PdgCode();
	            float k_endE=particle->EndE()*1000;
	      
	            if(particle && k_origin==1 && k_isPrimary && k_endE<500 && k_pdg==321){
	               k_match=true;
		       k_id=particle->TrackId(); 
		       n_k_true_rec_match++;
	            }
	      
	        } // pri_k_inex!=-1
	      
	        ///////////////////////////////////////////////////////////////////////////////////////////////////////////
	      
	      
	        /////////////////////////// Primary muon truth matching ///////////////////////////////////////////////////
	      
	        if(pri_mu_index!=-1){
		   std::vector< art::Ptr<recob::Hit> > allpriMuHits = fmht.at(pri_mu_index);
	           std::map<int,double> trk_primu_ide;
	           for(size_t h = 0; h < allpriMuHits.size(); ++h){
		       art::Ptr<recob::Hit> hit = allpriMuHits[h];
		       std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(hit);
		       for(size_t e = 0; e < TrackIDs.size(); ++e){
		           trk_primu_ide[TrackIDs[e].trackID] += TrackIDs[e].energy;
		       }
	           }
	      
	           double maxprimue = -1;
	           double totprimue = 0;
	           int Track_primu_id = 0;
	      
	           for(std::map<int,double>::iterator ii=trk_primu_ide.begin();ii!=trk_primu_ide.end(); ++ii){
		       totprimue += ii->second;
		       if((ii->second)>maxprimue){
		           maxprimue = ii->second;
		           Track_primu_id=ii->first;
		       }
	           }
	      
	           int primu_origin=bt->TrackIDToMCTruth(Track_primu_id)->Origin();
	           const simb::MCParticle* particle=bt->TrackIDToParticle(Track_primu_id);
	           std::string pri("primary");
		   bool primu_isPrimary=0;
	           primu_isPrimary=particle->Process()==pri;
	           int primu_pdg=particle->PdgCode();
	      
	           if(particle && primu_origin==1 && TMath::Abs(primu_pdg)==13 && primu_isPrimary){ 
	              pri_mu_match=true;
		      n_primu_true_rec_match++;
	           }
	      
	      } // pri mu !=-1
	      
	      /////////////////////////////////////////////////////////////////////////////////////////////////////////
	      
	      /////////////////////////////////// Decaying muon matching ////////////////////////////////////////////
	      
	      if(dec_mu_index!=-1){
	         std::vector< art::Ptr<recob::Hit> > alldecMuHits = fmht.at(dec_mu_index);
	         std::map<int,double> trk_decmu_ide;
	         for(size_t h = 0; h < alldecMuHits.size(); ++h){
		     art::Ptr<recob::Hit> hit = alldecMuHits[h];
		     std::vector<sim::TrackIDE> TrackIDs = bt->HitToEveID(hit);
		     for(size_t e = 0; e < TrackIDs.size(); ++e){
		         trk_decmu_ide[TrackIDs[e].trackID] += TrackIDs[e].energy;
		     }
	         }
	      
	         double maxdecmue = -1;
	         double totdecmue = 0;
	         int Track_decmu_id = 0;
	      
	         for(std::map<int,double>::iterator ii=trk_decmu_ide.begin();ii!=trk_decmu_ide.end(); ++ii){
		     totdecmue += ii->second;
		     if((ii->second)>maxdecmue){
		         maxdecmue = ii->second;
		         Track_decmu_id=ii->first;
		     }
	         }
	      
	         const simb::MCParticle* particle=bt->TrackIDToParticle(Track_decmu_id);
	         int decmu_pdg=particle->PdgCode();
	      
	         if(particle && decmu_pdg==-13 && particle->Process()=="Decay" && particle->Mother()==k_id){
	            dec_mu_match=true;
		    n_decmu_true_rec_match++;
	         }
	      
	      } // dec mu !=-1
	      
	     //////////////////////////////////// End of decaying muon matching //////////////////////////////////
	      
	    if(k_match==true && pri_mu_match==true && dec_mu_match==true) n_true_events++;
	 
	 } // flash info saved
    } // not data
     
    fEventTree->Fill();
} // end of analyze function
	   
/////////////////// Defintion of reset function ///////////
void Kaonfiltercali2::reset(){
     run=-99999;
     subrun=-99999;
     event=-99999;
     rec_k_mu=-9999;
     rec_pri_mu=-9999;
     n_cc_reco=-9999;
     n_k_true_rec_match=-9999;
     n_primu_true_rec_match=-9999;
     n_decmu_true_rec_match=-9999;
     n_true_events=-9999;
}
//////////////////////// End of definition ///////////////

void Kaonfiltercali2::endSubRun(const art::SubRun& sr){ 
     art::Handle< sumdata::POTSummary > potListHandle;
     if(sr.getByLabel(fPOTModuleLabel,potListHandle))
        potbnb=potListHandle->totpot;
     fPOT->Fill(); 
} 		  

}
DEFINE_ART_MODULE(microboone::Kaonfiltercali2)

