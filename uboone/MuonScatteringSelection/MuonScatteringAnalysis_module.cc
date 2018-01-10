#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include <vector>
#include <string>

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "canvas/Persistency/Common/FindMany.h"

#include "uboone/MuonScatteringSelection/MuonScatteringAlg.h"

//#include "larsim/MCCheater/BackTracker.h"

#include "lardataobj/AnalysisBase/T0.h"

#include "TH1F.h"

namespace andy {
  class MuonScatteringAnalysis;
}

class andy::MuonScatteringAnalysis : public art::EDAnalyzer {
public:
  explicit MuonScatteringAnalysis(fhicl::ParameterSet const & p);
  virtual ~MuonScatteringAnalysis();

  void analyze(art::Event const & e) override;
//  void fillTruthVariables(art::Event const & e);

  void beginJob() override;

private:
  
  //TH1D *FlashCentre;
  TTree *OutputTree;
  TTree *SelectionTree;
  int event;
  int truthpdg;
  std::string process; 
  bool isRecoNu;
  bool isRecoCosmic;

  bool isRecoTruth;
  bool manyReco;

  int run;
  int subrun;
  TVector3 *vtxpos_reco;
  TVector3 *muonMom_reco;
  TVector3 *protonMom_reco;
  
  TVector3 *vtxpos_true;
  double t0;
  TVector3 *trueMomentum;


  int motherPDG;
  TVector3 *motherStartPosTPC;

  int selectionBranch[2];
  double closestTrackDist;
  int closestTrackID;
  double closestTrackStartY;
  double closestTrackEndY;
  
  int topTrackID;
  int bottomTrackID;
  double topTrackEnterY;
  double bottomTrackExitY;
  double closestEndDistBottom;
  double closestEndDistTop;
  
  int processIsMuonNuclear; // 1 - muon elastic; 0 - anything else

  MuonScatteringAlg algs;

  double vertexAssocCut;

//  std::string _t0TagName;
};


andy::MuonScatteringAnalysis::MuonScatteringAnalysis(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
//  std::cout << "reading parameter set" << std::endl;
  vertexAssocCut = p.get<double>("vertexAssocCut",50);
//  _t0TagName = p.get<std::string>("t0TagName");
  //this->reconfigure(p); 
}

andy::MuonScatteringAnalysis::~MuonScatteringAnalysis()
{
  // Clean up dynamic memory and other resources here.
  //a+=1;
  delete vtxpos_reco;
  delete muonMom_reco;
  delete protonMom_reco;
  
  delete vtxpos_true;
  delete trueMomentum;

}

void andy::MuonScatteringAnalysis::analyze(art::Event const & e)
{
  
  ////get the flash data
  //art::Handle< std::vector<recob::OpFlash> > flashHandle;
  //std::string _flash_producer_name = "opflashBeam";
  //e.getByLabel(_flash_producer_name,flashHandle);
  //if(!flashHandle.isValid()) {
  //    std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve recob::OpFlash from "
  //    << _flash_producer_name << std::endl;
  //    throw std::exception();
  //}

  // get MCtruths objects
  art::Handle< std::vector<simb::MCTruth> > truthHandle;
  std::string _mctruth_name = "generator";
  e.getByLabel(_mctruth_name, truthHandle);
  if(!truthHandle.isValid()) {
      std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve simb::MCTruth from "
      << _mctruth_name << std::endl;
      throw std::exception();
  }

  // get reco nu tracks
  std::string _nuTrack_name = "pandoraNu";
  auto nuTrackHandle = e.getValidHandle< std::vector< recob::Track > >(_nuTrack_name);
//  art::Handle< std::vector<recob::Track> > nuTrackHandle;
//  e.getByLabel(_nuTrack_name, nuTrackHandle);
//  if(!nuTrackHandle.isValid()) {
//      std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve recob::Track from "
//      << _nuTrack_name << std::endl;
//      throw std::exception();
//  }

  // get reco cosmic tracks
  art::Handle< std::vector<recob::Track> > cosTrackHandle;
  std::string _cosTrack_name = "pandoraCosmic";
  e.getByLabel(_cosTrack_name, cosTrackHandle);
  if(!cosTrackHandle.isValid()) {
      std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve recob::Track from "
      << _cosTrack_name << std::endl;
      throw std::exception();
  }
  
  // get reco vertixes
  art::Handle< std::vector<recob::Vertex> > vertexHandle;
  e.getByLabel(_nuTrack_name, vertexHandle);
  if(!vertexHandle.isValid()) {
      std::cerr << "\033[93m[ERROR]\033[00m Could not retrieve recob::Vertex from "
      << _nuTrack_name << std::endl;
      throw std::exception();
  }
  
  // get MCParticles
  auto mcParticleHandle = e.getValidHandle< std::vector< simb::MCParticle> >("largeant"); // check the name here!

  // Truth matching objects
  //auto t0TagHandle = e.getValidHandle< std::vector< anab::T0 > >(_t0TagName);
  art::FindManyP<anab::T0> FindMatchT0(nuTrackHandle, e, "mctrutht0matching");
  art::FindManyP<anab::T0> FindMatchT0Cosmic(cosTrackHandle, e, "mctrutht0matchingCosmic");

//  std::vector<recob::OpFlash> const& flashVector(*flashHandle);
  std::vector<recob::Track> const& nuTrackVector(*nuTrackHandle);
  std::vector<recob::Track> const& cosTrackVector(*cosTrackHandle);
  std::vector<recob::Vertex> const& vertexVector(*vertexHandle);

  std::vector<int> vertexCandidateList(vertexVector.size(),0);
//  int maxTracks = std::max(nuTrackVector.size(), cosTrackVector.size());
  std::vector< std::vector<int> > nuTrackCandidateList(vertexVector.size(),std::vector<int>(nuTrackVector.size(),0));
  std::vector< std::vector<int> > cosTrackCandidateList(vertexVector.size(),std::vector<int>(cosTrackVector.size(),0));
  
  int mcpart_i(-1);
  for (auto const mcpart : *mcParticleHandle ){ // Loop through MCParticles bare
    mcpart_i++;
//    std::vector<art::Ptr<anab::T0>> const& T0tags = FindRecoFromTruth.at(mcpart_i);

    if (mcpart.PdgCode() == 2212){
//      if (T0tags.size() < 1){isRecoTruth=false;}
//      if (T0tags.size() > 1){isRecoTruth=true; manyReco=true;}
      
      truthpdg = 2212;
      isRecoNu = false;
      isRecoCosmic = false;
      processIsMuonNuclear = 0;
      if (mcpart.Process() == "muonNuclear"){
        processIsMuonNuclear=1;
      }
      TVector3 tmpmom = mcpart.Momentum().Vect();
      trueMomentum->SetXYZ(tmpmom.X(), tmpmom.Y(), tmpmom.Z());
      TVector3 tmppos = mcpart.Position().Vect();
      vtxpos_true->SetXYZ(tmppos.X(), tmppos.Y(), tmppos.Z());
      t0 = mcpart.Position().T();
      auto motherId = mcpart.Mother();
      for (auto const mcpartMum : *mcParticleHandle ){ // Loop through MCParticles for mother
        if (mcpartMum.TrackId() == motherId){
          motherPDG = mcpartMum.PdgCode();
          auto tmpVct = algs.TPCentryPoint(mcpartMum);
          motherStartPosTPC->SetXYZ(tmpVct.X(),tmpVct.Y(),tmpVct.Z() );

          double motherMomentumAtScatter = 0;
          break;
        }
      }

      OutputTree->Fill();
    }
  }

  for (unsigned int i(0); i < vertexVector.size(); ++i){
    double xyz[3];
    vertexVector.at(i).XYZ(xyz); // read the vertex position into the array (this is insane)
    vtxpos_reco->SetXYZ(xyz[0], xyz[1], xyz[2]);
    if (algs.inFV(xyz)){
      vertexCandidateList.at(i)=1;
    }
    for (unsigned int j(0); j < nuTrackVector.size(); ++j){
      if (nuTrackVector.at(j).NPoints() < 2){
        continue;}
      if (not algs.fullyContained(nuTrackVector.at(j))){
        continue;
      }
      if (algs.closestApproach(nuTrackVector.at(j), xyz) < vertexAssocCut){
        nuTrackCandidateList.at(i).at(j)=1;
        std::vector<art::Ptr<anab::T0>> const& T0tags = FindMatchT0.at(j);
        if (T0tags.size() > 1){
//          std::cout << "more than one match!" << std::endl;
          continue;
        }
        else if (T0tags.size() < 1){
          continue;
        }
        for (auto const mcpart : *mcParticleHandle ){
          if ((*(T0tags.at(0))).TriggerBits() == mcpart.TrackId() ){ //matched
            int pdg = mcpart.PdgCode();
            truthpdg = pdg;
            isRecoNu = true;
            isRecoCosmic = false;
            TVector3 tmpmom = mcpart.Momentum().Vect();
            trueMomentum->SetXYZ(tmpmom.X(), tmpmom.Y(), tmpmom.Z());
            TVector3 tmppos = mcpart.Position().Vect();
            vtxpos_true->SetXYZ(tmppos.X(), tmppos.Y(), tmppos.Z());
            t0 = mcpart.Position().T();
            if (mcpart.Process() == "muonNuclear"){
              processIsMuonNuclear=1;
            }
            OutputTree->Fill();
          }
        }
      }
    }
    for (unsigned int j(0); j < cosTrackVector.size(); ++j){
      if (cosTrackVector.at(j).NPoints() < 2){
        continue;}
      if (not algs.fullyContained(cosTrackVector.at(j))){
        continue;
      }
      if (algs.closestApproach(cosTrackVector.at(j), xyz) < vertexAssocCut){
        cosTrackCandidateList.at(i).at(j)=1;
        std::vector<art::Ptr<anab::T0>> const& T0tags = FindMatchT0Cosmic.at(j);
        if (T0tags.size() > 1){
//          std::cout << "more than one match!" << std::endl;
          continue;
        }
        else if (T0tags.size() < 1){
//          std::cout << "no matches" << std::endl;
          continue;
        }
        for (auto const mcpart : *mcParticleHandle ){
          if ((*(T0tags.at(0))).TriggerBits() == mcpart.TrackId() ){ //matched
            int pdg = mcpart.PdgCode();
            truthpdg = pdg;
            isRecoNu = false;
            isRecoCosmic = true;
            TVector3 tmpmom = mcpart.Momentum().Vect();
            trueMomentum->SetXYZ(tmpmom.X(), tmpmom.Y(), tmpmom.Z());
            TVector3 tmppos = mcpart.Position().Vect();
            t0 = mcpart.Position().T();
            vtxpos_true->SetXYZ(tmppos.X(), tmppos.Y(), tmppos.Z());
            if (mcpart.Process() == "muonNuclear"){
              processIsMuonNuclear=1;
            }
            OutputTree->Fill();
          }
        }
      }
    }
  }
//---------------------------------------------------------------------//
// Actual selection code
//---------------------------------------------------------------------//
//  closestTrackDist(9999);
//  closestTrackID(-9999);
//  std::string process;
  for (unsigned int i(0); i < nuTrackVector.size(); ++i){
    if (nuTrackVector.at(i).NPoints() < 2){
      continue;}
    if (not algs.fullyContained(nuTrackVector.at(i))){
      continue;
    }
    double startPos[3] = { nuTrackVector.at(i).Vertex().X(), nuTrackVector.at(i).Vertex().Y(), nuTrackVector.at(i).Vertex().Z() };
    double endPos[3] = { nuTrackVector.at(i).End().X(), nuTrackVector.at(i).End().Y(), nuTrackVector.at(i).End().Z() };

// Get true information about proton candidate here
    std::vector<art::Ptr<anab::T0>> const& T0tags = FindMatchT0.at(i);
    if (T0tags.size() < 1){std::cout << "proton candidate has no truth match" << std::endl;}
    else{
      int mcpartid = (*(T0tags.at(0))).TriggerBits();
      for (auto const mcpart : *mcParticleHandle ){ // Loop through MCParticles bare
        if ( mcpartid == mcpart.TrackId() ){ //matched
          truthpdg = mcpart.PdgCode();
          process = mcpart.Process();
          TVector3 tmpmom = mcpart.Momentum().Vect();
          trueMomentum->SetXYZ(tmpmom.X(), tmpmom.Y(), tmpmom.Z());
          TVector3 tmppos = mcpart.Position().Vect();
          vtxpos_true->SetXYZ(tmppos.X(), tmppos.Y(), tmppos.Z());
          vtxpos_reco->SetXYZ(startPos[0], startPos[1], startPos[2]);
        }
      }
    }
    processIsMuonNuclear = 0;
    if ( process == "muonNuclear"){processIsMuonNuclear = 1;}
// reset whether the track is tagged by one or other selection branch
    for (int branch(0); branch < 2; branch++){
      selectionBranch[branch]=0;
    } 

// Test event selection (0)
// Match this pandoraNu track start/end within X cm of pandoraCosmic track bulk
    closestTrackDist = 9999;
    closestTrackID = -9999;
    closestTrackStartY = -9999;
    closestTrackEndY = -9999;
    for (unsigned int j(0); j < cosTrackVector.size(); ++j){
      if (cosTrackVector.at(j).NPoints() < 2){
        continue;}
      double closestApproachStart = algs.closestApproach(cosTrackVector.at(j), startPos);
      double closestApproachEnd = algs.closestApproach(cosTrackVector.at(j), endPos);
      if ((closestApproachStart < vertexAssocCut) or (closestApproachEnd < vertexAssocCut)){
        double closestTrackDistThisNu = (closestApproachStart < closestApproachEnd ? closestApproachStart : closestApproachEnd);
        if (closestTrackDistThisNu < closestTrackDist){
          closestTrackDist = closestTrackDistThisNu;
          closestTrackID = cosTrackVector.at(j).ID();
          closestTrackStartY = cosTrackVector.at(j).Vertex().Y();
          closestTrackEndY = cosTrackVector.at(j).End().Y();
          selectionBranch[0] = 1;
        }
        else continue;
      } // end if close to track
    } // cosTrack loop end

// Test event selection (1)
// Match a pandoraNu track start/end within X cm of pandoraNu track start/end, which enters/leaves through top/bottom
// Then also require a second track going out the other top/bottom
    topTrackID = -9999;
    bottomTrackID = -9999;
    topTrackEnterY = -9999;
    bottomTrackExitY = -9999;
    closestEndDistBottom = 9999;
    closestEndDistTop = 9999;
    for (unsigned int j(0); j < nuTrackVector.size(); ++j){
      if (nuTrackVector.at(j).NPoints() < 2){
        continue;}
      if (i==j){continue;} // this would literally be matching the proton candidate to itself!
      double startPosCos[3] = {nuTrackVector.at(j).Vertex().X(), nuTrackVector.at(j).Vertex().Y(), nuTrackVector.at(j).Vertex().Z()};
      double endPosCos[3] = {nuTrackVector.at(j).End().X(), nuTrackVector.at(j).End().Y(), nuTrackVector.at(j).End().Z()};
      bool upwards = (startPosCos[1] < endPosCos[1] ? true : false);
      bool aboveProtonCand;
      if (upwards){aboveProtonCand = (endPosCos[1] > endPos[1] ? true : false);} // check if end position of muon segment candidate is above proton candidate (reco up)
      else {aboveProtonCand = (startPosCos[1] > endPos[1] ? true : false);} // check if start position of muon segment candidate is above proton candidate (reco down)
      
      double *vertexEnd;
      double *boundaryEnd;
      if (aboveProtonCand){
        if (upwards){vertexEnd = startPosCos; boundaryEnd = endPosCos;}
        else {vertexEnd = endPosCos; boundaryEnd = startPosCos;}
      }
      else{
        if (upwards) { vertexEnd = endPosCos; boundaryEnd = startPosCos;}
        else { vertexEnd = startPosCos; boundaryEnd = endPosCos;}
      } // worked out which end to compare to what
      // find if we want the proton forwars or backwards
      double thisClosestEndDist = (algs.absDistance(startPos, vertexEnd) < algs.absDistance(endPos, vertexEnd) ? algs.absDistance(startPos, vertexEnd) : algs.absDistance(endPos, vertexEnd));
      // tracking top/bottom muon candidate track distances from proton candidate
      if (aboveProtonCand and thisClosestEndDist > closestEndDistTop){continue;}
      if ((not aboveProtonCand) and thisClosestEndDist > closestEndDistBottom){continue;}
      if (aboveProtonCand){closestEndDistTop = thisClosestEndDist;}
      if (not aboveProtonCand){closestEndDistBottom = thisClosestEndDist;}
      // If the end is closest enough
      if ( thisClosestEndDist < vertexAssocCut){
        if (abs(boundaryEnd[1]) > 100){ // TODO - define y better
          if (aboveProtonCand){
            topTrackID = nuTrackVector.at(j).ID();
            topTrackEnterY = boundaryEnd[1];
          }
          else {
            bottomTrackID = nuTrackVector.at(j).ID();
            bottomTrackExitY = boundaryEnd[1];
          }

        } // other track end is close to a boundary
      } // track is close to proton candidate end
    } // nuTrack (cosmic candidate) loop end
    selectionBranch[1] = 1;
// No more ideas for selecting tracks!
  } // nuTrack loop end
  SelectionTree->Fill();
  
//  fillTruthVariables(art::Event e);

}

//void andy::MuonScatteringAnalysis::fillTruthVariables(art::Event const & e)
//{
//  
//}

void andy::MuonScatteringAnalysis::beginJob()
{
  //std::cout << "calling beginJob" << std::endl;
  art::ServiceHandle<art::TFileService> tfs;

  vtxpos_reco = new TVector3();
  muonMom_reco = new TVector3();
  protonMom_reco = new TVector3();
  
  vtxpos_true = new TVector3();
  trueMomentum = new TVector3();

  FlashCentre = tfs->make<TH1D>("FlashCentre"      , "FlashCentre"        , 30, -50, 100050);
  OutputTree = tfs->make<TTree>("truthTree","truthTree");
//  OutputTree->Branch("event",&event,"event/I");
  OutputTree->Branch("truthpdg",&truthpdg,"truthpdg/I");
  OutputTree->Branch("isRecoNu",&isRecoNu,"isRecoNu/B");
  OutputTree->Branch("isRecoCosmic",&isRecoCosmic,"isRecoCosmic/B");
  OutputTree->Branch("isRecoTruth",&isRecoTruth,"isRecoTruth/B");
  OutputTree->Branch("manyReco",&manyReco,"manyReco/B");

  OutputTree->Branch("processIsMuonNuclear",&processIsMuonNuclear,"processIsMuonNuclear/O");
  OutputTree->Branch("motherPDG",&motherPDG,"motherPDG/I");
  OutputTree->Branch("motherStartPosTPC","TVector3",&motherStartPosTPC);

  //OutputTree->Branch("subrun",&subrun,"subrun/I");
  //OutputTree->Branch("run",&run,"run/I");
  OutputTree->Branch("vtxpos_reco","TVector3",&vtxpos_reco);
  //OutputTree->Branch("muonMom_reco","TVector3",&muonMom_reco);
  OutputTree->Branch("protonMom_reco","TVector3",&protonMom_reco);
  
  OutputTree->Branch("vtxpos_true","TVector3",&vtxpos_true);
  OutputTree->Branch("t0",&t0,"t0/D");
  OutputTree->Branch("trueMomentum","TVector3",&trueMomentum);
  
  SelectionTree = tfs->make<TTree>("muonSelectionTree","muonSelectionTree");
  SelectionTree->Branch("truthpdg",&truthpdg,"truthpdg/I");
  SelectionTree->Branch("process", &process);
  SelectionTree->Branch("processIsMuonNuclear",&processIsMuonNuclear,"processIsMuonNuclear/O");
  SelectionTree->Branch("protonMom_reco","TVector3",&protonMom_reco);
  SelectionTree->Branch("trueMomentum","TVector3",&trueMomentum);
  SelectionTree->Branch("vtxpos_true","TVector3",&vtxpos_true);
  SelectionTree->Branch("vtxpos_reco","TVector3",&vtxpos_reco);
  SelectionTree->Branch("selectionBranch",&selectionBranch,"selectionBranch[2]/I");
  
  // selectionBranch[0] variables
  SelectionTree->Branch("closestTrackDist",&closestTrackDist,"closestTrackDist/D");
  SelectionTree->Branch("closestTrackID",&closestTrackID,"closestTrackID/I");
  SelectionTree->Branch("closestTrackStartY",&closestTrackStartY,"closestTrackStartY/D");
  SelectionTree->Branch("closestTrackEndY",&closestTrackEndY,"closestTrackEndY/D");

  // selectionBranch[1] variables
  SelectionTree->Branch("topTrackID",&topTrackID,"topTrackID/I");
  SelectionTree->Branch("topTrackEnterY",&topTrackEnterY,"topTrackEnterY/D");
  SelectionTree->Branch("bottomTrackID",&bottomTrackID,"bottomTrackID/I");
  SelectionTree->Branch("bottomTrackExitY",&bottomTrackExitY,"bottomTrackExitY/D");
  SelectionTree->Branch("closestEndDistTop",&closestEndDistTop,"closestEndDistTop/D");
  SelectionTree->Branch("closestEndDistBottom",&closestEndDistBottom,"closestEndDistBottom/D");
  

}

DEFINE_ART_MODULE(andy::MuonScatteringAnalysis)
