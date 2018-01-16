#ifndef TRACKFEATURES_CXX
#define TRACKFEATURES_CXX

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"


#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/OpDetGeo.h"
#include "uboone/Geometry/UBOpReadoutMap.h"

#include "larsim/MCCheater/BackTracker.h"

#include "TMath.h"

#include <iostream>

#include "uboone/DecisionTreeID/Algorithms/TrackFeatures.h"

namespace dtfeatures {

  TrackFeatures::TrackFeatures()
  {
    _nendhits      = 3;
    _beamwinstartt = 3.2;
    _beamwinendt   = 4.8;

  }

  void TrackFeatures::Configure(fhicl::ParameterSet const& p)
  {
    _nendhits      = p.get<unsigned int>("NEndHits");
    _beamwinstartt = p.get<float>("BeamWindowStartT");
    _beamwinendt   = p.get<float>("BeamWindowEndT");

    _pidassoclabel        = p.get<std::string>("PIDAssocLabel");
    _caloassoclabel       = p.get<std::string>("CaloAssocLabel");
    _cosmictagassoclabel  = p.get<std::string>("CosmicTagAssocLabel");
    _containtagassoclabel = p.get<std::string>("ContainTagAssocLabel");
    _hitassoclabel        = p.get<std::string>("HitAssocLabel");
    _flashmodulelabel     = p.get<std::string>("FlashModuleLabel");
  }

  std::vector<std::vector<float>> TrackFeatures::CreateFeatures(art::Event const& e, 
                                       art::Handle<std::vector<recob::Track>>& trackVecHandle )
  {
    // recover associations
    art::FindManyP<anab::ParticleID> pidAssns(trackVecHandle, e, _pidassoclabel);

    art::FindManyP<anab::Calorimetry> caloAssns(trackVecHandle, e, _caloassoclabel);
    art::FindManyP<anab::CosmicTag> cosmictagAssns(trackVecHandle, e, _cosmictagassoclabel); 
    art::FindManyP<anab::CosmicTag> containtagAssns(trackVecHandle, e, _containtagassoclabel);
    // recover flash handle
    art::Handle< std::vector<recob::OpFlash> > flashVecHandle;
    e.getByLabel(_flashmodulelabel, flashVecHandle);

    std::vector< std::vector<float> > evtdata;
    if(trackVecHandle.isValid())
    {
      // Loop over tracks in event and get features
      for(size_t trkIdx = 0; trkIdx < trackVecHandle->size(); trkIdx++)
      {
        art::Ptr<recob::Track> trk(trackVecHandle,trkIdx);

        std::vector<float> trkdata;
        
        float  _minflashpe = 6.5;
        
        std::vector<art::Ptr<anab::ParticleID> > pidVec    = pidAssns.at(trk.key());
        std::vector<art::Ptr<anab::Calorimetry> > caloVec  = caloAssns.at(trk.key());
        std::vector<art::Ptr<anab::CosmicTag> > cosmicVec  = cosmictagAssns.at(trk.key());
        std::vector<art::Ptr<anab::CosmicTag> > containVec = containtagAssns.at(trk.key());
        
        int   nhits         = 0;
        float length        = -9999.;
        float starty        = -9999.;
        float startz        = -9999.;
        float endy          = -9999.;
        float endz          = -9999.;
        float theta         = -9999.;
        float phi           = -9999.;
        float distlenratio  = -9999.;
        float startdqdx     = -9999.;
        float enddqdx       = -9999.;
        float dqdxdiff      = -9999.;
        float dqdxratio     = -9999.;
        float totaldqdx     = -9999.;
        float averagedqdx   = -9999.;
        float cosmicscore   = -9999.;
        float coscontscore  = -9999.;
        float pidpida       = -9999.;
        float pidchipr      = -9999.;
        float pidchika      = -9999.;
        float pidchipi      = -9999.;
        float pidchimu      = -9999.;
        float cfdistance    = 9999.;
       
        float tmpphi;
        float tmpstartdqdx;
        float tmpenddqdx;
        float startx;
        float endx;
        float dist;
        
        bool tflip          = false;

        size_t bstplane = 2;
        int bstplanehits = 0;

        if(!caloVec.empty())
        {
          
          std::vector<double> dqdx;
          nhits     = 0;
          totaldqdx = 0.;
          tmpstartdqdx = 0.;
          tmpenddqdx   = 0.;
          for(size_t iplane = 0; iplane < caloVec.size(); iplane++)
          {
            dqdx  = caloVec.at(iplane)->dQdx();
            if(dqdx.empty()) continue;

            int pnhits = dqdx.size();
            if(pnhits > bstplanehits)
            {
              bstplanehits = pnhits;
              bstplane     = caloVec.at(iplane)->PlaneID().Plane;
            }
            nhits += pnhits;
            for( auto& dq : dqdx) totaldqdx += dq;
            size_t endhits;
            if(dqdx.size() >= _nendhits)
            {
              endhits = _nendhits;
            } else {
              endhits = TMath::Floor((float)_nendhits/2.);
            }

            for(size_t ihd=0; ihd < endhits; ihd++)
            {
              tmpenddqdx   += dqdx.at(ihd);
              tmpstartdqdx += dqdx.at(pnhits - (1+ihd));
            }
          }
          averagedqdx = (float)totaldqdx / (float)nhits; 
          if(tmpenddqdx < tmpstartdqdx)
          {
            tflip     = true;
            enddqdx   = tmpstartdqdx;
            startdqdx = tmpenddqdx;
          } else {
            tflip     = false;
            enddqdx   = tmpenddqdx;
            startdqdx = tmpstartdqdx;
          }
          dqdxdiff  = enddqdx - startdqdx;
          dqdxratio = enddqdx/startdqdx;
        }
        
        length = trk->Length();
        if(!tflip)
        {
          startx = trk->Vertex().X();
          starty = trk->Vertex().Y();
          startz = trk->Vertex().Z();
          endx   = trk->End().X();
          endy   = trk->End().Y();
          endz   = trk->End().Z();
          theta  = trk->VertexDirection().Theta();
          phi    = trk->VertexDirection().Phi();
        } else {
          startx = trk->End().X();
          starty = trk->End().Y();
          startz = trk->End().Z();
          endx   = trk->Start().X();
          endy   = trk->Start().Y();
          endz   = trk->Start().Z();
          theta  = TMath::Pi() - trk->EndDirection().Theta();
          tmpphi = trk->EndDirection().Phi();
          phi    = (tmpphi > 0 ? tmpphi - TMath::Pi() :  tmpphi + TMath::Pi());
        }
        
        dist = TMath::Sqrt( TMath::Power(endx-startx,2) 
                          + TMath::Power(endy-starty,2) 
                          + TMath::Power(endz-startz,2));
        distlenratio = dist/length;
        
        if(!cosmicVec.empty()) cosmicscore = cosmicVec.at(0)->CosmicScore();
        if(!containVec.empty()) coscontscore = containVec.at(0)->CosmicScore();
        if(!pidVec.empty())
        {
          for(size_t iplane = 0; iplane < pidVec.size(); iplane++)
          {
            if(pidVec[iplane]->PlaneID().Plane == bstplane)
            {
              pidpida  = pidVec[iplane]->PIDA();
              pidchipr = pidVec[iplane]->Chi2Proton();
              pidchika = pidVec[iplane]->Chi2Kaon();
              pidchipi = pidVec[iplane]->Chi2Pion();
              pidchimu = pidVec[iplane]->Chi2Muon();
            }
          }
        }
        
        // find closest flash
        if(flashVecHandle.isValid())
        {
          float trkzcenter = (startz + endz)/2.;
          for(size_t flIdx = 0; flIdx < flashVecHandle->size(); flIdx++)
          {
            art::Ptr<recob::OpFlash> flash(flashVecHandle,flIdx);
            if(flash->Time() < _beamwinstartt || flash->Time() > _beamwinendt || flash->TotalPE() < _minflashpe) continue; 
            if(TMath::Abs(flash->ZCenter() - trkzcenter) < cfdistance) 
            {
              cfdistance = TMath::Abs(flash->ZCenter() - trkzcenter);
            }
          }
        }
        
        trkdata.clear();
        
        trkdata.emplace_back(nhits);
        trkdata.emplace_back(length);
        trkdata.emplace_back(starty);
        trkdata.emplace_back(startz);
        trkdata.emplace_back(endy);
        trkdata.emplace_back(endz);
        trkdata.emplace_back(theta);
        trkdata.emplace_back(phi);
        trkdata.emplace_back(distlenratio);
        trkdata.emplace_back(startdqdx);
        trkdata.emplace_back(enddqdx);
        trkdata.emplace_back(dqdxdiff);
        trkdata.emplace_back(dqdxratio);
        trkdata.emplace_back(totaldqdx);
        trkdata.emplace_back(averagedqdx);
        trkdata.emplace_back(cosmicscore);
        trkdata.emplace_back(coscontscore);
        trkdata.emplace_back(pidpida);
        trkdata.emplace_back(pidchipr);
        trkdata.emplace_back(pidchika);
        trkdata.emplace_back(pidchipi);
        trkdata.emplace_back(pidchimu);
        trkdata.emplace_back(cfdistance);

        evtdata.emplace_back(trkdata);
      }
    }
    
    return evtdata; 
  }

  std::vector<float> TrackFeatures::ClassLabel(art::Event const& e,
                                       art::Handle<std::vector<recob::Track>>& trackVecHandle )
  {
    std::vector<float> trklabelvec;

    // recover hits associated with tracks
    art::FindManyP<recob::Hit> hitAssns(trackVecHandle, e, _hitassoclabel);

    if(trackVecHandle.isValid())
    {
      if(e.isRealData())
      {
        size_t num_tracks = trackVecHandle->size();
        std::vector<float> tmplabelvec(num_tracks, -1.);
        return tmplabelvec;
      }

      art::ServiceHandle<cheat::BackTracker> bt;

      // loop over tracks and assign class label
      for(size_t trkIdx = 0; trkIdx < trackVecHandle->size(); trkIdx++)
      {
        art::Ptr<recob::Track> trk(trackVecHandle,trkIdx);
        std::vector<art::Ptr<recob::Hit>> hitVec = hitAssns.at(trk.key());
        
        std::vector<int> labelVec(5,0);
        
        for(auto const& hit : hitVec)
        {
          std::vector<sim::TrackIDE> g4ides;
          g4ides = bt->HitToTrackID(hit);
        
          if(g4ides.empty()) continue;    
          const simb::MCParticle* mcparticle = bt->TrackIDToParticle( g4ides.at(0).trackID );
          if(!mcparticle) continue;
          art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth( g4ides.at(0).trackID );
          if(mctruth.isNull()) continue;
        
          if(mcparticle->PdgCode() == 2212)
          { labelVec[0] = labelVec[0] + 1; }
          else if(mctruth->Origin() == 1 && TMath::Abs(mcparticle->PdgCode()) == 13)
          { labelVec[1] = labelVec[1] + 1; }
          else if(mctruth->Origin() == 1 && TMath::Abs(mcparticle->PdgCode()) == 211)
          { labelVec[2] = labelVec[2] + 1; }
          else if(mctruth->Origin() == 1) 
          { labelVec[3] = labelVec[3] + 1; }
          else if(mctruth->Origin() == 2 && mcparticle->PdgCode() != 2212)
          { labelVec[4] = labelVec[4] + 1; }
        
        }
        // PDG index with the most hits wins
        int maxidx = (int)(max_element(labelVec.begin(),labelVec.end()) - labelVec.begin());

        // if no mc tracks were found give label -1
        if(labelVec[maxidx] == 0) maxidx = -1;

        trklabelvec.emplace_back( (float)maxidx ); 
      }
    } 
    return trklabelvec;
  }

  std::vector<anab::CosmicTagID_t> TrackFeatures::CosmicTagLabel()
  {
    std::vector<anab::CosmicTagID_t> ctlabelvec;

    const anab::CosmicTagID_t TAGID_P  = anab::CosmicTagID_t::kGeometry_YY;
    const anab::CosmicTagID_t TAGID_MU = anab::CosmicTagID_t::kGeometry_YZ;
    const anab::CosmicTagID_t TAGID_PI = anab::CosmicTagID_t::kGeometry_ZZ;
    const anab::CosmicTagID_t TAGID_EM = anab::CosmicTagID_t::kGeometry_XX;
    const anab::CosmicTagID_t TAGID_CS = anab::CosmicTagID_t::kGeometry_XY;

    ctlabelvec.clear();
    ctlabelvec.emplace_back(TAGID_P);    
    ctlabelvec.emplace_back(TAGID_MU);    
    ctlabelvec.emplace_back(TAGID_PI);    
    ctlabelvec.emplace_back(TAGID_EM);    
    ctlabelvec.emplace_back(TAGID_CS);    

    return ctlabelvec;
  }
}

#endif
