////////////////////////////////////////////////////////////////////////
// Class:       TriggerPrimitiveHitMaker
// Plugin Type: producer (art v2_11_02)
// File:        TriggerPrimitiveHitMaker_module.cc
//
// Generated at Mon Jun  3 13:29:25 2019 by Claire Hinrichs using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "TH1D.h"
#include <string>
#include "TCanvas.h"
#include <iostream>
#include <memory>
#include <cmath>
#include <algorithm>
#include "TF1.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "TH2D.h"
#include "TF2.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "lardata/ArtDataHelper/HitCreator.h"
#include "larreco/HitFinder/HitFinderTools/ICandidateHitFinder.h"
#include "larreco/HitFinder/HitFinderTools/IPeakFitter.h"
#include "TGraphErrors.h"
#include "TDecompSVD.h"
#include "TMath.h"
#include "TStopwatch.h"

class TriggerPrimitiveHitMaker;


class TriggerPrimitiveHitMaker : public art::EDProducer {
public:
  explicit TriggerPrimitiveHitMaker(fhicl::ParameterSet const & pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TriggerPrimitiveHitMaker(TriggerPrimitiveHitMaker const &) = delete;
  TriggerPrimitiveHitMaker(TriggerPrimitiveHitMaker &&) = delete;
  TriggerPrimitiveHitMaker & operator = (TriggerPrimitiveHitMaker const &) = delete;
  TriggerPrimitiveHitMaker & operator = (TriggerPrimitiveHitMaker &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;
 
  // Selected optional functions.
  void beginJob() override;

  void endJob() override;

private:

  // Declare member data here.

int fPrimMode;
std::string fAllHitsInstanceName;
void reconfigure(fhicl::ParameterSet const& p) ;

/*
TH2D* fInthist;
TH2D* fMaxhist;
TH2D* fTothist;
TH1D* fInthist_u;
TH1D* fInthist_v;
TH1D* fInthist_y;
TH1D* fMaxhist_u;
TH1D* fMaxhist_v;
TH1D* fMaxhist_y;
TH1D* fTothist_u;
TH1D* fTothist_v;
TH1D* fTothist_y;


TH2D* fIntTothist_u;
TH2D* fIntTothist_v;
TH2D* fIntTothist_y;
TH2D* fIntMaxhist_u;
TH2D* fIntMaxhist_v;
TH2D* fIntMaxhist_y;
TH2D* fMaxTothist_u;
TH2D* fMaxTothist_v;
TH2D* fMaxTothist_y;
*/

};


TriggerPrimitiveHitMaker::TriggerPrimitiveHitMaker(fhicl::ParameterSet const & pset)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
	this->reconfigure(pset);
	recob::HitCollectionCreator::declare_products(*this,fAllHitsInstanceName);

	return;
}


void TriggerPrimitiveHitMaker::reconfigure(fhicl::ParameterSet const& p){
	fPrimMode = p.get<int>("PrimMode");
	fAllHitsInstanceName = "";

	return;

}



void TriggerPrimitiveHitMaker::beginJob()
{

art::ServiceHandle<art::TFileService const> tfs;

 /*
 fInthist               = tfs->make<TH2D>("IntegralAll", "Integral All;Channel;Integral", 8256, 0, 8256, 100, 0, 1000);
 fTothist               = tfs->make<TH2D>("TimeOverThresholdAll", "Time Over Threshold All;Channel;Tot",8256, 0, 8256, 100, 0, 100);
 fMaxhist               = tfs->make<TH2D>("MaxValueAll", "Max Value All;Channel;Max Value",8256, 0, 8256, 100, 0, 200);
 fInthist_u             = tfs->make<TH1D>("IntegralU","Integral U Plane;Integral;Count", 100,0,1000);
 fInthist_v             = tfs->make<TH1D>("IntegralV","Integral V Plane;Integral;Count", 100,0,1000);
 fInthist_y             = tfs->make<TH1D>("IntegralY","Integral Y Plane;Integral;Count", 100,0,1000);
 fTothist_u             = tfs->make<TH1D>("TimeOverThresholdU","Time Over Threshold U Plane;Tot;Count", 100,0,100);
 fTothist_v             = tfs->make<TH1D>("TimeOverThresholdV","Time Over Threshold V Plane;Tot;Count", 100,0,100);
 fTothist_y             = tfs->make<TH1D>("TimeOverThresholdY","Time Over Threshold Y Plane;Tot;Count", 100,0,100);
 fMaxhist_u             = tfs->make<TH1D>("MaxValueU","Max Peak U Plane;Peak Value;Count", 100,0,200);
 fMaxhist_v             = tfs->make<TH1D>("MaxValueV","Max Peak V Plane;Peak Value;Count", 100,0,200);
 fMaxhist_y             = tfs->make<TH1D>("MaxValueY","Max Peak Y Plane;Peak Value;Count", 100,0,200);
 fIntTothist_u          = tfs->make<TH2D>("IntTotU", "Tot vs. Integral U Plane;Tot;Integral", 100, 0, 100, 100, 0, 1000); 
 fIntTothist_v          = tfs->make<TH2D>("IntTotV", "Tot vs. Integral V Plane;Tot;Integral", 100, 0, 100, 100, 0, 1000);
 fIntTothist_y          = tfs->make<TH2D>("IntTotY", "Tot vs. Integral Y Plane;Tot;Integral", 100, 0, 100, 100, 0, 1000);
 fIntMaxhist_u          = tfs->make<TH2D>("IntMaxU", "Max vs. Integral U Plane;Max;Integral", 100, 0, 200, 100, 0, 1000);
 fIntMaxhist_v          = tfs->make<TH2D>("IntMaxV", "Max vs. Integral V Plane;Max;Integral", 100, 0, 200, 100, 0, 1000);
 fIntMaxhist_y          = tfs->make<TH2D>("IntMaxY", "Max vs. Integral Y Plane;Max;Integral", 100, 0, 200, 100, 0, 1000);
 fMaxTothist_u          = tfs->make<TH2D>("MaxTotU", "Tot vs. Max U Plane;Tot;Max", 100, 0, 100, 100, 0, 200);
 fMaxTothist_v          = tfs->make<TH2D>("MaxTotV", "Tot vs. Max V Plane;Tot;Max", 100, 0, 100, 100, 0, 200);
 fMaxTothist_y          = tfs->make<TH2D>("MaxTotY", "Tot vs. Max Y Plane;Tot;Max", 100, 0, 100, 100, 0, 200);
 */

 }


void TriggerPrimitiveHitMaker::produce(art::Event & e)
{
 // Implementation of required member function here.
    recob::HitCollectionCreator allHitCol(*this, e, fAllHitsInstanceName);
    //recob::HitCollectionCreator* filteredHitCol = 0;  
    int numHits = 0;

    art::Handle<std::vector<recob::Wire>> wiredata;
    e.getByLabel("snnodeco",wiredata);

    for(size_t rdIter = 0; rdIter < wiredata->size(); ++rdIter){
	// get the reference to the current non-deconvolved recob::Wire                                                                                               
	art::Ptr<recob::Wire> wireVec(wiredata, rdIter);
	//art::Ptr<raw::RawDigit> rawdigits = RawDigits.at(rdIter);
	auto channel = wireVec->Channel();
	auto zsROIs = wireVec->SignalROI();
	art::ServiceHandle<geo::Geometry const> geom;
	std::vector<geo::WireID> wids = geom->ChannelToWire(channel);
	geo::WireID wid  = wids[0];
	//geo::PlaneID::PlaneID_t plane = wid.Plane;

	std::vector<recob::Hit>  filteredHitVec;

	for (auto iROI = zsROIs.begin_range(); iROI != zsROIs.end_range(); ++iROI) {
	    auto ROI = *iROI;
	    const size_t firstTick = ROI.begin_index();
	    const size_t endTick = ROI.end_index();

	    //std::string titlestring = "roi_original"+ std::to_string(channel);
	    //TH1D horig(titlestring.c_str(), "roi_original;Tick;ADC", endTick + 1 - firstTick, firstTick, endTick + 1);
	    //horig.SetLineColor(kBlack);

	    float integralsum = 0; 
	    float maxpeak = 0;
	    float  tot = 0;
	    float peaktime = 0;
 

	    for (size_t iTick = ROI.begin_index(); iTick < ROI.end_index(); iTick++ ){
		//if (channel > 5500 && channel < 5505){
		//std::string titlestring = "roi_original" + std::to_string(channel);
		//TH1D horig(titlestring.c_str(), "roi_original;Tick:ADC", endTick+1 - firstTick, firstTick, endTick + 1);
		//horig.SetLineColor(kBlack);
		//horig.Fill((int)iTick,ROI[iTick]);
		//TCanvas c ("c", "c");
		//horig.Draw("hist ][");
		//c.Modified();
		//c.Update();
		//c.Print(".png");}//This Line
		integralsum +=  std::abs(ROI[iTick]);
		if (std::abs(ROI[iTick]) > maxpeak){
		    maxpeak = std::abs(ROI[iTick]);
		    peaktime = iTick;
		}
	    }


	    //TCanvas c ("c", "c");
	    //horig.Draw("hist ][");
	    //c.Modified();
	    //c.Update();
	    //c.Print(".png");



	    tot = endTick - firstTick;
	    float widthtot = tot/2;

	    float ucut = -1;
	    float vcut = -1;
	    float ycut = -1;
	    float primitive = -1;

	    if (fPrimMode == 0){
		primitive = integralsum;
		ucut = 200;
		vcut = 120;
		ycut = 200;
	    } else if (fPrimMode == 1 ) {
		primitive = maxpeak;
		ucut = 20;
		vcut = 13;
		ycut = 25;
	    } else if (fPrimMode == 2 ){
		primitive = tot;
		ucut = 21;
		vcut = 20;
		ycut = 24;
	    }

	    int pass  = 0;   //If the waveform passes the cut, pass = 1, if not pass = 0
	    if (channel <= 2399 && primitive >= ucut){
		pass = 1;
	    } else if (channel > 2399 && channel <= 4799 && primitive >= vcut){
		pass = 1; 
	    } else if (channel >= 4800 && primitive >= ycut){
		pass = 1;
	    }

	    if (pass == 1){
		recob::HitCreator hitcreator(
					(recob::Wire)*wireVec,
					(geo::WireID)wid,
					(float)widthtot,
					(float)peaktime, 
					(float)1,
					(float)primitive,
					(float)1,
					(float)primitive,
					(float)1.0,
					(float)primitive,
					(short int)1,
					(short int)0,
					(float)1.0,
					(int)1,
					(size_t)firstTick
					);
		filteredHitVec.push_back(hitcreator.copy());
	    
		art::Ptr<raw::RawDigit> dummyRawDigits;

		const recob::Hit hit(hitcreator.move()); 
		allHitCol.emplace_back(std::move(hit), wireVec, dummyRawDigits);

		numHits++;

		//for(const auto& filteredHit : filteredHitVec)
		//    filteredHitCol->emplace_back(filteredHit, wireVec, dummyRawDigits);
	    }
            

	    /* 
	    if (channel <=  2399){
		fInthist_u->Fill(integralsum);
		fMaxhist_u->Fill(maxpeak);
		fTothist_u->Fill(tot);
		fIntTothist_u->Fill(tot,integralsum);
		fIntMaxhist_u->Fill(maxpeak,integralsum);
		fMaxTothist_u->Fill(tot,maxpeak);
	    }
	    if (channel > 2399 && channel <= 4799){
		fInthist_v->Fill(integralsum);
		fMaxhist_v->Fill(maxpeak);
		fTothist_v->Fill(tot);
		fIntTothist_v->Fill(tot,integralsum);  
		fIntMaxhist_v->Fill(maxpeak,integralsum);
		fMaxTothist_v->Fill(tot,maxpeak);
	    }
	    if (channel >= 4800){
		fInthist_y->Fill(integralsum);
		fMaxhist_y->Fill(maxpeak);
		fTothist_y->Fill(tot);
		fIntTothist_y->Fill(tot,integralsum);
		fIntMaxhist_y->Fill(maxpeak,integralsum);
		fMaxTothist_y->Fill(tot,maxpeak);
	    }

	   fInthist->Fill(channel,integralsum); 

           fMaxhist->Fill(channel,maxpeak);

           fTothist->Fill(channel,tot);
	   */
	} // end of roi loop
    
    } // end of channel loop
         
/*
	 fIntTothist_u->Fit("pol1");
         fIntMaxhist_u->Fit("pol1");
         fMaxTothist_u->Fit("pol1");
         fIntTothist_v->Fit("pol1");
         fIntMaxhist_v->Fit("pol1");
         fMaxTothist_v->Fit("pol1");
         fIntTothist_y->Fit("pol1");
         fIntMaxhist_y->Fit("pol1");
         fMaxTothist_y->Fit("pol1");
*/

    allHitCol.put_into(e);
}


void TriggerPrimitiveHitMaker::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(TriggerPrimitiveHitMaker)
