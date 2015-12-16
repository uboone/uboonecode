#include "WeightCalcCreator.h"
#include "WeightCalc.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "NuReweight/art/NuReweight.h" //GENIEReweight.h"

#include "SimulationBase/MCFlux.h"
#include "SimulationBase/MCTruth.h"
#include "SimulationBase/GTruth.h"

namespace evwgh {
  class GenieWeightCalc : public WeightCalc
  {
  public:
    GenieWeightCalc();
    void Configure(fhicl::ParameterSet const& pset);
    std::vector<std::vector<double> > GetWeight(art::Event & e);
    
  private:
    // The reweighting utility class:
    std::vector<rwgt::NuReweight *> reweightVector;

    CLHEP::RandGaussQ *fGaussRandom;
    std::string fGenieModuleLabel;

    enum EReweight {kNCEL, kQEMA, kQEVec, kResGanged, kCCRes, kNCRes, 
		    kCoh, kNonResRvp1pi, kNonResRvbarp1pi, kNonResRvp2pi,	
		    kNonResRvbarp2pi, kResDecay, kNC, kDIS, kDISnucl, 	
		    kAGKY, kNReWeights};
    
    DECLARE_WEIGHTCALC(GenieWeightCalc)
  };
  GenieWeightCalc::GenieWeightCalc()
  {
  }

  void GenieWeightCalc::Configure(fhicl::ParameterSet const& p)
  {
    //global config
    fGenieModuleLabel= p.get< std::string > ("genie_module_label");

    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    //calc config
    std::vector<std::string> pars = pset.get< std::vector<std::string> > ("parameter_list");	
    std::vector<float> parsigmas = pset.get< std::vector<float> > ("parameter_sigma");	
    std::string mode             = pset.get<std::string>("mode");

    if (pars.size() != parsigmas.size() )
      throw cet::exception(__FUNCTION__) << GetName()<<"::Bad fcl configuration. parameter_list and parameter_sigma need to have same number of parameters."<<std::endl;

    int number_of_multisims = pset.get< int > ("number_of_multisims");
      
    std::vector<EReweight> erwgh;
    for( auto & s : pars){
      if (s == "NCEL") erwgh.push_back(kNCEL);
      else if (s == "QEMA") erwgh.push_back(kQEMA);
      else if (s == "QEVec") erwgh.push_back(kQEVec);
      else if (s == "ResGanged") erwgh.push_back(kResGanged);
      else if (s == "CCRes") erwgh.push_back(kCCRes);
      else if (s == "NCRes") erwgh.push_back(kNCRes);
      else if (s == "Coh") erwgh.push_back(kCoh);
      else if (s == "NonResRvp1pi") erwgh.push_back(kNonResRvp1pi);
      else if (s == "NonResRvbarp1pi") erwgh.push_back(kNonResRvbarp1pi);
      else if (s == "NonResRvp2pi") erwgh.push_back(kNonResRvp2pi);
      else if (s == "NonResRvbarp2pi") erwgh.push_back(kNonResRvbarp2pi);
      else if (s == "ResDecay") erwgh.push_back(kResDecay);
      else if (s == "NC") erwgh.push_back(kNC);
      else if (s == "DIS") erwgh.push_back(kDIS);
      else if (s == "DISnucl") erwgh.push_back(kDISnucl);
      else if (s == "AGKY") erwgh.push_back(kAGKY);
      else {
	throw cet::exception(__FUNCTION__) << GetName()<<"::Physical process "<<s<<" you requested is not available to reweight." << std::endl;
      }
    }
      
    //Prepare sigmas
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    fGaussRandom = new CLHEP::RandGaussQ(rng->getEngine(GetName()));

    std::vector<std::vector<float> > reweightingSigmas(erwgh.size());
    for (unsigned int i = 0; i < reweightingSigmas.size(); ++i) {
      reweightingSigmas[i].resize(number_of_multisims);
      for (int j = 0; j < number_of_multisims; j ++) {
	if (mode.find("multisim") != std::string::npos )
	  reweightingSigmas[i][j] = parsigmas[i]*fGaussRandom->shoot(&rng->getEngine(GetName()),0.,1.);
	else
	  reweightingSigmas[i][j] = parsigmas[i];
      }
    }

    reweightVector.resize(number_of_multisims);
    
    for (int weight_point = 0; 
	 weight_point < number_of_multisims;
	 weight_point++){
      
      reweightVector[weight_point] = new rwgt::NuReweight;
      
      for (unsigned int i_reweightingKnob=0;i_reweightingKnob<erwgh.size();i_reweightingKnob++) {
	std::cout<<GetName()<<"::Setting up rwgh "<<weight_point<<"\t"<<i_reweightingKnob<<"\t"<<erwgh[i_reweightingKnob]<<std::endl; 
	switch (erwgh[i_reweightingKnob]){
	case kNCEL:
	  // reweightVector[weight_point]
	  //   -> ReweightNCEL(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kQEMA:
	  reweightVector[weight_point]
	    -> ReweightQEMA(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kQEVec:
	  reweightVector[weight_point]
	    -> ReweightQEVec(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kResGanged:
	  reweightVector[weight_point]
	    -> ReweightResGanged(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kCCRes:
	  reweightVector[weight_point]
	    -> ReweightCCRes(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNCRes:
	  reweightVector[weight_point]
	    -> ReweightNCRes(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kCoh:
	  // reweightVector[weight_point]
	  //   -> ReweightCoh(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNonResRvp1pi:
	  reweightVector[weight_point]
	    -> ReweightNonResRvp1pi(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNonResRvbarp1pi:
	  reweightVector[weight_point]
	    -> ReweightNonResRvbarp1pi(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNonResRvp2pi:
	  reweightVector[weight_point]
	    -> ReweightNonResRvp2pi(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNonResRvbarp2pi:
	  reweightVector[weight_point]
	    -> ReweightNonResRvbarp2pi(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kResDecay:
	  // reweightVector[weight_point]
	  //   -> ReweightResDecay(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNC:
	  reweightVector[weight_point]
	    -> ReweightNC(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kDIS:
	  // reweightVector[weight_point]
	  //   -> ReweightDIS(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kDISnucl:
	  reweightVector[weight_point]
	    -> ReweightDISnucl(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kAGKY:
	  // reweightVector[weight_point]
	  //   -> ReweightAGKY(reweightingSigmas[i_reweightingKnob][weight_point]);
	  break;
	case kNReWeights:
	  break;
	}
      }
   
    } //loop over nWeights
    // Tell all of the reweight drivers to configure themselves:
    std::cout<< GetName()<<"::Setting up "<<reweightVector.size()<<" reweightcalcs"<<std::endl;
    for(auto & driver : reweightVector){
      driver -> Configure();
    }
  }

  std::vector<std::vector<double> > GenieWeightCalc::GetWeight(art::Event & e)
  { 
    //returns a vector of weights for each neutrino interaction in the event

    //get the MC generator information out of the event       
    //these are all handles to mc information.
    art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;  
    art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
    art::Handle< std::vector<simb::GTruth> > gTruthHandle;

    //actually go and get the stuff
    e.getByLabel(fGenieModuleLabel,mcTruthHandle);
    e.getByLabel(fGenieModuleLabel,mcFluxHandle);
    e.getByLabel(fGenieModuleLabel,gTruthHandle);

    std::vector<art::Ptr<simb::MCTruth> > mclist;
    art::fill_ptr_vector(mclist, mcTruthHandle);

    std::vector<art::Ptr<simb::GTruth > > glist;
    art::fill_ptr_vector(glist, gTruthHandle);

    //calculate weight(s) here 
    std::vector<std::vector<double> >weight(mclist.size());
    for ( unsigned int inu=0; inu<mclist.size();inu++) {
      weight[inu].resize(reweightVector.size());    
      for (unsigned int i_weight = 0; 
	   i_weight < reweightVector.size(); 
	   i_weight ++){
	weight[inu][i_weight]= reweightVector[i_weight]-> CalcWeight(*mclist[inu],*glist[inu]);
      }
    }
    return weight;

  }
  REGISTER_WEIGHTCALC(GenieWeightCalc)
}
