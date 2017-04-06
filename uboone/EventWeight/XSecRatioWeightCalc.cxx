#include "WeightCalcCreator.h"
#include "WeightCalc.h"
#include "uboone/EventWeight/IFDHFileTransfer.h"

#include <iostream>

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
//#include "artextensions/SeedService/SeedService.hh"
#include "nutools/RandomUtils/NuRandomService.h"

#include "CLHEP/Random/RandGaussQ.h"

#include "nusimdata/SimulationBase/MCFlux.h"
//#include "SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"
//#include "SimulationBase/MCTruth.h"
#include "TH2D.h"
#include "TFile.h"
#include "TH1F.h"

namespace evwgh {
  class XSecRatioWeightCalc : public WeightCalc
  {  
     evwgh::IFDHFileTransfer IFDH;
     public:
       XSecRatioWeightCalc();
       void Configure(fhicl::ParameterSet const& p);

       std::array<std::array<double, 101>, 1001> reweightingRatio; 
       std::array<std::array<double, 101>, 1001> reweightingSigmas;
       std::array<std::array<double, 101>, 1001> reweightingSigmas_e;
 
       std::array<double, 101> enumuon_nominal;
       std::array<double, 1001> q2muon_nominal;
       std::array<double, 101> binwidth_enu;
       std::array<double, 1001> binwidth_q2;
       std::array<std::array<double, 101>, 1001> xsecratio_nominal;
       std::array<std::array<double, 101>, 1001>  xsecratio_Fv3;
       std::array<std::array<double, 101>, 1001>  xsecratio_Fa3;
       std::array<std::array<double, 101>, 1001>  xsecratio_Fv3Fa3;

       std::array<std::array<double, 101>, 1001> xsecmuon_nominal;
       std::array<std::array<double, 101>, 1001>  xsecmuon_Fv3;
       std::array<std::array<double, 101>, 1001>  xsecmuon_Fa3;
       std::array<std::array<double, 101>, 1001>  xsecmuon_Fv3Fa3;
       std::array<std::array<double, 101>, 1001> xsece_nominal;
       std::array<std::array<double, 101>, 1001>  xsece_Fv3;
       std::array<std::array<double, 101>, 1001>  xsece_Fa3;
       std::array<std::array<double, 101>, 1001>  xsece_Fv3Fa3;
  
       enum XSecRatioRW{
         knominal,
         kFv3,
         kFa3,
         kFv3Fa3
       };
       std::vector<XSecRatioRW> xsecratiorwgh;

 
       void GetXSecRatio();

       std::vector<std::vector<double > > GetWeight(art::Event & e);
     private:

       CLHEP::RandGaussQ *fGaussRandom; 
       //The reweighting utility class
       //std::vector<rwgt::NuReweight *>  reweightVector;
       //std::vector<rwgt::NuReweight *> reweightVector_enu;
       //std::vector<rwgt::NuReweight *> reweightVector_q2;
     
       std::vector<double> fWeightArray;
       int fNmultisims;
       std::string fMode;
       std::string fGenieModuleLabel;

/*       enum XSecRatioRW{
         knominal,
         kFv3,
         kFa3,
         kFv3Fa3
       };
*/       
     DECLARE_WEIGHTCALC(XSecRatioWeightCalc)
  };
  XSecRatioWeightCalc::XSecRatioWeightCalc()
  {
  }

  void XSecRatioWeightCalc::Configure(fhicl::ParameterSet const& p)
  {
    //global config
    fGenieModuleLabel= p.get< std::string > ("genie_module_label");
        
    fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());
    //calc config
    fNmultisims=pset.get<int>("number_of_multisims");
    fMode      =pset.get<std::string>("mode");


    //GetXSecRatio();

    //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~~^~^~^~^~^~^
    //Prepare random generator
    art::ServiceHandle<art::RandomNumberGenerator> rng;
    fGaussRandom = new CLHEP::RandGaussQ(rng->getEngine(GetName()));


    //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^
    //calc config
    std::cout<<"tryingn to get the fcl parameters:   "<<std::endl;
    std::vector<std::string> pars = pset.get< std::vector<std::string> > ("parameter_list");
    std::vector<float> parsigmas = pset.get< std::vector<float> > ("parameter_sigma");
    fWeightArray.resize(fNmultisims);


    if(fMode.find("multisim") !=std::string::npos)
       for(int i=0; i<fNmultisims; i++) fWeightArray[i]=fGaussRandom->shoot(&rng->getEngine(GetName()),0,1.);
    else
       for (int i=0;i<fNmultisims;i++) fWeightArray[i]=1.;  

    
    //-------------------------------------------------------- 
//    std::vector<XSecRatioRW> xsecratiorwgh;

    if(pars.size() != parsigmas.size() )
       throw cet::exception(__FUNCTION__) << GetName()<<"::Bad fcl configuration. parameter_list and parameter_sigma need to have same number of parameters."<<std::endl;
    //int number_of_q2bins=pset.get< int > ("number_of_q2bins", 1000);
    //int number_of_ebins=pset.get< int > ("number_of_ebins", 100);

    for(auto & s: pars){
       if (s=="nominal") xsecratiorwgh.push_back(knominal);
       else if(s=="Fv3") xsecratiorwgh.push_back(kFv3);
       else if(s=="Fa3") xsecratiorwgh.push_back(kFa3);
       else if(s=="Fv3Fa3") xsecratiorwgh.push_back(kFv3Fa3);
       else {
       throw cet::exception(__FUNCTION__) << GetName()<<"::Physical process "<<s<<" you requested is not available to reweight." << std::endl;
       }
    }
    //--------------------------------------------------------------------
    for(unsigned int cc=0; cc<xsecratiorwgh.size(); cc++){
       std::cout<<"s= "<<xsecratiorwgh[cc]<<std::endl;
     
    }
    GetXSecRatio();



    //--------------------------------------------------------------------    
  }
  void XSecRatioWeightCalc::GetXSecRatio() 
  {

    //upload the root file here and get the ratios between nue xsec to numuxsec
    TFile *inputf = new TFile("2DhistogramsForJoseph.root");
    TH2D *h_muon_nominal;
    TH2D *h_muon_Fa3;
    TH2D *h_muon_Fv3;
    TH2D *h_muon_Fv3Fa3;
   
    TH2D *h_e_nominal;
    TH2D *h_e_Fa3;
    TH2D *h_e_Fv3;
    TH2D *h_e_Fv3Fa3;

    h_muon_nominal=(TH2D*)inputf->Get("muonNominal");
    h_muon_Fa3=(TH2D*)inputf->Get("muonFa3");
    h_muon_Fv3=(TH2D*)inputf->Get("muonFv3");
    h_muon_Fv3Fa3=(TH2D*)inputf->Get("muonFv3Fa3");

    h_e_nominal=(TH2D*)inputf->Get("eNominal");
    h_e_Fa3=(TH2D*)inputf->Get("eFa3");
    h_e_Fv3=(TH2D*)inputf->Get("eFv3");
    h_e_Fv3Fa3=(TH2D*)inputf->Get("eFv3Fa3");
 
    int NbinsX_muon_nominal=h_muon_nominal->GetNbinsX();
    int NbinsX_muon_Fa3=h_muon_Fa3->GetNbinsX();
    int NbinsX_muon_Fv3=h_muon_Fv3->GetNbinsX();
    int NbinsX_muon_Fv3Fa3=h_muon_Fv3Fa3->GetNbinsX();

    int NbinsY_muon_nominal=h_muon_nominal->GetNbinsY();
    int NbinsY_muon_Fa3=h_muon_Fa3->GetNbinsY();
    int NbinsY_muon_Fv3=h_muon_Fv3->GetNbinsY();
    int NbinsY_muon_Fv3Fa3=h_muon_Fv3Fa3->GetNbinsY();

    int NbinsX_e_nominal=h_e_nominal->GetNbinsX();
    int NbinsX_e_Fa3=h_e_Fa3->GetNbinsX();
    int NbinsX_e_Fv3=h_e_Fv3->GetNbinsX();
    int NbinsX_e_Fv3Fa3=h_e_Fv3Fa3->GetNbinsX();

    int NbinsY_e_nominal=h_e_nominal->GetNbinsY();
    int NbinsY_e_Fa3=h_e_Fa3->GetNbinsY();
    int NbinsY_e_Fv3=h_e_Fv3->GetNbinsY();
    int NbinsY_e_Fv3Fa3=h_e_Fv3Fa3->GetNbinsY();
 


    std::cout<<"NbinsX_muon_nominal "<<NbinsX_muon_nominal<<std::endl;
    std::cout<<"NbinsY_muon_nominal "<<NbinsY_muon_nominal<<std::endl;

    int number_of_q2bins=NbinsY_muon_nominal;
    int number_of_ebins=NbinsX_muon_nominal;

    
    std::cout<<number_of_q2bins<<std::endl;
    std::cout<<number_of_ebins<<std::endl;
    std::cout<<NbinsX_muon_nominal<<std::endl;
    std::cout<<NbinsX_muon_Fv3<<std::endl;
    std::cout<<NbinsX_muon_Fa3<<std::endl;
    std::cout<<NbinsX_muon_Fv3Fa3<<std::endl;

    std::cout<<NbinsY_muon_nominal<<std::endl;
    std::cout<<NbinsY_muon_Fv3<<std::endl;
    std::cout<<NbinsY_muon_Fa3<<std::endl;
    std::cout<<NbinsY_muon_Fv3Fa3<<std::endl;

    std::cout<<NbinsX_e_nominal<<std::endl;
    std::cout<<NbinsX_e_Fv3<<std::endl;
    std::cout<<NbinsX_e_Fa3<<std::endl;
    std::cout<<NbinsX_e_Fv3Fa3<<std::endl;

    std::cout<<NbinsY_e_nominal<<std::endl;
    std::cout<<NbinsY_e_Fv3<<std::endl;
    std::cout<<NbinsY_e_Fa3<<std::endl;
    std::cout<<NbinsY_e_Fv3Fa3<<std::endl;
    

/*    
*/

/*    
*/

    
    //get the bin content for each histogram
    for(int ind_enu=0; ind_enu<=100; ind_enu++){
      enumuon_nominal[ind_enu]=h_muon_nominal->GetXaxis()->GetBinCenter(ind_enu);
      binwidth_enu[ind_enu]=h_muon_nominal->GetXaxis()->GetBinWidth(ind_enu);
/*      enue_nominal[ind_enu]=h_e_nominal->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_Fv3[ind_enu]=h_muon_Fv3->GetXaxis()->GetBinCenter(ind_enu);
      enue_Fv3[ind_enu]=h_e_Fv3->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_Fa3[ind_enu]=h_muon_Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enue_Fa3[ind_enu]=h_e_Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_Fv3Fa3[ind_enu]=h_muon_Fv3Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enue_Fv3Fa3[ind_enu]=h_e_Fv3Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_nominal[ind_enu]=h_muon_nominal->GetXaxis()->GetBinCenter(ind_enu);
      enue_nominal[ind_enu]=h_e_nominal->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_Fv3[ind_enu]=h_muon_Fv3->GetXaxis()->GetBinCenter(ind_enu);
      enue_Fv3[ind_enu]=h_e_Fv3->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_Fa3[ind_enu]=h_muon_Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enue_Fa3[ind_enu]=h_e_Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_Fv3Fa3[ind_enu]=h_muon_Fv3Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enue_Fv3Fa3[ind_enu]=h_e_Fv3Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_nominal[ind_enu]=h_muon_nominal->GetXaxis()->GetBinCenter(ind_enu);
      enue_nominal[ind_enu]=h_e_nominal->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_Fv3[ind_enu]=h_muon_Fv3->GetXaxis()->GetBinCenter(ind_enu);
      enue_Fv3[ind_enu]=h_e_Fv3->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_Fa3[ind_enu]=h_muon_Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enue_Fa3[ind_enu]=h_e_Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_Fv3Fa3[ind_enu]=h_muon_Fv3Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enue_Fv3Fa3[ind_enu]=h_e_Fv3Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_nominal[ind_enu]=h_muon_nominal->GetXaxis()->GetBinCenter(ind_enu);
      enue_nominal[ind_enu]=h_e_nominal->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_Fv3[ind_enu]=h_muon_Fv3->GetXaxis()->GetBinCenter(ind_enu);
      enue_Fv3[ind_enu]=h_e_Fv3->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_Fa3[ind_enu]=h_muon_Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enue_Fa3[ind_enu]=h_e_Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enumuon_Fv3Fa3[ind_enu]=h_muon_Fv3Fa3->GetXaxis()->GetBinCenter(ind_enu);
      enue_Fv3Fa3[ind_enu]=h_e_Fv3Fa3->GetXaxis()->GetBinCenter(ind_enu);
*/    }
    for(auto ind_q2=0; ind_q2<=1000; ind_q2++){
      q2muon_nominal[ind_q2]=h_muon_nominal->GetYaxis()->GetBinCenter(ind_q2);
      binwidth_q2[ind_q2]=h_muon_nominal->GetYaxis()->GetBinWidth(ind_q2);
/*      q2e_nominal[ind_q2]=h_e_nominal->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_Fv3[ind_q2]=h_muon_Fv3->GetYaxis()->GetBinCenter(ind_q2);
      q2e_Fv3[ind_q2]=h_e_Fv3->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_Fa3[ind_q2]=h_muon_Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2e_Fa3[ind_q2]=h_e_Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_Fv3Fa3[ind_q2]=h_muon_Fv3Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2e_Fv3Fa3[ind_q2]=h_e_Fv3Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_nominal[ind_q2]=h_muon_nominal->GetYaxis()->GetBinCenter(ind_q2);
      q2e_nominal[ind_q2]=h_e_nominal->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_Fv3[ind_q2]=h_muon_Fv3->GetYaxis()->GetBinCenter(ind_q2);
      q2e_Fv3[ind_q2]=h_e_Fv3->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_Fa3[ind_q2]=h_muon_Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2e_Fa3[ind_q2]=h_e_Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_Fv3Fa3[ind_q2]=h_muon_Fv3Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2e_Fv3Fa3[ind_q2]=h_e_Fv3Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_nominal[ind_q2]=h_muon_nominal->GetYaxis()->GetBinCenter(ind_q2);
      q2e_nominal[ind_q2]=h_e_nominal->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_Fv3[ind_q2]=h_muon_Fv3->GetYaxis()->GetBinCenter(ind_q2);
      q2e_Fv3[ind_q2]=h_e_Fv3->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_Fa3[ind_q2]=h_muon_Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2e_Fa3[ind_q2]=h_e_Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_Fv3Fa3[ind_q2]=h_muon_Fv3Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2e_Fv3Fa3[ind_q2]=h_e_Fv3Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_nominal[ind_q2]=h_muon_nominal->GetYaxis()->GetBinCenter(ind_q2);
      q2e_nominal[ind_q2]=h_e_nominal->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_Fv3[ind_q2]=h_muon_Fv3->GetYaxis()->GetBinCenter(ind_q2);
      q2e_Fv3[ind_q2]=h_e_Fv3->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_Fa3[ind_q2]=h_muon_Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2e_Fa3[ind_q2]=h_e_Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2muon_Fv3Fa3[ind_q2]=h_muon_Fv3Fa3->GetYaxis()->GetBinCenter(ind_q2);
      q2e_Fv3Fa3[ind_q2]=h_e_Fv3Fa3->GetYaxis()->GetBinCenter(ind_q2);
*/    }  
     
    for(int ind_alg=0; ind_alg<4; ind_alg++){
     for(int ind_enu_nom=0; ind_enu_nom<=number_of_ebins; ind_enu_nom++){
       for(int ind_q2_nom=0; ind_q2_nom<=number_of_q2bins; ind_q2_nom++){
         if(ind_alg==0){
            xsecmuon_nominal[ind_enu_nom][ind_q2_nom]=h_muon_nominal->GetBinContent(ind_enu_nom,ind_q2_nom);
         }
         else if(ind_alg==1){
            xsecmuon_Fv3[ind_enu_nom][ind_q2_nom]=h_muon_Fv3->GetBinContent(ind_enu_nom,ind_q2_nom); 
         }
         else if(ind_alg==2){
            xsecmuon_Fa3[ind_enu_nom][ind_q2_nom]=h_muon_Fa3->GetBinContent(ind_enu_nom,ind_q2_nom); 
         }
         else if(ind_alg==3){
            xsecmuon_Fv3Fa3[ind_enu_nom][ind_q2_nom]=h_muon_Fv3Fa3->GetBinContent(ind_enu_nom,ind_q2_nom); 
         }       
       }     
     }
    }
    std::cout<<"xsecmuon_nominal[33][13]"<< xsecmuon_nominal[33][13]<<std::endl;
    std::cout<<"xsecmuon_Fv3[33][13]"<< xsecmuon_Fv3[33][13]<<std::endl;

     for(int ind_alg=0; ind_alg<4; ind_alg++){
     for(int ind_enu_nom=0; ind_enu_nom<=number_of_ebins; ind_enu_nom++){
       for(int ind_q2_nom=0; ind_q2_nom<=number_of_q2bins; ind_q2_nom++){
         if(ind_alg==0){
            xsece_nominal[ind_enu_nom][ind_q2_nom]=h_e_nominal->GetBinContent(ind_enu_nom,ind_q2_nom);
         }
         else if(ind_alg==1){
            xsece_Fv3[ind_enu_nom][ind_q2_nom]=h_e_Fv3->GetBinContent(ind_enu_nom,ind_q2_nom); 
         }
         else if(ind_alg==2){
            xsece_Fa3[ind_enu_nom][ind_q2_nom]=h_e_Fa3->GetBinContent(ind_enu_nom,ind_q2_nom); 
         }
         else if(ind_alg==3){
            xsece_Fv3Fa3[ind_enu_nom][ind_q2_nom]=h_e_Fv3Fa3->GetBinContent(ind_enu_nom,ind_q2_nom); 
         }       
       }     
     }
    }





    //
    //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^ 
    //std::vector<XSecRatioRW> xsecratiorwgh;

    for(unsigned int i=0; i<4; ++i){
       for(int j=0; j<= 100; j++){
         for(int k=0; k<= 1000; k++){
            if(i==0){xsecratio_nominal[j][k]=xsecmuon_nominal[j][k]/xsece_nominal[j][k];}
            else if(i==1){xsecratio_Fv3[j][k]=xsecmuon_Fv3[j][k]/xsece_Fv3[j][k];}
            else if(i==2){xsecratio_Fa3[j][k]=xsecmuon_Fa3[j][k]/xsece_Fa3[j][k];}
            else if(i==3){xsecratio_Fv3Fa3[j][k]=xsecmuon_Fv3Fa3[j][k]/xsece_Fv3Fa3[j][k];}
         }
       }
    }
    std::cout<<"xsecratio in the  ....."<<std::endl;
    std::cout<<xsecratio_nominal[70][50]<<" "<<xsecratio_Fv3[70][50]<<std::endl;
    std::cout<<xsecratio_Fa3[70][50]<<" "<<xsecratio_Fv3Fa3[70][50]<<std::endl;

    art::ServiceHandle<art::RandomNumberGenerator> rng;
    fGaussRandom = new CLHEP::RandGaussQ(rng->getEngine(GetName()));


    std::cout<<"start to get tthe reweighting sigmas and reweighting ratios:  "<<std::endl;
    for(unsigned int i_reweightingKnob=0; i_reweightingKnob<xsecratiorwgh.size(); i_reweightingKnob++){   
    //for(unsigned int i_reweightingKnob=2; i_reweightingKnob<3; i_reweightingKnob++){   

   
     //std::cout<<GetName()<<"::Setting up rwgh "<<"\t"<<i_reweightingKnob<<"\t" <<xsecratiorwgh[i_reweightingKnob]<<std::endl;
     switch (xsecratiorwgh[i_reweightingKnob]){
      case knominal:
      break;
      //--------------------------------------------------------------------------------------------
      case kFv3:
      for(int jj=0; jj<=100; jj++){
       for(int kk=0; kk<=1000; kk++){
          //reweightingSigmas[jj][kk]=abs(xsecratio_nominal[jj][kk]-xsecratio_Fv3[jj][kk])
                                   //*xsecratio_Fv3[jj][kk]/xsecratio_nominal[jj][kk];
                                   //*fGaussRandom->shoot(&rng->getEngine(GetName()), 0., 1);        
          reweightingRatio[jj][kk]=xsecratio_Fv3[jj][kk]; 
          reweightingSigmas[jj][kk]=xsecmuon_Fv3[jj][kk]-xsecmuon_nominal[jj][kk];
          reweightingSigmas_e[jj][kk]=xsece_Fv3[jj][kk]-xsece_nominal[jj][kk];
       }
      }
      break;
      //-------------------------------------------------------------------------
      case kFa3:
      for(int jj=0; jj<=100; jj++){
       for(int kk=0; kk<=1000; kk++){
          //reweightingSigmas[jj][kk]=abs(xsecratio_nominal[jj][kk]-xsecratio_Fa3[jj][kk])
                                   //*xsecratio_Fa3[jj][kk]/xsecratio_nominal[jj][kk];
                                   //*fGaussRandom->shoot(&rng->getEngine(GetName()), 0., 1);
          reweightingRatio[jj][kk]=xsecratio_Fa3[jj][kk]; 
          reweightingSigmas[jj][kk]=xsecmuon_Fa3[jj][kk]-xsecmuon_nominal[jj][kk];
          reweightingSigmas_e[jj][kk]=xsece_Fa3[jj][kk]-xsece_nominal[jj][kk];
        }
      }
      break;
      //------------------------------------------------------------------------------
      case kFv3Fa3:
      for(int jj=0; jj<=100; jj++){
       for(int kk=0; kk<=1000; kk++){
          //reweightingSigmas[jj][kk]=abs(xsecratio_nominal[jj][kk]-xsecratio_Fv3Fa3[jj][kk])
                                   //*xsecratio_Fv3Fa3[jj][kk]/xsecratio_nominal[jj][kk];
                                   //*fGaussRandom->shoot(&rng->getEngine(GetName()), 0., 1);
          reweightingRatio[jj][kk]=xsecratio_Fv3Fa3[jj][kk]; 
          reweightingSigmas[jj][kk]=xsecmuon_Fv3Fa3[jj][kk]-xsecmuon_nominal[jj][kk];
          reweightingSigmas_e[jj][kk]=xsece_Fv3Fa3[jj][kk]-xsece_nominal[jj][kk];
        }
      }
      break;
      //----------------------------------------------------------------------------
     } 
    }
    std::cout<<"end of getting reweighting sigmas and reweighting ratios  "<<std::endl;
    //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^
  }  
  std::vector<std::vector<double> > XSecRatioWeightCalc::GetWeight(art::Event & e)
  { 
    //returns a vector of weights for each neutrino interaction in the event

    //get the MC generator information out of the event       
    //these are all handles to mc information.
    art::Handle< std::vector<simb::MCTruth> > mcTruthHandle;  
    art::Handle< std::vector<simb::MCFlux> > mcFluxHandle;
    //art::Handle< std::vector<simb::GTruth> > gTruthHandle;

    //actually go and get the stuff
    e.getByLabel(fGenieModuleLabel,mcTruthHandle);
    e.getByLabel(fGenieModuleLabel,mcFluxHandle);
    //e.getByLabel(fGenieModuleLabel,gTruthHandle);

    std::vector<art::Ptr<simb::MCTruth> > mclist;
    art::fill_ptr_vector(mclist, mcTruthHandle);

    //std::vector<art::Ptr<simb::GTruth > > glist;
    //art::fill_ptr_vector(glist, gTruthHandle);
    //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
    //libo start
    //calculate weight(s) here 
    std::vector<std::vector<double> > weight;
    
    weight.resize(mclist.size());
    //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^
    //get the Q2 of GENIE truth and Enu of GENIE truth
        
    art::Ptr<simb::MCTruth> mctruth;
    double _fnuPDGtruth;
    double _fq2truth;
    double _fEnutruth;

    //get the incoming neutrino's PDG Code
    
 
    for ( unsigned int inu=0; inu<mclist.size();inu++) {
    
      mctruth=mclist[inu];
      weight[inu].resize(fNmultisims);

      _fnuPDGtruth=mctruth->GetNeutrino().Nu().PdgCode();
      _fq2truth=mctruth->GetNeutrino().QSqr();
      _fEnutruth=mctruth->GetNeutrino().Nu().E();

      std::cout<<"Q2= "<<_fq2truth<<std::endl; 
      std::cout<<"Enu= "<<_fEnutruth<<std::endl;

      std::array<double, 100>::iterator iter_enu;
      std::array<double, 1000>::iterator iter_q2;
      for(int ind_multisim=0; ind_multisim<fNmultisims; ind_multisim++){

      int Nbin_enu=int(_fEnutruth*1000/binwidth_enu[1])+1; 
      if(Nbin_enu>100) Nbin_enu=100;

      int Nbin_q2=int(_fq2truth/binwidth_q2[1])+1; 
      if(Nbin_q2>1000) Nbin_q2=1000;


      std::cout<<"Nbin_q2= "<<Nbin_q2<<std::endl;
      std::cout<<"Nbin_enu= "<<Nbin_enu<<std::endl; 
      if (_fnuPDGtruth ==14){

                    //  std::cout<<"test for the xsecratio reweight ..."<<std::endl;
                    //weight[inu][ind_multisim]=1-(1-reweightingRatio[Nbin_enu][Nbin_q2])*fWeightArray[ind_multisim];
                      weight[inu][ind_multisim]=1-reweightingSigmas[Nbin_enu][Nbin_q2]*fWeightArray[ind_multisim];
                    //  std::cout<<"fWeightArray[ind_multisim]=" <<fWeightArray[ind_multisim]<<std::endl;
                    //  std::cout<<"reweightingSigmas[Nbin_enu][Nbin_q2]:"<<reweightingSigmas[Nbin_enu][Nbin_q2]<<std::endl;               
                    //  std::cout<<"ind_multisim = "<<ind_multisim<<std::endl;
                    //  std::cout<<"reweighting factor= "<<weight[inu][ind_multisim]<<std::endl;
       }   //---------------------------------------------------------
      if (_fnuPDGtruth ==12){

                    //  std::cout<<"test for the xsecratio reweight ..."<<std::endl;
                    //weight[inu][ind_multisim]=1-(1-reweightingRatio[Nbin_enu][Nbin_q2])*fWeightArray[ind_multisim];
                      weight[inu][ind_multisim]=1-reweightingSigmas_e[Nbin_enu][Nbin_q2]*fWeightArray[ind_multisim];
                    //  std::cout<<"fWeightArray[ind_multisim]=" <<fWeightArray[ind_multisim]<<std::endl;
                    //  std::cout<<"reweightingSigmas[Nbin_enu][Nbin_q2]:"<<reweightingSigmas[Nbin_enu][Nbin_q2]<<std::endl;               
                    //  std::cout<<"ind_multisim = "<<ind_multisim<<std::endl;
                    //  std::cout<<"reweighting factor= "<<weight[inu][ind_multisim]<<std::endl;
       }   //---------------------------------------------------------
 



       } //end of loop over multisims

       //=========================================================================
       //--------------------------------------------------------
    } //end of loop over inu mctruth
    //libo end~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^^~^~^
    return weight;

  }

/*  std::vector<std::vector<double> > XSecRatioWeightCalc::GetWeight(art::Event & e)
  {
    //calculate weight(s) here 
    std::vector<std::vector<double> > weight;
    return weight;
  }
*/
  REGISTER_WEIGHTCALC(XSecRatioWeightCalc)
}
