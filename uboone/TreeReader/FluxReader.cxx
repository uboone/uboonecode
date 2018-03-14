///////////////////////////////////////////////////////////////////////
/// \file  FluxReader.cxx
/// \brief Source to read beam flux files
/// \author  zarko@fnal.gov
////////////////////////////////////////////////////////////////////////

//LArSoft 
#include "larcoreobj/SummaryData/POTSummary.h"

//ART, ...
#include "art/Framework/IO/Sources/put_product_in_principal.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

//root
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"

#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "uboone/FluxReader/FluxReader.h"
#include "GSimpleInterface.h"
#include "NTupleInterface.h"
#include "cetlib/exception.h"

namespace fluxr {

  FluxReader::FluxReader(fhicl::ParameterSet const & ps, 
      art::ProductRegistryHelper &helper,
      art::SourceHelper const &pm)
    :
      fSourceHelper(pm),
      fCurrentSubRunID()
  {
    helper.reconstitutes<std::vector<simb::MCFlux>, art::InEvent>("flux");
    helper.reconstitutes<std::vector<simb::MCTruth>, art::InEvent>("flux");
    helper.reconstitutes<sumdata::POTSummary, art::InSubRun >("flux");

    art::ServiceHandle<art::TFileService> tfs;
    //initialize beam histograms specified in fhicl file
    art::TFileDirectory tffluxdir = tfs->mkdir( "Flux" ); 
    string nutype[]    ={    "nue",        "nuebar",      "numu",         "numubar"};
    string nultx[]     ={"#nu_{e}", "#bar{#nu}_{e}", "#nu_{#mu}", "#bar{#nu}_{#mu}"};
    string pltx[]      ={"#mu^{#pm}","#pi^{#pm}","K^{0}_{L}","K^{#pm}"};
    string secltx[]  ={"pBe->#pi^{#pm}->...->#mu^{#pm}",
      "pBe->#pi^{#pm}->..(not #mu^{#pm})..",
      "pBe->K^{0}_{L}->...",
      "pBe->K^{#pm}->...",
      "pBe->(p or n)->..."};


    int nbins=ps.get<int>("nBins",0);
    int Elow=ps.get<float>("Elow",0);
    int Ehigh=ps.get<float>("Ehigh",0);

    for (int i=0;i<4;i++) {
      fHFlux[i]=tffluxdir.make<TH1D>(Form("h50%i",i+1),
          Form("%s (all);Energy %s (GeV);#phi(%s)/50MeV/POT",nultx[i].c_str(),nultx[i].c_str(),nultx[i].c_str()),
          nbins,Elow,Ehigh);
      fHFlux[i]->Sumw2();
    }

    for (int inu=0;inu<4;inu++) {
      for (int ipar=0;ipar<4;ipar++) {
        fHFluxParent[inu][ipar]=tffluxdir.make<TH1D>(Form("h5%i%i",ipar+1,inu+1),
            Form("...->%s->%s;Energy %s (GeV);#phi(%s)/50MeV/POT",pltx[ipar].c_str(),nultx[inu].c_str(),nultx[inu].c_str(),nultx[inu].c_str()),
            nbins,Elow,Ehigh);

        fHFluxParent[inu][ipar]->Sumw2();
      }
      for (int isec=0;isec<5;isec++) {
        fHFluxSec[inu][isec]=tffluxdir.make<TH1D>(Form("h7%i%i",isec+1,inu+1),
            Form("%s->%s;Energy %s (GeV);#phi(%s)/50MeV/POT",secltx[isec].c_str(),nultx[inu].c_str(),nultx[inu].c_str(),nultx[inu].c_str()),
            nbins,Elow,Ehigh);
        fHFluxSec[inu][isec]->Sumw2();
      }
    }


    //  TH1D* h=tffluxdir.make<TH1D>("numu","numu",240,0,120);
    fEventCounter=0; 
    fEntry=ps.get<uint32_t>("skipEvents", 0);
    fMaxEvents=ps.get<int>("maxEvents", -1);
    fInputType=ps.get<std::string>("inputType");
  }

  void FluxReader::closeCurrentFile()  
  {    
    mf::LogInfo(__FUNCTION__)<<"File boundary (processed "<<fEventCounter<<" events)"<<std::endl;
    fCurrentSubRunID.flushSubRun();
    fEventCounter=0;
    fFluxInputFile->Close();
    delete fFluxInputFile;
    fSkipEvents=0; //if crossing file boundary don't skip events in next file
    fEntry=0;
  }

  void FluxReader::readFile(std::string const &name,
      art::FileBlock* &fb)
  {
    // Fill and return a new Fileblock.
    fb = new art::FileBlock(art::FileFormatVersion(1, "FluxReader"), name);

    fFluxInputFile=new TFile(name.c_str());
    if (fFluxInputFile->IsZombie()) {
      throw cet::exception(__PRETTY_FUNCTION__) << "Failed to open "<<fFluxInputFile<<std::endl;
    } 

    if (fInputType=="gsimple") {
      fFluxDriver=new GSimpleInterface();
      ((GSimpleInterface*)fFluxDriver)->SetRootFile(fFluxInputFile);
    } 
    else if (fInputType=="ntuple"){
      fFluxDriver=new NTupleInterface();
      ((NTupleInterface*)fFluxDriver)->SetRootFile(fFluxInputFile);
    }
    else {
      throw cet::exception(__PRETTY_FUNCTION__) << "Ntuple format "<<fInputType<<" not supported"<<std::endl;
    }

    std::cout<<"File has "<<fFluxDriver->GetEntries()<<" entries"<<std::endl;
    std::cout<<"POT = "<<fFluxDriver->GetPOT()<<std::endl;
    std::cout<<"Run = "<<fFluxDriver->GetRun()<<std::endl;
    fPOT+=fFluxDriver->GetPOT();
  }


  bool FluxReader::readNext(art::RunPrincipal* const &/*inR*/,
      art::SubRunPrincipal* const &/*inSR*/,
      art::RunPrincipal* &outR,
      art::SubRunPrincipal* &outSR,
      art::EventPrincipal* &outE)
  {
    if (fMaxEvents > 0 && fEventCounter == unsigned(fMaxEvents))
      return false;

    if (fEventCounter%10000==0)
      mf::LogInfo(__FUNCTION__)<<"Attempting to read event: "<<fEventCounter<<std::endl;
    // Create empty result, then fill it from current file:
    std::unique_ptr< std::vector<simb::MCFlux>  > mcfluxvec(new std::vector<simb::MCFlux >);
    std::unique_ptr< std::vector<simb::MCTruth>  > mctruthvec(new std::vector<simb::MCTruth >);

    simb::MCFlux flux;
    if (!fFluxDriver->FillMCFlux(fEntry,flux))
      return false;

    //fake mctruth product to cheat eventweight that gets neutrino energy from it
    simb::MCTruth mctruth;
    simb::MCParticle mcpnu(0,flux.fntype,"Flux");
    mcpnu.AddTrajectoryPoint(fFluxDriver->GetNuPosition(), fFluxDriver->GetNuMomentum());
    mctruth.Add(mcpnu);
    mctruth.SetNeutrino(0,0,0,0,0,0,0,0,0,0);
    mctruthvec->push_back(mctruth);

    std::cout << "fNType: " << flux.fntype << " fPType: " << flux.fptype << std::endl;

    int ipdg=0;
    //select index matching order in nutype defined in constructor
    if (flux.fntype==12)       ipdg=0;
    else if (flux.fntype==-12) ipdg=1;
    else if (flux.fntype==14)  ipdg=2;
    else if (flux.fntype==-14) ipdg=3;

    double enu=flux.fnenergyn;
    double totwgh=flux.fnwtnear*flux.fnimpwt;//for gsimple files weight is actually 1 
    //    std::cout<<"tot wgh "<<totwgh<<"\t"<<flux.fnwtnear<<"\t"<<flux.fnimpwt<<std::endl;
    totwgh=1;
    fHFlux[ipdg]->Fill(enu,totwgh);

    std::cout << "filling histogram " << ipdg << " bin " << enu << " with " << totwgh << std::endl; 

    if (flux.fptype==13 || flux.fptype==-13) //mu+-
      fHFluxParent[ipdg][0]->Fill(enu,totwgh);
    else if (flux.fptype==211 || flux.fptype==-211) //pi+-
      fHFluxParent[ipdg][1]->Fill(enu,totwgh);
    else if (flux.fptype==130) //K0L
      fHFluxParent[ipdg][2]->Fill(enu,totwgh);
    else if (flux.fptype==321 || flux.fptype==-321) //K+-
      fHFluxParent[ipdg][3]->Fill(enu,totwgh);

    if (fabs(flux.ftptype==211) && fabs(flux.fptype)==13)
      fHFluxSec[ipdg][0]->Fill(enu,totwgh);
    else if (fabs(flux.ftptype)==211)
      fHFluxSec[ipdg][1]->Fill(enu,totwgh);
    else if (fabs(flux.ftptype)==130)
      fHFluxSec[ipdg][2]->Fill(enu,totwgh);
    else if (fabs(flux.ftptype)==321)
      fHFluxSec[ipdg][3]->Fill(enu,totwgh);
    else if (flux.ftptype==2212 || flux.ftptype==2112)
      fHFluxSec[ipdg][4]->Fill(enu,totwgh);

    mcfluxvec->push_back(flux);
    fEventCounter++;
    fEntry++;

    art::RunNumber_t rn = fFluxDriver->GetRun();
    if (rn==0) rn=999999;
    art::Timestamp tstamp(time(0));

    art::SubRunID newID(rn, 0); //subrun not used in flux files, so set to 0
    if (fCurrentSubRunID.runID() != newID.runID()) { // New Run
      outR = fSourceHelper.makeRunPrincipal(rn, tstamp);
    }
    if (fCurrentSubRunID != newID) { // New SubRun
      outSR = fSourceHelper.makeSubRunPrincipal(rn,0,tstamp);
      std::unique_ptr<sumdata::POTSummary> pot(new sumdata::POTSummary);    
      pot->totpot = fPOT;
      pot->totgoodpot = fPOT;
      fPOT=0;

      art::put_product_in_principal(std::move(pot),
          *outSR,
          "flux");

      fCurrentSubRunID = newID;        
    }


    outE = fSourceHelper.makeEventPrincipal(fCurrentSubRunID.run(),
        fCurrentSubRunID.subRun(),
        fEventCounter,
        tstamp);

    // Put products in the event.
    art::put_product_in_principal(std::move(mcfluxvec),
        *outE,
        "flux"); // Module label
    art::put_product_in_principal(std::move(mctruthvec),
        *outE,
        "flux"); // Module label
    return true;
  }

  void FluxReader::endJob()
  {
    std::cout<<"End job function"<<std::endl;
  }
}
