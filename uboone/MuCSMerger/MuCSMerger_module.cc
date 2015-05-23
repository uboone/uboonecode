
////////////////////////////////////////////////////////////////////////
// Class:       MuCSMerger
// Module Type: producer
// File:        MuCSMerger_module.cc
//
// Generated at Wed May 20 14:08:15 2015 by Leonidas N. Kalousis using artmod
// from cetpkgsupport v1_08_05.
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
//    trivia : The main routine developed to merge data from the TPC 
//             and the Muon Counter System (MuCS), May 2015
//    author : Odysseas Kanenas
//    e-mail : kalousis@vt.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef MUCSMERGER_H
#define MUCSMERGER_H 

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>
#include <iostream>
#include "RawData/AuxDetDigit.h"
#include "uboone/RawData/uboone_datatypes/MuCSData.h"

using namespace std;
// using namespace gov::fnal::uboone::datatypes;

class MuCSMerger;

class MuCSMerger : public art::EDProducer {
public:
  explicit MuCSMerger( fhicl::ParameterSet const &pset );
  virtual ~MuCSMerger();
  
  // void reconfigure( fhicl::ParameterSet const &pset ) override;
  void produce( art::Event &evt ) override;
      
  private:
  
};

MuCSMerger::MuCSMerger( fhicl::ParameterSet const &pset )
// :
// Initialize member data here.
{
  produces< std::vector<gov::fnal::uboone::datatypes::MuCSData> >();  
  
}

MuCSMerger::~MuCSMerger()
{}

void MuCSMerger::produce( art::Event &evt )
{
  cout << " ( 2 ) ----> " << endl; // getchar();
  
  std::unique_ptr< std::vector<gov::fnal::uboone::datatypes::MuCSData> > mucsdatacol(new std::vector<gov::fnal::uboone::datatypes::MuCSData>);
    
  Int_t fadc = 0;
  
  gov::fnal::uboone::datatypes::MuCSData mucsevt( fadc ); 
    
  mucsdatacol->push_back( mucsevt );

  evt.put( std::move( mucsdatacol ) );
  
  return;
  
}

DEFINE_ART_MODULE( MuCSMerger )

#endif

////////////////////////////////////////////////////////////////////////
//
//    The end !
//
////////////////////////////////////////////////////////////////////////
