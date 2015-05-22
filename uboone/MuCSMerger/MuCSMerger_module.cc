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

using namespace std;

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
  produces< std::vector<raw::AuxDetDigit> >();  
      
}

MuCSMerger::~MuCSMerger()
{}

void MuCSMerger::produce( art::Event &evt )
{
  cout << " ( 2 ) ----> " << endl; // getchar();
    
  std::unique_ptr< std::vector<raw::AuxDetDigit> > mucsdatacol( new std::vector<raw::AuxDetDigit> );
  
  std::vector< short > fADC; fADC.clear(); fADC.push_back( 6 );
    
  unsigned short fChannel = 1;
  
  std::string  fAuxDetName = "MUCS";
  
  unsigned long long fTimeStamp = 1000;
  
  raw::AuxDetDigit mucsevt( fChannel, fADC, fAuxDetName, fTimeStamp ); 
    
  mucsdatacol->push_back( mucsevt );
  
  evt.put( std::move( mucsdatacol ) );
  
  return;
  
}

DEFINE_ART_MODULE(MuCSMerger)

#endif

////////////////////////////////////////////////////////////////////////
//
//    The end !
//
////////////////////////////////////////////////////////////////////////
