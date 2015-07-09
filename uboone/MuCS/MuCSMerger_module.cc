
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
#include "MuCSData.h"

using namespace std;

class MuCSMerger;

class MuCSMerger : public art::EDProducer {
public:
  explicit MuCSMerger( fhicl::ParameterSet const &pset );
  virtual ~MuCSMerger();
  
  void reconfigure( fhicl::ParameterSet const &pset ); // override;
  void produce( art::Event &evt ) override;
      
  private:
  
  Int_t group;
  
};

void MuCSMerger::reconfigure( fhicl::ParameterSet const &pset )
  {
    group = pset.get< int >( "group" );
  
  }

MuCSMerger::MuCSMerger( fhicl::ParameterSet const &pset )
// :
// Initialize member data here.
{
  this->reconfigure( pset );
  
  produces< std::vector<MuCS::MuCSData> >();  
  
}

MuCSMerger::~MuCSMerger()
{}

void MuCSMerger::produce( art::Event &evt )
{
  cout << " ( 0 ) ----> " << group << endl; getchar();
  
  std::unique_ptr< std::vector<MuCS::MuCSData> > mucsdatacol(new std::vector<MuCS::MuCSData>);
  
  Float_t t0 = 0;
  
  Float_t adc1[24], adc2[24], adc3[24], adc7[24];
  
  std::vector<Int_t> hits1, hits2, hits3, hits7;
  
  t0 = -1.0;  
  for ( Int_t i=0; i<24; i++ ) 
    { 
      adc1[i]=-1.0; adc2[i]=-1.0;
      adc3[i]=-1.0; adc7[i]=-1.0;
      
    }
  hits1.push_back(-1.0);
  hits2.push_back(-1.0);hits2.push_back(-1.0);
  hits3.push_back(-1.0);hits3.push_back(-1.0);hits3.push_back(-1.0);
  hits7.push_back(-1.0);hits7.push_back(-1.0);hits7.push_back(-1.0);hits7.push_back(-1.0);
    
  MuCS::MuCSData mucsevt( t0, adc1, adc2, adc3, adc7, hits1, hits2, hits3, hits7 ); 
  
  
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
