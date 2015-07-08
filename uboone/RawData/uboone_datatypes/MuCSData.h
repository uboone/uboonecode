
////////////////////////////////////////////////////////////////////////
//
//    trivia : Class to hold the Muon Counter System (MuCS) data, 
//             a.k.a. the best class there's ever been, May 2015 
//    author : Odysseas Kanenas
//    e-mail : kalousis@vt.edu
//
////////////////////////////////////////////////////////////////////////

#ifndef _UBOONETYPES_MUCSDATA_H
#define _UBOONETYPES_MUCSDATA_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <iosfwd>

#include "TMath.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {
 
class MuCSData 
{

 public:
  MuCSData();
  virtual ~MuCSData();  
  
  MuCSData( Float_t t0, Float_t adc1[24], Float_t adc2[24], Float_t adc3[24], Float_t adc7[24], 
	    std::vector<Int_t> hits1, std::vector<Int_t> hits2, std::vector<Int_t> hits3, std::vector<Int_t> hits7 ); 
  
  Float_t T0();
    
 private:
  
  Int_t group;
  
  Float_t ft0;
  
  Float_t fadc1[24];
  Float_t fadc2[24];
  Float_t fadc3[24];
  Float_t fadc7[24];
  
  std::vector<Int_t> fhits1;
  std::vector<Int_t> fhits2;
  std::vector<Int_t> fhits3;
  std::vector<Int_t> fhits7;
    
};

}}}}

#endif 

////////////////////////////////////////////////////////////////////////
//
//    The end !
//
////////////////////////////////////////////////////////////////////////
