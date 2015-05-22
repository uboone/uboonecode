
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

#include "TMath.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {

using namespace gov::fnal::uboone;
 
class MuCSData 
{

 public:
  
  MuCSData();
  
 private:
  
  Int_t cha;
  
};

}}}}

#endif 



