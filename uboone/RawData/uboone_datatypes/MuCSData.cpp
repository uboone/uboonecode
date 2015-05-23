
////////////////////////////////////////////////////////////////////////
//
//    trivia : Class to hold the Muon Counter System (MuCS) data, 
//             a.k.a. the best class there's ever been, May 2015 
//    author : Odysseas Kanenas
//    e-mail : kalousis@vt.edu
//
////////////////////////////////////////////////////////////////////////

#include "MuCSData.h"

#include "cetlib/exception.h"

namespace gov {
namespace fnal {
namespace uboone {
namespace datatypes {


MuCSData::MuCSData()
{
  cha = 1;
    
}

MuCSData::~MuCSData()
{}

MuCSData::MuCSData( Int_t fadc )
{
  cha = fadc;
    
}

}}}}

////////////////////////////////////////////////////////////////////////
//
//    The end !
//
////////////////////////////////////////////////////////////////////////
