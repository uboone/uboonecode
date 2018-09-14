//
// Build a dictionary.
//
// $Id: classes.h,v 1.8 2010/04/12 18:12:28  Exp $
// $Author:  $
// $Date: 2010/04/12 18:12:28 $
// 
// Original author Rob Kutschke, modified by wes
//

#include "canvas/Persistency/Common/Wrapper.h"
#include "CRTSimData.hh"
#include "CRTHit.hh"
#include "CRTTrack.hh"
////#include "MSetCRTFrag.hh"
#include <utility>
#include <vector>
#include <artdaq-core/Data/Fragment.hh>
#include <map>

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

template class std::vector<crt::CRTSimData>;
template class art::Wrapper< std::vector<crt::CRTSimData> >;

template class std::vector< std::pair<int,float> >;
template class std::map< unsigned char, std::vector< std::pair<int,float> > >;

template class std::vector<crt::CRTHit>;
template class art::Wrapper< std::vector<crt::CRTHit> >;

template class std::vector<crt::CRTTrack>;
template class art::Wrapper< std::vector<crt::CRTTrack> >;

////template class std::vector<crt::MSetCRTFrag>;
////template class art::Wrapper< std::vector<crt::MSetCRTFrag> >;

//template class std::vector< std::vector< artdaq::Fragment> >;
//template class art::Wrapper< std::vector< std::vector< artdaq::Fragment> > >;

