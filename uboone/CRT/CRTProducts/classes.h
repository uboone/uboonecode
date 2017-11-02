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
#include "canvas/Persistency/Common/Assns.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "CRTSimData.hh"
#include "CRTHit.hh"
#include "CRTTrack.hh"
#include <utility>
#include <vector>
#include <map>

//
// Only include objects that we would like to be able to put into the event.
// Do not include the objects they contain internally.
//

template class std::vector<crt::CRTSimData>;
template class art::Wrapper< std::vector<crt::CRTSimData> >;

template class std::vector< std::pair<int,double> >;
template class std::map< unsigned char, std::vector< std::pair<int,double> > >;

template class std::vector<crt::CRTHit>;
template class art::Wrapper< std::vector<crt::CRTHit> >;

template class std::vector<crt::CRTTrack>;
template class art::Wrapper< std::vector<crt::CRTTrack> >;

// Added by SDP
template class art::Assns<crt::CRTSimData,sim::AuxDetSimChannel>;
template class art::Assns<sim::AuxDetSimChannel,crt::CRTSimData>;
template class art::Wrapper<art::Assns<crt::CRTSimData,sim::AuxDetSimChannel> >;
template class art::Wrapper<art::Assns<sim::AuxDetSimChannel,crt::CRTSimData> >;
