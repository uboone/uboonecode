/**
 * @file classes.h
 * @brief Dictionary selection for RecoBase
 * 
 * Original author Rob Kutschke, modified by klg
 * 
 * @note The system is not able to deal with
 * `art::Wrapper<std::vector<std::string>>`.
 * The problem is somewhere inside ROOT's REFLEX mechanism
 * and Philippe Canal says that it is (as of March 2010) a known problem.
 * He also says that they do not have any plans to fix it soon.
 * We can always work around it by putting the string inside another object.
 */

#include "art/Persistency/Common/PtrVector.h" 
#include "art/Persistency/Common/Wrapper.h"
#include "art/Persistency/Common/Assns.h"

#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/OpFlash.h"

template class std::pair< art::Ptr<recob::Track>,        art::Ptr<recob::OpFlash>    >;
template class std::pair< art::Ptr<recob::OpFlash>,      art::Ptr<recob::Track>      >;

template class art::Assns<recob::Track,      recob::OpFlash,    void>;
template class art::Assns<recob::OpFlash,    recob::Track,      void>;

template class art::Wrapper< art::Assns<recob::Track,      recob::OpFlash,    void> >;
template class art::Wrapper< art::Assns<recob::OpFlash,    recob::Track,      void> >;
