// PandoraNuTrack for CCQE analysis
#include "uboone/ErezCCQEana/MyObjects/PandoraNuTrack.h"
#include "uboone/ErezCCQEana/MyObjects/hit.h"
#include "uboone/ErezCCQEana/MyObjects/box.h"
#include "uboone/ErezCCQEana/MyObjects/flash.h"
#include "uboone/ErezCCQEana/MyObjects/GENIEinteraction.h"
#include "uboone/ErezCCQEana/MyObjects/pairVertex.h"
#include "uboone/ErezCCQEana/MyObjects/tripleVertex.h"
#include <vector>

template class std::vector<PandoraNuTrack>;
template class std::vector<hit>;
template class std::vector<box>;
template class std::vector<flash>;
template class std::vector<GENIEinteraction>;
template class std::vector<pairVertex>;
template class std::vector<tripleVertex>;

