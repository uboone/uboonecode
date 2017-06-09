////////////////////////////////////////////////////////////////////////
/// \file   NoBaseline.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "uboone/CalData/DeconTools/NoBaseline.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

#include <fstream>

namespace uboone_tool
{

//----------------------------------------------------------------------
// Constructor.
NoBaseline::NoBaseline(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
NoBaseline::~NoBaseline()
{
}
    
void NoBaseline::configure(const fhicl::ParameterSet& pset)
{
    return;
}

    
float NoBaseline::GetBaseline(std::vector<float> const& holder,
                              raw::ChannelID_t    channel,
                              size_t              roiStart,
                              size_t              roiLen) const
{
    return 0.;
}
    
void NoBaseline::outputHistograms(art::TFileDirectory& histDir) const
{
    return;
}
    
}
