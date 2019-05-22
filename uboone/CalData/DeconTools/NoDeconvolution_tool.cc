////////////////////////////////////////////////////////////////////////
/// \file   NoDeconvolution.cc
/// \author J.I. Crespo-Anadón (based on ROI_deconvolution.cc by T. Usher)
////////////////////////////////////////////////////////////////////////

// Tool to bypass deconvolution

#include "uboone/CalData/DeconTools/IDeconvolution.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
//#include "messagefacility/MessageLogger/MessageLogger.h"
//#include "cetlib_except/exception.h"


#include "art/Utilities/make_tool.h"

namespace uboone_tool
{

class NoDeconvolution : public IDeconvolution
{
public:
    explicit NoDeconvolution(const fhicl::ParameterSet& pset);
    
    ~NoDeconvolution();
    
    void configure(const fhicl::ParameterSet& pset)              override;
    void outputHistograms(art::TFileDirectory&)            const override;
    
    void Deconvolve(const IROIFinder::Waveform&,
                    raw::ChannelID_t,
                    IROIFinder::CandidateROIVec const&,
                    recob::Wire::RegionsOfInterest_t& )    const override;
    
private:

};
    
//----------------------------------------------------------------------
// Constructor.
NoDeconvolution::NoDeconvolution(const fhicl::ParameterSet& pset)
{
    configure(pset);
}
    
NoDeconvolution::~NoDeconvolution()
{
}
    
void NoDeconvolution::configure(const fhicl::ParameterSet& pset)
{
  return;
}

void NoDeconvolution::Deconvolve(const IROIFinder::Waveform&        waveform,
                                  raw::ChannelID_t                   channel,
                                  IROIFinder::CandidateROIVec const& roiVec,
                                  recob::Wire::RegionsOfInterest_t&  ROIVec) const
{
    for(auto const& roi : roiVec)
    {
      // add the range into ROIVec
      ROIVec.add_range(roi.first, waveform.begin() + roi.first, waveform.begin() + roi.second);
    } // loop over candidate roi's
    
    return;
}

void NoDeconvolution::outputHistograms(art::TFileDirectory& histDir) const
{
  return;
}
    
DEFINE_ART_CLASS_TOOL(NoDeconvolution)
}
