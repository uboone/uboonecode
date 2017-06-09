///////////////////////////////////////////////////////////////////////
///
/// \file   IBaseline.h
///
/// \brief  This provides an interface for tools which are tasked with
///         finding the baselines in input waveforms, primarily ROI's
///
/// \author T. Usher
///
////////////////////////////////////////////////////////////////////////

#ifndef IBaseline_H
#define IBaseline_H

#include "fhiclcpp/ParameterSet.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t

namespace art
{
    class TFileDirectory;
}

namespace uboone_tool
{
    class IBaseline
    {
    public:
        virtual ~IBaseline() noexcept = default;
        
        virtual void configure(const fhicl::ParameterSet& pset)                                = 0;
        virtual void outputHistograms(art::TFileDirectory&)                              const = 0;
        
        // Find the baseline
        virtual float GetBaseline(std::vector<float> const&, raw::ChannelID_t, size_t, size_t) const = 0;
    };
}

#endif
