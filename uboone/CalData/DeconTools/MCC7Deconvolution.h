////////////////////////////////////////////////////////////////////////
/// \file   MCC7Deconvolution.cc
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "uboone/CalData/DeconTools/IDeconvolution.h"
#include "uboone/CalData/DeconTools/IBaseline.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "lardata/Utilities/LArFFT.h"
#include "uboone/CalData/DeconTools/WaveformPropertiesAlg.h"

namespace uboone_tool
{

class MCC7Deconvolution : public IDeconvolution
{
public:
    explicit MCC7Deconvolution(const fhicl::ParameterSet& pset);
    
    ~MCC7Deconvolution();
    
    void configure(const fhicl::ParameterSet& pset)              override;
    void outputHistograms(art::TFileDirectory&)            const override;
    
    void Deconvolve(IROIFinder::Waveform const&,
                    raw::ChannelID_t,
                    IROIFinder::CandidateROIVec const&,
                    recob::Wire::RegionsOfInterest_t& )    const override;
    
private:

    bool                                                     fDodQdxCalib;                ///< Do we apply wire-by-wire calibration?
    std::string                                              fdQdxCalibFileName;          ///< Text file for constants to do wire-by-wire calibration
    std::map<unsigned int, float>                            fdQdxCalib;                  ///< Map to do wire-by-wire calibration, key is channel
                                                                                          ///< number, content is correction factor
    
    std::unique_ptr<uboone_tool::IBaseline>                  fBaseline;

    const geo::GeometryCore*                                 fGeometry = lar::providerFrom<geo::Geometry>();
    art::ServiceHandle<util::LArFFT>                         fFFT;
    art::ServiceHandle<util::SignalShapingServiceMicroBooNE> fSignalShaping;
};
    
}
