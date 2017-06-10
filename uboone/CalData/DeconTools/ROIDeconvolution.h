////////////////////////////////////////////////////////////////////////
/// \file   ROIDeconvolution.h
/// \author T. Usher
////////////////////////////////////////////////////////////////////////

#include "uboone/CalData/DeconTools/IDeconvolution.h"
#include "uboone/CalData/DeconTools/IBaseline.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/LArFFT.h"
#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"

namespace uboone_tool
{

class ROIDeconvolution : public IDeconvolution
{
public:
    explicit ROIDeconvolution(const fhicl::ParameterSet& pset);
    
    ~ROIDeconvolution();
    
    void configure(const fhicl::ParameterSet& pset)              override;
    void outputHistograms(art::TFileDirectory&)            const override;
    
    void Deconvolve(IROIFinder::Waveform const&,
                    raw::ChannelID_t,
                    IROIFinder::CandidateROIVec const&,
                    recob::Wire::RegionsOfInterest_t& )    const override;
    
private:
    // Member variables from the fhicl file
    size_t                                                   fFFTSize;                    ///< FFT size for ROI deconvolution
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
