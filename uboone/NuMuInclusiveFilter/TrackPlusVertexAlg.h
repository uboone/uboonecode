/**
 *  @file   TrackPlusVertexAlg.h
 * 
 *  @brief  This is an algorithm for finding neutrino candidates using tracks and vertices
 * 
 */
#ifndef TrackPlusVertexAlg_h
#define TrackPlusVertexAlg_h

#include "uboone/TPCNeutrinoIDFilter/NeutrinoIDAlgBase.h"

// LArSoft includes
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// Root includes
#include "TH1D.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace neutrinoid
{

/**
 *  @brief  TrackPlusVertexAlg class
 */
class TrackPlusVertexAlg : virtual public NeutrinoIDAlgBase
{
public:
    /**
     *  @brief  Constructor
     * 
     *  @param  pset
     */
    TrackPlusVertexAlg(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    virtual ~TrackPlusVertexAlg();
    
    /**
     *  @brief a handler for the case where the algorithm control parameters are to be reset
     */
    virtual void reconfigure(fhicl::ParameterSet const&);
    
    /**
     *  @brief Set up for "beginJob" phase if requested
     */
    virtual void beginJob(art::ServiceHandle<art::TFileService>&);
    
    /**
     *  @brief Each algorithm may have different objects it wants "produced" so use this to
     *         let the top level producer module "know" what it is outputting
     */
    virtual void produces(art::EDProducer*);

    /**
     *  @brief Given the list of hits this will search for candidate Seed objects and return them
     */
    virtual bool findNeutrinoCandidates(art::Event&) const;
    
    virtual bool inFV(Double_t, Double_t, Double_t) const;
    
    virtual double FlashTrackDist(double&, double, double) const;

private:
    
    /**
     *  @ brief FHICL parameters.
     */
    std::string                fTrackModuleLabel;        ///< Producer of input tracks
    std::string                fVertexModuleLabel;       ///< Producer of input vertices
    std::string                fCosmicModuleLabel;       ///< Producer of cosmic track tags
    std::string                fOpFlashModuleLabel;      ///< Producer of flashes
    double                     fCosmicScoreCut;          ///< Cut value for possible cosmic tag scores
    double                     fNeutrinoVtxTrackDistCut; ///< Cut to select neutrino candidate
    bool                       fDoHists;                 ///< Fill histograms
    
    TH1D*                      fMaxDistHists;            ///< maximum distance all triangles
    TH1D*                      fBestMaxDistHists;        ///< best max dist
    
    art::EDProducer*           fMyProducerModule;        ///< The producer module driving us
    
    /// @{
    /**
     *  @brief Standard useful properties
     */
    geo::GeometryCore const*             m_geometry;           ///< pointer to the Geometry service
    detinfo::DetectorProperties const* m_detector;           ///< Pointer to the detector properties
    /// @}
};

} // namespace lar_cluster3d
#endif
