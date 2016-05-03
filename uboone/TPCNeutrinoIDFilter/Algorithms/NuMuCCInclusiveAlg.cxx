/**
 *  @file   NuMuCCInclusiveAlg.cxx
 * 
 *  @brief  Implementation of the Track/Vertex Neutrino ID alg
 *          This module outputs associations between vertices
 *          and tracks that are found to be within the cut value
 *
 *  Original implementation September 20, 2015 by usher@slac.stanford.edu
 *  This is based on work by Anne Schukraft (aschu@fnal.gov) and her
 *  Neutrino ID task force
 */

// The main include
#include "uboone/TPCNeutrinoIDFilter/Algorithms/NuMuCCInclusiveAlg.h"

// Framework Includes
#include "art/Framework/Core/FindManyP.h"

// LArSoft includes
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/WireGeo.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/RecoBase/OpFlash.h"
#include "lardata/AnalysisBase/CosmicTag.h"

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace neutrinoid {

NuMuCCInclusiveAlg::NuMuCCInclusiveAlg(fhicl::ParameterSet const &pset) :
    fMyProducerModule(0),
    fGeometry(lar::providerFrom<geo::Geometry>()),
    fDetector(lar::providerFrom<detinfo::DetectorPropertiesService>())
{
    this->reconfigure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

NuMuCCInclusiveAlg::~NuMuCCInclusiveAlg()
{
}

void NuMuCCInclusiveAlg::reconfigure(fhicl::ParameterSet const &inputPset)
{
    // Assume we could be called externally with the top level module's complete parameter set
    const fhicl::ParameterSet& pset = inputPset.get<fhicl::ParameterSet>("NuMuCCInclusiveAlg");
    
    fTrackModuleLabel        = pset.get<std::string> ("TrackModuleLabel");
    fVertexModuleLabel       = pset.get<std::string> ("VertexModuleLabel");
    fOpFlashModuleLabel      = pset.get<std::string> ("OpFlashModuleLabel");
    
    fDistToEdgeX             = fGeometry->DetHalfWidth()   - pset.get<double>("DistToEdgeX",   10.);
    fDistToEdgeY             = fGeometry->DetHalfHeight()  - pset.get<double>("DistToEdgeY",   20.);
    fDistToEdgeZ             = fGeometry->DetLength() / 2. - pset.get<double>("DistToEdgeZ",   10.);
    
    fFlashWidth              = pset.get<double>      ("FlashWidth",                            80.);
    fBeamMin                 = pset.get<double>      ("BeamMin",                              3.55);
    fBeamMax                 = pset.get<double>      ("BeamMax",                              5.15);
    fPEThresh                = pset.get<double>      ("PEThresh",                              50.);
    fMinTrk2VtxDist          = pset.get<double>      ("MinTrk2VtxDist",                         5.);
    fMinTrackLen             = pset.get<double>      ("MinTrackLen",                           75.);
    
    fDoHists                 = pset.get<bool>        ("FillHistograms",                      false);
}
    
void NuMuCCInclusiveAlg::beginJob(art::ServiceHandle<art::TFileService>& tfs)
{
    if (fDoHists)
    {
        // Define the histograms. Putting semi-colons around the title
        // causes it to be displayed as the x-axis label if the histogram
        // is drawn.
        fMaxDistHists     = tfs->make<TH1D>("TriangleMaxDist", "Max distance for each triangle found",            2000, 0, 1000);
        fBestMaxDistHists = tfs->make<TH1D>("TriBestMaxDist",  "Max distance for the best triangle in the event", 2000, 0, 1000);
    }
    
    return;
}
    
void NuMuCCInclusiveAlg::produces(art::EDProducer* owner)
{
    fMyProducerModule = owner;
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::Track> >();
}

    
bool NuMuCCInclusiveAlg::findNeutrinoCandidates(art::Event & event) const
{
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<art::Assns<recob::Vertex, recob::Track> > vertexTrackAssociations(new art::Assns<recob::Vertex, recob::Track>);
    
    // Recover the hanles to the vertex and track collections we want to analyze.
    art::Handle<std::vector<recob::Vertex>>  vertexVecHandle;
    art::Handle<std::vector<recob::Track>>   trackVecHandle;
    art::Handle<std::vector<recob::OpFlash>> flashListHandle;
    
    event.getByLabel(fVertexModuleLabel,    vertexVecHandle);
    event.getByLabel(fTrackModuleLabel,     trackVecHandle);
    
    //----------------------------------------------------
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
    
    if (event.getByLabel(fOpFlashModuleLabel,flashListHandle))
        art::fill_ptr_vector(flashlist, flashListHandle);
    
    // Require valid handles, otherwise nothing to do
    if (vertexVecHandle.isValid() && vertexVecHandle->size() > 0 && trackVecHandle.isValid() && trackVecHandle->size() > 0)
    {
        //----loop over all the flashes and check if there are flashes within the beam
        //window and above the PE threshold
        const recob::OpFlash* flashPtr(0);
        double                flashmax(0);
        bool                  flashtag(false);
        
        for(const auto& opFlash : flashlist)
        {
            std::cout << "Flash time: " << opFlash->Time() << ", PE: " << opFlash->TotalPE() << std::endl;
            
            if (opFlash->Time() > fBeamMin && opFlash->Time() < fBeamMax && opFlash->TotalPE() > fPEThresh)
            {
                flashtag = true;
                
                // Keep track of the largest flash
                if (opFlash->TotalPE() > flashmax)
                {
                    flashPtr = opFlash.get();
                    flashmax = opFlash->TotalPE();
                }
            }
        }  //end of loop over all the flashes
        
        if(flashtag)
        {
            // We need to keep track of the best combination
            // Can we assign art ptrs? I don't think so...
            int    VertexCandidate=-1;
            int    TrackCandidate=-1;
            double TrackCandLength=0;
            double trackstartzcandidate=0;
            double trackstartxcandidate=0;
            double trackstartycandidate=0;
            double trackendzcandidate=0;
            double trackendxcandidate=0;
            double trackendycandidate=0;
            
            //-----------------------------------------------------------
            for(size_t vertexIdx = 0; vertexIdx < vertexVecHandle->size(); vertexIdx++)
            {
                // Recover art ptr to vertex
                art::Ptr<recob::Vertex> vertex(vertexVecHandle, vertexIdx);
                
                // Get the position of the vertex
                // Ultimately we really want the vertex position in a TVector3 object...
                double vertexXYZ[3];
                
                vertex->XYZ(vertexXYZ);
                
                TVector3 vertexPos(vertexXYZ[0],vertexXYZ[1],vertexXYZ[2]);
            	   
                if(inFV(vertexPos.X(),vertexPos.Y(),vertexPos.Z()))
                {
                    // For each vertex we loop over all tracks looking for matching pairs
                    // The outer loop here, then is over one less than all tracks
                    for(size_t trackIdx = 0; trackIdx < trackVecHandle->size(); trackIdx++)
                    {
                        // Work with an art Ptr here
                        art::Ptr<recob::Track> track(trackVecHandle,trackIdx);
                        
                        // so we need to get the track direction sorted out.
                        TVector3 trackPos = track->Vertex();
                        TVector3 trackEnd = track->End();
                        
                        // Take the closer end---------------------------------
                        double trackToVertexDist = (trackPos - vertexPos).Mag();
                        
                        if ((trackEnd - vertexPos).Mag() < trackToVertexDist)
                        {
                            trackPos          = track->End();
                            trackEnd          = track->Vertex();
                            trackToVertexDist = (trackPos - vertexPos).Mag();
                        }
                        
                        //--------------------------------------------------------------------------
                        if(trackToVertexDist<fMinTrk2VtxDist)
                        {
                            if((trackEnd-trackPos).Mag()>TrackCandLength)
                            {
                                TrackCandLength = (trackEnd-trackPos).Mag();
                                TrackCandidate=trackIdx;
                                VertexCandidate=vertexIdx;
                                trackstartzcandidate=trackPos.z();
                                trackstartxcandidate=trackPos.x();
                                trackstartycandidate=trackPos.y();
                                trackendzcandidate=trackEnd.z();
                                trackendxcandidate=trackEnd.x();
                                trackendycandidate=trackEnd.y();
                            }
                        } //end of if track distance is within 5cm
                    }  //end of loop over the tracks
                }  //end of if the vertex is contained
            } //end of loop over all the vertex
            
            if(TrackCandidate > -1)
            {
                bool trackContainedFlag = FlashTrackDist(flashPtr->ZCenter(), trackstartzcandidate, trackendzcandidate) < fFlashWidth;
                
                // Check to see if we think we have a candidate
                if(TrackCandLength>fMinTrackLen && trackContainedFlag && inFV(trackstartxcandidate, trackstartycandidate, trackstartzcandidate) && inFV(trackendxcandidate, trackendycandidate, trackendzcandidate) )
                {
                    // Make an association between the best vertex and the matching tracks
                    art::Ptr<recob::Vertex> vertex(vertexVecHandle,VertexCandidate);
                    art::Ptr<recob::Track>  track(trackVecHandle,TrackCandidate);
                    
                    util::CreateAssn(*fMyProducerModule, event, track, vertex, *vertexTrackAssociations);
                }
            }
        }  //end of if flag
    }
    
    // Add associations to event.
    event.put(std::move(vertexTrackAssociations));
    
    return true;
}
    
bool NuMuCCInclusiveAlg::inFV(double x, double y, double z) const
{
    double distInX = x - fGeometry->DetHalfWidth();
    double distInY = y - fGeometry->DetHalfHeight();
    double distInZ = z - 0.5 * fGeometry->DetLength();
    
    if (fabs(distInX) < fDistToEdgeX && fabs(distInY) < fDistToEdgeY && fabs(distInZ) < fDistToEdgeZ) return true;
    
    return false;
}

//This function returns the distance between a flash and
//a track (in one dimension, here used only for z direction)
double NuMuCCInclusiveAlg::FlashTrackDist(double flash, double start, double end) const
{
    if(end >= start) {
        if(flash < end && flash > start) return 0;
        else return std::min(fabs(flash-start), fabs(flash-end));
    }
    else {
        if(flash > end && flash < start) return 0;
        else return std::min(fabs(flash-start), fabs(flash-end));
    }
}

} // namespace
