////////////////////////////////////////////////////////////////////////
//
// @file FindPandoraNuVertexAna_module.cc
//
// @brief An analyzer that runs through the pandora hierarchy and tries to find a numu CC inclusive candidate event
//
// @authors aschu@fnal.gov (cloned from an example)
//
///////////////////////////////////////////////////////////////////////

#ifndef  FindPandoraNuVertexAna_Module
#define  FindPandoraNuVertexAna_Module

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Core/FileBlock.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"
#include "cetlib/cpu_timer.h"

// LArSoft includes
#include "lardata/RecoBase/PFParticle.h"
#include "lardata/RecoBase/Seed.h"
#include "lardata/RecoBase/Hit.h"
#include "lardata/RecoBase/Cluster.h"
#include "lardata/RecoBase/Track.h"
#include "lardata/RecoBase/Vertex.h"
#include "lardata/AnalysisBase/CosmicTag.h"
#include "larcore/SimpleTypesAndConstants/geo_types.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <http://root.cern.ch/root/html532/ClassIndex.html>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"

namespace lar_pandora{class LArPandoraHelper;}

namespace FindPandoraNuVertexAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition

class  FindPandoraNuVertexAna : public art::EDAnalyzer
{
public:
 
    // Standard constructor and destructor for an ART module.
    explicit  FindPandoraNuVertexAna(fhicl::ParameterSet const& pset);
    virtual ~ FindPandoraNuVertexAna();

    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun(const art::Run& run);

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);
    
    // Override the response to input and output fils so we can get the
    // fully qualified path for argo
    void respondToOpenInputFile(art::FileBlock const&);
    void respondToOpenOutputFile(art::FileBlock const&);

    // The analysis routine, called once per event. 
    void analyze (const art::Event& evt); 

private:

    // Need vectors here because we have have several instantiations
    // of the module running depending on matching
    std::string fInputFileName;
    
    // Pointers to the histograms we'll create. 

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;

    // Other variables that will be shared between different methods.
    geo::GeometryCore const*             fGeometry;           ///< pointer to the Geometry service
    detinfo::DetectorProperties const* fDetectorProperties; ///< Pointer to the detector properties

}; // class  FindPandoraNuVertexAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
 FindPandoraNuVertexAna:: FindPandoraNuVertexAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
 FindPandoraNuVertexAna::~ FindPandoraNuVertexAna()
{}
   
//-----------------------------------------------------------------------
void  FindPandoraNuVertexAna::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
//    fPDGCodeHist        = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
}
   
//-----------------------------------------------------------------------
void  FindPandoraNuVertexAna::beginRun(const art::Run& /*run*/)
{
}

//-----------------------------------------------------------------------
void  FindPandoraNuVertexAna::reconfigure(fhicl::ParameterSet const& pset)
{
    // For now require that we input the fully qualified input file name, including full path to file
    // **TODO** learn how to recover from art framework
    //fInputFileName = pset.get<std::string>("FullyQualifiedInputFile");
    
    fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    
    return;
}
    
void FindPandoraNuVertexAna::respondToOpenInputFile(art::FileBlock const& fileBlock)
{
    // Override the fhicl parameter for the input file name
    //fInputFileName = fileBlock.fileName();

    return;
}
    
void FindPandoraNuVertexAna::respondToOpenOutputFile(art::FileBlock const& fileBlock)
{
    // TODO
    return;
}

//-----------------------------------------------------------------------
void  FindPandoraNuVertexAna::analyze(const art::Event& event)
{
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "This is event " << fEvent << " in run " << fRun << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    lar_pandora::PFParticleVector particleVector;
    lar_pandora::LArPandoraHelper::CollectPFParticles(event, "pandoraNu", particleVector);

    // Get the reconstructed vertices
    // ==============================
    lar_pandora::VertexVector vertexVector;
    lar_pandora::PFParticlesToVertices particlesToVertices;
    lar_pandora::LArPandoraHelper::CollectVertices(event, "pandoraNu", vertexVector, particlesToVertices);

    // Get the reconstructed tracks
    // ============================
    lar_pandora::TrackVector trackVector, trackVector2;
    lar_pandora::PFParticlesToTracks particlesToTracks;
    lar_pandora::TracksToHits tracksToHits;
    lar_pandora::LArPandoraHelper::CollectTracks(event, "pandoraNu", trackVector, particlesToTracks);
    lar_pandora::LArPandoraHelper::CollectTracks(event, "pandoraNu", trackVector2, tracksToHits);

    for (unsigned int n = 0; n < particleVector.size(); ++n)
    {
        const art::Ptr<recob::PFParticle> particle = particleVector.at(n);
        std::cout << "Index: " << n << std::endl;
        std::cout << "PDG Code: " << particle->PdgCode() << std::endl;
        std::cout << "Is primary? " << particle->IsPrimary() << std::endl;
        std::cout << "No Daughters: " << particle->NumDaughters() << std::endl;

        if(particle->IsPrimary()) {
            lar_pandora::PFParticlesToVertices::const_iterator vIter = particlesToVertices.find(particle);
            if (particlesToVertices.end() != vIter)
            {
                const lar_pandora::VertexVector &vertexVector = vIter->second;
                if (!vertexVector.empty())
                {
                    if (vertexVector.size() !=1)
                        std::cout << " Warning: Found particle with more than one associated vertex " << std::endl;

                    const art::Ptr<recob::Vertex> vertex = *(vertexVector.begin());
                    double xyz[3] = {0.0, 0.0, 0.0} ;
                    vertex->XYZ(xyz);

                    std::cout << "Vertex: " << xyz[0] << "\t" << xyz[1] << "\t" << xyz[2] << std::endl;
                }
            }
        }
    }

    return;
}

DEFINE_ART_MODULE(FindPandoraNuVertexAna)

} // namespace  FindPandoraNuVertexAna

#endif //  FindPandoraNuVertexAna_Module
