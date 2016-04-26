/**
 *  @file   TrackPlusVertexAlg.cxx
 * 
 *  @brief  Implementation of the Track/Vertex Neutrino ID alg
 *          This module outputs associations between vertices
 *          and tracks that are found to be within the cut value
 *
 *  Original implementation September 20, 2015 by usher@slac.stanford.edu
 *  This is based on work by Anne Schukraft (aschu@fnal.gov) and her
 *  Neutrino ID task force
 */

/////
///   To do list...

// 1) "at least 1" track-vertex association instead of "at least 2"
// 2) Import christoph's cuts
// 3) ...

///// Casting of Christoph's cuts
/*
  
  
  
  
*/

// The main include
#include "uboone/NuMuInclusiveFilter/TrackPlusVertexAlg.h"

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
#include "lardata/AnalysisBase/FlashMatch.h"
#include "lardata/AnalysisBase/T0.h"

//------------------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows
//constexpr int kMaxFlashes    = 1000;   //maximum number of flashes

double flash_time[1000];      
double flash_pe[1000];         
double flash_ycenter[1000];    
double flash_zcenter[1000];    
double flash_ywidth[1000];  
double flash_zwidth[1000];  
double flash_timewidth[1000];



namespace neutrinoid {

TrackPlusVertexAlg::TrackPlusVertexAlg(fhicl::ParameterSet const &pset) : fMyProducerModule(0)
{
    this->reconfigure(pset);
    
    m_geometry = lar::providerFrom<geo::Geometry>();
    m_detector = lar::providerFrom<detinfo::DetectorPropertiesService>();
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackPlusVertexAlg::~TrackPlusVertexAlg()
{
}
    
void TrackPlusVertexAlg::reconfigure(fhicl::ParameterSet const &pset)
{
    // Assume we could be called externally with the top level module's complete parameter set
    const fhicl::ParameterSet& myPset = pset.get<fhicl::ParameterSet>("TPCTrackPlusVertexAlg");
    
    fTrackModuleLabel        = myPset.get<std::string> ("TrackModuleLabel",  "pandoraNu");
    fVertexModuleLabel       = myPset.get<std::string> ("VertexModuleLabel", "pandoraNu");    
    fCosmicModuleLabel       = myPset.get<std::string> ("CosmicModuleLabel", "pandoraNuTag");
    fOpFlashModuleLabel      = myPset.get<std::string> ("OpFlashModuleLabel","opflash");
    fCosmicScoreCut          = myPset.get<double>      ("CosmicScoreCut",    0.4);
    fNeutrinoVtxTrackDistCut = myPset.get<double>      ("NuVtxTrackDistCut", 4.5);
    fDoHists                 = myPset.get<bool>        ("FillHistograms",    true);
    
}
    
void TrackPlusVertexAlg::beginJob(art::ServiceHandle<art::TFileService>& tfs)
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
    
void TrackPlusVertexAlg::produces(art::EDProducer* owner)
{
    fMyProducerModule = owner;
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::Track> >();
    

}

//check is the vertex/track ends position with in the Fiducial Volume
bool neutrinoid::TrackPlusVertexAlg::inFV(Double_t x, Double_t y, Double_t z) const {
      art::ServiceHandle<geo::Geometry> geom;      
      //set borders of the FV
      double xmin=10.;
      double xmax=2.*geom->DetHalfWidth()-10.;
      double ymin=-geom->DetHalfHeight()+20.;
      double ymax=geom->DetHalfHeight()-20.;
      double zmin=10.;
      double zmax=geom->DetLength()-10.;

      if((x < xmax) && (x > xmin) && (y < ymax) && (y > ymin) && (z < zmax) && (z > zmin)) return true;
      else return false;
}

//This function returns the distance between a flash and 
//a track (in one dimension, here used only for z direction)
double neutrinoid::TrackPlusVertexAlg::FlashTrackDist(double& flash, double start, double end) const {
      if(end >= start) {
        if(flash < end && flash > start) return 0;
        else return TMath::Min(fabs(flash-start), fabs(flash-end));
      }
      else {
        if(flash > end && flash < start) return 0;
        else return TMath::Min(fabs(flash-start), fabs(flash-end));
      }
}

    
bool TrackPlusVertexAlg::findNeutrinoCandidates(art::Event & event) const
{
    art::ServiceHandle<geo::Geometry> geom;


    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection
    std::unique_ptr<std::vector<recob::Track> > neutrinoTracks(new std::vector<recob::Track>);
    std::unique_ptr<art::Assns<recob::Vertex, recob::Track> > vertexTrackAssociations(new art::Assns<recob::Vertex, recob::Track>);
    
    // Recover the hanles to the vertex and track collections we want to analyze.
    art::Handle< std::vector<recob::Vertex> >  vertexVecHandle;
    art::Handle< std::vector<recob::Track> >   trackVecHandle;
    art::Handle< std::vector<recob::OpFlash> > flashListHandle;
    //    art::Handle< std::vector<anab::T0> >       t0Handle;

    event.getByLabel(fVertexModuleLabel,    vertexVecHandle);
    event.getByLabel(fTrackModuleLabel,     trackVecHandle);

    //----------------------------------------------------
    std::vector<art::Ptr<recob::OpFlash> > flashlist;    

    //NFlashes-->total number of flashes in each event

    //int kMaxFlashes=1000;
    if (event.getByLabel(fOpFlashModuleLabel,flashListHandle))
      art::fill_ptr_vector(flashlist, flashListHandle);
    
    const int NFlashes = flashlist.size();
    
    for (int i = 0; i < NFlashes && i < 1000; ++i){//loop over hits
 
      //for (size_t i = 0; i < NFlashes && i < kMaxFlashes ; ++i){//loop over hits
      flash_time[i]       = flashlist[i]->Time();
      flash_pe[i]         = flashlist[i]->TotalPE();
      flash_ycenter[i]    = flashlist[i]->YCenter();
      flash_zcenter[i]    = flashlist[i]->ZCenter();
      flash_ywidth[i]     = flashlist[i]->YWidth();
      flash_zwidth[i]     = flashlist[i]->ZWidth();
      flash_timewidth[i]  = flashlist[i]->TimeWidth();
      }

    // Require valid handles, otherwise nothing to do
    // means?
    if (vertexVecHandle.isValid() && vertexVecHandle->size() > 0 && trackVecHandle.isValid() && trackVecHandle->size() > 0)
    {
        // Recover associations relating cosmic tags and track
        art::FindManyP<anab::CosmicTag> cosmicAssns(trackVecHandle, event, fCosmicModuleLabel);

        // We need to keep track of the best combination
        // Can we assign art ptrs? I don't think so...
        std::vector<art::Ptr<recob::Vertex>> bestVertexVec;
        std::vector<art::Ptr<recob::Track> > bestTrackVec;
        double bestDistance(100000.);
       

        //define cut variables
        double flashwidth = 80; //cm. Distance flash-track
        //double distcut = 5; //cm. Distance track start/end to vertex
        //double lengthcut = 75; //cm. Length of longest track
        double beammin = 3.55/*-0.36*/; //us. Beam window start
        double beammax = 5.15/*-0.36*/; //us. Beam window end
        double PEthresh = 50; //PE
                 
        double trklen_longest=0;      
        //double trktheta_longest=0;
        //int longtracktag=0;
        //int trackcandidate=0;
	// int vertexcandidate=0;
        double flashmax=0;
        int theflash=-1;
        bool flashtag=false;
 
        std::cout<<"Start looping over all the flashes and get the maximum and flash number"<<std::endl; 
        //----loop over all the flashes and check if there are flashes within the beam
        //window and above the PE threshold
        for(int f=0; f < NFlashes && f < 1000; f++){         
        //for(int f=0; f<NFlashes && f < kMaxFlashes; f++){ 
          //flash time within the beam window????
          if((flash_time[f]> beammin && flash_time[f]<beammax)&&flash_pe[f]>PEthresh) {
               flashtag=true;
               //get the flash index corresponding to the maximum PE
               if(flash_pe[f]>flashmax)
               {
                   theflash=f;
                   flashmax=flash_pe[f];
               }
          }          
        }  //end of loop over all the flashes f
        if(flashtag==false) return false;

        //
        bool flashmatchtag=false;
        //----------------------------------------------------
        // Outer loop is over the vertices in the collection
        for(size_t vertexIdx = 0; vertexIdx < vertexVecHandle->size(); vertexIdx++)
        {
            // Recover art ptr to vertex
            art::Ptr<recob::Vertex> vertex(vertexVecHandle, vertexIdx);
            
            // Get the position of the vertex
            // Ultimately we really want the vertex position in a TVector3 object...
            double vertexXYZ[3];
            
            vertex->XYZ(vertexXYZ);
            
            TVector3 vertexPos(vertexXYZ[0],vertexXYZ[1],vertexXYZ[2]);
	    
            //check if the vertex is within the FV;
            if(!inFV(vertexPos.X(),vertexPos.Y(),vertexPos.Z())) return false;
            
            // For each vertex we loop over all tracks looking for matching pairs
            // The outer loop here, then is over one less than all tracks
            for(size_t track1Idx = 0; track1Idx < trackVecHandle->size() - 1; track1Idx++)
            {
                // Work with an art Ptr here
                art::Ptr<recob::Track> track1(trackVecHandle,track1Idx);
                
                // Is there a CosmicTag associated with this track?
                // There are other/better ways to handle this but I think this covers worst case scenario
                if (cosmicAssns.isValid() && cosmicAssns.size() > 0)
                {
                    std::vector<art::Ptr<anab::CosmicTag> > cosmicVec = cosmicAssns.at(track1.key());
                    
                    if (!cosmicVec.empty())
                    {
                        art::Ptr<anab::CosmicTag>& cosmicTag(cosmicVec.front());
                        
                        if (cosmicTag->CosmicScore() > fCosmicScoreCut) continue;
                    }
                }
                
                // Currently we have the problem that tracks can be fit in the "wrong" direction
                // so we need to get the track direction sorted out.
                TVector3 track1Pos = track1->Vertex();
                TVector3 track1End = track1->End();

                //declare track length/track angle  will be used for the track selection
                //double trklen = track.Length();
                //double trktheta  = track1.Theta();
                // select the trklen and trktheta               
 
                // Take the closer end---------------------------------
                double track1ToVertexDist = (track1Pos - vertexPos).Mag();
                
                if ((track1End - vertexPos).Mag() < track1ToVertexDist)
                {
                    track1Pos          = track1->End();
                    track1End          = track1->Vertex();
                    track1ToVertexDist = (track1Pos - vertexPos).Mag();
                }
                // calculate the cos angle of track1
                double trktheta=fabs(track1End.z()-track1Pos.z())/(track1End-track1Pos).Mag();

                //-------------------------------------------------------
                //get the length of the longest forward going track and close to vertex
                if((track1End-track1Pos).Mag()> trklen_longest && trktheta>0.85 && track1ToVertexDist<5) {
                  trklen_longest=(track1End-track1Pos).Mag();
                  //trackcandidate=track1Idx;
		  //                  vertexcandidate=vertexIdx;
                  //check if the longest track if flash matched.
                  if (FlashTrackDist(flash_zcenter[theflash], track1Pos.z(), track1End.z() )< flashwidth ) {flashmatchtag=true;}

                }  
                // Is there a cut at this point?
                
                // Now loop over rest of tracks looking for best match
                for(size_t track2Idx = track1Idx+1; track2Idx < trackVecHandle->size(); track2Idx++)
                {
                    // Still working the art ptrs
                    art::Ptr<recob::Track> track2(trackVecHandle,track2Idx);
                    
                    // Check cosmic ray tags again
                    if (cosmicAssns.isValid() && cosmicAssns.size() > 0)
                    {
                        std::vector<art::Ptr<anab::CosmicTag> > cosmicVec = cosmicAssns.at(track2.key());
                        
                        if (!cosmicVec.empty())
                        {
                            art::Ptr<anab::CosmicTag>& cosmicTag(cosmicVec.front());
                            
                            if (cosmicTag->CosmicScore() > fCosmicScoreCut) continue;
                        }
                    }
                    
                    // Same dance for closest position
                    TVector3 track2Pos = track2->Vertex();
                    TVector3 track2End = track2->End();
                    
                    // Take the closer end
                    double track2ToVertexDist = (track2Pos - vertexPos).Mag();
                    
                    if ((track2End - vertexPos).Mag() < track2ToVertexDist)
                    {
                        track2Pos          = track2->End();
                        track2End          = track2->Vertex();
                        track2ToVertexDist = (track2Pos - vertexPos).Mag();
                    }
                    
                    // Which distance is larger?
                    double maxDist = std::max(track1ToVertexDist,track2ToVertexDist);
                    
                    // Now also get the distance between the start of the two tracks
                    double track1ToTrack2Dist = (track1Pos - track2Pos).Mag();
                    
                    // is it larger? what is this maxDist for??
                    maxDist = std::max(maxDist,track1ToTrack2Dist);
                    
                    // fill hist if asked
                    if (fDoHists) fMaxDistHists->Fill(maxDist, 1.);
                    
                    // Is this the best?
                    if (maxDist < bestDistance)
                    {
                        double length1 = (track1End - track1Pos).Mag();
                        double length2 = (track2End - track2Pos).Mag();
                        double cos = 0;

                        if(length1 > 10 && length1 > length2) {
                           double length1z = fabs(track1End.z() - track1Pos.z());
                           cos = length1z/length1;
                        }
                        else if(length2 > 10) {
                           double length2z = fabs(track2End.z() - track2Pos.z());
                           cos = length2z/length2;
                        }

                        if(cos > 0.85) {

			   // Clear out the old results
                           bestVertexVec.clear();
                           bestTrackVec.clear();
                        
                           //Now store away
                           bestVertexVec.push_back(vertex);
                           bestTrackVec.push_back(track1);
                           bestTrackVec.push_back(track2);
                           bestDistance = maxDist;
                       }
                    }
                }   //end of loop over the second track
            }  //end of loop over the first track
        } // loop over the vertex
        if (trklen_longest<75) return false;         
        if (flashmatchtag==false) return false;


        // Fill hists if asked
        if (fDoHists) fBestMaxDistHists->Fill(bestDistance, 1.);
        
        // Check to see if we think we have a candidate
        if (bestDistance < fNeutrinoVtxTrackDistCut)
        {
            // Make an association between the best vertex and the matching tracks
            util::CreateAssn(*fMyProducerModule, event, bestVertexVec[0], bestTrackVec, *vertexTrackAssociations);
        }
    }
    
    // Add tracks and associations to event.
    event.put(std::move(vertexTrackAssociations));
    
    return true;
}

} // namespace
