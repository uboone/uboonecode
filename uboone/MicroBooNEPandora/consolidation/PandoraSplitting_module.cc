////////////////////////////////////////////////////////////////////////
// Class:       PandoraSplitting
// Plugin Type: producer (art v2_07_03)
// File:        PandoraSplitting_module.cc
//
// Generated at Thu Sep 21 13:56:07 2017 by Andrew D. Smith using cetskelgen
// from cetlib version v3_00_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/PCAxis.h"

#include <memory>

/*!
 *  \breif   A module that splits collections from a single Pandora instance into multiple collections
 *
 *  This module is used to split the output of the consolidated Pandora reconstruction into separate collections of "neutrinos" and "cosmic-rays".
 *  It can be configured to produce a copy of all Pandora-produced objects associated with a top-level PFParticle that is:
 *
 *  1. ...       identified as a neutrino by Pandora   (set ShouldProduceNeutrinos to true)
 *  2. ... *not* identified as a neutrino by Pandora   (set ShouldProduceCosmics   to true)
 */
class PandoraSplitting;

class PandoraSplitting : public art::EDProducer {
public:
  explicit PandoraSplitting(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PandoraSplitting(PandoraSplitting const &) = delete;
  PandoraSplitting(PandoraSplitting &&) = delete;
  PandoraSplitting & operator = (PandoraSplitting const &) = delete;
  PandoraSplitting & operator = (PandoraSplitting &&) = delete;

  // Required functions.
  virtual void reconfigure(fhicl::ParameterSet const & p);
  void produce(art::Event & e) override;

private:

  /*!
   *  \brief  Filters primary PFParticles from an input list
   *
   *  \param  allPFParticles      input vector of all PFParticles in the event
   *  \param  primaryPFParticles  output vector of all primary PFParticles in the input vector
   */
  void GetPrimaryPFParticles( const art::Handle< std::vector< recob::PFParticle > > & allPFParticles, std::vector< art::Ptr< recob::PFParticle > > & primaryPFParticles );

  /*!
   *  \brief  Determines if the supplied PFParticle should be persisted into the new collection (given ShouldProduceNeutrinos and ShouldProduceCosmics)
   *
   *  \param  part  the PFParticle in question
   *
   *  \return whether the supplied PFParticle should be added to the new collection
   */
  bool ShouldPersist( const art::Ptr< recob::PFParticle > & part );

  /*!
   *  \brief  Filters the PFParticles that should be persisted from an input list
   *
   *  \param  inputPFParticles       input vector of PFParticles
   *  \param  persistentPFParticles  output vector of all PFParticles in the input vector that should be persisted
   */
  void GetPFParticlesToPersist( const std::vector< art::Ptr< recob::PFParticle > > & inputPFParticles, std::vector< art::Ptr< recob::PFParticle > > & persistentPFParticles );

  /*!
   *  \brief  Produce a mapping between PFParticles and their ID
   *
   *  \param  allPFParticles      input vector of all PFParticles in the event
   *  \param  idToPFParticleMap   output mapping between PFParticles and their IDs
   */
  void GetIdToPFParticleMap( const art::Handle< std::vector< recob::PFParticle > > & allPFParticles, std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap );

  /*!
   *  \brief  Collects all PFParticles downstream (children, grandchildren, ...) of a given PFParticle
   *
   *  \param  part                         input PFParticle
   *  \param  idToPFParticleMap            input mapping between PFParticles and their IDs
   *  \param  idToDownstreamPFParticleMap  output mapping between of all downstream PFParticles and their IDs
   */
  void GetDownstreamPFParticles( art::Ptr< recob::PFParticle > part, const std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap, std::map< size_t, art::Ptr< recob::PFParticle > > & idToDownstreamPFParticleMap );

  /*!
   *  \brief  Collects all PFParticles downstream (children, grandchildren, ...) of a given vector of PFParticle
   *
   *  \param  inputPFParticles             input vector of PFParticles
   *  \param  idToPFParticleMap            input mapping between PFParticles and their IDs
   *  \param  idToDownstreamPFParticleMap  output mapping between of all downstream PFParticles and their IDs
   */
  void GetDownstreamPFParticles( const std::vector< art::Ptr< recob::PFParticle > > & inputPFParticles, const std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap, std::map< size_t, art::Ptr< recob::PFParticle > > & idToDownstreamPFParticleMap );

  /*!
   *  \brief  Collects all objects of type U associated to a given object of type T
   *
   *  \param  anObject         an input object of type T with which we want to collect associated objects of type U
   *  \param  associationTtoU  the general input association between objects of type U and T
   *  \param  associatedU      output vector of objects of type U associated with anObject
   */
  template < class T, class U >
  void CollectAssociated( const art::Ptr< T > & anObject, const art::FindManyP< U > & associationTtoU, std::unique_ptr< std::vector< U > > & associatedU );

  /*!
   *  \brief  Makes associations between two collections (A and B) based on the associations between the input collections of the same types
   *
   *  Example usage: 
   *
   *  Suppose you have an input collection of SpacePoints (A) and PFParticles (B) from a single Pandora instance. 
   *  You also have associations between these collections (inputAssnAtoB)
   *
   *  You have already decided which objects in these collections you want to persist. 
   *  These are in your output collections (collectionA and collectionB). 
   *
   *  Now you need to produce output associations (outputAssnAtoB) only between the objects that are in collectionA and collectionB.
   *  First make an empty association, and then use this function as follows:
   *
   *  \code
   *  // Example usage
   *  std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint > > outputParticlesToSpacePoints( new art::Assns<recob::PFParticle, recob::SpacePoint> );
   *  this->MakeAssociation( outputSpacePoints, outputParticles, assnPFParticleSpacePoint, outputParticlesToSpacePoints);
   *  \code
   *
   *  \param  collectionA     an collection of type A
   *  \param  collectionB     an collection of type B
   *  \param  inputAssnAtoB   input assocations between the full collections of type A and B
   *  \param  outputAssnAtoB  outuput associations between collectionA and collectionB
   */
  template < class A, class B >
  void MakeAssociation( const std::unique_ptr< std::vector< A > > & collectionA,  std::unique_ptr< std::vector< B > > & collectionB, const art::FindManyP< A > & inputAssnAtoB, std::unique_ptr< art::Assns< B, A > > & outputAssnAtoB );

  // FHicL congifurable parameters
  std::string     fInputProducerLabel;           ///< Label for the Pandora instance that produced the collections we want to split up
  bool            fShouldProduceNeutrinos;       ///< If we should produce collections related to neutrino top-level PFParticles
  bool            fShouldProduceCosmics;         ///< If we should produce collections related to cosmic (== non-neutrino) top-level PFParticles


  
  // Useful PDG code for readability
  enum Pdg {
    nue   = 12,
    numu  = 14,
    nutau = 16
  };

};

// ---------------------------------------------------------------------------------------

PandoraSplitting::PandoraSplitting(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  reconfigure(p);
  
  // Define which types of collections this the module produces
  produces< std::vector<recob::PFParticle> >();
  produces< std::vector<recob::SpacePoint> >();
  produces< std::vector<recob::Cluster> >();
  produces< std::vector<recob::Seed> >();
  produces< std::vector<recob::Vertex> >();
  produces< std::vector<recob::Track> >(); 
  produces< std::vector<recob::Shower> >();
  produces< std::vector<recob::PCAxis> >();

  /*
  produces< art::Assns<recob::PFParticle, recob::SpacePoint> >();
  produces< art::Assns<recob::PFParticle, recob::Cluster> >();
  produces< art::Assns<recob::PFParticle, recob::Seed> >();
  produces< art::Assns<recob::PFParticle, recob::Vertex> >();
  produces< art::Assns<recob::PFParticle, recob::Track> >();
  produces< art::Assns<recob::PFParticle, recob::Shower> >();
  produces< art::Assns<recob::PFParticle, recob::PCAxis> >();
  produces< art::Assns<recob::Track, recob::Hit> >();
  produces< art::Assns<recob::Shower, recob::Hit> >();
  produces< art::Assns<recob::Shower, recob::PCAxis> >();
  produces< art::Assns<recob::SpacePoint, recob::Hit> >();
  produces< art::Assns<recob::Cluster, recob::Hit> >();
  produces< art::Assns<recob::Seed, recob::Hit> >();
  */
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::produce(art::Event & e)
{
  /// \todo  Error check for invalid handles etc.

  // ---------------------------------------------------------------------------------------
  // Get all of the input collections related to the Pandora instance specified
  // ---------------------------------------------------------------------------------------

  art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
  art::Handle< std::vector<recob::SpacePoint> > spacePointHandle;
  art::Handle< std::vector<recob::Cluster>    > clusterHandle;
  art::Handle< std::vector<recob::Seed>       > seedHandle;
  art::Handle< std::vector<recob::Vertex>     > vertexHandle;
  art::Handle< std::vector<recob::Track>      > trackHandle;
  art::Handle< std::vector<recob::Shower>     > showerHandle;
  art::Handle< std::vector<recob::PCAxis>     > pcAxisHandle;

  e.getByLabel(fInputProducerLabel, pfParticleHandle);
  e.getByLabel(fInputProducerLabel, spacePointHandle);
  e.getByLabel(fInputProducerLabel, clusterHandle);
  e.getByLabel(fInputProducerLabel, seedHandle);
  e.getByLabel(fInputProducerLabel, vertexHandle);
  e.getByLabel(fInputProducerLabel, trackHandle);
  e.getByLabel(fInputProducerLabel, showerHandle);
  e.getByLabel(fInputProducerLabel, pcAxisHandle);

  // Get the associations
  art::FindManyP< recob::SpacePoint > assnPFParticleSpacePoint( pfParticleHandle, e, fInputProducerLabel );
  art::FindManyP< recob::Cluster    > assnPFParticleCluster(    pfParticleHandle, e, fInputProducerLabel );
  art::FindManyP< recob::Seed       > assnPFParticleSeed(       pfParticleHandle, e, fInputProducerLabel );
  art::FindManyP< recob::Vertex     > assnPFParticleVertex(     pfParticleHandle, e, fInputProducerLabel );
  art::FindManyP< recob::Track      > assnPFParticleTrack(      pfParticleHandle, e, fInputProducerLabel );
  art::FindManyP< recob::Shower     > assnPFParticleShower(     pfParticleHandle, e, fInputProducerLabel );
  art::FindManyP< recob::PCAxis     > assnPFParticlePCAxis(     pfParticleHandle, e, fInputProducerLabel );

  // ---------------------------------------------------------------------------------------
  // Identify the PFParticles to persist
  // ---------------------------------------------------------------------------------------

  // Determine which top-level PFParticles we want to persist
  std::vector< art::Ptr< recob::PFParticle > > primaryPFParticles;
  this->GetPrimaryPFParticles( pfParticleHandle, primaryPFParticles );
  
  std::vector< art::Ptr< recob::PFParticle > > persistentPrimaryPFParticles;
  this->GetPFParticlesToPersist( primaryPFParticles, persistentPrimaryPFParticles );

  // Collect all daughter PFParticles to produce the final list of selected PFParticles to persist
  std::map< size_t, art::Ptr< recob::PFParticle > > idToPFParticleMap;
  this->GetIdToPFParticleMap( pfParticleHandle, idToPFParticleMap );

  std::map< size_t, art::Ptr< recob::PFParticle > > idToSelectedPFParticleMap;
  this->GetDownstreamPFParticles( persistentPrimaryPFParticles, idToPFParticleMap, idToSelectedPFParticleMap);

  // ---------------------------------------------------------------------------------------
  // Output the required collections
  // ---------------------------------------------------------------------------------------

  // Setup the output collections
  std::unique_ptr< std::vector<recob::PFParticle> > outputParticles(   new std::vector<recob::PFParticle> );
  std::unique_ptr< std::vector<recob::SpacePoint> > outputSpacePoints( new std::vector<recob::SpacePoint> );
  std::unique_ptr< std::vector<recob::Cluster>    > outputClusters(    new std::vector<recob::Cluster>    );
  std::unique_ptr< std::vector<recob::Seed>       > outputSeeds(       new std::vector<recob::Seed>       );
  std::unique_ptr< std::vector<recob::Vertex>     > outputVertices(    new std::vector<recob::Vertex>     );
  std::unique_ptr< std::vector<recob::Track>      > outputTracks(      new std::vector<recob::Track>      );
  std::unique_ptr< std::vector<recob::Shower>     > outputShowers(     new std::vector<recob::Shower>     ); 
  std::unique_ptr< std::vector<recob::PCAxis>     > outputPCAxes(      new std::vector<recob::PCAxis>     );

  // Setup the output associations
  /*
  std::unique_ptr< art::Assns<recob::PFParticle, recob::SpacePoint > > outputParticlesToSpacePoints( new art::Assns<recob::PFParticle, recob::SpacePoint> );
  std::unique_ptr< art::Assns<recob::PFParticle, recob::Cluster >    > outputParticlesToClusters(    new art::Assns<recob::PFParticle, recob::Cluster>    );
  std::unique_ptr< art::Assns<recob::PFParticle, recob::Seed >       > outputParticlesToSeeds(       new art::Assns<recob::PFParticle, recob::Seed>       );
  std::unique_ptr< art::Assns<recob::PFParticle, recob::Vertex >     > outputParticlesToVertices(    new art::Assns<recob::PFParticle, recob::Vertex>     );
  std::unique_ptr< art::Assns<recob::PFParticle, recob::Track >      > outputParticlesToTracks(      new art::Assns<recob::PFParticle, recob::Track>      );
  std::unique_ptr< art::Assns<recob::PFParticle, recob::Shower >     > outputParticlesToShowers(     new art::Assns<recob::PFParticle, recob::Shower>     );
  std::unique_ptr< art::Assns<recob::PFParticle, recob::PCAxis >     > outputParticlesToPCAxes(      new art::Assns<recob::PFParticle, recob::PCAxis>     );

  std::unique_ptr< art::Assns<recob::Track     , recob::Hit >        > outputTracksToHits(           new art::Assns<recob::Track,      recob::Hit>        );
  std::unique_ptr< art::Assns<recob::Shower    , recob::Hit >        > outputShowersToHits(          new art::Assns<recob::Shower,     recob::Hit>        );
  std::unique_ptr< art::Assns<recob::SpacePoint, recob::Hit >        > outputSpacePointsToHits(      new art::Assns<recob::SpacePoint, recob::Hit>        );
  std::unique_ptr< art::Assns<recob::Cluster   , recob::Hit >        > outputClustersToHits(         new art::Assns<recob::Cluster,    recob::Hit>        );
  std::unique_ptr< art::Assns<recob::Seed      , recob::Hit >        > outputSeedsToHits(            new art::Assns<recob::Seed,       recob::Hit>        );

  std::unique_ptr< art::Assns<recob::Shower    , recob::PCAxis >     > outputShowersToPCAxes(        new art::Assns<recob::Shower,     recob::PCAxis>     );
  */

  // Output the collections related to the selected PFParticles
  for ( std::map< size_t, art::Ptr< recob::PFParticle > >::const_iterator selectedParticleIt=idToSelectedPFParticleMap.begin(); selectedParticleIt != idToSelectedPFParticleMap.end(); ++selectedParticleIt ) {
      art::Ptr< recob::PFParticle > part = selectedParticleIt->second;
    
      // Add the particle to the ouput collection
      outputParticles->push_back( *part );

      // Collect all other associated objects
      this->CollectAssociated( part, assnPFParticleSpacePoint, outputSpacePoints );
      this->CollectAssociated( part, assnPFParticleCluster   , outputClusters    );
      this->CollectAssociated( part, assnPFParticleSeed      , outputSeeds       );
      this->CollectAssociated( part, assnPFParticleVertex    , outputVertices    );
      this->CollectAssociated( part, assnPFParticleTrack     , outputTracks      );
      this->CollectAssociated( part, assnPFParticleShower    , outputShowers     );
      this->CollectAssociated( part, assnPFParticlePCAxis    , outputPCAxes      );
  }

  // Make associations between the output objects
  //this->MakeAssociation( outputSpacePoints, outputParticles, assnPFParticleSpacePoint, outputParticlesToSpacePoints);

  // Put the new collections into the event
  e.put( std::move( outputParticles   ) );
  e.put( std::move( outputSpacePoints ) );
  e.put( std::move( outputClusters    ) );
  e.put( std::move( outputSeeds       ) );
  e.put( std::move( outputVertices    ) );
  e.put( std::move( outputTracks      ) );
  e.put( std::move( outputShowers     ) );
  e.put( std::move( outputPCAxes      ) );
  
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::reconfigure(fhicl::ParameterSet const & p)
{
  fInputProducerLabel     = p.get<std::string>("InputProducerLabel");
  fShouldProduceNeutrinos = p.get<bool>("ShouldProduceNeutrinos", false);
  fShouldProduceCosmics   = p.get<bool>("ShouldProduceCosmics"  , false);
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetPrimaryPFParticles( const art::Handle< std::vector< recob::PFParticle > > & allPFParticles, std::vector< art::Ptr< recob::PFParticle > > & primaryPFParticles )
{
  for( size_t pfPartIdx = 0; pfPartIdx != allPFParticles->size(); pfPartIdx++ ) {
      art::Ptr<recob::PFParticle> part(allPFParticles, pfPartIdx);
      if ( part->IsPrimary() ) 
        primaryPFParticles.push_back( part );
  }
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetPFParticlesToPersist( const std::vector< art::Ptr< recob::PFParticle > > & inputPFParticles, std::vector< art::Ptr< recob::PFParticle > > & persistentPFParticles )
{
  for ( art::Ptr< recob::PFParticle > part : inputPFParticles )
    if ( this->ShouldPersist( part ) )
      persistentPFParticles.push_back( part );
}

// ---------------------------------------------------------------------------------------

bool PandoraSplitting::ShouldPersist( const art::Ptr< recob::PFParticle > & part ) 
{
  unsigned int pdg = std::abs(part->PdgCode());
  bool isNeutrino = ( pdg == nue || pdg == numu || pdg == nutau );

  if (  isNeutrino && fShouldProduceNeutrinos ) return true;
  if ( !isNeutrino && fShouldProduceCosmics   ) return true;

  return false;
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetIdToPFParticleMap( const art::Handle< std::vector< recob::PFParticle > > & allPFParticles, std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap )
{
  for( size_t pfPartIdx = 0; pfPartIdx != allPFParticles->size(); pfPartIdx++ ) {
      art::Ptr<recob::PFParticle> part(allPFParticles, pfPartIdx);
      idToPFParticleMap.insert( std::map< size_t, art::Ptr< recob::PFParticle > >::value_type( part->Self(), part ) );
  }
}

// ---------------------------------------------------------------------------------------
//
void PandoraSplitting::GetDownstreamPFParticles( art::Ptr< recob::PFParticle > part, const std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap, std::map< size_t, art::Ptr< recob::PFParticle > > & idToDownstreamPFParticleMap )
{
    // Add part to the downstream map
    if ( idToDownstreamPFParticleMap.find( part->Self() ) == idToDownstreamPFParticleMap.end() )
        idToDownstreamPFParticleMap.insert( std::map< size_t, art::Ptr< recob::PFParticle > >::value_type( part->Self(), part ) );

    // And all of its daughters
    for ( size_t daughterId : part->Daughters() ) {

        std::map< size_t, art::Ptr< recob::PFParticle > >::const_iterator daughterIt = idToPFParticleMap.find( daughterId );

        /// \todo handle error properly
        if ( daughterIt == idToPFParticleMap.end() ) std::exit(1);

        this->GetDownstreamPFParticles( daughterIt->second, idToPFParticleMap, idToDownstreamPFParticleMap );
    }
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetDownstreamPFParticles( const std::vector< art::Ptr< recob::PFParticle > > & inputPFParticles, const std::map< size_t, art::Ptr< recob::PFParticle > > & idToPFParticleMap, std::map< size_t, art::Ptr< recob::PFParticle > > & idToDownstreamPFParticleMap )
{
  for ( art::Ptr< recob::PFParticle > part : inputPFParticles )
    this->GetDownstreamPFParticles( part, idToPFParticleMap, idToDownstreamPFParticleMap );
}

// ---------------------------------------------------------------------------------------

template < class T, class U >
void PandoraSplitting::CollectAssociated( const art::Ptr< T > & anObject, const art::FindManyP< U > & associationTtoU, std::unique_ptr< std::vector< U > > & associatedU )
{
  std::vector< art::Ptr< U > > associatedObjects = associationTtoU.at( anObject.key() );
  for ( art::Ptr< U > associatedObject : associatedObjects )
      associatedU->push_back( *associatedObject );
}

// ---------------------------------------------------------------------------------------

DEFINE_ART_MODULE(PandoraSplitting)
