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

#include "lardataobj/RecoBase/PFParticle.h"

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
  void GetPrimaryPFParticles( const std::vector< recob::PFParticle > & allPFParticles, std::vector< recob::PFParticle > & primaryPFParticles );

  /*!
   *  \brief  Determines if the supplied PFParticle should be persisted into the new collection (given ShouldProduceNeutrinos and ShouldProduceCosmics)
   *
   *  \param  part  the PFParticle in question
   *
   *  \return whether the supplied PFParticle should be added to the new collection
   */
  bool ShouldPersist( const recob::PFParticle & part );

  /*!
   *  \brief  Filters the PFParticles that should be persisted from an input list
   *
   *  \param  inputPFParticles       input vector of PFParticles
   *  \param  persistentPFParticles  output vector of all PFParticles in the input vector that should be persisted
   */
  void GetPFParticlesToPersist( const std::vector< recob::PFParticle > & inputPFParticles, std::vector< recob::PFParticle > & persistentPFParticles );

  /*!
   *  \brief  Produce a mapping between PFParticles and their ID
   *
   *  \param  allPFParticles      input vector of all PFParticles in the event
   *  \param  idToPFParticleMap   output mapping between PFParticles and their IDs
   */
  void GetIdToPFParticleMap( const std::vector< recob::PFParticle > & allPFParticles, std::map< size_t, recob::PFParticle > & idToPFParticleMap );

  /*!
   *  \brief  Collects all PFParticles downstream (children, grandchildren, ...) of a given PFParticle
   *
   *  \param  part                         input PFParticle
   *  \param  idToPFParticleMap            input mapping between PFParticles and their IDs
   *  \param  idToDownstreamPFParticleMap  output mapping between of all downstream PFParticles and their IDs
   */
  void GetDownstreamPFParticles( recob::PFParticle part, const std::map< size_t, recob::PFParticle > & idToPFParticleMap, std::map< size_t, recob::PFParticle > & idToDownstreamPFParticleMap );

  /*!
   *  \brief  Collects all PFParticles downstream (children, grandchildren, ...) of a given vector of PFParticle
   *
   *  \param  inputPFParticles             input vector of PFParticles
   *  \param  idToPFParticleMap            input mapping between PFParticles and their IDs
   *  \param  idToDownstreamPFParticleMap  output mapping between of all downstream PFParticles and their IDs
   */
  void GetDownstreamPFParticles( const std::vector< recob::PFParticle > & inputPFParticles, const std::map< size_t, recob::PFParticle > & idToPFParticleMap, std::map< size_t, recob::PFParticle > & idToDownstreamPFParticleMap );


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
  // Call appropriate produces<>() functions here.
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::produce(art::Event & e)
{
  /// \todo  Error check for invalid handles etc.

  // Get all of the input collections related to the Pandora instance specified
  art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
  e.getByLabel(fInputProducerLabel, pfParticleHandle);
 
  // Determine which top-level PFParticles we want to persist
  std::vector< recob::PFParticle > primaryPFParticles;
  this->GetPrimaryPFParticles( *pfParticleHandle, primaryPFParticles );
  
  std::vector< recob::PFParticle > persistentPrimaryPFParticles;
  this->GetPFParticlesToPersist( primaryPFParticles, persistentPrimaryPFParticles );

  // Collect all daughter PFParticles
  std::map< size_t, recob::PFParticle > idToPFParticleMap;
  this->GetIdToPFParticleMap( *pfParticleHandle, idToPFParticleMap );

  std::map< size_t, recob::PFParticle > idToSelectedPFParticleMap;
  this->GetDownstreamPFParticles( persistentPrimaryPFParticles, idToPFParticleMap, idToSelectedPFParticleMap);

  // Output the selected PFParticles
  std::unique_ptr< std::vector<recob::PFParticle> > outputParticles( new std::vector<recob::PFParticle> );
  for ( std::map< size_t, recob::PFParticle >::const_iterator it=idToSelectedPFParticleMap.begin(); it != idToSelectedPFParticleMap.end(); ++it )
      outputParticles->push_back(it->second);

  e.put(std::move(outputParticles));
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::reconfigure(fhicl::ParameterSet const & p)
{
  fInputProducerLabel     = p.get<std::string>("InputProducerLabel");
  fShouldProduceNeutrinos = p.get<bool>("ShouldProduceNeutrinos", false);
  fShouldProduceCosmics   = p.get<bool>("ShouldProduceCosmics"  , false);
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetPrimaryPFParticles( const std::vector< recob::PFParticle > & allPFParticles, std::vector< recob::PFParticle > & primaryPFParticles )
{
  for ( recob::PFParticle part : allPFParticles )
    if ( part.IsPrimary() )
      primaryPFParticles.push_back( part );
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetPFParticlesToPersist( const std::vector< recob::PFParticle > & inputPFParticles, std::vector< recob::PFParticle > & persistentPFParticles )
{
  for ( recob::PFParticle part : inputPFParticles )
    if ( this->ShouldPersist( part ) )
      persistentPFParticles.push_back( part );
}

// ---------------------------------------------------------------------------------------

bool PandoraSplitting::ShouldPersist( const recob::PFParticle & part ) 
{
  unsigned int pdg = std::abs(part.PdgCode());
  bool isNeutrino = ( pdg == nue || pdg == numu || pdg == nutau );

  if (  isNeutrino && fShouldProduceNeutrinos ) return true;
  if ( !isNeutrino && fShouldProduceCosmics   ) return true;

  return false;
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetIdToPFParticleMap( const std::vector< recob::PFParticle > & allPFParticles, std::map< size_t, recob::PFParticle > & idToPFParticleMap )
{
  for ( recob::PFParticle part : allPFParticles )
      idToPFParticleMap.insert( std::map< size_t, recob::PFParticle >::value_type( part.Self(), part ) );
}

// ---------------------------------------------------------------------------------------
//
void PandoraSplitting::GetDownstreamPFParticles( recob::PFParticle part, const std::map< size_t, recob::PFParticle > & idToPFParticleMap, std::map< size_t, recob::PFParticle > & idToDownstreamPFParticleMap )
{
    // Add part to the downstream map
    if ( idToDownstreamPFParticleMap.find( part.Self() ) == idToDownstreamPFParticleMap.end() )
        idToDownstreamPFParticleMap.insert( std::map< size_t, recob::PFParticle >::value_type( part.Self(), part ) );

    // And all of its daughters
    for ( size_t daughterId : part.Daughters() ) {

        std::map< size_t, recob::PFParticle >::const_iterator daughterIt = idToPFParticleMap.find( daughterId );

        /// \todo handle error properly
        if ( daughterIt == idToPFParticleMap.end() ) std::exit(1);

        this->GetDownstreamPFParticles( daughterIt->second, idToPFParticleMap, idToDownstreamPFParticleMap );
    }
}

// ---------------------------------------------------------------------------------------

void PandoraSplitting::GetDownstreamPFParticles( const std::vector< recob::PFParticle > & inputPFParticles, const std::map< size_t, recob::PFParticle > & idToPFParticleMap, std::map< size_t, recob::PFParticle > & idToDownstreamPFParticleMap )
{
  for ( recob::PFParticle part : inputPFParticles )
    this->GetDownstreamPFParticles( part, idToPFParticleMap, idToDownstreamPFParticleMap );
}

// ---------------------------------------------------------------------------------------
DEFINE_ART_MODULE(PandoraSplitting)
