////////////////////////////////////////////////////////////////////////
// Class:       ParticleId
// Plugin Type: producer (art v2_05_01)
// File:        ParticleId_module.cc
//
// Generated at Wed Jan 31 11:25:52 2018 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 *
 * \class ParticleId
 *
 * \brief ParticleId producer module
 *
 * \author Kirsty Duffy (kduffy@fnal.gov), Adam Lister (a.lister1@lancaster.ac.uk)
 *
 * \date 2018/04/18
 *
 */

// base includes
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

// data products
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larana/TruncatedMean/Algorithm/TruncMean.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// local includes
#include "uboone/ParticleID/Algorithms/GetDaughterTracksShowers.h"
#include "uboone/ParticleID/Algorithms/FiducialVolume.h"
#include "uboone/ParticleID/Algorithms/PIDA.h"
#include "uboone/ParticleID/Algorithms/Chisq.h"
#include "uboone/ParticleID/Algorithms/Bragg_Likelihood_Estimator.h"
#include "uboone/ParticleID/Algorithms/LandauGaussian.h"
#include "uboone/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

// root includes
#include "TVector3.h"

// cpp includes
#include <memory>
#include <bitset>

namespace UBPID{
  class ParticleId;
}

class UBPID::ParticleId : public art::EDProducer {
  public:
    explicit ParticleId(fhicl::ParameterSet const & p);

    ParticleId();
    ParticleId(ParticleId const &) = delete;
    ParticleId(ParticleId &&) = delete;
    ParticleId & operator = (ParticleId const &) = delete;
    ParticleId & operator = (ParticleId &&) = delete;

    // Required functions.
    void produce(art::Event & e) override;

    std::vector<double> fv;

  private:

    // fcl
    std::string fTrackLabel;
    std::string fCaloLabel;
    double fCutDistance;
    double fCutFraction;
    bool fIsSimSmear;

    // fidvol related
    fidvol::fiducialVolume fid;

    // for Chisq
    particleid::Chisquare chisq;

    // for PIDA
    particleid::PIDA pida;

    // for likelihood-based PID
    particleid::Bragg_Likelihood_Estimator braggcalc;

    // For truncated mean
    TruncMean trm;

    //other
    bool isData;

    fhicl::ParameterSet p_holder;
};


UBPID::ParticleId::ParticleId(fhicl::ParameterSet const & p)
{

  std::cout << "[ParticleID] Note: A plane ID of -1 is expected when a track has no calorimetry" << std::endl;
  std::cout << "                   objects in a plane. Because reconstruction demands hits in two" << std::endl;
  std::cout << "                   planes, this should happen a maximum of one time per track." << std::endl;
  std::cout << "                   If more than one plane has an ID of -1 for a single track," << std::endl;
  std::cout << "                   contact the authors." << std::endl;
  std::cout << "[ParticleID] Note: Note that each PID variable is provided on a per-plane basis," << std::endl;
  std::cout << "                   however we currently do not recommend using induction plane" << std::endl;
  std::cout << "                   Particle ID or calorimetry. Proceed with caution." << std::endl;

  fhicl::ParameterSet const p_fv     = p.get<fhicl::ParameterSet>("FiducialVolume");
  fhicl::ParameterSet const p_labels = p.get<fhicl::ParameterSet>("ProducerLabels");
  fhicl::ParameterSet const p_bragg  = p.get<fhicl::ParameterSet>("BraggAlgo");
  fhicl::ParameterSet const p_chi2pidalg = p.get<fhicl::ParameterSet>("Chi2PIDAlg");

  p_holder = p_chi2pidalg;

  // fcl parameters
  fTrackLabel = p_labels.get< std::string > ("TrackLabel");
  fCaloLabel = p_labels.get< std::string > ("CalorimetryLabel");
  fCutDistance  = p.get< double > ("DaughterFinderCutDistance");
  fCutFraction  = p.get< double > ("DaughterFinderCutFraction");

  fv = fid.setFiducialVolume(fv, p_fv);
  fid.printFiducialVolume(fv);
  braggcalc.configure(p_bragg);
  braggcalc.printConfiguration();

  // this module produces a anab::ParticleID object and
  // an association to the track which produced it
  produces< std::vector<anab::ParticleID> >();
  produces< art::Assns< recob::Track, anab::ParticleID> >();

  std::cout << "[ParticleID] >> Track Label: " << fTrackLabel << std::endl;
  std::cout << "[ParticleID] >> Calorimetry Label: " << fCaloLabel << std::endl;
  std::cout << "[ParticleID] >> Cut Distance : " << fCutDistance << std::endl;
  std::cout << "[ParticleID] >> Cut Fraction : " << fCutFraction << std::endl;

}

void UBPID::ParticleId::produce(art::Event & e)
{

  bool isData = e.isRealData();

  if (!isData) std::cout << "[ParticleID] Running simulated data." << std::endl;
  else std::cout << "[ParticleID] Running physics data." << std::endl;

  // produce collection of particleID objects
  std::unique_ptr< std::vector<anab::ParticleID> > particleIDCollection( new std::vector<anab::ParticleID> );
  std::unique_ptr< art::Assns <recob::Track, anab::ParticleID> > trackParticleIdAssn( new art::Assns<recob::Track, anab::ParticleID> );

  //
  // get handles to needed information
  //

  // tracks...
  art::Handle < std::vector<recob::Track> > trackHandle;
  e.getByLabel(fTrackLabel, trackHandle);
  std::vector< art::Ptr<recob::Track> > trackCollection;
  art::fill_ptr_vector(trackCollection, trackHandle);

  // calorimetry object...
  art::FindManyP<anab::Calorimetry> caloFromTracks(trackHandle, e, fCaloLabel);

  // reserve space for the particle ID collection
  // this seems to help with seg faults, but it's a little unclear why
  particleIDCollection->reserve(trackCollection.size());

  for (auto& track : trackCollection){

    // Skip tracks/events with no valid calorimetry object associated to it.
    // In practice this should be 0 tracks (0 hits = 0 tracks), but this is
    // just here for edge cases.
    if (!caloFromTracks.isValid()){
      std::cout << "[ParticleID] Calorimetry<->Track associations are not valid for this track. Skipping." << std::endl;
      continue;
    }

    std::vector< art::Ptr<anab::Calorimetry> > caloFromTrack = caloFromTracks.at(track->ID());

    std::vector<anab::sParticleIDAlgScores> AlgScoresVec;

    // Variables for ParticleID Class
    std::vector<anab::sParticleIDAlgScores> Bragg_fwd_mu    = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_fwd_p     = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_fwd_pi    = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_fwd_k     = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_bwd_mu    = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_bwd_p     = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_bwd_pi    = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_bwd_k     = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> noBragg_fwd_MIP = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;

    // Variables for storing residual range shift from ParticleID Class
    std::vector<anab::sParticleIDAlgScores> Bragg_fwd_mu_shift    = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_fwd_p_shift     = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_fwd_pi_shift    = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_fwd_k_shift     = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_bwd_mu_shift    = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_bwd_p_shift     = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_bwd_pi_shift    = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> Bragg_bwd_k_shift     = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;

    std::vector<anab::sParticleIDAlgScores> Chi2_mu         = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()};
    std::vector<anab::sParticleIDAlgScores> Chi2_p         = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()};
    std::vector<anab::sParticleIDAlgScores> Chi2_pi         = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()};
    std::vector<anab::sParticleIDAlgScores> Chi2_k         = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()};

    std::vector<anab::sParticleIDAlgScores> PIDAval_mean   = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> PIDAval_median = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> PIDAval_kde    = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;

    std::vector<anab::sParticleIDAlgScores> dEdxtruncmean  = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;
    std::vector<anab::sParticleIDAlgScores> dQdxtruncmean  = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;

    std::vector<anab::sParticleIDAlgScores> trk_depE       = {anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores(), anab::sParticleIDAlgScores()} ;

    anab::sParticleIDAlgScores trk_rangeE_mu;
    anab::sParticleIDAlgScores trk_rangeE_p;
    anab::sParticleIDAlgScores trklen; // only here for ease of validation

    art::Ptr< anab:: Calorimetry > calo;
    int planenum = -1;
    for (auto c : caloFromTrack){
      planenum = c->PlaneID().Plane;
      calo = c;
      // Check that caloFromTrack is a valid object
      if (!calo || planenum < 0 || planenum > 2){
        std::cout << "[ParticleID] Calorimetry on plane " << planenum << " is unavailable. Skipping." << std::endl;
        continue;
      }

      std::vector<double> dEdx = calo->dEdx();
      std::vector<double> dQdx = calo->dQdx();
      std::vector<double> resRange = calo->ResidualRange();
      std::vector<double> trkpitchvec = calo->TrkPitchVec();

      /**
       * Initially wanted to only perform particle ID on tracks which Bragged,
       * in the TPC, but changed direction. This is left here for testing purposes.
       */

      // int nDaughters = GetNDaughterTracks((*trackHandle), track->ID(), fCutDistance, fCutFraction);
      // std::cout << "[ParticleID]  Found track with " << nDaughters << " reconstructed daughters." << std::endl;

      // bool   isContained = false;

      // Evaluate PID only for fully-contained particles
      // TVector3 trackStart = track->Vertex();
      // TVector3 trackEnd = track->End();

      // if (fid.isInFiducialVolume(trackStart, fv) && fid.isInFiducialVolume(trackEnd, fv)){
      //   isContained = true;
      // }

      // if (isContained){

      // Check if particle has reconstructed "daughters" - if it does, there may be no Bragg peak
      // and PID might not be accurate
      // if (nDaughters == 0){

      // std::cout << "[ParticleID]  >> Track is fully contained and has no daughters " << std::endl;

      /**
       * Algorithm 1: Bragg_Likelihood
       * Uses B. Ballers theory, along with landau-gaussian distributions with
       * widths measured from data and simulation to estimate the likelihood for
       * each hit in a track to have come from each particle species.
       */
       double shift_fwd_mu = -999.;
       double shift_fwd_p  = -999.;
       double shift_fwd_pi = -999.;
       double shift_fwd_k  = -999.;
       double shift_bwd_mu = -999.;
       double shift_bwd_p  = -999.;
       double shift_bwd_pi = -999.;
       double shift_bwd_k  = -999.;

      Bragg_fwd_mu.at(planenum).fAlgName      = "BraggPeakLLH";
      Bragg_fwd_p.at(planenum).fAlgName       = "BraggPeakLLH";
      Bragg_fwd_pi.at(planenum).fAlgName      = "BraggPeakLLH";
      Bragg_fwd_k.at(planenum).fAlgName       = "BraggPeakLLH";
      Bragg_bwd_mu.at(planenum).fAlgName      = "BraggPeakLLH";
      Bragg_bwd_p.at(planenum).fAlgName       = "BraggPeakLLH";
      Bragg_bwd_pi.at(planenum).fAlgName      = "BraggPeakLLH";
      Bragg_bwd_k.at(planenum).fAlgName       = "BraggPeakLLH";
      Bragg_fwd_mu.at(planenum).fVariableType = anab::kLikelihood;
      Bragg_fwd_p.at(planenum).fVariableType  = anab::kLikelihood;
      Bragg_fwd_pi.at(planenum).fVariableType = anab::kLikelihood;
      Bragg_fwd_k.at(planenum).fVariableType  = anab::kLikelihood;
      Bragg_bwd_mu.at(planenum).fVariableType = anab::kLikelihood;
      Bragg_bwd_p.at(planenum).fVariableType  = anab::kLikelihood;
      Bragg_bwd_pi.at(planenum).fVariableType = anab::kLikelihood;
      Bragg_bwd_k.at(planenum).fVariableType  = anab::kLikelihood;
      Bragg_fwd_mu.at(planenum).fTrackDir = anab::kForward;
      Bragg_fwd_p.at(planenum).fTrackDir  = anab::kForward;
      Bragg_fwd_pi.at(planenum).fTrackDir = anab::kForward;
      Bragg_fwd_k.at(planenum).fTrackDir  = anab::kForward;
      Bragg_bwd_mu.at(planenum).fTrackDir = anab::kBackward;
      Bragg_bwd_p.at(planenum).fTrackDir  = anab::kBackward;
      Bragg_bwd_pi.at(planenum).fTrackDir = anab::kBackward;
      Bragg_bwd_k.at(planenum).fTrackDir  = anab::kBackward;
      Bragg_fwd_mu.at(planenum).fAssumedPdg   = 13;
      Bragg_fwd_p.at(planenum).fAssumedPdg    = 2212;
      Bragg_fwd_pi.at(planenum).fAssumedPdg   = 211;
      Bragg_fwd_k.at(planenum).fAssumedPdg    = 321;
      Bragg_bwd_mu.at(planenum).fAssumedPdg   = 13;
      Bragg_bwd_p.at(planenum).fAssumedPdg    = 2212;
      Bragg_bwd_pi.at(planenum).fAssumedPdg   = 211;
      Bragg_bwd_k.at(planenum).fAssumedPdg    = 321;
      Bragg_fwd_mu.at(planenum).fValue        = braggcalc.getLikelihood(dEdx, resRange, Bragg_fwd_mu.at(planenum).fAssumedPdg, true, planenum, shift_fwd_mu);
      Bragg_fwd_p.at(planenum).fValue         = braggcalc.getLikelihood(dEdx, resRange, Bragg_fwd_p.at(planenum).fAssumedPdg,  true, planenum, shift_fwd_p);
      Bragg_fwd_pi.at(planenum).fValue        = braggcalc.getLikelihood(dEdx, resRange, Bragg_fwd_pi.at(planenum).fAssumedPdg, true, planenum, shift_fwd_pi);
      Bragg_fwd_k.at(planenum).fValue         = braggcalc.getLikelihood(dEdx, resRange, Bragg_fwd_k.at(planenum).fAssumedPdg,  true, planenum, shift_fwd_k);
      Bragg_bwd_mu.at(planenum).fValue        = braggcalc.getLikelihood(dEdx, resRange, Bragg_bwd_mu.at(planenum).fAssumedPdg, false, planenum, shift_bwd_mu);
      Bragg_bwd_p.at(planenum).fValue         = braggcalc.getLikelihood(dEdx, resRange, Bragg_bwd_p.at(planenum).fAssumedPdg,  false, planenum, shift_bwd_p);
      Bragg_bwd_pi.at(planenum).fValue        = braggcalc.getLikelihood(dEdx, resRange, Bragg_bwd_pi.at(planenum).fAssumedPdg, false, planenum, shift_bwd_pi);
      Bragg_bwd_k.at(planenum).fValue         = braggcalc.getLikelihood(dEdx, resRange, Bragg_bwd_k.at(planenum).fAssumedPdg,  false, planenum, shift_bwd_k);
      Bragg_fwd_mu.at(planenum).fPlaneID      = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_fwd_p.at(planenum).fPlaneID       = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_fwd_pi.at(planenum).fPlaneID      = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_fwd_k.at(planenum).fPlaneID       = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_bwd_mu.at(planenum).fPlaneID      = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_bwd_p.at(planenum).fPlaneID       = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_bwd_pi.at(planenum).fPlaneID      = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_bwd_k.at(planenum).fPlaneID       = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      // Special case: MIP-like probability. fAssumedPdg == 0 tells the Bragg
      // algorithm to use the "No-Bragg" theory case
      noBragg_fwd_MIP.at(planenum).fAlgName = "BraggPeakLLH";
      noBragg_fwd_MIP.at(planenum).fVariableType = anab::kLikelihood;
      noBragg_fwd_MIP.at(planenum).fTrackDir = anab::kForward;
      noBragg_fwd_MIP.at(planenum).fAssumedPdg = 0;
      noBragg_fwd_MIP.at(planenum).fValue = braggcalc.getLikelihood(dEdx, resRange, noBragg_fwd_MIP.at(planenum).fAssumedPdg, true, planenum);
      noBragg_fwd_MIP.at(planenum).fPlaneID = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);

      AlgScoresVec.push_back(Bragg_fwd_mu.at(planenum));
      AlgScoresVec.push_back(Bragg_fwd_p.at(planenum));
      AlgScoresVec.push_back(Bragg_fwd_pi.at(planenum));
      AlgScoresVec.push_back(Bragg_fwd_k.at(planenum));
      AlgScoresVec.push_back(Bragg_bwd_mu.at(planenum));
      AlgScoresVec.push_back(Bragg_bwd_p.at(planenum));
      AlgScoresVec.push_back(Bragg_bwd_pi.at(planenum));
      AlgScoresVec.push_back(Bragg_bwd_k.at(planenum));
      AlgScoresVec.push_back(noBragg_fwd_MIP.at(planenum));

      // Objects to store the residual range shift favoured by the likelihood PID
      Bragg_fwd_mu_shift.at(planenum).fAlgName      = "BraggPeakLLH_shift";
      Bragg_fwd_p_shift.at(planenum).fAlgName       = "BraggPeakLLH_shift";
      Bragg_fwd_pi_shift.at(planenum).fAlgName      = "BraggPeakLLH_shift";
      Bragg_fwd_k_shift.at(planenum).fAlgName       = "BraggPeakLLH_shift";
      Bragg_bwd_mu_shift.at(planenum).fAlgName      = "BraggPeakLLH_shift";
      Bragg_bwd_p_shift.at(planenum).fAlgName       = "BraggPeakLLH_shift";
      Bragg_bwd_pi_shift.at(planenum).fAlgName      = "BraggPeakLLH_shift";
      Bragg_bwd_k_shift.at(planenum).fAlgName       = "BraggPeakLLH_shift";
      Bragg_fwd_mu_shift.at(planenum).fVariableType = anab::kLikelihood;
      Bragg_fwd_p_shift.at(planenum).fVariableType  = anab::kLikelihood;
      Bragg_fwd_pi_shift.at(planenum).fVariableType = anab::kLikelihood;
      Bragg_fwd_k_shift.at(planenum).fVariableType  = anab::kLikelihood;
      Bragg_bwd_mu_shift.at(planenum).fVariableType = anab::kLikelihood;
      Bragg_bwd_p_shift.at(planenum).fVariableType  = anab::kLikelihood;
      Bragg_bwd_pi_shift.at(planenum).fVariableType = anab::kLikelihood;
      Bragg_bwd_k_shift.at(planenum).fVariableType  = anab::kLikelihood;
      Bragg_fwd_mu_shift.at(planenum).fTrackDir = anab::kForward;
      Bragg_fwd_p_shift.at(planenum).fTrackDir  = anab::kForward;
      Bragg_fwd_pi_shift.at(planenum).fTrackDir = anab::kForward;
      Bragg_fwd_k_shift.at(planenum).fTrackDir  = anab::kForward;
      Bragg_bwd_mu_shift.at(planenum).fTrackDir = anab::kBackward;
      Bragg_bwd_p_shift.at(planenum).fTrackDir  = anab::kBackward;
      Bragg_bwd_pi_shift.at(planenum).fTrackDir = anab::kBackward;
      Bragg_bwd_k_shift.at(planenum).fTrackDir  = anab::kBackward;
      Bragg_fwd_mu_shift.at(planenum).fAssumedPdg   = 13;
      Bragg_fwd_p_shift.at(planenum).fAssumedPdg    = 2212;
      Bragg_fwd_pi_shift.at(planenum).fAssumedPdg   = 211;
      Bragg_fwd_k_shift.at(planenum).fAssumedPdg    = 321;
      Bragg_bwd_mu_shift.at(planenum).fAssumedPdg   = 13;
      Bragg_bwd_p_shift.at(planenum).fAssumedPdg    = 2212;
      Bragg_bwd_pi_shift.at(planenum).fAssumedPdg   = 211;
      Bragg_bwd_k_shift.at(planenum).fAssumedPdg    = 321;
      Bragg_fwd_mu_shift.at(planenum).fValue        = shift_fwd_mu;
      Bragg_fwd_p_shift.at(planenum).fValue         = shift_fwd_p;
      Bragg_fwd_pi_shift.at(planenum).fValue        = shift_fwd_pi;
      Bragg_fwd_k_shift.at(planenum).fValue         = shift_fwd_k;
      Bragg_bwd_mu_shift.at(planenum).fValue        = shift_bwd_mu;
      Bragg_bwd_p_shift.at(planenum).fValue         = shift_bwd_p;
      Bragg_bwd_pi_shift.at(planenum).fValue        = shift_bwd_pi;
      Bragg_bwd_k_shift.at(planenum).fValue         = shift_bwd_k;
      Bragg_fwd_mu_shift.at(planenum).fPlaneID      = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_fwd_p_shift.at(planenum).fPlaneID       = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_fwd_pi_shift.at(planenum).fPlaneID      = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_fwd_k_shift.at(planenum).fPlaneID       = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_bwd_mu_shift.at(planenum).fPlaneID      = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_bwd_p_shift.at(planenum).fPlaneID       = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_bwd_pi_shift.at(planenum).fPlaneID      = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      Bragg_bwd_k_shift.at(planenum).fPlaneID       = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);

      AlgScoresVec.push_back(Bragg_fwd_mu_shift.at(planenum));
      AlgScoresVec.push_back(Bragg_fwd_p_shift.at(planenum));
      AlgScoresVec.push_back(Bragg_fwd_pi_shift.at(planenum));
      AlgScoresVec.push_back(Bragg_fwd_k_shift.at(planenum));
      AlgScoresVec.push_back(Bragg_bwd_mu_shift.at(planenum));
      AlgScoresVec.push_back(Bragg_bwd_p_shift.at(planenum));
      AlgScoresVec.push_back(Bragg_bwd_pi_shift.at(planenum));
      AlgScoresVec.push_back(Bragg_bwd_k_shift.at(planenum));

      /**
       * Algorithm 2: Chi2
       */
      std::vector<double> chisqValues = chisq.getChisq(calo, p_holder);

      Chi2_mu.at(planenum).fAlgName = "Chi2";
      Chi2_mu.at(planenum).fVariableType = anab::kGOF;
      Chi2_mu.at(planenum).fAssumedPdg = 13;
      Chi2_mu.at(planenum).fTrackDir = anab::kForward;
      Chi2_mu.at(planenum).fValue = chisqValues.at(0);
      Chi2_mu.at(planenum).fPlaneID = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);

      Chi2_p.at(planenum).fAlgName = "Chi2";
      Chi2_p.at(planenum).fVariableType = anab::kGOF;
      Chi2_p.at(planenum).fAssumedPdg = 2212;
      Chi2_p.at(planenum).fTrackDir = anab::kForward;
      Chi2_p.at(planenum).fValue = chisqValues.at(1);
      Chi2_p.at(planenum).fPlaneID = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);

      Chi2_pi.at(planenum).fAlgName = "Chi2";
      Chi2_pi.at(planenum).fVariableType = anab::kGOF;
      Chi2_pi.at(planenum).fAssumedPdg = 211;
      Chi2_pi.at(planenum).fTrackDir = anab::kForward;
      Chi2_pi.at(planenum).fValue = chisqValues.at(2);
      Chi2_pi.at(planenum).fPlaneID = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);

      Chi2_k.at(planenum).fAlgName = "Chi2";
      Chi2_k.at(planenum).fVariableType = anab::kGOF;
      Chi2_k.at(planenum).fAssumedPdg = 321;
      Chi2_k.at(planenum).fTrackDir = anab::kForward;
      Chi2_k.at(planenum).fValue = chisqValues.at(3);
      Chi2_k.at(planenum).fPlaneID = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);

      AlgScoresVec.push_back(Chi2_mu.at(planenum));
      AlgScoresVec.push_back(Chi2_p.at(planenum));
      AlgScoresVec.push_back(Chi2_pi.at(planenum));
      AlgScoresVec.push_back(Chi2_k.at(planenum));

      /**
       * Algorithm 3: PIDA
       * This makes use of Bruce home-brewed PIDA calculation, which can be
       * calculated via three methods:
       * (1) mean (original implementation from B. Baller)
       * (2) median (T. Yang & V. Meddage)
       * (3) kernel density estimator (A. Lister)
       */
      // mean
      PIDAval_mean.at(planenum).fAlgName = "PIDA_mean";
      PIDAval_mean.at(planenum).fVariableType = anab::kPIDA;
      PIDAval_mean.at(planenum).fTrackDir = anab::kForward;
      PIDAval_mean.at(planenum).fValue = pida.getPida(dEdx, resRange, "mean");
      PIDAval_mean.at(planenum).fPlaneID = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      AlgScoresVec.push_back(PIDAval_mean.at(planenum));

      // median
      PIDAval_median.at(planenum).fAlgName = "PIDA_median";
      PIDAval_median.at(planenum).fVariableType = anab::kPIDA;
      PIDAval_median.at(planenum).fTrackDir = anab::kForward;
      PIDAval_median.at(planenum).fValue = pida.getPida(dEdx, resRange, "median");
      PIDAval_median.at(planenum).fPlaneID = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      AlgScoresVec.push_back(PIDAval_median.at(planenum));

      // median
      PIDAval_kde.at(planenum).fAlgName = "PIDA_kde";
      PIDAval_kde.at(planenum).fVariableType = anab::kPIDA;
      PIDAval_kde.at(planenum).fTrackDir = anab::kForward;
      PIDAval_kde.at(planenum).fValue = pida.getPida(dEdx, resRange, "kde");
      PIDAval_kde.at(planenum).fPlaneID = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);
      AlgScoresVec.push_back(PIDAval_kde.at(planenum));

      /**
       * Algorithm 4: Truncated mean dE/dx versus track length
       * Makes use of the "Truncated Mean" algorithm developed by D. Caratelli
       * to plot the truncated mean dE/dx  of a track
       * versus its length for separation.
       */
      size_t nmin = 1;
      size_t nmax = 1;
      const size_t currentiteration = 0;
      const size_t lmin = 1;
      const float convergencelimit = 0.1;
      const float nsigma = 1.0;

      dEdxtruncmean.at(planenum).fAlgName = "TruncatedMean";
      dEdxtruncmean.at(planenum).fVariableType = anab::kdEdxtruncmean;
      dEdxtruncmean.at(planenum).fTrackDir = anab::kForward;
      if (dEdx.size()>0)
      {
        // Convert dEdx vector from double to a float. This is a bad hack because the truncated mean algorithm expects a float as input but dEdx is stored as a vector of doubles
        std::vector<float> dEdx_float(dEdx.begin(),dEdx.end());
        // Now calculate truncated mean
        dEdxtruncmean.at(planenum).fValue = (double)trm.CalcIterativeTruncMean(dEdx_float, nmin, nmax, currentiteration, lmin, convergencelimit, nsigma);
      }
      dEdxtruncmean.at(planenum).fPlaneID = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);

      dQdxtruncmean.at(planenum).fAlgName = "TruncatedMean";
      dQdxtruncmean.at(planenum).fVariableType = anab::kdQdxtruncmean;
      if (dQdx.size()>0)
      {
        // Convert dQdx vector from double to a float. This is a bad hack because the truncated mean algorithm expects a float as input but dQdx is stored as a vector of doubles
        std::vector<float> dQdx_float(dQdx.begin(),dQdx.end());
        // Now calculate truncated mean
        dQdxtruncmean.at(planenum).fValue = (double)trm.CalcIterativeTruncMean(dQdx_float, nmin, nmax, currentiteration, lmin, convergencelimit, nsigma);
      }
      dQdxtruncmean.at(planenum).fPlaneID = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);

      AlgScoresVec.push_back(dEdxtruncmean.at(planenum));
      AlgScoresVec.push_back(dQdxtruncmean.at(planenum));


      /**
       * Algorithm 5: Deposited energy vs energy by range
       * Calculate deposited energy from product of dEdx and trkpitchvec vectors
       * (there is a KineticEnergy object in anab::Calorimetry that already
       * does this, but due to a bug it currently does not use the calibrated
       * dEdx)
       */
      double depE = 0;
      for (size_t i_hit=0; i_hit < dEdx.size(); i_hit++){
        depE += dEdx.at(i_hit)*trkpitchvec.at(i_hit);
      }

      trk_depE.at(planenum).fAlgName = "DepEvsRangeE";
      trk_depE.at(planenum).fVariableType = anab::kEdeposited;
      trk_depE.at(planenum).fTrackDir = anab::kNoDirection;
      trk_depE.at(planenum).fValue = depE;
      trk_depE.at(planenum).fPlaneID = UBPID::uB_SinglePlaneGetBitset(c->PlaneID().Plane);

      AlgScoresVec.push_back(trk_depE.at(planenum));

    } // loop calorimetry objects

    /**
     * Now get additional information which isn't on a per-plane basis
     */

    trklen.fAlgName = "TruncatedMean";
    trklen.fVariableType = anab::kTrackLength;
    trklen.fTrackDir = anab::kNoDirection;
    trklen.fPlaneID = UBPID::uB_SinglePlaneGetBitset(2); // dummy
    trklen.fValue = track->Length();
    AlgScoresVec.push_back(trklen);

    /**
     * Get energy estimation by range (code for momentum by range copied from
     * analysistree, then convert momentum to energy)
     * Calculations only exist in TrackMomentumCalculator for muons and protons
     * TrackMomentumCalculator returns GeV, multiply by 1000 to get MeV
     */
    trkf::TrackMomentumCalculator trkm;
    double track_rangeP_mu = trkm.GetTrackMomentum(track->Length(),13)*1000.;
    double track_rangeP_p = trkm.GetTrackMomentum(track->Length(),2212)*1000.;

    /**
     * Now convert P->E
     * From TrackMomentumCalculator::GetTrackMomentum:
     * P = TMath::Sqrt((KE*KE)+(2*M*KE))
     * P = TMath::Sqrt((E*E)-(M*M)) and E = KE+M
     * => KE = TMath::Sqrt((P*P)+(M*M))-M
     * TrackMometumCalculator uses
     * Muon_M = 105.7 MeV,
     * Proton_M = 938.272 MeV
     * so use these values here
     */

    trk_rangeE_mu.fAlgName = "DepEvsRangeE";
    trk_rangeE_mu.fVariableType = anab::kEbyRange;
    trk_rangeE_mu.fTrackDir = anab::kNoDirection;
    trk_rangeE_mu.fAssumedPdg = 13;
    trk_rangeE_mu.fPlaneID = UBPID::uB_SinglePlaneGetBitset(2); // dummy
    trk_rangeE_p.fAlgName = "DepEvsRangeE";
    trk_rangeE_p.fVariableType = anab::kEbyRange;
    trk_rangeE_p.fTrackDir = anab::kNoDirection;
    trk_rangeE_p.fAssumedPdg = 2212;
    trk_rangeE_mu.fPlaneID = UBPID::uB_SinglePlaneGetBitset(2); // dummy
    trk_rangeE_mu.fValue = TMath::Sqrt((track_rangeP_mu*track_rangeP_mu)+(105.7*105.7)) - 105.7;
    trk_rangeE_p.fValue = TMath::Sqrt((track_rangeP_p*track_rangeP_p)+(938.272*938.272)) - 938.272;
    AlgScoresVec.push_back(trk_rangeE_mu);
    AlgScoresVec.push_back(trk_rangeE_p);

    /**
     * Initially wanted to only perform particle ID on tracks which Bragged
     * in the TPC, but changed direction. This is left here for testing purposes.
     */

    /*  }

        } // end if(isContained)
        else{
    // If particle is *not* contained, assume it is a muon
    // Set pdg=13. No other PID variables will be set because those are only filled for contained particles
    pdg = 13;
    }*/

    /**
     * Fill ParticleID object and push back to event
     */

    anab::ParticleID PID_object(AlgScoresVec);
    particleIDCollection->push_back(PID_object);

    util::CreateAssn(*this, e, *particleIDCollection, track, *trackParticleIdAssn);
  }

  e.put(std::move(particleIDCollection));
  e.put(std::move(trackParticleIdAssn));
}

DEFINE_ART_MODULE(UBPID::ParticleId)
