/******************************************************************************
 * @file DrawTreeFiller.cxx
 * @brief Useful class for handling pseudo-vertices between two track/shower origins
 * @author salvatore.porzio@postgrad.manchester.ac.uk
 * @see  DrawTreeFiller.h
 * ****************************************************************************/

// Decay vertex header
#include "DrawTreeFiller.h"

namespace AuxEvent
{
  DrawTreeFiller::DrawTreeFiller()
  {}
  DrawTreeFiller::~DrawTreeFiller()
  {}

  void DrawTreeFiller::Initialize(AuxEvent::EventTreeFiller & etf, int i_hsnID, const AuxVertex::DecayVertex & decayVertex)
  {
    // General
    run = etf.run;
    subrun = etf.subrun;
    event = etf.event;
    nHsnCandidatesInSameEvent = etf.nHsnCandidates;
    hsnID = i_hsnID;

    // Get vector of hits
    std::vector<art::Ptr<recob::Hit>> prong1_hits = decayVertex.GetProngHits(0);
    std::vector<art::Ptr<recob::Hit>> prong2_hits = decayVertex.GetProngHits(1);
    std::vector<art::Ptr<recob::Hit>> thisTot_hits = decayVertex.GetTotHits();

    // Fill prong1
    prong1_hits_p0_wireCoordinates.clear();
    prong1_hits_p0_tickCoordinates.clear();
    prong1_hits_p1_wireCoordinates.clear();
    prong1_hits_p1_tickCoordinates.clear();
    prong1_hits_p2_wireCoordinates.clear();
    prong1_hits_p2_tickCoordinates.clear();
    for (auto hit : prong1_hits)
    {
      if (hit->View() == 0) {
        prong1_hits_p0_wireCoordinates.push_back(hit->Channel());
        prong1_hits_p0_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
      if (hit->View() == 1) {
        prong1_hits_p1_wireCoordinates.push_back(hit->Channel());
        prong1_hits_p1_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
      if (hit->View() == 2) {
        prong1_hits_p2_wireCoordinates.push_back(hit->Channel());
        prong1_hits_p2_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
    }

    // Fill prong2
    prong2_hits_p0_wireCoordinates.clear();
    prong2_hits_p0_tickCoordinates.clear();
    prong2_hits_p1_wireCoordinates.clear();
    prong2_hits_p1_tickCoordinates.clear();
    prong2_hits_p2_wireCoordinates.clear();
    prong2_hits_p2_tickCoordinates.clear();
    for (auto hit : prong2_hits)
    {
      if (hit->View() == 0) {
        prong2_hits_p0_wireCoordinates.push_back(hit->Channel());
        prong2_hits_p0_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
      if (hit->View() == 1) {
        prong2_hits_p1_wireCoordinates.push_back(hit->Channel());
        prong2_hits_p1_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
      if (hit->View() == 2) {
        prong2_hits_p2_wireCoordinates.push_back(hit->Channel());
        prong2_hits_p2_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
    }

    // Fill prong2
    tot_hits_p0_wireCoordinates.clear();
    tot_hits_p0_tickCoordinates.clear();
    tot_hits_p1_wireCoordinates.clear();
    tot_hits_p1_tickCoordinates.clear();
    tot_hits_p2_wireCoordinates.clear();
    tot_hits_p2_tickCoordinates.clear();
    for (auto hit : thisTot_hits)
    {
      if (hit->View() == 0) {
        tot_hits_p0_wireCoordinates.push_back(hit->Channel());
        tot_hits_p0_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
      if (hit->View() == 1) {
        tot_hits_p1_wireCoordinates.push_back(hit->Channel());
        tot_hits_p1_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
      if (hit->View() == 2) {
        tot_hits_p2_wireCoordinates.push_back(hit->Channel());
        tot_hits_p2_tickCoordinates.push_back((hit->StartTick() + hit->EndTick())/2.);
      }
    }

    // Get coordinates
    dv_p0_wireCoordinates = decayVertex.fChannelLoc[0];
    dv_p0_tickCoordinates = decayVertex.fTickLoc[0];
    dv_p1_wireCoordinates = decayVertex.fChannelLoc[1];
    dv_p1_tickCoordinates = decayVertex.fTickLoc[1];
    dv_p2_wireCoordinates = decayVertex.fChannelLoc[2];
    dv_p2_tickCoordinates = decayVertex.fTickLoc[2];

    prong1_p0_wireCoordinates = decayVertex.fProngChannelLoc[0][0];
    prong1_p0_tickCoordinates = decayVertex.fProngTickLoc[0][0];
    prong1_p1_wireCoordinates = decayVertex.fProngChannelLoc[0][1];
    prong1_p1_tickCoordinates = decayVertex.fProngTickLoc[0][1];
    prong1_p2_wireCoordinates = decayVertex.fProngChannelLoc[0][2];
    prong1_p2_tickCoordinates = decayVertex.fProngTickLoc[0][2];

    prong2_p0_wireCoordinates = decayVertex.fProngChannelLoc[1][0];
    prong2_p0_tickCoordinates = decayVertex.fProngTickLoc[1][0];
    prong2_p1_wireCoordinates = decayVertex.fProngChannelLoc[1][1];
    prong2_p1_tickCoordinates = decayVertex.fProngTickLoc[1][1];
    prong2_p2_wireCoordinates = decayVertex.fProngChannelLoc[1][2];
    prong2_p2_tickCoordinates = decayVertex.fProngTickLoc[1][2];

  } // END function initialize
} // END namespace DrawTreeFiller 
