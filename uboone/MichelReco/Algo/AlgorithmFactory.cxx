#ifndef ALGORITHMFACTORY_CXX
#define ALGORITHMFACTORY_CXX

#include "AlgorithmFactory.h"

namespace michel {

  AlgoDefault::AlgoDefault() {

    msg::MSGLevel_t debugmode = msg::kINFO;

    //MichelRecoManager* mgr = new MichelRecoManager();
    this->SetVerbosity(debugmode);

    FilterStraightLineClusters* filterAlgo = new FilterStraightLineClusters();
    filterAlgo->setMinRMS(0.5);
    this->AddFilteringAlgo(filterAlgo);
    EdgeMerger* mergeAlgo = new EdgeMerger();
    mergeAlgo->SetVerbosity(debugmode);
    this->AddMergingAlgo(mergeAlgo);

    // calculate various cluster parameters...
    CalcTruncated* ctrunk = new CalcTruncated();
    ctrunk->SetVerbosity(debugmode);
    this->AddAlgo(ctrunk);

    // Attach algorithm for boundary finding
    BoundaryFromTQMaxQ* boundaryalgo = new BoundaryFromTQMaxQ();
    boundaryalgo->SetMaxDistancesTruncatedQMaxQ(15);
    boundaryalgo->SetVerbosity(debugmode);
    this->AddAlgo(boundaryalgo);
    
    RequireBoundaryInLowCov* lowcovbound = new RequireBoundaryInLowCov();
    lowcovbound->SetMaxCovarianceAtStart(0.8);
    lowcovbound->SetVerbosity(debugmode);
    this->AddAlgo(lowcovbound);

    // Attach algorithm for finding michel cluster
    ForwardMichelID* findMichel = new ForwardMichelID();
    findMichel->SetMaxMichelHits(0);
    findMichel->SetVerbosity(debugmode);
    this->AddAlgo(findMichel);

    // MID finding algorithms
    DecideIfStoppingMuon* midalgo = new DecideIfStoppingMuon();
    midalgo->SetChiMin     ( 0.9 );
    midalgo->SetFracMinHits( 0.7 );
    midalgo->SetHitRadius  ( 30  );
    midalgo->SetMaxDist    ( 3.0 );
    midalgo->SetMinBadHits ( 10  );
    midalgo->SetVerbosity(debugmode);
    this->AddAlgo(midalgo);

    // BraggArea filter algo
    FindBraggPeak* braggalgo = new FindBraggPeak();
    braggalgo->SetMinBraggArea(1000.);
    braggalgo->SetVerbosity(debugmode);
    this->AddAlgo(braggalgo);

    // MINIMUM LENGTH REQUIREMENT
    CutOnMuonLength* minlength = new CutOnMuonLength();
    minlength->SetMinMuonLength(10);
    minlength->SetVerbosity(debugmode);
    this->AddAlgo(minlength);

    // MUON LINEARITY CUT
    CutOnMuonLinearity* minlinearity = new CutOnMuonLinearity();
    minlinearity->SetChiMin     ( 0.8 );
    minlinearity->SetFracMinHits( 0.5 );
    minlinearity->SetVerbosity(debugmode);
    this->AddAlgo(minlinearity);

    /*
    // MID filter that removes michels close to wire gaps/edges
    CutOnFiducialVolume* fidvolfilter = new CutOnFiducialVolume();

    wires_to_exclude_min, \
        wires_to_exclude_max, \
        times_to_exclude_min, \
        times_to_exclude_max = fidparser.list_wires_times_to_exclude();
    
    fidvolfilter.SetExcludedWireRanges(wires_to_exclude_min,wires_to_exclude_max);
    fidvolfilter.SetExcludedTimeRanges(times_to_exclude_min,times_to_exclude_max);
    this->AddAlgo(fidvolfilter);
    */

    // Attach algorithm to recluster michel
    SuperSonicClusterer* supersonic = new SuperSonicClusterer();
    supersonic->SetMergeTillConverge(true);
    supersonic->SetUseHitRadius     (true);
    supersonic->SetMaxRadius( 15 );
    supersonic->SetHitRadius(  3 );
    supersonic->SetVerbosity(debugmode);
    this->AddAlgo(supersonic);

    // cone-finding algorithm
    ConeHitFinder* conefinder = new ConeHitFinder();
    conefinder->SetMaxRadius               ( 20 );
    conefinder->SetMaxPerpendicularDistance(  3 );
    conefinder->SetVerbosity(debugmode);
    this->AddAlgo(conefinder);

    // final mid algo cutting on num of michel hits
    CutOnMichelNumHits* michelhits = new CutOnMichelNumHits();
    michelhits->SetMinMichelHits (  5 );
    michelhits->SetMaxMichelHits ( 35 );
    michelhits->SetVerbosity(debugmode);
    this->AddAlgo(michelhits);

    // remove weird horizontal tracks from PMT
    RemoveFakePMTSignals* pmtremoved = new RemoveFakePMTSignals();
    pmtremoved->SetMaxErrorTime(0.1);
    pmtremoved->SetVerbosity(debugmode);
    this->AddAlgo(pmtremoved);

    // require large angle between michel and muon
    RequireLargeAngle* largeangle = new RequireLargeAngle();
    largeangle->SetMinAngle(30.*3.14/180.);
    largeangle->SetMinStraightMichelHits(5);
    largeangle->SetMuonLengthUsed(10000);
    //largeangle.SetVerbosity(new msg.kDEBUG);
    largeangle->SetVerbosity(debugmode);
    this->AddAlgo(largeangle);

    // remove bragg peak hits
    RemoveBraggPeakHits* removeBragg = new RemoveBraggPeakHits();
    removeBragg->SetMaxRadius(1.);
    removeBragg->SetChargeFactor(2.);
    removeBragg->SetVerbosity(debugmode);
    this->AddAlgo(removeBragg);;

    // cut on Michels with high avg. Q/hit
    CutOnMeanHitCharge* cutonavgq = new CutOnMeanHitCharge();
    cutonavgq->SetMaxAvgQ(400.);
    cutonavgq->SetVerbosity(debugmode);
    this->AddAlgo(cutonavgq);

    // find photons
    PhotonFinder* photonFinder = new PhotonFinder();
    photonFinder->SetMaxRadius(80);
    photonFinder->SetMinDotProduct(0.9);
    photonFinder->SetVerbosity(debugmode);
    this->AddAlgo(photonFinder);
    
    // Separate photon clusters
    ClusterPhotons* clusterphotons = new ClusterPhotons();
    clusterphotons->SetMaxDist(3.);
    clusterphotons->SetVerbosity(debugmode);
    this->AddAlgo(clusterphotons);

    // remove bad photon clusters
    RemoveBadPhotonClusters* badphotons = new RemoveBadPhotonClusters();
    badphotons->SetMaxLinearity(0.8);
    badphotons->SetVerbosity(debugmode);
    this->AddAlgo(badphotons);

  }

}

#endif
