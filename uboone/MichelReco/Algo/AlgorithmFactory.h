#ifndef ALGORITHMFACTORY_H
#define ALGORITHMFACTORY_H

#include <string>
#include <cstdlib>
#include "uboone/MichelReco/Fmwk/MichelRecoManager.h"
#include "BoundaryFromTQMaxQ.h"
#include "CalcTruncated.h"
#include "ChiBoundary.h"
#include "ConeHitFinder.h"
#include "CovarianceFollowBoundary.h"
#include "CutOnFiducialVolume.h"
#include "CutOnMeanHitCharge.h"
#include "CutOnMichelNumHits.h"
#include "CutOnMuonLength.h"
#include "CutOnMuonLinearity.h"
#include "CutOnTotNumHits.h"
#include "DecideIfStoppingMuon.h"
#include "EdgeMerger.h"
#include "FilterStraightLineClusters.h"
#include "FindBraggPeak.h"
#include "ForwardMichelID.h"
#include "MatchBoundaries.h"
#include "MaxQBoundary.h"
#include "PhotonFinder.h"
#include "RadiusMichelCluster.h"
#include "RecoMichelDirection.h"
#include "RemoveBadPhotonClusters.h"
#include "RemoveBraggPeakHits.h"
#include "RemoveFakePMTSignals.h"
#include "RequireBoundaryInLowCov.h"
#include "RequireCloseTruncatedPeaks.h"
#include "RequireCovarianceDip.h"
#include "RequireLargeAngle.h"
#include "RequireSlopeSignFlip.h"
#include "StepAroundCluster.h"
#include "StepSuperSonicCluster.h"
#include "SuperSonicClusterer.h"
#include "TruncatedQBoundary.h"
#include "PhotonFinder.h"
#include "ClusterPhotons.h"
#include "RemoveBadPhotonClusters.h"

namespace michel {


  class AlgoDefault : public MichelRecoManager {

  public:
    AlgoDefault();

    ~AlgoDefault(){}

  };

}

#endif
