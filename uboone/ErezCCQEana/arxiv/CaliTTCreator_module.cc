//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <utility>      // std::pair, std::make_pair

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TStopwatch.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVectorBase.h"

//"larsoft" object includes
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/RecoAlg/TrajectoryMCSFitter.h"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larsim/MCCheater/BackTracker.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"

#include "uboone/MyClasses/BackTrackerTruthMatch.h"
#include "uboone/MyClasses/TruncMean.h"

namespace mynamespace {
	class TTreeCreator;
}


class mynamespace::TTreeCreator : public art::EDAnalyzer {

	public:
		
		explicit TTreeCreator(fhicl::ParameterSet const & p);
		TTreeCreator(TTreeCreator const &) = delete;
		TTreeCreator(TTreeCreator &&) = delete;
		TTreeCreator & operator = (TTreeCreator const &) = delete;
		TTreeCreator & operator = (TTreeCreator &&) = delete;

		void analyze(art::Event const & e) override;
		void beginJob() override;
		void endJob() override;

	private:

		art::ServiceHandle<art::TFileService>tfs;
		TTree* myTTree;
		int NEvents;
		int RunNumber;
		int SubRunNumber;
		int EventNumber;
		int EventPassedSwTrigger;

		//____________________________________________________________________________________________________________________________________________________________________________________

		// Tracks

		int NumberTracks;
		std::vector<int> Track_IsBeamEvent;
		std::vector<int> Track_IsCosmicEvent;
		std::vector<double> Track_StartX;
		std::vector<double> Track_StartY;
		std::vector<double> Track_StartZ;
		std::vector<double> Track_EndX;
		std::vector<double> Track_EndY;
		std::vector<double> Track_EndZ;
		std::vector<int> Track_ID;
		std::vector<int> Track_ParticleId;
		std::vector<double> Track_Length;
		std::vector<double> Track_Phi;
		std::vector<double> Track_Theta;
		std::vector<double> Track_AzimuthAngle; 
		std::vector<double> Track_ZenithAngle; 
		std::vector<double> Track_VertexX;
		std::vector<double> Track_VertexY;
		std::vector<double> Track_VertexZ;
		std::vector<double> Track_Chi2;
		std::vector<double> Track_Chi2PerNdof;
		std::vector<int> Track_CountValidPoints;
		std::vector<std::vector<double>> Track_UnitaryDirectionVectorCoordinateXAtPoint;
		std::vector<std::vector<double>> Track_UnitaryDirectionVectorCoordinateYAtPoint;
		std::vector<std::vector<double>> Track_UnitaryDirectionVectorCoordinateZAtPoint;
		std::vector<std::vector<double>> Track_ThetaAtPoint;
		std::vector<std::vector<double>> Track_PhiAtPoint;
		std::vector<double> Track_HasMomentum;
		std::vector<double> Track_EndMomentum;
		std::vector<double> Track_EndMomentumVectorX;
		std::vector<double> Track_EndMomentumVectorY;
		std::vector<double> Track_EndMomentumVectorZ;
		std::vector<std::vector<int>> Track_HasPoint;
		std::vector<std::vector<int>> Track_HasValidPoint;
		std::vector<double> Track_FirstPoint;
		std::vector<double> Track_FirstValidPoint;
		std::vector<double> Track_LastPoint;
		std::vector<double> Track_LastValidPoint;
		std::vector<std::vector<double>> Track_LocationVectorXAtPoint;
		std::vector<std::vector<double>> Track_LocationVectorYAtPoint;
		std::vector<std::vector<double>> Track_LocationVectorZAtPoint;
		std::vector<std::vector<double>> Track_MomentumAtPoint;
		std::vector<std::vector<double>> Track_MomentumVectorXAtPoint;
		std::vector<std::vector<double>> Track_MomentumVectorYAtPoint;
		std::vector<std::vector<double>> Track_MomentumVectorZAtPoint;
		std::vector<int> Track_Ndof;
		std::vector<std::vector<int>> Track_NextValidPoint;
		std::vector<std::vector<int>> Track_PreviousValidPoint;
		std::vector<int> Track_NPoints;
		std::vector<int> Track_NumberTrajectoryPoints;
		std::vector<double> Track_Calorimetry_Plane0_KineticEnergy;
		std::vector<double> Track_Calorimetry_Plane1_KineticEnergy;
		std::vector<double> Track_Calorimetry_Plane2_KineticEnergy;
		std::vector<double> Track_Momentum;
		std::vector<double> Track_Momentum_MCS;
		std::vector<std::vector<double>> Track_Calorimetry_Plane0_ResidualRange;
		std::vector<std::vector<double>> Track_Calorimetry_Plane1_ResidualRange;
		std::vector<std::vector<double>> Track_Calorimetry_Plane2_ResidualRange;
		std::vector<std::vector<double>> Track_Calorimetry_Plane0_dEdx;
		std::vector<std::vector<double>> Track_Calorimetry_Plane1_dEdx;
		std::vector<std::vector<double>> Track_Calorimetry_Plane2_dEdx;
		std::vector<std::vector<double>> Track_Calorimetry_Plane0_dQdx;
		std::vector<std::vector<double>> Track_Calorimetry_Plane1_dQdx;
		std::vector<std::vector<double>> Track_Calorimetry_Plane2_dQdx;
		std::vector<std::vector<double>> Track_Calorimetry_Plane0_TruncdEdx;
		std::vector<std::vector<double>> Track_Calorimetry_Plane1_TruncdEdx;
		std::vector<std::vector<double>> Track_Calorimetry_Plane2_TruncdEdx;
		std::vector<std::vector<double>> Track_Calorimetry_Plane0_TruncdQdx;
		std::vector<std::vector<double>> Track_Calorimetry_Plane1_TruncdQdx;
		std::vector<std::vector<double>> Track_Calorimetry_Plane2_TruncdQdx;
		std::vector<double> Track_MCParticle_E;
		std::vector<double> Track_MCParticle_Mass;
		std::vector<double> Track_MCParticle_Momentum_E;
		std::vector<double> Track_MCParticle_Momentum_Px;
		std::vector<double> Track_MCParticle_Momentum_Py;
		std::vector<double> Track_MCParticle_Momentum_Pz;
		std::vector<int> Track_MCParticle_PdgCode;
		std::vector<int> Track_MCParticle_TrackId;
		std::vector<double> Track_MCParticle_Purity;
		std::vector<double> Track_MCParticle_EndX;
		std::vector<double> Track_MCParticle_EndY;
		std::vector<double> Track_MCParticle_EndZ;
		std::vector<double> Track_MCParticle_Vx;
		std::vector<double> Track_MCParticle_Vy;
		std::vector<double> Track_MCParticle_Vz;

		//____________________________________________________________________________________________________________________________________________________________________________________

		// Hits

		int NumberHits;
		std::vector<int> Hit_IsBeamEvent;
		std::vector<int> Hit_IsCosmicEvent;
		std::vector<int> Hit_Channel;
		std::vector<double> Hit_Integral;
		std::vector<double> Hit_Multiplicity;
		std::vector<double> Hit_PeakTime;
		std::vector<double> Hit_PeakAmplitude;
		std::vector<double> Hit_RMS;
		std::vector<double> Hit_Wire;
		std::vector<double> Hit_SigmaPeakAmplitude;
		std::vector<double> Hit_SigmaPeakTime;
		std::vector<int> Hit_WireID_Plane;

		std::vector<int> Hit_MCParticle_PdgCode;
		std::vector<double> Hit_MCParticle_Purity;
		std::vector<double> Hit_MCParticle_EndX;
		std::vector<double> Hit_MCParticle_EndY;
		std::vector<double> Hit_MCParticle_EndZ;
		std::vector<double> Hit_MCParticle_Vx;
		std::vector<double> Hit_MCParticle_Vy;
		std::vector<double> Hit_MCParticle_Vz;

		//____________________________________________________________________________________________________________________________________________________________________________________

		// Beam Flashes

		int NumberFlashesBeam;
		std::vector<double> FlashesBeam_Ywidth;
		std::vector<double> FlashesBeam_Zwidth;
		std::vector<double> FlashesBeam_Twidth;
		std::vector<double> FlashesBeam_Ycenter;
		std::vector<double> FlashesBeam_Zcenter;
		std::vector<double> FlashesBeam_Time;
		std::vector<double> FlashesBeam_TotalPE;
		std::vector<std::vector<double>> FlashesBeam_PE_Per_PMT;

		//____________________________________________________________________________________________________________________________________________________________________________________

		// Cosmic Flashes

		int NumberFlashesCosmic;
		std::vector<double> FlashesCosmic_Ywidth;
		std::vector<double> FlashesCosmic_Zwidth;
		std::vector<double> FlashesCosmic_Twidth;
		std::vector<double> FlashesCosmic_Ycenter;
		std::vector<double> FlashesCosmic_Zcenter;
		std::vector<double> FlashesCosmic_Time;
		std::vector<double> FlashesCosmic_TotalPE;
		std::vector<std::vector<double>> FlashesCosmic_PE_Per_PMT;

		//____________________________________________________________________________________________________________________________________________________________________________________

		// Beam OpHits

		int NumberOpHitsBeam;
		std::vector<double> OpHitsBeam_Amplitude;
		std::vector<double> OpHitsBeam_Area;
		std::vector<double> OpHitsBeam_OpChannel;
		std::vector<double> OpHitsBeam_PE;
		std::vector<double> OpHitsBeam_PeakTime;
		std::vector<double> OpHitsBeam_Width;

		//____________________________________________________________________________________________________________________________________________________________________________________

		// Cosmic OpHits

		int NumberOpHitsCosmic;
		std::vector<double> OpHitsCosmic_Amplitude;
		std::vector<double> OpHitsCosmic_Area;
		std::vector<double> OpHitsCosmic_OpChannel;
		std::vector<double> OpHitsCosmic_PE;
		std::vector<double> OpHitsCosmic_PeakTime;
		std::vector<double> OpHitsCosmic_Width;

		//____________________________________________________________________________________________________________________________________________________________________________________

		// Vertices

		int NumberVertices;
		std::vector<int> Vertex_ID;
		std::vector<int> Vertex_IsBeamEvent;
		std::vector<int> Vertex_IsCosmicEvent;
		std::vector<double> Vertex_PositionX;
		std::vector<double> Vertex_PositionY;
		std::vector<double> Vertex_PositionZ;

		//____________________________________________________________________________________________________________________________________________________________________________________

		// MCTruth Info

		int NumberMCTruthEvents;
		std::vector<int> MCTruth_CCNC;
		std::vector<int> MCTruth_Mode;
		std::vector<double> MCTruth_Pt;
		std::vector<double> MCTruth_QSqr;
		std::vector<int> MCTruth_Target;
		std::vector<double> MCTruth_Theta;
		std::vector<double> MCTruth_X;
		std::vector<double> MCTruth_Y;
		std::vector<double> MCTruth_W;
		std::vector<int> MCTruth_NParticles;
		std::vector<int> MCTruth_Origin;
		std::vector<int> MCTruth_InteractionType;
		std::vector<std::vector<int> > MCTruth_Particle_Pdg;
		std::vector<std::vector<int> > MCTruth_Particle_TrackId;
		std::vector<std::vector<double> > MCTruth_Particle_EndX;
		std::vector<std::vector<double> > MCTruth_Particle_EndY;
		std::vector<std::vector<double> > MCTruth_Particle_EndZ;
		std::vector<std::vector<double> > MCTruth_Particle_Vx;
		std::vector<std::vector<double> > MCTruth_Particle_Vy;
		std::vector<std::vector<double> > MCTruth_Particle_Vz;
		std::vector<double> MCTruth_Particle_Nu_Vx;
		std::vector<double> MCTruth_Particle_Nu_Vy;
		std::vector<double> MCTruth_Particle_Nu_Vz;

		//____________________________________________________________________________________________________________________________________________________________________________________

		// PFParticle

		/*int NumberPFParticles;
		std::vector<int> PFParticle_IsBeamEvent;
		std::vector<int> PFParticle_IsCosmicEvent;
		std::vector<int> PFParticle_IsPrimary;
		std::vector<int> PFParticle_NumDaughters;
		std::vector<double> PFParticle_Track_Length;
		std::vector<double> PFParticle_Track_Theta;
		std::vector<double> PFParticle_Track_Phi;
		std::vector<double> PFParticle_Track_StartX;
		std::vector<double> PFParticle_Track_StartY;
		std::vector<double> PFParticle_Track_StartZ;
		std::vector<double> PFParticle_Track_EndX;
		std::vector<double> PFParticle_Track_EndY;
		std::vector<double> PFParticle_Track_EndZ;*/

};
// ____________________________________________________________________________________________________________________________________________________________________________________________________

mynamespace::TTreeCreator::TTreeCreator(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p) {

}
// ____________________________________________________________________________________________________________________________________________________________________________________________________

void mynamespace::TTreeCreator::beginJob()
{
	myTTree = tfs->make<TTree>("myTTree", "myTTree");
// ____________________________________________________________________________________________________________________________________________________________________________________________________

	myTTree->Branch("RunNumber", &RunNumber, "RunNumber/I");
	myTTree->Branch("SubRunNumber", &SubRunNumber, "SubRunNumber/I");
	myTTree->Branch("EventNumber", &EventNumber, "EventNumber/I");
	myTTree->Branch("EventPassedSwTrigger", &EventPassedSwTrigger, "EventPassedSwTrigger/I");

	// Tracks

	myTTree->Branch("NumberTracks", &NumberTracks, "NumberTracks/I");
	myTTree->Branch("Track_IsBeamEvent", &Track_IsBeamEvent);
	myTTree->Branch("Track_IsCosmicEvent", &Track_IsCosmicEvent);
	myTTree->Branch("Track_StartX", &Track_StartX);
	myTTree->Branch("Track_StartY", &Track_StartY);
	myTTree->Branch("Track_StartZ", &Track_StartZ);
	myTTree->Branch("Track_EndX", &Track_EndX);
	myTTree->Branch("Track_EndY", &Track_EndY);
	myTTree->Branch("Track_EndZ", &Track_EndZ);
	myTTree->Branch("Track_ID", &Track_ID);
	myTTree->Branch("Track_ParticleId", &Track_ParticleId);
	myTTree->Branch("Track_Length", &Track_Length);
	myTTree->Branch("Track_Phi", &Track_Phi);
	myTTree->Branch("Track_Theta", &Track_Theta);
	myTTree->Branch("Track_AzimuthAngle", &Track_AzimuthAngle);
	myTTree->Branch("Track_ZenithAngle", &Track_ZenithAngle);
	myTTree->Branch("Track_VertexX", &Track_VertexX);
	myTTree->Branch("Track_VertexY", &Track_VertexY);
	myTTree->Branch("Track_VertexZ", &Track_VertexZ);
	myTTree->Branch("Track_Chi2", &Track_Chi2);
	myTTree->Branch("Track_Chi2PerNdof", &Track_Chi2PerNdof);
	myTTree->Branch("Track_CountValidPoints", &Track_CountValidPoints);
	myTTree->Branch("Track_UnitaryDirectionVectorCoordinateXAtPoint", &Track_UnitaryDirectionVectorCoordinateXAtPoint);
	myTTree->Branch("Track_UnitaryDirectionVectorCoordinateYAtPoint", &Track_UnitaryDirectionVectorCoordinateYAtPoint);
	myTTree->Branch("Track_UnitaryDirectionVectorCoordinateZAtPoint", &Track_UnitaryDirectionVectorCoordinateZAtPoint);
	myTTree->Branch("Track_ThetaAtPoint", &Track_ThetaAtPoint);
	myTTree->Branch("Track_PhiAtPoint", &Track_PhiAtPoint);
	myTTree->Branch("Track_HasMomentum", &Track_HasMomentum);
	myTTree->Branch("Track_EndMomentum", &Track_EndMomentum);
	myTTree->Branch("Track_EndMomentumVectorX", &Track_EndMomentumVectorX);
	myTTree->Branch("Track_EndMomentumVectorY", &Track_EndMomentumVectorY);
	myTTree->Branch("Track_EndMomentumVectorZ", &Track_EndMomentumVectorZ);
	myTTree->Branch("Track_HasPoint", &Track_HasPoint);
	myTTree->Branch("Track_HasValidPoint", &Track_HasValidPoint);
	myTTree->Branch("Track_FirstPoint", &Track_FirstPoint);
	myTTree->Branch("Track_FirstValidPoint", &Track_FirstValidPoint);
	myTTree->Branch("Track_LastPoint", &Track_LastPoint);
	myTTree->Branch("Track_LastValidPoint", &Track_LastValidPoint);
	myTTree->Branch("Track_LocationVectorXAtPoint", &Track_LocationVectorXAtPoint);
	myTTree->Branch("Track_LocationVectorYAtPoint", &Track_LocationVectorYAtPoint);
	myTTree->Branch("Track_LocationVectorZAtPoint", &Track_LocationVectorZAtPoint);
	myTTree->Branch("Track_MomentumAtPoint", &Track_MomentumAtPoint);
	myTTree->Branch("Track_MomentumVectorXAtPoint", &Track_MomentumVectorXAtPoint);
	myTTree->Branch("Track_MomentumVectorYAtPoint", &Track_MomentumVectorYAtPoint);
	myTTree->Branch("Track_MomentumVectorZAtPoint", &Track_MomentumVectorZAtPoint);
	myTTree->Branch("Track_Ndof", &Track_Ndof);
	myTTree->Branch("Track_NextValidPoint", &Track_NextValidPoint);
	myTTree->Branch("Track_PreviousValidPoint", &Track_PreviousValidPoint);
	myTTree->Branch("Track_NPoints", &Track_NPoints);
	myTTree->Branch("Track_NumberTrajectoryPoints", &Track_NumberTrajectoryPoints);
	myTTree->Branch("Track_Calorimetry_Plane0_KineticEnergy", &Track_Calorimetry_Plane0_KineticEnergy);
	myTTree->Branch("Track_Calorimetry_Plane1_KineticEnergy", &Track_Calorimetry_Plane1_KineticEnergy);
	myTTree->Branch("Track_Calorimetry_Plane2_KineticEnergy", &Track_Calorimetry_Plane2_KineticEnergy);
	myTTree->Branch("Track_Momentum", &Track_Momentum);
	myTTree->Branch("Track_Momentum_MCS", &Track_Momentum_MCS);
	myTTree->Branch("Track_Calorimetry_Plane0_ResidualRange", &Track_Calorimetry_Plane0_ResidualRange);
	myTTree->Branch("Track_Calorimetry_Plane1_ResidualRange", &Track_Calorimetry_Plane1_ResidualRange);
	myTTree->Branch("Track_Calorimetry_Plane2_ResidualRange", &Track_Calorimetry_Plane2_ResidualRange);
	myTTree->Branch("Track_Calorimetry_Plane0_dEdx", &Track_Calorimetry_Plane0_dEdx);
	myTTree->Branch("Track_Calorimetry_Plane1_dEdx", &Track_Calorimetry_Plane1_dEdx);
	myTTree->Branch("Track_Calorimetry_Plane2_dEdx", &Track_Calorimetry_Plane2_dEdx);
	myTTree->Branch("Track_Calorimetry_Plane0_dQdx", &Track_Calorimetry_Plane0_dQdx);
	myTTree->Branch("Track_Calorimetry_Plane1_dQdx", &Track_Calorimetry_Plane1_dQdx);
	myTTree->Branch("Track_Calorimetry_Plane2_dQdx", &Track_Calorimetry_Plane2_dQdx);
	myTTree->Branch("Track_Calorimetry_Plane0_TruncdEdx", &Track_Calorimetry_Plane0_TruncdEdx);
	myTTree->Branch("Track_Calorimetry_Plane1_TruncdEdx", &Track_Calorimetry_Plane1_TruncdEdx);
	myTTree->Branch("Track_Calorimetry_Plane2_TruncdEdx", &Track_Calorimetry_Plane2_TruncdEdx);
	myTTree->Branch("Track_Calorimetry_Plane0_TruncdQdx", &Track_Calorimetry_Plane0_TruncdQdx);
	myTTree->Branch("Track_Calorimetry_Plane1_TruncdQdx", &Track_Calorimetry_Plane1_TruncdQdx);
	myTTree->Branch("Track_Calorimetry_Plane2_TruncdQdx", &Track_Calorimetry_Plane2_TruncdQdx);
	myTTree->Branch("Track_MCParticle_E", &Track_MCParticle_E);
	myTTree->Branch("Track_MCParticle_Mass", &Track_MCParticle_Mass);
	myTTree->Branch("Track_MCParticle_Momentum_E", &Track_MCParticle_Momentum_E);
	myTTree->Branch("Track_MCParticle_Momentum_Px", &Track_MCParticle_Momentum_Px);
	myTTree->Branch("Track_MCParticle_Momentum_Py", &Track_MCParticle_Momentum_Py);
	myTTree->Branch("Track_MCParticle_Momentum_Pz", &Track_MCParticle_Momentum_Pz);
	myTTree->Branch("Track_MCParticle_PdgCode", &Track_MCParticle_PdgCode);
	myTTree->Branch("Track_MCParticle_TrackId", &Track_MCParticle_TrackId);
	myTTree->Branch("Track_MCParticle_Purity", &Track_MCParticle_Purity);
	myTTree->Branch("Track_MCParticle_EndX", &Track_MCParticle_EndX);
	myTTree->Branch("Track_MCParticle_EndY", &Track_MCParticle_EndY);
	myTTree->Branch("Track_MCParticle_EndZ", &Track_MCParticle_EndZ);
	myTTree->Branch("Track_MCParticle_Vx", &Track_MCParticle_Vx);
	myTTree->Branch("Track_MCParticle_Vy", &Track_MCParticle_Vy);
	myTTree->Branch("Track_MCParticle_Vz", &Track_MCParticle_Vz);

	//____________________________________________________________________________________________________________________________________________________________________________________________

	myTTree->Branch("NumberHits", &NumberHits, "NumberHits/I");
	myTTree->Branch("Hit_IsBeamEvent", &Hit_IsBeamEvent);
	myTTree->Branch("Hit_IsCosmicEvent", &Hit_IsCosmicEvent);
	myTTree->Branch("Hit_Channel", &Hit_Channel);
	myTTree->Branch("Hit_Integral", &Hit_Integral);
	myTTree->Branch("Hit_Multiplicity", &Hit_Multiplicity);
	myTTree->Branch("Hit_PeakTime", &Hit_PeakTime);
	myTTree->Branch("Hit_PeakAmplitude", &Hit_PeakAmplitude);
	myTTree->Branch("Hit_RMS", &Hit_RMS);
	myTTree->Branch("Hit_Wire", &Hit_Wire);
	myTTree->Branch("Hit_SigmaPeakAmplitude", &Hit_SigmaPeakAmplitude);
	myTTree->Branch("Hit_SigmaPeakTime", &Hit_SigmaPeakTime);
	myTTree->Branch("Hit_WireID_Plane", &Hit_WireID_Plane);
	myTTree->Branch("Hit_MCParticle_PdgCode", &Hit_MCParticle_PdgCode);
	myTTree->Branch("Hit_MCParticle_Purity", &Hit_MCParticle_Purity);
	myTTree->Branch("Hit_MCParticle_EndX", &Hit_MCParticle_EndX);
	myTTree->Branch("Hit_MCParticle_EndY", &Hit_MCParticle_EndY);
	myTTree->Branch("Hit_MCParticle_EndZ", &Hit_MCParticle_EndZ);
	myTTree->Branch("Hit_MCParticle_Vx", &Hit_MCParticle_Vx);
	myTTree->Branch("Hit_MCParticle_Vy", &Hit_MCParticle_Vy);
	myTTree->Branch("Hit_MCParticle_Vz", &Hit_MCParticle_Vz);

	//_____________________________________________________________________________________________________________________________________________________________________________________________

	myTTree->Branch("NumberFlashesBeam", &NumberFlashesBeam, "NumberFlashesBeam/I");
	myTTree->Branch("FlashesBeam_Ywidth", &FlashesBeam_Ywidth);
	myTTree->Branch("FlashesBeam_Zwidth", &FlashesBeam_Zwidth);
	myTTree->Branch("FlashesBeam_Twidth", &FlashesBeam_Twidth);
	myTTree->Branch("FlashesBeam_Ycenter", &FlashesBeam_Ycenter);
	myTTree->Branch("FlashesBeam_Zcenter", &FlashesBeam_Zcenter);
	myTTree->Branch("FlashesBeam_Time", &FlashesBeam_Time);
	myTTree->Branch("FlashesBeam_TotalPE", &FlashesBeam_TotalPE);
	myTTree->Branch("FlashesBeam_PE_Per_PMT", &FlashesBeam_PE_Per_PMT);

	//_____________________________________________________________________________________________________________________________________________________________________________________________

	myTTree->Branch("NumberFlashesCosmic", &NumberFlashesCosmic, "NumberFlashesCosmic/I");
	myTTree->Branch("FlashesCosmic_Ywidth", &FlashesCosmic_Ywidth);
	myTTree->Branch("FlashesCosmic_Zwidth", &FlashesCosmic_Zwidth);
	myTTree->Branch("FlashesCosmic_Twidth", &FlashesCosmic_Twidth);
	myTTree->Branch("FlashesCosmic_Ycenter", &FlashesCosmic_Ycenter);
	myTTree->Branch("FlashesCosmic_Zcenter", &FlashesCosmic_Zcenter);
	myTTree->Branch("FlashesCosmic_Time", &FlashesCosmic_Time);
	myTTree->Branch("FlashesCosmic_TotalPE", &FlashesCosmic_TotalPE);
	myTTree->Branch("FlashesCosmic_PE_Per_PMT", &FlashesCosmic_PE_Per_PMT);

	//_____________________________________________________________________________________________________________________________________________________________________________________________

	myTTree->Branch("NumberOpHitsBeam", &NumberOpHitsBeam, "NumberOpHitsBeam/I");
	myTTree->Branch("OpHitsBeam_Amplitude", &OpHitsBeam_Amplitude);
	myTTree->Branch("OpHitsBeam_Area", &OpHitsBeam_Area);
	myTTree->Branch("OpHitsBeam_OpChannel", &OpHitsBeam_OpChannel);
	myTTree->Branch("OpHitsBeam_PE", &OpHitsBeam_PE);
	myTTree->Branch("OpHitsBeam_PeakTime", &OpHitsBeam_PeakTime);
	myTTree->Branch("OpHitsBeam_Width", &OpHitsBeam_Width);

	//_____________________________________________________________________________________________________________________________________________________________________________________________

	myTTree->Branch("NumberOpHitsCosmic", &NumberOpHitsCosmic, "NumberOpHitsCosmic/I");
	myTTree->Branch("OpHitsCosmic_Amplitude", &OpHitsCosmic_Amplitude);
	myTTree->Branch("OpHitsCosmic_Area", &OpHitsCosmic_Area);
	myTTree->Branch("OpHitsCosmic_OpChannel", &OpHitsCosmic_OpChannel);
	myTTree->Branch("OpHitsCosmic_PE", &OpHitsCosmic_PE);
	myTTree->Branch("OpHitsCosmic_PeakTime", &OpHitsCosmic_PeakTime);
	myTTree->Branch("OpHitsCosmic_Width", &OpHitsCosmic_Width);

	//_____________________________________________________________________________________________________________________________________________________________________________________________

	myTTree->Branch("NumberVertices", &NumberVertices, "NumberVertices/I");
	myTTree->Branch("Vertex_ID", &Vertex_ID);
	myTTree->Branch("Vertex_IsBeamEvent", &Vertex_IsBeamEvent);
	myTTree->Branch("Vertex_IsCosmicEvent", &Vertex_IsCosmicEvent);
	myTTree->Branch("Vertex_PositionX", &Vertex_PositionX);
	myTTree->Branch("Vertex_PositionY", &Vertex_PositionY);
	myTTree->Branch("Vertex_PositionZ", &Vertex_PositionZ);

	//_____________________________________________________________________________________________________________________________________________________________________________________________

	myTTree->Branch("NumberMCTruthEvents", &NumberMCTruthEvents, "NumberMCTruthEvents/I");
	myTTree->Branch("MCTruth_CCNC", &MCTruth_CCNC);
	myTTree->Branch("MCTruth_Mode", &MCTruth_Mode);
	myTTree->Branch("MCTruth_Pt", &MCTruth_Pt);
	myTTree->Branch("MCTruth_QSqr", &MCTruth_QSqr);
	myTTree->Branch("MCTruth_Target", &MCTruth_Target);
	myTTree->Branch("MCTruth_Theta", &MCTruth_Theta);
	myTTree->Branch("MCTruth_X", &MCTruth_X);
	myTTree->Branch("MCTruth_Y", &MCTruth_Y);
	myTTree->Branch("MCTruth_W", &MCTruth_W);
	myTTree->Branch("MCTruth_NParticles", &MCTruth_NParticles);
	myTTree->Branch("MCTruth_Origin", &MCTruth_Origin);
	myTTree->Branch("MCTruth_InteractionType", &MCTruth_InteractionType);
	myTTree->Branch("MCTruth_Particle_Pdg", &MCTruth_Particle_Pdg);
	myTTree->Branch("MCTruth_Particle_TrackId", &MCTruth_Particle_TrackId);
	myTTree->Branch("MCTruth_Particle_EndX", &MCTruth_Particle_EndX);
	myTTree->Branch("MCTruth_Particle_EndY", &MCTruth_Particle_EndY);
	myTTree->Branch("MCTruth_Particle_EndZ", &MCTruth_Particle_EndZ);
	myTTree->Branch("MCTruth_Particle_Vx", &MCTruth_Particle_Vx);
	myTTree->Branch("MCTruth_Particle_Vy", &MCTruth_Particle_Vy);
	myTTree->Branch("MCTruth_Particle_Vz", &MCTruth_Particle_Vz);
	myTTree->Branch("MCTruth_Particle_Nu_Vx", &MCTruth_Particle_Nu_Vx);
	myTTree->Branch("MCTruth_Particle_Nu_Vy", &MCTruth_Particle_Nu_Vy);
	myTTree->Branch("MCTruth_Particle_Nu_Vz", &MCTruth_Particle_Nu_Vz);

	//_____________________________________________________________________________________________________________________________________________________________________________________________

	// PFParticle

	/*myTTree->Branch("NumberPFParticles", &NumberPFParticles, "NumberPFParticles/I");
	myTTree->Branch("PFParticle_IsBeamEvent", &PFParticle_IsBeamEvent);
	myTTree->Branch("PFParticle_IsCosmicEvent", &PFParticle_IsCosmicEvent);
	myTTree->Branch("PFParticle_IsPrimary", &PFParticle_IsPrimary);
	myTTree->Branch("PFParticle_NumDaughters", &PFParticle_NumDaughters);
	myTTree->Branch("PFParticle_Track_StartX", &PFParticle_Track_StartX);
	myTTree->Branch("PFParticle_Track_StartY", &PFParticle_Track_StartY);
	myTTree->Branch("PFParticle_Track_StartZ", &PFParticle_Track_StartZ);
	myTTree->Branch("PFParticle_Track_EndX", &PFParticle_Track_EndX);
	myTTree->Branch("PFParticle_Track_EndY", &PFParticle_Track_EndY);
	myTTree->Branch("PFParticle_Track_EndZ", &PFParticle_Track_EndZ);
	myTTree->Branch("PFParticle_Track_Length", &PFParticle_Track_Length);
	myTTree->Branch("PFParticle_Track_Theta", &PFParticle_Track_Theta);
	myTTree->Branch("PFParticle_Track_Phi", &PFParticle_Track_Phi);*/

	//_____________________________________________________________________________________________________________________________________________________________________________________________

	NEvents = 0;

}
// ____________________________________________________________________________________________________________________________________________________________________________________________________

void mynamespace::TTreeCreator::analyze(art::Event const & e)
{

	NEvents++;
	std::cout << std::endl << std::endl << "Processing " << "Run " << e.run() << ", " << "SubRun " << e.subRun() << ", " << "Event " << e.event() << std::endl;
	std::cout << std::endl << "Processing event # " << NEvents << std::endl << std::endl; 

	//_____________________________________________________________________________________________________________________________________________________________________________________________

	// Handles

	// Tracks

	art::Handle<std::vector<recob::Track>> trk_handle; 
	e.getByLabel("pandoraNu",trk_handle); 
	std::vector<art::Ptr<recob::Track>> trk_vec; 
	art::fill_ptr_vector(trk_vec,trk_handle);

	// Showers

	art::Handle<std::vector<recob::Shower>> shower_handle; 
	e.getByLabel("pandoraNu",shower_handle); 
	std::vector<art::Ptr<recob::Shower>> shower_vec; 
	art::fill_ptr_vector(shower_vec,shower_handle);

	// Hits 

	art::Handle<std::vector<recob::Hit > > hit_handle; 
	e.getByLabel("pandoraCosmicHitRemoval",hit_handle); 
	std::vector<art::Ptr<recob::Hit> > hit_vec; 
	art::fill_ptr_vector(hit_vec,hit_handle);

	// Beam Flashes 

	art::Handle<std::vector<recob::OpFlash>> flash_handle_beam; 
	e.getByLabel("simpleFlashBeam",flash_handle_beam); 
	std::vector<art::Ptr<recob::OpFlash> > flash_beam_vec; 
	art::fill_ptr_vector(flash_beam_vec,flash_handle_beam);

	// Cosmic Flashes

	art::Handle<std::vector<recob::OpFlash>> flash_handle_cosmic; 
	e.getByLabel("simpleFlashCosmic",flash_handle_cosmic); 
	std::vector<art::Ptr<recob::OpFlash> > flash_cosmic_vec; 
	art::fill_ptr_vector(flash_cosmic_vec,flash_handle_cosmic);

	// Beam OpHits

	art::Handle<std::vector<recob::OpHit>> ophit_handle_beam; 
	e.getByLabel("ophitBeam",ophit_handle_beam); 
	std::vector<art::Ptr<recob::OpHit> > ophit_beam_vec; 
	art::fill_ptr_vector(ophit_beam_vec,ophit_handle_beam);

	// Cosmic OpHits

	art::Handle<std::vector<recob::OpHit>> ophit_handle_cosmic; 
	e.getByLabel("ophitCosmic",ophit_handle_cosmic); 
	std::vector<art::Ptr<recob::OpHit> > ophit_cosmic_vec; 
	art::fill_ptr_vector(ophit_cosmic_vec,ophit_handle_cosmic);

	// Vertices 

	art::Handle<std::vector<recob::Vertex>> vertex_handle; 
	e.getByLabel("pandoraNu",vertex_handle); 
	std::vector<art::Ptr<recob::Vertex> > vertex_vec; 
	art::fill_ptr_vector(vertex_vec,vertex_handle);

	// PFParticle 

	art::Handle<std::vector<recob::PFParticle>> pfparticle_handle; 
	e.getByLabel("pandoraNu",pfparticle_handle); 
	std::vector<art::Ptr<recob::PFParticle> > pfparticle_vec; 
	art::fill_ptr_vector(pfparticle_vec,pfparticle_handle);

	// MCParticle

	art::Handle<std::vector<simb::MCParticle>> mcparticle_handle; 
	e.getByLabel("largeant",mcparticle_handle); 
	std::vector<art::Ptr<simb::MCParticle> > largeant_vec; 
	art::fill_ptr_vector(largeant_vec,mcparticle_handle);

	// MCTruth

	art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
	e.getByLabel("generator",mctruthListHandle);
	std::vector<art::Ptr<simb::MCTruth> > mclist;
	art::fill_ptr_vector(mclist, mctruthListHandle);

	// Software Trigger

	art::Handle<raw::ubdaqSoftwareTriggerData> softwareTriggerHandle;
	e.getByLabel("swtrigger", softwareTriggerHandle);

	if (!softwareTriggerHandle.isValid()) {std::cout << "SW Handle doesn't exist!" << std::endl;}
	if (softwareTriggerHandle.isValid()) {
 
		std::cout << "Got the SW Trigger Handle!" << std::endl;
		std::cout << "Number of algorithms = " << softwareTriggerHandle->getNumberOfAlgorithms() << std::endl;	
		std::vector<std::string> algoNames = softwareTriggerHandle->getListOfAlgorithms();	

		for (int i = 0; i < int(algoNames.size()); i++) {

			std::cout << "algoNames["<< i << "] = " << algoNames[i] << std::endl;

			if (algoNames[i] == "BNB_FEMBeamTriggerAlgo") {

				EventPassedSwTrigger = softwareTriggerHandle->passedAlgo(algoNames[i]) ? 1 : 0;

			}

		}

	}

	// Multiple Coulomb Scattering (MCS)

art::ValidHandle<std::vector<recob::MCSFitResult> > MCSMu = e.getValidHandle<std::vector<recob::MCSFitResult> >("pandoraNu");

	// Associations

	art::FindManyP<recob::Hit> hits_per_track(trk_handle,e,"pandoraNu");
	art::FindManyP<recob::Hit> hits_per_shower(shower_handle,e,"pandoraNu");
	art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> mcps_from_hits(hit_handle,e,"crHitRemovalTruthMatch");
	art::FindOneP<recob::PFParticle> VertexToPFParticle(vertex_handle,e,"pandoraNu");
	art::FindOneP<recob::Vertex> PFParticleToVertex(pfparticle_handle,e,"pandoraNu");
	art::FindOneP<recob::Track> PFParticleToTrack(pfparticle_handle,e,"pandoraNu");
	art::FindOneP<recob::Shower> PFParticleToShower(pfparticle_handle,e,"pandoraNu");
	art::FindOneP<simb::MCTruth> MCParticleToMCTruth(mcparticle_handle,e,"largeant");

	//_____________________________________________________________________________________________________________________________________________________________________________________________

	RunNumber = e.run();
	SubRunNumber = e.subRun();
	EventNumber = e.event();

	//_____________________________________________________________________________________________________________________________________________________________________________________________

	// Tracks

	NumberTracks = trk_vec.size();

	// Clear The Vectors

	Track_IsBeamEvent.clear();
	Track_IsCosmicEvent.clear();
	Track_StartX.clear();
	Track_StartY.clear();
	Track_StartZ.clear();
	Track_EndX.clear();
	Track_EndY.clear();
	Track_EndZ.clear();
	Track_ID.clear();
	Track_ParticleId.clear();
	Track_Length.clear();
	Track_Phi.clear();
	Track_AzimuthAngle.clear();
	Track_ZenithAngle.clear();
	Track_Theta.clear();
	Track_VertexX.clear();
	Track_VertexY.clear();
	Track_VertexZ.clear();
	Track_Chi2.clear();
	Track_Chi2PerNdof.clear();
	Track_CountValidPoints.clear();
	Track_HasMomentum.clear();
	Track_EndMomentum.clear();
	Track_EndMomentumVectorX.clear();
	Track_EndMomentumVectorY.clear();
	Track_EndMomentumVectorZ.clear();
	Track_FirstPoint.clear();
	Track_FirstValidPoint.clear();
	Track_LastPoint.clear();
	Track_LastValidPoint.clear();
	Track_Ndof.clear();
	Track_NPoints.clear();
	Track_NumberTrajectoryPoints.clear();

	// Clear The Vectors of Vectors

	Track_UnitaryDirectionVectorCoordinateXAtPoint.clear();
	Track_UnitaryDirectionVectorCoordinateYAtPoint.clear();
	Track_UnitaryDirectionVectorCoordinateZAtPoint.clear();
	Track_ThetaAtPoint.clear();
	Track_PhiAtPoint.clear();
	Track_LocationVectorXAtPoint.clear();
	Track_LocationVectorYAtPoint.clear();
	Track_LocationVectorZAtPoint.clear();
	Track_MomentumAtPoint.clear();
	Track_MomentumVectorXAtPoint.clear();
	Track_MomentumVectorYAtPoint.clear();
	Track_MomentumVectorZAtPoint.clear();
	Track_HasPoint.clear();
	Track_HasValidPoint.clear();
	Track_NextValidPoint.clear();
	Track_PreviousValidPoint.clear();

	Track_Calorimetry_Plane0_KineticEnergy.clear();
	Track_Calorimetry_Plane0_ResidualRange.clear();
	Track_Calorimetry_Plane0_dEdx.clear();
	Track_Calorimetry_Plane0_dQdx.clear();
	Track_Calorimetry_Plane0_TruncdEdx.clear();
	Track_Calorimetry_Plane0_TruncdQdx.clear();
	Track_Calorimetry_Plane1_KineticEnergy.clear();
	Track_Calorimetry_Plane1_ResidualRange.clear();
	Track_Calorimetry_Plane1_dEdx.clear();
	Track_Calorimetry_Plane1_dQdx.clear();
	Track_Calorimetry_Plane1_TruncdEdx.clear();
	Track_Calorimetry_Plane1_TruncdQdx.clear();
	Track_Calorimetry_Plane2_KineticEnergy.clear();
	Track_Momentum.clear();
	Track_Momentum_MCS.clear();
	Track_Calorimetry_Plane2_ResidualRange.clear();
	Track_Calorimetry_Plane2_dEdx.clear();
	Track_Calorimetry_Plane2_dQdx.clear();
	Track_Calorimetry_Plane2_TruncdEdx.clear();
	Track_Calorimetry_Plane2_TruncdQdx.clear();

	Track_MCParticle_E.clear();
	Track_MCParticle_Mass.clear();
	Track_MCParticle_Momentum_E.clear();
	Track_MCParticle_Momentum_Px.clear();
	Track_MCParticle_Momentum_Py.clear();
	Track_MCParticle_Momentum_Pz.clear();
	Track_MCParticle_PdgCode.clear();
	Track_MCParticle_TrackId.clear();
	Track_MCParticle_Purity.clear();
	Track_MCParticle_EndX.clear();
	Track_MCParticle_EndY.clear();
	Track_MCParticle_EndZ.clear();
	Track_MCParticle_Vx.clear();
	Track_MCParticle_Vy.clear();
	Track_MCParticle_Vz.clear();

std::cout << std::endl << "Loop over the tracks" << std::endl << std::endl; 

	for (int i_t = 0; i_t < int(trk_vec.size()); ++i_t) {

		art::Ptr<recob::Track> CurrentTrack = trk_vec.at(i_t);

		Track_StartX.push_back(CurrentTrack->Start().X());
		Track_StartY.push_back(CurrentTrack->Start().Y());
		Track_StartZ.push_back(CurrentTrack->Start().Z());
		Track_EndX.push_back(CurrentTrack->End()[0]);
		Track_EndY.push_back(CurrentTrack->End()[1]);
		Track_EndZ.push_back(CurrentTrack->End()[2]);
		Track_ID.push_back(CurrentTrack->ID());
		Track_ParticleId.push_back(CurrentTrack->ParticleId());
		Track_Length.push_back(CurrentTrack->Length(i_t));
		Track_Phi.push_back(CurrentTrack->Phi());
		Track_AzimuthAngle.push_back(CurrentTrack->AzimuthAngle());
		Track_ZenithAngle.push_back(CurrentTrack->ZenithAngle());
		Track_Theta.push_back(CurrentTrack->Theta());
		Track_VertexX.push_back(CurrentTrack->Vertex().X());
		Track_VertexY.push_back(CurrentTrack->Vertex().Y());
		Track_VertexZ.push_back(CurrentTrack->Vertex().Z());
		Track_Chi2.push_back(CurrentTrack->Chi2());
		Track_Chi2PerNdof.push_back(CurrentTrack->Chi2PerNdof());
		Track_CountValidPoints.push_back(CurrentTrack->CountValidPoints());
		Track_HasMomentum.push_back(CurrentTrack->HasMomentum());
		Track_EndMomentum.push_back(CurrentTrack->EndMomentum());
		Track_EndMomentumVectorX.push_back(CurrentTrack->EndMomentumVector().X());
		Track_EndMomentumVectorY.push_back(CurrentTrack->EndMomentumVector().Y());
		Track_EndMomentumVectorZ.push_back(CurrentTrack->EndMomentumVector().Z());
		Track_FirstPoint.push_back(CurrentTrack->FirstPoint());
		Track_FirstValidPoint.push_back(CurrentTrack->FirstValidPoint());
		Track_LastPoint.push_back(CurrentTrack->LastPoint());
		Track_LastValidPoint.push_back(CurrentTrack->LastValidPoint());
		Track_Ndof.push_back(CurrentTrack->Ndof());
		Track_NPoints.push_back(CurrentTrack->NPoints());
		Track_NumberTrajectoryPoints.push_back(CurrentTrack->NumberTrajectoryPoints());

		int NPoints = CurrentTrack->NPoints();

		// Declare The Subvectors

		std::vector<double> CurrentTrack_UnitaryDirectionVectorCoordinateXAtPoint;
		std::vector<double> CurrentTrack_UnitaryDirectionVectorCoordinateYAtPoint;
		std::vector<double> CurrentTrack_UnitaryDirectionVectorCoordinateZAtPoint;
		std::vector<double> CurrentTrack_ThetaAtPoint;
		std::vector<double> CurrentTrack_PhiAtPoint;
		std::vector<int> CurrentTrack_HasPoint;
		std::vector<int> CurrentTrack_HasValidPoint;
		std::vector<double> CurrentTrack_LocationVectorXAtPoint;
		std::vector<double> CurrentTrack_LocationVectorYAtPoint;
		std::vector<double> CurrentTrack_LocationVectorZAtPoint;
		std::vector<double> CurrentTrack_MomentumAtPoint;
		std::vector<double> CurrentTrack_MomentumVectorXAtPoint;
		std::vector<double> CurrentTrack_MomentumVectorYAtPoint;
		std::vector<double> CurrentTrack_MomentumVectorZAtPoint;
		std::vector<int> CurrentTrack_NextValidPoint;
		std::vector<int> CurrentTrack_PreviousValidPoint;

		// Clear The Subvectors

		CurrentTrack_UnitaryDirectionVectorCoordinateXAtPoint.clear();
		CurrentTrack_UnitaryDirectionVectorCoordinateXAtPoint.clear();
		CurrentTrack_UnitaryDirectionVectorCoordinateYAtPoint.clear();
		CurrentTrack_UnitaryDirectionVectorCoordinateZAtPoint.clear();
		CurrentTrack_ThetaAtPoint.clear();
		CurrentTrack_PhiAtPoint.clear();
		CurrentTrack_LocationVectorXAtPoint.clear();
		CurrentTrack_LocationVectorYAtPoint.clear();
		CurrentTrack_LocationVectorZAtPoint.clear();
		CurrentTrack_MomentumAtPoint.clear();
		CurrentTrack_MomentumVectorXAtPoint.clear();
		CurrentTrack_MomentumVectorYAtPoint.clear();
		CurrentTrack_MomentumVectorZAtPoint.clear();
		CurrentTrack_HasPoint.clear();
		CurrentTrack_HasValidPoint.clear();
		CurrentTrack_NextValidPoint.clear();
		CurrentTrack_PreviousValidPoint.clear();

		for (int i_Point = 0; i_Point < NPoints; i_Point++) {

			TVector3 CurrentTrack_UnitaryDirectionVectorComponentsAt_i_Point = CurrentTrack->DirectionAtPoint(i_Point);

			CurrentTrack_UnitaryDirectionVectorCoordinateXAtPoint.push_back(CurrentTrack_UnitaryDirectionVectorComponentsAt_i_Point.X());
			CurrentTrack_UnitaryDirectionVectorCoordinateYAtPoint.push_back(CurrentTrack_UnitaryDirectionVectorComponentsAt_i_Point.Y());
			CurrentTrack_UnitaryDirectionVectorCoordinateZAtPoint.push_back(CurrentTrack_UnitaryDirectionVectorComponentsAt_i_Point.Z());
			CurrentTrack_ThetaAtPoint.push_back(CurrentTrack_UnitaryDirectionVectorComponentsAt_i_Point.Theta());
			CurrentTrack_PhiAtPoint.push_back(CurrentTrack_UnitaryDirectionVectorComponentsAt_i_Point.Phi());
			CurrentTrack_HasPoint.push_back(CurrentTrack->HasPoint(i_Point));
			CurrentTrack_HasValidPoint.push_back(CurrentTrack->HasValidPoint(i_Point));
			CurrentTrack_LocationVectorXAtPoint.push_back(CurrentTrack->LocationAtPoint(i_Point).X());
			CurrentTrack_LocationVectorYAtPoint.push_back(CurrentTrack->LocationAtPoint(i_Point).Y());
			CurrentTrack_LocationVectorZAtPoint.push_back(CurrentTrack->LocationAtPoint(i_Point).Z());
			CurrentTrack_MomentumAtPoint.push_back(CurrentTrack->MomentumAtPoint(i_Point));
			CurrentTrack_MomentumVectorXAtPoint.push_back(CurrentTrack->MomentumVectorAtPoint(i_Point).X());
			CurrentTrack_MomentumVectorYAtPoint.push_back(CurrentTrack->MomentumVectorAtPoint(i_Point).Y());
			CurrentTrack_MomentumVectorZAtPoint.push_back(CurrentTrack->MomentumVectorAtPoint(i_Point).Z());
			CurrentTrack_NextValidPoint.push_back(CurrentTrack->NextValidPoint(i_Point));
			CurrentTrack_PreviousValidPoint.push_back(CurrentTrack->PreviousValidPoint(i_Point));

		}

		Track_UnitaryDirectionVectorCoordinateXAtPoint.push_back(CurrentTrack_UnitaryDirectionVectorCoordinateXAtPoint);
		Track_UnitaryDirectionVectorCoordinateYAtPoint.push_back(CurrentTrack_UnitaryDirectionVectorCoordinateYAtPoint);
		Track_UnitaryDirectionVectorCoordinateZAtPoint.push_back(CurrentTrack_UnitaryDirectionVectorCoordinateZAtPoint);
		Track_ThetaAtPoint.push_back(CurrentTrack_ThetaAtPoint);
		Track_PhiAtPoint.push_back(CurrentTrack_PhiAtPoint);
		Track_HasPoint.push_back(CurrentTrack_HasPoint);
		Track_HasValidPoint.push_back(CurrentTrack_HasValidPoint);
		Track_LocationVectorXAtPoint.push_back(CurrentTrack_LocationVectorXAtPoint);
		Track_LocationVectorYAtPoint.push_back(CurrentTrack_LocationVectorYAtPoint);
		Track_LocationVectorZAtPoint.push_back(CurrentTrack_LocationVectorZAtPoint);
		Track_MomentumAtPoint.push_back(CurrentTrack_MomentumAtPoint);
		Track_MomentumVectorXAtPoint.push_back(CurrentTrack_MomentumVectorXAtPoint);
		Track_MomentumVectorYAtPoint.push_back(CurrentTrack_MomentumVectorYAtPoint);
		Track_MomentumVectorZAtPoint.push_back(CurrentTrack_MomentumVectorZAtPoint);
		Track_NextValidPoint.push_back(CurrentTrack_NextValidPoint);
		Track_PreviousValidPoint.push_back(CurrentTrack_PreviousValidPoint);

		//____________________________________________________________________________________________________________________________________________________________________________________

		// Multiple Coulomb Scattering

const recob::MCSFitResult& mcsMu = MCSMu->at(CurrentTrack);
double trkMom_MuFwd = mcsMu.fwdMomentum();
Track_Momentum_MCS.push_back(trkMom_MuFwd);

		//____________________________________________________________________________________________________________________________________________________________________________________

		// Calorimetry
std::cout << std::endl << "Doing Calorimetry" << std::endl << std::endl;
		art::FindManyP<anab::Calorimetry> CalorimetryFromTracks(trk_handle,e,"pandoraNucali");
		std::vector<art::Ptr<anab::Calorimetry>> Calorimetry_vec = CalorimetryFromTracks.at(i_t);

		TruncMean truncmean;

		// Plane 0

		art::Ptr<anab::Calorimetry> plane0 = Calorimetry_vec.at(0);
		Track_Calorimetry_Plane0_KineticEnergy.push_back(plane0->KineticEnergy());
		std::vector< double > dedxPlane0 = plane0->dEdx();
		std::vector< double > dqdxPlane0 = plane0->dQdx();
		std::vector< double > ResidualRangePlane0 = plane0->ResidualRange();
		int NHits_Plane0 = dedxPlane0.size();

		std::vector<double> CurrentTrack_Calorimetry_Plane0_ResidualRange;
		std::vector<double> CurrentTrack_Calorimetry_Plane0_dEdx;
		std::vector<double> CurrentTrack_Calorimetry_Plane0_dQdx;
		std::vector<double> CurrentTrack_Calorimetry_Plane0_TruncdEdx;
		std::vector<double> CurrentTrack_Calorimetry_Plane0_TruncdQdx;

		CurrentTrack_Calorimetry_Plane0_ResidualRange.clear();
		CurrentTrack_Calorimetry_Plane0_dEdx.clear();
		CurrentTrack_Calorimetry_Plane0_dQdx.clear();
		CurrentTrack_Calorimetry_Plane0_TruncdEdx.clear();
		CurrentTrack_Calorimetry_Plane0_TruncdQdx.clear();

		for (int i_Hit_Plane0 = 0; i_Hit_Plane0 < NHits_Plane0; i_Hit_Plane0++) {

			CurrentTrack_Calorimetry_Plane0_ResidualRange.push_back(ResidualRangePlane0.at(i_Hit_Plane0));
			CurrentTrack_Calorimetry_Plane0_dEdx.push_back(dedxPlane0.at(i_Hit_Plane0));
			CurrentTrack_Calorimetry_Plane0_dQdx.push_back(dqdxPlane0.at(i_Hit_Plane0));

		}

		Track_Calorimetry_Plane0_ResidualRange.push_back(CurrentTrack_Calorimetry_Plane0_ResidualRange);
		Track_Calorimetry_Plane0_dEdx.push_back(CurrentTrack_Calorimetry_Plane0_dEdx);
		Track_Calorimetry_Plane0_dQdx.push_back(CurrentTrack_Calorimetry_Plane0_dQdx);

		truncmean.CalcTruncMean(CurrentTrack_Calorimetry_Plane0_ResidualRange,CurrentTrack_Calorimetry_Plane0_dEdx,CurrentTrack_Calorimetry_Plane0_TruncdEdx);
		truncmean.CalcTruncMean(CurrentTrack_Calorimetry_Plane0_ResidualRange,CurrentTrack_Calorimetry_Plane0_dQdx,CurrentTrack_Calorimetry_Plane0_TruncdQdx);

		Track_Calorimetry_Plane0_TruncdEdx.push_back(CurrentTrack_Calorimetry_Plane0_TruncdEdx);
		Track_Calorimetry_Plane0_TruncdQdx.push_back(CurrentTrack_Calorimetry_Plane0_TruncdQdx);

		// Plane 1

		art::Ptr<anab::Calorimetry> plane1 = Calorimetry_vec.at(1);
		Track_Calorimetry_Plane1_KineticEnergy.push_back(plane1->KineticEnergy());
		std::vector< double > dedxPlane1 = plane1->dEdx();
		std::vector< double > dqdxPlane1 = plane1->dQdx();
		std::vector< double > ResidualRangePlane1 = plane1->ResidualRange();
		int NHits_Plane1 = dedxPlane1.size();

		std::vector<double> CurrentTrack_Calorimetry_Plane1_ResidualRange;
		std::vector<double> CurrentTrack_Calorimetry_Plane1_dEdx;
		std::vector<double> CurrentTrack_Calorimetry_Plane1_dQdx;
		std::vector<double> CurrentTrack_Calorimetry_Plane1_TruncdEdx;
		std::vector<double> CurrentTrack_Calorimetry_Plane1_TruncdQdx;

		CurrentTrack_Calorimetry_Plane1_ResidualRange.clear();
		CurrentTrack_Calorimetry_Plane1_dEdx.clear();
		CurrentTrack_Calorimetry_Plane1_dQdx.clear();
		CurrentTrack_Calorimetry_Plane1_TruncdEdx.clear();
		CurrentTrack_Calorimetry_Plane1_TruncdQdx.clear();

		for (int i_Hit_Plane1 = 0; i_Hit_Plane1 < NHits_Plane1; i_Hit_Plane1++) {

			CurrentTrack_Calorimetry_Plane1_ResidualRange.push_back(ResidualRangePlane1.at(i_Hit_Plane1));
			CurrentTrack_Calorimetry_Plane1_dEdx.push_back(dedxPlane1.at(i_Hit_Plane1));
			CurrentTrack_Calorimetry_Plane1_dQdx.push_back(dqdxPlane1.at(i_Hit_Plane1));		

		}

		Track_Calorimetry_Plane1_ResidualRange.push_back(CurrentTrack_Calorimetry_Plane1_ResidualRange);
		Track_Calorimetry_Plane1_dEdx.push_back(CurrentTrack_Calorimetry_Plane1_dEdx);
		Track_Calorimetry_Plane1_dQdx.push_back(CurrentTrack_Calorimetry_Plane1_dQdx);

		truncmean.CalcTruncMean(CurrentTrack_Calorimetry_Plane1_ResidualRange,CurrentTrack_Calorimetry_Plane1_dEdx,CurrentTrack_Calorimetry_Plane1_TruncdEdx);
		truncmean.CalcTruncMean(CurrentTrack_Calorimetry_Plane1_ResidualRange,CurrentTrack_Calorimetry_Plane1_dQdx,CurrentTrack_Calorimetry_Plane1_TruncdQdx);

		Track_Calorimetry_Plane1_TruncdEdx.push_back(CurrentTrack_Calorimetry_Plane1_TruncdEdx);
		Track_Calorimetry_Plane1_TruncdQdx.push_back(CurrentTrack_Calorimetry_Plane1_TruncdQdx);

		// Plane 2

		art::Ptr<anab::Calorimetry> plane2 = Calorimetry_vec.at(2);
		Track_Calorimetry_Plane2_KineticEnergy.push_back(plane2->KineticEnergy());	
		std::vector< double > dedxPlane2 = plane2->dEdx();
		std::vector< double > dqdxPlane2 = plane2->dQdx();
		std::vector< double > ResidualRangePlane2 = plane2->ResidualRange();
		int NHits_Plane2 = dedxPlane2.size();

		std::vector<double> CurrentTrack_Calorimetry_Plane2_ResidualRange;
		std::vector<double> CurrentTrack_Calorimetry_Plane2_dEdx;
		std::vector<double> CurrentTrack_Calorimetry_Plane2_dQdx;
		std::vector<double> CurrentTrack_Calorimetry_Plane2_TruncdEdx;
		std::vector<double> CurrentTrack_Calorimetry_Plane2_TruncdQdx;

		CurrentTrack_Calorimetry_Plane2_ResidualRange.clear();
		CurrentTrack_Calorimetry_Plane2_dEdx.clear();
		CurrentTrack_Calorimetry_Plane2_dQdx.clear();
		CurrentTrack_Calorimetry_Plane2_TruncdEdx.clear();
		CurrentTrack_Calorimetry_Plane2_TruncdQdx.clear();

		for (int i_Hit_Plane2 = 0; i_Hit_Plane2 < NHits_Plane2; i_Hit_Plane2++) {

			CurrentTrack_Calorimetry_Plane2_ResidualRange.push_back(ResidualRangePlane2.at(i_Hit_Plane2));
			CurrentTrack_Calorimetry_Plane2_dEdx.push_back(dedxPlane2.at(i_Hit_Plane2));
			CurrentTrack_Calorimetry_Plane2_dQdx.push_back(dqdxPlane2.at(i_Hit_Plane2));		

		}

		Track_Calorimetry_Plane2_ResidualRange.push_back(CurrentTrack_Calorimetry_Plane2_ResidualRange);
		Track_Calorimetry_Plane2_dEdx.push_back(CurrentTrack_Calorimetry_Plane2_dEdx);
		Track_Calorimetry_Plane2_dQdx.push_back(CurrentTrack_Calorimetry_Plane2_dQdx);

		truncmean.CalcTruncMean(CurrentTrack_Calorimetry_Plane2_ResidualRange,CurrentTrack_Calorimetry_Plane2_dEdx,CurrentTrack_Calorimetry_Plane2_TruncdEdx);
		truncmean.CalcTruncMean(CurrentTrack_Calorimetry_Plane2_ResidualRange,CurrentTrack_Calorimetry_Plane2_dQdx,CurrentTrack_Calorimetry_Plane2_TruncdQdx);

		Track_Calorimetry_Plane2_TruncdEdx.push_back(CurrentTrack_Calorimetry_Plane2_TruncdEdx);
		Track_Calorimetry_Plane2_TruncdQdx.push_back(CurrentTrack_Calorimetry_Plane2_TruncdQdx);

// ____________________________________________________________________________________________________________________________________________________________________________________________________

		// BackTracker

		std::vector< art::Ptr<recob::Hit> > trk_hits_ptrs = hits_per_track.at(i_t);
		BackTrackerTruthMatch backtrackertruthmatch; 
		backtrackertruthmatch.MatchToMCParticle(hit_handle,e,trk_hits_ptrs);
		art::Ptr< simb::MCParticle > maxp_me = backtrackertruthmatch.ReturnMCParticle();
std::cout << "Hello, I'm going to call the BackTracker!" << std::endl;

		if ( maxp_me.isNull() ) { 

			std::cout << "Event not matched by the BackTracker" << std::endl;
			Track_IsBeamEvent.push_back(0);
			Track_IsCosmicEvent.push_back(1);

			Track_Momentum.push_back(-99.);

			Track_MCParticle_E.push_back(-99);
			Track_MCParticle_Mass.push_back(-99);
			Track_MCParticle_Momentum_E.push_back(-99);
			Track_MCParticle_Momentum_Px.push_back(-99);
			Track_MCParticle_Momentum_Py.push_back(-99);
			Track_MCParticle_Momentum_Pz.push_back(-99);
			Track_MCParticle_PdgCode.push_back(-99);
			Track_MCParticle_TrackId.push_back(-99);
			Track_MCParticle_Purity.push_back(-99.);
			Track_MCParticle_EndX.push_back(-99.);
			Track_MCParticle_EndY.push_back(-99.);
			Track_MCParticle_EndZ.push_back(-99.);
			Track_MCParticle_Vx.push_back(-99.);
			Track_MCParticle_Vy.push_back(-99.);
			Track_MCParticle_Vz.push_back(-99.);

		} else {

			double purity = backtrackertruthmatch.ReturnPurity();

			std::cout << "Event matched by the BackTracker" << std::endl;
			const art::Ptr<simb::MCTruth> mctruth = MCParticleToMCTruth.at(maxp_me.key());


			// Beginning of the matched beam selection

			if (mctruth->Origin() == 1) { 

				std::cout << "Beam Event" << std::endl;
				Track_IsBeamEvent.push_back(1);
				Track_IsCosmicEvent.push_back(0);

				trkf::TrackMomentumCalculator TrackMomCalc;
				double TrackMomentum =  TrackMomCalc.GetTrackMomentum(CurrentTrack->Length(i_t),maxp_me->PdgCode()); 
				Track_Momentum.push_back(TrackMomentum);

				Track_MCParticle_E.push_back(maxp_me->E());
				Track_MCParticle_Mass.push_back(maxp_me->Mass());
				Track_MCParticle_Momentum_E.push_back(maxp_me->Momentum().E());
				Track_MCParticle_Momentum_Px.push_back(maxp_me->Momentum().Px());
				Track_MCParticle_Momentum_Py.push_back(maxp_me->Momentum().Py());
				Track_MCParticle_Momentum_Pz.push_back(maxp_me->Momentum().Pz());
				Track_MCParticle_PdgCode.push_back(maxp_me->PdgCode());
				Track_MCParticle_TrackId.push_back(maxp_me->TrackId());
				Track_MCParticle_Purity.push_back(purity);
				Track_MCParticle_EndX.push_back(maxp_me->EndX());
				Track_MCParticle_EndY.push_back(maxp_me->EndY());
				Track_MCParticle_EndZ.push_back(maxp_me->EndZ());
				Track_MCParticle_Vx.push_back(maxp_me->Vx());
				Track_MCParticle_Vy.push_back(maxp_me->Vy());
				Track_MCParticle_Vz.push_back(maxp_me->Vz());


			} // End of the matched beam selection


			// Beginning of the unmatched / cosmic selection

			if (mctruth->Origin() == 2) { 

				std::cout << "Cosmic Event" << std::endl;
				Track_IsBeamEvent.push_back(0);
				Track_IsCosmicEvent.push_back(1);

				trkf::TrackMomentumCalculator TrackMomCalc;
				double TrackMomentum =  TrackMomCalc.GetTrackMomentum(CurrentTrack->Length(i_t),maxp_me->PdgCode()); 
				Track_Momentum.push_back(TrackMomentum);

				Track_MCParticle_E.push_back(maxp_me->E());
				Track_MCParticle_Mass.push_back(maxp_me->Mass());
				Track_MCParticle_Momentum_E.push_back(maxp_me->Momentum().E());
				Track_MCParticle_Momentum_Px.push_back(maxp_me->Momentum().Px());
				Track_MCParticle_Momentum_Py.push_back(maxp_me->Momentum().Py());
				Track_MCParticle_Momentum_Pz.push_back(maxp_me->Momentum().Pz());
				Track_MCParticle_PdgCode.push_back(maxp_me->PdgCode());
				Track_MCParticle_TrackId.push_back(maxp_me->TrackId());
				Track_MCParticle_Purity.push_back(purity);
				Track_MCParticle_EndX.push_back(maxp_me->EndX());
				Track_MCParticle_EndY.push_back(maxp_me->EndY());
				Track_MCParticle_EndZ.push_back(maxp_me->EndZ());
				Track_MCParticle_Vx.push_back(maxp_me->Vx());
				Track_MCParticle_Vy.push_back(maxp_me->Vy());
				Track_MCParticle_Vz.push_back(maxp_me->Vz());

			} // End of the unmatched / cosmic selection

		} // End of the if-else statement for the BackTracker

	} // End of the loop over the tracks
	
	//_____________________________________________________________________________________________________________________________________________________________________________________________

	// Hits 
std::cout << std::endl << "Loop over the hits" << std::endl << std::endl;
	NumberHits = hit_vec.size();
	Hit_IsBeamEvent.clear();
	Hit_IsCosmicEvent.clear();
	Hit_Channel.clear();
	Hit_Integral.clear();
	Hit_Multiplicity.clear();
	Hit_PeakTime.clear();
	Hit_PeakAmplitude.clear();
	Hit_RMS.clear();
	Hit_Wire.clear();
	Hit_SigmaPeakAmplitude.clear();
	Hit_SigmaPeakTime.clear();
	Hit_WireID_Plane.clear();

	Hit_MCParticle_PdgCode.clear();
	Hit_MCParticle_Purity.clear();
	Hit_MCParticle_EndX.clear();
	Hit_MCParticle_EndY.clear();
	Hit_MCParticle_EndZ.clear();
	Hit_MCParticle_Vx.clear();
	Hit_MCParticle_Vy.clear();
	Hit_MCParticle_Vz.clear();

	for (int i_h = 0; i_h < int(hit_vec.size()); ++i_h) {

		art::Ptr<recob::Hit> CurrentHit = hit_vec.at(i_h);

		Hit_Channel.push_back(CurrentHit->Channel());
		Hit_Integral.push_back(CurrentHit->Integral());
		Hit_Multiplicity.push_back(CurrentHit->Multiplicity());
		Hit_PeakTime.push_back(CurrentHit->PeakTime());
		Hit_PeakAmplitude.push_back(CurrentHit->PeakAmplitude());
		Hit_RMS.push_back(CurrentHit->RMS());
		Hit_Wire.push_back(CurrentHit->WireID().Wire);
		Hit_SigmaPeakAmplitude.push_back(CurrentHit->SigmaPeakAmplitude());
		Hit_SigmaPeakTime.push_back(CurrentHit->SigmaPeakTime());
		Hit_WireID_Plane.push_back(CurrentHit->WireID().Plane);

		// BackTracker

		BackTrackerTruthMatch backtrackertruthmatchhit;
		std::vector<art::Ptr<recob::Hit> > hit_vec_ptr;
		hit_vec_ptr.push_back(CurrentHit);
		backtrackertruthmatchhit.MatchToMCParticle(hit_handle,e,hit_vec_ptr);

		art::Ptr< simb::MCParticle > maxp_me_hit = backtrackertruthmatchhit.ReturnMCParticle();

		if ( maxp_me_hit.isNull()) { 

			//std::cout << "Event not matched by the BackTracker" << std::endl;
			Hit_IsBeamEvent.push_back(0);
			Hit_IsCosmicEvent.push_back(1);

			Hit_MCParticle_PdgCode.push_back(-99);
			Hit_MCParticle_EndX.push_back(-99.);
			Hit_MCParticle_EndY.push_back(-99.);
			Hit_MCParticle_EndZ.push_back(-99.);
			Hit_MCParticle_Vx.push_back(-99.);
			Hit_MCParticle_Vy.push_back(-99.);
			Hit_MCParticle_Vz.push_back(-99.);

		} else {

			//std::cout << "Event matched by the BackTracker" << std::endl;
			const art::Ptr<simb::MCTruth> mctruth_hit = MCParticleToMCTruth.at(maxp_me_hit.key());


			// Beginning of the matched beam selection

			if (mctruth_hit->Origin() == 1) { 

				//std::cout << "Beam Event" << std::endl;
				Hit_IsBeamEvent.push_back(1);
				Hit_IsCosmicEvent.push_back(0);

				Hit_MCParticle_PdgCode.push_back(maxp_me_hit->PdgCode());
				Hit_MCParticle_EndX.push_back(maxp_me_hit->EndX());
				Hit_MCParticle_EndY.push_back(maxp_me_hit->EndY());
				Hit_MCParticle_EndZ.push_back(maxp_me_hit->EndZ());
				Hit_MCParticle_Vx.push_back(maxp_me_hit->Vx());
				Hit_MCParticle_Vy.push_back(maxp_me_hit->Vy());
				Hit_MCParticle_Vz.push_back(maxp_me_hit->Vz());


			} // End of the matched beam selection


			// Beginning of the unmatched / cosmic selection

			if (mctruth_hit->Origin() == 2) { 

				//std::cout << "Cosmic Event" << std::endl;
				Hit_IsBeamEvent.push_back(0);
				Hit_IsCosmicEvent.push_back(1);

				Hit_MCParticle_PdgCode.push_back(maxp_me_hit->PdgCode());
				Hit_MCParticle_EndX.push_back(maxp_me_hit->EndX());
				Hit_MCParticle_EndY.push_back(maxp_me_hit->EndY());
				Hit_MCParticle_EndZ.push_back(maxp_me_hit->EndZ());
				Hit_MCParticle_Vx.push_back(maxp_me_hit->Vx());
				Hit_MCParticle_Vy.push_back(maxp_me_hit->Vy());
				Hit_MCParticle_Vz.push_back(maxp_me_hit->Vz());

			} // End of the unmatched / cosmic selection

		} // End of the if-else statement for the BackTracker

	} // End of the loop over the hits
// ____________________________________________________________________________________________________________________________________________________________________________________________________

	// Beam Flashes 
std::cout << std::endl << "Loop over the beam flashes" << std::endl << std::endl;
	NumberFlashesBeam = flash_beam_vec.size();

	FlashesBeam_Ywidth.clear();
	FlashesBeam_Zwidth.clear();
	FlashesBeam_Twidth.clear();
	FlashesBeam_Ycenter.clear();
	FlashesBeam_Zcenter.clear();
	FlashesBeam_Time.clear();
	FlashesBeam_TotalPE.clear();
	FlashesBeam_PE_Per_PMT.clear();

	for (int i_fb = 0; i_fb < int(flash_beam_vec.size()); ++i_fb) {

		art::Ptr<recob::OpFlash> CurrentBeamFlash = flash_beam_vec.at(i_fb);

		FlashesBeam_Ywidth.push_back(CurrentBeamFlash->YWidth());
		FlashesBeam_Zwidth.push_back(CurrentBeamFlash->ZWidth());
		FlashesBeam_Twidth.push_back(CurrentBeamFlash->TimeWidth());
		FlashesBeam_Ycenter.push_back(CurrentBeamFlash->YCenter());
		FlashesBeam_Zcenter.push_back(CurrentBeamFlash->ZCenter());
		FlashesBeam_Time.push_back(CurrentBeamFlash->Time());
		FlashesBeam_TotalPE.push_back(CurrentBeamFlash->TotalPE());

		std::vector<double> CurrentBeamFlash_PE_Per_PMT;

		CurrentBeamFlash_PE_Per_PMT.clear();

		for(int i_pmt = 0; i_pmt < 32; i_pmt++) {

			const double optical_channel_flashes_beam = CurrentBeamFlash->PE(i_pmt);

			CurrentBeamFlash_PE_Per_PMT.push_back(optical_channel_flashes_beam);

		}

		FlashesBeam_PE_Per_PMT.push_back(CurrentBeamFlash_PE_Per_PMT);

	} // End of the loop over the beam flashes
// ____________________________________________________________________________________________________________________________________________________________________________________________________

	// Cosmic Flashes
std::cout << std::endl << "Loop over the cosmic flashes" << std::endl << std::endl;
	NumberFlashesCosmic = flash_cosmic_vec.size();
	FlashesCosmic_Ywidth.clear();
	FlashesCosmic_Zwidth.clear();
	FlashesCosmic_Twidth.clear();
	FlashesCosmic_Ycenter.clear();
	FlashesCosmic_Zcenter.clear();
	FlashesCosmic_Time.clear();
	FlashesCosmic_TotalPE.clear();

	for (int i_fc = 0; i_fc < int(flash_cosmic_vec.size()); ++i_fc) {

		art::Ptr<recob::OpFlash> CurrentCosmicFlash = flash_cosmic_vec.at(i_fc);

		FlashesCosmic_Ywidth.push_back(CurrentCosmicFlash->YWidth());
		FlashesCosmic_Zwidth.push_back(CurrentCosmicFlash->ZWidth());
		FlashesCosmic_Twidth.push_back(CurrentCosmicFlash->TimeWidth());
		FlashesCosmic_Ycenter.push_back(CurrentCosmicFlash->YCenter());
		FlashesCosmic_Zcenter.push_back(CurrentCosmicFlash->ZCenter());
		FlashesCosmic_Time.push_back(CurrentCosmicFlash->Time());
		FlashesCosmic_TotalPE.push_back(CurrentCosmicFlash->TotalPE());

		std::vector<double> CurrentCosmicFlash_PE_Per_PMT;

		CurrentCosmicFlash_PE_Per_PMT.clear();

		for(int i_pmt = 200; i_pmt < 232; i_pmt++) {

			const double optical_channel_flashes_cosmic = CurrentCosmicFlash->PE(i_pmt);

			CurrentCosmicFlash_PE_Per_PMT.push_back(optical_channel_flashes_cosmic);

		}

		FlashesCosmic_PE_Per_PMT.push_back(CurrentCosmicFlash_PE_Per_PMT);

	} // End of the loop over the cosmic flashes

// ____________________________________________________________________________________________________________________________________________________________________________________________________

	// Beam OpHits 
std::cout << std::endl << "Loop over the beam optical hits" << std::endl << std::endl;
	NumberOpHitsBeam = ophit_beam_vec.size();
	OpHitsBeam_Amplitude.clear();
	OpHitsBeam_Area.clear();
	OpHitsBeam_OpChannel.clear();
	OpHitsBeam_PE.clear();
	OpHitsBeam_PeakTime.clear();
	OpHitsBeam_Width.clear();

	for (int i_ohb = 0; i_ohb < int(ophit_beam_vec.size()); ++i_ohb) {

		art::Ptr<recob::OpHit> CurrentBeamOpHit = ophit_beam_vec.at(i_ohb);

		OpHitsBeam_Amplitude.push_back(CurrentBeamOpHit->Amplitude());
		OpHitsBeam_Area.push_back(CurrentBeamOpHit->Area());
		OpHitsBeam_OpChannel.push_back(CurrentBeamOpHit->OpChannel());
		OpHitsBeam_PE.push_back(CurrentBeamOpHit->PE());
		OpHitsBeam_PeakTime.push_back(CurrentBeamOpHit->PeakTime());
		OpHitsBeam_Width.push_back(CurrentBeamOpHit->Width());

	} // End of the loop over the beam optical hits

// ____________________________________________________________________________________________________________________________________________________________________________________________________

	// Cosmic OpHits 
std::cout << std::endl << "Loop over the cosmic optical hits" << std::endl << std::endl;
	NumberOpHitsCosmic = ophit_cosmic_vec.size();
	OpHitsCosmic_Amplitude.clear();
	OpHitsCosmic_Area.clear();
	OpHitsCosmic_OpChannel.clear();
	OpHitsCosmic_PE.clear();
	OpHitsCosmic_PeakTime.clear();
	OpHitsCosmic_Width.clear();

	for (int i_ohc = 0; i_ohc < int(ophit_cosmic_vec.size()); ++i_ohc) {

		art::Ptr<recob::OpHit> CurrentCosmicOpHit = ophit_cosmic_vec.at(i_ohc);

		OpHitsCosmic_Amplitude.push_back(CurrentCosmicOpHit->Amplitude());
		OpHitsCosmic_Area.push_back(CurrentCosmicOpHit->Area());
		OpHitsCosmic_OpChannel.push_back(CurrentCosmicOpHit->OpChannel());
		OpHitsCosmic_PE.push_back(CurrentCosmicOpHit->PE());
		OpHitsCosmic_PeakTime.push_back(CurrentCosmicOpHit->PeakTime());
		OpHitsCosmic_Width.push_back(CurrentCosmicOpHit->Width());

	} // End of the loop over the cosmic optical hits

// ____________________________________________________________________________________________________________________________________________________________________________________________________

	// Vertices 
std::cout << std::endl << "Loop over the vertices" << std::endl << std::endl;
	NumberVertices = vertex_vec.size();
	Vertex_ID.clear();
	Vertex_IsBeamEvent.clear();
	Vertex_IsCosmicEvent.clear();
	Vertex_PositionX.clear();
	Vertex_PositionY.clear();
	Vertex_PositionZ.clear();

	double VertexPosition[3] = {0,0,0};

	for (int i_v = 0; i_v < int(vertex_vec.size()); ++i_v) {

		art::Ptr<recob::Vertex> CurrentVertex = vertex_vec.at(i_v);

		Vertex_ID.push_back(CurrentVertex->ID());
		CurrentVertex->XYZ(VertexPosition);
		Vertex_PositionX.push_back(VertexPosition[0]);
		Vertex_PositionY.push_back(VertexPosition[1]);
		Vertex_PositionZ.push_back(VertexPosition[2]);

		art::Ptr<recob::PFParticle> CurrentPFParticle = VertexToPFParticle.at(CurrentVertex.key());

		if (CurrentPFParticle.isNull()) {

			Vertex_IsBeamEvent.push_back(0);
			Vertex_IsCosmicEvent.push_back(1);

		} else {

			art::Ptr<recob::Track> CurrentTrack = PFParticleToTrack.at(CurrentPFParticle.key());
			art::Ptr<recob::Shower> CurrentShower = PFParticleToShower.at(CurrentPFParticle.key());

			if (CurrentTrack.isNull() && CurrentShower.isNull()) {

				Vertex_IsBeamEvent.push_back(0);
				Vertex_IsCosmicEvent.push_back(1);

			} else {

				if (!CurrentTrack.isNull()) {

					std::vector< art::Ptr<recob::Hit> > trk_hits_ptrs_vertex = hits_per_track.at(CurrentTrack.key());
					BackTrackerTruthMatch backtrackertruthmatch_track_vertex; 
					backtrackertruthmatch_track_vertex.MatchToMCParticle(hit_handle,e,trk_hits_ptrs_vertex);
					art::Ptr< simb::MCParticle > maxp_me_track_vertex = backtrackertruthmatch_track_vertex.ReturnMCParticle();

					if (maxp_me_track_vertex.isNull()) {

						Vertex_IsBeamEvent.push_back(0);
						Vertex_IsCosmicEvent.push_back(1);

					} else {

						Vertex_IsBeamEvent.push_back(1);
						Vertex_IsCosmicEvent.push_back(0);

					}

				} // End of the track case

				if (!CurrentShower.isNull()) {

					std::vector< art::Ptr<recob::Hit> > shower_hits_ptrs_vertex = hits_per_shower.at(CurrentShower.key());
					BackTrackerTruthMatch backtrackertruthmatch_shower_vertex; 
					backtrackertruthmatch_shower_vertex.MatchToMCParticle(hit_handle,e,shower_hits_ptrs_vertex);
					art::Ptr< simb::MCParticle > maxp_me_shower_vertex = backtrackertruthmatch_shower_vertex.ReturnMCParticle();

					if (maxp_me_shower_vertex.isNull()) {

						Vertex_IsBeamEvent.push_back(0);
						Vertex_IsCosmicEvent.push_back(1);

					} else {

						Vertex_IsBeamEvent.push_back(1);
						Vertex_IsCosmicEvent.push_back(0);

					}					

				} // End of the shower case

			} // End of the if-statement whether there is a track or a shower

		} // End of the if-statement whether the VertexToPFParticle assosiation exists 

	} // End of the loop over the vertices

// ____________________________________________________________________________________________________________________________________________________________________________________________________

	// MCTruth
std::cout << std::endl << "Loop over the mctruth" << std::endl << std::endl;
	NumberMCTruthEvents = mclist.size();
	MCTruth_CCNC.clear();
	MCTruth_Mode.clear();
	MCTruth_Pt.clear();
	MCTruth_QSqr.clear();
	MCTruth_Target.clear();
	MCTruth_Theta.clear();
	MCTruth_X.clear();
	MCTruth_Y.clear();
	MCTruth_W.clear();
	MCTruth_NParticles.clear();
	MCTruth_Origin.clear();
	MCTruth_InteractionType.clear();
	MCTruth_Particle_Pdg.clear();
	MCTruth_Particle_TrackId.clear();
	MCTruth_Particle_EndX.clear();
	MCTruth_Particle_EndY.clear();
	MCTruth_Particle_EndZ.clear();
	MCTruth_Particle_Vx.clear();
	MCTruth_Particle_Vy.clear();
	MCTruth_Particle_Vz.clear();
	MCTruth_Particle_Nu_Vx.clear();
	MCTruth_Particle_Nu_Vy.clear();
	MCTruth_Particle_Nu_Vz.clear();

	for (int i_mclist = 0; i_mclist < NumberMCTruthEvents; i_mclist ++ ) {

		art::Ptr<simb::MCTruth> CurrentMCTruth = mclist.at(i_mclist);
		MCTruth_CCNC.push_back(CurrentMCTruth->GetNeutrino().CCNC());
		MCTruth_Mode.push_back(CurrentMCTruth->GetNeutrino().Mode());
		MCTruth_Pt.push_back(CurrentMCTruth->GetNeutrino().Pt());
		MCTruth_QSqr.push_back(CurrentMCTruth->GetNeutrino().QSqr());
		MCTruth_Target.push_back(CurrentMCTruth->GetNeutrino().Target());
		MCTruth_Theta.push_back(CurrentMCTruth->GetNeutrino().Theta());
		MCTruth_X.push_back(CurrentMCTruth->GetNeutrino().X());
		MCTruth_Y.push_back(CurrentMCTruth->GetNeutrino().Y());
		MCTruth_W.push_back(CurrentMCTruth->GetNeutrino().W());
		MCTruth_InteractionType.push_back(CurrentMCTruth->GetNeutrino().InteractionType());
		MCTruth_NParticles.push_back(CurrentMCTruth->NParticles());
		MCTruth_Origin.push_back(CurrentMCTruth->Origin());
		MCTruth_Particle_Nu_Vx.push_back(CurrentMCTruth->GetNeutrino().Nu().Vx());
		MCTruth_Particle_Nu_Vy.push_back(CurrentMCTruth->GetNeutrino().Nu().Vy());
		MCTruth_Particle_Nu_Vz.push_back(CurrentMCTruth->GetNeutrino().Nu().Vz());

		int NMCTruthParticles = CurrentMCTruth->NParticles();

		std::vector<int> CurrentMCTruth_Particle_Pdg;
		std::vector<int> CurrentMCTruth_Particle_TrackId;
		std::vector<double> CurrentMCTruth_Particle_EndX;
		std::vector<double> CurrentMCTruth_Particle_EndY;
		std::vector<double> CurrentMCTruth_Particle_EndZ;
		std::vector<double> CurrentMCTruth_Particle_Vx;
		std::vector<double> CurrentMCTruth_Particle_Vy;
		std::vector<double> CurrentMCTruth_Particle_Vz;

		CurrentMCTruth_Particle_Pdg.clear();
		CurrentMCTruth_Particle_TrackId.clear();
		CurrentMCTruth_Particle_EndX.clear();
		CurrentMCTruth_Particle_EndY.clear();
		CurrentMCTruth_Particle_EndZ.clear();
		CurrentMCTruth_Particle_Vx.clear();
		CurrentMCTruth_Particle_Vy.clear();
		CurrentMCTruth_Particle_Vz.clear();

		for (int WhichParticle = 0; WhichParticle < NMCTruthParticles; WhichParticle++) {

			const simb::MCParticle& mcparticle( CurrentMCTruth->GetParticle(WhichParticle) );
			CurrentMCTruth_Particle_Pdg.push_back(mcparticle.PdgCode());
			CurrentMCTruth_Particle_TrackId.push_back(mcparticle.TrackId());
			CurrentMCTruth_Particle_EndX.push_back(mcparticle.EndX());
			CurrentMCTruth_Particle_EndY.push_back(mcparticle.EndY());
			CurrentMCTruth_Particle_EndZ.push_back(mcparticle.EndZ());
			CurrentMCTruth_Particle_Vx.push_back(mcparticle.Vx());
			CurrentMCTruth_Particle_Vy.push_back(mcparticle.Vy());
			CurrentMCTruth_Particle_Vz.push_back(mcparticle.Vz());

		} // End of the loop over the MCParticles

		MCTruth_Particle_Pdg.push_back(CurrentMCTruth_Particle_Pdg);
		MCTruth_Particle_TrackId.push_back(CurrentMCTruth_Particle_TrackId);
		MCTruth_Particle_EndX.push_back(CurrentMCTruth_Particle_EndX);
		MCTruth_Particle_EndY.push_back(CurrentMCTruth_Particle_EndY);
		MCTruth_Particle_EndZ.push_back(CurrentMCTruth_Particle_EndZ);
		MCTruth_Particle_Vx.push_back(CurrentMCTruth_Particle_Vx);
		MCTruth_Particle_Vy.push_back(CurrentMCTruth_Particle_Vy);
		MCTruth_Particle_Vz.push_back(CurrentMCTruth_Particle_Vz);

	} // End of the loop over the mclist

// __________________________________________________________________________________________________________________________________________________________________________________________________

	// Tracking Efficiency

	// Loop over the PFParticles

	/*NumberPFParticles = pfparticle_vec.size();
	PFParticle_IsBeamEvent.clear();
	PFParticle_IsCosmicEvent.clear();
	PFParticle_IsPrimary.clear();
	PFParticle_NumDaughters.clear();

	for (int i_pf = 0; i_pf < NumberPFParticles; ++i_pf) {

		art::Ptr<recob::PFParticle> CurrentPFParticle = pfparticle_vec.at(i_pf);
		art::Ptr<recob::Vertex> CurrentVertex = PFParticleToVertex.at(CurrentPFParticle.key());
std::cout << "CurrentVertex->ID() = " << CurrentVertex->ID() << std::endl;
		art::Ptr<recob::Track> CurrentTrack = PFParticleToTrack.at(CurrentPFParticle.key());

		PFParticle_IsPrimary.push_back(CurrentPFParticle->IsPrimary());
		PFParticle_NumDaughters.push_back(CurrentPFParticle->NumDaughters());

		std::vector< art::Ptr<recob::Hit> > trk_hits_ptrs_pfpar = hits_per_track.at(CurrentTrack.key());
		BackTrackerTruthMatch backtrackertruthmatch_pfpar; 
		backtrackertruthmatch_pfpar.MatchToMCParticle(hit_handle,e,trk_hits_ptrs_pfpar);
		art::Ptr<simb::MCParticle> maxp_me_pfpar = backtrackertruthmatch_pfpar.ReturnMCParticle();

		if (maxp_me_pfpar.isNull()) { 

			PFParticle_IsBeamEvent.push_back(0);
			PFParticle_IsCosmicEvent.push_back(1);	

		} else {

			art::Ptr<simb::MCTruth> mctruth_pfpar = MCParticleToMCTruth.at(maxp_me_pfpar.key());

			if (mctruth_pfpar->Origin() == 1) { 

				PFParticle_IsBeamEvent.push_back(1);
				PFParticle_IsCosmicEvent.push_back(0);	
		
			} 

			if (mctruth_pfpar->Origin() == 2) { 

				PFParticle_IsBeamEvent.push_back(0);
				PFParticle_IsCosmicEvent.push_back(1);	
		
			}

		}

	}*/
// __________________________________________________________________________________________________________________________________________________________________________________________________

	std::cout << std::endl << "I just finished processing event # " << NEvents << std::endl << std::endl; 

	myTTree->Fill();

} // End of the analysis module
// __________________________________________________________________________________________________________________________________________________________________________________________________

void mynamespace::TTreeCreator::endJob()
{

	TFile& file = tfs->file();
	file.cd();

	std::cout << std::endl << "The TTree has been created and contains " << NEvents << " events." << std::endl << std::endl << std::cout << std::endl << std::endl; 
}
//___________________________________________________________________________________________________________________________________________________________________________________________________ 

DEFINE_ART_MODULE(mynamespace::TTreeCreator)
