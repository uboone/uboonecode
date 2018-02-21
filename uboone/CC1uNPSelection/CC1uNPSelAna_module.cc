////////////////////////////////////////////////////////////////////////
//
// @file CC1uNPSelAna_module.cc
//
// @brief A basic "skeleton" Ana module to serve as an example/basis
//        for the neutrino ID chain
//
// @authors usher@slac.stanford.edu (cloned from an example)
//
///////////////////////////////////////////////////////////////////////

#ifndef  CC1uNPSelAna_Module
#define  CC1uNPSelAna_Module

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTracker.h"
//include header files for new backtracker test
#include "uboone/AnalysisTree/MCTruth/IMCTruthMatching.h"
#include "uboone/AnalysisTree/MCTruth/AssociationsTruth_tool.h"
#include "uboone/AnalysisTree/MCTruth/BackTrackerTruth_tool.h"
//----------------------------------------------------------


#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Wire.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larreco/Deprecated/BezierTrack.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "uboone/EventWeight/MCEventWeight.h"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcore/Geometry/WireGeo.h"

#include "lardata/Utilities/AssociationUtil.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include <cstddef> // std::ptrdiff_t
#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <typeinfo>
#include <memory> // std::unique_ptr<>

#include "TTree.h"
#include "TTimeStamp.h"


constexpr int kmaxg4par= 90000;

constexpr int kNplanes       = 3;     //number of wire planes
constexpr int kMaxHits       = 40000; //maximum number of hits;
constexpr int kMaxTrackHits  = 2000;  //maximum number of hits on a track
constexpr int kMaxTrackers   = 15;    //number of trackers passed into fTrackModuleLabel
constexpr int kMaxVertices   = 500;    //max number of 3D vertices
constexpr int kMaxVertexAlgos = 10;    //max number of vertex algorithms
constexpr int kMaxFlashAlgos = 10;    //max number of flash algorithms
constexpr int kNOpDets = 32;          ///< number of optical detectors (PMTs)
constexpr unsigned short kMaxAuxDets = 4; ///< max number of auxiliary detector cells per MC particle
constexpr int kMaxFlashes    = 1000;   //maximum number of flashes
constexpr int kMaxShowerHits = 10000;  //maximum number of hits on a shower
constexpr int kMaxTruth      = 10;     //maximum number of neutrino truth interactions
constexpr int kMaxClusters   = 1000;   //maximum number of clusters;
constexpr int kMaxTicks   = 9600;   //maximum number of ticks (time samples)

constexpr int kMaxNDaughtersPerPFP = 100; //maximum number of daughters per PFParticle
constexpr int kMaxNClustersPerPFP  = 100; //maximum number of clusters per PFParticle
constexpr int kMaxNPFPNeutrinos    = 10;  //maximum number of reconstructed neutrino PFParticles

constexpr int kMaxSysts = 1000;
constexpr int kMaxWeights = 1000;

constexpr double pi = 3.1415;
constexpr float muonmass= 0.105658;
constexpr float protonmass=0.938;
/// total_extent\<T\>::value has the total number of elements of an array
template <typename T>
struct total_extent {
  using value_type = size_t;
  static constexpr value_type value
  = sizeof(T) / sizeof(typename std::remove_all_extents<T>::type);
}; // total_extent<>


namespace  CC1uNPSelAna
{
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class definition
  class CC1uNpSelAnaDataStruct {
  public:
    
    /// A wrapper to a C array (needed to embed an array into a vector)
    template <typename Array_t>
     class BoxedArray {
    protected:
      Array_t array; // actual data
      
    public:
      using This_t = BoxedArray<Array_t>;
      typedef typename std::remove_all_extents<Array_t>::type Data_t;
      
      BoxedArray() {} // no initialization
      BoxedArray(const This_t& from)
      { std::memcpy((char*) &(data()), (char*) &(from.data()), sizeof(Array_t)); }
      
      Array_t& data() { return array; }
      const Array_t& data() const { return array; }
      
      //@{
      /// begin/end interface
      static constexpr size_t size() { return total_extent<Array_t>::value; }
      Data_t* begin() { return reinterpret_cast<Data_t*>(&array); }
      const Data_t* begin() const { return reinterpret_cast<const Data_t*>(&array); }
      Data_t* end() { return begin() + size(); }
      const Data_t* end() const { return begin() + size(); }
      //@}
      
      //@{
      /// Array interface
      auto operator[] (size_t index) -> decltype(*array) { return array[index]; }
      auto operator[] (size_t index) const -> decltype(*array) { return array[index]; }
      auto operator+ (ptrdiff_t index) -> decltype(&*array) { return array + index; }
      auto operator+ (ptrdiff_t index) const -> decltype(&*array) { return array + index; }
      auto operator- (ptrdiff_t index) -> decltype(&*array) { return array - index; }
      auto operator- (ptrdiff_t index) const -> decltype(&*array) { return array - index; }
      auto operator* () -> decltype(*array) { return *array; }
      auto operator* () const -> decltype(*array) { return *array; }
      
      operator decltype(&array[0]) () { return &array[0]; }
      operator decltype(&array[0]) () const { return &array[0]; }
      //@}
      
    }; // BoxedArray
    class TrackDataStruct {
    public:
      /* Data structure size:
       *
       * TrackData_t<Short_t>                    :  2  bytes/track
       * TrackData_t<Float_t>                    :  4  bytes/track
       * PlaneData_t<Float_t>, PlaneData_t<Int_t>: 12  bytes/track
       * HitData_t<Float_t>                      : 24k bytes/track
       * HitCoordData_t<Float_t>                 : 72k bytes/track
       */
      template <typename T>
      using TrackData_t = std::vector<T>;
      template <typename T>
      using PlaneData_t = std::vector<BoxedArray<T[kNplanes]>>;
      template <typename T>
      using HitData_t = std::vector<BoxedArray<T[kNplanes][kMaxTrackHits]>>;
      template <typename T>
      using HitCoordData_t = std::vector<BoxedArray<T[kNplanes][kMaxTrackHits][3]>>;
      
      size_t MaxTracks; ///< maximum number of storable tracks
      
      Short_t  ntracks;             //number of reconstructed tracks
      PlaneData_t<Float_t>    trkke;
      PlaneData_t<Float_t>    trkrange;
      PlaneData_t<Int_t>      trkidtruth;  //true geant trackid
      PlaneData_t<Short_t>    trkorigin;   //_ev_origin 0: unknown, 1: neutrino, 2: cosmic, 3: supernova, 4: singles
      PlaneData_t<Int_t>      trkpdgtruth; //true pdg code
      PlaneData_t<Float_t>    trkefftruth; //completeness
      PlaneData_t<Float_t>    trkpurtruth; //purity of track
      PlaneData_t<Float_t>    trkpitchc;
      PlaneData_t<Short_t>    ntrkhits;
      HitData_t<Float_t>      trkdedx;
      HitData_t<Float_t>      trkdqdx;
      HitData_t<Float_t>      trkresrg;
      HitCoordData_t<Float_t> trkxyz;

      // more track info
      TrackData_t<Short_t> trkId;
      TrackData_t<Short_t> trkncosmictags_tagger;
      TrackData_t<Float_t> trkcosmicscore_tagger;
      TrackData_t<Short_t> trkcosmictype_tagger;
      TrackData_t<Short_t> trkncosmictags_containmenttagger;
      TrackData_t<Float_t> trkcosmicscore_containmenttagger;
      TrackData_t<Short_t> trkcosmictype_containmenttagger;
      TrackData_t<Short_t> trkncosmictags_flashmatch;
      TrackData_t<Float_t> trkcosmicscore_flashmatch;
      TrackData_t<Short_t> trkcosmictype_flashmatch;
      TrackData_t<Float_t> trkstartx;     // starting x position.
      TrackData_t<Float_t> trkstarty;     // starting y position.
      TrackData_t<Float_t> trkstartz;     // starting z position. 
      TrackData_t<Float_t> trkstartd;     // starting distance to boundary.
      TrackData_t<Float_t> trkendx;       // ending x position.
      TrackData_t<Float_t> trkendy;       // ending y position.
      TrackData_t<Float_t> trkendz;       // ending z position.
      TrackData_t<Float_t> trkendd;       // ending distance to boundary.
      TrackData_t<Float_t> trkACpierceT0;   // t0 per track from anode or cathode piercing tracks (in ns)     
      TrackData_t<Float_t> trkflashT0;   // t0 per track from matching tracks to flashes (in ns)
      TrackData_t<Float_t> trktrueT0;    // t0 per track from truth information (in ns)
      TrackData_t<Float_t> trkpurity;    // track purity based on hit information
      TrackData_t<Float_t> trkcompleteness; //track completeness based on hit information
      TrackData_t<int> trkg4id;        //true g4 track id for the reconstructed track
      TrackData_t<int> trkorig;        //origin of the track 
      TrackData_t<Float_t> trktheta;      // theta.
      TrackData_t<Float_t> trkphi;        // phi.
      TrackData_t<Float_t> trkstartdcosx;
      TrackData_t<Float_t> trkstartdcosy;
      TrackData_t<Float_t> trkstartdcosz;
      TrackData_t<Float_t> trkenddcosx;
      TrackData_t<Float_t> trkenddcosy;
      TrackData_t<Float_t> trkenddcosz;
      TrackData_t<Float_t> trkthetaxz;    // theta_xz.
      TrackData_t<Float_t> trkthetayz;    // theta_yz.
      TrackData_t<Float_t> trkmom;        // momentum.
      TrackData_t<Float_t> trklen;        // length.
      TrackData_t<Float_t> trkmomrange;    // track momentum from range using CSDA tables
      TrackData_t<Float_t> trkmommschi2;   // track momentum from multiple scattering Chi2 method
      TrackData_t<Float_t> trkmommsllhd;   // track momentum from multiple scattering LLHD method
      TrackData_t<Short_t> trksvtxid;     // Vertex ID associated with the track start
      TrackData_t<Short_t> trkevtxid;     // Vertex ID associated with the track end
      PlaneData_t<Int_t> trkpidpdg;       // particle PID pdg code
      PlaneData_t<Float_t> trkpidchi;
      PlaneData_t<Float_t> trkpidchipr;   // particle PID chisq for proton
      PlaneData_t<Float_t> trkpidchika;   // particle PID chisq for kaon
      PlaneData_t<Float_t> trkpidchipi;   // particle PID chisq for pion
      PlaneData_t<Float_t> trkpidchimu;   // particle PID chisq for muon
      PlaneData_t<Float_t> trkpidpida;    // particle PIDA
      TrackData_t<Short_t> trkpidbestplane; // this is defined as the plane with most hits   
	
      TrackData_t<Short_t> trkhasPFParticle; // whether this belongs to a PFParticle 
      TrackData_t<Short_t> trkPFParticleID;  // if hasPFParticle, its ID
       
      /// Creates an empty tracker data structure
      TrackDataStruct(): MaxTracks(0) { Clear(); }
      /// Creates a tracker data structure allowing up to maxTracks tracks
      TrackDataStruct(size_t maxTracks): MaxTracks(maxTracks) { Clear(); }
      void Clear();
      void SetMaxTracks(size_t maxTracks)
      { MaxTracks = maxTracks; Resize(MaxTracks); }
      void Resize(size_t nTracks);
      void SetAddresses(TTree* pTree, std::string tracker, bool isCosmics);
      size_t GetMaxTracks() const { return MaxTracks; }
      size_t GetMaxPlanesPerTrack(int /* iTrack */ = 0) const
      { return (size_t) kNplanes; }
      size_t GetMaxHitsPerTrack(int /* iTrack */ = 0, int /* ipl */ = 0) const
      { return (size_t) kMaxTrackHits; }
      
    }; // class TrackDataStruct






    //Vertex data struct
    class VertexDataStruct {
    public:
      template <typename T>
      using VertexData_t = std::vector<T>;

      size_t MaxVertices; ///< maximum number of storable vertices

      Short_t  nvtx;             //number of reconstructed vertices
      VertexData_t<Short_t> vtxId;    // the vertex ID.
      VertexData_t<Float_t> vtxx;     // x position.
      VertexData_t<Float_t> vtxy;     // y position.
      VertexData_t<Float_t> vtxz;     // z position.
	  
	  VertexData_t<Short_t> vtxhasPFParticle; // whether this belongs to a PFParticle 
	  VertexData_t<Short_t> vtxPFParticleID;  // if hasPFParticle, its ID

      VertexDataStruct(): MaxVertices(0) { Clear(); }
      VertexDataStruct(size_t maxVertices): MaxVertices(maxVertices) { Clear(); }
      void Clear();
      void SetMaxVertices(size_t maxVertices)
      { MaxVertices = maxVertices; Resize(MaxVertices); }
      void Resize(size_t nVertices);
      void SetAddresses(TTree* pTree, std::string tracker, bool isCosmics);

      size_t GetMaxVertices() const { return MaxVertices; }
    }; // class VertexDataStruct 

    //Neutrino vertex data struct
    class NeutrinoVertexDataStruct {
    public:
      template <typename T>
      using VertexData_t = std::vector<T>;

      size_t MaxVertices;             ///< maximum number of storable vertices

      Short_t  nvtx;                  ///< number of neutrino reconstructed vertices
      VertexData_t<Short_t> vtxId;    ///< the vertex ID.
      VertexData_t<Float_t> vtxx;     ///< x position.
      VertexData_t<Float_t> vtxy;     ///< y position.
      VertexData_t<Float_t> vtxz;     ///< z position.
      VertexData_t<Int_t>   vtxpdg;   ///< pdg of pfp

      VertexData_t<Short_t> vtxhasPFParticle; // whether this belongs to a PFParticle 
      VertexData_t<Short_t> vtxPFParticleID;  // if hasPFParticle, its ID

      NeutrinoVertexDataStruct(): MaxVertices(0) { Clear(); }
      NeutrinoVertexDataStruct(size_t maxVertices): MaxVertices(maxVertices) { Clear(); }
      void Clear();
      void SetMaxVertices(size_t maxVertices)   
      { MaxVertices = maxVertices; Resize(MaxVertices); }
      void Resize(size_t nVertices);
      void SetAddresses(TTree* pTree, std::string tracker, bool isCosmics);

      size_t GetMaxVertices() const { return MaxVertices; }
    }; // class VertexDataStruct 


    // Flash data struct
    class FlashDataStruct {
    public:
      template <typename T>
      using FlashData_t     = std::vector<T>;
      template <typename T>
      using FlashDataSpec_t = std::vector<BoxedArray<T[kNOpDets]>>;

      size_t MaxFlashes;                                ///< Maximum number of storable flashes

      Short_t nfls;                                     ///< Number of reconstructed flashes
      FlashData_t<Float_t> flsTime;                     ///< Flash time (us)
      FlashData_t<Float_t> flsPe;                       ///< Flash total PE
      FlashDataSpec_t<Double_t> flsPePerOpDet;          ///< Flash PE per optical detector
      FlashData_t<Float_t> flsXcenter;                  ///< Flash X center (cm)
      FlashData_t<Float_t> flsYcenter;                  ///< Flash Y center (cm)
      FlashData_t<Float_t> flsZcenter;                  ///< Flash Z center (cm)
      FlashData_t<Float_t> flsYwidth;                   ///< Flash Y width (cm)
      FlashData_t<Float_t> flsZwidth;                   ///< Flash Z width (cm)
      FlashData_t<Float_t> flsTwidth;                   ///< Flash time width (us)

      void Clear();
      void SetMaxFlashes(size_t maxFlashes)
      { MaxFlashes = maxFlashes; Resize(MaxFlashes); }
      void Resize(size_t nFlashes);
      void SetAddresses(TTree* pTree, std::string tracker, bool isCosmics);

      size_t GetMaxFlashes() const { return MaxFlashes; }

      size_t GetNOpDet() const { auto const* geom = lar::providerFrom<geo::Geometry>(); return geom->NOpDets(); }
    }; // class FlashDataStruct

    /// Shower algorithm result
    /// 
    /// Can connect to a tree, clear its fields and resize its data.
    class ShowerDataStruct {
    public:
      /* Data structure size:
       *
       * ShowerData_t<Short_t>                   :  2  bytes/shower
       * ShowerData_t<Float_t>                   :  4  bytes/shower
       * PlaneData_t<Float_t>, PlaneData_t<Int_t>: 12  bytes/shower
       * HitData_t<Float_t>                      : 24k bytes/shower
       * HitCoordData_t<Float_t>                 : 72k bytes/shower
       */
      template <typename T>
      using ShowerData_t = std::vector<T>;
      template <typename T>
      using PlaneData_t = std::vector<BoxedArray<T[kNplanes]>>;
      template <typename T>
      using HitData_t = std::vector<BoxedArray<T[kNplanes][kMaxShowerHits]>>;
      template <typename T>
      using HitCoordData_t = std::vector<BoxedArray<T[kNplanes][kMaxShowerHits][3]>>;
      
      std::string name; ///< name of the shower algorithm (for branch names)
      
      size_t MaxShowers; ///< maximum number of storable showers
      
      /// @{
      /// @name Branch data structures
      Short_t  nshowers;                      ///< number of showers
      ShowerData_t<Short_t>  showerID;        ///< Shower ID
      ShowerData_t<Short_t>  shwr_bestplane;  ///< Shower best plane
      ShowerData_t<Float_t>  shwr_length;     ///< Shower length
      ShowerData_t<Float_t>  shwr_startdcosx; ///< X directional cosine at start of shower
      ShowerData_t<Float_t>  shwr_startdcosy; ///< Y directional cosine at start of shower
      ShowerData_t<Float_t>  shwr_startdcosz; ///< Z directional cosine at start of shower
      ShowerData_t<Float_t>  shwr_startx;     ///< startx of shower
      ShowerData_t<Float_t>  shwr_starty;     ///< starty of shower
      ShowerData_t<Float_t>  shwr_startz;     ///< startz of shower
      PlaneData_t<Float_t>   shwr_totEng;     ///< Total energy of the shower per plane
      PlaneData_t<Float_t>   shwr_dedx;       ///< dE/dx of the shower per plane
      PlaneData_t<Float_t>   shwr_mipEng;     ///< Total MIP energy of the shower per plane
	  
	  ShowerData_t<Short_t>  shwr_hasPFParticle; // whether this belongs to a PFParticle 
	  ShowerData_t<Short_t>  shwr_PFParticleID;  // if hasPFParticle, its ID
      /// @}
      
      /// Creates a shower data structure allowing up to maxShowers showers
      ShowerDataStruct(std::string new_name = "", size_t maxShowers = 0):
        name(new_name), MaxShowers(maxShowers) { Clear(); }
      
      std::string Name() const { return name; }
      
      void Clear();    
      /// Applies a special prescription to mark shower information as missing
      void MarkMissing(TTree* pTree);
      void SetName(std::string new_name) { name = new_name; }
      void SetMaxShowers(size_t maxShowers)
      { MaxShowers = maxShowers; Resize(MaxShowers); }
      void Resize(size_t nShowers);
      void SetAddresses(TTree* pTree);
      
      size_t GetMaxShowers() const { return MaxShowers; }
      size_t GetMaxPlanesPerShower(int /* iShower */ = 0) const
      { return (size_t) kNplanes; }
      size_t GetMaxHitsPerShower(int /* iShower */ = 0, int /* ipl */ = 0) const
      { return (size_t) kMaxShowerHits; }
      
    }; // class ShowerDataStruct
    class PFParticleDataStruct {
    public:
      /* Data structure size:
       *
       * PFParticleData_t<Short_t>   :  2  bytes/PFParticle
       * PFParticleData_t<Int_t>     :  4  bytes/PFParticle
       * DaughterData_t<Short_t>     :  20 bytes/PFParticle
       * ClusterData_t<Short_t>      :  20 bytes/PFParticle
       * Short_t [kMaxNPFPNeutrinos] :  10 bytes in total
       */
      template <typename T>
      using PFParticleData_t = std::vector<T>;
      template <typename T>
      using DaughterData_t = std::vector<BoxedArray<T[kMaxNDaughtersPerPFP]>>;
      template <typename T>
      using ClusterData_t = std::vector<BoxedArray<T[kMaxNClustersPerPFP]>>;
      
      size_t MaxPFParticles; ///< maximum number of storable PFParticles
      
      /// @{
      /// @name Branch data structures
      Short_t                   nPFParticles;     ///< the total number of PFParticles
      PFParticleData_t<Short_t> pfp_selfID;       ///< the PFParticles' own IDs
      PFParticleData_t<Short_t> pfp_isPrimary;    ///< whether the PFParticle is a primary particle
      
      PFParticleData_t<Short_t> pfp_numDaughters; ///< the number of daughters belonging to this PFParticle
      DaughterData_t<Short_t>   pfp_daughterIDs;  ///< the IDs of the daughter PFParticles
      PFParticleData_t<Short_t> pfp_parentID;     ///< the ID of this PFParticle's immediate parent
      
      PFParticleData_t<Short_t> pfp_vertexID;     ///< the ID of the vertex belonging to this PFParticle
      PFParticleData_t<Short_t> pfp_isShower;     ///< whether this PFParticle corresponds to a shower
      PFParticleData_t<Short_t> pfp_isTrack;      ///< whether this PFParticle corresponds to a track
      PFParticleData_t<Short_t> pfp_trackID;      ///< the ID of the track object corresponding to this PFParticle, if !isShower
      PFParticleData_t<Short_t> pfp_showerID;     ///< the ID of the shower object corresponding to this PFParticle, if isShower
      
      PFParticleData_t<Short_t> pfp_isNeutrino;   ///< whether this PFParticle is a neutrino
      PFParticleData_t<Int_t>   pfp_pdgCode;      ///< the preliminary estimate of the PFParticle type using the PDG code
      
      PFParticleData_t<Short_t> pfp_numClusters;  ///< the number of associated clusters
      ClusterData_t<Short_t>    pfp_clusterIDs;   ///< the IDs of any associated clusters
      
      Short_t                   pfp_numNeutrinos; ///< the number of reconstructed neutrinos
      Short_t pfp_neutrinoIDs[kMaxNPFPNeutrinos]; ///< the PFParticle IDs of the neutrinos
      /// @}
      
      /// Creates a PFParticle data structure allowing up to maxPFParticles PFParticles
      PFParticleDataStruct(size_t maxPFParticles = 0):
        MaxPFParticles(maxPFParticles) { Clear(); }
      
      void Clear(); 
      void SetMaxPFParticles(size_t maxPFParticles)
        { MaxPFParticles = maxPFParticles; Resize(MaxPFParticles); }
      void Resize(size_t numPFParticles);
      void SetAddresses(TTree* pTree);
      
      size_t GetMaxPFParticles() const { return MaxPFParticles; }
      size_t GetMaxDaughtersPerPFParticle(int /* iPFParticle */ = 0) const
        { return (size_t) kMaxNDaughtersPerPFP; }
      size_t GetMaxClustersPerPFParticle(int /* iPFParticle */ = 0) const
        { return (size_t) kMaxNClustersPerPFP; }
      
    }; // class PFParticleDataStruct





    // Flash information
    Char_t kNFlashAlgos;
    std::vector<FlashDataStruct> FlashData;

    //track information
    Char_t   kNTracker;
    std::vector<TrackDataStruct> TrackData;

    //vertex information
    Char_t   kNVertexAlgos;
    std::vector<VertexDataStruct> VertexData;
   
    //neutrino vertex information
    Char_t   kNNeutrinoVertexAlgos;
    std::vector<NeutrinoVertexDataStruct> NeutrinoVertexData;

    // shower information
    Char_t   kNShowerAlgos;
    std::vector<ShowerDataStruct> ShowerData;

    // PFParticle information
    PFParticleDataStruct PFParticleData;
/*    
    // Raw Waveform information
    RawDataStruct RawData;
    
    // Calibration Waveform information
    CalibDataStruct CalibData;
*/


}; //end of CC1uNPSelAnaDataStruct 

class  CC1uNPSelAna : public art::EDAnalyzer
{
public:
 
    // Standard constructor and destructor for an ART module.
    explicit  CC1uNPSelAna(fhicl::ParameterSet const& pset);
    virtual ~ CC1uNPSelAna();

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

    bool inFV(double x, double y, double z) const;
    bool inTPC(double x, double y, double z) const;// starts in active TPC *** need to be done
    bool containmentFV(double startx, double starty, double startz, double endx, double endy, double endz) const;// starts & ends in FV *** need to be done
    bool containmentTPC(double startx, double starty, double startz, double endx, double endy, double endz) const;// starts & end in active TPC *** need to be done


    //declare the function to calculate truncated dQdx
    double GetDqDxTruncatedMean(std::vector<art::Ptr<anab::Calorimetry>> calos);
   
    double GetMean(std::vector<double> dqdx_v);

    double GetMedian(std::vector<double> dqdx_v);

    double GetVariance(std::vector<double> dqdx_v);

    double GetSTD(std::vector<double> dqdx_v);

    bool MIPConsistency(double dqds, double length);

   
    double GetFlashTrackDist(double flash, double start, double end) const;
    double GetFlashTrackDistMod(double flash, double start, double end) const;


    virtual void produces(art::EDProducer*);
 
    double GetTrackRange(art::Ptr<recob::Track>  InputTrackPtr) const;  
    double GetTrackLength(art::Ptr<recob::Track>  InputTrackPtr) const;  
    int GetNhitsTrk(art::Ptr<recob::Track> InputTrackPtr); 


    void truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet);
    double driftedLength(const simb::MCParticle& part, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi);
 
    double driftedLength(const sim::MCTrack& mctrack, TLorentzVector& tpcstart, TLorentzVector& tpcend, TLorentzVector& tpcmom);
 
    double Simlength(const simb::MCParticle& part, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi);

    double Recolength(const recob::Track& track);
    
    int Topology(int nmuons, int nelectrons, int npions, int npi0, int nprotons, bool cosmicflag, bool OOFVflag);
    //int Reaction(int mode);/// *** need to be finish                                                                                                        
    double GetDistTracktoVtx(art::Ptr<recob::Track> InputTrackPtr, TVector3 InputvertexPos);

    double GetTrackDirection(art::Ptr<recob::Track> InputTrackPtr, TVector3 InputvertexPos);/// *** need to finish/revisit                                  

    //double GetTruthNuEnergy();/// truth nu energy according to GENIE                                                                                        
    //double GetTruthInvMass();/// truth hadronic invariant mass according to GENIE                                                                           
    //double GetTruthQ2();/// truth Q^2 according to GENIE                                                                                                    
    //double GetTruthQ3();/// truth |Q_3| according to GENIE*** hacer                                                                                       

    /// **** need to be done                                                                                                                                
    /*                                                                                                                                                      
    double GetTruthNuEnergyRecoFormula();/// reconstructed nu energy using MC truth kinematics *** hacer                                                    
    double GetTruthInvMassRecoFormula();/// reconstructed invariant mass using MC truth kinematics *** hacer                                                
    double GetTruthQ2RecoFormula();/// reconstructed Q^2 using MC truth kinematics *** hacer                                                                
    double GetTruthQ3RecoFormula();/// reconstructed |Q_3| using MC truth kinematics *** hacer                                                              
                                                                                                                                                          
    double GetRecoNuEnergyRecoFormula();/// reconstructed nu energy using reconstructed kinematics *** hacer                                                
    double GetRecohInvMassRecoFormula();/// reconstructed invariant mass using reconstructed kinematics *** hacer                                           
    double GetRecoQ2RecoFormula();/// reconstructed Q^2 using reconstructed kinematics *** hacer                                                            
    double GetRecoQ3RecoFormula();/// reconstructed |Q_3| using reconstructed kinematics *** hacer                                                          
    */
    // Clear all fields if this object (not the tracker algorithm data)
    void ClearLocalData();
  

    // Override the response to input and output fils so we can get the
    // fully qualified path for argo
    //void respondToOpenInputFile(art::FileBlock const&);
    //void respontToOpenOutputFile(art::FileBlock const&);
    
    

    // The analysis routine, called once per event. 
    void analyze (const art::Event& evt); 
    TTree*    fMC_Truth;
    TTree*    fMC_Geant;

    TTree*   fMC_allsel;
    TTree*   fMC_flashwin;
    TTree*   fMC_flashtag;
    TTree*   fMC_vtxinFV;
    TTree*   fMC_ntrks;
    TTree*   fMC_noshwr;
    TTree*   fMC_trkfls;
    TTree*   fMC_NoExTrk;
    TTree*   fMC_mupinFV;
    TTree*   fMC_TrunMean;
private:


    //==================================================================
    // For keeping track of the replacement backtracker
    std::unique_ptr<truth::IMCTruthMatching> fMCTruthMatching;
 
    //------------------------------------------------------
    int   _fTrueccnc;
    int   _fTruemode;
    int   _fTrueinttype;
    int   _fTruenupdg;
    float _fTrueenu;
    float _fTrueq2truth;
    float _fTrueWtruth;
    float _fTrueXtruth;
    float _fTrueYtruth;
    float _fTruenuvrtxx;
    float _fTruenuvrtxy;
    float _fTruenuvrtxz;
    float _fTruenuvrtxx_SCE;
    float _fTruenuvrtxy_SCE;
    float _fTruenuvrtxz_SCE;
    int Nevt_truth; 
    //------------------------------------------------------
    int   _fmparpdg;
    int   _fmparStatusCode;
    float _fmparTheta;
    float _fmparCosTheta;	
    float _fmparSinTheta;
    float _fmparPhi;
    float _fmparCosPhi;
    float _fmparSinPhi;    
    float _fmparE;
    float _fmparMass;
    float _fmparKE;
    int   _fmctrue_origin;
    float _fmparEndE;
    float _fmparPx;
    float _fmparPy;
    float _fmparPz;
    float _fmparP;
    float _fmparStartx;
    float _fmparStarty;
    float _fmparStartz;
    float _fmparEndx;
    float _fmparEndy;
    float _fmparEndz;
    //arrays declared for the efficiency calculation
    
    std::vector<int> *fhg4parpdg;
    std::vector<int> *fhg4parstatus;
    std::vector<float> *fhg4parpx;
    std::vector<float> *fhg4parpy;
    std::vector<float> *fhg4parpz;
    std::vector<float> *fhg4partheta;
    std::vector<float> *fhg4parphi;
    std::vector<float> *fhg4parp;
    float vectex[kmaxg4par];
    //===================================================
    int nGEANTparticles;
    std::vector<int> _fg4parpdg;
    std::vector<int> _fg4parstatus; //status code of secondary particles assigned by geant
    std::vector<float> _fg4EngS; //energy of particles at started in GeV;
    std::vector<float> _fg4EngE; //energy of particles at ended in GeV;
    std::vector<float> _fg4px;
    std::vector<float> _fg4py;
    std::vector<float> _fg4pz;
    std::vector<float> _fg4P;
    std::vector<float> _fg4Mass;
    std::vector<float> _fg4StartPointx;
    std::vector<float> _fg4StartPointy;
    std::vector<float> _fg4StartPointz;
    std::vector<float> _fg4EndPointx;
    std::vector<float> _fg4EndPointy;
    std::vector<float> _fg4EndPointz;

    std::vector<float> _fg4SCEcorrStartPointx;
    std::vector<float> _fg4SCEcorrStartPointy;
    std::vector<float> _fg4SCEcorrStartPointz;



    std::vector<float> _fg4SCEcorrEndPointx;
    std::vector<float> _fg4SCEcorrEndPointy;
    std::vector<float> _fg4SCEcorrEndPointz;

    std::vector<float> _fg4StartT;
    std::vector<float> _fg4EndT;
    std::vector<float> _fg4theta;
    std::vector<float> _fg4phi;
    std::vector<float> _fg4theta_xz;
    std::vector<float> _fg4theta_yz;
    std::vector<float> _fg4pathlen;
    std::vector<int>  _fg4inTPCActive;
    std::vector<float> _fg4StartPointx_tpcAV;
    std::vector<float> _fg4StartPointy_tpcAV;
    std::vector<float> _fg4StartPointz_tpcAV;
    std::vector<float> _fg4EndPointx_tpcAV;
    std::vector<float> _fg4EndPointy_tpcAV;
    std::vector<float> _fg4EndPointz_tpcAV;
    std::vector<float> _fg4pathlen_drifted;
    std::vector<float> _fg4inTPCDrifted; 

    std::vector<int>   _fg4NumberDaughters;
    std::vector<int>   _fg4TrackId;
    std::vector<int>   _fg4Mother;
    std::vector<string>   _fg4processname; //Physics process by which the particle was created
    std::vector<int>   _fg4MergeId; //Geant track segments which belong to the same particle, get the same merge id
    std::vector<simb::Origin_t>   _fg4origin;  //1 cosmic, 2 neutrino 3 supernova  4 singles
    std::vector<int>   _fg4MCTruthIndex;
    //std::vector<int>   _fg4MCTruthIndex_new;
    std::vector<int> NuVertexID;
    std::vector<float> NuVertexx;
    std::vector<float> NuVertexy;
    std::vector<float> NuVertexz;
    

    //------------------------------------------------------------------
    // Need vectors here because we have have several instantiations
    // of the module running depending on matching
    std::vector<std::string> fVertexModuleLabelVec;
    std::vector<std::string> fPandoraNuVertexModuleLabelVec;
    std::vector<std::string> fVtxTrackAssnsModuleLabelVec;
    std::string              fInputFileName;
    //std::vector<std::string> fOpFlashModuleLabel; 
    // Pointers to the histograms we'll create. 

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;
    int truthtop; // true topology
    int truthtop_200thresh; // true topology assuming a 200MeV/c proton threshold
    int truthtop_300thresh; // true topology assuming a 300MeV/c proton threshold
    int truthtop_400thresh; // true topology assuming a 400MeV/c proton threshold

    double trueMuonTrueMomentum;
    double trueMuonTrueTheta;
    double trueMuonTruePhi;
    std::vector<double> *trueProtonsTrueMomentum;
    std::vector<double> *trueProtonsTrueTheta;
    std::vector<double> *trueProtonsTruePhi;

    float fLlep;
    float fLhad;
    float fPlep;
    float fPhad;
    float fThetaLep;
    float fThetaHad;
    float fCosThetaLep;
    float fCosThetaHad;
    float fPhiLep;
    float fPhiHad;
    int flashmax;
    float fopflashtime;
    float fopflashmax;
    float fvtxx, fvtxy,fvtxz;
    float flstrkdist;
    float vershwrdist; 
    float ftrklenmuoncand, ftrklenprotoncand;
    int Nhits_muoncand, Nhits_protoncand;
    //muon candidate info
    float trackendxcandidate=-999.0;
    float trackendycandidate=-999.0;
    float trackendzcandidate=-999.0;
    float trackstartxcandidate=-999.0;
    float trackstartycandidate=-999.0;
    float trackstartzcandidate=-999.0;
    float trackmomcandidate=-999.0;
    float trackmomcandidate_mcs=-999.0;
    std::vector<float> trackdedxcandidate;
    std::vector<float> trackresrgcandidate;

 
    //true information of track candidate
    int trackcand_origin=-999;
    int trackcand_nuset=-999;     
 
    int trackcand_parPDG=-999;
    int trackcand_parStatusCode=-999;
    float trackcand_parTheta=-999.0;                  
    float trackcand_parCosTheta=-999.0;
    float trackcand_parSinTheta=-999.0;                  
    float trackcand_parE=-999.0;        
    float trackcand_parMass=-999.0;
    float trackcand_parKE=-999.0;
    float trackcand_parEndE=-999.0;
    float trackcand_parPx=-999.0;
    float trackcand_parPy=-999.0;
    float trackcand_parPz=-999.0;
    float trackcand_parPhi=-999.0;
    float trackcand_parCosPhi=-999.0;
    float trackcand_parSinPhi=-999.0;

    //proton candidate info
    float trackendxprotoncandidate=-999.0;
    float trackendyprotoncandidate=-999.0;
    float trackendzprotoncandidate=-999.0;
    float trackstartxprotoncandidate=-999.0;
    float trackstartyprotoncandidate=-999.0;
    float trackstartzprotoncandidate=-999.0;
    float trackmomprotoncandidate=-999.0;
    std::vector<float> trackdedxprotoncandidate;
    std::vector<float> trackresrgprotoncandidate;

    int trackpcand_origin=-999;
    int trackpcand_nuset=-999;     
     
    int trackpcand_parPDG=-999;
    int trackpcand_parStatusCode=-999;
    float trackpcand_parTheta=-999.0;                  
    float trackpcand_parCosTheta=-999.0;
    float trackpcand_parSinTheta=-999.0;                  
    float trackpcand_parE=-999.0;        
    float trackpcand_parMass=-999.0;
    float trackpcand_parKE=-999.0;
    float trackpcand_parEndE=-999.0;
    float trackpcand_parPx=-999.0;
    float trackpcand_parPy=-999.0;
    float trackpcand_parPz=-999.0;
    float trackpcand_parPhi=-999.0;
    float trackpcand_parCosPhi=-999.0;
    float trackpcand_parSinPhi=-999.0;
    
    //=========================================
    float Evis=-999.0;
    float Q2cal=-999.0;
    float Wcal=-999.0;

    TLorentzVector *fHitNucP4;

    //=============================================================
    //declare the vectors for all the proton candidates
    std::vector<int> *trackidpcand;
    std::vector<float> *trackstartxpcand;
    std::vector<float> *trackstartypcand;    
    std::vector<float> *trackstartzpcand;
    std::vector<float> *trackendxpcand;
    std::vector<float> *trackendypcand;    
    std::vector<float> *trackendzpcand;
    
    std::vector<double> *trackmompcand;// = new std::vector<double>;
    std::vector<double> *trackthetapcand;// = new std::vector<double>;
    std::vector<double> *tracklengthpcand;// = new std::vector<double>;
    std::vector<double> *trackphipcand;// = new std::vector<double>;;
    std::vector<double> *tracktrunmeanpcand;// = new std::vector<double>;;
    std::vector<double> *trackpidapcand;
    //----------------------------------------
    int fNRecoTrks=-999;  //total number of tracks including muon, proton and others
    int fNRecoPTrks=-999; //total number of reco proton tracks
    int fNTruePTrks=-999; //total number of true proton tracks   


    float TrunMean_cand;
    float TrunMean_pcand;


    // Other variables that will be shared between different methods.
    detinfo::DetectorProperties const* fDetectorProperties; ///< Pointer to the detector properties
    //declare the fhicl parameters here
    std::string              fHitsModuleLabel ;    
    std::string              fTrackModuleLabel;    
    std::string              fVertexModuleLabel;   
    std::string              fPandoraNuVertexModuleLabel;
    std::string              fGenieGenModuleLabel;   
    std::string              fG4ModuleLabel; 
    std::string              fOpFlashModuleLabel;    
    std::string              fCalorimetryModuleLabel; 
    std::string              fShowerModuleLabel; 
    std::string              fTrackMCSFitLabel;



    //std::string fPFParticleModuleLabel;
    //what does this producer do??????? 
    art::EDProducer*           fMyProducerModule;        ///< The producer module driving us 

    float fG4minE;

    
    double fFlashWidth;      
    double fBeamMin;    
    double fBeamMax;      
    double fPEThresh;     
    double fMinTrk2VtxDist;     
    double fMinTrackLen;
    double fDistToEdgeX;
    double fDistToEdgeY;
    double fDistToEdgeZ;

    std::vector<double> _svm_x;

}; // class  CC1uNPSelAna


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// class implementation

//-----------------------------------------------------------------------
// Constructor
 CC1uNPSelAna:: CC1uNPSelAna(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet)
{
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
}

//-----------------------------------------------------------------------
// Destructor
 CC1uNPSelAna::~ CC1uNPSelAna()
{
  delete fhg4parpdg;
  delete fhg4parstatus;
  delete fhg4parpx;
  delete fhg4parpy; 
  delete fhg4parpz;
  delete fhg4parp;
  delete fhg4partheta;
  delete fhg4parphi; 

  delete trackidpcand;
  delete trackstartxpcand;
  delete trackstartypcand;
  delete trackstartzpcand;
  delete trackendxpcand;
  delete trackendypcand;
  delete trackendzpcand; 
  delete trackmompcand;
  delete trackphipcand;
  delete trackthetapcand;
  delete tracklengthpcand;
  delete tracktrunmeanpcand;

  delete trackpidapcand;
  delete fHitNucP4;
}
   
//-----------------------------------------------------------------------
void  CC1uNPSelAna::beginJob()
{
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
    //    fPDGCodeHist        = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
    // declar histograms or trees here
    //-----------------------------------------------------------
    //fill the ntuples here
    fMC_Truth=tfs->make<TTree>("fMC_Truth", "MC Holder"); 
    fMC_Truth->Branch("fRun",&fRun,"fRun/I");
    fMC_Truth->Branch("fSubRun",&fSubRun,"fSubRun/I");
    fMC_Truth->Branch("fEvent",&fEvent,"fEvent/I");
    fMC_Truth->Branch("_fTrueccnc",&_fTrueccnc,"_fTrueccnc/I");
    fMC_Truth->Branch("_fTruemode",&_fTruemode,"_fTruemode/I");
    fMC_Truth->Branch("_fTrueinttype",&_fTrueinttype,"_fTrueinttype/I");
    fMC_Truth->Branch("_fTruenupdg",&_fTruenupdg,"_fTruenupdg/I");
    fMC_Truth->Branch("_fTrueenu",&_fTrueenu,"_fTrueenu/F");
    fMC_Truth->Branch("_fTrueq2truth",&_fTrueq2truth,"_fTrueq2truth/F");
    fMC_Truth->Branch("_fTrueWtruth",&_fTrueWtruth,"_fTrueWtruth/F");
    fMC_Truth->Branch("_fTrueXtruth",&_fTrueXtruth,"_fTrueXtruth/F");
    fMC_Truth->Branch("_fTrueYtruth",&_fTrueYtruth,"_fTrueYtruth/F");
    fMC_Truth->Branch("_fTruenuvrtxx",&_fTruenuvrtxx,"_fTruenuvrtxx/F");
    fMC_Truth->Branch("_fTruenuvrtxy",&_fTruenuvrtxy,"_fTruenuvrtxy/F");
    fMC_Truth->Branch("_fTruenuvrtxz",&_fTruenuvrtxz,"_fTruenuvrtxz/F");
    
    fMC_Truth->Branch("_fTruenuvrtxx_SCE",&_fTruenuvrtxx_SCE,"_fTruenuvrtxx_SCE/F");
    fMC_Truth->Branch("_fTruenuvrtxy_SCE",&_fTruenuvrtxy_SCE,"_fTruenuvrtxy_SCE/F");
    fMC_Truth->Branch("_fTruenuvrtxz_SCE",&_fTruenuvrtxz_SCE,"_fTruenuvrtxz_SCE/F");
    //-------------------------------------------------------------------------------
    fMC_Geant=tfs->make<TTree>("fMC_Geant", "MC Holder"); 
    fMC_Geant->Branch("fRun",&fRun,"fRun/I");
    fMC_Geant->Branch("fSubRun",&fSubRun,"fSubRun/I");
    fMC_Geant->Branch("fEvent",&fEvent,"fEvent/I");
    fMC_Geant->Branch("nGEANTparticles", &nGEANTparticles , "nGEANTparticles/I");
    fMC_Geant->Branch("fhg4parpdg", "std::vector<int>",&fhg4parpdg);
    fMC_Geant->Branch("fhg4parstatus", "std::vector<int>", &fhg4parstatus);
    fMC_Geant->Branch("fhg4parphi", "std::vector<float>", &fhg4parphi);
    fMC_Geant->Branch("fhg4partheta","std::vector<float>", &fhg4partheta);
    fMC_Geant->Branch("fhg4parpx", "std::vector<float>", &fhg4parpx);
    fMC_Geant->Branch("fhg4parpy", "std::vector<float>", &fhg4parpy);
    fMC_Geant->Branch("fhg4parpz", "std::vector<float>", &fhg4parpz);
    fMC_Geant->Branch("fhg4parp",  "std::vector<float>", &fhg4parp);
    fMC_Geant->Branch("truthtop", &truthtop, "truthtop/I");
    fMC_Geant->Branch("truthtop_200thresh", &truthtop_200thresh, "truthtop_200thresh/I");
    fMC_Geant->Branch("truthtop_300thresh", &truthtop_300thresh, "truthtop_300thresh/I");
    fMC_Geant->Branch("truthtop_400thresh", &truthtop_400thresh, "truthtop_400thresh/I");
     
    //========================================================================================= 


    fMC_allsel=tfs->make<TTree>("fMC_allsel","Data Holder");    
    fMC_allsel->Branch("fRun",&fRun,"fRun/I");
    fMC_allsel->Branch("fSubRun",&fSubRun,"fSubRun/I");
    fMC_allsel->Branch("fEvent",&fEvent,"fEvent/I");
    fMC_allsel->Branch("fopflashtime", &fopflashtime, "fopflashtime/F");
    fMC_allsel->Branch("truthtop", &truthtop, "truthtop/I");
    fMC_allsel->Branch("truthtop_200thresh", &truthtop_200thresh, "truthtop_200thresh/I");
    fMC_allsel->Branch("truthtop_300thresh", &truthtop_300thresh, "truthtop_300thresh/I");
    fMC_allsel->Branch("truthtop_400thresh", &truthtop_400thresh, "truthtop_400thresh/I");

    //============================================================================================
 
    fMC_flashwin=tfs->make<TTree>("fMC_flashwin","Data Holder");    
    fMC_flashwin->Branch("fRun",&fRun,"fRun/I");
    fMC_flashwin->Branch("fSubRun",&fSubRun,"fSubRun/I");
    fMC_flashwin->Branch("fEvent",&fEvent,"fEvent/I");
    fMC_flashwin->Branch("fopflashtime", &fopflashtime, "fopflashtime/F");
    fMC_flashwin->Branch("fopflashmax",  &fopflashmax,  "fopflashmax/F");
    fMC_flashwin->Branch("truthtop", &truthtop, "truthtop/I");
    fMC_flashwin->Branch("truthtop_200thresh", &truthtop_200thresh, "truthtop_200thresh/I");
    fMC_flashwin->Branch("truthtop_300thresh", &truthtop_300thresh, "truthtop_300thresh/I");
    fMC_flashwin->Branch("truthtop_400thresh", &truthtop_400thresh, "truthtop_400thresh/I");



    //=====================================================================================
    fMC_flashtag=tfs->make<TTree>("fMC_flashtag","Data Holder");    
    fMC_flashtag->Branch("fRun",&fRun,"fRun/I");
    fMC_flashtag->Branch("fSubRun",&fSubRun,"fSubRun/I");
    fMC_flashtag->Branch("fEvent",&fEvent,"fEvent/I");
    fMC_flashtag->Branch("fopflashtime", &fopflashtime, "fopflashtime/F");
    fMC_flashtag->Branch("fopflashmax",  &fopflashmax,  "fopflashmax/F");
    fMC_flashtag->Branch("truthtop", &truthtop, "truthtop/I");
    fMC_flashtag->Branch("truthtop_200thresh", &truthtop_200thresh, "truthtop_200thresh/I");
    fMC_flashtag->Branch("truthtop_300thresh", &truthtop_300thresh, "truthtop_300thresh/I");
    fMC_flashtag->Branch("truthtop_400thresh", &truthtop_400thresh, "truthtop_400thresh/I");






    fMC_vtxinFV=tfs->make<TTree>("fMC_vtxinFV","Data Holder");    
    fMC_vtxinFV->Branch("fRun",&fRun,"fRun/I");
    fMC_vtxinFV->Branch("fSubRun",&fSubRun,"fSubRun/I");
    fMC_vtxinFV->Branch("fEvent",&fEvent,"fEvent/I");
    //branch for vertex position is here---------------------!
    fMC_vtxinFV->Branch("fvtxx", &fvtxx, "fvtxx/F");         
    fMC_vtxinFV->Branch("fvtxy", &fvtxy, "fvtxy/F");
    fMC_vtxinFV->Branch("fvtxz", &fvtxz, "fvtxz/F");
    //-------------------------------------------------------!
    fMC_vtxinFV->Branch("fopflashtime", &fopflashtime, "fopflashtime/F");
    fMC_vtxinFV->Branch("fopflashmax",  &fopflashmax,  "fopflashmax/F");
    fMC_vtxinFV->Branch("truthtop", &truthtop, "truthtop/I");
    fMC_vtxinFV->Branch("truthtop_200thresh", &truthtop_200thresh, "truthtop_200thresh/I");
    fMC_vtxinFV->Branch("truthtop_300thresh", &truthtop_300thresh, "truthtop_300thresh/I");
    fMC_vtxinFV->Branch("truthtop_400thresh", &truthtop_400thresh, "truthtop_400thresh/I");


    //===================================================================================

    fMC_ntrks=tfs->make<TTree>("fMC_ntrks","Data Holder");    
    fMC_ntrks->Branch("fRun",&fRun,"fRun/I");
    fMC_ntrks->Branch("fSubRun",&fSubRun,"fSubRun/I");
    fMC_ntrks->Branch("fEvent",&fEvent,"fEvent/I");
    fMC_ntrks->Branch("fvtxx", &fvtxx, "fvtxx/F");         
    fMC_ntrks->Branch("fvtxy", &fvtxy, "fvtxy/F");
    fMC_ntrks->Branch("fvtxz", &fvtxz, "fvtxz/F");
    fMC_ntrks->Branch("fopflashtime", &fopflashtime, "fopflashtime/F");
    fMC_ntrks->Branch("fopflashmax",  &fopflashmax,  "fopflashmax/F");
    fMC_ntrks->Branch("vershwrdist",  &vershwrdist,  "vershwrdist/F");
    fMC_ntrks->Branch("truthtop", &truthtop, "truthtop/I");
    fMC_ntrks->Branch("truthtop_200thresh", &truthtop_200thresh, "truthtop_200thresh/I");
    fMC_ntrks->Branch("truthtop_300thresh", &truthtop_300thresh, "truthtop_300thresh/I");
    fMC_ntrks->Branch("truthtop_400thresh", &truthtop_400thresh, "truthtop_400thresh/I");

    //===================================================================================
    fMC_noshwr=tfs->make<TTree>("fMC_noshwr","Data Holder");    
    fMC_noshwr->Branch("fRun",&fRun,"fRun/I");
    fMC_noshwr->Branch("fSubRun",&fSubRun,"fSubRun/I");
    fMC_noshwr->Branch("fEvent",&fEvent,"fEvent/I");
    fMC_noshwr->Branch("fvtxx", &fvtxx, "fvtxx/F");         
    fMC_noshwr->Branch("fvtxy", &fvtxy, "fvtxy/F");
    fMC_noshwr->Branch("fvtxz", &fvtxz, "fvtxz/F");
    fMC_noshwr->Branch("fopflashtime", &fopflashtime, "fopflashtime/F");
    fMC_noshwr->Branch("fopflashmax",  &fopflashmax,  "fopflashmax/F");
    fMC_noshwr->Branch("vershwrdist",  &vershwrdist,  "vershwrdist/F"); 
    fMC_noshwr->Branch("truthtop", &truthtop, "truthtop/I");
    fMC_noshwr->Branch("truthtop_200thresh", &truthtop_200thresh, "truthtop_200thresh/I");
    fMC_noshwr->Branch("truthtop_300thresh", &truthtop_300thresh, "truthtop_300thresh/I");
    fMC_noshwr->Branch("truthtop_400thresh", &truthtop_400thresh, "truthtop_400thresh/I");

    //===============================================================================
    fMC_trkfls=tfs->make<TTree>("fMC_trkfls","Data Holder");    
    fMC_trkfls->Branch("fRun",&fRun,"fRun/I");
    fMC_trkfls->Branch("fSubRun",&fSubRun,"fSubRun/I");
    fMC_trkfls->Branch("fEvent",&fEvent,"fEvent/I");

    fMC_trkfls->Branch("fvtxx", &fvtxx, "fvtxx/F");         
    fMC_trkfls->Branch("fvtxy", &fvtxy, "fvtxy/F");
    fMC_trkfls->Branch("fvtxz", &fvtxz, "fvtxz/F");

    fMC_trkfls->Branch("fopflashtime", &fopflashtime, "fopflashtime/F");
    fMC_trkfls->Branch("fopflashmax",  &fopflashmax,  "fopflashmax/F");
    fMC_trkfls->Branch("flstrkdist" ,  &flstrkdist,   "flstrkdist/F");
    fMC_trkfls->Branch("truthtop", &truthtop, "truthtop/I");
    fMC_trkfls->Branch("truthtop_200thresh", &truthtop_200thresh, "truthtop_200thresh/I");
    fMC_trkfls->Branch("truthtop_300thresh", &truthtop_300thresh, "truthtop_300thresh/I");
    fMC_trkfls->Branch("truthtop_400thresh", &truthtop_400thresh, "truthtop_400thresh/I");

    //================================================================= 
    fMC_NoExTrk=tfs->make<TTree>("fMC_NoExTrk","Data Holder");    
    fMC_NoExTrk->Branch("fRun",&fRun,"fRun/I");
    fMC_NoExTrk->Branch("fSubRun",&fSubRun,"fSubRun/I");
    fMC_NoExTrk->Branch("fEvent",&fEvent,"fEvent/I");

    fMC_NoExTrk->Branch("fvtxx", &fvtxx, "fvtxx/F");         
    fMC_NoExTrk->Branch("fvtxy", &fvtxy, "fvtxy/F");
    fMC_NoExTrk->Branch("fvtxz", &fvtxz, "fvtxz/F");

    fMC_NoExTrk->Branch("fopflashtime", &fopflashtime, "fopflashtime/F");
    fMC_NoExTrk->Branch("fopflashmax",  &fopflashmax,  "fopflashmax/F");
    fMC_NoExTrk->Branch("flstrkdist" ,  &flstrkdist,   "flstrkdist/F");
    fMC_NoExTrk->Branch("truthtop", &truthtop, "truthtop/I");
    fMC_NoExTrk->Branch("truthtop_200thresh", &truthtop_200thresh, "truthtop_200thresh/I");
    fMC_NoExTrk->Branch("truthtop_300thresh", &truthtop_300thresh, "truthtop_300thresh/I");
    fMC_NoExTrk->Branch("truthtop_400thresh", &truthtop_400thresh, "truthtop_400thresh/I");



    //===============================================================
    fMC_mupinFV=tfs->make<TTree>("fMC_mupinFV","Data Holder");    
    fMC_mupinFV->Branch("fRun",&fRun,"fRun/I");
    fMC_mupinFV->Branch("fSubRun",&fSubRun,"fSubRun/I");
    fMC_mupinFV->Branch("fEvent",&fEvent,"fEvent/I");

    fMC_mupinFV->Branch("fvtxx", &fvtxx, "fvtxx/F");         
    fMC_mupinFV->Branch("fvtxy", &fvtxy, "fvtxy/F");
    fMC_mupinFV->Branch("fvtxz", &fvtxz, "fvtxz/F");
 
    fMC_mupinFV->Branch("fLlep", &fLlep, "fLlep/F");
    fMC_mupinFV->Branch("fLhad", &fLhad, "fLhad/F");
    fMC_mupinFV->Branch("fPlep", &fPlep, "fPlep/F");
    fMC_mupinFV->Branch("fPhad", &fPhad, "fPhad/F");
    fMC_mupinFV->Branch("fThetaLep",&fThetaLep, "fThetaLep/F");
    fMC_mupinFV->Branch("fThetaHad",&fThetaHad, "fThetaHad/F");
    fMC_mupinFV->Branch("fPhiLep",&fPhiLep, "fPhiLep/F");
    fMC_mupinFV->Branch("fPhiHad",&fPhiHad, "fPhiHad/F");
    fMC_mupinFV->Branch("fCosThetaLep",&fCosThetaLep, "fCosThetaLep/F");
    fMC_mupinFV->Branch("fCosThetaHad",&fCosThetaHad, "fCosThetaHad/F");
    fMC_mupinFV->Branch("trackendxcandidate", &trackendxcandidate, "trackendxcandidate/F");
    fMC_mupinFV->Branch("trackendycandidate", &trackendycandidate, "trackendycandidate/F");
    fMC_mupinFV->Branch("trackendzcandidate", &trackendzcandidate, "trackendzcandidate/F");
    fMC_mupinFV->Branch("trackstartxcandidate", &trackstartxcandidate, "trackstartxcandidate/F");
    fMC_mupinFV->Branch("trackstartycandidate", &trackstartycandidate, "trackstartycandidate/F");
    fMC_mupinFV->Branch("trackstartzcandidate", &trackstartzcandidate, "trackstartzcandidate/F");
    fMC_mupinFV->Branch("fNRecoTrks", &fNRecoTrks, "fNRecoTrks/I");
    fMC_mupinFV->Branch("fNRecoPTrks", &fNRecoPTrks, "fNRecoPTrks/I");
    fMC_mupinFV->Branch("fNTruePTrks", &fNTruePTrks, "fNTruePTrks/I");


    fMC_mupinFV->Branch("ftrklenmuoncand", &ftrklenmuoncand, "ftrklenmuoncand/F");
    fMC_mupinFV->Branch("ftrklenprotoncand", &ftrklenprotoncand, "ftrklenprotoncand/F");
    fMC_mupinFV->Branch("trackmomcandidate", &trackmomcandidate, "trackmomcandidate/F");
    fMC_mupinFV->Branch("trackmomcandidate_mcs", &trackmomcandidate_mcs, "trackmomcandidate_mcs/F");
    fMC_mupinFV->Branch("trackmomprotoncandidate", &trackmomprotoncandidate, "trackmomprotoncandidate/F");
 



    fMC_mupinFV->Branch("fopflashtime", &fopflashtime, "fopflashtime/F");
    fMC_mupinFV->Branch("fopflashmax",  &fopflashmax,  "fopflashmax/F");
    fMC_mupinFV->Branch("flstrkdist" ,  &flstrkdist,   "flstrkdist/F");


    fMC_mupinFV->Branch("truthtop", &truthtop, "truthtop/I");
    fMC_mupinFV->Branch("truthtop_200thresh", &truthtop_200thresh, "truthtop_200thresh/I");
    fMC_mupinFV->Branch("truthtop_300thresh", &truthtop_300thresh, "truthtop_300thresh/I");
    fMC_mupinFV->Branch("truthtop_400thresh", &truthtop_400thresh, "truthtop_400thresh/I");


    //=================================================================
    fMC_TrunMean=tfs->make<TTree>("fMC_TrunMean","Data Holder");    
    fMC_TrunMean->Branch("fRun",&fRun,"fRun/I");
    fMC_TrunMean->Branch("fSubRun",&fSubRun,"fSubRun/I");
    fMC_TrunMean->Branch("fEvent",&fEvent,"fEvent/I");
    //--------------------------------------------------
    fMC_TrunMean->Branch("_fTrueccnc",&_fTrueccnc,"_fTrueccnc/I");
    fMC_TrunMean->Branch("_fTruemode",&_fTruemode,"_fTruemode/I");
    fMC_TrunMean->Branch("_fTrueinttype",&_fTrueinttype,"_fTrueinttype/I");
    fMC_TrunMean->Branch("_fTruenupdg",&_fTruenupdg,"_fTruenupdg/I");
    fMC_TrunMean->Branch("_fTrueenu",&_fTrueenu,"_fTrueenu/F");
    fMC_TrunMean->Branch("_fTrueq2truth",&_fTrueq2truth,"_fTrueq2truth/F");
    fMC_TrunMean->Branch("_fTrueWtruth",&_fTrueWtruth,"_fTrueWtruth/F");
    fMC_TrunMean->Branch("_fTrueXtruth",&_fTrueXtruth,"_fTrueXtruth/F");
    fMC_TrunMean->Branch("_fTrueYtruth",&_fTrueYtruth,"_fTrueYtruth/F");
  
    fMC_TrunMean->Branch("trueMuonTrueMomentum",&trueMuonTrueMomentum,"trueMuonTrueMomentum/D");
    fMC_TrunMean->Branch("trueMuonTrueTheta",&trueMuonTrueTheta,"trueMuonTrueTheta/D");
    fMC_TrunMean->Branch("trueMuonTruePhi",&trueMuonTruePhi,"trueMuonTruePhi/D");

    fMC_TrunMean->Branch("trueProtonsTrueMomentum","std::vector<double>",&trueProtonsTrueMomentum);
    fMC_TrunMean->Branch("trueProtonsTrueTheta","std::vector<double>",&trueProtonsTrueTheta);
    fMC_TrunMean->Branch("trueProtonsTruePhi","std::vector<double>",&trueProtonsTruePhi);


    //---------------------------------------------------------

    fMC_TrunMean->Branch("fvtxx", &fvtxx, "fvtxx/F");         
    fMC_TrunMean->Branch("fvtxy", &fvtxy, "fvtxy/F");
    fMC_TrunMean->Branch("fvtxz", &fvtxz, "fvtxz/F");
 
    fMC_TrunMean->Branch("fLlep", &fLlep, "fLlep/F");
    fMC_TrunMean->Branch("fLhad", &fLhad, "fLhad/F");
    fMC_TrunMean->Branch("fPlep", &fPlep, "fPlep/F");
    fMC_TrunMean->Branch("fPhad", &fPhad, "fPhad/F");
    fMC_TrunMean->Branch("fThetaLep",&fThetaLep, "fThetaLep/F");
    fMC_TrunMean->Branch("fThetaHad",&fThetaHad, "fThetaHad/F");
    fMC_TrunMean->Branch("fPhiLep",&fPhiLep, "fPhiLep/F");
    fMC_TrunMean->Branch("fPhiHad",&fPhiHad, "fPhiHad/F");
    fMC_TrunMean->Branch("fCosThetaLep",&fCosThetaLep, "fCosThetaLep/F");
    fMC_TrunMean->Branch("fCosThetaHad",&fCosThetaHad, "fCosThetaHad/F");
    fMC_TrunMean->Branch("trackendxcandidate", &trackendxcandidate, "trackendxcandidate/F");
    fMC_TrunMean->Branch("trackendycandidate", &trackendycandidate, "trackendycandidate/F");
    fMC_TrunMean->Branch("trackendzcandidate", &trackendzcandidate, "trackendzcandidate/F");
    fMC_TrunMean->Branch("trackstartxcandidate", &trackstartxcandidate, "trackstartxcandidate/F");
    fMC_TrunMean->Branch("trackstartycandidate", &trackstartycandidate, "trackstartycandidate/F");
    fMC_TrunMean->Branch("trackstartzcandidate", &trackstartzcandidate, "trackstartzcandidate/F");
    fMC_TrunMean->Branch("fNRecoTrks", &fNRecoTrks, "fNRecoTrks/I");
    
    fMC_TrunMean->Branch("TrunMean_cand", &TrunMean_cand, "TrunMean_cand/F");
    fMC_TrunMean->Branch("TrunMean_pcand", &TrunMean_pcand, "TrunMean_pcand/F");

    fMC_TrunMean->Branch("fNRecoPTrks", &fNRecoPTrks, "fNRecoPTrks/I");
    fMC_TrunMean->Branch("fNTruePTrks", &fNTruePTrks, "fNTruePTrks/I");

    fMC_TrunMean->Branch("ftrklenmuoncand", &ftrklenmuoncand, "ftrklenmuoncand/F");
    fMC_TrunMean->Branch("ftrklenprotoncand", &ftrklenprotoncand, "ftrklenprotoncand/F");

    fMC_TrunMean->Branch("trackmomcandidate", &trackmomcandidate, "trackmomcandidate/F");
    fMC_TrunMean->Branch("trackmomcandidate_mcs", &trackmomcandidate_mcs, "trackmomcandidate_mcs/F");
    fMC_TrunMean->Branch("trackmomprotoncandidate", &trackmomprotoncandidate, "trackmomprotoncandidate/F");


    fMC_TrunMean->Branch("fopflashtime", &fopflashtime, "fopflashtime/F");
    fMC_TrunMean->Branch("fopflashmax",  &fopflashmax,  "fopflashmax/F");
    fMC_TrunMean->Branch("flstrkdist" ,  &flstrkdist,   "flstrkdist/F");


    fMC_TrunMean->Branch("Nhits_muoncand", &Nhits_muoncand, "Nhits_muoncand/I");
    fMC_TrunMean->Branch("Nhits_protoncand", &Nhits_protoncand, "Nhits_protoncand/I");
    fMC_TrunMean->Branch("truthtop", &truthtop, "truthtop/I");
    fMC_TrunMean->Branch("truthtop_200thresh", &truthtop_200thresh, "truthtop_200thresh/I");
    fMC_TrunMean->Branch("truthtop_300thresh", &truthtop_300thresh, "truthtop_300thresh/I");
    fMC_TrunMean->Branch("truthtop_400thresh", &truthtop_400thresh, "truthtop_400thresh/I");

    fMC_TrunMean->Branch("trackcand_origin", &trackcand_origin, "trackcand_origin/I");
    fMC_TrunMean->Branch("trackcand_nuset", &trackcand_nuset, "trackcand_nuset/I");
 
    fMC_TrunMean->Branch("trackcand_parPDG", &trackcand_parPDG, "trackcand_parPDG/I");
    fMC_TrunMean->Branch("trackcand_parStatusCode",&trackcand_parStatusCode, "trackcand_parStatusCode/I");
    fMC_TrunMean->Branch("trackcand_parTheta",&trackcand_parTheta, "trackcand_parTheta/F");                  
    fMC_TrunMean->Branch("trackcand_parCosTheta",&trackcand_parCosTheta, "trackcand_parCosTheta/F");
    fMC_TrunMean->Branch("trackcand_parSinTheta",&trackcand_parSinTheta, "trackcand_parSinTheta/F");                  
    fMC_TrunMean->Branch("trackcand_parE", &trackcand_parE, "trackcand_parE/F");        
    fMC_TrunMean->Branch("trackcand_parMass",&trackcand_parMass, "trackcand_parMass/F");
    fMC_TrunMean->Branch("trackcand_parKE", &trackcand_parKE, "trackcand_parKE/F");
    fMC_TrunMean->Branch("trackcand_parEndE", &trackcand_parEndE, "trackcand_parEndE/F");
    fMC_TrunMean->Branch("trackcand_parPx", &trackcand_parPx, "trackcand_parPx/F");
    fMC_TrunMean->Branch("trackcand_parPy", &trackcand_parPy, "trackcand_parPy/F");
    fMC_TrunMean->Branch("trackcand_parPz", &trackcand_parPz, "trackcand_parPz/F");
    fMC_TrunMean->Branch("trackcand_parPhi", &trackcand_parPhi, "trackcand_parPhi/F");
    fMC_TrunMean->Branch("trackcand_parCosPhi",&trackcand_parCosPhi, "trackcand_parCosPhi/F");
    fMC_TrunMean->Branch("trackcand_parSinPhi",&trackcand_parSinPhi, "trackcand_parSinPhi/F");

    fMC_TrunMean->Branch("trackpcand_origin", &trackpcand_origin, "trackpcand_origin/I");
    fMC_TrunMean->Branch("trackpcand_nuset", &trackpcand_nuset, "trackpcand_nuset/I");
 
    fMC_TrunMean->Branch("trackpcand_parPDG", &trackpcand_parPDG, "trackpcand_parPDG/I");
    fMC_TrunMean->Branch("trackpcand_parStatusCode",&trackpcand_parStatusCode, "trackpcand_parStatusCode/I");
    fMC_TrunMean->Branch("trackpcand_parTheta",&trackpcand_parTheta, "trackpcand_parTheta/F");                  
    fMC_TrunMean->Branch("trackpcand_parCosTheta",&trackpcand_parCosTheta, "trackpcand_parCosTheta/F");
    fMC_TrunMean->Branch("trackpcand_parSinTheta",&trackpcand_parSinTheta, "trackpcand_parSinTheta/F");                  
    fMC_TrunMean->Branch("trackpcand_parE", &trackpcand_parE, "trackpcand_parE/F");        
    fMC_TrunMean->Branch("trackpcand_parMass",&trackpcand_parMass, "trackpcand_parMass/F");
    fMC_TrunMean->Branch("trackpcand_parKE", &trackpcand_parKE, "trackpcand_parKE/F");
    fMC_TrunMean->Branch("trackpcand_parEndE", &trackpcand_parEndE, "trackpcand_parEndE/F");
    fMC_TrunMean->Branch("trackpcand_parPx", &trackpcand_parPx, "trackpcand_parPx/F");
    fMC_TrunMean->Branch("trackpcand_parPy", &trackpcand_parPy, "trackpcand_parPy/F");
    fMC_TrunMean->Branch("trackpcand_parPz", &trackpcand_parPz, "trackpcand_parPz/F");
    fMC_TrunMean->Branch("trackpcand_parPhi", &trackpcand_parPhi, "trackpcand_parPhi/F");
    fMC_TrunMean->Branch("trackpcand_parCosPhi",&trackpcand_parCosPhi, "trackpcand_parCosPhi/F");
    fMC_TrunMean->Branch("trackpcand_parSinPhi",&trackpcand_parSinPhi, "trackpcand_parSinPhi/F");

    fMC_TrunMean->Branch("Evis", &Evis, "Evis/F");
    fMC_TrunMean->Branch("Q2cal", &Q2cal, "Q2cal/F");
    fMC_TrunMean->Branch("Wcal", &Wcal, "Wcal/F");

    fMC_TrunMean->Branch("fHitNucP4","TLorentzVector",&fHitNucP4);
     
    //example of vector branch adding
    fMC_TrunMean->Branch("protoncandidate_id", "std::vector<int>", &trackidpcand);
    fMC_TrunMean->Branch("protoncandidate_startx", "std::vector<float>", &trackstartxpcand);
    fMC_TrunMean->Branch("protoncandidate_starty", "std::vector<float>", &trackstartypcand);
    fMC_TrunMean->Branch("protoncandidate_startz", "std::vector<float>", &trackstartzpcand);
    fMC_TrunMean->Branch("protoncandidate_endx", "std::vector<float>", &trackendxpcand);
    fMC_TrunMean->Branch("protoncandidate_endy", "std::vector<float>", &trackendypcand);
    fMC_TrunMean->Branch("protoncandidate_endz", "std::vector<float>", &trackendzpcand);

    fMC_TrunMean->Branch("protoncandidate_momentum","std::vector<double>",&trackmompcand);
    fMC_TrunMean->Branch("protoncandidate_length","std::vector<double>",&tracklengthpcand);
    fMC_TrunMean->Branch("protoncandidate_theta","std::vector<double>",&trackthetapcand);
    fMC_TrunMean->Branch("protoncandidate_phi","std::vector<double>",&trackphipcand);
    fMC_TrunMean->Branch("protoncandidate_trunmeandqdx","std::vector<double>",&tracktrunmeanpcand);
    fMC_TrunMean->Branch("protoncandidate_pida", "std::vector<double>", &trackpidapcand);

    //-----------------------------------------------------------
}
   
//-----------------------------------------------------------------------
void  CC1uNPSelAna::beginRun(const art::Run& /*run*/)
{
}

//-----------------------------------------------------------------------
void  CC1uNPSelAna::reconfigure(fhicl::ParameterSet const& pset)
{
    // Read parameters from the .fcl file.
    //fVertexModuleLabelVec        = pset.get< std::vector<std::string> >("VertexModuleLabelVec",       std::vector<std::string>() ={"pandoraNu"});

    fVertexModuleLabelVec        = pset.get< std::vector<std::string> >("VertexModuleLabelVec",       std::vector<std::string>() ={"vertex3d"});
    fPandoraNuVertexModuleLabelVec  = pset.get< std::vector<std::string> >("PandoraNuVertexModuleLabelVec", std::vector<std::string>()={"pandoraNu"}); 
    //fPandoraNuVertexModuleLabelVec  = pset.get< std::vector<std::string> >("PandoraNuVertexModuleLabel", std::vector<std::string>()={"pandoraNu"}); 
  
    fVtxTrackAssnsModuleLabelVec = pset.get< std::vector<std::string> >("VtxTrackAssnModuleLabelVec", std::vector<std::string>() ={"neutrinoID"});
    //fOpFlashModuleLabel          = pset.get< std::vector<std::string> >("OpFlashModuleLabel",         std::vector<std::string>() ={"simpleFlashBeam"});

    if (fVertexModuleLabelVec.size() != fVtxTrackAssnsModuleLabelVec.size())
    {
        mf::LogError("TPCNeutrinoIDFilter") << "Mismatch between string vector lengths input from fhicl!" << std::endl;
    }
    
    // For now require that we input the fully qualified input file name, including full path to file
    // **TODO** learn how to recover from art framework
    //fInputFileName = pset.get<std::string>("FullyQualifiedInputFile");
    auto const* geom = lar::providerFrom<geo::Geometry>(); // geometry is needed to go from OpChannel to OpDet 
    //fGeometry = lar::providerFrom<geo::Geometry>();
    fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
    //set the fhicl parameter here----------------------------------------------
    //fOpFlashModuleLabel      = pset.get<std::string> ("OpFlashModuleLabel");
    fOpFlashModuleLabel      = pset.get<std::string  > ("OpFlashModuleLabel", "simpleFlashBeam");
    //fOpFlashModuleLabel      = pset.get ("OpFlashModuleLabel", "opflash"); //In MCC7    
    fHitsModuleLabel         = pset.get< std::string > ("HitsModuleLabel", "gaushit");    
    fTrackModuleLabel        = pset.get< std::string > ("TrackModuleLabel",  "pandoraNu");    
    fVertexModuleLabel       = pset.get< std::string > ("VertexModuleLabel",  "pandoraNu");    
    fPandoraNuVertexModuleLabel=pset.get<std::string > ("PandoraNuVertexModuleLabel", "pandoraNu");
    fGenieGenModuleLabel     = pset.get< std::string > ("GenieGenModuleLabel", "generator");        
    fG4ModuleLabel           = pset.get< std::string > ("G4ModuleLabel", "largeant");
    //fCalorimetryModuleLabel  = pset.get< std::vector<std::string> >("CalorimetryModuleLabel");
    fCalorimetryModuleLabel  = pset.get< std::string > ("CalorimetryModuleLabel", "pandoraNucalo");

    //fShowerModuleLabel       = pset.get< std::string > ("ShowerModuleLabel", "mcreco");
    fShowerModuleLabel       = pset.get< std::string > ("ShowerModuleLabel", "showerrecopandora");


    fTrackMCSFitLabel        = pset.get< std::string > ("TrackMCSFitLabel", "pandoraNuMCSMu" );



    fDistToEdgeX             = geom->DetHalfWidth()   - pset.get("DistToEdgeX",   10.);    
    fDistToEdgeY             = geom->DetHalfHeight()  - pset.get("DistToEdgeY",   20.);    
    fDistToEdgeZ             = geom->DetLength() / 2. - pset.get("DistToEdgeZ",   10.);        
    
    fFlashWidth              = pset.get      ("FlashWidth", 80.);    

    fBeamMin                 = pset.get      ("BeamMin", 3.2);   //BNB+COSMIC
    fBeamMax                 = pset.get      ("BeamMax", 4.8);   //BNB+COSMIC



    //fBeamMin                 = pset.get      ("BeamMin", 3.65);   //extbnb 
    //fBeamMax                 = pset.get      ("BeamMax", 5.25);   //extbnb

    //fBeamMin                 = pset.get      ("BeamMin", 3.3);   //bnb 
    //fBeamMax                 = pset.get      ("BeamMax", 4.9);   //bnb
  
    fPEThresh                = pset.get      ("PEThresh", 50.);    
    fMinTrk2VtxDist          = pset.get      ("MinTrk2VtxDist", 5.);    
    fMinTrackLen             = pset.get      ("MinTrackLen", 75.);
    fG4minE                  = pset.get      ("G4minE",0.01); 

    _svm_x = pset.get<std::vector<double>> ("SVM_X");
 

    //Get the tool for MC truth Matching
    const fhicl::ParameterSet& truthParams = pset.get<fhicl::ParameterSet>("MCTruthMatching");

    if (truthParams.get<std::string>("tool_type") == "AssociationsTruth")
    {
        fMCTruthMatching = std::unique_ptr<truth::IMCTruthMatching>(new truth::AssociationsTruth(truthParams));
    }
    else
    {
        fMCTruthMatching = std::unique_ptr<truth::IMCTruthMatching>(new truth::BackTrackerTruth(truthParams));
    }
    /*
    */
    
    fhg4parpdg = new std::vector<int>;
    fhg4parstatus= new std::vector<int>;
    fhg4parpx= new std::vector<float>;
    fhg4parpy= new std::vector<float>;
    fhg4parpz= new std::vector<float>;
    fhg4parp = new std::vector<float>;
    fhg4partheta= new std::vector<float>;
    fhg4parphi= new std::vector<float>;
    
    trackidpcand = new std::vector<int>;
    trackstartxpcand = new std::vector<float>;
    trackstartypcand = new std::vector<float>;
    trackstartzpcand = new std::vector<float>;
    trackendxpcand = new std::vector<float>;
    trackendypcand = new std::vector<float>;
    trackendzpcand = new std::vector<float>;


    trackmompcand = new std::vector<double>;
    trackphipcand = new std::vector<double>;
    trackthetapcand = new std::vector<double>;
    tracklengthpcand = new std::vector<double>;
    tracktrunmeanpcand = new std::vector<double>;

    trackpidapcand= new std::vector<double>;
    
    fHitNucP4 = new TLorentzVector(-999,-999,-999,-999);


    //--------------------------------------------------------------------------------
    return;
}
bool CC1uNPSelAna::inFV(double x, double y, double z) const
{
    auto const* geom = lar::providerFrom<geo::Geometry>(); // geometry is needed to go from OpChannel to OpDet 


    double distInX = x - geom->DetHalfWidth();
    double distInY = y;
    double distInZ = z - 0.5 * geom->DetLength();
    
    if (fabs(distInX) < fDistToEdgeX && fabs(distInY) < fDistToEdgeY && fabs(distInZ) < fDistToEdgeZ) return true;
    
    return false;
}



double CC1uNPSelAna::GetDqDxTruncatedMean(std::vector<art::Ptr<anab::Calorimetry>> calos) {
  double result = -9999;
  double n = 1.;

  for (auto c : calos) {
    if (!c) continue;
    if (!c->PlaneID().isValid) continue;
    int planenum = c->PlaneID().Plane;
    if (planenum != 2) continue;
   
    std::vector<double> dqdx_v = c->dQdx(); 

    if (dqdx_v.size() == 0)
      return result;

    //for (auto q : dqdx_v) {
    //  std::cout << "dqdx before trim: " << q << std::endl;
    //}
    
    double median = GetMedian(dqdx_v);
    double std    = GetSTD(dqdx_v);
    
    std::cout << "median " << median << std::endl;
    std::cout << "std    " << std << std::endl;
    
    std::vector<double> dqdx_v_trimmed;
    dqdx_v_trimmed.clear();
    
    for (auto q : dqdx_v) {
      if (q > median - n * std && 
          q < median + n * std) {
         dqdx_v_trimmed.emplace_back(q);
         //std::cout << "dqdx after trim: " << q << std::endl;
      }
    }
    
    result = GetMean(dqdx_v_trimmed);
  }
    
  return result;




}

double CC1uNPSelAna::GetMean(std::vector<double> dqdx_v) {

  double mean = -9999;
  size_t size = dqdx_v.size();

  if (size == 0)
    return mean;

  double sum = 0;

  for (auto v : dqdx_v) {
    sum += v;
  }

  mean = sum/(double)size;

  return mean;
}

double CC1uNPSelAna::GetMedian(std::vector<double> dqdx_v) {

  double median = -9999;
  size_t size = dqdx_v.size();

  if (size == 0)
    return median;

  std::sort(dqdx_v.begin(), dqdx_v.end());
  if (size % 2 == 0){
    median = (dqdx_v[size/2 - 1] + dqdx_v[size/2]) / 2;
  }
  else{
    median = dqdx_v[size/2];
  }

  return median;
}
double CC1uNPSelAna::GetVariance(std::vector<double> dqdx_v) {

  double variance = -1;

  double sum = 0;
  double sum2 = 0;
  size_t size = dqdx_v.size();

  if (size == 0)
    return variance;

  for (auto value : dqdx_v) {

    sum  += value;
    sum2 += value*value;

  }  

  variance = sum2/(double)size - (sum/(double)size)*(sum/(double)size);

  return variance;

}

double CC1uNPSelAna::GetSTD(std::vector<double> dqdx_v) {

  if (dqdx_v.size() == 0)
    return -9999;

  double variance = GetVariance(dqdx_v);
  if (variance > 0)
    return std::sqrt(variance);
  else 
    return -9999;

}


bool CC1uNPSelAna::MIPConsistency(double dqds, double length) {

    if (length > 1000)
      return true;

    if (length < 0) {
      std::cout << "[MuonCandidateFinder] Track length < 0?!" << std::endl;
      return false;
    }

    int l = std::round(length);
    double dqds_cut = _svm_x.at(l);

    std::cout << "[MuonCandidateFinder] Track length is " << length << ", dqds_cut is " << dqds_cut << ", dqds value is " << dqds << std::endl;
 
    if (dqds*198 <= dqds_cut)
      return true;
  

    return false;

}


double CC1uNPSelAna::GetTrackRange(art::Ptr<recob::Track> InputTrackPtr) const
{
    return std::sqrt( std::pow(InputTrackPtr->Vertex().x()-InputTrackPtr->End().x(),2) + 
                      std::pow(InputTrackPtr->Vertex().y()-InputTrackPtr->End().y(),2) +
                      std::pow(InputTrackPtr->Vertex().z()-InputTrackPtr->End().z(),2) );    
}
double CC1uNPSelAna::GetTrackLength(art::Ptr<recob::Track> InputTrackPtr) const
{   
    return InputTrackPtr->Length();
}
int GetNhitsTrk(art::Ptr<recob::Track> InputTrackPtr) 
{   int NhitsTrk=-1;
    //art::ServiceHandle<cheat::BackTracker> bt;  //need this or not?
      
    //NhitsTrk=Hits.size();

    return NhitsTrk;    
}
double CC1uNPSelAna::GetFlashTrackDist(double flash, double start, double end) const
{
    if (end >= start) {
        if (flash < end && flash > start) return 0;
        else return std::min(fabs(flash-start), fabs(flash-end));
    }
    else {
        if (flash > end && flash < start) return 0;
        else return std::min(fabs(flash-start), fabs(flash-end));
    }
}
double CC1uNPSelAna::GetFlashTrackDistMod(double flash, double start, double end) const
{
     return std::min(fabs(flash-start), fabs(flash-end));
}
void CC1uNPSelAna::truthMatcher( std::vector<art::Ptr<recob::Hit>> all_hits, std::vector<art::Ptr<recob::Hit>> track_hits, const simb::MCParticle *&MCparticle, double &Efrac, double &Ecomplet)
{
  //art::ServiceHandle<cheat::BackTracker> bt;	
  std::map<int,double> trkID_E;	
  for(size_t j = 0; j < track_hits.size(); ++j)
  {	
    art::Ptr<recob::Hit> hit = track_hits[j];
    //const auto& hit = *track_hits[j];
    //std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);
    std::vector<sim::TrackIDE> TrackIDs = fMCTruthMatching->HitToTrackID(hit);
    for(size_t k = 0; k < TrackIDs.size(); k++)
    {	
      trkID_E[TrackIDs[k].trackID] += TrackIDs[k].energy;
    }
  }

  double E_em =0.0;	
  double max_E = -999.0;	
  double total_E = 0.0;	
  int TrackID = -999;	
  double partial_E =0.0; // amount of energy deposited by the particle that deposited more energy... tomato potato... blabla	
  //!if the collection of hits have more than one particle associate save the particle w/ the highest energy deposition 	
  //!since we are looking for muons/pions/protons this should be enough 	
  if( !trkID_E.size() ) 
  {
    MCparticle = 0;	
    return; //Ghost track???	
  }	
  for(std::map<int,double>::iterator ii = trkID_E.begin(); ii!=trkID_E.end(); ++ii)
  {	
    total_E += ii->second;	
    if((ii->second)>max_E)
    {	
      partial_E = ii->second;
      max_E = ii->second;
      TrackID = ii->first;
      if( TrackID < 0 ) E_em += ii->second;
    }	
  }
   	
  //MCparticle = bt->TrackIDToParticle(TrackID);		
  MCparticle = fMCTruthMatching->TrackIDToParticle(TrackID);		
  //In the current simulation, we do not save EM Shower daughters in GEANT. But we do save the energy deposition in TrackIDEs. If the energy deposition is from a particle that is the daughter of 	
  //an EM particle, the negative of the parent track ID is saved in TrackIDE for the daughter particle	
  //we don't want to track gammas or any other EM activity 	
  if( TrackID < 0 ) return;		
  //Efrac = (partial_E+E_em)/total_E;	
  Efrac = (partial_E)/total_E;		
  //completeness	
  double totenergy =0;	
  for(size_t k = 0; k < all_hits.size(); ++k)
  {	  
    art::Ptr<recob::Hit> hit = all_hits[k];	
    //std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackID(hit);	
    std::vector<sim::TrackIDE> TrackIDs = fMCTruthMatching->HitToTrackID(hit);	
    for(size_t l = 0; l < TrackIDs.size(); ++l)
    {	
      if(TrackIDs[l].trackID==TrackID) totenergy += TrackIDs[l].energy;	
    }	
  } 	
  Ecomplet = partial_E/totenergy;
}
double CC1uNPSelAna::driftedLength(const sim::MCTrack& mctrack, TLorentzVector& tpcstart, TLorentzVector& tpcend, TLorentzVector& tpcmom){
  // Get geometry.
  auto const* geom = lar::providerFrom<geo::Geometry>();
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

  //compute the drift x range
  double vDrift = detprop->DriftVelocity()*1e-3; //cm/ns
  double xrange[2] = {detprop->ConvertTicksToX(0,0,0,0),detprop->ConvertTicksToX(detprop->NumberTimeSamples(),0,0,0)};
  
  // Get active volume boundary.
  double bnd[6] = {0.,2.*geom->DetHalfWidth(),-geom->DetHalfHeight(),geom->DetHalfHeight(),0.,geom->DetLength()};

  double result = 0.;
  TVector3 disp;
  bool first = true;

  for(auto step: mctrack) {
    // check if the particle is inside a TPC
    if (step.X() >= bnd[0] && step.X() <= bnd[1] && step.Y() >= bnd[2] && step.Y() <= bnd[3] && step.Z() >= bnd[4] && step.Z() <= bnd[5]){
      // Doing some manual shifting to account for
      // an interaction not occuring with the beam dump
      // we will reconstruct an x distance different from
      // where the particle actually passed to to the time
      // being different from in-spill interactions
      double newX = step.X()+(step.T()*vDrift);
      if (newX < xrange[0] || newX > xrange[1]) continue;
     
      TLorentzVector pos(newX,step.Y(),step.Z(),step.T());
      if(first){
	tpcstart = pos;
	tpcmom = step.Momentum();
	first = false;
      }
      else {
	disp -= pos.Vect();
	result += disp.Mag();
      }
      disp = pos.Vect();
      tpcend = pos;
    }
  }
  return result;
}
// Length of MC particle, trajectory by trajectory (with the manual shifting for x correction)
double CC1uNPSelAna::driftedLength(const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi)
{
  // Get geometry.
  auto const* geom = lar::providerFrom<geo::Geometry>();
  auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  
  //compute the drift x range
  double vDrift = detprop->DriftVelocity()*1e-3; //cm/ns
  double xrange[2] = {detprop->ConvertTicksToX(0,0,0,0),detprop->ConvertTicksToX(detprop->NumberTimeSamples(),0,0,0)};
  
  // Get active volume boundary.
  double bnd[6] = {0.,2.*geom->DetHalfWidth(),-geom->DetHalfHeight(),geom->DetHalfHeight(),0.,geom->DetLength()};

  double result = 0.;
  TVector3 disp;
  bool first = true;

  for(unsigned int i = 0; i < p.NumberTrajectoryPoints(); ++i) {
    // check if the particle is inside a TPC
    if (p.Vx(i) >= bnd[0] && p.Vx(i) <= bnd[1] && p.Vy(i) >= bnd[2] && p.Vy(i) <= bnd[3] && p.Vz(i) >= bnd[4] && p.Vz(i) <= bnd[5]){
      // Doing some manual shifting to account for
      // an interaction not occuring with the beam dump
      // we will reconstruct an x distance different from
      // where the particle actually passed to to the time
      // being different from in-spill interactions
      double newX = p.Vx(i)+(p.T(i)*vDrift);
      if (newX < xrange[0] || newX > xrange[1]) continue;
      TLorentzVector pos(newX,p.Vy(i),p.Vz(i),p.T());
      if(first){
	start = pos;
	starti=i;
	first = false;
      }
      else {
	disp -= pos.Vect();
	result += disp.Mag();
      }
      disp = pos.Vect();
      end = pos;
      endi = i;
    }
  }
  return result;
}



// Length of MC particle, trajectory by trajectory (with out the manual shifting for x correction)
double CC1uNPSelAna::Simlength(const simb::MCParticle& p, TLorentzVector& start, TLorentzVector& end, unsigned int &starti, unsigned int &endi)
{
  // Get geometry.
  art::ServiceHandle<geo::Geometry> geom;
  
  // Get active volume boundary.
  double bnd[6] = {0.,2.*geom->DetHalfWidth(),-geom->DetHalfHeight(),geom->DetHalfHeight(),0.,geom->DetLength()};
  double result = 0.;
  TVector3 disp;
  bool first = true;

  for(unsigned int i = 0; i < p.NumberTrajectoryPoints(); ++i) {
    // check if the particle is inside a TPC
    if (p.Vx(i) >= bnd[0] && p.Vx(i) <= bnd[1] && p.Vy(i) >= bnd[2] && p.Vy(i) <= bnd[3] && p.Vz(i) >= bnd[4] && p.Vz(i) <= bnd[5]){
      if(first){
	start = p.Position(i);
	first = false;
	starti = i;
      }else{
	disp -= p.Position(i).Vect();
	result += disp.Mag();
      }
      disp = p.Position(i).Vect();
      end = p.Position(i);
      endi = i;
    }
  }
  return result;
}



double CC1uNPSelAna::Recolength(const recob::Track& track)
{
  double result = 0.;
  TVector3 disp = track.LocationAtPoint(0);
  int n = track.NumberTrajectoryPoints();

  for(int i = 1; i < n; ++i) {
    const TVector3& pos = track.LocationAtPoint(i);
    disp -= pos;
    result += disp.Mag();
    disp = pos;
  }
  return result;
}




int CC1uNPSelAna::Topology(int nmuons, int nelectrons, int npions, int npi0, int nprotons, bool cosmicflag, bool OOFVflag)
{
  //// This function return the true topology of the event, numu & anti-numu                                                        
  ////1. CC0Pion0Proton                                                                                                             
  ////2. CC0Pion1Proton                                  
  ////3. CC0Pion2Proton                                                                                                             
  ////4. CC0PionNProton                                                                                                             
  ////5. CC1PionNProton (1 Pion= 1 charged pion || 1 neutral pion)
  ////6. CCNPionNProton    
  ////7. CCnue-antinue   
  ////8. NC                                                                                            
  ////9. OOFV (nu event out of FV)                                                                                                  
  ////10. Cosmic                                                                                                       
  ////11. Other (just in case, let's check!)                                                                        
  /// 12. 2&3&4 -> CC0PinProton (N>0) /// *** add this one!                                                                         
  /// e.g. numu CC inclusive= Topology >0 && Topology < 7 
  
  
  if (nmuons >0 && (nelectrons + npions + npi0 + cosmicflag + OOFVflag) == 0 && nprotons ==0 ) return 1;
  if (nmuons >0 && (nelectrons + npions + npi0 + cosmicflag + OOFVflag) == 0 && nprotons ==1 ) return 2;
  if (nmuons >0 && (nelectrons + npions + npi0 + cosmicflag + OOFVflag) == 0 && nprotons ==2 ) return 3;
  if (nmuons >0 && (nelectrons + npions + npi0 + cosmicflag + OOFVflag) == 0 && nprotons >2 ) return 4;
  if (nmuons >0 && (nelectrons + cosmicflag + OOFVflag) == 0 && (npions + npi0) == 1 ) return 5;
  if (nmuons >0 && (nelectrons + cosmicflag + OOFVflag) == 0 && (npions + npi0) > 1 ) return 6;
  if (nmuons == 0 && (cosmicflag + OOFVflag) == 0 && nelectrons > 0 ) return 7;
  if (nmuons==0 && (cosmicflag+ OOFVflag==0) && nelectrons ==0) return 8;
  if (OOFVflag) return 9;
  if (cosmicflag) return 10;  //check with colton how to select cosmic event

  else return 11;


}


/*int CC1uNPSelAna::Reaction(std::vector<art::Ptr<simb::MCTruth> > mclist){
  //// detail here the reaction numbering ****                                                                                                            

  int mode_truth = mclist[0]->GetNeutrino().Mode();/// only for the first neutrino                                                                        

  return mode_truth;
}
*/
double CC1uNPSelAna::GetDistTracktoVtx(art::Ptr<recob::Track> InputTrackPtr, TVector3 InputvertexPos){

  TVector3 trackPos = InputTrackPtr->Vertex();
  TVector3 trackEnd = InputTrackPtr->End();

  // Take the closer end---------------------------------                                                                                                 
  double diststart = (trackPos - InputvertexPos).Mag();
  double distend = (trackPos - InputvertexPos).Mag();

  if(diststart> distend ) return distend;
  else return diststart;

}

double CC1uNPSelAna::GetTrackDirection(art::Ptr<recob::Track> InputTrackPtr, TVector3 InputvertexPos){
// revisit with Libo's algorithm
// Gianluca                                                                                                                                                 

   // art::Ptr<recob::Track> OutputTrackPtr;                                                                                                               
   TVector3 trackPos = InputTrackPtr->Vertex();
   TVector3 trackEnd = InputTrackPtr->End();
   double trackTheta = InputTrackPtr->Theta();

   // Take the closer end---------------------------------                                                                                                 
   double diststart = (trackPos - InputvertexPos).Mag();
   double distend = (trackPos - InputvertexPos).Mag();

   if(diststart> distend ) {

     trackPos          = InputTrackPtr->End();
     trackEnd          = InputTrackPtr->Vertex();

     // *** re-calculate phi here /// *** missing phi!                                                                                                     
     trackTheta =  pi - InputTrackPtr->Theta();/// Libo's code                                                                                             

    }
    
 return trackTheta;
 //return trackPos, trackEnd, trackTheta;/// *** revisit and add phi!!!                                                                                    
}

void CC1uNPSelAna::ClearLocalData(){
  //std::fill(fhg4parpdg, fhg4parpdg + sizeof(fhg4parpdg)/sizeof(fhg4parpdg[0]), -99999.);
  //std::fill(fhg4parstatus, fhg4parstatus + sizeof(fhg4parstatus)/sizeof(fhg4parstatus[0]), -99999.);
  //std::fill(fhg4parp, fhg4parp + sizeof(fhg4parp)/sizeof(fhg4parp[0]), -99999.);
  //std::fill(fhg4parpx, fhg4parpx + sizeof(fhg4parpx)/sizeof(fhg4parpx[0]), -99999.);
  //std::fill(fhg4parpy, fhg4parpy + sizeof(fhg4parpy)/sizeof(fhg4parpy[0]), -99999.);
  //std::fill(fhg4parpz, fhg4parpz + sizeof(fhg4parpz)/sizeof(fhg4parpz[0]), -99999.);
  //std::fill(fhg4partheta, fhg4partheta + sizeof(fhg4partheta)/sizeof(fhg4partheta[0]), -99999.);
  //std::fill(fhg4parphi, fhg4parphi + sizeof(fhg4parphi)/sizeof(fhg4parphi[0]), -99999.);
  std::fill(vectex,       vectex+sizeof(vectex)/sizeof(vectex[0]), -9999.);
}

/*
double CC1uNPSelAna::GetTruthNuEnergy(std::vector<art::Ptr<simb::MCTruth> > mclist){/// *** need to be check                                              
  double enu_truth = mclist[0]->GetNeutrino().Nu().E();/// only for the first neutrino                                                                    
  return enu_truth;
}

double CC1uNPSelAna::GetTruthInvMass(std::vector<art::Ptr<simb::MCTruth> > mclist){
  double W_truth = mclist[0]->GetNeutrino().W();/// only for the first neutrino                                                                           
  return W_truth;
}

double CC1uNPSelAna::GetTruthQ2(std::vector<art::Ptr<simb::MCTruth> > mclist){
  double Q2_truth = mclist[0]->GetNeutrino().QSqr();/// only for the first neutrino                                                                       
  return Q2_truth;
}
*/
bool CC1uNPSelAna::inTPC(double x, double y, double z) const // *** need to revisit values                                                                
{  //TPC dimensions neeed to be checked

  if(x<-50 || x>300) return false;
  if(y<-106 || y>116) return false;
  if(z<0.01 || z>1050) return false;
  else return true;

}




void CC1uNPSelAna::produces(art::EDProducer* owner)
{
    fMyProducerModule = owner;
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::Track> >();
    fMyProducerModule->produces< art::Assns<recob::Vertex, recob::PFParticle> >();
}


    
/*void CC1uNPSelAna::respondToOpenInputFile(art::FileBlock const& fileBlock)
{
    // Override the fhicl parameter for the input file name
    fInputFileName = fileBlock.fileName();

    return;
}
    
void CC1uNPSelAna::respontToOpenOutputFile(art::FileBlock const& fileBlock)
{
    // TODO
    return;
}
*/
//-----------------------------------------------------------------------
void CC1uNpSelAnaDataStruct::FlashDataStruct::Resize(size_t nFlashes)
{
  MaxFlashes = nFlashes;

  flsTime.resize(MaxFlashes);
  flsPe.resize(MaxFlashes);
  flsPePerOpDet.resize(MaxFlashes);
  flsXcenter.resize(MaxFlashes);
  flsYcenter.resize(MaxFlashes);
  flsZcenter.resize(MaxFlashes);
  flsYwidth.resize(MaxFlashes);
  flsZwidth.resize(MaxFlashes);
  flsTwidth.resize(MaxFlashes);

  size_t nOpDets = this->GetNOpDet();
  if (kNOpDets != nOpDets) {
    mf::LogError("AnalysisTree") 
        << "Number of optical detectors from geometry services is" << nOpDets 
        << ", which is different from the one expected of " << kNOpDets 
        << ". Check variale kNOpDets in analysis tree module.";
  }

  //for (size_t opfls = 0; opfls < MaxFlashes; opfls++) {
  //  FillWith(flsPePerOpDet[opfls]   , -99999.);
  //}
}
/*
void CC1uNPSelAnaDataStruct::FlashDataStruct::Clear() {
  Resize(MaxFlashes);
  nfls = -9999;

  FillWith(flsTime        , -9999  );
  FillWith(flsPe          , -9999  );
  FillWith(flsXcenter     , -9999  );
  FillWith(flsYcenter     , -9999  );
  FillWith(flsZcenter     , -9999  );
  FillWith(flsYwidth      , -9999  );
  FillWith(flsZwidth      , -9999  );
  FillWith(flsTwidth      , -9999  );

  for (size_t opfls = 0; opfls < MaxFlashes; opfls++) {
    FillWith(flsPePerOpDet[opfls]  , -9999  );
  }
  
}

void CC1uNPSelAna::FlashDataStruct::Clear() {
  Resize(MaxFlashes);
  nfls = -9999;

  FillWith(flsTime        , -9999  );
  FillWith(flsPe          , -9999  );
  FillWith(flsXcenter     , -9999  );
  FillWith(flsYcenter     , -9999  );
  FillWith(flsZcenter     , -9999  );
  FillWith(flsYwidth      , -9999  );
  FillWith(flsZwidth      , -9999  );
  FillWith(flsTwidth      , -9999  );

  for (size_t opfls = 0; opfls < MaxFlashes; opfls++) {
    FillWith(flsPePerOpDet[opfls]  , -9999  );
  }
}
*/






//---------------------------------------------------------------------
void  CC1uNPSelAna::analyze(const art::Event& event)
{

    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    //if this is MC truth add more variables=======================
    /* 
    */
    //=============================================================
    fMC_allsel->Fill();
    //recover the back tracker
    //art::ServiceHandle<cheat::BackTracker> bt;

    // "Rebuild" the maps used by the parallel backtracker
    //fMCTruthMatching->Rebuild(event);
    
    
    std::cout<<"The Event Number is: "<<fEvent<<std::endl;
    std::cout<<"The Run Number is: "<<fRun<<std::endl;
    std::cout<<"The SubRun Number is: "<<fSubRun<<std::endl;


    bool isMC = !event.isRealData();

    //std::cout<<"is this an MC event or not?????????   "<<isMC<<std::endl;

    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     
    //use handles to connect the sim::Particles
    Int_t nGeniePrimaries=0;
    //Int_t nGeantPrimaries=0;
    Int_t nGeantParticles=0;
    if(isMC){ 

    // "Rebuild" the maps used by the parallel backtracker
    fMCTruthMatching->Rebuild(event);
 
    //comment out for test of this SCE
    //auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
    auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    //GENIE-----------------------------------------------   
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    std::vector<art::Ptr<simb::MCTruth> > mclist;
    if (event.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

    Nevt_truth=mclist.size();
    art::Ptr<simb::MCTruth> mctruth; 
    mctruth = mclist[0];

    if (mctruth->NeutrinoSet()) nGeniePrimaries = mctruth->NParticles();
    std::cout<<"<<<<total number of GENIE particle is :"<<nGeniePrimaries<<std::endl;
    //get the MC truth here 
    if(mctruth->NeutrinoSet()&& mctruth->Origin()==simb::kBeamNeutrino) 
    {
    _fTrueccnc=mctruth->GetNeutrino().CCNC();  //0->CC, 1->NC
    _fTruemode=mctruth->GetNeutrino().Mode();   //0=Quasi-elastic or Elastic, 1=Resonant (RES), 2=DIS, 3=Coherent productio
    //neutrino PDG code (nue=12; anti-nue=-12; numu=14; anti-numu=-14; nutau=16; anti-nutau=-16)_fTrueinttype=mctruth->GetNeutrino().InteractionType();
    _fTrueinttype=mctruth->GetNeutrino().InteractionType();
    _fTruenupdg=mctruth->GetNeutrino().Nu().PdgCode();
    _fTrueenu=mctruth->GetNeutrino().Nu().E(); //true neutrino energy
    _fTrueq2truth=mctruth->GetNeutrino().QSqr(); //true Q2
    _fTrueWtruth =mctruth->GetNeutrino().W();
    _fTrueXtruth =mctruth->GetNeutrino().X();
    _fTrueYtruth =mctruth->GetNeutrino().Y();

    trueMuonTrueMomentum = mctruth->GetNeutrino().Lepton().P();
    trueMuonTrueTheta = mctruth->GetNeutrino().Lepton().Momentum().Theta();
    trueMuonTruePhi = mctruth->GetNeutrino().Lepton().Momentum().Phi();
    
    trueProtonsTrueMomentum->clear();
    trueProtonsTrueTheta->clear();
    trueProtonsTruePhi->clear();
    for (int igeniepart(0); igeniepart<nGeniePrimaries; igeniepart++){
      simb::MCParticle part = mctruth->GetParticle(igeniepart);
      if (part.PdgCode()==2212 && part.StatusCode()==1){
        trueProtonsTrueMomentum->push_back(part.P());
        trueProtonsTrueTheta->push_back(part.Momentum().Theta());
        trueProtonsTruePhi->push_back(part.Momentum().Phi());
      }
    }
    /// Also here we should get things like the true struck neutron momentum - I think we need a GTruth object for this
    art::ValidHandle< std::vector<simb::GTruth> > gtruth = event.getValidHandle< std::vector<simb::GTruth> >("generator");
    if (gtruth->size() <1){
      std::cout << "WARNING NO GTRUTH OBJECT" << std::endl;
    }
    else{
      TLorentzVector v_tmp = gtruth->at(0).fHitNucP4;
      fHitNucP4->SetXYZT(v_tmp.X(), v_tmp.X(), v_tmp.Y(), v_tmp.E());
    }

    _fTruenuvrtxx=mctruth->GetNeutrino().Nu().Vx(); //true vertex x
    _fTruenuvrtxy=mctruth->GetNeutrino().Nu().Vy(); //true vertex y
    _fTruenuvrtxz=mctruth->GetNeutrino().Nu().Vz(); //true vertex z
    //get the vertex after the SCE correction
    _fTruenuvrtxx_SCE=  mctruth->GetNeutrino().Nu().Vx() - SCE->GetPosOffsets(mctruth->GetNeutrino().Nu().Vx(),mctruth->GetNeutrino().Nu().Vy(),mctruth->GetNeutrino().Nu().Vz())[0];
    _fTruenuvrtxy_SCE=  mctruth->GetNeutrino().Nu().Vy() + SCE->GetPosOffsets(mctruth->GetNeutrino().Nu().Vx(),mctruth->GetNeutrino().Nu().Vy(),mctruth->GetNeutrino().Nu().Vz())[1];
    _fTruenuvrtxz_SCE = mctruth->GetNeutrino().Nu().Vz() + SCE->GetPosOffsets(mctruth->GetNeutrino().Nu().Vx(),mctruth->GetNeutrino().Nu().Vy(),mctruth->GetNeutrino().Nu().Vz())[2];

    //std::cout<<"the true neutrino vertex position is : "<<_fTruenuvrtxx<<" "<<_fTruenuvrtxy<<" "<<_fTruenuvrtxz<<std::endl;
    //std::cout<<"test space charge effect correction x " <<SCE->GetPosOffsets(mctruth->GetNeutrino().Nu().Vx(),mctruth->GetNeutrino().Nu().Vy(),mctruth->GetNeutrino().Nu().Vz())[0]<<std::endl;
    //std::cout<<"test space charge effect correction y " <<SCE->GetPosOffsets(mctruth->GetNeutrino().Nu().Vx(),mctruth->GetNeutrino().Nu().Vy(),mctruth->GetNeutrino().Nu().Vz())[1]<<std::endl;
    //std::cout<<"test space charge effect correction z " <<SCE->GetPosOffsets(mctruth->GetNeutrino().Nu().Vx(),mctruth->GetNeutrino().Nu().Vy(),mctruth->GetNeutrino().Nu().Vz())[2]<<std::endl;


    //loop over all the particles from GENIE Stage and try to connect to GEANT4------------------
    //Need to do the space charge effect correction 






       //-------------------------------------------------------------------------------------------
    }

    

    std::cout<<"<<<<<<<<<<<<<<<libo test after GENIE Stage<<<<<<<<<<<<<<<<<"<<std::endl;





    //const art::Ptr<simb::MCTruth> mctruth = bt->TrackIDToMCTruth(mcparticle.TrackId());

    //std::cout<<"CC or NC"<<_fTrueccnc<<std::endl;
    //std::cout<<"interaction mode is "<<_fTruemode<<std::endl;
    //std::cout<<"interration type "<<_fTrueinttype<<std::endl;
    //std::cout<<"neutrino pdg "<<_fTruenupdg<<std::endl;
    fMC_Truth->Fill();
    //GEANT-------------------------------------------------
    art::Handle< std::vector<simb::MCParticle> > mcparticleListHandle;
    //art::ValidHandle< std::vector<simb::MCParticle> > mcparticleListHandle =
    //event.getValidHandle< std::vector<simb::MCParticle> >(fG4ModuleLabel);
    std::vector< art::Ptr<simb::MCParticle> > mcparticlelist;
    if (event.getByLabel(fG4ModuleLabel,mcparticleListHandle))
    art::fill_ptr_vector(mcparticlelist, mcparticleListHandle);


    //check if the truth vertex within the FV or not

   
    //helper to map the track id index
    std::map<int, size_t> TrackIDtoIndex;
    std::vector<int> gpdg;   //vector to store the pdg id of geant4
    std::vector<int> gmother; //vector to store the mother id of geant4 


    Int_t nmuons=0;
    Int_t npions=0;
    Int_t npi0=0;
    Int_t nprotons=0;
    Int_t nprotons_200thresh=0;
    Int_t nprotons_300thresh=0;
    Int_t nprotons_400thresh=0;
    Int_t nelectrons=0;
    Bool_t cosmicflag=false; //cosmic event
    Bool_t OOFVflag=false; //nu event outside FV
    //auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

    fhg4parpdg->clear();
    fhg4parstatus->clear();
    fhg4parpx->clear();
    fhg4parpy->clear();
    fhg4parpz->clear();
    fhg4parp->clear();
    fhg4partheta->clear();
    fhg4parphi->clear();

    std::cout<<"start looping over all the geant particles and get the topology"<<std::endl;

    std::string pri("primary");


    nGeantParticles=mcparticleListHandle->size();
    nGEANTparticles=nGeantParticles;
    std::cout<<"total number of geant particle is "<<nGeantParticles<<std::endl;
    for(size_t g4pt=0; g4pt<mcparticleListHandle->size(); g4pt++){
      if(mcparticleListHandle.isValid() && mcparticleListHandle->size()>0){
      //const art::Ptr<simb::MCParticle> pPart(mcparticleListHandle,g4pt);
      simb::MCParticle const & pPart = mcparticleListHandle->at(g4pt);

      //use a back tracker to get the mc truth of a certain G4 particles
      //const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth(pPart.TrackId());
      //new Truth matching test      
      const art::Ptr<simb::MCTruth> mc_truth=fMCTruthMatching->TrackIDToMCTruth(pPart.TrackId());   
       



      //++geant_particle;
      gpdg.push_back(pPart.PdgCode());          //get the PDG ID of all the particles
      gmother.push_back(pPart.Mother());        //get the mother particles ID
      //if (pPart->E()<fG4minE &&(!isPrimary)) continue;  //select the particles with a minimum GEANT 4 energy cut and primary particles
      if (pPart.E()<fG4minE ) continue;       


      Bool_t isPrimary = pPart.Process() == pri;  
      _fg4processname.push_back(pPart.Process());


      _fg4Mother.push_back(pPart.Mother());

      _fg4TrackId.push_back(pPart.TrackId());

      _fg4parpdg.push_back(pPart.PdgCode());
      _fg4parstatus.push_back(pPart.StatusCode());
      _fg4EngS.push_back(pPart.E());
      _fg4EngE.push_back(pPart.EndE());
      _fg4Mass.push_back(pPart.Mass());
      _fg4px.push_back(pPart.Px());
      _fg4py.push_back(pPart.Py());
      _fg4pz.push_back(pPart.Pz());

      _fg4P.push_back(pPart.Momentum().Vect().Mag());
      _fg4StartPointx.push_back(pPart.Vx());
      _fg4StartPointy.push_back(pPart.Vy());
      _fg4StartPointz.push_back(pPart.Vz());




      _fg4EndPointx.push_back(pPart.EndPosition()[0]);
      _fg4EndPointy.push_back(pPart.EndPosition()[1]);
      _fg4EndPointz.push_back(pPart.EndPosition()[2]);
      //do the space charge effect correction with the track start, end point position
      _fg4SCEcorrStartPointx.push_back(pPart.Vx() - SCE->GetPosOffsets(pPart.Vx(),pPart.Vy(),pPart.Vz())[0]);
      _fg4SCEcorrStartPointy.push_back(pPart.Vy() + SCE->GetPosOffsets(pPart.Vx(),pPart.Vy(),pPart.Vz())[1]);
      _fg4SCEcorrStartPointz.push_back(pPart.Vz() + SCE->GetPosOffsets(pPart.Vx(),pPart.Vy(),pPart.Vz())[2]);



      _fg4SCEcorrEndPointx.push_back(pPart.EndPosition()[0] - SCE->GetPosOffsets(pPart.EndPosition()[0],pPart.EndPosition()[1],pPart.EndPosition()[2])[0]);
      _fg4SCEcorrEndPointy.push_back(pPart.EndPosition()[1] + SCE->GetPosOffsets(pPart.EndPosition()[0],pPart.EndPosition()[1],pPart.EndPosition()[2])[1]);
      _fg4SCEcorrEndPointz.push_back(pPart.EndPosition()[2] + SCE->GetPosOffsets(pPart.EndPosition()[0],pPart.EndPosition()[1],pPart.EndPosition()[2])[2]);
      


 
     //std::cout<<"libo test 0 "<<std::endl;



      _fg4EndT.push_back(  pPart.EndT());
      _fg4theta.push_back( pPart.Momentum().Theta());
      _fg4phi.push_back( pPart.Momentum().Phi());
      _fg4theta_xz.push_back( std::atan2(pPart.Px(), pPart.Pz()));
      _fg4theta_yz.push_back( std::atan2(pPart.Py(), pPart.Pz()));
      _fg4origin.push_back(mc_truth->Origin());  //what created particle
      _fg4MCTruthIndex.push_back(mc_truth.key());
      _fg4NumberDaughters.push_back(pPart.NumberDaughters());
      //_fg4MCTruthIndex_new.push_back(mc_truth_new.key());
      //------------------------------------------------------------------------
      /*
      fhg4parpdg[g4pt]=pPart.PdgCode();
      fhg4parstatus[g4pt]=pPart.StatusCode();
      fhg4parpx[g4pt]=pPart.Px();
      fhg4parpy[g4pt]=pPart.Py();
      fhg4parpz[g4pt]=pPart.Pz();
      fhg4partheta[g4pt]=pPart.Momentum().Theta();
      fhg4parphi[g4pt]=pPart.Momentum().Phi();
      fhg4parp[g4pt]=pPart.Momentum().Vect().Mag();
      */
      
      
      //=========================================================================== 
      if( isPrimary && mctruth->NeutrinoSet() && pPart.StatusCode()==1 && pPart.Mother()==0 && mc_truth->Origin()== simb::kBeamNeutrino)  
      //primary tells you if the particle is from Michel electron or decay of other particle
      { 
       
 
        fhg4parpdg->push_back(pPart.PdgCode());
        fhg4parstatus->push_back(pPart.StatusCode());
        fhg4parpx->push_back(pPart.Px());
        fhg4parpy->push_back(pPart.Py());
        fhg4parpz->push_back(pPart.Pz());
        fhg4partheta->push_back(pPart.Momentum().Theta());
        fhg4parphi->push_back(pPart.Momentum().Phi());
        fhg4parp->push_back(pPart.Momentum().Vect().Mag());

 

        if(_fTrueccnc==0 && abs(pPart.PdgCode())==13) {
          nmuons=nmuons+1;
          //of all the muons select the primary muon from neutrino interaction
          //std::cout<<"Mother of the muon is "<<pPart.Mother()<<std::endl;
          //std::cout<<"TrackId of the muon is "<<pPart.TrackId()<<std::endl;
          //std::cout<<"Origin of this muon: "<<mc_truth->Origin()<<std::endl;

          //------------------------------------------------------------------
        }
        if(_fTrueccnc==0 && abs(pPart.PdgCode())==211) {npions=npions+1;}
        if(_fTrueccnc==0 && abs(pPart.PdgCode())==111) {npi0=npi0+1;}
        if(_fTrueccnc==0 && abs(pPart.PdgCode())==11) {nelectrons=nelectrons+1;}
        if(_fTrueccnc==0 && abs(pPart.PdgCode())==2212) {nprotons=nprotons+1;}
        if(_fTrueccnc==0 && abs(pPart.PdgCode())==2212 && pPart.Momentum().Vect().Mag() > 0.2 ) {nprotons_200thresh=nprotons_200thresh+1;}
        if(_fTrueccnc==0 && abs(pPart.PdgCode())==2212 && pPart.Momentum().Vect().Mag() > 0.3 ) {nprotons_300thresh=nprotons_300thresh+1;}
        if(_fTrueccnc==0 && abs(pPart.PdgCode())==2212 && pPart.Momentum().Vect().Mag() > 0.4 ) {nprotons_400thresh=nprotons_400thresh+1;}
      }
      //std::cout<<"libo test 1 "<<std::endl;
    
    }//end of is the g4 handle is valid and the size is greater than 0;
   }//end of loop over all the geant 4 particles
   //redefine the OOFV here
   if(!inFV(_fTruenuvrtxx, _fTruenuvrtxy, _fTruenuvrtxz) && (nmuons!=0 || nelectrons!=0 || npions!=0 || npi0!=0 || nprotons!=0)) {OOFVflag=true;}


   //check the cosmic flag before include the Topology function
   if(!inFV(_fTruenuvrtxx, _fTruenuvrtxy, _fTruenuvrtxz) && (nmuons==0 && nelectrons==0 && npions==0 && npi0==0 && nprotons==0)) {cosmicflag=true;} 

   Int_t TopFlag=Topology(nmuons, nelectrons, npions, npi0, nprotons, cosmicflag, OOFVflag);
   Int_t TopFlag200=Topology(nmuons, nelectrons, npions, npi0, nprotons_200thresh, cosmicflag, OOFVflag);
   Int_t TopFlag300=Topology(nmuons, nelectrons, npions, npi0, nprotons_300thresh, cosmicflag, OOFVflag);
   Int_t TopFlag400=Topology(nmuons, nelectrons, npions, npi0, nprotons_400thresh, cosmicflag, OOFVflag);


   truthtop=TopFlag;
   truthtop_200thresh=TopFlag200;
   truthtop_300thresh=TopFlag300;
   truthtop_400thresh=TopFlag400;
   //std::cout<<"Topology flag is: "<<TopFlag <<" OOFVflag= "<<OOFVflag <<std::endl;
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
   fNTruePTrks=nprotons;
   fMC_Geant->Fill();
   }//end of if MC


    // In principle we can have several producers running over various configurations of vertices and tracks.
    // The output associations we want to check are then encapsuated in the input vectors of strings
    // So the outer loop is over the indices
    // fVertexModuelLabelVec is the vector of vertex??????????????
    
    //auto const* geom = lar::providerFrom<geo::Geometry>(); // geometry is needed to go from OpChannel to OpDet 

    //Check flash tag-------------------------------------------------------------------------

    //------------------------------------------------------------------------------------ 
    // Agreed convention is to ALWAYS output to the event store so get a pointer to our collection 
    std::unique_ptr<art::Assns<recob::Vertex, recob::Track>>      vertexTrackAssociations(new art::Assns<recob::Vertex, recob::Track>);       
    std::unique_ptr<art::Assns<recob::Vertex, recob::PFParticle>> vertexPFParticleAssociations(new art::Assns<recob::Vertex, recob::PFParticle>);

    // Recover the hanles to the vertex and track collections we want to analyze.
    art::Handle<std::vector<recob::Vertex>>  vertexVecHandle;
    art::Handle<std::vector<recob::Track>>   trackVecHandle;
    art::Handle<std::vector<recob::MCSFitResult>>   MCSFitHandle;

    art::Handle<std::vector<recob::OpFlash>> flashListHandle;
    //art::Handle< std::vector<recob::Shower> > showerHandle; 
    //art::Handle< std::vector<sim::MCShower> > showerHandle;
                                        
    event.getByLabel(fVertexModuleLabel,    vertexVecHandle);
    event.getByLabel(fTrackModuleLabel,     trackVecHandle);
    event.getByLabel(fTrackMCSFitLabel,     MCSFitHandle);

    std::vector<art::Ptr<recob::MCSFitResult> >  mcsfitlist;
    art::fill_ptr_vector(mcsfitlist, MCSFitHandle);

    //std::vector<std::vector<recob::Shower> const*> showerList;
    art::ValidHandle< std::vector<recob::Shower> > showerHandle =
    event.getValidHandle< std::vector<recob::Shower> >(fShowerModuleLabel);
    std::vector< art::Ptr<recob::Shower> > showerList;
    //if(event.getByLabel( fShowerModuleLabel,  showerHandle))
    art::fill_ptr_vector(showerList, showerHandle);

    //std::vector<std::vector<recob::Shower> const*> showerList;
    //std::vector< art::Ptr<sim::MCShower> > showerList;
    //event.getByLabel( fShowerModuleLabel , showerHandle);
    //art::fill_ptr_vector(showerList, showerHandle);

    art::Handle< std::vector<recob::Hit> > hitListHandle;
    std::vector<art::Ptr<recob::Hit> > hitlist;

    if(event.getByLabel(fHitsModuleLabel,hitListHandle))
    art::fill_ptr_vector(hitlist, hitListHandle);

    //-----------------------------------------------------------------------------------------
    std::vector< std::vector<art::Ptr<recob::Vertex> > >     nuvertexlist(1);
    std::vector< std::vector<art::Ptr<recob::PFParticle> > > nuvertexlistToPfp(1);

    //get the neutrino vertex, but not the vertex of the tracks
    lar_pandora::VertexVector vertexVector;
    lar_pandora::PFParticlesToVertices particlesToVertices;
    //lar_pandora::LArPandoraHelper::CollectVertices(event, fPandoraNuVertexModuleLabelVec[0], vertexVector, particlesToVertices);  
    lar_pandora::LArPandoraHelper::CollectVertices(event, fPandoraNuVertexModuleLabel, vertexVector, particlesToVertices);  



    for (lar_pandora::PFParticlesToVertices::iterator it = particlesToVertices.begin(); it != particlesToVertices.end(); ++it){
      
      art::Ptr<recob::PFParticle> pfParticle = it->first;
      lar_pandora::VertexVector   vertex_v   = it->second;
      if (vertex_v.empty()) continue;
      if (lar_pandora::LArPandoraHelper::IsNeutrino(pfParticle)) {
        if (vertex_v.size() == 1) { // require 1 vtx associated to the neutrino PFP
          nuvertexlist[0].emplace_back(vertex_v[0]);
          nuvertexlistToPfp[0].emplace_back(pfParticle);
        }
      }
    }
    //art::Handle<std::vector<recob::Vertex>> nuvertexHandle;
    
    //if(event.getByLabel(fPandoraNuVertexModuleLabel, nuvertexHandle))

    //art::fill_ptr_vector(nuvertexlist[0], nuvertexHandle); 
    //---------------------------------------------------------------------------------------



    /*    //GENIE-----------------------------------------------   
    art::Handle< std::vector<simb::MCTruth> > mctruthListHandle;
    if (evt.getByLabel(fGenieGenModuleLabel,mctruthListHandle))
    art::fill_ptr_vector(mclist, mctruthListHandle);

    //GEANT-------------------------------------------------
    art::Handle< std::vector<simb::MCParticle> > mcparticleListHandle;
    std::vector< art::Ptr<simb::MCParticle> > mcparticlelist;
    //if (evt.getByLabel(fG4ModuleLabel,mcparticleListHandle))
    art::fill_ptr_vector(mcparticlelist, mcparticleListHandle);
   */



    //----------------------------------------------------
    std::vector<art::Ptr<recob::OpFlash> > flashlist;
                                                                
                                                                    
    if (event.getByLabel(fOpFlashModuleLabel,flashListHandle)) 
    art::fill_ptr_vector(flashlist, flashListHandle);  

    bool                  flashtag(false);   
    bool                  flashwintag(false);
   
    const   recob::OpFlash* flashPtr(0);
    flashmax = 0;
    fopflashmax=0;
 
    for(const auto& opFlash : flashlist)   //loop over all the flashes
    {
       //if (opFlash->Time() > fBeamMin && opFlash->Time() < fBeamMax && opFlash->TotalPE() > fPEThresh)
       if (opFlash->Time() > fBeamMin && opFlash->Time() < fBeamMax)
       {
         flashwintag = true;
         if(opFlash->TotalPE() > fPEThresh){
         flashtag=true;
         //std::cout<<"event number= "<<fEvent<<"flashtime is "<<opFlash->Time()<< std::endl;    
         // Keep track of the largest flash
         if (opFlash->TotalPE() > flashmax)              
         {

           fopflashtime=opFlash->Time(); //get the flash time corresponding to the maximum flash
           fopflashmax= opFlash->TotalPE();

           flashPtr = opFlash.get();
           flashmax = opFlash->TotalPE();
         }
         

           //fFlashPE->Fill(opFlash->TotalPE(), 1.);
           //fFlashTime->Fill(opFlash->Time(), 1.);
         }
       }
    } //end of loop over all the flash  
    
    if(flashwintag){
        fMC_flashwin->Fill(); 
    }
    
    //if(Nevt_truth==1)  //only for MC event
    {
    //std::cout<<"only one neutrino event"<<std::endl; 
    if(flashwintag && flashtag) { // if the flashtag=true
        fMC_flashtag->Fill();  
        //check the true vertex position see if they locate in the FV
        //std::cout<<"True Vertex Position"<<_fTruenuvrtxx<<" "<< _fTruenuvrtxy<<" "<< _fTruenuvrtxz<<std::endl;        
 
        //--------------------------------------------------------------------------------------------
        // Require valid handles, otherwise nothing to do
        if (vertexVecHandle.isValid() && vertexVecHandle->size() > 0 && trackVecHandle.isValid() && trackVecHandle->size() > 0)
        {
           //std::cout<<"vertex handle is valide"<<std::endl;
           //Recover associations to PFparticles...
           art::FindManyP<recob::PFParticle> trackToPFPartAssns(trackVecHandle,  event, fTrackModuleLabel);
           //art::FindManyP<recob::PFParticle> pfp_from_track(track_h, e, _pfp_producer);
           //----loop over all the flashes and check if there are flashes within the beam
           //window and above the PE threshold
           //const recob::OpFlash* flashPtr(0);
           //get the vertex position and check if this vertex is within Fiducial Volume

            int    VertexCandidate=-1;
            int    TrackCandidate=-1;
            int    TrackProtonCandidate=-1;
            
  

            double TrackCandLength=0;
            double TrackProtonCandLength=0;
            Nhits_muoncand=-1;
            Nhits_protoncand=-1;

            TrunMean_cand=0.0;
            TrunMean_pcand=0.0;
            //need to be replaced with a vector for all the protons

            //double vertexXYZ[3]; //only for old use??????????????????????????????????????

            // Initialize a vertex and associated track collection
            std::map< int,std::vector<int> > VertexTrackCollection;
            std::map< int,std::vector<int> > VertexTrackCollection2;

            //int NumTracksNearVertex=0;
            int NumTracksNearVertex2=0;
            //std::cout<<"event number is "<<fEvent<<"total number of vertex is :"<<nuvertexlist[0].size()<<std::endl;

            float nuvtxxtest[500];
            float nuvtxytest[500];
            float nuvtxztest[500];



            //check if there are a track within 5cm from neutrino vertex
            //loop over all the neutrino vertex and check the vertex position
            for(size_t nuvtxIdx=0 ; nuvtxIdx<nuvertexlist[0].size(); nuvtxIdx++){
              NuVertexID.push_back(nuvertexlist[0][nuvtxIdx]->ID());
              Double_t xyz[3]={};
              nuvertexlist[0][nuvtxIdx]->XYZ(xyz);
              NuVertexx.push_back(xyz[0]);
              NuVertexy.push_back(xyz[1]);
              NuVertexz.push_back(xyz[2]);
              nuvtxxtest[nuvtxIdx]=xyz[0];
              nuvtxytest[nuvtxIdx]=xyz[1];
              nuvtxztest[nuvtxIdx]=xyz[2];
    
              TVector3 nuvertexPos(xyz[0],xyz[1],xyz[2]);
  
                //std::cout<<"vertex position is"<< nuvtxIdx <<" "<< xyz[0]<<" "<<xyz[1]<<" "<<xyz[2]<<std::endl;
                //std::cout<<"vertex position is"<< nuvtxIdx <<" "<<nuvtxxtest[nuvtxIdx]<<" "<<nuvtxytest[nuvtxIdx]<<" "<<nuvtxztest[nuvtxIdx]<<std::endl;


                // For each vertex we loop over all tracks looking for matching pairs
                // The outer loop here, then is over one less than all tracks
                unsigned int TrackCountAtVertex2=0;

                unsigned int NTrackflip=0;
                for (size_t trackIdx = 0; trackIdx < trackVecHandle->size(); trackIdx++)
                {
                    // Work with an art Ptr here
                    art::Ptr<recob::Track> track(trackVecHandle,trackIdx);

                    // so we need to get the track direction sorted out.
                    TVector3 trackPos = track->Vertex();
                    TVector3 trackEnd = track->End();

                    // Take the closer end---------------------------------
                    double trackToVertexDist = (trackPos - nuvertexPos).Mag();

                    if ((trackEnd - nuvertexPos).Mag() < trackToVertexDist)
                    {
                        trackPos          = track->End();
                        trackEnd          = track->Vertex();
                        trackToVertexDist = (trackPos - nuvertexPos).Mag();
                        //recalculate the theta and phi here
                        NTrackflip=NTrackflip+1;

                    }

                    //--------------------------------------------------------------------------
                    if (trackToVertexDist<fMinTrk2VtxDist)
                    {
                        NumTracksNearVertex2=NumTracksNearVertex2+1;
                        //if ((trackEnd-trackPos).Mag()>TrackCandLength)  Dec01
                        {
                            // If we are looking at the first track which fulfills the distance to vertex cut
                            
                            //if (!TrackCountAtVertex)  Dec01
                            {
                                // Fill vertex ID into the collection map
                                VertexTrackCollection2.insert(std::pair< int,std::vector<int> >(nuvtxIdx,std::vector<int>()));
                            }

                            // Push back track ID for vertex v
                            VertexTrackCollection2.at(nuvtxIdx).push_back(trackIdx);
                            
                            //VertexTrackCollection2[nuvtxIdx].push_back(trackIdx);

                            // Increase track at vertex count
                            TrackCountAtVertex2++;
                        }
                    } //end of if track distance is within 5cm
                }  //end of loop over the tracks
                //std::cout<<fEvent<<" Ratio of Total number of Tracks flipped divide by the Total number of Tracks "<<float(NTrackflip*1.0)/(trackVecHandle->size()*1.0)<<std::endl;
              
              
            } //end of loop over all the reco neutrino vertex

            //check if there is vertex in FV
            //if(NumTracksNearVertex2>0){
            //std::cout<<"event number 6"<<fEvent<<" total number of vertex is"<<nuvertexlist[0].size()<<endl;
            //}

            //-----------------------------------------------------------
            /*for(size_t vertexIdx = 0; vertexIdx < vertexVecHandle->size(); vertexIdx++)
            {
                // Recover art ptr to vertex
                art::Ptr<recob::Vertex> vertex(vertexVecHandle, vertexIdx);
                
                // Get the position of the vertex
                // Ultimately we really want the vertex position in a TVector3 object...

                unsigned int TrackCountAtVertex=0;
                //---------------------------------------------------------------
                vertex->XYZ(vertexXYZ);
                //fvtxx=vertexXYZ[0]; fvtxy=vertexXYZ[1], fvtxz=vertexXYZ[2];
                
                TVector3 vertexPos(vertexXYZ[0],vertexXYZ[1],vertexXYZ[2]);
                //----------------------------------------------------------------
                // For each vertex we loop over all tracks looking for matching pairs
                // The outer loop here, then is over one less than all tracks
                for (size_t trackIdx = 0; trackIdx < trackVecHandle->size(); trackIdx++)
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
                        //recalculate the theta and phi here


                    }

                    //--------------------------------------------------------------------------
                    if (trackToVertexDist<fMinTrk2VtxDist)
                    {
                        NumTracksNearVertex=NumTracksNearVertex+1;
                        //if ((trackEnd-trackPos).Mag()>TrackCandLength)  Dec01
                        {
                            // If we are looking at the first track which fulfills the distance to vertex cut
                            //if (!TrackCountAtVertex)  Dec01
                            {
                                // Fill vertex ID into the collection map
                                VertexTrackCollection.insert(std::pair< int,std::vector<int> >(vertexIdx,std::vector<int>()));
                            }

                            // Push back track ID for vertex v
                            VertexTrackCollection.at(vertexIdx).push_back(trackIdx);

                            // Increase track at vertex count
                            TrackCountAtVertex++;
                        }
                    } //end of if track distance is within 5cm
                }  //end of loop over the tracks
            } //end of loop over all the vertex
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            */
            // Vertex candidate properties
            //if(NumTracksNearVertex>0){std::cout<<"event number with tracks near vertex"<<std::endl;}

            bool Ntrack_sig=false;
            int NumofRecoTrack=0;
            bool VertexContained=false;
            VertexCandidate = -1;
            float VertexCosTheta = 0.0;
            //check  with Marco about the finding of the primary vertex
            // Loop over the collection of vertices
            if(NumTracksNearVertex2>0){
	      for (auto const& VtxTrack : VertexTrackCollection2)
	    //for (auto const& VtxTrack : VertexTrackCollection)
            {
                // Get vertexID
                int VertexID = VtxTrack.first;

                // Weighted cos theta average
                float WeightedCosTheta = 0.0;

                // Normalization factor
                float NormFactor = 0.0;
                NumofRecoTrack=0;
                //get the vertex position here and check if this is within the FV
                //std::cout<<"event number 8="<<fEvent<<"vertexID="<<VertexID<<" Vertex Position= "<<nuvtxxtest[VertexID]<<" "<<nuvtxytest[VertexID]<<" "<< nuvtxztest[VertexID]<<std::endl;

                if(!inFV(nuvtxxtest[VertexID],nuvtxytest[VertexID], nuvtxztest[VertexID])) continue;
                VertexContained=true;
                //std::cout<<"event number 9="<<fEvent<<"vertexID="<<VertexID<<" Vertex Position= "<<nuvtxxtest[VertexID]<<" "<<nuvtxytest[VertexID]<<" "<< nuvtxztest[VertexID]<<std::endl;
                // Loop over all associated track IDs of this vertex
		for (auto const& TrackID : VtxTrack.second)
                {
                    art::Ptr<recob::Track> track(trackVecHandle,TrackID);
		    
		    // Calculate track range (trackstart - trackend)
		    //double TrackRange = GetTrackRange(track);
                    double TrackRange= GetTrackLength(track);
                    if(TrackRange==0){
                    std::cout<<"Event number is "<<fEvent<<" Track Length is 0  "<<TrackID<<std::endl;
                    }
                    // Add all track range of tracks close to vertex (Normalization factor of the weighted average)
                    NormFactor += TrackRange;
                    // Add cos(theta) weighted by track range
                    WeightedCosTheta += TrackRange*cos(track->Theta());
                    NumofRecoTrack++;
                }// track ID loop

                // Make average
                WeightedCosTheta /= NormFactor;
                //std::cout<<"Event number is "<<fEvent<<" NormFactor is "<<NormFactor<<std::endl;
                // Check for flatest angle (also backwards pointing)
                if (fabs(WeightedCosTheta) > VertexCosTheta)
                {
                    VertexCandidate = VertexID;
                    VertexCosTheta = fabs(WeightedCosTheta);
                }
                if(NumofRecoTrack>=2){
                  Ntrack_sig=true;
                } 
            }// vertex collection loop
            //std::cout<<"event number 11"<<fEvent<<"VertexCandidate= "<<VertexCandidate<<std::endl;
            }
       
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            /* if(NumTracksNearVertex2>0 ){
                  if( VertexContained==true)
                  {
                     fMC_vtxinFV->Fill();                       
                     if(Ntrack_sig==true)
                     {
                       fMC_ntrks->Fill();
                     }
                  }
            }
            */ 
            //std::cout<<"event number 7= "<<fEvent<<"VertexCandidate= "<<VertexCandidate<<std::endl; 
 
            /*{// Create the candidate vertex begin of shower check
            art::Ptr<recob::Vertex> vertex(vertexVecHandle,VertexCandidate);
	    
	    // Fill vertex candidate coordinates
	    vertex->XYZ(vertexXYZ);
            fvtxx=vertexXYZ[0]; fvtxy=vertexXYZ[1], fvtxz=vertexXYZ[2];
            //std::cout<<"event number= "<<fEvent<<"vertex in FV"<<fvtxx<<" "<<fvtxy<<"  "<<fvtxz<<std::endl;
            */
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            //get the vertex candidate position first


            //get the shower location and make sure there is no shower close to the primary vertex
            //distance from primary vertex to start point of shower or the perpendicular distance
            bool noshowerFlag=true;
            /* check total number of showers for each event
            int no_showers=0;
            //loop over all the showers
            //std::vector<recob::Shower> const & showers;
  	    //for(std::vector<sim::MCShower>::const_iterator imcshwr = showerHandle->begin();    
	    //imcshwr != showerHandle->end(); ++imcshwr) {
            if(showerHandle.isValid()){
                   no_showers=showerHandle->size();
    	    }     
            //}
            */
            //for(recob::Shower const& shwr: *showerHandle) {  //another way to loop over all the showers
            for(size_t shwrIdx=0; shwrIdx<  showerHandle->size(); shwrIdx++){
                  if(showerHandle.isValid()&& showerHandle->size()>0){
                  //handle allows to access to the showerr at ***
                  //sim::MCShower const & shwr = showerHandle->at(shwrIdx);
                  recob::Shower const & shwr = showerHandle->at(shwrIdx);
 
                  //get the start and end point position of the shower
                  //TVector3 const& shwrPos = {shwr.Start().X(), shwr.Start().Y(), shwr.Start().Z()};
                  TVector3 const & shwrPos = shwr.ShowerStart();
                  //TVector3 const& shwrEnd = shwr.ShowerEnd();
                  //TVector3 const& shwrDir = shwr.Direction ();                  
                  TVector3 const& temp={nuvtxxtest[VertexCandidate], nuvtxytest[VertexCandidate],nuvtxztest[VertexCandidate]};
                  //TVector3 const& temp = {vertexXYZ[0],vertexXYZ[1],vertexXYZ[2]};
                  //vershwrdist=(shwrPos-temp).Angle(shwrDir);
                  vershwrdist=(shwrPos-temp).Mag();
                  
               
                  if(vershwrdist<=50)  noshowerFlag=false;
                  }
	    }
            fvtxx=nuvtxxtest[VertexCandidate];
            fvtxy=nuvtxytest[VertexCandidate];
            fvtxz=nuvtxztest[VertexCandidate];
    
            //std::cout<<"check if there is a shower near the primary vertex"<<noshowerFlag<<std::endl;
            // Check if the wertex candidate is contained and pick the longest track
            // check with Marco and others with this
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
            if(NumTracksNearVertex2>0 ){
                  if( VertexContained==true)
                  {
                     fMC_vtxinFV->Fill();                       
                     if(Ntrack_sig==true)
                     {
                       fMC_ntrks->Fill();
                       if(noshowerFlag==true){
                          fMC_noshwr->Fill();
                       } 
                     }
                  }
            }
             
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            //call the track momentum algorithm that gives you momentum based on track range
            trkf::TrackMomentumCalculator trkm;
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            float momentum=0.0;
            if(NumTracksNearVertex2>0 ){
            if( VertexContained==true)
            {
            if(Ntrack_sig==true)
            {
            if(noshowerFlag==true)
            {   //std::cout<<"event number 12 "<<fEvent<<"the vertex candidate is: "<<VertexCandidate<<std::endl;
                // Looping over track IDs of tracks associated with the vertex candidate
		for (auto const& TrackID : VertexTrackCollection2.find(VertexCandidate)->second)
                {
                    art::Ptr<recob::Track> track(trackVecHandle,TrackID);
                    //recob::MCSFitResult const & mcsfitresult = MCSFitHandle.at(TrackID);
		    momentum = mcsfitlist[TrackID]->bestMomentum();
            bool mcs_isBestFwd = mcsfitlist[TrackID]->isBestFwd();
            if (!mcs_isBestFwd){std::cout << "muon candidate MCS fits better backwards" << std::endl;}
                    /*
                    * mcsfitlist[TrackID]->
                    *
                    */
                    //if (!mcsfitlist[TrackID]->isBestFwd()) continue;
                    //{
                      //We should ignore this because it's probably a cosmic
                      //std::cout<<"momentum"<<momentum<<std::endl;
                    //}

		    // Calculate track range (trackstart - trackend)
		    double TrackRange = GetTrackRange(track);
                    double TrackLength= GetTrackLength(track);
                    //check again if need to swap the start and end point position here

                    TVector3 trackPos = track->Vertex();
                    TVector3 trackEnd = track->End();
                    float trackMom=track->VertexMomentum();
                    float trackTheta=track->Theta();
                    float trackPhi=track->Phi();         		


                    //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^    
                    // Check for if track is longer
                    if( TrackLength > TrackCandLength && (inFV(trackPos.x(), trackPos.y(),trackPos.z()) || inFV(trackEnd.x(),trackEnd.y(),trackEnd.z())))
                    {
			// Pick the numbers of the longest track
                        TrackCandidate = TrackID;
                        //TrackCandLength = TrackRange;
                        TrackCandLength = TrackLength;
                        trackmomcandidate=trackMom;
                        trackmomcandidate_mcs=momentum;    
                        trackstartzcandidate=trackPos.z();
                        trackstartxcandidate=trackPos.x();
                        trackstartycandidate=trackPos.y();
                        trackendzcandidate=trackEnd.z();
                        trackendxcandidate=trackEnd.x();
                        trackendycandidate=trackEnd.y();
                        fThetaLep=trackTheta;
                        fPhiLep=trackPhi;
                        fLlep=TrackRange;
                        /*
                        */
                    }
                } //end of loop over all the tracks associated to the vertex candidate
            }
            }
            }
            }  //end of if NumofRecoTracks>0, vertex contained and Ntrack_sig=true
            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            //Set the value for fPlep for the lepton candidate
            //if the muon is within the FV
            if(inFV(trackstartxcandidate, trackstartycandidate, trackstartzcandidate) &&inFV(trackendxcandidate, trackendycandidate, trackendzcandidate)){
              fPlep=trkm.GetTrackMomentum(TrackCandLength,13);
            }           
            //if it is an existing muon
            else {fPlep=trackmomcandidate_mcs;}
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<, 
             
            bool trackFlashFlag=false;
            if(NumTracksNearVertex2>0 ){
            if( VertexContained==true)
            {
            if(Ntrack_sig==true)
            {
            if(noshowerFlag==true)
            {   
            if(TrackCandidate>-1 && VertexCandidate>-1){
               // Create the candidate track
               art::Ptr<recob::Track>  track(trackVecHandle,TrackCandidate);
                                           //(handle, index/key)
               art::FindManyP<recob::Hit> fmh(trackVecHandle, event, fTrackModuleLabel);
               std::vector<art::Ptr<recob::Hit>> hits = fmh.at(track.key());
               Nhits_muoncand=hits.size(); 
               //std::cout<<"total number of hits of muon candidate is : "<<hits.size()<<std::endl;

               art::FindManyP<anab::Calorimetry> calos_from_track(trackVecHandle, event, fCalorimetryModuleLabel);
 
               std::vector<art::Ptr<anab::Calorimetry>> calos = calos_from_track.at(track.key());
               //std::cout<<"<<<<<<<<<<<<<<<<<,libotest before calculating truncated dQdx of muon2 "<<std::endl;
                            
               TrunMean_cand=GetDqDxTruncatedMean(calos);     
               //loop over all the hits of muon and get dedx and residual range for each hit 
               


               ftrklenmuoncand=TrackCandLength;
               //use the backtracker to get the true PDG of a track/particle

               //------------------------------------------------------------------------------------------
               // Check if track and flash are matched
               flstrkdist=GetFlashTrackDist(flashPtr->ZCenter(), track->Vertex().z(), track->End().z());
               if(GetFlashTrackDist(flashPtr->ZCenter(), track->Vertex().z(), track->End().z()) < fFlashWidth) {
                trackFlashFlag=true; 
               }
               //Get the truth information of the muon candidate including true PDG and Momentum
               if(isMC){
               double tmpEfracm=0;
               double tmpEcompletm=0;
               const simb::MCParticle *mparticle;
               truthMatcher(hitlist, hits, mparticle, tmpEfracm, tmpEcompletm );

               const art::Ptr<simb::MCTruth> MCtruth = fMCTruthMatching->ParticleToMCTruth(mparticle);
               //std:string pri("primary");
               trackcand_origin=MCtruth->Origin();
               trackcand_nuset=MCtruth->NeutrinoSet(); 
 

               if(mparticle){
                  trackcand_parPDG=mparticle->PdgCode();
                  trackcand_parStatusCode=mparticle->StatusCode(); 
                  trackcand_parTheta=mparticle->Momentum().Theta(); 
                  trackcand_parCosTheta=TMath::Cos(mparticle->Momentum().Theta());
                  trackcand_parSinTheta=TMath::Sin(mparticle->Momentum().Theta());
                  trackcand_parE=mparticle->E();
                  trackcand_parMass=mparticle->Mass();
                  trackcand_parKE=(mparticle->E())-(mparticle->Mass());
                  trackcand_parEndE=mparticle->EndE();
                  trackcand_parPx=mparticle->Px();
                  trackcand_parPy=mparticle->Py();
                  trackcand_parPz=mparticle->Pz();
                  trackcand_parPhi=mparticle->Momentum().Phi();
                  trackcand_parCosPhi=TMath::Cos(mparticle->Momentum().Phi());
                  trackcand_parSinPhi=TMath::Cos(mparticle->Momentum().Phi());
               }
               }
               //----------------------------------------------------------------------------------------- 
            }
            //std::cout<<"Event number= "<<fEvent<<" true PDG id of the track candidate is "<<trackcand_parPDG<<std::endl;
            if (trackFlashFlag==true){
              fMC_trkfls->Fill();
            }//end of if trackFlagFlag is true
            }//end of if noshowerFlag
            }//end of Ntrack_sig=true
            }//end of if vertex contained
            }//end of if the number rof tracks near vertex2 is greater than 0
            //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>.
            TrackProtonCandidate=-1;
            TrackProtonCandLength=0.0;	
            //initialize all the vectors for all the proton candidates


            if(NumTracksNearVertex2>0 ){
            if( VertexContained==true)
            {
            if(Ntrack_sig==true)
            {
            if(noshowerFlag==true)
            {   
            if(trackFlashFlag==true)
            {
            //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
                fNRecoPTrks=0;
                fNRecoTrks=0;
                int ntrkstovtx=0;
                trackidpcand->clear();
                trackstartxpcand->clear();
                trackstartypcand->clear();
                trackstartzpcand->clear();
                trackendxpcand->clear();
                trackendypcand->clear();
                trackendzpcand->clear();
                trackmompcand->clear();
                trackthetapcand->clear();
                tracklengthpcand->clear();
                trackphipcand->clear();
                tracktrunmeanpcand->clear();
                trackpidapcand->clear();
		for (auto const& TrackID : VertexTrackCollection2.find(VertexCandidate)->second)
                {
                    art::Ptr<recob::Track> track(trackVecHandle,TrackID);

                    art::FindMany<anab::ParticleID> fmpid(trackVecHandle, event, "pandoraNuPID");
                    double PIDA(-1.);
                    if(fmpid.isValid()) {
                      std::vector<const anab::ParticleID*> pids = fmpid.at(TrackID);
                      for (size_t ipid = 0; ipid < pids.size(); ++ipid){
                        if (!pids[ipid]->PlaneID().isValid) continue;
                        int planenum = pids[ipid]->PlaneID().Plane;
                        if (planenum != 2){continue;}
                        PIDA = pids[ipid]->PIDA();
                      }
                    }

		    
		    // Calculate track range (trackstart - trackend)
		    double TrackRange = GetTrackRange(track);
                    double TrackLength= GetTrackLength(track);

                    TVector3 trackPos = track->Vertex();
                    TVector3 trackEnd = track->End();
                    float trackMom=track->VertexMomentum();
                    float trackTheta=track->Theta();
                    float trackPhi=track->Phi();         		    
                    // Check for if track is longer
                    //make sure the proton candidate is also within 5cm from the vertex
                    //check the number rof track here       
                    ntrkstovtx=ntrkstovtx+1;
                    if(TrackID ==TrackCandidate)   continue;           
                    if(inFV(trackPos.X(), trackPos.Y(), trackPos.Z())&& inFV(trackEnd.X(), trackEnd.Y(), trackEnd.Z())){
                      trackidpcand->push_back(TrackID);
                      trackstartxpcand->push_back(trackPos.X());
                      trackstartypcand->push_back(trackPos.Y());    
                      trackstartzpcand->push_back(trackPos.Z());
                      trackendxpcand->push_back(trackEnd.X());
                      trackendypcand->push_back(trackEnd.Y());    
                      trackendzpcand->push_back(trackEnd.Z());

                      //trackmompcand.push_back(track->VertexMomentum());
                      trackmompcand->push_back(trkm.GetTrackMomentum(TrackLength, 2212));
                      trackthetapcand->push_back(trackTheta);
                      tracklengthpcand->push_back(TrackLength);
                      trackphipcand->push_back(trackPhi);
                      trackpidapcand->push_back(PIDA);
                      //get the truncated mean dqdx of all the proton candidates-------------
                      art::FindManyP<recob::Hit> fmh(trackVecHandle, event, fTrackModuleLabel);
                      std::vector<art::Ptr<recob::Hit>> hits =fmh.at(track.key());
                      
                      art::FindManyP<anab::Calorimetry> calos_from_track(trackVecHandle, event, fCalorimetryModuleLabel) ;        
                      std:: vector<art::Ptr<anab::Calorimetry>> calos=calos_from_track.at(track.key());
                      tracktrunmeanpcand->push_back(GetDqDxTruncatedMean(calos));
                      //---------------------------------------------------------------------- 
                    }





                    //if (TrackRange > TrackProtonCandLength && TrackID !=TrackCandidate)
                    if (TrackLength > TrackProtonCandLength && inFV(trackPos.X(), trackPos.Y(), trackPos.Z())&& inFV(trackEnd.X(), trackEnd.Y(), trackEnd.Z()))
                    {
			// Pick the numbers of the longest track
                        TrackProtonCandidate = TrackID;
                        //TrackProtonCandLength = TrackRange;
                        TrackProtonCandLength =TrackLength;
                        trackmomprotoncandidate=trackMom;
                        trackstartzprotoncandidate=trackPos.z();
                        trackstartxprotoncandidate=trackPos.x();
                        trackstartyprotoncandidate=trackPos.y();
                        trackendzprotoncandidate=trackEnd.z();
                        trackendxprotoncandidate=trackEnd.x();
                        trackendyprotoncandidate=trackEnd.y();
                        fThetaHad=trackTheta;
                        fPhiHad=trackPhi;
                        fLhad=TrackRange;
                        fPhad=trkm.GetTrackMomentum(TrackProtonCandLength,2212);;
  
                    }
                } //end of loop over all the tracks associated to track candidate
                 
                //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<,
                ftrklenprotoncand=TrackProtonCandLength;
                                
                //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^
                bool ProtonTag=true;
                fNRecoTrks=trackidpcand->size()+1;
                fNRecoPTrks=trackidpcand->size();
                Evis=0.0;
                float Eptot=0.0;
                float Pxptot=0.0;
                float Pyptot=0.0;
                float Pzptot=0.0;  
                for(int pcand=0; pcand<fNRecoPTrks; pcand++){
                       //std::cout<<"event number= "<<fEvent<<" total number of protons are "<<trackidpcand.size()<<" "<<tracklengthpcand[pcand]<<" "<<tracktrunmeanpcand[pcand]<<std::endl;
                       /*if((tracktrunmeanpcand[pcand]<300 &&tracklengthpcand[pcand]<(-16*tracktrunmeanpcand[pcand]+5000)) ||
                          (tracktrunmeanpcand[pcand]>300 &&tracklengthpcand[pcand]<(20000*exp(-0.015*tracktrunmeanpcand[pcand])))) {
                            ProtonTag=false;              
                       }
                       */
                       
                       //std::cout<<"pcand= "<<pcand<<" dqdx of pcand is "<<tracktrunmeanpcand[pcand]<<" trklen of pcand is "<<tracklengthpcand[pcand]<<std::endl;
                       if(MIPConsistency(tracktrunmeanpcand->at(pcand), tracklengthpcand->at(pcand))){
                             ProtonTag=false;
                       }
                       //Evis=fPlep+trackmompcand[pcand];
                       //Eptot=Eptot+sqrt(trackmompcand[pcand]*trackmompcand[pcand]+protonmass*protonmass);
                       //Pxptot=Pxptot+trackmompcand[pcand]*TMath::Sin(trackthetapcand[pcand])*TMath::Cos(trackphipcand[pcand]);
                       //Pyptot=Pyptot+trackmompcand[pcand]*TMath::Sin(trackthetapcand[pcand])*TMath::Sin(trackphipcand[pcand]);
                       //Pzptot=Pzptot+trackmompcand[pcand]*TMath::Cos(trackthetapcand[pcand]);


                       Evis=fPlep+trackmompcand->at(pcand);
                       Eptot=Eptot+TMath::Sqrt(trackmompcand->at(pcand)*trackmompcand->at(pcand)+protonmass*protonmass);
                       Pxptot=Pxptot+trackmompcand->at(pcand)*TMath::Sin(trackthetapcand->at(pcand))*TMath::Cos(trackphipcand->at(pcand));
                       Pyptot=Pyptot+trackmompcand->at(pcand)*TMath::Sin(trackthetapcand->at(pcand))*TMath::Sin(trackphipcand->at(pcand));
                       Pzptot=Pzptot+trackmompcand->at(pcand)*TMath::Cos(trackthetapcand->at(pcand));

                }
                

                //calculate Q2 and W from here
                 
                float Emuoncand=TMath::Sqrt(muonmass*muonmass+fPlep*fPlep);
                Q2cal=-(Evis-Emuoncand)*(Evis-Emuoncand)+fPlep*fPlep; 
                Wcal=TMath::Sqrt(Eptot*Eptot-Pxptot*Pxptot-Pyptot*Pyptot-Pzptot*Pzptot);                
                    

                if(ntrkstovtx==fNRecoTrks){
                   fMC_NoExTrk->Fill();
              
                if(TrackProtonCandidate>-1 && TrackCandidate>-1){
                     //trackidmuoncand=TrackCandidate;
                     //trackidprotoncand=TrackProtonCandidate;
                     fMC_mupinFV->Fill();
                     //get the truncated dQdx of proton candidate
                     art::Ptr<recob::Track>  track(trackVecHandle,TrackProtonCandidate);

                     art::FindManyP<recob::Hit> fmh(trackVecHandle, event, fTrackModuleLabel);
                     std::vector<art::Ptr<recob::Hit> > hits = fmh.at(track.key());
                     Nhits_protoncand=hits.size();
                     //std::cout<<"event number= 17  "<<fEvent<<" VertexCandidate is "<<VertexCandidate<<std::endl;
 
                     art::FindManyP<anab::Calorimetry> calos_from_track(trackVecHandle, event, fCalorimetryModuleLabel);
 
                     std::vector<art::Ptr<anab::Calorimetry>> calos = calos_from_track.at(track.key());
                        
                     TrunMean_pcand=GetDqDxTruncatedMean(calos);    
                     //loop over all the hits of proton candidate and get the dedx and residual range here
                  
 



                     //perform the cut to TrunMean_cand and TrunMean_pcand and Fill the tree 
                     //std::cout<<"event number= "<<fEvent<<"Track Candidate is "<<TrackCandidate<<" Track Proton Candidate is "<<TrackProtonCandidate<<" TrunMean cand = "<<TrunMean_cand<<" TrunMean_pcand = "<<TrunMean_pcand<<std::endl; 
                     //get the MC truth of the proton candidate if MC events 
                     //---------------------------------------------------------------------
                     if(isMC){
                     double tmpEfracm=0;
                     double tmpEcompletm=0;
                     const simb::MCParticle *mparticle;
                     truthMatcher(hitlist, hits, mparticle, tmpEfracm, tmpEcompletm );

                     //check if this particle is from neutrino interaction
                     const art::Ptr<simb::MCTruth> MCtruth = fMCTruthMatching->ParticleToMCTruth(mparticle);
                     trackpcand_origin=MCtruth->Origin();
                     trackpcand_nuset=MCtruth->NeutrinoSet(); 
 
                     if(mparticle){
                         trackpcand_parPDG=mparticle->PdgCode();
                         trackpcand_parStatusCode=mparticle->StatusCode(); 
                         trackpcand_parTheta=mparticle->Momentum().Theta(); 
                         trackpcand_parCosTheta=TMath::Cos(mparticle->Momentum().Theta());
                         trackpcand_parSinTheta=TMath::Sin(mparticle->Momentum().Theta());
                         trackpcand_parE=mparticle->E();
                         trackpcand_parMass=mparticle->Mass();
                         trackpcand_parKE=(mparticle->E())-(mparticle->Mass());
                         trackpcand_parEndE=mparticle->EndE();
                         trackpcand_parPx=mparticle->Px();
                         trackpcand_parPy=mparticle->Py();
                         trackpcand_parPz=mparticle->Pz();
                         trackpcand_parPhi=mparticle->Momentum().Phi();
                         trackpcand_parCosPhi=TMath::Cos(mparticle->Momentum().Phi());
                         trackpcand_parSinPhi=TMath::Cos(mparticle->Momentum().Phi());
                     }
                     } //end of ifMC          
                      //std::cout<<"Event number= "<<fEvent<<" true PDG id of the track proton candidate is "<<trackpcand_parPDG<<std::endl;



                     //-----------------------------------------------------------------------
                     if(TrunMean_cand !=-9999 && TrunMean_pcand !=-9999){
                       /* 
                       if((TrunMean_cand<300 &&ftrklenmuoncand<(-16*TrunMean_cand+5000)) ||
                          (TrunMean_cand>300 &&ftrklenmuoncand<(20000*exp(-0.015*TrunMean_cand)))) {
                          if(ProtonTag==true){
                                fMC_TrunMean->Fill();
                          }  
                       }*/
                       //std::cout<<"is this a muon>>>"<<MIPConsistency(TrunMean_cand, ftrklenmuoncand)<<std::endl;
                       
                       if(MIPConsistency(TrunMean_cand, ftrklenmuoncand)){
                           //std::cout<<"ProtonTag is "<<ProtonTag<<std::endl;
                           if(ProtonTag==true){
                            fMC_TrunMean->Fill();
                           }
                       } 

                     } //end of if trun mean not equal to -9999 
                     //-------------------------------------------------------------------------
                } //end of if TrackCandidate and TrackProtonCandidate 
                }  //end of if there is no extra tracks around vertex candidate    
                //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
            //~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~^~
            } //end of trkfls cut
            } //end of no shwr cut 
            } //end of require Ntrk_sig=2
            } //end of vertex contained
            } //end of NumTrksVertex2>0
            //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
       
        } //end of if the vertex is valid
     //-------------------------------------------------------
    } //end of if flashtag is true
    } //end of if there is only one neutrino event
    ClearLocalData();
   //-----------------------------------------------------------------------------------------------
    return;
}

DEFINE_ART_MODULE( CC1uNPSelAna)

} // namespace  CC1uNPSelAna

#endif //  CC1uNPSelAna_Module
