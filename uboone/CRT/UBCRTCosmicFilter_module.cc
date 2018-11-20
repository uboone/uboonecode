////////////////////////////////////////////////////////////////////////
// Class:       UBCRTCosmicFilter
// Plugin Type: filter (art v2_11_03)
// File:        UBCRTCosmicFilter_module.cc
//
// Generated at Fri Sep 28 13:12:47 2018 by Christopher Barnes using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// declare the package with the CRT hit information.
#include "uboone/CRT/CRTProducts/CRTHit.hh"
#include "lardataobj/RecoBase/OpFlash.h"
#include "uboone/RawData/utils/DAQHeaderTimeUBooNE.h"

// Declare the file for the associations to be made.
#include "lardata/Utilities/AssociationUtil.h"

// C++
#include <memory>

// ROOT
#include <TTree.h>

class UBCRTCosmicFilter;

class UBCRTCosmicFilter : public art::EDFilter
{
public:
  explicit UBCRTCosmicFilter(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UBCRTCosmicFilter(UBCRTCosmicFilter const &) = delete;
  UBCRTCosmicFilter(UBCRTCosmicFilter &&) = delete;
  UBCRTCosmicFilter &operator=(UBCRTCosmicFilter const &) = delete;
  UBCRTCosmicFilter &operator=(UBCRTCosmicFilter &&) = delete;

  // Use the 'beginJob()' function to save information in the TTree.
  void beginJob() override;

  // Required functions.
  bool filter(art::Event &e) override;

private:
  // Declare member data here.
  std::string fBeamFlashProducer;
  std::string fCRTHitProducer;
  std::string fDAQHeaderProducer;
  double fBeamStart;
  double fBeamEnd;
  double fPEMin;
  double fDTOffset;
  double fResolution;
  bool fuseAsFilter;
  bool verbose;
  bool fAnode, fTop, fBottom, fCathode;

  art::ServiceHandle<art::TFileService> tfs;

  // TTree Declaration.
  TTree *_tree;

  // Event info.
  int _run;
  int _subrun;
  int _event;

  // Beam flash info.
  int _nflashes_in_beamgate;
  int _nflashes_in_beamgate_passing_beamspill_and_PE_cuts;
  double _beam_flash_time;
  double _beam_flash_PE;

  // CRT hit info.
  int _nCRThits_in_event;
  // Fields storing information about the CRT hit closest in time to the Flash
  double _CRT_hit_time;
  double _CRT_hit_PE;
  double _CRT_hit_x;
  double _CRT_hit_y;
  double _CRT_hit_z;
  // Fields storing information about all CRT hits in event
  std::vector<float> _CRT_hits_time;
  std::vector<float> _CRT_hits_PE;
  std::vector<float> _CRT_hits_x;
  std::vector<float> _CRT_hits_y;
  std::vector<float> _CRT_hits_z;

  // Comparison between CRT hit time & CRT hit PE.
  double _dt;

  // Is the closest CRT hit within resolution?  Binary flag.
  int _within_resolution;

  // Counter for the events.
  int event_counter;
};

UBCRTCosmicFilter::UBCRTCosmicFilter(fhicl::ParameterSet const &p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  produces<art::Assns<crt::CRTHit, recob::OpFlash>>();

  fBeamFlashProducer = p.get<std::string>("BeamFlashProducer");
  fCRTHitProducer = p.get<std::string>("CRTHitProducer");
  fDAQHeaderProducer = p.get<std::string>("DAQHeaderProducer");
  fBeamStart = p.get<double>("BeamStart");
  fBeamEnd = p.get<double>("BeamEnd");
  fPEMin = p.get<double>("PEMin");
  fDTOffset = p.get<double>("DTOffset");
  fResolution = p.get<double>("Resolution");
  fuseAsFilter = p.get<bool>("useAsFilter");
  verbose = p.get<bool>("verbose");
  fTop     = p.get<bool>("Top");
  fBottom  = p.get<bool>("Bottom");
  fAnode   = p.get<bool>("Anode");
  fCathode = p.get<bool>("Cathode");
}

void UBCRTCosmicFilter::beginJob()
{

  art::ServiceHandle<art::TFileService> tfs;
  _tree = tfs->make<TTree>("_tree", "Flash/CRT Matching Information For All Events");
  _tree->Branch("_run", &_run, "run/I");
  _tree->Branch("_subrun", &_subrun, "subrun/I");
  _tree->Branch("_event", &_event, "event/I");
  _tree->Branch("_nflashes_in_beamgate", &_nflashes_in_beamgate, "nflashes_in_beamgate/I");
  _tree->Branch("_nflashes_in_beamgate_passing_beamspill_and_PE_cuts", &_nflashes_in_beamgate_passing_beamspill_and_PE_cuts, "nflashes_in_beamgate_passing_beamspill_and_PE_cuts/I");
  _tree->Branch("_beam_flash_time", &_beam_flash_time, "beam_flash_time/D");
  _tree->Branch("_beam_flash_PE", &_beam_flash_PE, "beam_flash_PE/D");
  _tree->Branch("_nCRThits_in_event", &_nCRThits_in_event, "_nCRThits_in_event/I");
  _tree->Branch("_CRT_hit_time", &_CRT_hit_time, "CRT_hit_time/D");
  _tree->Branch("_CRT_hit_PE", &_CRT_hit_PE, "CRT_hit_PE/D");
  _tree->Branch("_CRT_hit_x", &_CRT_hit_x, "CRT_hit_x/D");
  _tree->Branch("_CRT_hit_y", &_CRT_hit_y, "CRT_hit_y/D");
  _tree->Branch("_CRT_hit_z", &_CRT_hit_z, "CRT_hit_z/D");
  _tree->Branch("CRT_hits_time", "std::vector< float >", &_CRT_hits_time);
  _tree->Branch("CRT_hits_PE", "std::vector< float >", &_CRT_hits_PE);
  _tree->Branch("CRT_hits_x", "std::vector< float >", &_CRT_hits_x);
  _tree->Branch("CRT_hits_y", "std::vector< float >", &_CRT_hits_y);
  _tree->Branch("CRT_hits_z", "std::vector< float >", &_CRT_hits_z);
  _tree->Branch("_dt", &_dt, "dt/D");
  _tree->Branch("_within_resolution", &_within_resolution, "within_resolution/I");
}

bool UBCRTCosmicFilter::filter(art::Event &e)
{
  // Implementation of required member function here.

  // Fill the event info.
  _event = e.event();
  _subrun = e.subRun();
  _run = e.run();

  if (verbose)
    std::cout << "Now looping over event #" << event_counter << "." << std::endl;

  // Iterate the event counter.
  event_counter++;

  // Save an association for CRT hits/PMT flashes.
  std::unique_ptr<art::Assns<crt::CRTHit, recob::OpFlash>> crthit_flash_assn_v(new art::Assns<crt::CRTHit, recob::OpFlash>);

  // Declare an object for the GPS timestamp of the event so that you can offset the cosmic t0 times.
  art::Handle<raw::DAQHeaderTimeUBooNE> rawHandle_DAQHeader;
  e.getByLabel(fDAQHeaderProducer, rawHandle_DAQHeader);
  if(!rawHandle_DAQHeader.isValid()) {
    e.put(std::move(crthit_flash_assn_v));
    return !fuseAsFilter;
  }

  raw::DAQHeaderTimeUBooNE const &my_DAQHeader(*rawHandle_DAQHeader);
  art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
  double evt_timeGPS_nsec = evtTimeGPS.timeLow();

  // load CRT hits.
  art::Handle<std::vector<crt::CRTHit>> crthit_h;
  e.getByLabel(fCRTHitProducer, crthit_h);

  // make sure CRT hits look good.
  if (!crthit_h.isValid()) {
    e.put(std::move(crthit_flash_assn_v));
    return !fuseAsFilter;
  }

  // load beam flashes here.
  art::Handle<std::vector<recob::OpFlash>> beamflash_h;
  e.getByLabel(fBeamFlashProducer, beamflash_h);

  // make sure beam flashes look good.
  if (!beamflash_h.isValid())
  {
    std::cerr << "\033[93m[ERROR]\033[00m ... could not locate beam flash!" << std::endl;
    throw std::exception();
  }

  // Set the variable for the number of flashes reconstructed in the beamgate for each event.
  _nflashes_in_beamgate = beamflash_h->size();
  _nflashes_in_beamgate_passing_beamspill_and_PE_cuts = 0;

  // Loop through the beam flashes and find the closest crt hit to the beam flashes.
  if (verbose)
    std::cout << "Number of beam flashes in this event = " << _nflashes_in_beamgate << "." << std::endl;

  // Set the variables for the closest CRT hit time.
  double _dt_abs = 100000.0;
  size_t flash_idx = 0;
  _dt = 0.;
  _CRT_hit_time = 0.;

  // Set the variable for if the CRT hit is within resolution of the flash.
  _within_resolution = 0;

  // Set the value for the '_beam_flash_PE' to be negative so that it will be overwritten by the information for the actual flash.
  _beam_flash_time = -10000.0;
  _beam_flash_PE = -1.0;

  // Loop over the flashes to find the one greatest in intensity.
  for (int i = 0; i < _nflashes_in_beamgate; i++)
  {

    // Throw the flash out if it is outside of the beamgate or too low in intensity.
    if (beamflash_h->at(i).Time() < fBeamStart || beamflash_h->at(i).Time() > fBeamEnd || beamflash_h->at(i).TotalPE() < fPEMin)
    {

      if (verbose)
        std::cout << "Flash is outside of beam gate or too low in intensity.  Skipping!" << std::endl;

      continue;
    }

    // Increment '_nflashes_in_beamgate_passing_beamspill_and_PE_cuts'.
    _nflashes_in_beamgate_passing_beamspill_and_PE_cuts++;

    if (beamflash_h->at(i).TotalPE() > _beam_flash_PE)
    {
      _beam_flash_time = beamflash_h->at(i).Time();
      _beam_flash_PE = beamflash_h->at(i).TotalPE();
      flash_idx = size_t(i);

    }

  } // End of the loop over the beam flashes.

  if (verbose)
  {
    std::cout << "Number of flashes in the event passing the beamspill and PE cuts = " << _nflashes_in_beamgate_passing_beamspill_and_PE_cuts << "." << std::endl;
    std::cout << "Time of greatest flash = " << _beam_flash_time << " us." << std::endl;
    std::cout << "PEs of greatest flash = " << _beam_flash_PE << " PEs." << std::endl;
  }

  // Set the variable for the number of CRT hits in the event.
  _nCRThits_in_event = crthit_h->size();
  // Reset the vectors to store the CRThit information
  _CRT_hits_time.clear();
  _CRT_hits_PE.clear();
  _CRT_hits_x.clear();
  _CRT_hits_y.clear();
  _CRT_hits_z.clear();

  // Loop over the CRT hits.
  for (int j = 0; j < _nCRThits_in_event; j++)
  {

    // figure out what plane this hit comes from
    // 3 -> top
    // 0 -> bottom
    // 1 -> anode
    // 2 -> cathode
    auto pl = crthit_h->at(j).plane;

    if ( (pl == 3) && !fTop)
      continue;

    if ( (pl == 2) && !fCathode)
      continue;

    if ( (pl == 1) && !fAnode)
      continue;

    if ( (pl == 0) && !fBottom)
      continue;

    //if (verbose)
    //std::cout << "Time of the CRT Hit wrt the event timestamp = " << ((crthit_h->at(j).ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.) << " us." << std::endl;

    double _crt_time_temp = ((crthit_h->at(j).ts0_ns - evt_timeGPS_nsec + fDTOffset) / 1000.);
    // Fill the vector variables.
    _CRT_hits_time.push_back(_crt_time_temp);
    _CRT_hits_PE.push_back(crthit_h->at(j).peshit);
    _CRT_hits_x.push_back(crthit_h->at(j).x_pos);
    _CRT_hits_y.push_back(crthit_h->at(j).y_pos);
    _CRT_hits_z.push_back(crthit_h->at(j).z_pos);
    
    if (fabs(_beam_flash_time - _crt_time_temp) < _dt_abs)
    {
      _dt_abs = fabs(_beam_flash_time - _crt_time_temp);
      _dt = _beam_flash_time - _crt_time_temp;
      _CRT_hit_time = _crt_time_temp;

      // set 'within_resolution' to 'true' and break the loop if 'closest_crt_diff' is less than fResolution.
      if (_dt_abs < fResolution)
      {
        _within_resolution = 1;

        // Set the position information and the intensity of the CRT hit.
        _CRT_hit_PE = crthit_h->at(j).peshit;
        _CRT_hit_x = crthit_h->at(j).x_pos;
        _CRT_hit_y = crthit_h->at(j).y_pos;
        _CRT_hit_z = crthit_h->at(j).z_pos;

        // Convert the CRT iterator to type 'size_t'.
        size_t crt_idx = size_t(j);

        // Make the two pointers.
        const art::Ptr<crt::CRTHit> crt_ptr(crthit_h, crt_idx);
        const art::Ptr<recob::OpFlash> flash_ptr(beamflash_h, flash_idx);

        // Make the association between the CRT hit and the OpFlash.
        crthit_flash_assn_v->addSingle(crt_ptr, flash_ptr);

        if (verbose)
        {
          std::cout << "CRT hit PE = " << _CRT_hit_PE << " PEs." << std::endl;
          std::cout << "CRT hit x = " << _CRT_hit_x << " cm." << std::endl;
          std::cout << "CRT hit y = " << _CRT_hit_y << " cm." << std::endl;
          std::cout << "CRT hit z = " << _CRT_hit_z << " cm." << std::endl;
        }

        break;
      }

    } // End of conditional for closest CRT hit time.

  } // End of loop over CRT hits.

  // Fill the tree.
  _tree->Fill();

  // Add the data products to the event.
  e.put(std::move(crthit_flash_assn_v));

  // Return true if you are within resolution.
  if (fuseAsFilter && _within_resolution == 1)
    return false;

  // Otherwise return false.
  return true;

} // End of the filter module


DEFINE_ART_MODULE(UBCRTCosmicFilter)
