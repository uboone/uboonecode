////////////////////////////////////////////////////////////////////////
// Class:       CRTFlashAna
// Plugin Type: analyzer (art v2_05_01)
// File:        CRTFlashAna_module.cc
//
// Generated at Fri Jul  6 13:52:52 2018 by Herbert Greenlee using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include <sstream>
#include "lardataobj/RawData/TriggerData.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "uboone/RawData/utils/DAQHeaderTimeUBooNE.h"
#include "uboone/CRT/CRTProducts/CRTHit.hh"
#include "uboone/RawData/utils/ubdaqSoftwareTriggerData.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h" 
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"

namespace crt {
  class CRTFlashAna;
}


class crt::CRTFlashAna : public art::EDAnalyzer {
public:
  explicit CRTFlashAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTFlashAna(CRTFlashAna const &) = delete;
  CRTFlashAna(CRTFlashAna &&) = delete;
  CRTFlashAna & operator = (CRTFlashAna const &) = delete;
  CRTFlashAna & operator = (CRTFlashAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

  // Methods.

  void add_algorithm(const std::string& algo);

private:

  // Declare member data here.

  // Fcl parameters.

  std::string fDAQTimeModuleLabel;
  std::string fTriggerModuleLabel;
  std::string fSWTriggerModuleLabel;
  std::string fFlashModuleLabel;
  std::string fCRTModuleLabel;
  double fFlashMinPE;

  // Event time descriptions.

  std::vector<std::string> fTimeDescrip;

  // Directories.

  art::TFileDirectory fTopDir;
  std::set<std::string> fAlgos;

  // Histograms.

  TH1F* fHgpsntp;                                         // GPS-NTP time difference.
  std::map<std::string, TH1F*> fHPE;                      // Flash PE.
  std::vector<std::map<std::string, TH1F*> > fHcrt0;      // Flash vs. CRT t0 time difference.
  std::vector<std::map<std::string, TH1F*> > fHcrt0x;     // Flash vs. CRT t0 time difference expanded.
  std::vector<std::map<std::string, TH1F*> > fHcrt0d;     // Flash vs. CRT t0 time difference detail.
  std::vector<std::map<std::string, TH1F*> > fHcrt0dd;    // Flash vs. CRT t0 time difference fine detail.
  std::vector<std::map<std::string, TH2F*> > fHcrt02d;    // Flash vs. CRT t0 time difference vs. Remainder.
  std::vector<std::map<std::string, TProfile*> > fHcrt0pr; // Flash vs. CRT t0 time difference vs. Remainder.
  std::map<std::string, TH1F*> fHcrt1;                    // Flash vs. CRT t1 time difference.
  std::map<std::string, TH1F*> fHcrt1x;                   // Flash vs. CRT t1 time difference expanded.
  std::map<std::string, TH1F*> fHcrt1d;                   // Flash vs. CRT t1 time difference detail.
  std::map<std::string, TH1F*> fHcrt1dd;                  // Flash vs. CRT t1 time difference fine detail.
  std::vector<std::map<std::string, TH1F*> > fHcrt10;     // CRT t1 vs. CRT t0 time difference.
  std::vector<std::map<std::string, TH1F*> > fHcrt10x;    // CRT t1 vs. CRT t0 time difference expanded.
  std::vector<std::map<std::string, TH1F*> > fHcrt10d;    // CRT t1 vs. CRT t0 time difference detail.
  std::vector<std::map<std::string, TH1F*> > fHcrt10dd;   // CRT t1 vs. CRT t0 time difference fine detail.
  std::vector<std::map<std::string, TH2F*> > fHcrt102d;   // CRT t1 vs. CRT t0 time difference vs. Remainder.
  std::vector<std::map<std::string, TProfile*> > fHcrt10pr; // CRT t1 vs. CRT t0 time difference vs. Remainder.
};

// Constructor.

crt::CRTFlashAna::CRTFlashAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p),
  fDAQTimeModuleLabel(p.get<std::string>("DAQTimeModuleLabel")),
  fTriggerModuleLabel(p.get<std::string>("TriggerModuleLabel")),
  fSWTriggerModuleLabel(p.get<std::string>("SWTriggerModuleLabel")),
  fFlashModuleLabel(p.get<std::string>("FlashModuleLabel")),
  fCRTModuleLabel(p.get<std::string>("CRTModuleLabel")),
  fFlashMinPE(p.get<double>("FlashMinPE")),
  fTopDir(art::ServiceHandle<art::TFileService>()->mkdir("CRTFlashAna", "CRT vs. Flash Timing"))
{
  // Add time descriptions.

  fTimeDescrip.push_back(std::string("NTP time from swizzler"));
  fTimeDescrip.push_back(std::string("GPS time from swizzler"));
  fTimeDescrip.push_back(std::string("Recalculated GPS time (no PPS round, no freq tune)"));
  fTimeDescrip.push_back(std::string("Recalculated GPS time (PPS round, no freq tune)"));
  fTimeDescrip.push_back(std::string("Recalculated GPS time (PPS round, freq tune)"));

  // Book histograms that don't depend on trigger algorithm.

  fHgpsntp = fTopDir.make<TH1F>("gpsntp", "GPS-NTP Time Difference", 200, -1., 1.);
  fHgpsntp->GetXaxis()->SetTitle("GPS-NTP Time Difference (sec)");

  // Add pseudoalgorithm "All."

  std::string all("All");
  add_algorithm(all);
}

void crt::CRTFlashAna::analyze(art::Event const & evt)
{
  // Get daq time header.

  art::Handle<raw::DAQHeaderTimeUBooNE> htime;
  evt.getByLabel(fDAQTimeModuleLabel, htime);
  if(!htime.isValid()) {
    std::cout << "No DAQHeaderTimeUBooNE data product with module label " 
	      << fDAQTimeModuleLabel << std::endl;
    return;
  }

  // Calculate various absolute event times.
  //
  // 0 - NTP time from swizzler.
  // 1 - GPS time from swizzler.
  // 2 - GPS time without PPS rounding to nearest second.
  // 3 - GPS time with PPS rounding, but without daq clock frequency tuning.
  // 4 - GPS time with PPS rounding and daq clock frequency tuning using hardwired parameter.

  std::vector<time_t> times;

  // NTP time from swizzler.

  times.push_back(htime->ntp_time());

  // GPS time from swizzler.

  times.push_back(htime->gps_time());

  // Recalculated GPS time (no PPS rounding or frequency tuning).

  long double gps0 = htime->pps_sec() +
                    1.e-6L * htime->pps_micro() +
                    1.e-9L * htime->pps_nano() +
                    1.6e-3L * (int(htime->trig_frame()) - int(htime->trig_pps_frame())) +
                    0.5e-6L * (int(htime->trig_sample()) - int(htime->trig_pps_sample())) +
                    0.0625e-6L * (int(htime->trig_div()) - int(htime->trig_pps_div()));
  uint64_t gps0_sec = gps0;
  uint64_t gps0_nsec = 1.e9L * (gps0 - gps0_sec);
  time_t gps0_time = (gps0_sec << 32) | gps0_nsec;
  times.push_back(gps0_time);

  // Recalculated GPS time (with PPS rounding, no frequency tuning).

  long double gps1 = htime->pps_sec() +
                    (htime->pps_micro() > 500000 ? 1.L : 0.L) +
                    1.6e-3L * (int(htime->trig_frame()) - int(htime->trig_pps_frame())) +
                    0.5e-6L * (int(htime->trig_sample()) - int(htime->trig_pps_sample())) +
                    0.0625e-6L * (int(htime->trig_div()) - int(htime->trig_pps_div()));
  uint64_t gps1_sec = gps1;
  uint64_t gps1_nsec = 1.e9L * (gps1 - gps1_sec);
  time_t gps1_time = (gps1_sec << 32) | gps1_nsec;
  times.push_back(gps1_time);

  // Recalculated GPS time (with PPS rounding and frequency tuning).

  long double tunefac = 1.000006882L;
  long double gps2 = htime->pps_sec() +
                    (htime->pps_micro() > 500000 ? 1.L : 0.L) +
                    tunefac * (1.6e-3L * (int(htime->trig_frame()) - int(htime->trig_pps_frame())) +
			       0.5e-6L * (int(htime->trig_sample()) - int(htime->trig_pps_sample())) +
			       0.0625e-6L * (int(htime->trig_div()) - int(htime->trig_pps_div())));
  uint64_t gps2_sec = gps2;
  uint64_t gps2_nsec = 1.e9L * (gps2 - gps2_sec);
  time_t gps2_time = (gps2_sec << 32) | gps2_nsec;
  times.push_back(gps2_time);

  if(times.size() != fTimeDescrip.size()) {
    throw cet::exception("CRTFlashAna") << "Number of times mismatched.";
  }

  // Calculate remainder time.

  double t_remainder = 1.6e-3L * (int(htime->trig_frame()) - int(htime->trig_pps_frame())) +
                       0.5e-6L * (int(htime->trig_sample()) - int(htime->trig_pps_sample())) +
                       0.0625e-6L * (int(htime->trig_div()) - int(htime->trig_pps_div()));
  std::cout << "Remainder = " << t_remainder << std::endl;

  // Print out the times we calculated.

  int n = 0;
  for(auto const& t : times) {
    unsigned int t_sec = (t >> 32);
    unsigned int t_nsec = (t & 0xffffffff);
    std::cout << "\n" << fTimeDescrip[n++] << "\n"
	      << "time = "<< t_sec << " seconds, " << t_nsec << " nanoseconds." << std::endl;
  }

  // Fill GPS-NTP time difference histogram.

  time_t gps = htime->gps_time();
  if(gps != 0) {
    time_t ntp = htime->ntp_time();
    double dt = (double(gps>>32) - double(ntp>>32)) +
      1.e-9 * (double(gps & 0xffffffff) - double(ntp & 0xffffffff));
    fHgpsntp->Fill(dt);
  }

  // Get hardware trigger information (not used except to print out).

  art::Handle<std::vector<raw::Trigger> > htrig;
  evt.getByLabel(fTriggerModuleLabel, htrig);
  if(!htrig.isValid()) {
    std::cout << "No Trigger data product with module label " 
	      << fTriggerModuleLabel << std::endl;
    return;
  }
  for(auto const& trig : *htrig) {
    double trigger_time = trig.TriggerTime();
    std::cout << "\nTrigger time = " << trigger_time << std::endl;
  }

  // Get Flashes.

  art::Handle<std::vector<recob::OpFlash> > hflash;
  evt.getByLabel(fFlashModuleLabel, hflash);

  // Get CRT hits.

  art::Handle<std::vector<crt::CRTHit> > hcrthit;
  evt.getByLabel(fCRTModuleLabel, hcrthit);
  if(!hcrthit.isValid())
    return;

  // Get software trigger information.

  std::string all("All");
  art::Handle<raw::ubdaqSoftwareTriggerData> hswtrig;
  evt.getByLabel(fSWTriggerModuleLabel, hswtrig);
  if(!hswtrig.isValid()) {
    std::cout << "No ubdaqSoftwareTriggerData data product with module label " 
	      << fSWTriggerModuleLabel << std::endl;
    return;
  }

  // Get software algorithms.

  auto const& swtrig = *hswtrig; 
  std::vector<std::string> algos = swtrig.getListOfAlgorithms();

  // Make a list of passed algorithms.
  // Add pseudoalgorithm "All" at the front of the list of passed algorithms.

  std::vector<std::string> passed_algos;
  passed_algos.push_back("All");
  for(auto const& algo : algos) {
    add_algorithm(algo);
    bool pass_algo = swtrig.passedAlgo(algo);
    bool pass_prescale = swtrig.passedPrescaleAlgo(algo);
    std::cout << algo << " " << pass_algo << " " << pass_prescale << std::endl;
    if(pass_algo && pass_prescale)
      passed_algos.push_back(algo);
  }

  // Loop over passed algorithms.

  for(auto const& algo : passed_algos) {

    // Loop over event times.

    for(unsigned int i=0; i<times.size(); ++i) {

      time_t tev = times[i];
      unsigned int tev_sec = (tev >> 32);
      unsigned int tev_nsec = (tev & 0xffffffff);

      // Make a double loop over CRT hits and flashes.
      // Outer loop over flashes.

      if(hflash.isValid()) {
	for(auto const& flash : *hflash) {
	  double flash_time = flash.Time();   // microseconds.
	  double flash_pe = flash.TotalPE();  // Total PE.

	  fHPE[algo]->Fill(flash_pe);

	  // Cut on flash pulse height.

	  if(flash_pe >= fFlashMinPE) {

	    // Inner loop over CRT hits.
  
	    for(auto const& crthit : *hcrthit) {

	    // Calculate CRT hit time relative to event time (units microseconds).

	      double crthit_t0 =
		1.e6 * (double(crthit.ts0_s) - double(tev_sec)) +
		1.e-3 * (double(crthit.ts0_ns) - double(tev_nsec));
	      double crthit_t1 = 1.e-3 * crthit.ts1_ns;

	      // Fill histograms.

	      fHcrt0[i][algo]->Fill(crthit_t0 - flash_time);
	      fHcrt0x[i][algo]->Fill(crthit_t0 - flash_time);
	      fHcrt0d[i][algo]->Fill(crthit_t0 - flash_time);
	      fHcrt0dd[i][algo]->Fill(crthit_t0 - flash_time);
	      fHcrt02d[i][algo]->Fill(t_remainder, crthit_t0 - flash_time);
	      fHcrt0pr[i][algo]->Fill(t_remainder, crthit_t0 - flash_time);

	      fHcrt1[algo]->Fill(crthit_t1 - flash_time);
	      fHcrt1x[algo]->Fill(crthit_t1 - flash_time);
	      fHcrt1d[algo]->Fill(crthit_t1 - flash_time);
	      fHcrt1dd[algo]->Fill(crthit_t1 - flash_time);
	    }
	  }
	}
      }

      // Make a single loop over crt hits to record t1-t0 time difference.

      for(auto const& crthit : *hcrthit) {

	// Calculate CRT hit time relative to event time (units microseconds).

	double crthit_t0 =
	  1.e6 * (double(crthit.ts0_s) - double(tev_sec)) +
	  1.e-3 * (double(crthit.ts0_ns) - double(tev_nsec));
	double crthit_t1 = 1.e-3 * crthit.ts1_ns;

	// Fill histograms.

	fHcrt10[i][algo]->Fill(crthit_t1 - crthit_t0);
	fHcrt10x[i][algo]->Fill(crthit_t1 - crthit_t0);
	fHcrt10d[i][algo]->Fill(crthit_t1 - crthit_t0);
	fHcrt10dd[i][algo]->Fill(crthit_t1 - crthit_t0);
	fHcrt102d[i][algo]->Fill(t_remainder, crthit_t1 - crthit_t0);
	fHcrt10pr[i][algo]->Fill(t_remainder, crthit_t1 - crthit_t0);
      }
    }
  }
}

void crt::CRTFlashAna::add_algorithm(const std::string& algo)
{
  if(fAlgos.count(algo) == 0) {
    fAlgos.insert(algo);
    art::TFileDirectory dir = fTopDir.mkdir(algo);

    // Flash pulse height.

    fHPE[algo] = dir.make<TH1F>("FlashPE", "Flash PE", 100, 0., 1000.);
    fHPE[algo]->GetXaxis()->SetTitle("Flash PE (ADC)");

    // Add t1 vs. flash time histograms.
    // These histograms don't depend on the absolute event time.
    // The only difference of these plots is the time scale.

    fHcrt1[algo] = dir.make<TH1F>("crt1", "CRT (t1) vs. Flash Time Difference",
				  1000, -5000., 5000.);
    fHcrt1[algo]->GetXaxis()->SetTitle("CRT (t1) Flash Time Difference (us)");

    fHcrt1x[algo] = dir.make<TH1F>("crt1x", "CRT (t1) vs. Flash Time Difference Expanded",
				  500, -500., 0.);
    fHcrt1x[algo]->GetXaxis()->SetTitle("CRT (t1) Flash Time Difference (us)");

    fHcrt1d[algo] = dir.make<TH1F>("crt1d", "CRT (t1) vs. Flash Time Difference Detail",
				  200, -75., -25.);
    fHcrt1d[algo]->GetXaxis()->SetTitle("CRT (t1) Flash Time Difference (us)");

    fHcrt1dd[algo] = dir.make<TH1F>("crt1dd", "CRT (t1) vs. Flash Time Difference Fine Detail",
				  200, -42.5, -37.5);
    fHcrt1dd[algo]->GetXaxis()->SetTitle("CRT (t1) Flash Time Difference (us)");

    // Loop over types of event times.

    for(unsigned int i=0; i<fTimeDescrip.size(); ++i) {

      // Add t0 vs. flash time histograms.
      // The only difference of these plots is the time scale.

      if(fHcrt0.size() <= i)
	fHcrt0.push_back(std::map<std::string, TH1F*>());
      {
	std::ostringstream name;
	std::ostringstream title;
	name << "crt0_t" << i;
	title << "CRT vs. Flash Time Difference, " << fTimeDescrip[i];
	fHcrt0[i][algo] = dir.make<TH1F>(name.str().c_str(), title.str().c_str(),
					 1000, -5000., 5000.);
	fHcrt0[i][algo]->GetXaxis()->SetTitle("CRT Flash Time Difference (us)");
      }

      if(fHcrt0x.size() <= i)
	fHcrt0x.push_back(std::map<std::string, TH1F*>());
      {
	std::ostringstream name;
	std::ostringstream title;
	name << "crt0x_t" << i;
	title << "CRT vs. Flash Time Difference Expanded, " << fTimeDescrip[i];
	fHcrt0x[i][algo] = dir.make<TH1F>(name.str().c_str(), title.str().c_str(), 500, -500., 0.);
	fHcrt0x[i][algo]->GetXaxis()->SetTitle("CRT Flash Time Difference (us)");
      }

      if(fHcrt0d.size() <= i)
	fHcrt0d.push_back(std::map<std::string, TH1F*>());
      {
	std::ostringstream name;
	std::ostringstream title;
	name << "crt0d_t" << i;
	title << "CRT vs. Flash Time Difference Detail, " << fTimeDescrip[i];
	fHcrt0d[i][algo] = dir.make<TH1F>(name.str().c_str(), title.str().c_str(),
					  200, -100., -50.);
	fHcrt0d[i][algo]->GetXaxis()->SetTitle("CRT Flash Time Difference (us)");
      }

      if(fHcrt0dd.size() <= i)
	fHcrt0dd.push_back(std::map<std::string, TH1F*>());
      {
	std::ostringstream name;
	std::ostringstream title;
	name << "crt0dd_t" << i;
	title << "CRT vs. Flash Time Difference Fine Detail, " << fTimeDescrip[i];
	fHcrt0dd[i][algo] = dir.make<TH1F>(name.str().c_str(), title.str().c_str(),
					   200, -70., -65.);
	fHcrt0dd[i][algo]->GetXaxis()->SetTitle("CRT Flash Time Difference (us)");
      }

       if(fHcrt02d.size() <= i)
	fHcrt02d.push_back(std::map<std::string, TH2F*>());
      {
	std::ostringstream name;
	std::ostringstream title;
	name << "crt02d_t" << i;
	title << "CRT vs. Flash Time Difference vs. Remainder, " << fTimeDescrip[i];
	fHcrt02d[i][algo] = dir.make<TH2F>(name.str().c_str(), title.str().c_str(),
					   20, 0., 1., 1000, -100., -50.);
	fHcrt02d[i][algo]->GetXaxis()->SetTitle("Remainder Time (sec)");
	fHcrt02d[i][algo]->GetYaxis()->SetTitle("CRT Flash Time Difference (us)");
      }

      if(fHcrt0pr.size() <= i)
	fHcrt0pr.push_back(std::map<std::string, TProfile*>());
      {
	std::ostringstream name;
	std::ostringstream title;
	name << "crt0pr_t" << i;
	title << "CRT vs. Flash Time Difference vs. remainder, " << fTimeDescrip[i];
	fHcrt0pr[i][algo] = dir.make<TProfile>(name.str().c_str(), title.str().c_str(),
					       20, 0., 1., -100., 100.);
	fHcrt0pr[i][algo]->GetXaxis()->SetTitle("Remainder Time (sec)");
	fHcrt0pr[i][algo]->GetYaxis()->SetTitle("CRT Flash Time Difference (us)");
      }

      // Add t0 vs. t1 histograms.
      // The only difference of these plots is the time scale.

      if(fHcrt10.size() <= i)
	fHcrt10.push_back(std::map<std::string, TH1F*>());
      {
	std::ostringstream name;
	std::ostringstream title;
	name << "crt10_t" << i;
	title << "CRT Hit t1 vs. t0 Time Difference, " << fTimeDescrip[i];
	fHcrt10[i][algo] = dir.make<TH1F>(name.str().c_str(), title.str().c_str(),
					  1000, -5000., 5000.);
	fHcrt10[i][algo]->GetXaxis()->SetTitle("CRT Hit Time Difference (us)");
      }

      if(fHcrt10x.size() <= i)
	fHcrt10x.push_back(std::map<std::string, TH1F*>());
      {
	std::ostringstream name;
	std::ostringstream title;
	name << "crt10x_t" << i;
	title << "CRT Hit t1 vs. t0 Time Difference Expanded, " << fTimeDescrip[i];
	fHcrt10x[i][algo] = dir.make<TH1F>(name.str().c_str(), title.str().c_str(),
					   200, -500., 500.);
	fHcrt10x[i][algo]->GetXaxis()->SetTitle("CRT Hit Time Difference (us)");
      }

      if(fHcrt10d.size() <= i)
	fHcrt10d.push_back(std::map<std::string, TH1F*>());
      {
	std::ostringstream name;
	std::ostringstream title;
	name << "crt10d_t" << i;
	title << "CRT Hit t1 vs. t0 Time Difference Detail, " << fTimeDescrip[i];
	fHcrt10d[i][algo] = dir.make<TH1F>(name.str().c_str(), title.str().c_str(), 200, 0., 50.);
	fHcrt10d[i][algo]->GetXaxis()->SetTitle("CRT Hit Time Difference (us)");
      }

      if(fHcrt10dd.size() <= i)
	fHcrt10dd.push_back(std::map<std::string, TH1F*>());
      {
	std::ostringstream name;
	std::ostringstream title;
	name << "crt10dd_t" << i;
	title << "CRT Hit t1 vs. t0 Time Difference FineDetail, " << fTimeDescrip[i];
	fHcrt10dd[i][algo] = dir.make<TH1F>(name.str().c_str(), title.str().c_str(), 200, 25., 30.);
	fHcrt10dd[i][algo]->GetXaxis()->SetTitle("CRT Hit Time Difference (us)");
      }

      if(fHcrt102d.size() <= i)
	fHcrt102d.push_back(std::map<std::string, TH2F*>());
      {
	std::ostringstream name;
	std::ostringstream title;
	name << "crt102d_t" << i;
	title << "CRT Hit t1 vs. t0 Time Difference vs. Remainder, " << fTimeDescrip[i];
	fHcrt102d[i][algo] = dir.make<TH2F>(name.str().c_str(), title.str().c_str(),
					    20, 0., 1., 1000, 0., 50.);
	fHcrt102d[i][algo]->GetXaxis()->SetTitle("Remainder Time (sec)");
	fHcrt102d[i][algo]->GetYaxis()->SetTitle("CRT Hit Time Difference (us)");
      }

      if(fHcrt10pr.size() <= i)
	fHcrt10pr.push_back(std::map<std::string, TProfile*>());
      {
	std::ostringstream name;
	std::ostringstream title;
	name << "crt10pr_t" << i;
	title << "CRT Hit t1 vs. t0 Time Difference vs. Remainder, " << fTimeDescrip[i];
	fHcrt10pr[i][algo] = dir.make<TProfile>(name.str().c_str(), title.str().c_str(),
						20, 0., 1., -100., 100.);
	fHcrt10pr[i][algo]->GetXaxis()->SetTitle("Remainder Time (sec)");
	fHcrt10pr[i][algo]->GetYaxis()->SetTitle("CRT Hit Time Difference (us)");
      }
    }
  }
}

DEFINE_ART_MODULE(crt::CRTFlashAna)
