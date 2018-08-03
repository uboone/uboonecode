#ifndef UBOPTICALADC_CXX
#define UBOPTICALADC_CXX

#include "UBOpticalADC.h"

namespace opdet {

  //----------------------------------------
  UBOpticalADC::UBOpticalADC() : UBADCBase()
  //----------------------------------------
  {
    Reset();
  }

  //------------------------
  void UBOpticalADC::Reset()
  //------------------------
  {
    UBADCBase::Reset();
    fSPE.Reset();
    fPED.Reset();
    fInputPhotonTime.clear();
    fDarkPhotonTime.clear();
    fPhotonTime.clear();
  }

  //---------------------------------------------------------------------
  void UBOpticalADC::SetPhotonsStatus(const std::vector<bool>& setinsd_v)
  //---------------------------------------------------------------------
  {
    fInputPhotonStatus.clear();
    fInputPhotonStatus.reserve(setinsd_v.size());
    for (auto const &v : setinsd_v) fInputPhotonStatus.push_back(v);
  }

  //---------------------------------------------------------------------
  void UBOpticalADC::SetPhotonsEnergy(const std::vector<float>& energy_v)
  //---------------------------------------------------------------------
  {
    fInputPhotonEnergy.clear();
    fInputPhotonEnergy.reserve(energy_v.size());
    for (auto const &v : energy_v) fInputPhotonEnergy.push_back(v);
  }

  //--------------------------------------------------------------
  void UBOpticalADC::SetPhotons(const std::vector<double>& g4time)
  //--------------------------------------------------------------
  {
    fInputPhotonTime.clear();
    fInputPhotonTime.reserve(g4time.size());
    for(auto const &v : g4time) fInputPhotonTime.push_back(v);
    /*
    double tmin=1e12;
    double tmax=0;
    size_t before_npe=0;
    size_t after_npe=0;
    size_t early_npe=0;
    size_t late_npe=0;
    size_t total_npe=0;
    for(auto const& v : g4time) {
      if(tmin > v) tmin = v;
      if(tmax < v) tmax = v;
    }
    tmin -=1;
    tmax +=1;
    for(auto const& v : g4time) {
      if(v < tmin) before_npe ++;
      if(v > tmin && v < tmin+1000) early_npe++;
      if(v > tmin+1000 && v < tmax) late_npe++;
      if(v > tmin && v < tmax) total_npe++;
      if(v < tmax ) after_npe++;
    }
    std::cout<< "Photons: " << total_npe << " ( "<<before_npe<<" => "<<early_npe<<" / "<<late_npe<<" => " <<after_npe<<" ) ... "<<tmin<<" => "<<tmax<<std::endl;
    */
  }

  //--------------------------------------------------------------------------
  //void UBOpticalADC::GenDarkNoise(double dark_rate, double period)
  void UBOpticalADC::GenDarkNoise(const unsigned int pmtid, const double g4start)
  //--------------------------------------------------------------------------
  {
    fDarkPhotonTime.clear();

    art::ServiceHandle<opdet::UBOpticalChConfig> ch_conf;
    art::ServiceHandle<geo::Geometry> geom;
    unsigned int ch = geom->OpChannel( pmtid, 0 ); // get channel reading out that PMT

    double dark_rate = ch_conf->GetFloat(kDarkRate,ch);

    unsigned int dark_count = RandomServer::GetME().Poisson(dark_rate * fDuration);

    fDarkPhotonTime.reserve(dark_count);

    for(size_t i=0; i<dark_count; ++i)

      fDarkPhotonTime.push_back(RandomServer::GetME().Uniform(fDuration*1.e3) + g4start);

  }

  //-------------------------------------------------------------------
  void UBOpticalADC::GenWaveform(const unsigned int ch, optdata::ChannelData& wf )
  //-------------------------------------------------------------------
  {
    //
    // Initialize, zero
    //
    size_t nticks = fTimeInfo.Ticks(fDuration);
    std::vector< float > wfm_tmp( nticks, 0.0 );
    if ( wf.size()!=nticks )
      wf.resize( nticks, 0 );
    else
      wf.assign( nticks, 0 );

    // services
    art::ServiceHandle<opdet::UBOpticalChConfig> ch_conf;
    //art::ServiceHandle<detinfo::LArPropertiesService> lar_prop;
auto const* lar_prop = lar::providerFrom<detinfo::LArPropertiesService>();
    //
    // Generate Signal & DarkNoise
    //
    // Configure to generate high gain SPE
    fSPE.Reset();

    fSPE.SetT0(ch_conf->GetFloat(kT0,ch),
	       ch_conf->GetFloat(kT0Spread,ch));
    /*
    if(ch<32) 
      std::cout<<"Gain: "<<ch_conf->GetFloat(kPMTGain,ch)<< " +/- "<< ch_conf->GetFloat(kGainSpread,ch)<<std::endl;
    */    
    fSPE.SetGain(ch_conf->GetFloat(kPMTGain,ch),
		 ch_conf->GetFloat(kGainSpread,ch));
    
    // Create combined photon time with QE applied on signal photons
    /*
    if(ch<32)
      std::cout<<"Channel: "<<ch<<" #photon: "<<fInputPhotonTime.size()<<std::endl;
    */
    fPhotonTime.clear();
    fPhotonTime.reserve(fInputPhotonTime.size() + fDarkPhotonTime.size());
    //const double qe = ch_conf->GetFloat(kQE,ch);

    for(auto const &v : fDarkPhotonTime) fPhotonTime.push_back(v);
    //std::cout << "[UBOpticalADC] >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> qe is set to " << qe  << " for ch " << ch  << std::endl;

    double mycounter = 0.;

    for(size_t i = 0; i < fInputPhotonTime.size(); i++) {

      double v = fInputPhotonTime.at(i);               // The photon time
      double set_in_sd = fInputPhotonStatus.at(i);     // False if it comes from photon library, True otherwise
      double energy = fInputPhotonEnergy.at(i) * 1.e6; // The photon enrgy [eV]

      //double qe = ch_conf->GetFloat(kQE,ch);
      //if (set_in_sd) {
      //  qe *= lar_prop->ScintPreScale();
      //}

      //std::cout << "Photon " << (set_in_sd ? "set" : "not set") << " in SD. Energy is " << energy << ". QE is " << ch_conf->GetFloatAtSpectrumValue(kQESpec, ch, energy) << std::endl;

      // Get the quantum efficiency for this channel (ch) and this photon energy
      double qe = ch_conf->GetFloatAtSpectrumValue(kQESpec, ch, energy);

      if (!set_in_sd) {
        // This is a scintillation photon propagated with optical library
        // A QE=ScintPreScale() has already been applied to these photons
        qe /= lar_prop->ScintPreScale();
      }

      if (v > 3911.87 - 100 && v < 3911.87 + 8000)
        std::cout << "[UBOpticalADC] >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> at our photon, qe is  " << qe << std::endl;

      // Apply quantum efficiency
      if(RandomServer::GetME().Uniform(1.) < qe) {
        if (v > 3911.87 - 100 && v < 3911.87 + 8000) {
          std::cout << "[UBOpticalADC] >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> >> passed " << std::endl;
          mycounter++;
        }
        fPhotonTime.push_back(v);
      }
    }

    std::cout << "[UBOpticalADC] >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> >> At the end we have mycounter = "<< mycounter << std::endl;

    fSPE.SetPhotons(fPhotonTime);
    fSPE.Process(wfm_tmp,fTimeInfo);
    // convert from pe waveform to adc
    /*
    double gain_ratio = ch_conf->GetFloat(kSplitterGain,ch);
    for(auto &v : wfm_tmp) 
      v *= gain_ratio;
    */

    //
    // Simulate pedestal
    //
    fPED.Reset();
    //if(ch<32) std::cout<<"Pedestal: "<<ch_conf->GetFloat(kPedestalMean,ch)<<" +/- "<<ch_conf->GetFloat(kPedestalSpread,ch)<<std::endl;
    fPED.SetPedestal(ch_conf->GetFloat(kPedestalMean,ch),
		     ch_conf->GetFloat(kPedestalSpread,ch));
    fPED.Process(wfm_tmp,fTimeInfo);
    
    // Make sure algorithms did not alter the waveform size

    if(wf.size()!=nticks || wfm_tmp.size()!=nticks)
      throw UBOpticalException("Waveform length changed (prohibited)!");

    //
    // Digitize amplitude
    //
    wf.clear();
    Digitize(wfm_tmp,wf);
    
  }
  
}

#endif
