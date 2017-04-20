////////////////////////////////////////////////////////////////////////
// Class:       ShowWire
// Module Type: analyzer
// File:        ShowWire_module.cc
//
// Generated at Thu Jul 31 15:07:10 2014 by Wesley Ketchum using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "larevt/CalibrationDBI/Interface/DetPedestalService.h"
#include "larevt/CalibrationDBI/Interface/DetPedestalProvider.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // geo::SigType_t
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcore/Geometry/Geometry.h"

#include "uboone/Utilities/SignalShapingServiceMicroBooNE.h"
#include "lardata/Utilities/LArFFT.h"


#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"


#include "lardata/ArtDataHelper/HitCreator.h"

#include "larreco/HitFinder/GausHitFinderAlg.h"


#include <string>
#include <sstream>
#include <iostream>


class VResponseCorrector : public art::EDProducer
{
public:
  // create calibrated signals on wires. this class runs
  // an fft to remove the electronics shaping.
  explicit VResponseCorrector(fhicl::ParameterSet const& pset);
  virtual ~VResponseCorrector();

  void produce(art::Event& evt);
  void beginJob();
  void endJob();
  void reconfigure(fhicl::ParameterSet const& p);


private:

  std::string    fHitModuleLabel;
  std::string    fWireModuleLabel;

  double         fSquaredDiffCutoff;

  double getSquareDiff(const std::vector<const recob::Hit *> & _roi_hits,
                       const lar::sparse_vector<float>::datarange_t & _datarange);

  float getPedestal(art::Ptr<raw::RawDigit> _raw_digit);

protected:

}; // class VResponseCorrector


DEFINE_ART_MODULE(VResponseCorrector)

VResponseCorrector::VResponseCorrector(fhicl::ParameterSet const& pset) {
  this->reconfigure(pset);

  // produces<art::Assns<raw::RawDigit, recob::Wire>>();
  // recob::HitCollectionCreator::declare_products(*this);
  // produces< std::vector<recob::Wire> >(fSpillName);
  // produces<art::Assns<raw::RawDigit, recob::Wire>>(fSpillName);
}

VResponseCorrector::~VResponseCorrector() {

}

void VResponseCorrector::reconfigure(fhicl::ParameterSet const& p) {

  fHitModuleLabel     = p.get< std::string > ("HitModuleLable", "gaushit");
  fWireModuleLabel    = p.get< std::string > ("WireModuleLable", "caldata");
  fSquaredDiffCutoff  = p.get<double>        ("SquaredDiffCutoff", 1.0);
}



void VResponseCorrector::produce(art::Event& evt) {

  //#############################################################
  // This module will try to correct regions in the V plane where
  // spurious hits are created due to deconvolution with the wrong
  // response function.  Generally, this is because the signals
  // have unipolar pulses and are deconvolved with bipolar response
  // functions.
  //
  // The procedure is:
  // 1) For every wire in the V plane, get the associated hits
  //    Find the regions of hits that have many hits in a row, or a very
  //    bad Chi^2 compared to the signal underneath.
  // 2) For the regions where there are determined to be bad hits, go
  //    back and get the raw digits and determine where the real hits are.
  // 3) Generally, the deconvolution failed on these regions which is
  //    why spurious hits are introduced.  So, re-deconvolve these regions
  //    with the unipolar response function, and refit them with hits.
  // 3a) For right now, just make a new collection of hits without the
  //     hits that are considered junk.
  //#############################################################

  // ##########################################
  // ### Container for the new hits produced ##
  // ##########################################
  // recob::HitCollectionCreator hcol(*this, evt);

  // ##########################################
  // ### Reading in the Wire List object(s) ###
  // ##########################################
  art::Handle< std::vector<recob::Wire> > wireVecHandle;
  evt.getByLabel(fWireModuleLabel, wireVecHandle);


  //#########################################
  // Get the hits associated to these wires #
  //#########################################

  art::InputTag assn_tag(fHitModuleLabel);
  art::FindMany<recob::Hit>
  hits_for_wire(wireVecHandle, evt, assn_tag);

  // ##############################################
  // Get the raw digits associated to these wires #
  // ##############################################
  art::FindOneP<raw::RawDigit> RawDigits
  (wireVecHandle, evt, fWireModuleLabel);

  // ################################
  // Get the signal shaping service #
  // ################################
  art::ServiceHandle<util::SignalShapingServiceMicroBooNE> sss;
  // re-initialize the FFT service for the request size
  art::ServiceHandle<util::LArFFT> fFFT;

  // ###################################
  // Get a place to store output wires #
  // ###################################
  // make a collection of Wires
  std::unique_ptr<std::vector<recob::Wire> > wirecol(new std::vector<recob::Wire>);
  // ... and an association set
  std::unique_ptr<art::Assns<raw::RawDigit, recob::Wire> >
  WireDigitAssn(new art::Assns<raw::RawDigit, recob::Wire>);




  // Loop over all the wires in the:
  for (size_t wireIter = 0; wireIter < wireVecHandle->size(); wireIter ++) {

    // Get the wire:
    auto wire = wireVecHandle -> at(wireIter);
    auto rawdigits = RawDigits.at(wireIter);

    // Get the plane of this wire.  If not v, continue:
    if (wire.View() != geo::kV) {
      // Make sure to add all the hits on
      // this wire to the new hit collection:
      continue;
    }

    // Got a V wire.  Get the hits associated with this wire:
    std::vector< const recob::Hit * > hits;
    hits_for_wire.get(wireIter, hits);

    if (hits.size() == 0) {
      continue;
    }

    // // Get the raw digit for this wire:
    // std::vector<const raw::RawDigit * > raw_digit;
    // digits_for_wire.get(wireIter, raw_digit);


    // Loop over the regions of interest, and find all the hits
    // that are associated with that region of interest.
    const recob::Wire::RegionsOfInterest_t& signalROI = wire.SignalROI();
    size_t i_roi = 0;
    std::vector<size_t> _bad_rois;
    std::vector<const recob::Hit *> _good_hits;
    for (const auto& range : signalROI.get_ranges())
    {
      const int _start_tick = range.begin_index();
      const int _end_tick = range.end_index();

      // Loop over the hits associated with this wire and find the ones
      // that are within the start and end tick of this ROI:
      std::vector< const recob::Hit *> _roi_hits;
      for (auto hit : hits) {
        if (hit ->StartTick() > _start_tick &&
            hit -> EndTick() < _end_tick) {
          _roi_hits.push_back(hit);
        }
      }

      // Next, compute the chi2 for the ROI and the hits within
      // that region:
      double _chi_sq = getSquareDiff(_roi_hits, range);

      if (_chi_sq >= fSquaredDiffCutoff) {
        std::cout << "ROI on wire " << wire.Channel() - 2400
                  << " over region [" << _start_tick << ", "
                  << _end_tick << "] has " << _roi_hits.size()
                  << " with a squared difference of "
                  << _chi_sq << "." << std::endl;

        _bad_rois.push_back(i_roi);

      }
      else {
        for (auto hit : _roi_hits) {
          _good_hits.push_back(hit);
        }
      }
      i_roi ++;

    } // for over ROI regions

    // Now we've found the bad regions of the wire.  If there are
    // none, copy the wire and hits and raw digits to
    // a new producer
    //
    // If there are bad regions, deconvolve the wire with a unipolar
    // response function, and re-find hits but only on the bad regions

    // std::vector < recob::Hit * > bad_hits = findBadHits(hits, wire);

    // if (bad_hits.size() != 0){
    //   std::vector<recob::Hit> remade_hits = replaceBadHits(bad_hits, raw_digit);
    // }

    recob::Wire::RegionsOfInterest_t ROIVec;

    if (_bad_rois.size() > 0) {


      size_t dataSize = rawdigits->Samples();
      auto channel = wire.Channel();
      size_t transformSize = 0;
      if (!transformSize)
      {
        sss->SetDecon(dataSize, channel);
        transformSize = fFFT->FFTSize();
      }

      std::vector<float> rawAdcLessPedVec;

      float pedestal = getPedestal(rawdigits);

      // uncompress the data
      std::vector<short> rawadc(dataSize);
      raw::Uncompress(rawdigits->ADCs(), rawadc, rawdigits->Compression());
      rawAdcLessPedVec.resize(transformSize, 0.);

      size_t binOffset(0);

      if (transformSize > dataSize) binOffset = (transformSize - dataSize) / 2;

      size_t startBin(binOffset);

      // Get the pedestal subtracted data, centered in the deconvolution vector
      std::transform(rawadc.begin(),
                     rawadc.end(),
                     rawAdcLessPedVec.begin() + startBin,
      [pedestal](const short & adc) {
        return std::round(float(adc) - pedestal);
      });
      std::fill(rawAdcLessPedVec.begin(),
                rawAdcLessPedVec.begin() + startBin,
                0.);    //rawAdcLessPedVec.at(startBin));
      std::fill(rawAdcLessPedVec.begin() + startBin + dataSize,
                rawAdcLessPedVec.end(),
                0.); //rawAdcLessPedVec.at(startBin+dataSize-1));


      std::cout << "Transform size: " << transformSize << std::endl;
      // Force the channel to be on the collection plane:
      sss->Deconvolute(5800, rawAdcLessPedVec);

      // Now, copy the rois onto the new wire collection
      // New roi ranges:
      for (_i_roi = 0; _i_roi < signalROI.get_ranges().size(); _i_roi ++){
        bool found = false;
        for (auto & _bad_roi : _bad_rois ) {
          if (_i_roi == _bad_roi) {
            std::vector<float> holder;

            holder.resize(roiLen);

            std::copy(rawAdcLessPedVec.begin() + binOffset + roi.first, rawAdcLessPedVec.begin() + binOffset + roi.second, holder.begin());
            std::transform(holder.begin(), holder.end(), holder.begin(), [deconNorm](float & deconVal) {return deconVal / deconNorm;});


            ROIVec.add_range(roi.first, std::move(holder));
            found = true;

          }
        }
        if (! found) {

        }
      }

    }

    // Finally, put the new hits and wires onto the event:
    evt.put(std::move(wirecol), fSpillName);
    evt.put(std::move(WireDigitAssn), fSpillName);


  } // end loop over wires

// hcol.put_into(evt);

}
void VResponseCorrector::beginJob() {

}
void VResponseCorrector::endJob() {

}

double VResponseCorrector::getSquareDiff(
  const std::vector<const recob::Hit *> & _roi_hits,
  const lar::sparse_vector<float>::datarange_t & _datarange) {
  // hit.peak_amplitude() * np.exp( - 0.5 * (xPts - hit.peak_time())**2 / hit.rms()**2  )

  // Initialize a vector to contain the waveform of the hits:
  std::vector<float> _hit_waveform(_datarange.size(), 0.0);
  for (auto hit : _roi_hits) {
    // Add the exponential waveform for this particular hit across
    // the entire region of this roi:
    for (size_t i = 0; i < _hit_waveform.size(); i ++) {
      float val = hit ->PeakAmplitude();
      val *= exp(-0.5 * pow(i - ( hit->PeakTime() - _datarange.begin_index() ), 2) / pow(hit->RMS(), 2) );
      _hit_waveform[i] +=  val;
    }
  }

  // Now, compare the two waveforms:
  double _squared_diff = 0.0;
  for (size_t i = 0; i < _hit_waveform.size(); i ++) {
    _squared_diff += pow(_hit_waveform[i] - _datarange.data()[i], 2);
  }
  return _squared_diff / _hit_waveform.size();
}

float VResponseCorrector::getPedestal(art::Ptr<raw::RawDigit> _raw_digit) {

  std::vector<float> _vals;

  _vals.reserve(101);


  for (size_t i = 0; i < _raw_digit -> Samples(); i += 80) {
    _vals.push_back(_raw_digit->ADC(i));
    if (_vals.size() == 101) {
      break;
    }
  }


  // This way gives a 25% decrease in execution time:
  std::nth_element(_vals.begin(), _vals.begin() + _vals.size() / 2, _vals.end());
  return _vals[_vals.size() / 2];


}