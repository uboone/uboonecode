//======================================================================
//
// Name: SparseRawDigitTest_module.cc
//
// Purpose: Analyzer module to test SparseRawDigit data product.
//          This module accepts a single fcl parameter, which is the
//          label of a Wire collection data product.  This module
//          fetches the Wire collection and then uses associations
//          to fetch associated SparseRawDigits and regular RawDigits.
//          It performs the following checks.
//
//          1.  Checks that ROIs are the same between Wire and
//              SparseRawdigit.
//
//         2.   Checks that the SparseRawDigit waveform is the same
//              as the RawDigit waveform.
//
// FCL parameters:
//
// CalDataModuleLabel - CalData module label.
//
// Created: 19-Nov-2017  H. Greenlee
//
//======================================================================

#include <iostream>
#include <cassert>
#include <string>

#include "art/Framework/Core/ModuleMacros.h" 
#include "art/Framework/Core/EDAnalyzer.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "cetlib_except/exception.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RawData/RawDigit.h"
#include "ubooneobj/SparseRawDigit.h"

namespace raw{

  class SparseRawDigitTest : public art::EDAnalyzer
  {
  public:

    // Constructor, destructor.

    SparseRawDigitTest(fhicl::ParameterSet const& pset);
    ~SparseRawDigitTest();

    // Overrides.

    void analyze(const art::Event& evt);

  private:

    // Data members.

    std::string fCalDataModuleLabel;

  };

  DEFINE_ART_MODULE(SparseRawDigitTest)

  // Constructor.

  SparseRawDigitTest::SparseRawDigitTest(fhicl::ParameterSet const& pset) :
    art::EDAnalyzer(pset),
    fCalDataModuleLabel(pset.get<std::string>("CalDataModuleLabel"))
  {
    std::cout << "Constructed SparseRawDigitTest.\n"
	      << "CalDataModuleLabel: " << fCalDataModuleLabel << std::endl;
  }

  // Destructor.

  SparseRawDigitTest::~SparseRawDigitTest()
  {}

  // Analyze event.

  void SparseRawDigitTest::analyze(const art::Event& evt)
  {
    // Fetch Wire data product.

    art::Handle< std::vector<recob::Wire> > wireh;
    evt.getByLabel(fCalDataModuleLabel, wireh);
    if(!wireh.isValid()) {
      throw cet::exception("SparseRawDigitTest") << "Wire data product not found.";
    }

    std::cout << "Number of wires: " << wireh->size() << std::endl;

    // Construct FindOneP objects for RawDigit and SparseRawDigit.

    art::FindOneP<raw::RawDigit> find_raw_digit(wireh, evt, fCalDataModuleLabel);
    art::FindOneP<raw::SparseRawDigit> find_sparse_raw_digit(wireh, evt, fCalDataModuleLabel);

    // Loop over wires.

    unsigned int iwire = 0;
    for(auto const& wire : *wireh) {

      unsigned int channel = wire.Channel();
      std::cout << "Wire Channel " << channel << std::endl;

      // Get associated RawDigit.

      const art::Ptr<raw::RawDigit>& digit_ptr = find_raw_digit.at(iwire);
      if(digit_ptr.isNull()) {
	throw cet::exception("SparseRawDigitTest") << "Associated RawDigit not found.";
      }
      if(digit_ptr->Channel() != channel) {
	throw cet::exception("SparseRawDigitTest") << "RawDigit channel does not match Wire.";
      }

      // Get associated SparseRawDigit.

      const art::Ptr<raw::SparseRawDigit>& sparse_ptr = find_sparse_raw_digit.at(iwire);
      if(sparse_ptr.isNull()) {
	throw cet::exception("SparseRawDigitTest") << "Associated SparseRawDigit not found.";
      }
      if(sparse_ptr->Channel() != channel) {
	throw cet::exception("SparseRawDigitTest") << "SparseRawDigit channel does not match Wire.";
      }

      // Check ROIs.

      auto const& wire_ranges = wire.SignalROI().get_ranges();
      auto const& sparse_ranges = sparse_ptr->ADCs().get_ranges();

      std::cout << "Number of Wire ROIs = " << wire_ranges.size() << std::endl;
      std::cout << "Number of SparseRawDigit ROIs = " << sparse_ranges.size() << std::endl;
      if(wire_ranges.size() != sparse_ranges.size()) {
	throw cet::exception("SparseRawDigitTest") << "Mismatch number of ROIs";
      }

      // Loop over ROIs.

      unsigned int irange = 0;
      for(auto const& wire_range : wire_ranges) {
	auto const& sparse_range = sparse_ranges[irange];

	std::cout << "Range " << irange << std::endl;
	std::cout << "  Wire begin index = " << wire_range.begin_index() << std::endl;
	std::cout << "  Wire end index = " << wire_range.end_index() << std::endl;
	std::cout << "  SparseRawDigit begin index = " << sparse_range.begin_index() << std::endl;
	std::cout << "  SparseRawDigit end index = " << sparse_range.end_index() << std::endl;
	if(wire_range.begin_index() != sparse_range.begin_index()) {
	  throw cet::exception("SparseRawDigitTest") << "Mismatch begin index";
	}
	if(wire_range.end_index() != sparse_range.end_index()) {
	  throw cet::exception("SparseRawDigitTest") << "Mismatch end index";
	}

	// Compare waveforms in this ROI.

	for(unsigned int i = sparse_range.begin_index(); i < sparse_range.end_index(); ++i) {
	  std::cout << "SparseRawDigit adc = " << sparse_ptr->ADC(i) << std::endl;
	  std::cout << "RawDigit adc = " << digit_ptr->ADC(i) << std::endl;
	  if(sparse_ptr->ADC(i) != digit_ptr->ADC(i)) {
	    throw cet::exception("SparseRawDigitTest") << "Mismatch waveform";
	  }
	}
	++irange;
      }
      ++iwire;
    }
  }

} // namespace raw
