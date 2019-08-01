#include <cassert>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"

void check_rw(const std::string& tweaked_file, const std::string& rw_file) {

  std::vector<std::string> tweaked_files { tweaked_file };
  std::vector<std::string> rw_files { rw_file };

  TH1D* tweak_histogram = new TH1D("Ev_tweaked", "E_{#nu} ; E_{#nu} (GeV); counts", 20, 0., 7.);
  tweak_histogram->SetLineColor(kBlack);
  tweak_histogram->SetLineWidth(2);

  TH1D* rw_histogram = new TH1D("Ev_rw", "E_{#nu} ; E_{#nu} (GeV); counts", 20, 0., 7.);
  rw_histogram->SetLineColor(kBlue);
  rw_histogram->SetLineWidth(2);

  // Set up input tags
  art::InputTag mc_truths_tag { "generator" };
  art::InputTag mc_weights_tag { "weighttweak" };

  // Tweaked event loop
  for ( gallery::Event ev( tweaked_files ) ; !ev.atEnd(); ev.next() ) {
    auto const& mc_truths = *ev.getValidHandle< std::vector<simb::MCTruth> >(
      mc_truths_tag);

    for ( size_t i = 0; i < mc_truths.size(); ++i ) {
      auto const& truth = mc_truths.at(i);

      double Ev = truth.GetNeutrino().Nu().E();

      std::cout << "Ev = " << Ev << " GeV\n";
      tweak_histogram->Fill( Ev );
    }
  }

  // Reweighted event loop
  for ( gallery::Event ev( rw_files ) ; !ev.atEnd(); ev.next() ) {
    auto const& mc_truths = *ev.getValidHandle< std::vector<simb::MCTruth> >(
      mc_truths_tag);

    auto const& mc_weights = *ev.getValidHandle< std::vector<evwgh::MCEventWeight> >(
      mc_weights_tag);

    assert( mc_truths.size() == mc_weights.size() );

    for ( size_t i = 0; i < mc_truths.size(); ++i ) {
      auto const& truth = mc_truths.at(i);
      auto const& mc_weight = mc_weights.at(i);

      double Ev = truth.GetNeutrino().Nu().E();

      std::map<std::string, std::vector<double> > weight_map = mc_weight.fWeight;

      const std::vector<double>& weights = weight_map.at("my_tweak_Genie");

      assert( weights.size() == 1u );
      double weight = weights.front();

      std::cout << "Ev = " << Ev << " GeV, weight = " << weight << '\n';

      rw_histogram->Fill(Ev, weight);
    }
  }

  double ymax = std::max( tweak_histogram->GetMaximum(), rw_histogram->GetMaximum() )*1.1;
  tweak_histogram->GetYaxis()->SetRangeUser(0., ymax);

  TCanvas* c = new TCanvas;
  tweak_histogram->Draw("hist");
  rw_histogram->Draw("same hist");

}
