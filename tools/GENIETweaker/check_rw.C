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

constexpr int NUM_BINS = 20;
constexpr double MIN_ENERGY = 0.; // GeV
constexpr double MAX_ENERGY = 7.; // GeV

// Perform a two-sample chi2 test using two histograms.
// See section 3 of inspirehep.net/record/775602/files/ACAT_060.pdf for details.
void chi2_test(TH1* hist, TH1* weighted_hist, bool exclude_problem_bins = false) {
  // Do a chi-squared test to see if the reweighted results are statistically consistent
  // with the tweaked ones
  int num_bins = hist->GetNbinsX();
  assert( num_bins == weighted_hist->GetNbinsX() );

  double N = hist->Integral();
  double W = weighted_hist->Integral();
  double W2 = std::pow(W, 2);

  int num_dof = num_bins - 1.;
  double chi2 = 0.;
  for ( int i = 1; i <= num_bins; ++i ) {
    double n_i = hist->GetBinContent( i );
    double w_i = weighted_hist->GetBinContent( i );
    double s2_i = std::pow( weighted_hist->GetBinError(i), 2 );

    double p_i = (W*w_i - N*s2_i + std::sqrt( std::pow(W*w_i - N*s2_i, 2)
      + 4.0*W2*s2_i*n_i )) / ( 2.0*W2 );

    // Protect against NaNs for empty bins by skipping them and
    // adjusting the degrees of freedom accordingly
    if ( w_i + n_i <= 0. ) {
      --num_dof;
      continue;
    }

    if ( exclude_problem_bins ) {
      if ( N*p_i  < 1.0 ) {
        std::cout << "Excluding bin #" << i << " because the unweighted histogram"
          << " has an estimated expected number of events less than 1.\n";
        --num_dof;
        continue;
      }
      if ( std::pow(w_i, 2) / s2_i < 10. ) {
        std::cout << "Excluding bin #" << i << " because the weighted histogram"
          << " has an effective number of events less than 10.\n";
        --num_dof;
        continue;
      }
    }

    chi2 += std::pow(n_i - N*p_i, 2) / (N * p_i) + std::pow(w_i - W*p_i, 2) / s2_i;
  }

  double p_value = ROOT::Math::chisquared_cdf_c( chi2, num_dof );
  std::cout << "chi2 / dof = " << chi2 << " / " << num_dof << '\n';
  std::cout << "p-value: " << p_value << '\n';
}

void check_rw(const std::string& tweaked_file, const std::string& rw_file)
{
  std::vector<std::string> tweaked_files { tweaked_file };
  std::vector<std::string> rw_files { rw_file };

  TH1D* tweak_histogram = new TH1D("Ev_tweaked", "E_{#nu} ; E_{#nu} (GeV); counts",
    NUM_BINS, MIN_ENERGY, MAX_ENERGY);
  tweak_histogram->SetLineColor(kBlack);
  tweak_histogram->SetLineWidth(2);

  TH1D* untweaked_histogram = new TH1D("Ev_untweaked", "E_{#nu} ; E_{#nu} (GeV); counts",
    NUM_BINS, MIN_ENERGY, MAX_ENERGY);
  untweaked_histogram->SetLineColor(kRed);
  untweaked_histogram->SetLineWidth(2);

  TH1D* rw_histogram = new TH1D("Ev_rw", "E_{#nu} ; E_{#nu} (GeV); counts",
    NUM_BINS, MIN_ENERGY, MAX_ENERGY);
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
      untweaked_histogram->Fill( Ev );
    }
  }

  chi2_test(tweak_histogram, rw_histogram, true);
  tweak_histogram->Chi2Test(rw_histogram, "uw p");

  // For plotting, normalize the untweaked and reweighted histograms to the
  // same number of events as the tweaked histogram
  rw_histogram->Scale( tweak_histogram->GetEntries() / rw_histogram->Integral() );
  untweaked_histogram->Scale( tweak_histogram->GetEntries() / untweaked_histogram->GetEntries() );

  // Set the maximum y-axis value high enough so that we can see all three histograms
  double ymax = 1.1 * std::max( { tweak_histogram->GetMaximum(), rw_histogram->GetMaximum(),
    untweaked_histogram->GetMaximum() } );
  tweak_histogram->GetYaxis()->SetRangeUser(0., ymax);

  // Draw the histograms on the same plot
  TCanvas* c = new TCanvas;
  tweak_histogram->Draw("hist");
  untweaked_histogram->Draw("same hist");
  rw_histogram->Draw("same hist");
}
