//Misc. Headers
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <queue>

//Pythia headers
#include "Pythia.h"

//Fastjet headers
#include "FastJet3.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"

//Root headers
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TNtuple.h"
#include "TArc.h"
#include "TVirtualFitter.h"
#include "TLine.h"

using namespace std;
using namespace Pythia8;
using namespace fastjet;
using namespace fastjet::contrib;

// class GeneralERecombiner : public fastjet::JetDefinition::Recombiner {
// public:
// // Constructor to choose value of alpha (defaulted to 1 for normal pT sum)
//   GeneralERecombiner(double delta) : _delta(delta) {}
  
//   std::string description() const {
//      return "General E-scheme recombination";
//   }
  
//   void recombine(const fastjet::PseudoJet & pa, const fastjet::PseudoJet & pb, fastjet::PseudoJet & pab) const {

//     double weighta = pow(pa.perp(), _delta);
//     double weightb = pow(pb.perp(), _delta);

//     double perp_ab = pa.perp() + pb.perp();
//     if (perp_ab != 0.0) { // weights also non-zero...
//       double y_ab = (weighta * pa.rap() + weightb * pb.rap())/(weighta+weightb);
      
//       // take care with periodicity in phi...
//       double phi_a = pa.phi(), phi_b = pb.phi();
//       if (phi_a - phi_b > pi)  phi_b += twopi;
//       if (phi_a - phi_b < -pi) phi_b -= twopi;
//       double phi_ab = (weighta * phi_a + weightb * phi_b)/(weighta+weightb);
      
//       pab.reset_PtYPhiM(perp_ab, y_ab, phi_ab);

//     }
//     else { // weights are zero
//       pab.reset(0.0,0.0,0.0,0.0);
//     }

//   }

// private:
//   double _delta;
// };

int next_comb(int comb[], int k, int n) {
  int i = k - 1;
  ++comb[i];
  while ((i >= 0) && (comb[i] >= n - k + 1 + i)) {
    --i;
    ++comb[i];
  }

  if (comb[0] > n - k) 
    return 0; 

  for (i = i + 1; i < k; ++i)
    comb[i] = comb[i - 1] + 1;

  return 1;
}

double triangle_perimeter(PseudoJet axis1, PseudoJet axis2, PseudoJet axis3) {
  double line1 = TMath::Sqrt(TMath::Power(axis1.eta() - axis2.eta(), 2) + TMath::Power(axis1.phi() - axis2.phi(), 2));
  double line2 = TMath::Sqrt(TMath::Power(axis1.eta() - axis3.eta(), 2) + TMath::Power(axis1.phi() - axis3.phi(), 2));
  double line3 = TMath::Sqrt(TMath::Power(axis2.eta() - axis3.eta(), 2) + TMath::Power(axis2.phi() - axis3.phi(), 2));
  return line1 + line2 + line3;
}

struct JetInformation {
  double tau;
  vector<PseudoJet> axes;
  vector<PseudoJet> jets;
};

JetInformation findMinAxes(vector<PseudoJet> input_particles, vector<PseudoJet> starting_axes, int njettiness, double beta, double Rcutoff) {

  vector<PseudoJet> min_manual_axes;
  JetInformation min_manual_set;

  std::string bitmask(njettiness, 1); // K leading 1's
  bitmask.resize(starting_axes.size(), 0); // N-K trailing 0's

  NjettinessPlugin njet_plugin_manual(njettiness, Manual_Axes(), UnnormalizedCutoffMeasure(beta, Rcutoff));
  JetDefinition njet_def_manual(&njet_plugin_manual);

  double min_manual_tau = std::numeric_limits<double>::max();
  do {
    vector<int> axis_indices;
    vector<PseudoJet> six_axes;
    for (int i = 0; i < starting_axes.size(); ++i) {
      if (bitmask[i]) axis_indices.push_back(i);
    }
    for (int j = 0; j < axis_indices.size(); j++) {
      six_axes.push_back(starting_axes[axis_indices[j]]);
    }

    double manual_tau = std::numeric_limits<double>::max();
    if (six_axes.size() == njettiness) {
      njet_plugin_manual.setAxes(six_axes);
      ClusterSequence njet_cluster_manual(input_particles, njet_def_manual);
      const NjettinessExtras *extras_manual = njettiness_extras(njet_cluster_manual);
      manual_tau = extras_manual->totalTau();

      if (manual_tau < min_manual_tau) {
        min_manual_tau = manual_tau;
        min_manual_set.tau = manual_tau;
        min_manual_set.axes = extras_manual->axes();
        min_manual_set.jets = extras_manual->jets();
      }
    }
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  return min_manual_set;
}

// vector<PseudoJet> findMinTriangle(vector<PseudoJet> starting_axes) {
//   // code for finding two triangles with minimum perimeter 

//   double min_perimeter = 1000;
//   PseudoJet min_ax1(0,0,0,0);
//   PseudoJet min_ax2(0,0,0,0);
//   PseudoJet min_ax3(0,0,0,0);

//   int n = 6; 
//   int k = 3; 
//   int comb[3]; 

//   for (int i = 0; i < k; ++i) {
//     comb[i] = i;
//   }

//   do {
//     double temp_perimeter = triangle_perimeter(starting_axes[comb[0]], starting_axes[comb[1]], starting_axes[comb[2]]);
//     if (temp_perimeter < min_perimeter) {
//       min_ax1 = starting_axes[comb[0]]; 
//       min_ax2 = starting_axes[comb[1]]; 
//       min_ax3 = starting_axes[comb[2]]; 
//       min_perimeter = temp_perimeter; 
//     }
//   } while (next_comb(comb, k, n));

//   vector<PseudoJet> min_axes;
//   min_axes.push_back(min_ax1);
//   min_axes.push_back(min_ax2);
//   min_axes.push_back(min_ax3);

//   return min_axes;
// }

// vector<PseudoJet> findSecondTriangle(vector<PseudoJet> starting_axes, vector<PseudoJet> first_triangle) {
//   vector<PseudoJet> second_triangle;
//   double epsilon = 0.0001;

//   // find the second triangle
//   for (int j = 0; j < 6; j++) {
//     if (starting_axes[j].delta_R(first_triangle[0]) > epsilon
//       && starting_axes[j].delta_R(first_triangle[1]) > epsilon
//       && starting_axes[j].delta_R(first_triangle[2]) > epsilon) {
//       second_triangle.push_back(starting_axes[j]);
//     }
//   }

//   return second_triangle;
// }

double calcMedian(vector<double> values)
{
  double median;
  int size = values.size();

  sort(values.begin(), values.end());

  if (size  % 2 == 0)
  {
      median = (values[size / 2 - 1] + values[size / 2]) / 2;
  }
  else 
  {
      median = values[size / 2];
  }

  return median;
}

int main(int argc, char* argv[]){

  // TFile out("njettiness_testing_10choose6_wtap2_beta1.root", "RECREATE");
  TFile out("axestest_antiktjets.root", "RECREATE");

  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW"); // Input File Name
  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW"); // Input File Name

  TStyle *plain  = new TStyle("Plain","plain");
  plain->SetLegendBorderSize(0);
  plain->SetLegendFillColor(0);
  plain->SetTitleFont(132, "a");
  plain->SetTitleFont(132, "xy");
  plain->SetLegendFont(132);
  // plain->SetTextSize(0.05);
  plain->SetLabelSize(0.04, "xy");
  // plain->SetLabelOffset(0.003, "xy");
  plain->SetTitleSize(0.06, "a");
  plain->SetTitleSize(0.05, "xy");
  plain->SetPadLeftMargin(0.12);
  plain->SetPadBottomMargin(0.12);
  plain->SetTitleOffset(1.25, "y");
  plain->SetHistLineWidth(3);

  // plain->SetLegendSize(12);
  plain->SetTitleBorderSize(0);
  plain->SetTitleX(0.1f);
  plain->SetTitleW(0.8f);
  plain->SetOptStat(0);
  plain->cd();

  gROOT->SetStyle("plain");

  double Rcutoff = 0.6;
  double power;
  double delta;

  int njettiness = 6;
  int initial_axes_num = njettiness;

  vector<double> mass_ratios;
  vector<double> mean_values;
  vector<double> rms_values;

  const int n_event = 1000; // # of events to be analyzed   
  const int n_powers = 32;
  const int n_deltas = 25;
  // const int n_powers = 5;
  // const int n_deltas = 5;
  const int n_betas = 8;

  const double infinity = std::numeric_limits<double>::max();

  // double optimal_power_values[n_betas] = {4.0, 2.0, 1.0, 0.67, 0.5, 0.33, 0.25, 0.1};
  // double optimal_delta_values[n_betas] = {100.0, 100.0, 100.0, 2.0, 1.0, 0.5, 0.33, 0.11};
  double beta_values[n_betas] = {0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 10.0};
  double optimal_njet_values_mean[n_betas];
  double optimal_njet_values_meansq[n_betas];
  double optimal_njet_values_rootmean[n_betas];
  vector<vector<double>> optimal_njet_values_median;

  double min_njet_values[n_betas];
  double min_njet_tau3_values[n_betas];
  for (int i_val = 0; i_val < n_betas; i_val++) {min_njet_values[i_val] = infinity;}
  for (int i_val = 0; i_val < n_betas; i_val++) {min_njet_tau3_values[i_val] = infinity;}
  double min_delta_values[n_betas];
  double min_delta_tau3_values[n_betas];
  double min_power_values[n_betas];
  double min_power_tau3_values[n_betas];
  double min_beta_values[n_betas];
  double min_beta_tau3_values[n_betas];

  double min_njet_values_meansq[n_betas];
  double min_njet_tau3_values_meansq[n_betas];
  for (int i_val = 0; i_val < n_betas; i_val++) {min_njet_values_meansq[i_val] = infinity;}
  for (int i_val = 0; i_val < n_betas; i_val++) {min_njet_tau3_values_meansq[i_val] = infinity;}
  double min_delta_values_meansq[n_betas];
  double min_delta_tau3_values_meansq[n_betas];
  double min_power_values_meansq[n_betas];
  double min_power_tau3_values_meansq[n_betas];
  double min_beta_values_meansq[n_betas];
  double min_beta_tau3_values_meansq[n_betas];

  double min_njet_values_rootmean[n_betas];
  double min_njet_tau3_values_rootmean[n_betas];
  for (int i_val = 0; i_val < n_betas; i_val++) {min_njet_values_rootmean[i_val] = infinity;}
  for (int i_val = 0; i_val < n_betas; i_val++) {min_njet_tau3_values_rootmean[i_val] = infinity;}
  double min_delta_values_rootmean[n_betas];
  double min_delta_tau3_values_rootmean[n_betas];
  double min_power_values_rootmean[n_betas];
  double min_power_tau3_values_rootmean[n_betas];
  double min_beta_values_rootmean[n_betas];
  double min_beta_tau3_values_rootmean[n_betas];

  double min_njet_values_median[n_betas];
  double min_njet_tau3_values_median[n_betas];
  for (int i_val = 0; i_val < n_betas; i_val++) {min_njet_values_median[i_val] = infinity;}
  for (int i_val = 0; i_val < n_betas; i_val++) {min_njet_tau3_values_median[i_val] = infinity;}
  double min_delta_values_median[n_betas];
  double min_delta_tau3_values_median[n_betas];
  double min_power_values_median[n_betas];
  double min_power_tau3_values_median[n_betas];
  double min_beta_values_median[n_betas];
  double min_beta_tau3_values_median[n_betas];


  for (int i_power = 0; i_power < (n_powers); i_power++) {
  // for (int i_power = 0; i_power < (n_betas); i_power++) {
    // power = optimal_power_values[i_power];
    power = (double)(i_power + 1)/8;

  for (int i_delta = 0; i_delta < (n_deltas); i_delta++) {
  // for (int i_delta = 0; i_delta < (n_betas); i_delta++) {
    // delta = optimal_delta_values[i_delta];
    delta = (double)(i_delta + 1)/8;
    if (i_delta == (n_deltas - 1)) delta = 100;

  int total_ttbar_jets = 0;
  int total_dijets_jets = 0;

  // TH1F *tau6_ttbar_manual_beta0p5 = new TH1F("tau6_ttbar_manual_beta0p5", "#tau_{6} for ttbar (Manual, #beta = 0.5)", 200, 0, 1000);
  // TH1F *tau6_dijets_manual_beta0p5 = new TH1F("tau6_dijets_manual_beta0p5", "#tau_{6} for dijets (Manual, #beta = 0.5)", 200, 0, 1000);

  // TH1F *tau6_ttbar_manual_beta1 = new TH1F("tau6_ttbar_manual_beta1", "#tau_{6} for ttbar (Manual, #beta = 1)", 200, 0, 1000);
  // TH1F *tau6_dijets_manual_beta1 = new TH1F("tau6_dijets_manual_beta1", "#tau_{6} for dijets (Manual, #beta = 1)", 200, 0, 1000);

  // TH1F *tau6_ttbar_manual_beta2 = new TH1F("tau6_ttbar_manual_beta2", "#tau_{6} for ttbar (Manual, #beta = 2)", 200, 0, 1000);
  // TH1F *tau6_dijets_manual_beta2 = new TH1F("tau6_dijets_manual_beta2", "#tau_{6} for dijets (Manual, #beta = 2)", 200, 0, 1000);

  // TH1F *tau6_ttbar_manual_beta3 = new TH1F("tau6_ttbar_manual_beta3", "#tau_{6} for ttbar (Manual, #beta = 3)", 200, 0, 1000);
  // TH1F *tau6_dijets_manual_beta3 = new TH1F("tau6_dijets_manual_beta3", "#tau_{6} for dijets (Manual, #beta = 3)", 200, 0, 1000);

  // TH1F *tau6_ttbar_manual_beta4 = new TH1F("tau6_ttbar_manual_beta4", "#tau_{6} for ttbar (Manual, #beta = 4)", 200, 0, 1000);
  // TH1F *tau6_dijets_manual_beta4 = new TH1F("tau6_dijets_manual_beta4", "#tau_{6} for dijets (Manual, #beta = 4)", 200, 0, 1000);

  cout << power << " " << delta << endl;

  vector<double> mean_tau6_ttbar_values;
  vector<double> mean_tau3_ttbar_values;

  vector<double> meansq_tau6_ttbar_values;
  vector<double> meansq_tau3_ttbar_values;

  vector<double> rootmean_tau6_ttbar_values;
  vector<double> rootmean_tau3_ttbar_values;

  vector<vector<double>> median_tau6_ttbar_values(n_betas);
  vector<vector<double>> median_tau3_ttbar_values(n_betas);

  for (int i = 0; i < n_betas; i++) {
    mean_tau6_ttbar_values.push_back(0.0);
    mean_tau3_ttbar_values.push_back(0.0);

    meansq_tau6_ttbar_values.push_back(0.0);
    meansq_tau3_ttbar_values.push_back(0.0);

    rootmean_tau6_ttbar_values.push_back(0.0);
    rootmean_tau3_ttbar_values.push_back(0.0);
  }

  TH1F *tau32_ttbar_wta = new TH1F("tau32_ttbar_wta", "#tau_{32} for ttbar (WTA)", 50, 0, 1);
  TH1F *tau32_dijets_wta = new TH1F("tau32_dijets_wta", "#tau_{32} for dijets (WTA)", 50, 0, 1);

  TH1F *tau32_ttbar_manual = new TH1F("tau32_ttbar_manual", "#tau_{32} for ttbar (manual)", 50, 0, 1);
  TH1F *tau32_dijets_manual = new TH1F("tau32_dijets_manual", "#tau_{32} for dijets (manual)", 50, 0, 1);

  TH1F *mass_ttbar_manual = new TH1F("mass_ttbar_manual", "mass ttbar (Manual)", 250, 0, 500);
  TH1F *mass_dijets_manual = new TH1F("mass_dijets_manual", "mass dijets (Manual)", 250, 0, 500);

  TH1F *mass_ttbar_antikt = new TH1F("mass_ttbar_antikt", "mass ttbar (antikt)", 250, 0, 500);
  TH1F *mass_dijets_antikt = new TH1F("mass_dijets_antikt", "mass dijets (antikt)", 250, 0, 500);

  std::string samples[2] = {"/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW", "/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW"};
  // std::string samples[2] = {"/home/tjwilk/Physics/thaler_urop/BOOST_samples/pythia64-tuneDW-lhc7-ttbar2hadrons-pt0500-0600.UW", "/home/tjwilk/Physics/thaler_urop/BOOST_samples/pythia64-tuneDW-lhc7-dijets-pt0500-0600.UW"};

  for (int i_sample = 0; i_sample < 1; i_sample++) {

    const char* current_data = samples[i_sample].c_str();
    ifstream inputStream(current_data); // Input File Name

    ////////// Set up Input //////////////////////////////////////////////////////////////////////////////
    vector<fastjet::PseudoJet> input_particles, partons;
    int particle_number, particle_code;
    double px, py, pz, E;
    int i_event = 1; bool is_first_event = true;
     
    ////////// Read in loop //////////////////////////////////////////////////////////////////////////////
    while (inputStream >> particle_number >> px >> py >> pz >> E >> particle_code){
       
      if (i_event > n_event) {break;}
      if (is_first_event) {is_first_event = false; continue;}
      if (particle_number == 0) {

        TH2F *event_display;

        if (i_sample == 0) event_display = new TH2F("event_display", "TTbar event (500 < pT < 600)", 60, -5, 5, 60, 0, 6.4);
        if (i_sample == 1) event_display = new TH2F("event_display", "Dijets event (500 < pT < 600)", 60, -5, 5, 60, 0, 6.4);
        TH2F *axes_manual_display = new TH2F("axes_manual_display", "Axes Plot", 300, -5, 5, 300, 0, 6.4);
        TH2F *axes_genrecomb_display = new TH2F("axes_genrecomb_display", "Axes Plot", 300, -5, 5, 300, 0, 6.4);
        TH2F *subjet1_display = new TH2F("subjet1_display", "Subjet display", 1000, -5, 5, 1000, 0, 6.4);
        TH2F *subjet2_display = new TH2F("subjet2_display", "Subjet display", 1000, -5, 5, 1000, 0, 6.4);
        TH2F *subjet3_display = new TH2F("subjet3_display", "Subjet display", 1000, -5, 5, 1000, 0, 6.4);
        TH2F *subjet4_display = new TH2F("subjet4_display", "Subjet display", 1000, -5, 5, 1000, 0, 6.4);
        TH2F *subjet5_display = new TH2F("subjet5_display", "Subjet display", 1000, -5, 5, 1000, 0, 6.4);
        TH2F *subjet6_display = new TH2F("subjet6_display", "Subjet display", 1000, -5, 5, 1000, 0, 6.4);
        TH2F *subjet7_display = new TH2F("subjet7_display", "Subjet display", 1000, -5, 5, 1000, 0, 6.4);

        vector<PseudoJet> jet_constituents = input_particles;

        for (int j = 0; j < jet_constituents.size(); j++) {
          event_display->Fill(jet_constituents[j].eta(), jet_constituents[j].phi(), jet_constituents[j].perp()*100);
        }

        // Double_t ghost_perp = 0.001;
        // Double_t n_ghosts = 50;
        // for (int i_ghosts = 0; i_ghosts < n_ghosts; i_ghosts++) {
        //   for (int j_ghosts = 0; j_ghosts < n_ghosts; j_ghosts++) {
        //     Double_t ghost_eta = -5 + (double)(10/n_ghosts)*i_ghosts;
        //     Double_t ghost_phi = (double)(2*TMath::Pi()/n_ghosts)*j_ghosts;
        //     PseudoJet ghost(ghost_perp*TMath::Cos(ghost_phi),ghost_perp*TMath::Sin(ghost_phi),
        //        ghost_perp*TMath::SinH(ghost_eta),ghost_perp*TMath::CosH(ghost_eta));
        //     jet_constituents.push_back(ghost);
        //   }
        // } 


        // Manual Axis Finding for general KT and general Recombiner
        // double Rparam = JetDefinition::max_allowable_R;
        double Rparam = Rcutoff;
        Strategy strategy = Best;
        const JetDefinition::Recombiner *recombScheme = new GeneralERecombiner(delta);
        // const JetDefinition::Recombiner *recombScheme = new WinnerTakeAllRecombiner();
        JetDefinition *jetDef = new JetDefinition(genkt_algorithm, Rparam, power, recombScheme, strategy);
        ClusterSequence clustSeq(jet_constituents, *jetDef);
        // vector<PseudoJet> exclusiveJets = clustSeq.exclusive_jets(njettiness);
        vector<PseudoJet> exclusiveJets = clustSeq.exclusive_jets(initial_axes_num);

        double beta;
        for (int i_beta = 0; i_beta < n_betas; i_beta++) {
          beta = beta_values[i_beta];

          UnnormalizedCutoffMeasure measure_function = UnnormalizedCutoffMeasure(beta, Rcutoff);
          // GeometricCutoffMeasure measure_function = GeometricCutoffMeasure(beta, Rcutoff);

          // NjettinessPlugin njet_plugin_genrecomb(initial_axes_num, Manual_Axes(), measure_function);
          // JetDefinition njet_def_genrecomb(&njet_plugin_genrecomb);
          // njet_plugin_genrecomb.setAxes(exclusiveJets);
          // ClusterSequence njet_cluster_genrecomb(jet_constituents, njet_def_genrecomb);
          // const NjettinessExtras *extras_genrecomb = njettiness_extras(njet_cluster_genrecomb);
          // vector<PseudoJet> axes_genrecomb = extras_genrecomb->axes();
          // vector<PseudoJet> jets_genrecomb = extras_genrecomb->jets();

          // JetInformation min_manual_set = findMinAxes(jet_constituents, axes_genrecomb, njettiness, beta, Rcutoff);
          // vector<PseudoJet> min_manual_axes = min_manual_set.axes;
          // double min_manual_tau = min_manual_set.tau;

          NjettinessPlugin njet_plugin_min_manual(njettiness, Manual_Axes(), measure_function);
          JetDefinition njet_def_min_manual(&njet_plugin_min_manual);
          // njet_plugin_min_manual.setAxes(min_manual_axes);
          njet_plugin_min_manual.setAxes(exclusiveJets);
          ClusterSequence njet_cluster_min_manual(jet_constituents, njet_def_min_manual);
          const NjettinessExtras *extras_min_manual = njettiness_extras(njet_cluster_min_manual);
          vector<PseudoJet> min_manual_jets = extras_min_manual->jets();
          double min_manual_tau = extras_min_manual->totalTau();

          if (i_sample == 0) {
            // cout << min_manual_tau << endl;

            mean_tau6_ttbar_values[i_beta] += min_manual_tau;
            meansq_tau6_ttbar_values[i_beta] += min_manual_tau*min_manual_tau;
            rootmean_tau6_ttbar_values[i_beta] += TMath::Sqrt(min_manual_tau);
            median_tau6_ttbar_values[i_beta].push_back(min_manual_tau);
          }
        
          // for (int i = 0; i < njettiness; i++) {
          //   // vector<PseudoJet> constituents = jets_genrecomb[i].constituents();
          //   vector<PseudoJet> constituents = min_manual_jets[i].constituents();
          //   // vector<PseudoJet> constituents = jets_wta[i].constituents();
          //   // vector<PseudoJet> constituents = jets_onepass[i].constituents();
          //   for (int i_const = 0; i_const < constituents.size(); i_const++) {
          //     if (i == 0) subjet1_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
          //     if (i == 1) subjet2_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
          //     if (i == 2) subjet3_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
          //     if (i == 3) subjet4_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
          //     if (i == 4) subjet5_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
          //     if (i == 5) subjet6_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
          //     if (i == 6) subjet7_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
          //   }          
          // }

          //LOOK AT INDIVIDUAL ANTI-KT JETS
          vector <PseudoJet> inclusiveJets_std, sortedJets_std, hardestJets_std;

          double Rparam_std = 1.0;
          Strategy strategy_std = Best;
          RecombinationScheme recombScheme_std = E_scheme;
          JetDefinition *jetDef_std = new JetDefinition(antikt_algorithm, Rparam_std, recombScheme_std, strategy_std);
          ClusterSequence clustSeq_std(input_particles, *jetDef_std);
          inclusiveJets_std = clustSeq_std.inclusive_jets(200.0);
          sortedJets_std = sorted_by_pt(inclusiveJets_std);
          Selector jet_selector = SelectorNHardest(2);
          hardestJets_std = jet_selector(sortedJets_std);

          for (int i = 0; i < hardestJets_std.size(); i++) {

            // if (i_sample == 1) {
            //   total_dijets_jets++;
            //   mass_dijets_antikt->Fill(hardestJets_std[i].m());
            // }

            ClusterSequence clustSeq_subjet(hardestJets_std[i].constituents(), *jetDef);
            vector<PseudoJet> starting_subjets = clustSeq_subjet.exclusive_jets(3);
            JetInformation subjet_set3 = findMinAxes(hardestJets_std[i].constituents(), starting_subjets, 3, beta, Rcutoff);
            JetInformation subjet_set2 = findMinAxes(hardestJets_std[i].constituents(), starting_subjets, 2, beta, Rcutoff);
            vector<PseudoJet> subjet_axes3 = subjet_set3.axes;
            vector<PseudoJet> subjet_axes2 = subjet_set2.axes;

            // NsubjettinessRatio nsub_32_wta(3, 2, WTA_KT_Axes(), UnnormalizedMeasure(1.0));

            if (i_sample == 0) {
              if (i_beta == 0) total_ttbar_jets++;
              // mass_ttbar_antikt->Fill(hardestJets_std[i].m());
              mean_tau3_ttbar_values[i_beta] += subjet_set3.tau;
              meansq_tau3_ttbar_values[i_beta] += pow(subjet_set3.tau, 2);
              rootmean_tau3_ttbar_values[i_beta] += pow(subjet_set3.tau, 0.5);
              median_tau3_ttbar_values[i_beta].push_back(subjet_set3.tau);
            }

            // if (hardestJets_std[i].m() > 160 && hardestJets_std[i].m() < 240) {
            //   if (i_sample == 0) {
            //     tau32_ttbar_manual->Fill(subjet_set3.tau/subjet_set2.tau);
            //     tau32_ttbar_wta->Fill(nsub_32_wta.result(hardestJets_std[i]));
            //   }
            //   if (i_sample == 1) {
            //     tau32_dijets_manual->Fill(subjet_set3.tau/subjet_set2.tau);
            //     tau32_dijets_wta->Fill(nsub_32_wta.result(hardestJets_std[i]));
            //   }
            // }
          }
        }

        event_display->SetStats(0);
        axes_genrecomb_display->SetStats(0);
        axes_manual_display->SetStats(0);
        subjet1_display->SetStats(0);
        subjet2_display->SetStats(0);
        subjet3_display->SetStats(0);
        subjet4_display->SetStats(0);
        subjet5_display->SetStats(0);
        subjet6_display->SetStats(0);
        subjet7_display->SetStats(0);

        TCanvas *display = new TCanvas("display", "Event Display", 1093, 700);
        display->cd();
        display->SetFixedAspectRatio();
        event_display->GetXaxis()->SetTitle("#eta");
        event_display->GetYaxis()->SetTitle("#phi");
        event_display->SetFillColor(kBlack);
        event_display->SetLineColor(kBlack);
        event_display->SetLineWidth(1);
        subjet1_display->SetMarkerColor(kRed);
        subjet2_display->SetMarkerColor(kBlue);
        subjet3_display->SetMarkerColor(kGreen);
        subjet4_display->SetMarkerColor(kYellow);
        subjet5_display->SetMarkerColor(kOrange);
        subjet6_display->SetMarkerColor(6);
        subjet7_display->SetMarkerColor(20);
        subjet1_display->SetMarkerStyle(21);
        subjet2_display->SetMarkerStyle(21);
        subjet3_display->SetMarkerStyle(21);
        subjet4_display->SetMarkerStyle(21);
        subjet5_display->SetMarkerStyle(21);
        subjet6_display->SetMarkerStyle(21);
        subjet7_display->SetMarkerStyle(21);
        subjet1_display->SetMarkerSize(0.5);
        subjet2_display->SetMarkerSize(0.5);
        subjet3_display->SetMarkerSize(0.5);
        subjet4_display->SetMarkerSize(0.5);
        subjet5_display->SetMarkerSize(0.5);
        subjet6_display->SetMarkerSize(0.5);
        subjet7_display->SetMarkerSize(0.5);

        event_display->Draw("box");

        // for (int a = 0; a < min_manual_axes.size(); a++) {
        //   axes_manual_display->Fill(min_manual_axes[a].eta(), min_manual_axes[a].phi());
        // }
        // axes_manual_display->SetMarkerStyle(3);
        // axes_manual_display->SetMarkerSize(3);
        // axes_manual_display->SetMarkerColor(kGreen);
        // axes_manual_display->Draw("SAMES");

        // for (int a = 0; a < axes_genrecomb.size(); a++) {
        //   axes_genrecomb_display->Fill(axes_genrecomb[a].eta(), axes_genrecomb[a].phi());
        // }
        // axes_genrecomb_display->SetMarkerStyle(3);
        // axes_genrecomb_display->SetMarkerSize(3);
        // axes_genrecomb_display->SetMarkerColor(kGreen);
        // axes_genrecomb_display->Draw("SAMES");

        // subjet1_display->Draw("SAMES");
        // subjet2_display->Draw("SAMES");
        // subjet3_display->Draw("SAMES");
        // subjet4_display->Draw("SAMES");
        // subjet5_display->Draw("SAMES");
        // subjet6_display->Draw("SAMES");
        // subjet7_display->Draw("SAMES");

        // ghost_wta_display->SetMarkerStyle(20);
        // ghost_wta_display->SetMarkerSize(1);
        // ghost_wta_display->Draw("SAMES");
        // ghost_onepass_display->Draw("SAMES");

        // display->Write();

        // ostringstream ss;
        // ss << i_event;
        // TString title;
        // if (i_sample == 0) title = "ttbar_display" + ss.str() + "_wtap8_beta0p125.pdf";
        // if (i_sample == 1) title = "dijets_display" + ss.str() + "_wtap8_beta0p125.pdf";

        // display->Print(title, "pdf");

        delete display;
        delete event_display;
        delete axes_genrecomb_display;
        delete axes_manual_display;
        delete subjet1_display;
        delete subjet2_display;
        delete subjet3_display;
        delete subjet4_display;
        delete subjet5_display;
        delete subjet6_display;
        delete subjet7_display;

        // Reset pseudojets for next event
        i_event = particle_code;
        input_particles.resize(0);
        partons.resize(0);
      }
      else if (particle_number == -1){
        partons.push_back(fastjet::PseudoJet(px,py,pz,E));
      }
      else if (abs(particle_code) == 12 || abs(particle_code) == 13 || abs(particle_code) == 14 || abs(particle_code) == 16) {
        // Do nothing
      }
      else {
        fastjet::PseudoJet temp = fastjet::PseudoJet(px,py,pz,E);
        if (fabs(temp.eta()) < 5) {
          input_particles.push_back(temp);
        }
      }
    }
  }

  // tau6_ttbar_manual->SetStats(0);
  // tau6_dijets_manual->SetStats(0);
  mass_ttbar_manual->SetStats(0);
  mass_dijets_manual->SetStats(0);
  tau32_ttbar_manual->SetStats(0);
  tau32_dijets_manual->SetStats(0);
  tau32_ttbar_wta->SetStats(0);
  tau32_dijets_wta->SetStats(0);

  // tau6_ttbar_manual->Write();
  // tau6_dijets_manual->Write();
  // mass_ttbar_manual->Write();
  // mass_dijets_manual->Write();
  // tau32_ttbar_manual->Write();
  // tau32_dijets_manual->Write();
  // tau32_ttbar_wta->Write();
  // tau32_dijets_wta->Write();

  // TCanvas *tau6_dijets_compare = new TCanvas("tau6_dijets_compare", "Mass Comparison", 600, 600);
  // tau6_dijets_compare->cd();
  // tau6_dijets_wta->SetLineColor(kRed);
  // tau6_dijets_onepass->SetLineColor(kBlue);
  // tau6_dijets_wta->Draw();
  // tau6_dijets_onepass->Draw("same");
  // TLegend *leg_tau6_dijets = new TLegend(0.1, 0.7, 0.3, 0.9);
  // leg_tau6_dijets->AddEntry(tau6_dijets_wta, "WTA", "L");
  // leg_tau6_dijets->AddEntry(tau6_dijets_onepass, "OnePass WTA", "L");
  // leg_tau6_dijets->Draw();
  // tau6_dijets_compare->Write();

  // int n_wtakt = tau32_ttbar_wta->GetSize() - 2;
  // double integral_dijets_wtakt[n_wtakt], integral_ttbar_wtakt[n_wtakt];
  // for (int i = 0; i < n_wtakt; i++) {
  //   integral_ttbar_wtakt[i] = tau32_ttbar_wta->Integral(0,i)/total_ttbar_jets;
  //   integral_dijets_wtakt[i] = tau32_dijets_wta->Integral(0,i)/total_dijets_jets;
  // }
  // TGraph* ROC_tau32_std = new TGraph(n_wtakt, integral_ttbar_wtakt, integral_dijets_wtakt);
  // ROC_tau32_std->GetXaxis()->SetLimits(0, 1);
  // ROC_tau32_std->GetYaxis()->SetLimits(0, 1);
  // ROC_tau32_std->SetLineColor(kBlack);
  // ROC_tau32_std->SetLineWidth(2);
  // ROC_tau32_std->SetMarkerStyle(5);
  // ROC_tau32_std->SetMarkerSize(2);
  // ROC_tau32_std->Write();

  // int n_size = tau32_ttbar_manual->GetSize() - 2;
  // double integral_dijets[n_size], integral_ttbar[n_size];
  // for (int i = 0; i < n_size; i++) {
  //   integral_ttbar[i] = tau32_ttbar_manual->Integral(0,i)/total_ttbar_jets;
  //   integral_dijets[i] = tau32_dijets_manual->Integral(0,i)/total_dijets_jets;
  // }
  // TGraph* ROC_tau32 = new TGraph(n_size, integral_ttbar, integral_dijets);
  // ROC_tau32->GetXaxis()->SetLimits(0, 1);
  // ROC_tau32->GetYaxis()->SetLimits(0, 1);
  // ROC_tau32->SetLineColor(kRed);
  // ROC_tau32->SetLineWidth(2);
  // ROC_tau32->SetMarkerStyle(5);
  // ROC_tau32->SetMarkerSize(2);
  // ROC_tau32->Write();

  // TCanvas *ROC_compare = new TCanvas("ROC_compare", "ROC Curve Comparison", 600, 600);
  // ROC_compare->cd();
  // ROC_compare->SetLogy();
  // TMultiGraph *ROC_multigraph = new TMultiGraph("ROC_multigraph", "ROC Comparison for #tau_{3}/#tau_{2} (500 < p_{T} < 600)");
  // ROC_multigraph->Add(ROC_tau32_std);
  // ROC_multigraph->Add(ROC_tau32);
  // ROC_multigraph->Draw("AL");
  // TLegend *leg_ROC = new TLegend(0.1, 0.7, 0.3, 0.9);
  // leg_ROC->AddEntry(ROC_tau32_std, "WTA kT", "L");
  // leg_ROC->AddEntry(ROC_tau32, "Manual Axes", "L");
  // leg_ROC->Draw();
  // ROC_compare->Write();

  for (int i_beta = 0; i_beta < n_betas; i_beta++) {

    double mean_value = (double)mean_tau6_ttbar_values[i_beta]/n_event;
    double mean_tau3_value = (double)mean_tau3_ttbar_values[i_beta]/total_ttbar_jets;
    if (mean_value < min_njet_values[i_beta]) {
      min_njet_values[i_beta] = mean_value;
      min_beta_values[i_beta] = beta_values[i_beta];
      min_delta_values[i_beta] = delta;
      min_power_values[i_beta] = power;
    }
    if (mean_tau3_value < min_njet_tau3_values[i_beta]) {
      min_njet_tau3_values[i_beta] = mean_tau3_value;
      min_beta_tau3_values[i_beta] = beta_values[i_beta];
      min_delta_tau3_values[i_beta] = delta;
      min_power_tau3_values[i_beta] = power;
    }

    double meansq_value = (double)meansq_tau6_ttbar_values[i_beta]/n_event;
    double meansq_tau3_value = (double)meansq_tau3_ttbar_values[i_beta]/total_ttbar_jets;
    if (meansq_value < min_njet_values_meansq[i_beta]) {
      min_njet_values_meansq[i_beta] = meansq_value;
      min_beta_values_meansq[i_beta] = beta_values[i_beta];
      min_delta_values_meansq[i_beta] = delta;
      min_power_values_meansq[i_beta] = power;
    }
    if (meansq_tau3_value < min_njet_tau3_values_meansq[i_beta]) {
      min_njet_tau3_values_meansq[i_beta] = meansq_tau3_value;
      min_beta_tau3_values_meansq[i_beta] = beta_values[i_beta];
      min_delta_tau3_values_meansq[i_beta] = delta;
      min_power_tau3_values_meansq[i_beta] = power;
    }

    double rootmean_value = (double)rootmean_tau6_ttbar_values[i_beta]/n_event;
    double rootmean_tau3_value = (double)rootmean_tau3_ttbar_values[i_beta]/total_ttbar_jets;
    if (rootmean_value < min_njet_values_rootmean[i_beta]) {
      min_njet_values_rootmean[i_beta] = rootmean_value;
      min_beta_values_rootmean[i_beta] = beta_values[i_beta];
      min_delta_values_rootmean[i_beta] = delta;
      min_power_values_rootmean[i_beta] = power;
    }
    if (rootmean_tau3_value < min_njet_tau3_values_rootmean[i_beta]) {
      min_njet_tau3_values_rootmean[i_beta] = rootmean_tau3_value;
      min_beta_tau3_values_rootmean[i_beta] = beta_values[i_beta];
      min_delta_tau3_values_rootmean[i_beta] = delta;
      min_power_tau3_values_rootmean[i_beta] = power;
    }

    double median_value = calcMedian(median_tau6_ttbar_values[i_beta]);
    double median_tau3_value = calcMedian(median_tau3_ttbar_values[i_beta]);
    if (median_value < min_njet_values_median[i_beta]) {
      min_njet_values_median[i_beta] = median_value;
      min_beta_values_median[i_beta] = beta_values[i_beta];
      min_delta_values_median[i_beta] = delta;
      min_power_values_median[i_beta] = power;
    }
    if (median_tau3_value < min_njet_tau3_values_median[i_beta]) {
      min_njet_tau3_values_median[i_beta] = median_tau3_value;
      min_beta_tau3_values_median[i_beta] = beta_values[i_beta];
      min_delta_tau3_values_median[i_beta] = delta;
      min_power_tau3_values_median[i_beta] = power;
    }

  }

  // cout << "mean manual tau6: " << tau6_ttbar_manual->GetMean() << endl;
  // cout << "sqrt(mean^2) manual tau6: " << TMath::Sqrt(TMath::Power(tau6_ttbar_manual->GetRMS(),2) + TMath::Power(tau6_ttbar_manual->GetMean(),2)) << endl;
  // cout << endl;

  // mean_values.push_back(tau6_ttbar_manual->GetMean());
  // rms_values.push_back(TMath::Sqrt(TMath::Power(tau6_ttbar_manual->GetRMS(),2) + TMath::Power(tau6_ttbar_manual->GetMean(),2)));

  // double tau6_ttbar_manual_scale = 1/tau6_ttbar_manual->Integral();
  // double tau6_dijets_manual_scale = 1/tau6_dijets_manual->Integral();
  double mass_ttbar_manual_scale = 1/mass_ttbar_manual->Integral();
  double mass_dijets_manual_scale = 1/mass_dijets_manual->Integral();

  // tau6_ttbar_manual->Scale(tau6_ttbar_manual_scale);
  // tau6_dijets_manual->Scale(tau6_dijets_manual_scale);
  mass_ttbar_manual->Scale(mass_ttbar_manual_scale);
  mass_dijets_manual->Scale(mass_dijets_manual_scale);

  // TCanvas *tau6_ttbar_compare = new TCanvas("tau6_ttbar_compare", "Tau32 Comparison", 600, 600);
  // tau6_ttbar_compare->cd();
  // tau6_ttbar_wta->SetLineColor(kRed);
  // tau6_ttbar_manual->SetLineColor(kBlue);
  // tau6_ttbar_wta->Draw();
  // tau6_ttbar_manual->Draw("same");
  // TLegend *leg_tau6_ttbar = new TLegend(0.1, 0.7, 0.3, 0.9);
  // leg_tau6_ttbar->AddEntry(tau6_ttbar_wta, "WTA", "L");
  // leg_tau6_ttbar->AddEntry(tau6_ttbar_manual, "Manual", "L");
  // leg_tau6_ttbar->Draw();
  // tau6_ttbar_compare->Write();

  // TCanvas *mass_ttbar_compare = new TCanvas("mass_ttbar_compare", "Mass Comparison", 600, 600);
  // mass_ttbar_compare->cd();
  // mass_ttbar_wta->SetLineColor(kRed);
  // mass_ttbar_manual->SetLineColor(kBlue);
  // mass_ttbar_wta->Draw();
  // mass_ttbar_manual->Draw("same");
  // TLegend *leg_mass_ttbar = new TLegend(0.1, 0.7, 0.3, 0.9);
  // leg_mass_ttbar->AddEntry(mass_ttbar_wta, "WTA", "L");
  // leg_mass_ttbar->AddEntry(mass_ttbar_manual, "Manual", "L");
  // leg_mass_ttbar->Draw();
  // mass_ttbar_compare->Write();

  // delete tau6_ttbar_manual;
  // delete tau6_dijets_manual;
  delete tau32_ttbar_wta;
  delete tau32_ttbar_manual;
  delete tau32_dijets_manual;
  delete tau32_dijets_wta;
  delete mass_ttbar_manual;
  delete mass_dijets_manual;
  delete mass_ttbar_antikt;
  delete mass_dijets_antikt;
  // delete ROC_compare;
  // delete mass_ttbar_compare;
  // delete tau6_ttbar_compare;

// }
// }
// }
}
}  


  // TH2F *min_beta_delta = new TH2F("min_beta_delta", "#beta vs #delta (#tau_{6})", );
  // TH2F *min_beta_delta = new TH2F("min_beta_delta", "#beta vs #delta (#tau_{6})");

  double optimal_njet_tau6_mean[n_betas] = {397.256, 248.407, 130.238, 81.5934, 55.476, 28.8354, 16.3192, 0.689889};
  double optimal_njet_tau6_meansq[n_betas] = {169827, 68395.2, 20011.8, 8210.13, 3906.22, 1100.11, 361.704, 0.679387};
  double optimal_njet_tau6_rootmean[n_betas] = {19.6949, 15.5159, 11.1409, 8.76259, 7.19512, 5.1565, 3.86299, 0.786335};
  double optimal_njet_tau6_median[n_betas] = {391.352, 240.156, 123.458, 76.5458, 51.413, 26.7442, 14.7982, 0.615171};

  double percent_diff_mean[n_betas];
  double percent_diff_meansq[n_betas];
  double percent_diff_rootmean[n_betas];
  double percent_diff_median[n_betas];

  for (int i = 0; i < n_betas; i++) {
    percent_diff_mean[i] = -(min_njet_values[i] - optimal_njet_tau6_mean[i])/min_njet_values[i];
    percent_diff_meansq[i] = -(min_njet_values_meansq[i] - optimal_njet_tau6_meansq[i])/min_njet_values_meansq[i];
    percent_diff_rootmean[i] = -(min_njet_values_rootmean[i] - optimal_njet_tau6_rootmean[i])/min_njet_values_rootmean[i];
    percent_diff_median[i] = -(min_njet_values_median[i] - optimal_njet_tau6_median[i])/min_njet_values_median[i];
  }

  TGraph *percent_diff_mean_plot = new TGraph(n_betas, beta_values, percent_diff_mean);
  TGraph *percent_diff_meansq_plot = new TGraph(n_betas, beta_values, percent_diff_meansq);
  TGraph *percent_diff_rootmean_plot = new TGraph(n_betas, beta_values, percent_diff_rootmean);
  TGraph *percent_diff_median_plot = new TGraph(n_betas, beta_values, percent_diff_median);

  TMultiGraph *percent_diff_hists = new TMultiGraph();
  TCanvas *percent_diff_hists_can = new TCanvas("percent_diff_hists_can", "percent_diff_hists", 600, 600);
  percent_diff_hists_can->cd();
  // percent_diff_mean_plot->SetMarkerStyle(3);
  percent_diff_mean_plot->SetMarkerSize(2);
  // percent_diff_mean_plot->SetMarkerColor(kBlack);
  percent_diff_mean_plot->SetLineColor(kRed);
  percent_diff_mean_plot->SetLineWidth(3);
  // percent_diff_meansq_plot->SetMarkerStyle(3);
  percent_diff_meansq_plot->SetMarkerSize(2);
  // percent_diff_meansq_plot->SetMarkerColor(kBlue);
  percent_diff_meansq_plot->SetLineColor(kOrange);
  percent_diff_meansq_plot->SetLineWidth(3);
  // percent_diff_rootmean_plot->SetMarkerStyle(3);
  percent_diff_rootmean_plot->SetMarkerSize(2);
  // percent_diff_rootmean_plot->SetMarkerColor(kGreen);
  percent_diff_rootmean_plot->SetLineColor(kGreen);
  percent_diff_rootmean_plot->SetLineWidth(3);
  // percent_diff_median_plot->SetMarkerStyle(3);
  percent_diff_median_plot->SetMarkerSize(2);
  // percent_diff_median_plot->SetMarkerColor(kRed);
  percent_diff_median_plot->SetLineColor(kBlue);
  percent_diff_median_plot->SetLineWidth(3);
  // percent_diff_mean_plot->GetHistogram()->SetMinimum(-1.0);
  // percent_diff_mean_plot->GetHistogram()->SetMaximum(1.0);
  // percent_diff_meansq_plot->GetHistogram()->SetMinimum(-1.0);
  // percent_diff_meansq_plot->GetHistogram()->SetMaximum(1.0);
  // percent_diff_rootmean_plot->GetHistogram()->SetMinimum(-1.0);
  // percent_diff_rootmean_plot->GetHistogram()->SetMaximum(1.0);
  // percent_diff_median_plot->GetHistogram()->SetMinimum(-1.0);
  // percent_diff_median_plot->GetHistogram()->SetMaximum(1.0);
  percent_diff_hists->Add(percent_diff_mean_plot);
  percent_diff_hists->Add(percent_diff_median_plot);
  // percent_diff_hists->Add(percent_diff_meansq_plot);
  // percent_diff_hists->Add(percent_diff_rootmean_plot);
  percent_diff_hists->SetMinimum(-0.5);
  percent_diff_hists->SetMaximum(0.5);
  percent_diff_hists->Draw("AL");
  percent_diff_hists->SetTitle("% Difference Between Minimum and ``Optimal'' #tau_{6}");
  percent_diff_hists->GetXaxis()->SetTitle("#beta");
  percent_diff_hists->GetYaxis()->SetTitle("% difference");
  TLegend *leg_percent_diff = new TLegend(0.14, 0.7, 0.4, 0.88);
  leg_percent_diff->SetFillColor(kWhite);
  leg_percent_diff->SetLineColor(kWhite);
  leg_percent_diff->AddEntry(percent_diff_mean_plot, "Mean", "L");
  // leg_percent_diff->AddEntry(percent_diff_meansq_plot, "Mean^{2}", "L");
  // leg_percent_diff->AddEntry(percent_diff_rootmean_plot, "#sqrt{Mean}", "L");
  leg_percent_diff->AddEntry(percent_diff_median_plot, "Median", "L");
  leg_percent_diff->Draw();
  percent_diff_hists_can->Write();
  percent_diff_hists_can->Print("percentdiff_tau6.eps", "eps");


  TF1 *beta_delta_func = new TF1("beta_delta_func", "(x <= 1)*100 + (x > 1)*1/(x-1)", 0, 10);
  TF1 *beta_power_func = new TF1("beta_power_func", "1/x", 0, 10);
  beta_delta_func->SetLineColor(kBlack);
  beta_power_func->SetLineColor(kBlack);

  // TGraph *min_beta_delta = new TGraph(n_betas, min_beta_values, min_delta_values);
  // TGraph *min_beta_power = new TGraph(n_betas, min_beta_values, min_power_values);
  // TGraph *min_beta_delta_meansq = new TGraph(n_betas, min_beta_values_meansq, min_delta_values_meansq);
  // TGraph *min_beta_power_meansq = new TGraph(n_betas, min_beta_values_meansq, min_power_values_meansq);
  // TGraph *min_beta_delta_rootmean = new TGraph(n_betas, min_beta_values_rootmean, min_delta_values_rootmean);
  // TGraph *min_beta_power_rootmean = new TGraph(n_betas, min_beta_values_rootmean, min_power_values_rootmean);
  // TGraph *min_beta_delta_median = new TGraph(n_betas, min_beta_values_median, min_delta_values_median);
  // TGraph *min_beta_power_median = new TGraph(n_betas, min_beta_values_median, min_power_values_median);

  TGraph *min_beta_delta = new TGraph(n_betas, beta_values, min_delta_values);
  TGraph *min_beta_power = new TGraph(n_betas, beta_values, min_power_values);
  TGraph *min_beta_delta_meansq = new TGraph(n_betas, beta_values, min_delta_values_meansq);
  TGraph *min_beta_power_meansq = new TGraph(n_betas, beta_values, min_power_values_meansq);
  TGraph *min_beta_delta_rootmean = new TGraph(n_betas, beta_values, min_delta_values_rootmean);
  TGraph *min_beta_power_rootmean = new TGraph(n_betas, beta_values, min_power_values_rootmean);
  TGraph *min_beta_delta_median = new TGraph(n_betas, beta_values, min_delta_values_median);
  TGraph *min_beta_power_median = new TGraph(n_betas, beta_values, min_power_values_median);

  TMultiGraph *min_tau6_values_delta = new TMultiGraph(); 
  TMultiGraph *min_tau6_values_power = new TMultiGraph(); 
  min_tau6_values_delta->SetTitle("#beta and #delta for Minimum #tau_{6}");
  min_tau6_values_power->SetTitle("#beta and p for Minimum #tau_{6}");
  // min_beta_delta->SetMarkerStyle(3);
  // min_beta_power->SetMarkerStyle(3);
  min_beta_delta->SetMarkerSize(2);
  min_beta_power->SetMarkerSize(2);
  min_beta_delta->SetLineColor(kRed);
  min_beta_power->SetLineColor(kRed);
  min_beta_delta->SetLineWidth(3);
  min_beta_power->SetLineWidth(3);


  // min_beta_delta_median->SetMarkerStyle(3);
  // min_beta_power_median->SetMarkerStyle(3);
  min_beta_delta_median->SetMarkerSize(2);
  min_beta_power_median->SetMarkerSize(2);
  // min_beta_delta_median->SetMarkerColor(kRed);
  // min_beta_power_median->SetMarkerColor(kRed);
  min_beta_delta_median->SetLineColor(kBlue);
  min_beta_power_median->SetLineColor(kBlue);
  min_beta_delta_median->SetLineWidth(3);
  min_beta_power_median->SetLineWidth(3);

  // min_beta_delta_meansq->SetMarkerStyle(3);
  // min_beta_power_meansq->SetMarkerStyle(3);
  min_beta_delta_meansq->SetMarkerSize(2);
  min_beta_power_meansq->SetMarkerSize(2);
  // min_beta_delta_meansq->SetMarkerColor(kBlue);
  // min_beta_power_meansq->SetMarkerColor(kBlue);
  min_beta_delta_meansq->SetLineColor(kOrange);
  min_beta_power_meansq->SetLineColor(kOrange);
  min_beta_delta_meansq->SetLineWidth(3);
  min_beta_power_meansq->SetLineWidth(3);
  // min_beta_delta_rootmean->SetMarkerStyle(3);
  // min_beta_power_rootmean->SetMarkerStyle(3);
  min_beta_delta_rootmean->SetMarkerSize(2);
  min_beta_power_rootmean->SetMarkerSize(2);
  // min_beta_delta_rootmean->SetMarkerColor(kGreen);
  // min_beta_power_rootmean->SetMarkerColor(kGreen);
  min_beta_delta_rootmean->SetLineColor(kGreen);
  min_beta_power_rootmean->SetLineColor(kGreen);
  min_beta_delta_rootmean->SetLineWidth(3);
  min_beta_power_rootmean->SetLineWidth(3);
  // min_beta_delta->Write();
  // min_beta_power->Write();

  TCanvas* beta_delta_can = new TCanvas("beta_delta_can", "beta_delta_can", 600, 600);
  beta_delta_can->cd();
  min_tau6_values_delta->Add(min_beta_delta);
  min_tau6_values_delta->Add(min_beta_delta_median);
  // min_tau6_values_delta->Add(min_beta_delta_meansq);
  // min_tau6_values_delta->Add(min_beta_delta_rootmean);
  min_tau6_values_delta->Draw("AL");
  beta_delta_func->Draw("SAME");
  min_tau6_values_delta->GetXaxis()->SetTitle("#beta");
  min_tau6_values_delta->GetYaxis()->SetTitle("#delta");
  min_tau6_values_delta->SetMaximum(10);
  TLegend *leg_beta_delta = new TLegend(0.55, 0.55, 0.88, 0.88);
  leg_beta_delta->SetFillColor(kWhite);
  leg_beta_delta->SetLineColor(kWhite);
  leg_beta_delta->AddEntry(min_beta_delta, "Mean", "L");
  leg_beta_delta->AddEntry(min_beta_delta_median, "Median", "L");
  // leg_beta_delta->AddEntry(min_beta_delta_meansq, "Mean^{2}", "L");
  // leg_beta_delta->AddEntry(min_beta_delta_rootmean, "#sqrt{Mean}", "L");
  leg_beta_delta->AddEntry(beta_delta_func, "Prediction", "L");
  leg_beta_delta->Draw();
  beta_delta_can->Write();
  beta_delta_can->Print("betavdelta_tau6.eps", "eps");

  TCanvas* beta_power_can = new TCanvas("beta_power_can", "beta_power_can", 600, 600);
  beta_power_can->cd();
  min_tau6_values_power->Add(min_beta_power);
  min_tau6_values_power->Add(min_beta_power_median);
  // min_tau6_values_power->Add(min_beta_power_meansq);
  // min_tau6_values_power->Add(min_beta_power_rootmean);
  min_tau6_values_power->Draw("AL");
  beta_power_func->Draw("SAME");
  min_tau6_values_power->GetXaxis()->SetTitle("#beta");
  min_tau6_values_power->GetYaxis()->SetTitle("power");
  TLegend *leg_beta_power = new TLegend(0.55, 0.55, 0.88, 0.88);
  leg_beta_power->SetFillColor(kWhite);
  leg_beta_power->SetLineColor(kWhite);
  leg_beta_power->AddEntry(min_beta_power, "Mean", "L");
  leg_beta_power->AddEntry(min_beta_power_median, "Median", "L");
  // leg_beta_power->AddEntry(min_beta_power_meansq, "Mean^{2}", "L");
  // leg_beta_power->AddEntry(min_beta_power_rootmean, "#sqrt{Mean}", "L");
  leg_beta_power->AddEntry(beta_power_func, "Prediction", "L");  
  leg_beta_power->Draw();
  beta_power_can->Write();
  beta_power_can->Print("betavpower_tau6.eps", "eps");

  double optimal_njet_tau3_mean[n_betas] = {135.389, 71.2266, 26.2456, 11.924, 6.28232, 2.31148, 1.05438, 0.0302973};
  double optimal_njet_tau3_meansq[n_betas] = {20401.3, 5864.69, 896.665, 206.831, 63.0371, 10.1063, 2.34966, 0.00333089};
  double optimal_njet_tau3_rootmean[n_betas] = {11.4668, 8.2811, 4.95462, 3.29348, 2.36186, 1.40286, 0.929607, 0.142458};
  double optimal_njet_tau3_median[n_betas] = {131.956, 67.2988, 22.952, 9.82452, 4.99128, 1.71188, 0.732627, 0.0152763};

  double percent_diff_tau3_mean[n_betas];
  double percent_diff_tau3_meansq[n_betas];
  double percent_diff_tau3_rootmean[n_betas];
  double percent_diff_tau3_median[n_betas];

  for (int i = 0; i < n_betas; i++) {
    percent_diff_tau3_mean[i] = -(min_njet_tau3_values[i] - optimal_njet_tau3_mean[i])/min_njet_tau3_values[i];
    percent_diff_tau3_meansq[i] = -(min_njet_tau3_values_meansq[i] - optimal_njet_tau3_meansq[i])/min_njet_tau3_values_meansq[i];
    percent_diff_tau3_rootmean[i] = -(min_njet_tau3_values_rootmean[i] - optimal_njet_tau3_rootmean[i])/min_njet_tau3_values_rootmean[i];
    percent_diff_tau3_median[i] = -(min_njet_tau3_values_median[i] - optimal_njet_tau3_median[i])/min_njet_tau3_values_median[i];
  }

  TGraph *percent_diff_tau3_mean_plot = new TGraph(n_betas, beta_values, percent_diff_tau3_mean);
  TGraph *percent_diff_tau3_meansq_plot = new TGraph(n_betas, beta_values, percent_diff_tau3_meansq);
  TGraph *percent_diff_tau3_rootmean_plot = new TGraph(n_betas, beta_values, percent_diff_tau3_rootmean);
  TGraph *percent_diff_tau3_median_plot = new TGraph(n_betas, beta_values, percent_diff_tau3_median);

  TMultiGraph *percent_diff_tau3_hists = new TMultiGraph();
  TCanvas *percent_diff_tau3_hists_can = new TCanvas("percent_diff_tau3_hists_can", "percent_diff_tau3_hists", 600, 600);
  percent_diff_tau3_hists_can->cd();
  // percent_diff_tau3_mean_plot->SetMarkerStyle(3);
  percent_diff_tau3_mean_plot->SetMarkerSize(2);
  // percent_diff_tau3_mean_plot->SetMarkerColor(kBlack);
  percent_diff_tau3_mean_plot->SetLineColor(kRed);
  percent_diff_tau3_mean_plot->SetLineWidth(3);
  percent_diff_tau3_median_plot->SetMarkerSize(2);
  // percent_diff_tau3_median_plot->SetMarkerColor(kRed);
  // percent_diff_tau3_median_plot->SetLineColor(kOrange);
  percent_diff_tau3_median_plot->SetLineColor(kBlue);
  percent_diff_tau3_median_plot->SetLineWidth(3);
  // percent_diff_tau3_meansq_plot->SetMarkerStyle(3);
  percent_diff_tau3_meansq_plot->SetMarkerSize(2);
  // percent_diff_tau3_meansq_plot->SetMarkerColor(kBlue);
  percent_diff_tau3_meansq_plot->SetLineColor(kOrange);
  percent_diff_tau3_meansq_plot->SetLineWidth(3);
  // percent_diff_tau3_rootmean_plot->SetMarkerStyle(3);
  percent_diff_tau3_rootmean_plot->SetMarkerSize(2);
  // percent_diff_tau3_rootmean_plot->SetMarkerColor(kGreen);
  percent_diff_tau3_rootmean_plot->SetLineColor(kGreen);
  percent_diff_tau3_rootmean_plot->SetLineWidth(3);
  // percent_diff_tau3_median_plot->SetMarkerStyle(3);

  percent_diff_tau3_hists->Add(percent_diff_tau3_mean_plot);
  percent_diff_tau3_hists->Add(percent_diff_tau3_median_plot);
  // percent_diff_tau3_hists->Add(percent_diff_tau3_meansq_plot);
  // percent_diff_tau3_hists->Add(percent_diff_tau3_rootmean_plot);
  percent_diff_tau3_hists->SetMinimum(-0.5);
  percent_diff_tau3_hists->SetMaximum(0.5);
  percent_diff_tau3_hists->Draw("AL");
  percent_diff_tau3_hists->SetTitle("% Difference Between Minimum and ``Optimal'' #tau_{3}");
  percent_diff_tau3_hists->GetXaxis()->SetTitle("#beta");
  percent_diff_tau3_hists->GetYaxis()->SetTitle("% difference");
  percent_diff_tau3_hists->GetYaxis()->SetLimits(-0.5, 0.5);
  TLegend *leg_percent_diff_tau3 = new TLegend(0.14, 0.7, 0.4, 0.88);
  leg_percent_diff_tau3->SetFillColor(kWhite);
  leg_percent_diff_tau3->SetLineColor(kWhite);
  leg_percent_diff_tau3->AddEntry(percent_diff_tau3_mean_plot, "Mean", "L");
  leg_percent_diff_tau3->AddEntry(percent_diff_tau3_median_plot, "Median", "L");
  // leg_percent_diff_tau3->AddEntry(percent_diff_tau3_meansq_plot, "Mean^{2}", "L");
  // leg_percent_diff_tau3->AddEntry(percent_diff_tau3_rootmean_plot, "#sqrt{Mean}", "L");
  leg_percent_diff_tau3->Draw();
  percent_diff_tau3_hists_can->Write();
  percent_diff_tau3_hists_can->Print("percentdiff_tau3.eps", "eps");

  TMultiGraph *min_tau3_values_delta = new TMultiGraph(); 
  TMultiGraph *min_tau3_values_power = new TMultiGraph(); 
  TGraph *min_beta_tau3_delta = new TGraph(n_betas, min_beta_tau3_values, min_delta_tau3_values);
  TGraph *min_beta_tau3_power = new TGraph(n_betas, min_beta_tau3_values, min_power_tau3_values);
  TGraph *min_beta_tau3_delta_meansq = new TGraph(n_betas, min_beta_tau3_values_meansq, min_delta_values_meansq);
  TGraph *min_beta_tau3_power_meansq = new TGraph(n_betas, min_beta_tau3_values_meansq, min_power_values_meansq);
  TGraph *min_beta_tau3_delta_rootmean = new TGraph(n_betas, min_beta_tau3_values_rootmean, min_delta_values_rootmean);
  TGraph *min_beta_tau3_power_rootmean = new TGraph(n_betas, min_beta_tau3_values_rootmean, min_power_values_rootmean);
  TGraph *min_beta_tau3_delta_median = new TGraph(n_betas, min_beta_tau3_values_median, min_delta_values_median);
  TGraph *min_beta_tau3_power_median = new TGraph(n_betas, min_beta_tau3_values_median, min_power_values_median);

  min_tau3_values_delta->SetTitle("#beta and #delta for Minimum #tau_{3}");
  min_tau3_values_power->SetTitle("#beta and p for Minimum #tau_{3}");
  // min_beta_tau3_delta->SetMarkerStyle(3);
  // min_beta_tau3_power->SetMarkerStyle(3);
  min_beta_tau3_delta->SetMarkerSize(2);
  min_beta_tau3_power->SetMarkerSize(2);
  min_beta_tau3_delta->SetLineColor(kRed);
  min_beta_tau3_power->SetLineColor(kRed);
  min_beta_tau3_delta->SetLineWidth(3);
  min_beta_tau3_power->SetLineWidth(3);
  // min_beta_tau3_delta_meansq->SetMarkerStyle(3);
  // min_beta_tau3_power_meansq->SetMarkerStyle(3);
  min_beta_tau3_delta_meansq->SetMarkerSize(2);
  min_beta_tau3_power_meansq->SetMarkerSize(2);
  // min_beta_tau3_delta_meansq->SetMarkerColor(kBlue);
  // min_beta_tau3_power_meansq->SetMarkerColor(kBlue);
  min_beta_tau3_delta_meansq->SetLineColor(kOrange);
  min_beta_tau3_power_meansq->SetLineColor(kOrange);
  min_beta_tau3_delta_meansq->SetLineWidth(3);
  min_beta_tau3_power_meansq->SetLineWidth(3);
  // min_beta_tau3_delta_rootmean->SetMarkerStyle(3);
  // min_beta_tau3_power_rootmean->SetMarkerStyle(3);
  min_beta_tau3_delta_rootmean->SetMarkerSize(2);
  min_beta_tau3_power_rootmean->SetMarkerSize(2);
  // min_beta_tau3_delta_rootmean->SetMarkerColor(kGreen);
  // min_beta_tau3_power_rootmean->SetMarkerColor(kGreen);
  min_beta_tau3_delta_rootmean->SetLineColor(kGreen);
  min_beta_tau3_power_rootmean->SetLineColor(kGreen);
  min_beta_tau3_delta_rootmean->SetLineWidth(3);
  min_beta_tau3_power_rootmean->SetLineWidth(3);
  // min_beta_tau3_delta_median->SetMarkerStyle(3);
  // min_beta_tau3_power_median->SetMarkerStyle(3);
  min_beta_tau3_delta_median->SetMarkerSize(2);
  min_beta_tau3_power_median->SetMarkerSize(2);
  // min_beta_tau3_delta_median->SetMarkerColor(kRed);
  // min_beta_tau3_power_median->SetMarkerColor(kRed);
  min_beta_tau3_delta_median->SetLineColor(kBlue);
  min_beta_tau3_power_median->SetLineColor(kBlue);
  min_beta_tau3_delta_median->SetLineWidth(3);
  min_beta_tau3_power_median->SetLineWidth(3);

  TCanvas* beta_tau3_delta_can = new TCanvas("beta_tau3_delta_can", "beta_tau3_delta_can", 600, 600);
  beta_tau3_delta_can->cd();
  min_tau3_values_delta->Add(min_beta_tau3_delta);
  min_tau3_values_delta->Add(min_beta_tau3_delta_median);
  // min_tau3_values_delta->Add(min_beta_tau3_delta_meansq);
  // min_tau3_values_delta->Add(min_beta_tau3_delta_rootmean);
  min_tau3_values_delta->Draw("AL");
  beta_delta_func->Draw("SAME");
  min_tau3_values_delta->GetXaxis()->SetTitle("#beta");
  min_tau3_values_delta->GetYaxis()->SetTitle("#delta");
  min_tau3_values_delta->SetMaximum(10);
  TLegend *leg_beta_tau3_delta = new TLegend(0.55, 0.55, 0.88, 0.88);
  leg_beta_tau3_delta->SetFillColor(kWhite);
  leg_beta_tau3_delta->SetLineColor(kWhite);
  leg_beta_tau3_delta->AddEntry(min_beta_tau3_delta, "Mean", "L");
  leg_beta_tau3_delta->AddEntry(min_beta_tau3_delta_median, "Median", "L");
  // leg_beta_tau3_delta->AddEntry(min_beta_tau3_delta_meansq, "Mean^{2}", "L");
  // leg_beta_tau3_delta->AddEntry(min_beta_tau3_delta_rootmean, "#sqrt{Mean}", "L");
  leg_beta_tau3_delta->AddEntry(beta_delta_func, "Prediction", "L");
  leg_beta_tau3_delta->Draw();
  beta_tau3_delta_can->Write();
  beta_tau3_delta_can->Print("betavdelta_tau3.eps", "eps");

  TCanvas* beta_tau3_power_can = new TCanvas("beta_tau3_power_can", "beta_tau3_power_can", 600, 600);
  beta_tau3_power_can->cd();
  min_tau3_values_power->Add(min_beta_tau3_power);
  min_tau3_values_power->Add(min_beta_tau3_power_median);
  // min_tau3_values_power->Add(min_beta_tau3_power_meansq);
  // min_tau3_values_power->Add(min_beta_tau3_power_rootmean);
  min_tau3_values_power->Draw("AL");
  beta_power_func->Draw("SAME");
  min_tau3_values_power->GetXaxis()->SetTitle("#beta");
  min_tau3_values_power->GetYaxis()->SetTitle("power");
  TLegend *leg_beta_tau3_power = new TLegend(0.55, 0.55, 0.88, 0.88);
  leg_beta_tau3_power->SetFillColor(kWhite);
  leg_beta_tau3_power->SetLineColor(kWhite);
  leg_beta_tau3_power->AddEntry(min_beta_tau3_power, "Mean", "L");
  leg_beta_tau3_power->AddEntry(min_beta_tau3_power_median, "Median", "L");
  // leg_beta_tau3_power->AddEntry(min_beta_tau3_power_meansq, "Mean^{2}", "L");
  // leg_beta_tau3_power->AddEntry(min_beta_tau3_power_rootmean, "#sqrt{Mean}", "L");
  leg_beta_tau3_power->AddEntry(beta_power_func, "Prediction", "L");
  leg_beta_tau3_power->Draw();
  beta_tau3_power_can->Write();
  beta_tau3_power_can->Print("betavpower_tau3.eps", "eps");

  // for (int i = 0; i < mean_values.size(); i++) {
  //   cout << mean_values[i] << " & ";
  // }
  // cout << endl << endl;
  // for (int i = 0; i < rms_values.size(); i++) {
  //   cout << rms_values[i] << " & ";
  // }
  // cout << endl << endl;
  // for (int i = 0; i < mass_ratios.size(); i++) {
  //   cout << mass_ratios[i] << " & ";
  // }
  // cout << endl;

  out.Close();

}