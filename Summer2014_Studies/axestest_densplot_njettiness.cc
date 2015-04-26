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
#include "TColor.h"
#include "TString.h"
#include "TPaletteAxis.h"

using namespace std;
using namespace Pythia8;
using namespace fastjet;
using namespace fastjet::contrib;

class AxesStruct {

private:
  // Shared Ptr so it handles memory management
  SharedPtr<AxesDefinition> _axes_def;

public:
  AxesStruct(const AxesDefinition & axes_def)
  : _axes_def(axes_def.create()) {}
  
  // Need special copy constructor to make it possible to put in a std::vector
  AxesStruct(const AxesStruct& myStruct)
  : _axes_def(myStruct._axes_def->create()) {}

  const AxesDefinition & def() const {return *_axes_def;}
  string description() const {return _axes_def->description();}
  string short_description() const {return _axes_def->short_description();}
  
};

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

  TFile out("axestest_densplot_colortest.root", "RECREATE");

  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW"); // Input File Name
  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW"); // Input File Name


  // Int_t palette_red[10];
  // Int_t palette_blue[5];
  // palette_red[0] = kWhite;
  // // palette_blue[0] = kWhite;
  // for (int i = 0; i < 10; i++) {
  //     palette_red[10-i] = kRed - i;
  // }

  Int_t palette_red[100];
  Double_t r_red[]    = {1., 1.0};
  Double_t g_red[]    = {1., 0.0};
  Double_t b_red[]    = {1., 0.0};
  Double_t stop_red[] = {0., 1.0};
  Int_t FI_red = TColor::CreateGradientColorTable(2, stop_red, r_red, g_red, b_red, 100);

  Int_t palette_blue[100];
  Double_t r_blue[]    = {1., 0.0};
  Double_t g_blue[]    = {1., 0.0};
  Double_t b_blue[]    = {1., 1.0};
  Double_t stop_blue[] = {0., 1.0};
  Int_t FI_blue = TColor::CreateGradientColorTable(2, stop_blue, r_blue, g_blue, b_blue, 100);
  
  for (int i=0;i<100;i++) {
    palette_blue[i] = FI_blue+i;
    palette_red[i] = FI_red+i;
  }

  TStyle *plain  = new TStyle("Plain","plain");

  plain->SetLegendBorderSize(0);
  plain->SetLegendFillColor(0);
  plain->SetTitleFont(132, "a");
  plain->SetTitleFont(132, "xy");
  plain->SetLegendFont(132);
  // plain->SetTextSize(0.05);
  plain->SetLabelSize(0.05, "xy");
  // plain->SetLabelOffset(0.003, "xy");
  plain->SetTitleSize(0.08, "a");
  plain->SetTitleSize(0.07, "xy");
  plain->SetPadLeftMargin(0.12);
  plain->SetPadBottomMargin(0.12);
  plain->SetTitleOffset(0.8, "xy");
  plain->SetHistLineWidth(3);

  // plain->SetLegendSize(12);
  plain->SetTitleBorderSize(0);
  plain->SetTitleX(0.1f);
  plain->SetTitleW(0.8f);
  plain->SetOptStat(0);
  plain->cd();

  // gROOT->ForceStyle();

  gROOT->SetStyle("plain");

  double Rcutoff = 0.5;
  int njettiness = 6;

  vector<double> mass_ratios;
  vector<double> mean_values;
  vector<double> rms_values;

  const int n_event = 100; // # of events to be analyzed   
  const int n_powers = 16;
  const int n_deltas = 17;
  // const int n_betas = 6;
  const int n_betas = 2;
  double epsilon = 0.001;

  const double infinity = std::numeric_limits<double>::max();

  // double optimal_power_values[n_betas] = {4.0, 2.0, 1.0, 0.67, 0.5, 0.33, 0.25, 0.1};
  // double optimal_delta_values[n_betas] = {100.0, 100.0, 100.0, 2.0, 1.0, 0.5, 0.33, 0.11};
  // double beta_values[n_betas] = {0.25, 0.5, 1.0, 1.5, 2.0, 4.0};
  double beta_values[n_betas] = {1.0, 2.0};
  // double beta_values[n_betas] = {0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0};
  double optimal_njet_values_mean[n_betas];
  double optimal_njet_values_meansq[n_betas];
  double optimal_njet_values_rootmean[n_betas];
  vector<vector<double>> optimal_njet_values_median;

  // 1st component is jets vs event, second component is no min vs with min, last component is for the different betas
  int percent_overlap[4][2][n_betas];

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < n_betas; k++) {
        percent_overlap[i][j][k] = 0;
      }
    }
  }

  TObjArray delta_power_density_plots;
  TObjArray delta_power_density_nopass_plots;

  for (int i_beta = 0; i_beta < n_betas; i_beta++) {
    // TH2* delta_power_density = new TH2F("delta_power_density", "", 20, 0, 5, 20, 0, 5);
    // TH2* delta_power_density_nopass = new TH2F("delta_power_density_nopass", "", 20, 0, 5, 20, 0, 5);
    TH2* delta_power_density = new TH2F("delta_power_density", "", 16, 0, 4, 16, 0, 4);
    TH2* delta_power_density_nopass = new TH2F("delta_power_density_nopass", "", 16, 0, 4, 16, 0, 4);
    delta_power_density_plots.Add(delta_power_density);
    delta_power_density_nopass_plots.Add(delta_power_density_nopass);
  }


  int total_ttbar_jets = 0;
  int total_dijets_jets = 0;

  // int overlap_power[n_powers][n_betas] = {0};
  // int overlap_delta[n_deltas][n_betas] = {0};
  int overlap_total[n_powers][n_deltas][n_betas] = {0};
  int overlap_minjets[n_powers][n_deltas][n_betas] = {0};

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
          event_display->Fill(jet_constituents[j].eta(), jet_constituents[j].phi(), jet_constituents[j].perp());
        }

        cout << i_event << endl;

        double beta;
        for (int i_beta = 0; i_beta < n_betas; i_beta++) {
          beta = beta_values[i_beta];

          // FIRST ESTABLISH THE "OTIMAL" AXIS CONFIGURATION
          double optimal_power = (double)1/(beta);
          double optimal_delta;
          if (beta > 1) optimal_delta = (double)1/(beta - 1);
          else optimal_delta = std::numeric_limits<int>::max();

          AxesStruct *axes_finder_optimal = new AxesStruct(GenRecomb_GenKT_Axes(optimal_delta, optimal_power, Rcutoff));

          AxesStruct *axes_finder_onepass_optimal;
          if (beta >= 1 && beta <= 3) axes_finder_onepass_optimal = new AxesStruct(OnePass_GenRecomb_GenKT_Axes(optimal_delta, optimal_power, Rcutoff));
          else axes_finder_onepass_optimal = axes_finder_optimal;

          // UnnormalizedCutoffMeasure measure_function = UnnormalizedCutoffMeasure(beta, Rcutoff);
          XConeCutoffMeasure measure_function = XConeCutoffMeasure(beta, Rcutoff);

          NjettinessPlugin njet_plugin_optimal(njettiness, axes_finder_onepass_optimal->def(), measure_function);
          JetDefinition njet_def_optimal(&njet_plugin_optimal);
          ClusterSequence njet_cluster_optimal(jet_constituents, njet_def_optimal);
          // vector<PseudoJet> optimal_jets = njet_cluster_optimal.inclusive_jets();
          const NjettinessExtras *extras_optimal = njettiness_extras(njet_cluster_optimal);
          vector<PseudoJet> optimal_jets = extras_optimal->jets();
          vector<PseudoJet> optimal_axes = extras_optimal->axes();


          NjettinessPlugin njet_plugin_optimal_nopass(njettiness, axes_finder_optimal->def(), measure_function);
          JetDefinition njet_def_optimal_nopass(&njet_plugin_optimal_nopass);
          ClusterSequence njet_cluster_optimal_nopass(jet_constituents, njet_def_optimal_nopass);
          // vector<PseudoJet> optimal_jets = njet_cluster_optimal.inclusive_jets();
          const NjettinessExtras *extras_optimal_nopass = njettiness_extras(njet_cluster_optimal_nopass);
          vector<PseudoJet> optimal_jets_nopass = extras_optimal_nopass->jets();
          vector<PseudoJet> optimal_axes_nopass = extras_optimal_nopass->axes();


          vector<PseudoJet> min_jets;
          vector<PseudoJet> min_axes;
          double min_tau = std::numeric_limits<int>::max();
          vector<PseudoJet> total_min_jets[n_powers][n_deltas];

          vector<PseudoJet> min_jets_nopass;
          vector<PseudoJet> min_axes_nopass;
          double min_tau_nopass = std::numeric_limits<int>::max();
          vector<PseudoJet> total_min_jets_nopass[n_powers][n_deltas];

          for (int i_power = 0; i_power < (n_powers); i_power++) {
            // double power = (double)(i_power + 1)/8;
            double power = (double)(i_power + 1)/4 - epsilon;
            double delta;

            for (int i_delta = 0; i_delta < (n_deltas); i_delta++) {
              // delta = (double)(i_delta + 1)/8;
              delta = (double)(i_delta + 1)/4 - epsilon;
              if (i_delta == (n_deltas - 1)) delta = infinity;

              AxesStruct *axes_finder = new AxesStruct(GenRecomb_GenKT_Axes(delta, power, Rcutoff));

              AxesStruct *axes_finder_onepass;
              if (beta >= 1 && beta <= 3) axes_finder_onepass = new AxesStruct(OnePass_GenRecomb_GenKT_Axes(delta, power, Rcutoff));
              else axes_finder_onepass = axes_finder;

              NjettinessPlugin njet_plugin_min_manual(njettiness, axes_finder_onepass->def(), measure_function);
              JetDefinition njet_def_min_manual(&njet_plugin_min_manual);
              ClusterSequence njet_cluster_min_manual(jet_constituents, njet_def_min_manual);
              // vector<PseudoJet> min_manual_jets = njet_cluster_min_manual.inclusive_jets();
              const NjettinessExtras *extras_min_manual = njettiness_extras(njet_cluster_min_manual);
              vector<PseudoJet> min_manual_jets = extras_min_manual->jets();
              vector<PseudoJet> min_manual_axes = extras_min_manual->axes();
              double min_manual_tau = extras_min_manual->totalTau();

              total_min_jets[i_power][i_delta] = min_manual_axes;

              if (min_manual_tau < min_tau) {
                min_tau = min_manual_tau;
                min_jets = min_manual_jets;
                min_axes = min_manual_axes;
              }


              NjettinessPlugin njet_plugin_min_manual_nopass(njettiness, axes_finder->def(), measure_function);
              JetDefinition njet_def_min_manual_nopass(&njet_plugin_min_manual_nopass);
              ClusterSequence njet_cluster_min_manual_nopass(jet_constituents, njet_def_min_manual_nopass);
              // vector<PseudoJet> min_manual_nopass_jets = njet_cluster_min_manual_nopass.inclusive_jets();
              const NjettinessExtras *extras_min_manual_nopass = njettiness_extras(njet_cluster_min_manual_nopass);
              vector<PseudoJet> min_manual_nopass_jets = extras_min_manual_nopass->jets();
              vector<PseudoJet> min_manual_nopass_axes = extras_min_manual_nopass->axes();
              double min_manual_nopass_tau = extras_min_manual_nopass->totalTau();

              total_min_jets_nopass[i_power][i_delta] = min_manual_nopass_axes;

              if (min_manual_nopass_tau < min_tau_nopass) {
                min_tau_nopass = min_manual_nopass_tau;
                min_jets_nopass = min_manual_nopass_jets;
                min_axes_nopass = min_manual_nopass_axes;
              }

              // mean_tau6_ttbar_values[i_beta] += min_manual_tau;
              // meansq_tau6_ttbar_values[i_beta] += min_manual_tau*min_manual_tau;
              // rootmean_tau6_ttbar_values[i_beta] += TMath::Sqrt(min_manual_tau);
              // median_tau6_ttbar_values[i_beta].push_back(min_manual_tau);

              // for (int i_jet = 0; i_jet < min_manual_axes.size(); i_jet++) {
              //   double min_distance = std::numeric_limits<int>::max();
              //   PseudoJet closest_jet;
              //   for (int j_jet = 0; j_jet < optimal_axes.size(); j_jet++) {
              //     if (min_manual_axes[i_jet].delta_R(optimal_axes[j_jet]) < min_distance) {
              //       min_distance = min_manual_axes[i_jet].delta_R(optimal_axes[j_jet]);
              //       closest_jet = optimal_axes[j_jet];
              //     }
              //   }
              //   if (min_manual_axes[i_jet].delta_R(closest_jet) < Rcutoff/2.0) overlap_total[i_power][i_delta][i_beta]++;
              // }

              // //LOOK AT INDIVIDUAL ANTI-KT JETS
              // vector <PseudoJet> inclusiveJets_std, sortedJets_std, hardestJets_std;

              // double Rparam_std = 1.0;
              // Strategy strategy_std = Best;
              // RecombinationScheme recombScheme_std = E_scheme;
              // JetDefinition *jetDef_std = new JetDefinition(antikt_algorithm, Rparam_std, recombScheme_std, strategy_std);
              // ClusterSequence clustSeq_std(input_particles, *jetDef_std);
              // inclusiveJets_std = clustSeq_std.inclusive_jets(200.0);
              // sortedJets_std = sorted_by_pt(inclusiveJets_std);
              // Selector jet_selector = SelectorNHardest(2);
              // hardestJets_std = jet_selector(sortedJets_std);

              // for (int i = 0; i < hardestJets_std.size(); i++) {

              //   ClusterSequence clustSeq_subjet(hardestJets_std[i].constituents(), *jetDef);
              //   vector<PseudoJet> starting_subjets = clustSeq_subjet.exclusive_jets(3);
              //   JetInformation subjet_set3 = findMinAxes(hardestJets_std[i].constituents(), starting_subjets, 3, beta, Rcutoff);
              //   JetInformation subjet_set2 = findMinAxes(hardestJets_std[i].constituents(), starting_subjets, 2, beta, Rcutoff);
              //   vector<PseudoJet> subjet_axes3 = subjet_set3.axes;
              //   vector<PseudoJet> subjet_axes2 = subjet_set2.axes;
              //   mean_tau3_ttbar_values[i_beta] += subjet_set3.tau;
              //   meansq_tau3_ttbar_values[i_beta] += pow(subjet_set3.tau, 2);
              //   rootmean_tau3_ttbar_values[i_beta] += pow(subjet_set3.tau, 0.5);
              //   median_tau3_ttbar_values[i_beta].push_back(subjet_set3.tau);

            }
          }

          double smallest_dR = 0.1;

          int jet_counter_nopass = 0;
          int jet_counter_onepass = 0;

          vector<PseudoJet> temp_optimal_axes = optimal_axes;
          vector<PseudoJet> temp_min_axes = min_axes;

          if (min_axes.size() > 0 && optimal_axes.size() > 0) {

            while (temp_optimal_axes.size() > 0) {
              double min_distance = std::numeric_limits<int>::max();
              PseudoJet closest_min_jet, closest_temp_jet;
              int min_index, temp_index;
              for (int i_jet = 0; i_jet < temp_optimal_axes.size(); i_jet++) {
                for (int j_jet = 0; j_jet < temp_min_axes.size(); j_jet++) {

                  if (temp_optimal_axes[i_jet].delta_R(temp_min_axes[j_jet]) < min_distance) {

                    min_distance = temp_optimal_axes[i_jet].delta_R(temp_min_axes[j_jet]);

                    closest_temp_jet = temp_optimal_axes[i_jet];
                    closest_min_jet = temp_min_axes[j_jet];
                    temp_index = i_jet;
                    min_index = j_jet;
                  }
                }
              }

              if (closest_min_jet.delta_R(closest_temp_jet) < smallest_dR) jet_counter_onepass++;

              std::swap(temp_optimal_axes.at(temp_index), temp_optimal_axes.back());
              std::swap(temp_min_axes.at(min_index), temp_min_axes.back());
              temp_optimal_axes.pop_back();
              temp_min_axes.pop_back();

            } 
          }


          vector<PseudoJet> temp_optimal_axes_nopass = optimal_axes_nopass;
          vector<PseudoJet> temp_min_axes_nopass = min_axes;

          if (min_axes.size() > 0 && optimal_axes_nopass.size() > 0) {

            while (temp_optimal_axes_nopass.size() > 0) {
              double min_distance = std::numeric_limits<int>::max();
              PseudoJet closest_min_jet, closest_temp_jet;
              int min_index, temp_index;
              for (int i_jet = 0; i_jet < temp_optimal_axes_nopass.size(); i_jet++) {
                for (int j_jet = 0; j_jet < temp_min_axes_nopass.size(); j_jet++) {

                  if (temp_optimal_axes_nopass[i_jet].delta_R(temp_min_axes_nopass[j_jet]) < min_distance) {

                    min_distance = temp_optimal_axes_nopass[i_jet].delta_R(temp_min_axes_nopass[j_jet]);

                    closest_temp_jet = temp_optimal_axes_nopass[i_jet];
                    closest_min_jet = temp_min_axes_nopass[j_jet];
                    temp_index = i_jet;
                    min_index = j_jet;
                  }
                }
              }

              if (closest_min_jet.delta_R(closest_temp_jet) < smallest_dR) jet_counter_nopass++;

              std::swap(temp_optimal_axes_nopass.at(temp_index), temp_optimal_axes_nopass.back());
              std::swap(temp_min_axes_nopass.at(min_index), temp_min_axes_nopass.back());
              temp_optimal_axes_nopass.pop_back();
              temp_min_axes_nopass.pop_back();

            } 
          }

          percent_overlap[0][0][i_beta] += jet_counter_nopass;
          if (jet_counter_nopass >= (optimal_axes.size() - 2)) percent_overlap[1][0][i_beta]++;
          if (jet_counter_nopass >= (optimal_axes.size() - 1)) percent_overlap[2][0][i_beta]++;
          if (jet_counter_nopass == optimal_axes.size()) percent_overlap[3][0][i_beta]++;

          percent_overlap[0][1][i_beta] += jet_counter_onepass;
          if (jet_counter_onepass >= (optimal_axes.size() - 2)) percent_overlap[1][1][i_beta]++;
          if (jet_counter_onepass >= (optimal_axes.size() - 1)) percent_overlap[2][1][i_beta]++;
          if (jet_counter_onepass == optimal_axes.size()) percent_overlap[3][1][i_beta]++;

          total_ttbar_jets += min_axes.size();

          TH2* delta_power_density = (TH2F*)delta_power_density_plots.At(i_beta);
          TH2* delta_power_density_nopass = (TH2F*)delta_power_density_nopass_plots.At(i_beta);

          for (int i_power = 0; i_power < (n_powers); i_power++) {
            double power = (double)(i_power + 1)/4 - epsilon;
            double delta;

            for (int i_delta = 0; i_delta < (n_deltas); i_delta++) {
              // delta = (double)(i_delta + 1)/8;
              delta = (double)(i_delta + 1)/4 - epsilon;
              if (i_delta == (n_deltas - 1)) delta = infinity;
  
              vector<PseudoJet> delta_power_axes = total_min_jets[i_power][i_delta];
              vector<PseudoJet> temp_delta_power_axes = delta_power_axes;
              vector<PseudoJet> temp_min_axes = min_axes;

              int event_counter_onepass = 0;

              if (min_axes.size() > 0 && delta_power_axes.size() > 0) {

                while (temp_delta_power_axes.size() > 0) {
                  double min_distance = std::numeric_limits<int>::max();
                  PseudoJet closest_min_jet, closest_temp_jet;
                  int min_index, temp_index;
                  for (int i_jet = 0; i_jet < temp_delta_power_axes.size(); i_jet++) {
                    for (int j_jet = 0; j_jet < temp_min_axes.size(); j_jet++) {

                      if (temp_delta_power_axes[i_jet].delta_R(temp_min_axes[j_jet]) < min_distance) {

                        min_distance = temp_delta_power_axes[i_jet].delta_R(temp_min_axes[j_jet]);

                        closest_temp_jet = temp_delta_power_axes[i_jet];
                        closest_min_jet = temp_min_axes[j_jet];
                        temp_index = i_jet;
                        min_index = j_jet;
                      }
                    }
                  }

                  if (closest_min_jet.delta_R(closest_temp_jet) < smallest_dR) delta_power_density->Fill(delta, power);

                  std::swap(temp_delta_power_axes.at(temp_index), temp_delta_power_axes.back());
                  std::swap(temp_min_axes.at(min_index), temp_min_axes.back());
                  temp_delta_power_axes.pop_back();
                  temp_min_axes.pop_back();

                  // temp_delta_power_axes.erase(temp_delta_power_axes.begin()+temp_index);
                  // temp_min_axes.erase(temp_min_axes.begin()+min_index);

                } 
              }

              vector<PseudoJet> delta_power_axes_nopass = total_min_jets_nopass[i_power][i_delta];
              vector<PseudoJet> temp_delta_power_axes_nopass = delta_power_axes_nopass;
              vector<PseudoJet> temp_min_axes_nopass = min_axes;

              int event_counter_nopass = 0;

              if (min_axes.size() > 0 && delta_power_axes_nopass.size() > 0) {

                while (temp_delta_power_axes_nopass.size() > 0) {
                  double min_distance = std::numeric_limits<int>::max();
                  PseudoJet closest_min_jet, closest_temp_jet;
                  int min_index, temp_index;
                  for (int i_jet = 0; i_jet < temp_delta_power_axes_nopass.size(); i_jet++) {
                    for (int j_jet = 0; j_jet < temp_min_axes_nopass.size(); j_jet++) {

                      if (temp_delta_power_axes_nopass[i_jet].delta_R(temp_min_axes_nopass[j_jet]) < min_distance) {

                        min_distance = temp_delta_power_axes_nopass[i_jet].delta_R(temp_min_axes_nopass[j_jet]);

                        closest_temp_jet = temp_delta_power_axes_nopass[i_jet];
                        closest_min_jet = temp_min_axes_nopass[j_jet];
                        temp_index = i_jet;
                        min_index = j_jet;
                      }
                    }
                  }

                  if (closest_min_jet.delta_R(closest_temp_jet) < smallest_dR) delta_power_density_nopass->Fill(delta, power);

                  std::swap(temp_delta_power_axes_nopass.at(temp_index), temp_delta_power_axes_nopass.back());
                  std::swap(temp_min_axes_nopass.at(min_index), temp_min_axes_nopass.back());
                  temp_delta_power_axes_nopass.pop_back();
                  temp_min_axes_nopass.pop_back();

                  // temp_delta_power_axes_nopass.erase(temp_delta_power_axes_nopass.begin()+temp_index);
                  // temp_min_axes_nopass.erase(temp_min_axes_nopass.begin()+min_index);

                } 
              }

              // for (int i_jet = 0; i_jet < delta_power_axes.size(); i_jet++) {
              //   double min_distance = std::numeric_limits<int>::max();
              //   PseudoJet closest_jet;
              //   for (int j_jet = 0; j_jet < min_axes.size(); j_jet++) {
              //     if (delta_power_axes[i_jet].delta_R(min_axes[j_jet]) < min_distance) {
              //       min_distance = delta_power_axes[i_jet].delta_R(min_axes[j_jet]);
              //       closest_jet = min_axes[j_jet];
              //     }
              //   }
              //   if (delta_power_axes[i_jet].delta_R(closest_jet) < Rcutoff/2.0) {
              //     delta_power_density->Fill(delta, power);
              //   }
              // }
            }
          }

            // for (int i_jet = 0; i_jet < min_jets.size(); i_jet++) {
            //   double min_distance = std::numeric_limits<int>::max();
            //   PseudoJet closest_jet;
            //   for (int j_jet = 0; j_jet < optimal_jets.size(); j_jet++) {
            //     if (min_jets[i_jet].delta_R(optimal_jets[j_jet]) < min_distance) {
            //       min_distance = min_jets[i_jet].delta_R(optimal_jets[j_jet]);
            //       closest_jet = optimal_jets[j_jet];
            //     }
            //   }
            //   if (min_jets[i_jet].delta_R(closest_jet) < Rcutoff/2.0) overlap_minjets[min_power][min_delta][i_beta]++;
            // }

            // do {

            //   double min_distance = std::numeric_limits<int>::max();
            //   PseudoJet closest_min_jet, closest_optimal_jet;
            //   int min_index, optimal_index;
            //   for (int i_jet = 0; i_jet < temp_min_axes.size(); i_jet++) {
            //     for (int j_jet = 0; j_jet < temp_optimal_axes.size(); j_jet++) {
            //       if (temp_min_axes[i_jet].delta_R(temp_optimal_axes[j_jet]) < min_distance) {
            //         min_distance = temp_min_axes[i_jet].delta_R(temp_optimal_axes[j_jet]);
            //         closest_min_jet = temp_min_axes[i_jet];
            //         min_index = i_jet;
            //         closest_optimal_jet = temp_optimal_axes[j_jet];
            //         optimal_index = j_jet;
            //       }
            //     }
            //   }

            //   if (closest_min_jet.delta_R(closest_optimal_jet) < Rcutoff/2.0) overlap_minjets[min_power][min_delta][i_beta]++;

            //   final_min_axes.push_back(closest_min_jet);
            //   final_optimal_axes.push_back(closest_optimal_jet);

            //   temp_min_axes.erase(temp_min_axes.begin()+min_index);
            //   temp_optimal_axes.erase(temp_optimal_axes.begin()+optimal_index);

            // } while (temp_min_axes.size() > 0);


        // event_display->SetStats(0);
        // axes_genrecomb_display->SetStats(0);
        // axes_manual_display->SetStats(0);
        // subjet1_display->SetStats(0);
        // subjet2_display->SetStats(0);
        // subjet3_display->SetStats(0);
        // subjet4_display->SetStats(0);
        // subjet5_display->SetStats(0);
        // subjet6_display->SetStats(0);
        // subjet7_display->SetStats(0);

        // TCanvas *display = new TCanvas("display", "Event Display", 1093, 700);
        // display->cd();
        // display->SetFixedAspectRatio();
        // event_display->GetXaxis()->SetTitle("#eta");
        // event_display->GetYaxis()->SetTitle("#phi");
        // event_display->SetFillColor(kBlack);
        // event_display->SetLineColor(kBlack);
        // event_display->SetLineWidth(1);
        // subjet1_display->SetMarkerColor(kRed);
        // subjet2_display->SetMarkerColor(kBlue);
        // subjet3_display->SetMarkerColor(kGreen);
        // subjet4_display->SetMarkerColor(kYellow);
        // subjet5_display->SetMarkerColor(kOrange);
        // subjet6_display->SetMarkerColor(6);
        // subjet7_display->SetMarkerColor(20);
        // subjet1_display->SetMarkerStyle(21);
        // subjet2_display->SetMarkerStyle(21);
        // subjet3_display->SetMarkerStyle(21);
        // subjet4_display->SetMarkerStyle(21);
        // subjet5_display->SetMarkerStyle(21);
        // subjet6_display->SetMarkerStyle(21);
        // subjet7_display->SetMarkerStyle(21);
        // subjet1_display->SetMarkerSize(0.5);
        // subjet2_display->SetMarkerSize(0.5);
        // subjet3_display->SetMarkerSize(0.5);
        // subjet4_display->SetMarkerSize(0.5);
        // subjet5_display->SetMarkerSize(0.5);
        // subjet6_display->SetMarkerSize(0.5);
        // subjet7_display->SetMarkerSize(0.5);

        // event_display->Draw("box");

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
        // if (i_sample == 0) title = "ttbar_display" + ss.str() + "_wtap8_beta0p125.eps";
        // if (i_sample == 1) title = "dijets_display" + ss.str() + "_wtap8_beta0p125.eps";

        // display->Print(title, "eps");
          }
        // delete display;
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

          // }
        // }

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

  total_ttbar_jets /= (n_betas);

  TH2* overlap_density_delta = new TH2F("overlap_density_delta", "", 20, 0, 5, 20, 0, 5);
  TH2* overlap_density_power = new TH2F("overlap_density_power", "", 20, 0, 5, 20, 0, 5);

  TH2* overlap_density_minjets_delta = new TH2F("overlap_density_minjets_delta", "", 20, 0, 5, 20, 0, 5);
  TH2* overlap_density_minjets_power = new TH2F("overlap_density_minjets_power", "", 20, 0, 5, 20, 0, 5);

  for (int i_beta = 0; i_beta < n_betas; i_beta++) {
  double beta = beta_values[i_beta];
  double power;
  double delta;

    for (int i_power = 0; i_power < (n_powers); i_power++) {
      power = (double)(i_power + 1)/4 - epsilon;

      for (int i_delta = 0; i_delta < (n_deltas - 1); i_delta++) {
        delta = (double)(i_delta + 1)/4 - epsilon;

        overlap_total[i_power][i_delta][i_beta] = (double)overlap_total[i_power][i_delta][i_beta]/(double)(total_ttbar_jets);
        // overlap_minjets[i_power][i_delta][i_beta] = (double)overlap_minjets[i_power][i_delta][i_beta]/(double)(total_ttbar_jets);

        overlap_density_delta->Fill(beta, delta, overlap_total[i_power][i_delta][i_beta]);
        overlap_density_power->Fill(beta, power, overlap_total[i_power][i_delta][i_beta]);
        overlap_density_minjets_delta->Fill(beta, delta, overlap_minjets[i_power][i_delta][i_beta]);
        overlap_density_minjets_power->Fill(beta, power, overlap_minjets[i_power][i_delta][i_beta]);
      }
    }

  }

  // overlap_density_power->SetStats(0);
  // overlap_density_delta->SetStats(0);

  // TCanvas *density_delta = new TCanvas("density_delta", "", 600, 600);
  // density_delta->cd();
  // overlap_density_delta->SetTitle("Density of True/Heuristic Minimum Overlap");
  // overlap_density_delta->GetXaxis()->SetTitle("#beta");
  // overlap_density_delta->GetYaxis()->SetTitle("#delta");
  // overlap_density_delta->Draw("COLZ");
  // density_delta->Write();

  // TCanvas *density_power = new TCanvas("density_power", "", 600, 600);
  // density_power->cd();
  // overlap_density_power->SetTitle("Density of True/Heuristic Minimum Overlap");
  // overlap_density_power->GetXaxis()->SetTitle("#beta");
  // overlap_density_power->GetYaxis()->SetTitle("p");
  // overlap_density_power->Draw("COLZ");
  // density_power->Write();


  // overlap_density_minjets_power->SetStats(0);
  // overlap_density_minjets_delta->SetStats(0);

  // TCanvas *density_minjets_delta = new TCanvas("density_minjets_delta", "", 600, 600);
  // density_minjets_delta->cd();
  // overlap_density_minjets_delta->SetTitle("Density of True/Heuristic Minimum Overlap");
  // overlap_density_minjets_delta->GetXaxis()->SetTitle("#beta");
  // overlap_density_minjets_delta->GetYaxis()->SetTitle("#delta");
  // overlap_density_minjets_delta->Draw("COLZ");
  // density_minjets_delta->Write();

  // TCanvas *density_minjets_power = new TCanvas("density_minjets_power", "", 600, 600);
  // density_minjets_power->cd();
  // overlap_density_minjets_power->SetTitle("Density of True/Heuristic Minimum Overlap");
  // overlap_density_minjets_power->GetXaxis()->SetTitle("#beta");
  // overlap_density_minjets_power->GetYaxis()->SetTitle("p");
  // overlap_density_minjets_power->Draw("COLZ");
  // density_minjets_power->Write();

  for (int i_beta = 0; i_beta < n_betas; i_beta++) {

    double beta = beta_values[i_beta];
    string beta_string;
    ostringstream ss;
    ss << beta;
    beta_string = ss.str();

    cout << "beta = " << beta << endl;
    cout << (double)percent_overlap[0][0][i_beta]/total_ttbar_jets << " " << (double)percent_overlap[0][1][i_beta]/total_ttbar_jets << endl;
    cout << (double)percent_overlap[1][0][i_beta]/n_event << " " << (double)percent_overlap[1][1][i_beta]/n_event << endl;
    cout << (double)percent_overlap[2][0][i_beta]/n_event << " " << (double)percent_overlap[2][1][i_beta]/n_event << endl;
    cout << (double)percent_overlap[3][0][i_beta]/n_event << " " << (double)percent_overlap[3][1][i_beta]/n_event << endl;
    cout << endl;

    TH2* delta_power_density = (TH2F*)delta_power_density_plots.At(i_beta);
    delta_power_density->SetStats(0);
    delta_power_density->Scale((double)1/total_ttbar_jets);

    if (i_beta == 0) gStyle->SetPalette(100, palette_red);
    else if (i_beta == 1) gStyle->SetPalette(100, palette_blue);

    TCanvas *density_minjets = new TCanvas("density_minjets", "", 600, 600);
    density_minjets->cd(); 
    // density_minjets->SetLogz();

    delta_power_density->SetTitle("Min. Axes for #beta = " + (TString)ss.str() + " (One pass)");
    delta_power_density->GetXaxis()->SetTitle("#delta");
    delta_power_density->GetYaxis()->SetTitle("p");
    delta_power_density->SetMinimum(0.75);
    delta_power_density->SetMaximum(1.0);
    delta_power_density->SetContour(100);
    delta_power_density->Draw("COLZ");
    density_minjets->SetRightMargin(0.15);
    gPad->Update();
    TPaletteAxis *palette = (TPaletteAxis*)delta_power_density->GetListOfFunctions()->FindObject("palette");
    palette->SetLabelSize(0.05);
    density_minjets->Write();
    // if (i_beta == 2 || i_beta == 3 || i_beta == 4) density_minjets->Print("delta_power_density_beta" + (TString)beta_string + "_color.eps", "eps");
    density_minjets->Print("delta_power_density_beta" + (TString)beta_string + "_color.eps", "eps");

    TH2* delta_power_density_nopass = (TH2F*)delta_power_density_nopass_plots.At(i_beta);
    delta_power_density_nopass->SetStats(0);
    delta_power_density_nopass->Scale((double)1/total_ttbar_jets);

    TCanvas *density_minjets_nopass = new TCanvas("density_minjets_nopass", "", 600, 600);
    density_minjets_nopass->cd();
    delta_power_density_nopass->SetTitle("Min. Axes for #beta = " + (TString)ss.str() + " (No pass)");
    delta_power_density_nopass->GetXaxis()->SetTitle("#delta");
    delta_power_density_nopass->GetYaxis()->SetTitle("p");
    delta_power_density_nopass->SetMinimum(0.75);
    delta_power_density_nopass->SetMaximum(1.0);
    delta_power_density_nopass->SetContour(100);
    delta_power_density_nopass->Draw("COLZ");
    density_minjets_nopass->SetRightMargin(0.15);
    gPad->Update();
    TPaletteAxis *palette_nopass = (TPaletteAxis*)delta_power_density_nopass->GetListOfFunctions()->FindObject("palette");
    palette_nopass->SetLabelSize(0.05);
    density_minjets_nopass->Write();
    // if (i_beta == 2 || i_beta == 3 || i_beta == 4) density_minjets_nopass->Print("delta_power_density_beta" + (TString)beta_string + "_nopass_color.eps", "eps");
    density_minjets_nopass->Print("delta_power_density_beta" + (TString)beta_string + "_nopass_color.eps", "eps");


    delete density_minjets;
    delete density_minjets_nopass;
  }

  out.Close();

}