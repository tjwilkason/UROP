//Misc. Headers
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <queue>
#include <unordered_map>

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
#include "TRint.h"
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

// double triangle_perimeter(PseudoJet axis1, PseudoJet axis2, PseudoJet axis3) {
// }

vector<PseudoJet> findMinAxes(vector<PseudoJet> input_particles, vector<PseudoJet> starting_axes, int njettiness, double beta, double Rcutoff) {

  vector<PseudoJet> min_manual_axes;

  std::string bitmask(njettiness, 1); // K leading 1's
  bitmask.resize(starting_axes.size(), 0); // N-K trailing 0's

  NjettinessPlugin njet_plugin_manual(njettiness, Manual_Axes(), UnnormalizedCutoffMeasure(beta, Rcutoff));
  JetDefinition njet_def_manual(&njet_plugin_manual);

  double min_manual_tau = std::numeric_limits<double>::max();
  do {
    vector<int> axis_indices;
    vector<PseudoJet> temp_axes;
    for (int i = 0; i < starting_axes.size(); ++i) {
      if (bitmask[i]) axis_indices.push_back(i);
    }
    for (int j = 0; j < axis_indices.size(); j++) {
      temp_axes.push_back(starting_axes[axis_indices[j]]);
    }

    double manual_tau = std::numeric_limits<int>::max();
    if (temp_axes.size() == njettiness) {
      njet_plugin_manual.setAxes(temp_axes);
      ClusterSequence njet_cluster_manual(input_particles, njet_def_manual);
      const NjettinessExtras *extras_manual = njettiness_extras(njet_cluster_manual);
      manual_tau = extras_manual->totalTau();

      if (manual_tau < min_manual_tau) {
        min_manual_tau = manual_tau;
        min_manual_axes = extras_manual->axes();
      }
    }
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  return min_manual_axes;
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

struct groupedJetsinfo {
  vector<PseudoJet> groupedJets;
  vector<vector<PseudoJet>> groupedJets_constituents;
  vector<int> numGroups;
  int triplet_counter;
};

groupedJetsinfo groupJets(vector<PseudoJet> inputJets, double Rgroup, double Rmin) {

  groupedJetsinfo jet_information;
  vector<PseudoJet> temp_bigjets(inputJets.size());
  int index_references[inputJets.size()];

  vector< vector<PseudoJet> > temp_bigjets_constituents;

  // initialize the index references to point to their own index
  for (int i_jets = 0; i_jets < inputJets.size(); i_jets++) {
    index_references[i_jets] = i_jets;
    temp_bigjets_constituents.push_back(vector<PseudoJet>());
  }

  // run through every combination of jets and reset index reference to the index of the closest jet with the lowest index
  for (int j_jets = 0; j_jets < inputJets.size(); j_jets++) {
    for (int k_jets = j_jets + 1; k_jets < inputJets.size(); k_jets++) {
      if ((inputJets[j_jets].delta_R(inputJets[k_jets]) < Rgroup) && (inputJets[j_jets].delta_R(inputJets[k_jets]) > Rmin)) {
        if (index_references[k_jets] < index_references[j_jets]) index_references[j_jets] = index_references[k_jets];
        else if (index_references[j_jets] < index_references[k_jets]) index_references[k_jets] = index_references[j_jets];
      }
    }
  }

  // cout << "hello first" << endl;
  // group the N subjets according to the index references from above
  vector<int> subjet_counter(inputJets.size());
  for (int i = 0; i < inputJets.size(); i++) {
    if (!temp_bigjets[index_references[i]].has_structure()) {
      temp_bigjets[index_references[i]] = inputJets[i];
      temp_bigjets_constituents[index_references[i]].push_back(inputJets[i]);
      subjet_counter[index_references[i]] = 1;
    }
    else {
      temp_bigjets[index_references[i]] = join(temp_bigjets[index_references[i]], inputJets[i]);
      temp_bigjets_constituents[index_references[i]].push_back(inputJets[i]);
      subjet_counter[index_references[i]]++;
    }
  }
  // cout << "hello!";

  // pick out only the combined jets that have mass (i.e. remove empty jets from the vector)
  vector<PseudoJet> massive_jets;
  vector<int> massive_subjet_counter;
  vector<vector<PseudoJet>> massive_jets_constituents;
  for (int i = 0; i < temp_bigjets.size(); i++) {
    if (temp_bigjets[i].m() > 0) {
      massive_jets.push_back(temp_bigjets[i]);
      massive_subjet_counter.push_back(subjet_counter[i]);
      massive_jets_constituents.push_back(temp_bigjets_constituents[i]);
    }
  }

  //count how many of the jets have 3 subjets
  int triplet_counter = 0;
  for (int i = 0; i < massive_subjet_counter.size(); i++) {
    if (massive_subjet_counter[i] >= 3) triplet_counter++;
  }

  jet_information.groupedJets = massive_jets;
  jet_information.groupedJets_constituents = massive_jets_constituents;
  jet_information.numGroups = massive_subjet_counter;
  jet_information.triplet_counter = triplet_counter;
  return jet_information;
}


vector<PseudoJet> findMinMass(vector<PseudoJet> initial_jets, int n_min) {

  PseudoJet minmass_jet;
  vector<PseudoJet> minmass_jets_const;

  std::string bitmask(n_min, 1); // K leading 1's
  bitmask.resize(initial_jets.size(), 0); // N-K trailing 0's

  double minmass = std::numeric_limits<int>::max();
  do {
    vector<int> axis_indices;
    for (int i = 0; i < initial_jets.size(); ++i) {
      if (bitmask[i]) axis_indices.push_back(i);
    }
    PseudoJet temp_jet;
    vector<PseudoJet> temp_jets_const;
    for (int j = 0; j < axis_indices.size(); j++) {
      temp_jet = join(temp_jet,initial_jets[axis_indices[j]]);
      temp_jets_const.push_back(initial_jets[axis_indices[j]]);
    }

    if (temp_jet.m() < minmass) {
      minmass = temp_jet.m();
      minmass_jet = temp_jet;
      minmass_jets_const = temp_jets_const;
    }
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  return minmass_jets_const;
}


vector<PseudoJet> findMaxPerp(vector<PseudoJet> initial_jets, int n_min) {

  PseudoJet maxperp_jet;
  vector<PseudoJet> maxperp_jets_const;

  std::string bitmask(n_min, 1); // K leading 1's
  bitmask.resize(initial_jets.size(), 0); // N-K trailing 0's

  double maxperp = 0;
  do {
    vector<int> axis_indices;
    for (int i = 0; i < initial_jets.size(); ++i) {
      if (bitmask[i]) axis_indices.push_back(i);
    }
    PseudoJet temp_jet;
    vector<PseudoJet> temp_jets_const;
    for (int j = 0; j < axis_indices.size(); j++) {
      temp_jet = join(temp_jet,initial_jets[axis_indices[j]]);
      temp_jets_const.push_back(initial_jets[axis_indices[j]]);
    }

    if (temp_jet.perp() > maxperp) {
      maxperp = temp_jet.perp();
      maxperp_jet = temp_jet;
      maxperp_jets_const = temp_jets_const;
    }
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  return maxperp_jets_const;
}

int main(int argc, char* argv[]){

  TFile out("ttbarstudy.root", "RECREATE");
  // TFile out("njettiness_testing_garbage.root", "RECREATE");

  TStyle *plain  = new TStyle("Plain","plain");
  plain->SetLegendBorderSize(0);
  plain->SetLegendFillColor(0);
  plain->SetTitleFont(132, "a");
  plain->SetTitleFont(132, "xy");
  plain->SetLegendFont(132);
  // plain->SetTextSize(0.05);
  plain->SetLabelSize(0.04, "x");
  plain->SetLabelSize(0.03, "y");
  plain->SetLabelOffset(0.003, "xy");
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
  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW"); // Input File Name
  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW"); // Input File Name

  double epsilon = 0.0001;
  double Rparam = 0.5;

  int initial_axes_num = 6;
  int antikt_initial_axes_num = 3;
  double ghost_maxrap = 5.0; // e.g. if particles go up to y=5

  vector<double> mass_ratios;
  vector<double> mean_values;
  vector<double> rms_values;

  //create list of various values of beta
  vector<double> betalist;
  // betalist.push_back(0.25);
  // betalist.push_back(0.5);
  betalist.push_back(1.0);
  betalist.push_back(2.0);
  // betalist.push_back(3.0);
  int n_betas = betalist.size();

  vector<Color_t> colorlist;
  colorlist.push_back(kRed);
  colorlist.push_back(kBlue);
  colorlist.push_back(kYellow);
  colorlist.push_back(8);

  // vector<int> ptlist;
  // ptlist.push_back(240);
  // ptlist.push_back(300);
  // ptlist.push_back(400);
  // ptlist.push_back(500);
  // ptlist.push_back(600);
  // ptlist.push_back(700);

  int n_perps = 6;
  // int n_perps = 1;
  double ptlist[n_perps];

  double top_efficiency_rawmass[n_betas + 2][n_perps];
  double top_efficiency_mass[n_betas + 2][n_perps];
  double top_efficiency_mass_2groups[n_betas + 2][n_perps];
  double top_efficiency_mass_2groups_improved[n_betas + 2][n_perps];
  double top_efficiency_mass_2groups_Wmasscut[n_betas + 2][n_perps];
  double top_efficiency_mass_2groups_improved_Wmasscut[n_betas + 2][n_perps];
  double top_efficiency_mass_3jet4jet_2groups[n_betas + 2][n_perps];
  double top_efficiency_mass_3jet4jet_2groups_Wmasscut[n_betas + 2][n_perps];

  double top_efficiency_ratio[n_betas][n_perps];
  double top_efficiency_ratio_2groups[n_betas][n_perps];

  double qcd_efficiency_rawmass[n_betas + 2][n_perps];
  double qcd_efficiency_mass[n_betas + 2][n_perps];
  double qcd_efficiency_mass_2groups[n_betas + 2][n_perps];
  double qcd_efficiency_mass_2groups_improved[n_betas + 2][n_perps];
  double qcd_efficiency_mass_2groups_Wmasscut[n_betas + 2][n_perps];
  double qcd_efficiency_mass_2groups_improved_Wmasscut[n_betas + 2][n_perps];
  double qcd_efficiency_mass_3jet4jet_2groups[n_betas + 2][n_perps];
  double qcd_efficiency_mass_3jet4jet_2groups_Wmasscut[n_betas + 2][n_perps];

  double qcd_efficiency_ratio[n_betas][n_perps];
  double qcd_efficiency_ratio_2groups[n_betas][n_perps];

  // double ROC_mass_compare[2][n_betas + 1][n_perps];
  TObjArray tau32_ttbar_manual_6jet_total_hists;
  TObjArray tau32_ttbar_manual_total_hists;
  TObjArray tau32_ttbar_manual_2groups_total_hists;
  TObjArray tau32_ttbar_manual_2groups_improved_total_hists;
  TObjArray tau32_ttbar_manual_2groups_bigcone_total_hists;
  TObjArray tau32_ttbar_antikt_total_hists;
  TObjArray tau32_ttbar_antikt_wta_total_hists;
  TObjArray tau32_dijets_manual_6jet_total_hists;
  TObjArray tau32_dijets_manual_total_hists;
  TObjArray tau32_dijets_manual_2groups_total_hists;
  TObjArray tau32_dijets_manual_2groups_improved_total_hists;
  TObjArray tau32_dijets_manual_2groups_bigcone_total_hists;
  TObjArray tau32_dijets_antikt_total_hists;
  TObjArray tau32_dijets_antikt_wta_total_hists;

  for (int i_beta = 0; i_beta < n_betas; i_beta++) {
    TH1* tau32_ttbar_manual_6jet_total = new TH1F("tau32_ttbar_manual_6jet_total", "#tau_{32} of Tops", 25, 0, 1);
    TH1* tau32_ttbar_manual_total = new TH1F("tau32_ttbar_manual_total", "#tau_{32} of Tops", 25, 0, 1);
    TH1* tau32_ttbar_manual_2groups_total = new TH1F("tau32_ttbar_manual_2groups_total", "#tau_{32} of Tops", 25, 0, 1);
    TH1* tau32_ttbar_manual_2groups_improved_total = new TH1F("tau32_ttbar_manual_2groups_improved_total", "#tau_{32} of Tops", 25, 0, 1);
    TH1* tau32_ttbar_manual_2groups_bigcone_total = new TH1F("tau32_ttbar_manual_2groups_bigcone_total", "", 25, 0, 1);
    TH1* tau32_ttbar_antikt_total = new TH1F("tau32_ttbar_antikt_total", "#tau_{32} of Tops", 25, 0, 1);
    TH1* tau32_ttbar_antikt_wta_total = new TH1F("tau32_ttbar_antikt_wta_total", "#tau_{32} of Tops", 25, 0, 1);
    TH1* tau32_dijets_manual_6jet_total = new TH1F("tau32_dijets_manual_6jet_total", "#tau_{32} of QCD", 25, 0, 1);
    TH1* tau32_dijets_manual_total = new TH1F("tau32_dijets_manual_total", "#tau_{32} of QCD", 25, 0, 1);
    TH1* tau32_dijets_manual_2groups_total = new TH1F("tau32_dijets_manual_2groups_total", "#tau_{32} of QCD", 25, 0, 1);
    TH1* tau32_dijets_manual_2groups_improved_total = new TH1F("tau32_dijets_manual_2groups_improved_total", "#tau_{32} of QCD", 25, 0, 1);
    TH1* tau32_dijets_manual_2groups_bigcone_total = new TH1F("tau32_dijets_manual_2groups_bigcone_total", "", 25, 0, 1);
    TH1* tau32_dijets_antikt_total = new TH1F("tau32_dijets_antikt_total", "#tau_{32} of QCD", 25, 0, 1);
    TH1* tau32_dijets_antikt_wta_total = new TH1F("tau32_dijets_antikt_wta_total", "#tau_{32} of QCD", 25, 0, 1);

    tau32_ttbar_manual_6jet_total_hists.Add(tau32_ttbar_manual_6jet_total);
    tau32_ttbar_manual_total_hists.Add(tau32_ttbar_manual_total);
    tau32_ttbar_manual_2groups_total_hists.Add(tau32_ttbar_manual_2groups_total);
    tau32_ttbar_manual_2groups_improved_total_hists.Add(tau32_ttbar_manual_2groups_improved_total);
    tau32_ttbar_manual_2groups_bigcone_total_hists.Add(tau32_ttbar_manual_2groups_bigcone_total);
    tau32_ttbar_antikt_total_hists.Add(tau32_ttbar_antikt_total);
    tau32_ttbar_antikt_wta_total_hists.Add(tau32_ttbar_antikt_wta_total);
    tau32_dijets_manual_6jet_total_hists.Add(tau32_dijets_manual_6jet_total);
    tau32_dijets_manual_total_hists.Add(tau32_dijets_manual_total);
    tau32_dijets_manual_2groups_total_hists.Add(tau32_dijets_manual_2groups_total);
    tau32_dijets_manual_2groups_improved_total_hists.Add(tau32_dijets_manual_2groups_improved_total);
    tau32_dijets_manual_2groups_bigcone_total_hists.Add(tau32_dijets_manual_2groups_bigcone_total);
    tau32_dijets_antikt_total_hists.Add(tau32_dijets_antikt_total);
    tau32_dijets_antikt_wta_total_hists.Add(tau32_dijets_antikt_wta_total);
  }

  int total_total_ttbar_jets = 0;
  int total_total_dijets_jets = 0;

  for (int i_pt = 0; i_pt < n_perps; i_pt++) {

    // double pt = 500;
    // double pt_2 = 600;
    double pt = (i_pt + 2)*100; 
    double pt_2 = (i_pt + 3)*100; 
    double pt_avg = (pt + pt_2)*0.5;

    ptlist[i_pt] = pt;

    ostringstream ss_pt;
    ostringstream ss_pt2;
    ss_pt << pt;
    ss_pt2 << pt_2;

    string filetitle_ttbar = "/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0" + ss_pt.str() + "-0" + ss_pt2.str() + ".UW";
    string filetitle_dijets = "/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0" + ss_pt.str() + "-0" + ss_pt2.str() + ".UW";

    int total_ttbar_jets = 0;
    int total_dijets_jets = 0;

    vector<int> total_ttbar_events;
    vector<int> total_ttbar_Wmasscut_events;
    vector<int> total_dijets_events;
    vector<int> total_dijets_Wmasscut_events;
    // int total_ttbar_events = 0;
    // int total_ttbar_Wmasscut_events = 0;
    // int total_dijets_events = 0;
    // int total_dijets_Wmasscut_events = 0;

    TObjArray area_ttbar_manual_hists;
    TObjArray area_dijets_manual_hists;

    TObjArray rawmass_ttbar_manual_2jets_hists;
    TObjArray rawmass_dijets_manual_2jets_hists;

    TObjArray rawmass_ttbar_manual_hists;
    TObjArray rawmass_dijets_manual_hists;

    TObjArray mass_ttbar_manual_hists;
    TObjArray mass_dijets_manual_hists;

    TObjArray mass_ttbar_manual_2groups_hists;
    TObjArray mass_dijets_manual_2groups_hists;

    TObjArray mass_ttbar_manual_2groups_improved_hists;
    TObjArray mass_dijets_manual_2groups_improved_hists;

    TObjArray minmass_ttbar_manual_2groups_hists;
    TObjArray minmass_dijets_manual_2groups_hists;

    TObjArray minmass_ttbar_manual_2groups_improved_hists;
    TObjArray minmass_dijets_manual_2groups_improved_hists;

    TObjArray mass_ttbar_manual_2groups_Wmasscut_hists;
    TObjArray mass_dijets_manual_2groups_Wmasscut_hists;

    TObjArray mass_ttbar_manual_2groups_improved_Wmasscut_hists;
    TObjArray mass_dijets_manual_2groups_improved_Wmasscut_hists;

    TObjArray mass_ttbar_manual_3jet4jet_2groups_compare_hists;
    TObjArray mass_dijets_manual_3jet4jet_2groups_compare_hists;

    TObjArray minmass_ttbar_manual_2jet3jet_2groups_compare_hists;
    TObjArray minmass_dijets_manual_2jet3jet_2groups_compare_hists;

    TObjArray mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare_hists;
    TObjArray mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare_hists;

    TObjArray tau32_ttbar_manual_6jet_hists;
    TObjArray tau32_dijets_manual_6jet_hists;

    TObjArray tau32_ttbar_manual_hists;
    TObjArray tau32_dijets_manual_hists;

    TObjArray tau32_ttbar_manual_2groups_hists;
    TObjArray tau32_dijets_manual_2groups_hists;

    TObjArray tau32_ttbar_manual_2groups_improved_hists;
    TObjArray tau32_dijets_manual_2groups_improved_hists;

    TObjArray tau32_ttbar_manual_2groups_bigcone_hists;
    TObjArray tau32_dijets_manual_2groups_bigcone_hists;

    TObjArray max_njet_ttbar_manual_2groups_improved_hists;
    TObjArray max_njet_dijets_manual_2groups_improved_hists;

    TObjArray max_njet_ttbar_manual_hists;
    TObjArray max_njet_dijets_manual_hists;

    TObjArray tau32_ttbar_antikt_hists;
    TObjArray tau32_dijets_antikt_hists;

    TObjArray volatility_ttbar_manual_hists;
    TObjArray volatility_dijets_manual_hists;

    for (int B = 0; B < n_betas; B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      total_ttbar_events.push_back(0);
      total_ttbar_Wmasscut_events.push_back(0);
      total_dijets_events.push_back(0);
      total_dijets_Wmasscut_events.push_back(0);

      TH1* area_ttbar_manual = new TH1F("area_ttbar_manual", "Area of top jets", 50, 0, 1.5);
      TH1* area_dijets_manual = new TH1F("area_dijets_manual", "Area of QCD jet", 50, 0, 1.5);
      
      TH1* rawmass_ttbar_manual_2jets = new TH1F("rawmass_ttbar_manual_2jets", "Raw Mass of reconstructed tops", 50, 0, 500);
      TH1* rawmass_dijets_manual_2jets = new TH1F("rawmass_dijets_manual_2jets", "Mass of QCD jet", 50, 0, 500);

      TH1* rawmass_ttbar_manual = new TH1F("rawmass_ttbar_manual", "Raw Mass of reconstructed tops", 50, 0, 500);
      TH1* rawmass_dijets_manual = new TH1F("rawmass_dijets_manual", "Mass of QCD jet", 50, 0, 500);
      TH1* mass_ttbar_manual_2groups = new TH1F("mass_ttbar_manual_2groups", "Mass of reconstructed tops from 2 groups", 50, 0, 500);
      TH1* mass_dijets_manual_2groups = new TH1F("mass_dijets_manual_2groups", "Mass of QCD jet from 2 groups", 50, 0, 500);
      TH1* mass_ttbar_manual_2groups_improved = new TH1F("mass_ttbar_manual_2groups_improved", "Mass of reconstructed tops from 2 groups", 50, 0, 500);
      TH1* mass_dijets_manual_2groups_improved = new TH1F("mass_dijets_manual_2groups_improved", "Mass of QCD jet from 2 groups", 50, 0, 500);
      TH1* mass_ttbar_manual = new TH1F("mass_ttbar_manual", "Mass of reconstructed tops from 3+3 grouping", 50, 0, 500);
      TH1* mass_dijets_manual = new TH1F("mass_dijets_manual", "Mass of QCD jet from 3+3 grouping", 50, 0, 500);
      TH1* minmass_ttbar_manual_2groups = new TH1F("minmass_ttbar_manual_2groups", "Mass of reconstructed tops from 2 groups", 50, 0, 200);
      TH1* minmass_dijets_manual_2groups = new TH1F("minmass_dijets_manual_2groups", "Mass of QCD jet from 2 groups", 50, 0, 200);
      TH1* minmass_ttbar_manual_2groups_improved = new TH1F("minmass_ttbar_manual_2groups_improved", "Mass of reconstructed tops from 2 groups", 50, 0, 200);
      TH1* minmass_dijets_manual_2groups_improved = new TH1F("minmass_dijets_manual_2groups_improved", "Mass of QCD jet from 2 groups", 50, 0, 200);
      TH1* mass_ttbar_manual_2groups_Wmasscut = new TH1F("mass_ttbar_manual_2groups_Wmasscut", "", 50, 0, 500);
      TH1* mass_dijets_manual_2groups_Wmasscut = new TH1F("mass_dijets_manual_2groups_Wmasscut", "", 50, 0, 500);
      TH1* mass_ttbar_manual_2groups_improved_Wmasscut = new TH1F("mass_ttbar_manual_2groups_improved_Wmasscut", "", 50, 0, 500);
      TH1* mass_dijets_manual_2groups_improved_Wmasscut = new TH1F("mass_dijets_manual_2groups_improved_Wmasscut", "", 50, 0, 500);

      TH2* mass_ttbar_manual_3jet4jet_2groups_compare = new TH2F("mass_ttbar_manual_3jet4jet_2groups_compare", "", 50, 0, 500, 50, 0, 500);
      TH2* mass_dijets_manual_3jet4jet_2groups_compare = new TH2F("mass_dijets_manual_3jet4jet_2groups_compare", "", 50, 0, 500, 50, 0, 500);
      TH2* minmass_ttbar_manual_2jet3jet_2groups_compare = new TH2F("minmass_ttbar_manual_2jet3jet_2groups_compare", "", 50, 0, 200, 50, 0, 200);
      TH2* minmass_dijets_manual_2jet3jet_2groups_compare = new TH2F("minmass_dijets_manual_2jet3jet_2groups_compare", "", 50, 0, 200, 50, 0, 200);
      TH2* mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare = new TH2F("mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare", "", 50, 0, 500, 50, 0, 500);
      TH2* mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare = new TH2F("mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare", "", 50, 0, 500, 50, 0, 500);

      TH1* tau32_ttbar_manual_6jet = new TH1F("tau32_ttbar_manual_6jet", "#tau_{32} of Tops", 25, 0, 1);
      TH1* tau32_dijets_manual_6jet = new TH1F("tau32_dijets_manual_6jet", "#tau_{32} of QCD", 25, 0, 1);      
      TH1* tau32_ttbar_manual = new TH1F("tau32_ttbar_manual", "#tau_{32} of Tops", 25, 0, 1);
      TH1* tau32_dijets_manual = new TH1F("tau32_dijets_manual", "#tau_{32} of QCD", 25, 0, 1);
      TH1* tau32_ttbar_manual_2groups = new TH1F("tau32_ttbar_manual_2groups", "#tau_{32} of Tops", 25, 0, 1);
      TH1* tau32_dijets_manual_2groups = new TH1F("tau32_dijets_manual_2groups", "#tau_{32} of QCD", 25, 0, 1);
      TH1* tau32_ttbar_manual_2groups_improved = new TH1F("tau32_ttbar_manual_2groups_improved", "#tau_{32} of Tops", 25, 0, 1);
      TH1* tau32_dijets_manual_2groups_improved = new TH1F("tau32_dijets_manual_2groups_improved", "#tau_{32} of QCD", 25, 0, 1);
      TH1* tau32_ttbar_manual_2groups_bigcone = new TH1F("tau32_ttbar_manual_2groups_bigcone", "", 25, 0, 1);
      TH1* tau32_dijets_manual_2groups_bigcone = new TH1F("tau32_dijets_manual_2groups_bigcone", "", 25, 0, 1);

      TH1* max_njet_ttbar_manual_2groups_improved = new TH1F("max_njet_ttbar_manual_2groups_improved", "Maximum N-jettiness of Tops from 2 groups", 10, 0, 10);
      TH1* max_njet_dijets_manual_2groups_improved = new TH1F("max_njet_dijets_manual_2groups_improved", "Maximum N-jettiness of QCD from 2 groups", 10, 0, 10);
      TH1* max_njet_ttbar_manual = new TH1F("max_njet_ttbar_manual", "Maximum N-jettiness of Tops", 20, 0, 20);
      TH1* max_njet_dijets_manual = new TH1F("max_njet_dijets_manual", "Maximum N-jettiness of QCD", 20, 0, 20);

      TH1F *tau32_ttbar_antikt = new TH1F("tau32_ttbar_antikt", "#tau_{32} for ttbar (akt)", 25, 0, 1);
      TH1F *tau32_dijets_antikt = new TH1F("tau32_dijets_antikt", "#tau_{32} for dijets (akt)", 25, 0, 1);

      TH1F* volatility_ttbar_manual = new TH1F("volatility_ttbar_manual", "", 25, 0, 1);
      TH1F* volatility_dijets_manual = new TH1F("volatility_dijets_manual", "", 25, 0, 1);

      area_ttbar_manual_hists.Add(area_ttbar_manual);
      area_dijets_manual_hists.Add(area_dijets_manual);
      
      rawmass_ttbar_manual_2jets_hists.Add(rawmass_ttbar_manual_2jets);
      rawmass_dijets_manual_2jets_hists.Add(rawmass_dijets_manual_2jets);
      rawmass_ttbar_manual_hists.Add(rawmass_ttbar_manual);
      rawmass_dijets_manual_hists.Add(rawmass_dijets_manual);
      
      mass_ttbar_manual_2groups_hists.Add(mass_ttbar_manual_2groups);
      mass_dijets_manual_2groups_hists.Add(mass_dijets_manual_2groups);
      mass_ttbar_manual_2groups_improved_hists.Add(mass_ttbar_manual_2groups_improved);
      mass_dijets_manual_2groups_improved_hists.Add(mass_dijets_manual_2groups_improved);
      mass_ttbar_manual_hists.Add(mass_ttbar_manual);
      mass_dijets_manual_hists.Add(mass_dijets_manual);
      minmass_ttbar_manual_2groups_hists.Add(minmass_ttbar_manual_2groups);
      minmass_dijets_manual_2groups_hists.Add(minmass_dijets_manual_2groups);
      minmass_ttbar_manual_2groups_improved_hists.Add(minmass_ttbar_manual_2groups_improved);
      minmass_dijets_manual_2groups_improved_hists.Add(minmass_dijets_manual_2groups_improved);
      mass_ttbar_manual_2groups_Wmasscut_hists.Add(mass_ttbar_manual_2groups_Wmasscut);
      mass_dijets_manual_2groups_Wmasscut_hists.Add(mass_dijets_manual_2groups_Wmasscut);
      mass_ttbar_manual_2groups_improved_Wmasscut_hists.Add(mass_ttbar_manual_2groups_improved_Wmasscut);
      mass_dijets_manual_2groups_improved_Wmasscut_hists.Add(mass_dijets_manual_2groups_improved_Wmasscut);

      mass_ttbar_manual_3jet4jet_2groups_compare_hists.Add(mass_ttbar_manual_3jet4jet_2groups_compare);
      mass_dijets_manual_3jet4jet_2groups_compare_hists.Add(mass_dijets_manual_3jet4jet_2groups_compare);

      minmass_ttbar_manual_2jet3jet_2groups_compare_hists.Add(minmass_ttbar_manual_2jet3jet_2groups_compare);
      minmass_dijets_manual_2jet3jet_2groups_compare_hists.Add(minmass_dijets_manual_2jet3jet_2groups_compare);
      mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare_hists.Add(mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare);
      mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare_hists.Add(mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare);

      tau32_ttbar_manual_6jet_hists.Add(tau32_ttbar_manual_6jet);
      tau32_dijets_manual_6jet_hists.Add(tau32_dijets_manual_6jet);
      tau32_ttbar_manual_hists.Add(tau32_ttbar_manual);
      tau32_dijets_manual_hists.Add(tau32_dijets_manual);
      tau32_ttbar_manual_2groups_hists.Add(tau32_ttbar_manual_2groups);
      tau32_dijets_manual_2groups_hists.Add(tau32_dijets_manual_2groups);
      tau32_ttbar_manual_2groups_improved_hists.Add(tau32_ttbar_manual_2groups_improved);
      tau32_dijets_manual_2groups_improved_hists.Add(tau32_dijets_manual_2groups_improved);
      tau32_ttbar_manual_2groups_bigcone_hists.Add(tau32_ttbar_manual_2groups_bigcone);
      tau32_dijets_manual_2groups_bigcone_hists.Add(tau32_dijets_manual_2groups_bigcone);
      
      max_njet_ttbar_manual_2groups_improved_hists.Add(max_njet_ttbar_manual_2groups_improved);
      max_njet_dijets_manual_2groups_improved_hists.Add(max_njet_dijets_manual_2groups_improved);
      max_njet_ttbar_manual_hists.Add(max_njet_ttbar_manual);
      max_njet_dijets_manual_hists.Add(max_njet_dijets_manual);

      tau32_ttbar_antikt_hists.Add(tau32_ttbar_antikt);
      tau32_dijets_antikt_hists.Add(tau32_dijets_antikt);

      volatility_ttbar_manual_hists.Add(volatility_ttbar_manual);
      volatility_dijets_manual_hists.Add(volatility_dijets_manual);
    } 


    total_ttbar_events.push_back(0);
    total_ttbar_Wmasscut_events.push_back(0);
    total_dijets_events.push_back(0);
    total_dijets_Wmasscut_events.push_back(0);

    total_ttbar_events.push_back(0);
    total_ttbar_Wmasscut_events.push_back(0);
    total_dijets_events.push_back(0);
    total_dijets_Wmasscut_events.push_back(0);

    // TH1F *tau6_ttbar_manual = new TH1F("tau6_ttbar_manual", "#tau_{6} for ttbar (Manual)", 240, 0, 1000);
    // TH1F *tau6_dijets_manual = new TH1F("tau6_dijets_manual", "#tau_{6} for dijets (Manual)", 240, 0, 1000);

    // TH1F *tau3_ttbar_manual_antikt = new TH1F("tau3_ttbar_manual_antikt", "#tau_{3} for ttbar (Manual)", 240, 0, 1000);
    // TH1F *tau3_dijets_manual_antikt = new TH1F("tau3_dijets_manual_antikt", "#tau_{3} for dijets (Manual)", 240, 0, 1000);

    TH1F* area_ttbar_antikt = new TH1F("area_ttbar_antikt", "area ttbar", 50, 0, 1.5);
    TH1F* area_dijets_antikt = new TH1F("area_dijets_antikt", "area dijets", 50, 0, 1.5);

    TH1F *mass_ttbar_antikt_6jets = new TH1F("mass_ttbar_antikt_6jets", "mass ttbar", 50, 0, 500);
    TH1F *mass_dijets_antikt_6jets = new TH1F("mass_dijets_antikt_6jets", "mass dijets", 50, 0, 500);

    TH1F *mass_ttbar_antikt_6jets_improved = new TH1F("mass_ttbar_antikt_6jets", "mass ttbar", 50, 0, 500);
    TH1F *mass_dijets_antikt_6jets_improved = new TH1F("mass_dijets_antikt_6jets", "mass dijets", 50, 0, 500);

    TH1F *mass_ttbar_antikt = new TH1F("mass_ttbar_antikt", "mass ttbar", 50, 0, 500);
    TH1F *mass_dijets_antikt = new TH1F("mass_dijets_antikt", "mass dijets", 50, 0, 500);

    TH1F *mass_ttbar_antikt_2groups = new TH1F("mass_ttbar_antikt_2groups", "mass ttbar", 50, 0, 500);
    TH1F *mass_dijets_antikt_2groups = new TH1F("mass_dijets_antikt_2groups", "mass dijets", 50, 0, 500);

    TH1F *mass_ttbar_antikt_2groups_improved = new TH1F("mass_ttbar_antikt_2groups_improved", "mass ttbar", 50, 0, 500);
    TH1F *mass_dijets_antikt_2groups_improved = new TH1F("mass_dijets_antikt_2groups_improved", "mass dijets", 50, 0, 500);

    TH1F *minmass_ttbar_antikt_2groups = new TH1F("minmass_ttbar_antikt_2groups", "", 50, 0, 200);
    TH1F *minmass_dijets_antikt_2groups = new TH1F("minmass_dijets_antikt_2groups", "", 50, 0, 200);

    TH1F* mass_ttbar_antikt_2groups_Wmasscut = new TH1F("mass_ttbar_antikt_2groups_Wmasscut", "", 50, 0, 500);
    TH1F* mass_dijets_antikt_2groups_Wmasscut = new TH1F("mass_dijets_antikt_2groups_Wmasscut", "", 50, 0, 500);

    TH1F *minmass_ttbar_antikt_2groups_improved = new TH1F("minmass_ttbar_antikt_2groups_improved", "", 50, 0, 200);
    TH1F *minmass_dijets_antikt_2groups_improved = new TH1F("minmass_dijets_antikt_2groups_improved", "", 50, 0, 200);

    TH1F* mass_ttbar_antikt_2groups_improved_Wmasscut = new TH1F("mass_ttbar_antikt_2groups_improved_Wmasscut", "", 50, 0, 500);
    TH1F* mass_dijets_antikt_2groups_improved_Wmasscut = new TH1F("mass_dijets_antikt_2groups_improved_Wmasscut", "", 50, 0, 500);

    TH2* mass_ttbar_antikt_3jet4jet_2groups_compare = new TH2F("mass_ttbar_antikt_3jet4jet_2groups_compare", "", 50, 0, 500, 50, 0, 500);
    TH2* mass_dijets_antikt_3jet4jet_2groups_compare = new TH2F("mass_dijets_antikt_3jet4jet_2groups_compare", "", 50, 0, 500, 50, 0, 500);

    TH2* minmass_ttbar_antikt_2jet3jet_2groups_compare = new TH2F("minmass_ttbar_antikt_2jet3jet_2groups_compare", "", 50, 0, 200, 50, 0, 200);
    TH2* minmass_dijets_antikt_2jet3jet_2groups_compare = new TH2F("minmass_dijets_antikt_2jet3jet_2groups_compare", "", 50, 0, 200, 50, 0, 200);

    TH2* mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare = new TH2F("mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare", "", 50, 0, 500, 50, 0, 500);
    TH2* mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare = new TH2F("mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare", "", 50, 0, 500, 50, 0, 500);

    TH1F* tau32_ttbar_antikt_wta = new TH1F("tau32_ttbar_antikt_wta", "tau32 for tops", 25, 0, 1);
    TH1F* tau32_dijets_antikt_wta = new TH1F("tau32_dijets_antikt_wta", "tau32 for QCD", 25, 0, 1);

    string samples[2] = {filetitle_ttbar, filetitle_dijets};
    // std::string samples[2] = {"/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW", "/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW"};
    // std::string samples[2] = {"/home/tjwilk/Physics/thaler_urop/BOOST_samples/pythia64-tuneDW-lhc7-ttbar2hadrons-pt0500-0600.UW", "/home/tjwilk/Physics/thaler_urop/BOOST_samples/pythia64-tuneDW-lhc7-dijets-pt0500-0600.UW"};
    const int nEvent = 100; // # of events to be analyzed   

    for (int i_sample = 0; i_sample < 2; i_sample++) {

      const char* current_data = samples[i_sample].c_str();
      ifstream inputStream(current_data); // Input File Name

    ////////// Set up Input //////////////////////////////////////////////////////////////////////////////
      vector<fastjet::PseudoJet> input_particles, partons;
      int particle_number, particle_code;
      double px, py, pz, E;
      int i_event = 1; bool is_first_event = true;

    ////////// Read in loop //////////////////////////////////////////////////////////////////////////////
      while (inputStream >> particle_number >> px >> py >> pz >> E >> particle_code){

        if (i_event > nEvent) {break;}
        if (is_first_event) {is_first_event = false; continue;}
        if (particle_number == 0) {

        // cout << i_event << endl;

          int max_njettiness;

          TH2F *event_display;
          if (i_sample == 0) event_display = new TH2F("event_display", "TTbar event (500 < pT < 600)", 60, -5, 5, 60, 0, 6.4);
          if (i_sample == 1) event_display = new TH2F("event_display", "QCD event (500 < pT < 600)", 60, -5, 5, 60, 0, 6.4);
          TH2F *axes_manual_display = new TH2F("axes_manual_display", "Axes Plot", 300, -5, 5, 300, 0, 6.4);
          TH2F *axes_manual_display_2 = new TH2F("axes_manual_display_2", "Axes Plot", 300, -5, 5, 300, 0, 6.4);
          TH2F *axes_manual_display_3 = new TH2F("axes_manual_display_3", "Axes Plot", 300, -5, 5, 300, 0, 6.4);
          TH2F *axes_manual_display_4 = new TH2F("axes_manual_display_4", "Axes Plot", 300, -5, 5, 300, 0, 6.4);
          TH2F *subjet_display = new TH2F("subjet_display", "Subjet display", 1000, -5, 5, 1000, 0, 6.4);
          TH2F *subjet_display_2 = new TH2F("subjet_display_2", "Subjet display", 1000, -5, 5, 1000, 0, 6.4);
          TH2F *subjet_display_3 = new TH2F("subjet_display_3", "Subjet display", 1000, -5, 5, 1000, 0, 6.4);
          TH2F *subjet_display_4 = new TH2F("subjet_display_4", "Subjet display", 1000, -5, 5, 1000, 0, 6.4);

          vector<PseudoJet> jet_constituents = input_particles;

          for (int j = 0; j < jet_constituents.size(); j++) {
            event_display->Fill(jet_constituents[j].eta(), jet_constituents[j].phi(), jet_constituents[j].perp()*100);
          }

          vector <PseudoJet> inclusiveJets_akt_big, sortedJets_akt_big, hardestJets_akt_big;
          vector <PseudoJet> inclusiveJets_akt, sortedJets_akt, hardestJets_akt;
          AreaDefinition area_def(active_area, GhostedAreaSpec(ghost_maxrap));

        // find 6 jets from anti-KT, group them, and fill the area and mass histogram

          Strategy strategy_akt = Best;
          RecombinationScheme recombScheme_akt = E_scheme;

          double grouping_radius = 2*Rparam;
        // double grouping_radius = 1.5;
        // double grouping_radius = 1.0; 
          double min_group_radius = 0.0;

          JetDefinition *jetDef_akt = new JetDefinition(antikt_algorithm, Rparam, recombScheme_akt, strategy_akt);
          ClusterSequenceArea clustSeq_akt(input_particles, *jetDef_akt, area_def);
          inclusiveJets_akt = clustSeq_akt.inclusive_jets(0.0);
          sortedJets_akt = sorted_by_pt(inclusiveJets_akt);
          Selector jet_selector = SelectorNHardest(6);
          hardestJets_akt = jet_selector(sortedJets_akt);

          for (int i_jet = 0; i_jet < hardestJets_akt.size(); i_jet++) {
            if (i_sample == 0) area_ttbar_antikt->Fill(hardestJets_akt[i_jet].area());
            if (i_sample == 1) area_dijets_antikt->Fill(hardestJets_akt[i_jet].area());
          }

          groupedJetsinfo bigjets_akt_info = groupJets(hardestJets_akt, grouping_radius, min_group_radius);
          vector<PseudoJet> bigjets_akt = bigjets_akt_info.groupedJets;
          Selector bigjet_selector = SelectorNHardest(2);
          vector<PseudoJet> hardest_bigjets_akt = bigjet_selector(bigjets_akt);

          for (int i_bigjet = 0; i_bigjet < hardest_bigjets_akt.size(); i_bigjet++) {
            if (i_sample == 0) mass_ttbar_antikt_6jets->Fill(hardest_bigjets_akt[i_bigjet].m());
            if (i_sample == 1) mass_dijets_antikt_6jets->Fill(hardest_bigjets_akt[i_bigjet].m());
          }

          for (unsigned int B = 0; B < betalist.size(); B++) {

            bool will_display = false;

            // TH1F *volatility_event = new TH1F("volatility_event", "", 50, 0, 500);

            double beta = betalist[B];
            double power = (double)1/beta;
          // double delta;

            const JetDefinition::Recombiner *recombScheme;
            if (beta > 1) recombScheme = new GeneralERecombiner((double)1/(beta - 1));
            else recombScheme = new WinnerTakeAllRecombiner();

            UnnormalizedCutoffMeasure measure_function = UnnormalizedCutoffMeasure(beta, Rparam);
            OnePass_Manual_Axes axes_finder = OnePass_Manual_Axes();

            vector<PseudoJet> final_bigjets;
            int njettiness = 6;

            Strategy strategy = Best;
            JetDefinition *jetDef = new JetDefinition(genkt_algorithm, Rparam, power, recombScheme, strategy);

            ClusterSequence clustSeq(jet_constituents, *jetDef);

            UnnormalizedCutoffMeasure measure_function_2jets = UnnormalizedCutoffMeasure(beta, grouping_radius);
            NjettinessPlugin njet_plugin_genrecomb_2jets(2, axes_finder, measure_function_2jets);
            JetDefinition njet_def_genrecomb_2jets(&njet_plugin_genrecomb_2jets);
            njet_plugin_genrecomb_2jets.setAxes(clustSeq.exclusive_jets(2));
            ClusterSequence njet_cluster_genrecomb_2jets(jet_constituents, njet_def_genrecomb_2jets);
            const NjettinessExtras *extras_genrecomb_2jets = njettiness_extras(njet_cluster_genrecomb_2jets);
            vector<PseudoJet> jets_genrecomb_2jets = extras_genrecomb_2jets->jets();

            for (int i_jet = 0; i_jet < jets_genrecomb_2jets.size(); i_jet++) {

              if (jets_genrecomb_2jets[i_jet].constituents().size() > 2) {

                ClusterSequence clustSeq_subjet(jets_genrecomb_2jets[i_jet].constituents(), *jetDef);
                vector<PseudoJet> subjet_axes3 = clustSeq_subjet.exclusive_jets(3);
                vector<PseudoJet> subjet_axes2 = clustSeq_subjet.exclusive_jets(2);

                Nsubjettiness nsub3_manual(3, axes_finder, UnnormalizedMeasure(beta));
                Nsubjettiness nsub2_manual(2, axes_finder, UnnormalizedMeasure(beta));

                nsub3_manual.setAxes(subjet_axes3);
                nsub2_manual.setAxes(subjet_axes2);

                if (jets_genrecomb_2jets[i_jet].m() > 160 && jets_genrecomb_2jets[i_jet].m() < 240) {
                  if (i_sample == 0) {
                    TH1* tau32_ttbar_manual_6jet = (TH1*)tau32_ttbar_manual_6jet_hists.At(B);
                    tau32_ttbar_manual_6jet->Fill((double)nsub3_manual.result(jets_genrecomb_2jets[i_jet])/nsub2_manual.result(jets_genrecomb_2jets[i_jet]));
                  }
                  if (i_sample == 1) {
                    TH1* tau32_dijets_manual_6jet = (TH1*)tau32_dijets_manual_6jet_hists.At(B);
                    tau32_dijets_manual_6jet->Fill((double)nsub3_manual.result(jets_genrecomb_2jets[i_jet])/nsub2_manual.result(jets_genrecomb_2jets[i_jet]));
                  }
                }
              }
    
              TH1* rawmass_ttbar_manual_2jets = (TH1*)rawmass_ttbar_manual_2jets_hists.At(B);
              TH1* rawmass_dijets_manual_2jets = (TH1*)rawmass_dijets_manual_2jets_hists.At(B);

              if (i_sample == 0) rawmass_ttbar_manual_2jets->Fill(jets_genrecomb_2jets[i_jet].m());
              if (i_sample == 1) rawmass_dijets_manual_2jets->Fill(jets_genrecomb_2jets[i_jet].m());
            }

            int initial_axes_num = njettiness;

            vector<PseudoJet> exclusiveJets = clustSeq.exclusive_jets(initial_axes_num);

            NjettinessPlugin njet_plugin_genrecomb(initial_axes_num, axes_finder, measure_function);
            JetDefinition njet_def_genrecomb(&njet_plugin_genrecomb);
            njet_plugin_genrecomb.setAxes(exclusiveJets);
            ClusterSequence njet_cluster_genrecomb(jet_constituents, njet_def_genrecomb);
            const NjettinessExtras *extras_genrecomb = njettiness_extras(njet_cluster_genrecomb);
            vector<PseudoJet> axes_genrecomb = extras_genrecomb->axes();
            vector<PseudoJet> jets_genrecomb = extras_genrecomb->jets();

            vector<PseudoJet> min_manual_axes = findMinAxes(jet_constituents, axes_genrecomb, njettiness, beta, Rparam);
          // vector<PseudoJet> min_manual_axes = axes_genrecomb;

            NjettinessPlugin njet_plugin_min_manual(njettiness, axes_finder, measure_function);
            JetDefinition njet_def_min_manual(&njet_plugin_min_manual);
            njet_plugin_min_manual.setAxes(min_manual_axes);

            ClusterSequenceArea njet_cluster_min_manual_area(jet_constituents, njet_def_min_manual, area_def);
            ClusterSequence njet_cluster_min_manual(jet_constituents, njet_def_min_manual);
            const NjettinessExtras *extras_min_manual = njettiness_extras(njet_cluster_min_manual);
            vector<PseudoJet> min_manual_jets = njet_cluster_min_manual_area.inclusive_jets();
          // vector<PseudoJet> min_manual_jets = extras_min_manual->jets();

          // double grouping_radius = 1.5;
            groupedJetsinfo jet_information = groupJets(min_manual_jets, grouping_radius, min_group_radius);
            int triplet_counter = jet_information.triplet_counter;
            vector<PseudoJet> massive_jets = jet_information.groupedJets;
            vector<int> massive_subjet_counter = jet_information.numGroups;

            for (int i_jet = 0; i_jet < min_manual_jets.size(); i_jet++) {
              double area_jet = min_manual_jets[i_jet].area();
              if (i_sample == 0) {
                TH1* area_ttbar_manual = (TH1*)area_ttbar_manual_hists.At(B);
                area_ttbar_manual->Fill(area_jet);
              }
              if (i_sample == 1) {
                TH1* area_dijets_manual = (TH1*)area_dijets_manual_hists.At(B);
                area_dijets_manual->Fill(area_jet);
              }
            }

            Selector hard_selector = SelectorNHardest(2);
            final_bigjets = hard_selector(massive_jets);
            for (int i_jet = 0; i_jet < final_bigjets.size(); i_jet++) {
            // if (final_bigjets[i_jet].perp() > 200) {
              // volatility_event->Fill(final_bigjets[i_jet].m());
              if (i_sample == 0) {
                TH1* rawmass_ttbar_manual = (TH1*)rawmass_ttbar_manual_hists.At(B);
                rawmass_ttbar_manual->Fill(final_bigjets[i_jet].m());
              }
              if (i_sample == 1) {
                TH1* rawmass_dijets_manual = (TH1*)rawmass_dijets_manual_hists.At(B);
                rawmass_dijets_manual->Fill(final_bigjets[i_jet].m());
              }
            // }
            }

          // METHOD #2: 2X3 Jettiness

            // vector<PseudoJet> final_manual_axes_2;
            // vector<PseudoJet> final_jets_2_improved_constituents;

            // double initial_bigradius = 1.8;
            UnnormalizedMeasure measure_function_infinite(beta);
            // UnnormalizedCutoffMeasure measure_function_infinite(beta, initial_bigradius);
            JetDefinition::Recombiner *recombScheme_big = new WinnerTakeAllRecombiner();
            JetDefinition *jetDef_big = new JetDefinition(genkt_algorithm, Rparam, power, recombScheme_big, strategy);
            ClusterSequence clustSeq_big(jet_constituents, *jetDef_big);

            NjettinessPlugin njet_plugin_2(2, axes_finder, measure_function_infinite);
            JetDefinition njet_def_2(&njet_plugin_2);
            njet_plugin_2.setAxes(clustSeq_big.exclusive_jets(2));
            ClusterSequence njet_cluster_2(jet_constituents, njet_def_2);
            const NjettinessExtras *extras_2 = njettiness_extras(njet_cluster_2);
            vector<PseudoJet> jets_2 = extras_2->jets();

            vector<PseudoJet> inclusiveJets_2, sortedJets_2, hardestJets_2;
            JetDefinition *jetDef_2 = new JetDefinition(antikt_algorithm, JetDefinition::max_allowable_R, recombScheme_akt, strategy_akt);
            ClusterSequence clustSeq_2(jet_constituents, *jetDef_2);
            inclusiveJets_2 = clustSeq_2.inclusive_jets(0.0);
            sortedJets_2 = sorted_by_pt(inclusiveJets_2);
            hardestJets_2 = jet_selector(sortedJets_2);

            for (int i_bigjet = 0; i_bigjet < jets_2.size(); i_bigjet++) {

              // vector<PseudoJet> final_bigjets_2;
              int njettiness_2 = 3;

              ClusterSequence clustSeq_2(jets_2[i_bigjet].constituents(), *jetDef);
              NjettinessPlugin njet_plugin_2times3(njettiness_2, axes_finder, measure_function);
              JetDefinition njet_def_2times3(&njet_plugin_2times3);
              njet_plugin_2times3.setAxes(clustSeq_2.exclusive_jets(njettiness_2));
              ClusterSequence njet_cluster_2times3(jets_2[i_bigjet].constituents(), njet_def_2times3);
              const NjettinessExtras *extras_2times3 = njettiness_extras(njet_cluster_2times3);
              vector<PseudoJet> jets_2times3 = extras_2times3->jets();

              PseudoJet bigjet_2(0,0,0,0);
              PseudoJet smalljet_2(0,0,0,0);
              for (int i_jet = 0; i_jet < jets_2times3.size(); i_jet++) {
                bigjet_2 = join(bigjet_2, jets_2times3[i_jet]);
              }
              vector<PseudoJet> jets_2times3_minmass = findMinMass(jets_2times3, 2);
              for (int i_jet = 0; i_jet < jets_2times3_minmass.size(); i_jet++) {
                smalljet_2 = join(smalljet_2, jets_2times3_minmass[i_jet]);
              }
              double min_mass = smalljet_2.m();

              // volatility_event->Fill(bigjet_2.m());

              TH1* mass_ttbar_manual_2groups = (TH1*)mass_ttbar_manual_2groups_hists.At(B);
              TH1* minmass_ttbar_manual_2groups = (TH1*)minmass_ttbar_manual_2groups_hists.At(B);
              TH1* mass_dijets_manual_2groups = (TH1*)mass_dijets_manual_2groups_hists.At(B);
              TH1* minmass_dijets_manual_2groups = (TH1*)minmass_dijets_manual_2groups_hists.At(B);
              TH1* mass_ttbar_manual_2groups_Wmasscut = (TH1*)mass_ttbar_manual_2groups_Wmasscut_hists.At(B);
              TH1* mass_dijets_manual_2groups_Wmasscut = (TH1*)mass_dijets_manual_2groups_Wmasscut_hists.At(B);
              
              if (i_sample == 0) {
                mass_ttbar_manual_2groups->Fill(bigjet_2.m());
                minmass_ttbar_manual_2groups->Fill(min_mass);
                if (min_mass > 50) mass_ttbar_manual_2groups_Wmasscut->Fill(bigjet_2.m());
              }
              if (i_sample == 1) {
                mass_dijets_manual_2groups->Fill(bigjet_2.m());
                minmass_dijets_manual_2groups->Fill(min_mass);
                if (min_mass > 50) mass_dijets_manual_2groups_Wmasscut->Fill(bigjet_2.m());
              }

            //now for the akt jets
              vector<PseudoJet> inclusiveJets_2times3, sortedJets_2times3, hardestJets_2times3;
              JetDefinition *jetDef_2times3 = new JetDefinition(antikt_algorithm, Rparam, recombScheme_akt, strategy_akt);
              ClusterSequence clustSeq_2times3(jet_constituents, *jetDef_2times3);
              inclusiveJets_2times3 = clustSeq_2times3.inclusive_jets(0.0);
              sortedJets_2times3 = sorted_by_pt(inclusiveJets_2times3);
              Selector threejet_selector = SelectorNHardest(njettiness_2);
              hardestJets_2times3 = threejet_selector(sortedJets_2times3);

              PseudoJet bigjet_2_akt(0,0,0,0);
              PseudoJet smalljet_2_akt(0,0,0,0);
              for (int i_jet = 0; i_jet < hardestJets_2times3.size(); i_jet++) {
                bigjet_2_akt = join(bigjet_2_akt, hardestJets_2times3[i_jet]);
              }
              vector<PseudoJet> hardestJets_2times3_minmass = findMinMass(hardestJets_2times3, 2);
              for (int i_jet = 0; i_jet < hardestJets_2times3_minmass.size(); i_jet++) {
                smalljet_2_akt = join(smalljet_2_akt, hardestJets_2times3_minmass[i_jet]);
              }
              double min_mass_akt = smalljet_2_akt.m();

              if (i_sample == 0) {
                mass_ttbar_antikt_2groups->Fill(bigjet_2_akt.m());
                minmass_ttbar_antikt_2groups->Fill(min_mass_akt);
                if (min_mass_akt > 50) mass_ttbar_antikt_2groups_Wmasscut->Fill(bigjet_2_akt.m());
              }
              if (i_sample == 1) {
                mass_dijets_antikt_2groups->Fill(bigjet_2_akt.m());
                minmass_dijets_antikt_2groups->Fill(min_mass_akt);
                if (min_mass_akt > 50) mass_dijets_antikt_2groups_Wmasscut->Fill(bigjet_2_akt.m());
              }

              if (bigjet_2.constituents().size() > 2) {

                ClusterSequence clustSeq_subjet(bigjet_2.constituents(), *jetDef);
                vector<PseudoJet> subjet_axes3 = clustSeq_subjet.exclusive_jets(3);
                vector<PseudoJet> subjet_axes2 = clustSeq_subjet.exclusive_jets(2);

                Nsubjettiness nsub3_manual(3, axes_finder, UnnormalizedMeasure(beta));
                Nsubjettiness nsub2_manual(2, axes_finder, UnnormalizedMeasure(beta));

                nsub3_manual.setAxes(subjet_axes3);
                nsub2_manual.setAxes(subjet_axes2);
                
                if (bigjet_2.m() > 160 && bigjet_2.m() < 240) {
                  if (i_sample == 0) {
                    TH1* tau32_ttbar_manual_2groups = (TH1*)tau32_ttbar_manual_2groups_hists.At(B);
                    tau32_ttbar_manual_2groups->Fill((double)nsub3_manual.result(bigjet_2)/nsub2_manual.result(bigjet_2));
                  }
                  if (i_sample == 1) {
                    TH1* tau32_dijets_manual_2groups = (TH1*)tau32_dijets_manual_2groups_hists.At(B);
                    tau32_dijets_manual_2groups->Fill((double)nsub3_manual.result(bigjet_2)/nsub2_manual.result(bigjet_2));
                  }
                }
              }


            // NOW FOR THE 4 CHOOSE 3 METHOD

              NjettinessPlugin njet_plugin_2times4(njettiness_2 + 1, axes_finder, measure_function);
              JetDefinition njet_def_2times4(&njet_plugin_2times4);
              njet_plugin_2times4.setAxes(clustSeq_2.exclusive_jets(njettiness_2 + 1));
              ClusterSequence njet_cluster_2times4(jets_2[i_bigjet].constituents(), njet_def_2times4);
              const NjettinessExtras *extras_2times4 = njettiness_extras(njet_cluster_2times4);
              vector<PseudoJet> jets_2times4 = extras_2times4->jets();

              PseudoJet bigjet_2_improved(0,0,0,0);
              PseudoJet smalljet_2_improved(0,0,0,0);
              vector<PseudoJet> jets_2times3_improved = findMinMass(jets_2times4, 3);
              // vector<PseudoJet> jets_2times3_improved = findMaxPerp(jets_2times4, 3);
              vector<PseudoJet> jets_2times3_improved_minmass = findMinMass(jets_2times3_improved, 2);
              for (int i_jet = 0; i_jet < jets_2times3_improved.size(); i_jet++) {
                bigjet_2_improved = join(bigjet_2_improved, jets_2times3_improved[i_jet]);
              }
              for (int i_jet = 0; i_jet < jets_2times3_improved_minmass.size(); i_jet++) {
                smalljet_2_improved = join(smalljet_2_improved, jets_2times3_improved_minmass[i_jet]);
              }
              double min_mass_improved = smalljet_2_improved.m();

              // volatility_event->Fill(final_bigjets_2[i_jet].m());
              TH1* mass_ttbar_manual_2groups_improved = (TH1*)mass_ttbar_manual_2groups_improved_hists.At(B);
              TH1* minmass_ttbar_manual_2groups_improved = (TH1*)minmass_ttbar_manual_2groups_improved_hists.At(B);
              TH1* mass_dijets_manual_2groups_improved = (TH1*)mass_dijets_manual_2groups_improved_hists.At(B);
              TH1* minmass_dijets_manual_2groups_improved = (TH1*)minmass_dijets_manual_2groups_improved_hists.At(B);
              TH1* mass_ttbar_manual_2groups_improved_Wmasscut = (TH1*)mass_ttbar_manual_2groups_improved_Wmasscut_hists.At(B);
              TH1* mass_dijets_manual_2groups_improved_Wmasscut = (TH1*)mass_dijets_manual_2groups_improved_Wmasscut_hists.At(B);
              
              TH2* mass_ttbar_manual_3jet4jet_2groups_compare = (TH2*)mass_ttbar_manual_3jet4jet_2groups_compare_hists.At(B);
              TH2* mass_dijets_manual_3jet4jet_2groups_compare = (TH2*)mass_dijets_manual_3jet4jet_2groups_compare_hists.At(B);
              TH2* minmass_ttbar_manual_2jet3jet_2groups_compare = (TH2*)minmass_ttbar_manual_2jet3jet_2groups_compare_hists.At(B);
              TH2* minmass_dijets_manual_2jet3jet_2groups_compare = (TH2*)minmass_dijets_manual_2jet3jet_2groups_compare_hists.At(B);
              TH2* mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare = (TH2*)mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare_hists.At(B);
              TH2* mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare = (TH2*)mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare_hists.At(B);

              if (i_sample == 0) {
                mass_ttbar_manual_2groups_improved->Fill(bigjet_2_improved.m());
                minmass_ttbar_manual_2groups_improved->Fill(min_mass_improved);
                mass_ttbar_manual_3jet4jet_2groups_compare->Fill(bigjet_2.m(), bigjet_2_improved.m());
                minmass_ttbar_manual_2jet3jet_2groups_compare->Fill(min_mass, min_mass_improved);
                if (min_mass_improved > 50) mass_ttbar_manual_2groups_improved_Wmasscut->Fill(bigjet_2_improved.m());
                if (min_mass > 50 && min_mass_improved > 50) mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare->Fill(bigjet_2.m(), bigjet_2_improved.m());

                if ((bigjet_2.m() > 160 && bigjet_2.m() < 240) || (bigjet_2_improved.m() > 160 && bigjet_2_improved.m() < 240)) {
                  total_ttbar_events[B + 2]++;
                  if (min_mass > 50 && min_mass_improved > 50) total_ttbar_Wmasscut_events[B + 2]++;
                }
              }
              if (i_sample == 1) {
                mass_dijets_manual_2groups_improved->Fill(bigjet_2_improved.m());
                minmass_dijets_manual_2groups_improved->Fill(min_mass_improved);
                mass_dijets_manual_3jet4jet_2groups_compare->Fill(bigjet_2.m(), bigjet_2_improved.m());
                minmass_dijets_manual_2jet3jet_2groups_compare->Fill(min_mass, min_mass_improved);
                if (min_mass_improved > 50) mass_dijets_manual_2groups_improved_Wmasscut->Fill(bigjet_2_improved.m());
                if (min_mass > 50 && min_mass_improved > 50) mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare->Fill(bigjet_2.m(), bigjet_2_improved.m());

                if ((bigjet_2.m() > 160 && bigjet_2.m() < 240) || (bigjet_2_improved.m() > 160 && bigjet_2_improved.m() < 240)) {
                  total_dijets_events[B + 2]++;
                  if (min_mass > 50 && min_mass_improved > 50) total_dijets_Wmasscut_events[B + 2]++;
                }
              }

            //ANTIKT 4 CHOOSE 3
              vector<PseudoJet> inclusiveJets_2times4, sortedJets_2times4, hardestJets_2times4;
              JetDefinition *jetDef_2times4 = new JetDefinition(antikt_algorithm, Rparam, recombScheme_akt, strategy_akt);
              ClusterSequence clustSeq_2times4(jet_constituents, *jetDef_2times4);
              inclusiveJets_2times4 = clustSeq_2times4.inclusive_jets(0.0);
              sortedJets_2times4 = sorted_by_pt(inclusiveJets_2times4);
              Selector fourjet_selector = SelectorNHardest(njettiness_2 + 1);
              hardestJets_2times4 = fourjet_selector(sortedJets_2times4);

              PseudoJet bigjet_2_akt_improved(0,0,0,0);
              PseudoJet smalljet_2_akt_improved(0,0,0,0);
              vector<PseudoJet> hardestJets_2times3_improved = findMinMass(hardestJets_2times4, 3);
              // vector<PseudoJet> hardestJets_2times3_improved = findMaxPerp(hardestJets_2times4, 3);
              vector<PseudoJet> hardestJets_2times3_improved_minmass = findMinMass(hardestJets_2times3_improved, 2);

              for (int i_jet = 0; i_jet < hardestJets_2times3_improved.size(); i_jet++) {
                bigjet_2_akt_improved = join(bigjet_2_akt_improved, hardestJets_2times3_improved[i_jet]);
              }
              for (int i_jet = 0; i_jet < hardestJets_2times3_improved_minmass.size(); i_jet++) {
                smalljet_2_akt_improved = join(smalljet_2_akt_improved, hardestJets_2times3_improved_minmass[i_jet]);
              }
              double min_mass_akt_improved = smalljet_2_akt_improved.m();

              if (i_sample == 0) {
                mass_ttbar_antikt_2groups_improved->Fill(bigjet_2_akt_improved.m());
                minmass_ttbar_antikt_2groups_improved->Fill(min_mass_akt_improved);
                mass_ttbar_antikt_3jet4jet_2groups_compare->Fill(bigjet_2_akt.m(), bigjet_2_akt_improved.m());
                minmass_ttbar_antikt_2jet3jet_2groups_compare->Fill(min_mass_akt, min_mass_akt_improved);
                if (min_mass_akt_improved > 50) mass_ttbar_antikt_2groups_improved_Wmasscut->Fill(bigjet_2_akt_improved.m());
                if (min_mass_akt > 50 && min_mass_akt_improved > 50) mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare->Fill(bigjet_2_akt.m(), bigjet_2_akt_improved.m());

                if ((bigjet_2_akt.m() > 160 && bigjet_2_akt.m() < 240) || (bigjet_2_akt_improved.m() > 160 && bigjet_2_akt_improved.m() < 240)) {
                  total_ttbar_events[1]++;
                  if (min_mass_akt > 50 && min_mass_akt_improved > 50) total_ttbar_Wmasscut_events[1]++;
                }
              }
              if (i_sample == 1) {
                mass_dijets_antikt_2groups_improved->Fill(bigjet_2_akt_improved.m());
                minmass_dijets_antikt_2groups_improved->Fill(min_mass_akt_improved);
                mass_dijets_antikt_3jet4jet_2groups_compare->Fill(bigjet_2_akt.m(), bigjet_2_akt_improved.m());
                minmass_dijets_antikt_2jet3jet_2groups_compare->Fill(min_mass_akt, min_mass_akt_improved);
                if (min_mass_akt_improved > 50) mass_dijets_antikt_2groups_improved_Wmasscut->Fill(bigjet_2_akt_improved.m());
                if (min_mass_akt > 50 && min_mass_akt_improved > 50) mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare->Fill(bigjet_2_akt.m(), bigjet_2_akt_improved.m());

                if ((bigjet_2_akt.m() > 160 && bigjet_2_akt.m() < 240) || (bigjet_2_akt_improved.m() > 160 && bigjet_2_akt_improved.m() < 240)) {
                  total_dijets_events[1]++;
                  if (min_mass_akt > 50 && min_mass_akt_improved > 50) total_dijets_Wmasscut_events[1]++;
                }
              }

              if (bigjet_2_improved.constituents().size() > 2) {

                ClusterSequence clustSeq_subjet(bigjet_2_improved.constituents(), *jetDef);
                vector<PseudoJet> subjet_axes3 = clustSeq_subjet.exclusive_jets(3);
                vector<PseudoJet> subjet_axes2 = clustSeq_subjet.exclusive_jets(2);

                Nsubjettiness nsub3_manual(3, axes_finder, UnnormalizedMeasure(beta));
                Nsubjettiness nsub2_manual(2, axes_finder, UnnormalizedMeasure(beta));

                nsub3_manual.setAxes(subjet_axes3);
                nsub2_manual.setAxes(subjet_axes2);

                if (bigjet_2_improved.m() > 160 && bigjet_2_improved.m() < 240) {
                  if (i_sample == 0) {
                    TH1* tau32_ttbar_manual_2groups_improved = (TH1*)tau32_ttbar_manual_2groups_improved_hists.At(B);
                    tau32_ttbar_manual_2groups_improved->Fill((double)nsub3_manual.result(bigjet_2_improved)/nsub2_manual.result(bigjet_2_improved));
                  }
                  if (i_sample == 1) {
                    TH1* tau32_dijets_manual_2groups_improved = (TH1*)tau32_dijets_manual_2groups_improved_hists.At(B);
                    tau32_dijets_manual_2groups_improved->Fill((double)nsub3_manual.result(bigjet_2_improved)/nsub2_manual.result(bigjet_2_improved));
                  }
                }
              }

              // vector<PseudoJet> final_bigjet_bigcone_constituents;
              // for (int i_part = 0; i_part < jet_constituents.size(); i_part++) {
              //   if (temp_jet_improved.delta_R(jet_constituents[i_part]) < 1.0) final_bigjet_bigcone_constituents.push_back(jet_constituents[i_part]);
              // }

              // PseudoJet final_bigjet_bigcone(0,0,0,0);
              // for (int i_const = 0; i_const < final_bigjet_bigcone_constituents.size(); i_const++) {
              //   final_bigjet_bigcone = join(final_bigjet_bigcone, final_bigjet_bigcone_constituents[i_const]);
              // }

              // ClusterSequence clustSeq_subjet(final_bigjet_bigcone.constituents(), *jetDef);
              // cout << final_bigjet_bigcone.constituents().size() << endl;
              // vector<PseudoJet> subjet_axes3 = clustSeq_subjet.exclusive_jets(3);
              // vector<PseudoJet> subjet_axes2 = clustSeq_subjet.exclusive_jets(2);

              // Nsubjettiness nsub3_manual(3, axes_finder, UnnormalizedMeasure(beta));
              // Nsubjettiness nsub2_manual(2, axes_finder, UnnormalizedMeasure(beta));

              // nsub3_manual.setAxes(subjet_axes3);
              // nsub2_manual.setAxes(subjet_axes2);

              // if (temp_jet_improved.m() > 160 && temp_jet_improved.m() < 240) {
              //   if (i_sample == 0) {
              //     TH1* tau32_ttbar_manual_2groups_bigcone = (TH1*) tau32_ttbar_manual_2groups_bigcone_hists.At(B);
              //     tau32_ttbar_manual_2groups_bigcone->Fill((double)nsub3_manual.result(final_bigjet_bigcone)/nsub2_manual.result(final_bigjet_bigcone));
              //   }
              //   if (i_sample == 1) {
              //     TH1* tau32_dijets_manual_2groups_bigcone = (TH1*) tau32_dijets_manual_2groups_bigcone_hists.At(B);
              //     tau32_dijets_manual_2groups_bigcone->Fill((double)nsub3_manual.result(final_bigjet_bigcone)/nsub2_manual.result(final_bigjet_bigcone));
              //   }
              // }
              // }

              // }
              // final_manual_axes = min_manual_axes;
              // max_njettiness_2 = njettiness_2;
              // break;
              // }
              // njettiness_2++;
              // }
            // }
            }

            // TH1F* volatility_ttbar_manual = (TH1F*) volatility_ttbar_manual_hists.At(B);
            // TH1F* volatility_dijets_manual = (TH1F*) volatility_dijets_manual_hists.At(B);

            // double stdev = volatility_event->GetRMS();
            // double mean = volatility_event->GetMean();

            // if (mean > 160 && mean < 240) { 
            //   if (i_sample == 0) volatility_ttbar_manual->Fill(stdev/mean);
            //   if (i_sample == 1) volatility_dijets_manual->Fill(stdev/mean);
            // }

            // delete volatility_event;
          // // ALL OF THIS STUFF IS SIMPLY FOR DISPLAY
          // if (abs(beta - 1) < epsilon && will_display) {
          //   Double_t ghost_perp = 0.001;
          //   Double_t n_ghosts = 50;
          //   for (int i_ghosts = 0; i_ghosts < n_ghosts; i_ghosts++) {
          //     for (int j_ghosts = 0; j_ghosts < n_ghosts; j_ghosts++) {
          //       Double_t ghost_eta = -5 + (double)(10/n_ghosts)*i_ghosts;
          //       Double_t ghost_phi = (double)(2*TMath::Pi()/n_ghosts)*j_ghosts;
          //       PseudoJet ghost(ghost_perp*TMath::Cos(ghost_phi),ghost_perp*TMath::Sin(ghost_phi),
          //         ghost_perp*TMath::SinH(ghost_eta),ghost_perp*TMath::CosH(ghost_eta));
          //       jet_constituents.push_back(ghost);
          //       // final_jets_2_improved_constituents.push_back(ghost);
          //     }
          //   } 

          //   //first method
          //   NjettinessPlugin njet_plugin_min_manual(6, axes_finder, measure_function);
          //   JetDefinition njet_def_min_manual(&njet_plugin_min_manual);
          //   njet_plugin_min_manual.setAxes(final_manual_axes);
          //   final_manual_axes.resize(0);
          //   ClusterSequence njet_cluster_min_manual(jet_constituents, njet_def_min_manual);
          //   const NjettinessExtras *extras_min_manual = njettiness_extras(njet_cluster_min_manual);
          //   vector<PseudoJet> final_manual_jets = extras_min_manual->jets();
          //   final_manual_axes = extras_min_manual->axes();

          //   //second method
          //   NjettinessPlugin njet_plugin_min_manual_2(6, axes_finder, measure_function);
          //   JetDefinition njet_def_min_manual_2(&njet_plugin_min_manual_2);
          //   njet_plugin_min_manual_2.setAxes(final_manual_axes_2);
          //   final_manual_axes_2.resize(0);
          //   ClusterSequence njet_cluster_min_manual_2(jet_constituents, njet_def_min_manual_2);
          //   const NjettinessExtras *extras_min_manual_2 = njettiness_extras(njet_cluster_min_manual_2);
          //   vector<PseudoJet> final_manual_jets_2 = extras_min_manual_2->jets();
          //   final_manual_axes_2 = extras_min_manual_2->axes();

          //   ClusterSequence clustSeq_akt(jet_constituents, *jetDef_akt);
          //   vector<PseudoJet> inclusiveJets_akt_disp = clustSeq_akt.inclusive_jets(0.0);
          //   vector<PseudoJet> sortedJets_akt_disp = sorted_by_pt(inclusiveJets_akt_disp);
          //   vector<PseudoJet> hardestJets_akt_disp = jet_selector(sortedJets_akt_disp);

          //   // //third method
          //   // NjettinessPlugin njet_plugin_min_manual(6, axes_finder, measure_function);
          //   // JetDefinition njet_def_min_manual(&njet_plugin_min_manual);
          //   // njet_plugin_min_manual.setAxes(final_manual_axes);
          //   // final_manual_axes.resize(0);
          //   // ClusterSequence njet_cluster_min_manual(jet_constituents, njet_def_min_manual);
          //   // const NjettinessExtras *extras_min_manual = njettiness_extras(njet_cluster_min_manual);
          //   // vector<PseudoJet> final_manual_jets = extras_min_manual->jets();
          //   // final_manual_axes = extras_min_manual->axes();

          //   // //fourth method
          //   // NjettinessPlugin njet_plugin_min_manual(6, axes_finder, measure_function);
          //   // JetDefinition njet_def_min_manual(&njet_plugin_min_manual);
          //   // njet_plugin_min_manual.setAxes(final_manual_axes);
          //   // final_manual_axes.resize(0);
          //   // ClusterSequence njet_cluster_min_manual(jet_constituents, njet_def_min_manual);
          //   // const NjettinessExtras *extras_min_manual = njettiness_extras(njet_cluster_min_manual);
          //   // vector<PseudoJet> final_manual_jets = extras_min_manual->jets();
          //   // final_manual_axes = extras_min_manual->axes();

          //   for (int a = 0; a < final_manual_axes.size(); a++) {
          //     axes_manual_display->Fill(final_manual_axes[a].eta(), final_manual_axes[a].phi());
          //   }

          //   for (int a = 0; a < final_manual_axes_2.size(); a++) {
          //     axes_manual_display_2->Fill(final_manual_axes_2[a].eta(), final_manual_axes_2[a].phi());
          //   }

          //   for (int i = 0; i < final_manual_jets.size(); i++) {
          //     vector<PseudoJet> constituents = final_manual_jets[i].constituents();
          //     for (int i_const = 0; i_const < constituents.size(); i_const++) {
          //       subjet_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
          //     }
          //   }

          //   for (int i = 0; i < final_manual_jets_2.size(); i++) {
          //     vector<PseudoJet> constituents = final_manual_jets_2[i].constituents();
          //     for (int i_const = 0; i_const < constituents.size(); i_const++) {
          //       subjet_display_2->Fill(constituents[i_const].eta(), constituents[i_const].phi());
          //     }
          //   }

          //   for (int i = 0; i < hardestJets_akt_disp.size(); i++) {
          //     vector<PseudoJet> constituents = hardestJets_akt_disp[i].constituents();
          //     for (int i_const = 0; i_const < constituents.size(); i_const++) {
          //       subjet_display_3->Fill(constituents[i_const].eta(), constituents[i_const].phi());
          //     }
          //   }

          //   vector<PseudoJet> constituents = final_jets_2_improved_constituents;
          //   for (int i_const = 0; i_const < constituents.size(); i_const++) {
          //     subjet_display_4->Fill(constituents[i_const].eta(), constituents[i_const].phi());
          //   }

          // }

            vector <PseudoJet> inclusiveJets_std, sortedJets_std, hardestJets_std;

            double Rparam_std = 1.0;
            Strategy strategy_std = Best;
            RecombinationScheme recombScheme_std = E_scheme;
          // JetDefinition::Recombiner *recombScheme_std = new WinnerTakeAllRecombiner();
            JetDefinition *jetDef_std = new JetDefinition(antikt_algorithm, Rparam_std, recombScheme_std, strategy_std);
            ClusterSequence clustSeq_std(input_particles, *jetDef_std);
            inclusiveJets_std = clustSeq_std.inclusive_jets(200.0);
            sortedJets_std = sorted_by_pt(inclusiveJets_std);
            Selector jet_selector = SelectorNHardest(2);
            hardestJets_std = jet_selector(sortedJets_std);

            for (int i = 0; i < hardestJets_std.size(); i++) {

              ClusterSequence clustSeq_subjet(hardestJets_std[i].constituents(), *jetDef);
            // vector<PseudoJet> subjet_axes3_start = clustSeq_subjet.exclusive_jets(4);
            // vector<PseudoJet> subjet_axes2_start = clustSeq_subjet.exclusive_jets(3);
            // vector<PseudoJet> subjet_axes3 = findMinAxes(hardestJets_std[i].constituents(), subjet_axes3_start, 3, beta, Rparam);
            // vector<PseudoJet> subjet_axes2 = findMinAxes(hardestJets_std[i].constituents(), subjet_axes2_start, 2, beta, Rparam);
              vector<PseudoJet> subjet_axes3 = clustSeq_subjet.exclusive_jets(3);
              vector<PseudoJet> subjet_axes2 = clustSeq_subjet.exclusive_jets(2);

              Nsubjettiness nsub3_manual(3, axes_finder, UnnormalizedMeasure(beta));
              Nsubjettiness nsub2_manual(2, axes_finder, UnnormalizedMeasure(beta));

              nsub3_manual.setAxes(subjet_axes3);
              nsub2_manual.setAxes(subjet_axes2);

              if (i_sample == 0) {
                total_ttbar_jets++;
                mass_ttbar_antikt->Fill(hardestJets_std[i].m());
              // tau3_ttbar_manual_antikt->Fill(nsub3_manual.result(hardestJets_std[i]));
              }
              if (i_sample == 1) {
                total_dijets_jets++;
                mass_dijets_antikt->Fill(hardestJets_std[i].m());
              // tau3_dijets_manual_antikt->Fill(nsub3_manual.result(hardestJets_std[i]));
              }

              NsubjettinessRatio nsub_32_wta(3, 2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(1.0));
              if (hardestJets_std[i].m() > 160 && hardestJets_std[i].m() < 240) {
                if (i_sample == 0) {
                  TH1* tau32_ttbar_antikt = (TH1*)tau32_ttbar_antikt_hists.At(B);
                  tau32_ttbar_antikt->Fill(nsub3_manual.result(hardestJets_std[i])/nsub2_manual.result(hardestJets_std[i]));
                  tau32_ttbar_antikt_wta->Fill(nsub_32_wta.result(hardestJets_std[i]));
                }
                if (i_sample == 1) {
                  TH1* tau32_dijets_antikt = (TH1*)tau32_dijets_antikt_hists.At(B);
                  tau32_dijets_antikt->Fill(nsub3_manual.result(hardestJets_std[i])/nsub2_manual.result(hardestJets_std[i]));
                  tau32_dijets_antikt_wta->Fill(nsub_32_wta.result(hardestJets_std[i]));
                }
              }
            }
          }

          event_display->SetStats(0);
          axes_manual_display->SetStats(0);
          axes_manual_display_2->SetStats(0);
          subjet_display->SetStats(0);
          subjet_display_2->SetStats(0);
          subjet_display_3->SetStats(0);
          subjet_display_4->SetStats(0);
        // subjet5_display->SetStats(0);
        // subjet6_display->SetStats(0);
        // subjet7_display->SetStats(0);

          TCanvas *display = new TCanvas("display", "Event Display", 1093, 700);
          display->cd();
          display->SetFixedAspectRatio();
          event_display->GetXaxis()->SetTitle("#eta");
          event_display->GetYaxis()->SetTitle("#phi");
          event_display->SetFillColor(kBlack);
          event_display->SetLineColor(kBlack);
        // event_display->SetLineWidth(1);
          subjet_display->SetMarkerColor(kBlue);
          subjet_display_2->SetMarkerColor(kRed);
          subjet_display_3->SetMarkerColor(kOrange);
          subjet_display_4->SetMarkerColor(8);
        // subjet5_display->SetMarkerColor(kBlue);
        // subjet6_display->SetMarkerColor(6);
        // subjet7_display->SetMarkerColor(20);
          subjet_display->SetMarkerStyle(21);
          subjet_display_2->SetMarkerStyle(21);
          subjet_display_3->SetMarkerStyle(21);
          subjet_display_4->SetMarkerStyle(21);
        // subjet5_display->SetMarkerStyle(21);
        // subjet6_display->SetMarkerStyle(21);
        // subjet7_display->SetMarkerStyle(21);
          subjet_display->SetMarkerSize(0.5);
          subjet_display_2->SetMarkerSize(0.5);
          subjet_display_3->SetMarkerSize(0.5);
          subjet_display_4->SetMarkerSize(0.5);
        // subjet5_display->SetMarkerSize(0.5);
        // subjet6_display->SetMarkerSize(0.5);
        // subjet7_display->SetMarkerSize(0.5);
          event_display->Draw("box");

          axes_manual_display->SetMarkerStyle(3);
          axes_manual_display->SetMarkerSize(3);
          axes_manual_display->SetMarkerColor(kBlue);
          axes_manual_display->Draw("SAMES");

          axes_manual_display_2->SetMarkerStyle(3);
          axes_manual_display_2->SetMarkerSize(3);
          axes_manual_display_2->SetMarkerColor(kRed);
          axes_manual_display_2->Draw("SAMES");

        // axes_manual_display_3->SetMarkerStyle(3);
        // axes_manual_display_3->SetMarkerSize(3);
        // axes_manual_display_3->SetMarkerColor(kRed);
        // axes_manual_display_3->Draw("SAMES");

        // axes_manual_display_4->SetMarkerStyle(3);
        // axes_manual_display_4->SetMarkerSize(3);
        // axes_manual_display_4->SetMarkerColor(kRed);
        // axes_manual_display_4->Draw("SAMES");

          subjet_display->Draw("SAMES");
          subjet_display_2->Draw("SAMES");
          subjet_display_3->Draw("SAMES");
          subjet_display_4->Draw("SAMES");
        // subjet5_display->Draw("SAMES");
        // subjet6_display->Draw("SAMES");
        // subjet7_display->Draw("SAMES");

          if (subjet_display->Integral() != 0) display->Write();

          ostringstream ss;
          ss << i_event;
          TString title;
          if (i_sample == 0) title = "ttbar_display" + (TString)ss_pt.str() + "_6jets_test_pt" + (TString)ss_pt.str() + ".eps";
          if (i_sample == 1) title = "dijets_display" + (TString)ss_pt.str() + "_6jets_test_pt" + (TString)ss_pt.str() + ".eps";

        //   // display->Print(title, "eps");
        // }

          delete display;
          delete event_display;
          delete axes_manual_display;
          delete axes_manual_display_2;
          delete axes_manual_display_3;
          delete axes_manual_display_4;
          delete subjet_display;
          delete subjet_display_2;
          delete subjet_display_3;
          delete subjet_display_4;

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

    total_ttbar_jets = total_ttbar_jets/n_betas;
    total_dijets_jets = total_dijets_jets/n_betas;

    total_total_ttbar_jets += total_ttbar_jets;
    total_total_dijets_jets += total_dijets_jets;

    area_ttbar_antikt->SetStats(0);
    area_dijets_antikt->SetStats(0);
    mass_ttbar_antikt->SetStats(0);
    mass_dijets_antikt->SetStats(0);
    mass_ttbar_antikt_6jets->SetStats(0);
    mass_dijets_antikt_6jets->SetStats(0);

    area_ttbar_antikt->Write();
    area_dijets_antikt->Write();
    mass_ttbar_antikt->Write();
    mass_dijets_antikt->Write();
    mass_ttbar_antikt_6jets->Write();
    mass_dijets_antikt_6jets->Write();


    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];
      ostringstream ss;
      ss << beta;

      TCanvas *ROC_compare = new TCanvas("ROC_compare", "ROC Curve Comparison", 600, 600);
      ROC_compare->cd();
      ROC_compare->SetLogy();
      TMultiGraph *ROC_multigraph = new TMultiGraph("ROC_multigraph", "ROC Comparison for #tau_{3}/#tau_{2} (" + (TString)ss_pt.str() + " < p_{T} < " + (TString)ss_pt2.str() + ", 160 < m < 240, #beta = " + (TString)ss.str() + ")");
      TLegend *leg_ROC = new TLegend(0.14, 0.6, 0.4, 0.88);
      leg_ROC->SetLineColor(0);
      leg_ROC->SetFillColor(0);

      TH1* tau32_ttbar_antikt_wta_total = (TH1F*)tau32_ttbar_antikt_wta_total_hists.At(B);
      TH1* tau32_dijets_antikt_wta_total = (TH1F*)tau32_dijets_antikt_wta_total_hists.At(B);
      tau32_ttbar_antikt_wta_total->Add(tau32_ttbar_antikt_wta);
      tau32_dijets_antikt_wta_total->Add(tau32_dijets_antikt_wta);

      tau32_ttbar_antikt_wta->SetStats(0);
      tau32_dijets_antikt_wta->SetStats(0);
      tau32_ttbar_antikt_wta->Write();
      tau32_dijets_antikt_wta->Write();

      int n_antikt_wta = tau32_ttbar_antikt_wta->GetSize() - 2;
      double integral_dijets_antikt_wta[n_antikt_wta], integral_ttbar_antikt_wta[n_antikt_wta];
      for (int i = 0; i < n_antikt_wta; i++) {
        integral_ttbar_antikt_wta[i] = tau32_ttbar_antikt_wta->Integral(0,i)/(total_ttbar_jets*n_betas);
        integral_dijets_antikt_wta[i] = tau32_dijets_antikt_wta->Integral(0,i)/(total_dijets_jets*n_betas);
      }
      TGraph* ROC_tau32_antikt_wta = new TGraph(n_antikt_wta, integral_ttbar_antikt_wta, integral_dijets_antikt_wta);
      ROC_tau32_antikt_wta->GetXaxis()->SetLimits(0, 1);
      ROC_tau32_antikt_wta->GetYaxis()->SetLimits(0, 1);
      ROC_tau32_antikt_wta->SetLineColor(kBlack);
      ROC_tau32_antikt_wta->SetLineWidth(2);
      ROC_tau32_antikt_wta->SetMarkerStyle(5);
      ROC_tau32_antikt_wta->SetMarkerSize(2);
      ROC_tau32_antikt_wta->Write();
      ROC_multigraph->Add(ROC_tau32_antikt_wta);
      leg_ROC->AddEntry(ROC_tau32_antikt_wta, "ak_{T}", "L");

      TH1F* tau32_ttbar_antikt = (TH1F*)tau32_ttbar_antikt_hists.At(B);
      TH1F* tau32_dijets_antikt = (TH1F*)tau32_dijets_antikt_hists.At(B);
      tau32_ttbar_antikt->SetStats(0);
      tau32_dijets_antikt->SetStats(0);
      tau32_ttbar_antikt->Write();
      tau32_dijets_antikt->Write();

      TH1F *tau32_ttbar_antikt_total = (TH1F*)tau32_ttbar_antikt_total_hists.At(B);
      tau32_ttbar_antikt_total->Add(tau32_ttbar_antikt);
      TH1F *tau32_dijets_antikt_total = (TH1F*)tau32_dijets_antikt_total_hists.At(B);
      tau32_dijets_antikt_total->Add(tau32_dijets_antikt);

      int n_antikt = tau32_ttbar_antikt->GetSize() - 2;
      double integral_dijets_antikt[n_antikt], integral_ttbar_antikt[n_antikt];
      for (int i = 0; i < n_antikt; i++) {
        integral_ttbar_antikt[i] = tau32_ttbar_antikt->Integral(0,i)/total_ttbar_jets;
        integral_dijets_antikt[i] = tau32_dijets_antikt->Integral(0,i)/total_dijets_jets;
      }
      TGraph* ROC_tau32_std = new TGraph(n_antikt, integral_ttbar_antikt, integral_dijets_antikt);
      ROC_tau32_std->GetXaxis()->SetLimits(0, 1);
      ROC_tau32_std->GetYaxis()->SetLimits(0, 1);
      ROC_tau32_std->SetLineColor(8);
      ROC_tau32_std->SetLineWidth(2);
      ROC_tau32_std->SetMarkerStyle(5);
      ROC_tau32_std->SetMarkerSize(2);
      ROC_tau32_std->Write();
    // ROC_multigraph->Add(ROC_tau32_std);
    // leg_ROC->AddEntry(ROC_tau32_std, "Genk_{T}", "L");

      TH1F* tau32_ttbar_manual_6jet = (TH1F*)tau32_ttbar_manual_6jet_hists.At(B);
      TH1F* tau32_dijets_manual_6jet = (TH1F*)tau32_dijets_manual_6jet_hists.At(B);
      TH1F* rawmass_ttbar_manual = (TH1F*)rawmass_ttbar_manual_hists.At(B);
      TH1F* rawmass_dijets_manual = (TH1F*)rawmass_dijets_manual_hists.At(B);

      TH1F *tau32_ttbar_manual_6jet_total = (TH1F*)tau32_ttbar_manual_6jet_total_hists.At(B);
      tau32_ttbar_manual_6jet_total->Add(tau32_ttbar_manual_6jet);
      TH1F *tau32_dijets_manual_6jet_total = (TH1F*)tau32_dijets_manual_6jet_total_hists.At(B);
      tau32_dijets_manual_6jet_total->Add(tau32_dijets_manual_6jet);

      tau32_ttbar_manual_6jet->SetStats(0);
      tau32_dijets_manual_6jet->SetStats(0);
      tau32_ttbar_manual_6jet->Write();
      tau32_dijets_manual_6jet->Write();

      int n_size_6jet = tau32_ttbar_manual_6jet->GetSize() - 2;
      double integral_dijets_6jet[n_size_6jet], integral_ttbar_6jet[n_size_6jet];
      for (int i = 0; i < n_size_6jet; i++) {
      // integral_ttbar_6jet[i] = tau32_ttbar_manual_6jet->Integral(0,i)/total_ttbar_jets;
      // integral_dijets_6jet[i] = tau32_dijets_manual_6jet->Integral(0,i)/total_dijets_jets;
        integral_ttbar_6jet[i] = tau32_ttbar_manual_6jet->Integral(0,i)/rawmass_ttbar_manual->Integral(0, rawmass_ttbar_manual->GetNbinsX() + 1);
        integral_dijets_6jet[i] = tau32_dijets_manual_6jet->Integral(0,i)/rawmass_dijets_manual->Integral(0, rawmass_dijets_manual->GetNbinsX() + 1);
      }
      TGraph* ROC_tau32_manual_6jet = new TGraph(n_size_6jet, integral_ttbar_6jet, integral_dijets_6jet);
      ROC_tau32_manual_6jet->GetXaxis()->SetLimits(0, 1);
      ROC_tau32_manual_6jet->GetYaxis()->SetLimits(0, 1);
      ROC_tau32_manual_6jet->SetLineColor(kOrange);
      ROC_tau32_manual_6jet->SetLineWidth(2);
      ROC_tau32_manual_6jet->SetMarkerStyle(5);
      ROC_tau32_manual_6jet->SetMarkerSize(2);
    // ROC_multigraph->Add(ROC_tau32_manual_6jet);
    // leg_ROC->AddEntry(ROC_tau32_manual_6jet, "manual_6jet 6 Axes (beta = " + (TString)ss.str() + ")", "L");
    // leg_ROC->AddEntry(ROC_tau32_manual_6jet, "6-jet", "L");
    // ROC_tau32_manual_6jet->Write();


      TH1F* tau32_ttbar_manual = (TH1F*)tau32_ttbar_manual_hists.At(B);
      TH1F* tau32_dijets_manual = (TH1F*)tau32_dijets_manual_hists.At(B);
      TH1F* mass_ttbar_manual = (TH1F*)mass_ttbar_manual_hists.At(B);
      TH1F* mass_dijets_manual = (TH1F*)mass_dijets_manual_hists.At(B);

      TH1F *tau32_ttbar_manual_total = (TH1F*)tau32_ttbar_manual_total_hists.At(B);
      tau32_ttbar_manual_total->Add(tau32_ttbar_manual);
      TH1F *tau32_dijets_manual_total = (TH1F*)tau32_dijets_manual_total_hists.At(B);
      tau32_dijets_manual_total->Add(tau32_dijets_manual);

      tau32_ttbar_manual->SetStats(0);
      tau32_dijets_manual->SetStats(0);
      tau32_ttbar_manual->Write();
      tau32_dijets_manual->Write();

      int n_size = tau32_ttbar_manual->GetSize() - 2;
      double integral_dijets[n_size], integral_ttbar[n_size];
      for (int i = 0; i < n_size; i++) {
      // integral_ttbar[i] = tau32_ttbar_manual->Integral(0,i)/total_ttbar_jets;
      // integral_dijets[i] = tau32_dijets_manual->Integral(0,i)/total_dijets_jets;
        integral_ttbar[i] = tau32_ttbar_manual->Integral(0,i)/mass_ttbar_manual->Integral(0, mass_ttbar_manual->GetNbinsX() + 1);
        integral_dijets[i] = tau32_dijets_manual->Integral(0,i)/mass_dijets_manual->Integral(0, mass_dijets_manual->GetNbinsX() + 1);
      }
      TGraph* ROC_tau32_manual = new TGraph(n_size, integral_ttbar, integral_dijets);
      ROC_tau32_manual->GetXaxis()->SetLimits(0, 1);
      ROC_tau32_manual->GetYaxis()->SetLimits(0, 1);
      ROC_tau32_manual->SetLineColor(kBlue);
      ROC_tau32_manual->SetLineWidth(2);
      ROC_tau32_manual->SetMarkerStyle(5);
      ROC_tau32_manual->SetMarkerSize(2);
      ROC_multigraph->Add(ROC_tau32_manual);
    // leg_ROC->AddEntry(ROC_tau32_manual, "Manual 6 Axes (beta = " + (TString)ss.str() + ")", "L");
      leg_ROC->AddEntry(ROC_tau32_manual, "6-jet imp", "L");
    // ROC_tau32_manual->Write();


      TH1F* tau32_ttbar_manual_2groups = (TH1F*)tau32_ttbar_manual_2groups_hists.At(B);
      TH1F* tau32_dijets_manual_2groups = (TH1F*)tau32_dijets_manual_2groups_hists.At(B);
      TH1F* mass_ttbar_manual_2groups = (TH1F*)mass_ttbar_manual_2groups_hists.At(B);
      TH1F* mass_dijets_manual_2groups = (TH1F*)mass_dijets_manual_2groups_hists.At(B);

      TH1F *tau32_ttbar_manual_2groups_total = (TH1F*)tau32_ttbar_manual_2groups_total_hists.At(B);
      tau32_ttbar_manual_2groups_total->Add(tau32_ttbar_manual_2groups);
      TH1F *tau32_dijets_manual_2groups_total = (TH1F*)tau32_dijets_manual_2groups_total_hists.At(B);
      tau32_dijets_manual_2groups_total->Add(tau32_dijets_manual_2groups);

      tau32_ttbar_manual_2groups->SetStats(0);
      tau32_dijets_manual_2groups->SetStats(0);
      tau32_ttbar_manual_2groups->Write();
      tau32_dijets_manual_2groups->Write();

      int n_size_manual_2groups = tau32_ttbar_manual_2groups->GetSize() - 2;
      double integral_dijets_manual_2groups[n_size_manual_2groups], integral_ttbar_manual_2groups[n_size_manual_2groups];
      for (int i = 0; i < n_size_manual_2groups; i++) {
      // integral_ttbar_manual_2groups[i] = tau32_ttbar_manual_2groups->Integral(0,i)/total_ttbar_jets;
      // integral_dijets_manual_2groups[i] = tau32_dijets_manual_2groups->Integral(0,i)/total_dijets_jets;
        integral_ttbar_manual_2groups[i] = tau32_ttbar_manual_2groups->Integral(0,i)/mass_ttbar_manual_2groups->Integral(0, mass_ttbar_manual_2groups->GetNbinsX() + 1);
        integral_dijets_manual_2groups[i] = tau32_dijets_manual_2groups->Integral(0,i)/mass_dijets_manual_2groups->Integral(0, mass_dijets_manual_2groups->GetNbinsX() + 1);
      }
      TGraph* ROC_tau32_manual_2groups = new TGraph(n_size_manual_2groups, integral_ttbar_manual_2groups, integral_dijets_manual_2groups);
      ROC_tau32_manual_2groups->GetXaxis()->SetLimits(0, 1);
      ROC_tau32_manual_2groups->GetYaxis()->SetLimits(0, 1);
      ROC_tau32_manual_2groups->SetLineColor(kRed);
      ROC_tau32_manual_2groups->SetLineWidth(2);
      ROC_tau32_manual_2groups->SetMarkerStyle(5);
      ROC_tau32_manual_2groups->SetMarkerSize(2);
      ROC_multigraph->Add(ROC_tau32_manual_2groups);
    // leg_ROC->AddEntry(ROC_tau32_manual_2groups, "Manual 2X3 Axes (beta = " + (TString)ss.str() + ")", "L");
      leg_ROC->AddEntry(ROC_tau32_manual_2groups, "2X3-jet", "L");
      ROC_tau32_manual_2groups->Write();


      TH1F* tau32_ttbar_manual_2groups_improved = (TH1F*)tau32_ttbar_manual_2groups_improved_hists.At(B);
      TH1F* tau32_dijets_manual_2groups_improved = (TH1F*)tau32_dijets_manual_2groups_improved_hists.At(B);
      TH1F* mass_ttbar_manual_2groups_improved = (TH1F*)mass_ttbar_manual_2groups_improved_hists.At(B);
      TH1F* mass_dijets_manual_2groups_improved = (TH1F*)mass_dijets_manual_2groups_improved_hists.At(B);

      TH1F *tau32_ttbar_manual_2groups_improved_total = (TH1F*)tau32_ttbar_manual_2groups_improved_total_hists.At(B);
      tau32_ttbar_manual_2groups_improved_total->Add(tau32_ttbar_manual_2groups_improved);
      TH1F *tau32_dijets_manual_2groups_improved_total = (TH1F*)tau32_dijets_manual_2groups_improved_total_hists.At(B);
      tau32_dijets_manual_2groups_improved_total->Add(tau32_dijets_manual_2groups_improved);

      tau32_ttbar_manual_2groups_improved->SetStats(0);
      tau32_dijets_manual_2groups_improved->SetStats(0);
      tau32_ttbar_manual_2groups_improved->Write();
      tau32_dijets_manual_2groups_improved->Write();

      int n_size_manual_2groups_improved = tau32_ttbar_manual_2groups_improved->GetSize() - 2;
      double integral_dijets_manual_2groups_improved[n_size_manual_2groups_improved], integral_ttbar_manual_2groups_improved[n_size_manual_2groups_improved];
      for (int i = 0; i < n_size_manual_2groups_improved; i++) {
      // integral_ttbar_manual_2groups_improved[i] = tau32_ttbar_manual_2groups_improved->Integral(0,i)/total_ttbar_jets;
      // integral_dijets_manual_2groups_improved[i] = tau32_dijets_manual_2groups_improved->Integral(0,i)/total_dijets_jets;
        integral_ttbar_manual_2groups_improved[i] = tau32_ttbar_manual_2groups_improved->Integral(0,i)/mass_ttbar_manual_2groups_improved->Integral(0, mass_ttbar_manual_2groups_improved->GetNbinsX() + 1);
        integral_dijets_manual_2groups_improved[i] = tau32_dijets_manual_2groups_improved->Integral(0,i)/mass_dijets_manual_2groups_improved->Integral(0, mass_dijets_manual_2groups_improved->GetNbinsX() + 1);
      }
      TGraph* ROC_tau32_manual_2groups_improved = new TGraph(n_size_manual_2groups_improved, integral_ttbar_manual_2groups_improved, integral_dijets_manual_2groups_improved);
      ROC_tau32_manual_2groups_improved->GetXaxis()->SetLimits(0, 1);
      ROC_tau32_manual_2groups_improved->GetYaxis()->SetLimits(0, 1);
      ROC_tau32_manual_2groups_improved->SetLineColor(8);
      ROC_tau32_manual_2groups_improved->SetLineWidth(2);
      ROC_tau32_manual_2groups_improved->SetMarkerStyle(5);
      ROC_tau32_manual_2groups_improved->SetMarkerSize(2);
      ROC_multigraph->Add(ROC_tau32_manual_2groups_improved);
    // leg_ROC->AddEntry(ROC_tau32_manual_2groups_improved, "Manual 2X3 Axes (beta = " + (TString)ss.str() + ")", "L");
      leg_ROC->AddEntry(ROC_tau32_manual_2groups_improved, "2X3-jet imp", "L");
      ROC_tau32_manual_2groups_improved->Write();

      TH1F* tau32_ttbar_manual_2groups_bigcone = (TH1F*)tau32_ttbar_manual_2groups_bigcone_hists.At(B);
      TH1F* tau32_dijets_manual_2groups_bigcone = (TH1F*)tau32_dijets_manual_2groups_bigcone_hists.At(B);

      TH1F *tau32_ttbar_manual_2groups_bigcone_total = (TH1F*)tau32_ttbar_manual_2groups_bigcone_total_hists.At(B);
      tau32_ttbar_manual_2groups_bigcone_total->Add(tau32_ttbar_manual_2groups_bigcone);
      TH1F *tau32_dijets_manual_2groups_bigcone_total = (TH1F*)tau32_dijets_manual_2groups_bigcone_total_hists.At(B);
      tau32_dijets_manual_2groups_bigcone_total->Add(tau32_dijets_manual_2groups_bigcone);

      tau32_ttbar_manual_2groups_bigcone->SetStats(0);
      tau32_dijets_manual_2groups_bigcone->SetStats(0);
      tau32_ttbar_manual_2groups_bigcone->Write();
      tau32_dijets_manual_2groups_bigcone->Write();

      int n_size_manual_2groups_bigcone = tau32_ttbar_manual_2groups_bigcone->GetSize() - 2;
      double integral_dijets_manual_2groups_bigcone[n_size_manual_2groups_bigcone], integral_ttbar_manual_2groups_bigcone[n_size_manual_2groups_bigcone];
      for (int i = 0; i < n_size_manual_2groups_bigcone; i++) {
        integral_ttbar_manual_2groups_bigcone[i] = tau32_ttbar_manual_2groups_bigcone->Integral(0,i)/total_ttbar_jets;
        integral_dijets_manual_2groups_bigcone[i] = tau32_dijets_manual_2groups_bigcone->Integral(0,i)/total_dijets_jets;
      }
      TGraph* ROC_tau32_manual_2groups_bigcone = new TGraph(n_size_manual_2groups_bigcone, integral_ttbar_manual_2groups_bigcone, integral_dijets_manual_2groups_bigcone);
      ROC_tau32_manual_2groups_bigcone->GetXaxis()->SetLimits(0, 1);
      ROC_tau32_manual_2groups_bigcone->GetYaxis()->SetLimits(0, 1);
      ROC_tau32_manual_2groups_bigcone->SetLineColor(kYellow);
      ROC_tau32_manual_2groups_bigcone->SetLineWidth(2);
      ROC_tau32_manual_2groups_bigcone->SetMarkerStyle(5);
      ROC_tau32_manual_2groups_bigcone->SetMarkerSize(2);
      ROC_multigraph->Add(ROC_tau32_manual_2groups_bigcone);
    // leg_ROC->AddEntry(ROC_tau32_manual_2groups_bigcone, "Manual 2X3 Axes (beta = " + (TString)ss.str() + ")", "L");
      leg_ROC->AddEntry(ROC_tau32_manual_2groups_bigcone, "2X3-jet bigcone", "L");
      ROC_tau32_manual_2groups_bigcone->Write();

      TH1F* volatility_ttbar_manual = (TH1F*)volatility_ttbar_manual_hists.At(B);
      TH1F* volatility_dijets_manual = (TH1F*)volatility_dijets_manual_hists.At(B);

    // TH1F *volatility_ttbar_manual_total = (TH1F*)volatility_ttbar_manual_total_hists.At(B);
    // volatility_ttbar_manual_total->Add(volatility_ttbar_manual);
    // TH1F *volatility_dijets_manual_total = (TH1F*)volatility_dijets_manual_total_hists.At(B);
    // volatility_dijets_manual_total->Add(volatility_dijets_manual);

      volatility_ttbar_manual->SetStats(0);
      volatility_dijets_manual->SetStats(0);
      volatility_ttbar_manual->Write();
      volatility_dijets_manual->Write();

      int n_size_vol = volatility_ttbar_manual->GetSize() - 2;
      double integral_dijets_vol[n_size], integral_ttbar_vol[n_size];
      for (int i = 0; i < n_size_vol; i++) {
        integral_ttbar_vol[i] = volatility_ttbar_manual->Integral(0,i)/nEvent;
        integral_dijets_vol[i] = volatility_dijets_manual->Integral(0,i)/nEvent;
      }
      TGraph* ROC_volatility_manual = new TGraph(n_size_vol, integral_ttbar_vol, integral_dijets_vol);
      ROC_volatility_manual->GetXaxis()->SetLimits(0, 1);
      ROC_volatility_manual->GetYaxis()->SetLimits(0, 1);
      ROC_volatility_manual->SetLineColor(kOrange);
      ROC_volatility_manual->SetLineWidth(2);
      ROC_volatility_manual->SetMarkerStyle(5);
      ROC_volatility_manual->SetMarkerSize(2);
    // ROC_multigraph->Add(ROC_volatility_manual);
    // leg_ROC->AddEntry(ROC_volatility_manual, "Manual 6 Axes (beta = " + (TString)ss.str() + ")", "L");
    // leg_ROC->AddEntry(ROC_volatility_manual, "volatility", "L");
    // ROC_volatility_manual->Write();

      ROC_multigraph->Draw("AL");
      ROC_multigraph->GetXaxis()->SetTitle("Tagging Efficiency");
      ROC_multigraph->GetYaxis()->SetTitle("Mistag Rate");
      ROC_multigraph->GetYaxis()->SetTitleOffset(0.85);
      leg_ROC->Draw();
      ROC_compare->Write();
      ROC_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ROC_compare_pt" + (TString)ss_pt.str() + "_beta" + (TString)ss.str() + ".eps", "eps");

    }

    TCanvas *mass_ttbar_compare = new TCanvas("mass_ttbar_compare", "Mass Comparison", 600, 600);
    mass_ttbar_compare->cd();
    double mass_ttbar_antikt_scale = 1/mass_ttbar_antikt->Integral(0, mass_ttbar_antikt->GetNbinsX() + 1);
    mass_ttbar_antikt->Scale(mass_ttbar_antikt_scale);
    mass_ttbar_antikt->SetLineColor(kBlack);
    mass_ttbar_antikt->SetTitle("Mass of Top Jets");
    mass_ttbar_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_ttbar_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_ttbar_antikt->Draw();

    double mass_ttbar_antikt_6jets_scale = 1/mass_ttbar_antikt_6jets->Integral(0, mass_ttbar_antikt_6jets->GetNbinsX() + 1);
    mass_ttbar_antikt_6jets->Scale(mass_ttbar_antikt_6jets_scale);
    mass_ttbar_antikt_6jets->SetLineColor(8);
    mass_ttbar_antikt_6jets->Draw("SAMES");

    TLegend *leg_mass_ttbar = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_mass_ttbar->SetLineColor(0);
    leg_mass_ttbar->SetFillColor(0);
    leg_mass_ttbar->AddEntry(mass_ttbar_antikt, "ak_{T} (2-jet)", "L");
    leg_mass_ttbar->AddEntry(mass_ttbar_antikt_6jets, "ak_{T} (6-jet)", "L");

    double max_val_ttbar_akt = (mass_ttbar_antikt->GetMaximum() > mass_ttbar_antikt_6jets->GetMaximum()) ? mass_ttbar_antikt->GetMaximum() : mass_ttbar_antikt_6jets->GetMaximum();
    double max_val_ttbar = max_val_ttbar_akt;

    top_efficiency_mass[0][i_pt] = mass_ttbar_antikt->Integral((double)160/500*mass_ttbar_antikt->GetNbinsX(), (double)240/500*mass_ttbar_antikt->GetNbinsX());
    top_efficiency_mass[1][i_pt] = mass_ttbar_antikt_6jets->Integral((double)160/500*mass_ttbar_antikt_6jets->GetNbinsX(), (double)240/500*mass_ttbar_antikt_6jets->GetNbinsX());

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_ttbar_manual = (TH1*)mass_ttbar_manual_hists.At(B);
      TH1* mass_dijets_manual = (TH1*)mass_dijets_manual_hists.At(B);
      double mass_ttbar_manual_scale = 1/mass_ttbar_manual->Integral(0, mass_ttbar_manual->GetNbinsX() + 1);
      double mass_dijets_manual_scale = 1/mass_dijets_manual->Integral(0, mass_dijets_manual->GetNbinsX() + 1);
      mass_ttbar_manual->Scale(mass_ttbar_manual_scale);
      mass_dijets_manual->Scale(mass_dijets_manual_scale);

      mass_ttbar_manual->SetStats(0);
      mass_dijets_manual->SetStats(0);

      top_efficiency_mass[B+2][i_pt] = mass_ttbar_manual->Integral((double)160/500*mass_ttbar_manual->GetNbinsX(), (double)240/500*mass_ttbar_manual->GetNbinsX());

      mass_ttbar_manual->SetLineColor(colorlist[B]);
    // if (B == 0) mass_ttbar_manual->SetLineColor(kRed);
    // if (B == 1) mass_ttbar_manual->SetLineColor(kBlue);
    // if (B == 2) mass_ttbar_manual->SetLineColor(kRed);
    // if (B == 3) mass_ttbar_manual->SetLineColor(kBlue);

      mass_ttbar_manual->Draw("SAMES");
      leg_mass_ttbar->AddEntry(mass_ttbar_manual, "#beta = " + (TString)ss.str());
      mass_ttbar_manual->Write();

      if (mass_ttbar_manual->GetMaximum() > max_val_ttbar) {
        max_val_ttbar = mass_ttbar_manual->GetMaximum();
      // mass_ttbar_manual->SetMaximum(1.2*max_val_ttbar);
      // mass_ttbar_compare->Update();
      }

    // cout << "manual ttbar mass ratio: " << mass_ttbar_manual->Integral(75, 110) << "(beta = " << beta << ")" << endl;
    }

    mass_ttbar_antikt->SetMaximum(1.2*max_val_ttbar);
    leg_mass_ttbar->Draw();
    mass_ttbar_compare->Write();
    mass_ttbar_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_mass_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *mass_ttbar_compare_2jets = new TCanvas("mass_ttbar_compare_2jets", "Mass Comparison", 600, 600);
    mass_ttbar_compare_2jets->cd();
    mass_ttbar_antikt->SetLineColor(kBlack);
    mass_ttbar_antikt->SetTitle("Mass of Top Jets");
    mass_ttbar_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_ttbar_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_ttbar_antikt->Draw();

    TLegend *leg_mass_ttbar_2jets = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_mass_ttbar_2jets->SetLineColor(0);
    leg_mass_ttbar_2jets->SetFillColor(0);
    leg_mass_ttbar_2jets->AddEntry(mass_ttbar_antikt, "ak_{T} (2-jet)", "L");

    double max_val_ttbar_2jets = max_val_ttbar_akt;

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* rawmass_ttbar_manual_2jets = (TH1*)rawmass_ttbar_manual_2jets_hists.At(B);
      double rawmass_ttbar_manual_2jets_scale = 1/rawmass_ttbar_manual_2jets->Integral(0, rawmass_ttbar_manual_2jets->GetNbinsX() + 1);
      rawmass_ttbar_manual_2jets->Scale(rawmass_ttbar_manual_2jets_scale);

      rawmass_ttbar_manual_2jets->SetStats(0);

      rawmass_ttbar_manual_2jets->SetLineColor(colorlist[B]);
    // if (B == 0) rawmass_ttbar_manual_2jets->SetLineColor(kRed);
    // if (B == 1) rawmass_ttbar_manual_2jets->SetLineColor(kBlue);
    // if (B == 2) rawmass_ttbar_manual_2jets->SetLineColor(kRed);
    // if (B == 3) rawmass_ttbar_manual_2jets->SetLineColor(kBlue);

      rawmass_ttbar_manual_2jets->Draw("SAMES");
      leg_mass_ttbar_2jets->AddEntry(rawmass_ttbar_manual_2jets, "#beta = " + (TString)ss.str());
      rawmass_ttbar_manual_2jets->Write();

      if (rawmass_ttbar_manual_2jets->GetMaximum() > max_val_ttbar_2jets) {
        max_val_ttbar_2jets = rawmass_ttbar_manual_2jets->GetMaximum();
      // rawmass_ttbar_manual_2jets->SetMaximum(1.2*max_val_ttbar);
      // mass_ttbar_compare_2jets->Update();
      }

    // cout << "manual ttbar mass ratio: " << mass_ttbar_manual->Integral(75, 110) << "(beta = " << beta << ")" << endl;
    }

    mass_ttbar_antikt->SetMaximum(1.2*max_val_ttbar_2jets);
    leg_mass_ttbar_2jets->Draw();
    mass_ttbar_compare_2jets->Write();
    mass_ttbar_compare_2jets->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_mass_2jets_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *mass_dijets_compare = new TCanvas("mass_dijets_compare", "Mass Comparison", 600, 600);
    mass_dijets_compare->cd();
    double mass_dijets_antikt_scale = 1/mass_dijets_antikt->Integral(0, mass_dijets_antikt->GetNbinsX() + 1);
    mass_dijets_antikt->Scale(mass_dijets_antikt_scale);
    mass_dijets_antikt->SetLineColor(kBlack);
    mass_dijets_antikt->SetTitle("Mass of QCD jets");
    mass_dijets_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_dijets_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_dijets_antikt->Draw();

    double mass_dijets_antikt_6jets_scale = 1/mass_dijets_antikt_6jets->Integral(0, mass_dijets_antikt_6jets->GetNbinsX() + 1);
    mass_dijets_antikt_6jets->Scale(mass_dijets_antikt_6jets_scale);
    mass_dijets_antikt_6jets->SetLineColor(8);
    mass_dijets_antikt_6jets->Draw("SAMES");

    TLegend *leg_mass_dijets = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_mass_dijets->SetLineColor(0);
    leg_mass_dijets->SetFillColor(0);
    leg_mass_dijets->AddEntry(mass_dijets_antikt, "ak_{T} (2-jet)", "L");
    leg_mass_dijets->AddEntry(mass_dijets_antikt_6jets, "ak_{T} (6-jet)", "L");
  // double max_val_dijets = (mass_dijets_antikt->GetMaximum() > mass_dijets_antikt_6jets->GetMaximum()) ? mass_dijets_antikt->GetMaximum() : mass_dijets_antikt_6jets->GetMaximum();

    double max_val_dijets_akt = (mass_dijets_antikt->GetMaximum() > mass_dijets_antikt_6jets->GetMaximum()) ? mass_dijets_antikt->GetMaximum() : mass_dijets_antikt_6jets->GetMaximum();
    double max_val_dijets = max_val_dijets_akt;

    qcd_efficiency_mass[0][i_pt] = mass_dijets_antikt->Integral((double)160/500*mass_dijets_antikt->GetNbinsX(), (double)240/500*mass_dijets_antikt->GetNbinsX());
    qcd_efficiency_mass[1][i_pt] = mass_dijets_antikt_6jets->Integral((double)160/500*mass_dijets_antikt_6jets->GetNbinsX(), (double)240/500*mass_dijets_antikt_6jets->GetNbinsX());


    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_dijets_manual = (TH1*)mass_dijets_manual_hists.At(B);

      double mass_dijets_manual_scale = 1/mass_dijets_manual->Integral(0, mass_dijets_manual->GetNbinsX() + 1);
      mass_dijets_manual->Scale(mass_dijets_manual_scale);

      mass_dijets_manual->SetStats(0);

      qcd_efficiency_mass[B+2][i_pt] = mass_dijets_manual->Integral((double)160/500*mass_dijets_manual->GetNbinsX(), (double)240/500*mass_dijets_manual->GetNbinsX());

      mass_dijets_manual->SetLineColor(colorlist[B]);
    // if (B == 0) mass_dijets_manual->SetLineColor(kRed);
    // if (B == 1) mass_dijets_manual->SetLineColor(kBlue);
    // if (B == 2) mass_dijets_manual->SetLineColor(kRed);
    // if (B == 3) mass_dijets_manual->SetLineColor(kBlue);

      mass_dijets_manual->Draw("SAMES");
      leg_mass_dijets->AddEntry(mass_dijets_manual, "#beta = " + (TString)ss.str());
      mass_dijets_manual->Write();

      if (mass_dijets_manual->GetMaximum() > max_val_dijets) {
        max_val_dijets = mass_dijets_manual->GetMaximum();
      // mass_dijets_manual->SetMaximum(1.2*max_val_dijets);
      // mass_dijets_compare->Update();
      }

    // cout << "manual QCD mass ratio: " << mass_dijets_manual->Integral(75, 110) << "(beta = " << beta << ")" << endl;
    }

    mass_dijets_antikt->SetMaximum(1.2*max_val_dijets);
    leg_mass_dijets->Draw();
    mass_dijets_compare->Write();
    mass_dijets_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_mass_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *mass_dijets_compare_2jets = new TCanvas("mass_dijets_compare_2jets", "Mass Comparison", 600, 600);
    mass_dijets_compare_2jets->cd();
    mass_dijets_antikt->SetTitle("Mass of Top Jets");
    mass_dijets_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_dijets_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_dijets_antikt->Draw();

    TLegend *leg_mass_dijets_2jets = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_mass_dijets_2jets->SetLineColor(0);
    leg_mass_dijets_2jets->SetFillColor(0);
    leg_mass_dijets_2jets->AddEntry(mass_dijets_antikt, "ak_{T} (2-jet)", "L");

    double max_val_dijets_2jets = max_val_dijets_akt;

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* rawmass_dijets_manual_2jets = (TH1*)rawmass_dijets_manual_2jets_hists.At(B);
      double rawmass_dijets_manual_2jets_scale = 1/rawmass_dijets_manual_2jets->Integral(0, rawmass_dijets_manual_2jets->GetNbinsX() + 1);
      rawmass_dijets_manual_2jets->Scale(rawmass_dijets_manual_2jets_scale);

      rawmass_dijets_manual_2jets->SetStats(0);

      rawmass_dijets_manual_2jets->SetLineColor(colorlist[B]);
    // if (B == 0) rawmass_dijets_manual_2jets->SetLineColor(kRed);
    // if (B == 1) rawmass_dijets_manual_2jets->SetLineColor(kBlue);
    // if (B == 2) rawmass_dijets_manual_2jets->SetLineColor(kRed);
    // if (B == 3) rawmass_dijets_manual_2jets->SetLineColor(kBlue);

      rawmass_dijets_manual_2jets->Draw("SAMES");
      leg_mass_dijets_2jets->AddEntry(rawmass_dijets_manual_2jets, "#beta = " + (TString)ss.str());
      rawmass_dijets_manual_2jets->Write();

      if (rawmass_dijets_manual_2jets->GetMaximum() > max_val_dijets_2jets) {
        max_val_dijets_2jets = rawmass_dijets_manual_2jets->GetMaximum();
      // rawmass_dijets_manual_2jets->SetMaximum(1.2*max_val_dijets);
      // mass_dijets_compare_2jets->Update();
      }

    // cout << "manual dijets mass ratio: " << mass_dijets_manual->Integral(75, 110) << "(beta = " << beta << ")" << endl;
    }

    mass_dijets_antikt->SetMaximum(1.2*max_val_dijets_2jets);
    leg_mass_dijets_2jets->Draw();
    mass_dijets_compare_2jets->Write();
    mass_dijets_compare_2jets->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_mass_2jets_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");

    TCanvas *area_ttbar_compare = new TCanvas("area_ttbar_compare", "area Comparison", 600, 600);
    area_ttbar_compare->cd();
    double area_ttbar_antikt_scale = 1/area_ttbar_antikt->Integral(0, area_ttbar_antikt->GetNbinsX() + 1);
    area_ttbar_antikt->Scale(area_ttbar_antikt_scale);
    area_ttbar_antikt->SetStats(0);
    area_ttbar_antikt->SetLineColor(8);
    area_ttbar_antikt->SetTitle("Area of 6 Top Jets");
    area_ttbar_antikt->GetXaxis()->SetTitle("Area");
    area_ttbar_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    area_ttbar_antikt->Draw();

    TLegend *leg_area_ttbar = new TLegend(0.14, 0.65, 0.35, 0.88);
    leg_area_ttbar->SetFillColor(0);
    leg_area_ttbar->SetLineColor(0);
    leg_area_ttbar->AddEntry(area_ttbar_antikt, "ak_{T}", "L");
    double max_val_ttbar_area = area_ttbar_antikt->GetMaximum();

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* area_ttbar_manual = (TH1*)area_ttbar_manual_hists.At(B);
      TH1* area_dijets_manual = (TH1*)area_dijets_manual_hists.At(B);

      double area_ttbar_manual_scale = 1/area_ttbar_manual->Integral(0, area_ttbar_manual->GetNbinsX() + 1);
      double area_dijets_manual_scale = 1/area_dijets_manual->Integral(0, area_dijets_manual->GetNbinsX() + 1);
      area_ttbar_manual->Scale(area_ttbar_manual_scale);
      area_dijets_manual->Scale(area_dijets_manual_scale);

      area_ttbar_manual->SetStats(0);
      area_dijets_manual->SetStats(0);
    // area_ttbar_manual->SetTitle("Area of 6 Top Jets");
    // area_ttbar_manual->GetXaxis()->SetTitle("Area");
    // area_ttbar_manual->GetYaxis()->SetTitle("Relative Occurrence");

      area_ttbar_manual->SetLineColor(colorlist[B]);
    // if (B == 0) area_ttbar_manual->SetLineColor(kRed);
    // if (B == 1) area_ttbar_manual->SetLineColor(kBlue);
    // if (B == 2) area_ttbar_manual->SetLineColor(kRed);
    // if (B == 3) area_ttbar_manual->SetLineColor(kBlue);

    // if (B == 0) area_ttbar_manual->Draw();
      area_ttbar_manual->Draw("SAMES");
      leg_area_ttbar->AddEntry(area_ttbar_manual, "#beta = " + (TString)ss.str());
      area_ttbar_manual->Write();

      if (area_ttbar_manual->GetMaximum() > max_val_ttbar_area) {
        max_val_ttbar_area = area_ttbar_manual->GetMaximum();
      // area_ttbar_manual->SetMaximum(1.2*max_val_ttbar_area);
      // area_ttbar_compare->Update();
      }
    }

    TLine *true_area = new TLine();
    true_area->SetX1(TMath::Pi()*Rparam*Rparam);
    true_area->SetVertical(true);
    true_area->SetY1(0);
    true_area->SetY2(1.2*max_val_ttbar_area);
    true_area->SetLineColor(kRed);
    true_area->SetLineWidth(3);
    true_area->SetLineStyle(7);
    true_area->Draw();

    area_ttbar_antikt->SetMaximum(1.2*max_val_ttbar_area);
    leg_area_ttbar->Draw();
    area_ttbar_compare->Write();
    area_ttbar_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_area_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");

    TCanvas *area_dijets_compare = new TCanvas("area_dijets_compare", "area Comparison", 600, 600);
    area_dijets_compare->cd();
    double area_dijets_antikt_scale = 1/area_dijets_antikt->Integral(0, area_dijets_antikt->GetNbinsX() + 1);
    area_dijets_antikt->Scale(area_dijets_antikt_scale);
    area_dijets_antikt->SetStats(0);
    area_dijets_antikt->SetLineColor(8);
    area_dijets_antikt->SetTitle("Area of 6 QCD Jets");
    area_dijets_antikt->GetXaxis()->SetTitle("Area");
    area_dijets_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    area_dijets_antikt->Draw();

    TLegend *leg_area_dijets = new TLegend(0.14, 0.65, 0.35, 0.88);
    leg_area_dijets->SetLineColor(0);
    leg_area_dijets->SetFillColor(0);
    leg_area_dijets->AddEntry(area_dijets_antikt, "ak_{T}", "L");
    double max_val_dijets_area = area_dijets_antikt->GetMaximum();

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* area_dijets_manual = (TH1*)area_dijets_manual_hists.At(B);

      double area_dijets_manual_scale = 1/area_dijets_manual->Integral(0, area_dijets_manual->GetNbinsX() + 1);
      area_dijets_manual->Scale(area_dijets_manual_scale);

      area_dijets_manual->SetStats(0);
    // area_dijets_manual->SetTitle("Area of 6 QCD Jets");
    // area_dijets_manual->GetXaxis()->SetTitle("Area");
    // area_dijets_manual->GetYaxis()->SetTitle("Relative Occurrence");

      area_dijets_manual->SetLineColor(colorlist[B]);
    // if (B == 0) area_dijets_manual->SetLineColor(kRed);
    // if (B == 1) area_dijets_manual->SetLineColor(kBlue);
    // if (B == 2) area_dijets_manual->SetLineColor(kRed);
    // if (B == 3) area_dijets_manual->SetLineColor(kBlue);

    // if (B == 0) area_dijets_manual->Draw();
      area_dijets_manual->Draw("SAMES");
      leg_area_dijets->AddEntry(area_dijets_manual, "#beta = " + (TString)ss.str());
      area_dijets_manual->Write();

      if (area_dijets_manual->GetMaximum() > max_val_dijets_area) {
        max_val_dijets_area = area_dijets_manual->GetMaximum();
      // area_dijets_manual->SetMaximum(1.2*max_val_dijets_area);
      // area_dijets_compare->Update();
      }
    }

    true_area->SetY2(1.2*max_val_dijets_area);
    true_area->Draw();

    area_dijets_antikt->SetMaximum(1.2*max_val_dijets_area);
    leg_area_dijets->Draw();
    area_dijets_compare->Write();
    area_dijets_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_area_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *rawmass_ttbar_compare = new TCanvas("rawmass_ttbar_compare", "rawmass Comparison", 600, 600);
    rawmass_ttbar_compare->cd();
    mass_ttbar_antikt->SetLineColor(kBlack);
    mass_ttbar_antikt->SetTitle("Mass of Top Jets");
    mass_ttbar_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_ttbar_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_ttbar_antikt->Draw();
    mass_ttbar_antikt_6jets->Draw("SAMES");

    TLegend *leg_rawmass_ttbar = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_rawmass_ttbar->SetLineColor(0);
    leg_rawmass_ttbar->SetFillColor(0);
    leg_rawmass_ttbar->AddEntry(mass_ttbar_antikt, "ak_{T} (2-jet)", "L");
    leg_rawmass_ttbar->AddEntry(mass_ttbar_antikt_6jets, "ak_{T} (6-jet)", "L");
  // double max_val_ttbar_raw = (mass_ttbar_antikt->GetMaximum() > mass_ttbar_antikt_6jets->GetMaximum()) ? mass_ttbar_antikt->GetMaximum() : mass_ttbar_antikt_6jets->GetMaximum();
    double max_val_ttbar_raw = max_val_ttbar_akt;

    top_efficiency_rawmass[0][i_pt] = mass_ttbar_antikt->Integral((double)160/500*mass_ttbar_antikt->GetNbinsX(), (double)240/500*mass_ttbar_antikt->GetNbinsX());
    top_efficiency_rawmass[1][i_pt] = mass_ttbar_antikt_6jets->Integral((double)160/500*mass_ttbar_antikt_6jets->GetNbinsX(), (double)240/500*mass_ttbar_antikt_6jets->GetNbinsX());

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* rawmass_ttbar_manual = (TH1*)rawmass_ttbar_manual_hists.At(B);
      TH1* rawmass_dijets_manual = (TH1*)rawmass_dijets_manual_hists.At(B);
      double rawmass_ttbar_manual_scale = 1/rawmass_ttbar_manual->Integral(0, rawmass_ttbar_manual->GetNbinsX() + 1);
      double rawmass_dijets_manual_scale = 1/rawmass_dijets_manual->Integral(0, rawmass_dijets_manual->GetNbinsX() + 1);
      rawmass_ttbar_manual->Scale(rawmass_ttbar_manual_scale);
      rawmass_dijets_manual->Scale(rawmass_dijets_manual_scale);

      TH1* mass_ttbar_manual = (TH1*)mass_ttbar_manual_hists.At(B);
      TH1* mass_dijets_manual = (TH1*)mass_dijets_manual_hists.At(B);
      double mass_ttbar_manual_scale = 1/mass_ttbar_manual->Integral(0, mass_ttbar_manual->GetNbinsX() + 1);
      double mass_dijets_manual_scale = 1/mass_dijets_manual->Integral(0, mass_dijets_manual->GetNbinsX() + 1);
      mass_ttbar_manual->Scale(mass_ttbar_manual_scale);
      mass_dijets_manual->Scale(mass_dijets_manual_scale);

      top_efficiency_ratio[B][i_pt] = (double)(mass_ttbar_manual->Integral((double)160/500*mass_ttbar_manual->GetNbinsX(), (double)240/500*mass_ttbar_manual->GetNbinsX()))/(rawmass_ttbar_manual->Integral((double)160/500*rawmass_ttbar_manual->GetNbinsX(), (double)240/500*rawmass_ttbar_manual->GetNbinsX()));

      qcd_efficiency_ratio[B][i_pt] = (double)(mass_dijets_manual->Integral((double)160/500*mass_dijets_manual->GetNbinsX(), (double)240/500*mass_dijets_manual->GetNbinsX()))/(rawmass_dijets_manual->Integral((double)160/500*rawmass_dijets_manual->GetNbinsX(), (double)240/500*rawmass_dijets_manual->GetNbinsX()));

      rawmass_ttbar_manual->SetStats(0);
      rawmass_dijets_manual->SetStats(0);

      top_efficiency_rawmass[B+2][i_pt] = rawmass_ttbar_manual->Integral((double)160/500*rawmass_ttbar_manual->GetNbinsX(), (double)240/500*rawmass_ttbar_manual->GetNbinsX());

      rawmass_ttbar_manual->SetLineColor(colorlist[B]);
    // if (B == 0) rawmass_ttbar_manual->SetLineColor(kRed);
    // if (B == 1) rawmass_ttbar_manual->SetLineColor(kBlue);
    // if (B == 2) rawmass_ttbar_manual->SetLineColor(kRed);
    // if (B == 3) rawmass_ttbar_manual->SetLineColor(kBlue);

      rawmass_ttbar_manual->Draw("SAMES");
      leg_rawmass_ttbar->AddEntry(rawmass_ttbar_manual, "#beta = " + (TString)ss.str());
      rawmass_ttbar_manual->Write();

      if (rawmass_ttbar_manual->GetMaximum() > max_val_ttbar_raw) {
        max_val_ttbar_raw = rawmass_ttbar_manual->GetMaximum();
      // rawmass_ttbar_manual->SetMaximum(1.2*max_val_ttbar_raw);
      // rawmass_ttbar_compare->Update();
      }
    }

    mass_ttbar_antikt->SetMaximum(1.2*max_val_ttbar_raw);
    leg_rawmass_ttbar->Draw();
    rawmass_ttbar_compare->Write();
    rawmass_ttbar_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_rawmass_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");

    TCanvas *rawmass_dijets_compare = new TCanvas("rawmass_dijets_compare", "rawmass Comparison", 600, 600);
    rawmass_dijets_compare->cd();
    mass_dijets_antikt->SetLineColor(kBlack);
    mass_dijets_antikt->SetTitle("Mass of QCD Jets");
    mass_dijets_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_dijets_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_dijets_antikt->Draw();
    mass_dijets_antikt_6jets->Draw("SAMES");

    TLegend *leg_rawmass_dijets = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_rawmass_dijets->SetLineColor(0);
    leg_rawmass_dijets->SetFillColor(0);
    leg_rawmass_dijets->AddEntry(mass_dijets_antikt, "ak_{T} (2-jet)", "L");
    leg_rawmass_dijets->AddEntry(mass_dijets_antikt_6jets, "ak_{T} (6-jet)", "L");
  // double max_val_dijets_raw = (mass_dijets_antikt->GetMaximum() > mass_dijets_antikt_6jets->GetMaximum()) ? mass_dijets_antikt->GetMaximum() : mass_dijets_antikt_6jets->GetMaximum();
    double max_val_dijets_raw = max_val_dijets_akt;

    qcd_efficiency_rawmass[0][i_pt] = mass_dijets_antikt->Integral((double)160/500*mass_dijets_antikt->GetNbinsX(), (double)240/500*mass_dijets_antikt->GetNbinsX());
    qcd_efficiency_rawmass[1][i_pt] = mass_dijets_antikt_6jets->Integral((double)160/500*mass_dijets_antikt_6jets->GetNbinsX(), (double)240/500*mass_dijets_antikt_6jets->GetNbinsX());

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* rawmass_dijets_manual = (TH1*)rawmass_dijets_manual_hists.At(B);

      double rawmass_dijets_manual_scale = 1/rawmass_dijets_manual->Integral(0, rawmass_dijets_manual->GetNbinsX() + 1);
      rawmass_dijets_manual->Scale(rawmass_dijets_manual_scale);

      rawmass_dijets_manual->SetStats(0);

      qcd_efficiency_rawmass[B+2][i_pt] = rawmass_dijets_manual->Integral((double)160/500*rawmass_dijets_manual->GetNbinsX(), (double)240/500*rawmass_dijets_manual->GetNbinsX());

      rawmass_dijets_manual->SetLineColor(colorlist[B]);
      // if (B == 0) rawmass_dijets_manual->SetLineColor(kRed);
      // if (B == 1) rawmass_dijets_manual->SetLineColor(kBlue);
      // if (B == 2) rawmass_dijets_manual->SetLineColor(kRed);
      // if (B == 3) rawmass_dijets_manual->SetLineColor(kBlue);

      rawmass_dijets_manual->Draw("SAMES");
      leg_rawmass_dijets->AddEntry(rawmass_dijets_manual, "#beta = " + (TString)ss.str());
      rawmass_dijets_manual->Write();

      if (rawmass_dijets_manual->GetMaximum() > max_val_dijets_raw) {
        max_val_dijets_raw = rawmass_dijets_manual->GetMaximum();
        // rawmass_dijets_manual->SetMaximum(1.2*max_val_dijets_raw);
        // rawmass_dijets_compare->Update();
      }
    }

    mass_dijets_antikt->SetMaximum(1.2*max_val_dijets_raw);
    leg_rawmass_dijets->Draw();
    rawmass_dijets_compare->Write();
    rawmass_dijets_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_rawmass_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *mass_2groups_ttbar_compare = new TCanvas("mass_2groups_ttbar_compare", "mass_2groups Comparison", 600, 600);
    mass_2groups_ttbar_compare->cd();
    mass_ttbar_antikt->SetLineColor(kBlack);
    mass_ttbar_antikt->SetTitle("Mass of Top Jets");
    mass_ttbar_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_ttbar_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_ttbar_antikt->Draw();

    double mass_ttbar_antikt_2groups_scale = 1/mass_ttbar_antikt_2groups->Integral(0, mass_ttbar_antikt_2groups->GetNbinsX() + 1);
    // double mass_ttbar_antikt_2groups_scale = 1/total_ttbar_jets;
    mass_ttbar_antikt_2groups->Scale(mass_ttbar_antikt_2groups_scale);
    mass_ttbar_antikt_2groups->SetLineColor(8);
    mass_ttbar_antikt_2groups->Draw("SAMES");

    TLegend *leg_mass_2groups_ttbar = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_mass_2groups_ttbar->SetLineColor(0);
    leg_mass_2groups_ttbar->SetFillColor(0);
    leg_mass_2groups_ttbar->AddEntry(mass_ttbar_antikt, "ak_{T} (2-jet)", "L");
    leg_mass_2groups_ttbar->AddEntry(mass_ttbar_antikt_2groups, "ak_{T} (2 groups)", "L");
    // double max_val_ttbar_2groups = (mass_ttbar_antikt->GetMaximum() > mass_ttbar_antikt_6jets->GetMaximum()) ? mass_ttbar_antikt->GetMaximum() : mass_ttbar_antikt_6jets->GetMaximum();
    double max_val_ttbar_2groups = max_val_ttbar_akt;

    top_efficiency_mass_2groups[0][i_pt] = mass_ttbar_antikt->Integral((double)160/500*mass_ttbar_antikt->GetNbinsX(), (double)240/500*mass_ttbar_antikt->GetNbinsX());
    top_efficiency_mass_2groups[1][i_pt] = mass_ttbar_antikt_2groups->Integral((double)160/500*mass_ttbar_antikt_2groups->GetNbinsX(), (double)240/500*mass_ttbar_antikt_2groups->GetNbinsX());

    double mass_ttbar_manual_2groups_scale;
    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_ttbar_manual_2groups = (TH1*)mass_ttbar_manual_2groups_hists.At(B);
      // TH1* mass_dijets_manual_2groups = (TH1*)mass_dijets_manual_2groups_hists.At(B);

      mass_ttbar_manual_2groups_scale = 1/mass_ttbar_manual_2groups->Integral(0, mass_ttbar_manual_2groups->GetNbinsX() + 1);
      // double mass_dijets_manual_2groups_scale = 1/mass_dijets_manual_2groups->Integral(0, mass_dijets_manual_2groups->GetNbinsX() + 1);

      // double mass_ttbar_manual_2groups_scale = 1/total_ttbar_jets;
      // double mass_dijets_manual_2groups_scale = 1/total_dijets_jets;
      mass_ttbar_manual_2groups->Scale(mass_ttbar_manual_2groups_scale);
      // mass_dijets_manual_2groups->Scale(mass_dijets_manual_2groups_scale);

      mass_ttbar_manual_2groups->SetStats(0);
      // mass_dijets_manual_2groups->SetStats(0);

      top_efficiency_mass_2groups[B+2][i_pt] = mass_ttbar_manual_2groups->Integral((double)160/500*mass_ttbar_manual_2groups->GetNbinsX(), (double)240/500*mass_ttbar_manual_2groups->GetNbinsX());

      mass_ttbar_manual_2groups->SetLineColor(colorlist[B]);
      // if (B == 0) mass_ttbar_manual_2groups->SetLineColor(kRed);
      // if (B == 1) mass_ttbar_manual_2groups->SetLineColor(kBlue);
      // if (B == 2) mass_ttbar_manual_2groups->SetLineColor(kRed);
      // if (B == 3) mass_ttbar_manual_2groups->SetLineColor(kBlue);

      mass_ttbar_manual_2groups->Draw("SAMES");
      leg_mass_2groups_ttbar->AddEntry(mass_ttbar_manual_2groups, "#beta = " + (TString)ss.str());
      mass_ttbar_manual_2groups->Write();

      if (mass_ttbar_manual_2groups->GetMaximum() > max_val_ttbar_2groups) {
        max_val_ttbar_2groups = mass_ttbar_manual_2groups->GetMaximum();
        // mass_ttbar_manual_2groups->SetMaximum(1.2*max_val_ttbar_2groups);
        // mass_2groups_ttbar_compare->Update();
      }
    }

    mass_ttbar_antikt->SetMaximum(1.2*max_val_ttbar_2groups);
    leg_mass_2groups_ttbar->Draw();
    mass_2groups_ttbar_compare->Write();
    mass_2groups_ttbar_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_mass_2groups_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");

    TCanvas *mass_2groups_dijets_compare = new TCanvas("mass_2groups_dijets_compare", "mass_2groups Comparison", 600, 600);
    mass_2groups_dijets_compare->cd();
    mass_dijets_antikt->SetLineColor(kBlack);
    mass_dijets_antikt->SetTitle("Mass of QCD Jets");
    mass_dijets_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_dijets_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_dijets_antikt->Draw();

    double mass_dijets_antikt_2groups_scale = 1/mass_dijets_antikt_2groups->Integral(0, mass_dijets_antikt_2groups->GetNbinsX() + 1);
    // double mass_dijets_antikt_2groups_scale = 1/total_dijets_jets;
    mass_dijets_antikt_2groups->Scale(mass_dijets_antikt_2groups_scale);
    mass_dijets_antikt_2groups->SetLineColor(8);
    mass_dijets_antikt_2groups->Draw("SAMES");

    TLegend *leg_mass_2groups_dijets = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_mass_2groups_dijets->SetLineColor(0);
    leg_mass_2groups_dijets->SetFillColor(0);
    leg_mass_2groups_dijets->AddEntry(mass_dijets_antikt, "ak_{T} (2-jet)", "L");
    leg_mass_2groups_dijets->AddEntry(mass_dijets_antikt_2groups, "ak_{T} (2 groups)", "L");
    // double max_val_dijets_2groups = (mass_dijets_antikt->GetMaximum() > mass_dijets_antikt_6jets->GetMaximum()) ? mass_dijets_antikt->GetMaximum() : mass_dijets_antikt_6jets->GetMaximum();
    double max_val_dijets_2groups = max_val_dijets_akt;

    qcd_efficiency_mass_2groups[0][i_pt] = mass_dijets_antikt->Integral((double)160/500*mass_dijets_antikt->GetNbinsX(), (double)240/500*mass_dijets_antikt->GetNbinsX());
    qcd_efficiency_mass_2groups[1][i_pt] = mass_dijets_antikt_2groups->Integral((double)160/500*mass_dijets_antikt_2groups->GetNbinsX(), (double)240/500*mass_dijets_antikt_2groups->GetNbinsX());

    double mass_dijets_manual_2groups_scale;
    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_dijets_manual_2groups = (TH1*)mass_dijets_manual_2groups_hists.At(B);

      mass_dijets_manual_2groups_scale = 1/mass_dijets_manual_2groups->Integral(0, mass_dijets_manual_2groups->GetNbinsX() + 1);
      // double mass_dijets_manual_2groups_scale = 1/total_dijets_jets;
      mass_dijets_manual_2groups->Scale(mass_dijets_manual_2groups_scale);

      mass_dijets_manual_2groups->SetStats(0);

      qcd_efficiency_mass_2groups[B+2][i_pt] = mass_dijets_manual_2groups->Integral((double)160/500*mass_dijets_manual_2groups->GetNbinsX(), (double)240/500*mass_dijets_manual_2groups->GetNbinsX());

      mass_dijets_manual_2groups->SetLineColor(colorlist[B]);
      // if (B == 0) mass_dijets_manual_2groups->SetLineColor(kRed);
      // if (B == 1) mass_dijets_manual_2groups->SetLineColor(kBlue);
      // if (B == 2) mass_dijets_manual_2groups->SetLineColor(kRed);
      // if (B == 3) mass_dijets_manual_2groups->SetLineColor(kBlue);

      mass_dijets_manual_2groups->Draw("SAMES");
      leg_mass_2groups_dijets->AddEntry(mass_dijets_manual_2groups, "#beta = " + (TString)ss.str());
      mass_dijets_manual_2groups->Write();

      if (mass_dijets_manual_2groups->GetMaximum() > max_val_dijets_2groups) {
        max_val_dijets_2groups = mass_dijets_manual_2groups->GetMaximum();
        // mass_dijets_manual_2groups->SetMaximum(1.2*max_val_dijets_2groups);
        // mass_2groups_dijets_compare->Update();
      }
    }

    mass_dijets_antikt->SetMaximum(1.2*max_val_dijets_2groups);
    leg_mass_2groups_dijets->Draw();
    mass_2groups_dijets_compare->Write();
    mass_2groups_dijets_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_mass_2groups_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *mass_2groups_improved_ttbar_compare = new TCanvas("mass_2groups_improved_ttbar_compare", "mass_2groups_improved Comparison", 600, 600);
    mass_2groups_improved_ttbar_compare->cd();
    mass_ttbar_antikt->SetLineColor(kBlack);
    mass_ttbar_antikt->SetTitle("Mass of Top Jets");
    mass_ttbar_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_ttbar_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_ttbar_antikt->Draw();
    mass_ttbar_antikt_6jets->Draw("SAMES");

    TLegend *leg_mass_2groups_improved_ttbar = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_mass_2groups_improved_ttbar->SetLineColor(0);
    leg_mass_2groups_improved_ttbar->SetFillColor(0);
    leg_mass_2groups_improved_ttbar->AddEntry(mass_ttbar_antikt, "ak_{T} (2-jet)", "L");
    leg_mass_2groups_improved_ttbar->AddEntry(mass_ttbar_antikt_6jets, "ak_{T} (6-jet)", "L");
    // double max_val_ttbar_2groups_improved = (mass_ttbar_antikt->GetMaximum() > mass_ttbar_antikt_6jets->GetMaximum()) ? mass_ttbar_antikt->GetMaximum() : mass_ttbar_antikt_6jets->GetMaximum();
    double max_val_ttbar_2groups_improved = max_val_ttbar_akt;

    top_efficiency_mass_2groups_improved[0][i_pt] = mass_ttbar_antikt->Integral((double)160/500*mass_ttbar_antikt->GetNbinsX(), (double)240/500*mass_ttbar_antikt->GetNbinsX());
    top_efficiency_mass_2groups_improved[1][i_pt] = mass_ttbar_antikt_6jets->Integral((double)160/500*mass_ttbar_antikt_6jets->GetNbinsX(), (double)240/500*mass_ttbar_antikt_6jets->GetNbinsX());

    double mass_ttbar_manual_2groups_improved_scale;
    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_ttbar_manual_2groups_improved = (TH1*)mass_ttbar_manual_2groups_improved_hists.At(B);
      // TH1* mass_dijets_manual_2groups_improved = (TH1*)mass_dijets_manual_2groups_improved_hists.At(B);

      mass_ttbar_manual_2groups_improved_scale = 1/mass_ttbar_manual_2groups_improved->Integral(0, mass_ttbar_manual_2groups_improved->GetNbinsX() + 1);
      // double mass_dijets_manual_2groups_improved_scale = 1/mass_dijets_manual_2groups_improved->Integral(0, mass_dijets_manual_2groups_improved->GetNbinsX() + 1);
      mass_ttbar_manual_2groups_improved->Scale(mass_ttbar_manual_2groups_improved_scale);
      // mass_dijets_manual_2groups_improved->Scale(mass_dijets_manual_2groups_improved_scale);

      mass_ttbar_manual_2groups_improved->SetStats(0);
      // mass_dijets_manual_2groups_improved->SetStats(0);

      top_efficiency_mass_2groups_improved[B+2][i_pt] = mass_ttbar_manual_2groups_improved->Integral((double)160/500*mass_ttbar_manual_2groups_improved->GetNbinsX(), (double)240/500*mass_ttbar_manual_2groups_improved->GetNbinsX());

      mass_ttbar_manual_2groups_improved->SetLineColor(colorlist[B]);
      // if (B == 0) mass_ttbar_manual_2groups_improved->SetLineColor(kRed);
      // if (B == 1) mass_ttbar_manual_2groups_improved->SetLineColor(kBlue);
      // if (B == 2) mass_ttbar_manual_2groups_improved->SetLineColor(kRed);
      // if (B == 3) mass_ttbar_manual_2groups_improved->SetLineColor(kBlue);

      mass_ttbar_manual_2groups_improved->Draw("SAMES");
      leg_mass_2groups_improved_ttbar->AddEntry(mass_ttbar_manual_2groups_improved, "#beta = " + (TString)ss.str());
      mass_ttbar_manual_2groups_improved->Write();

      if (mass_ttbar_manual_2groups_improved->GetMaximum() > max_val_ttbar_2groups_improved) {
        max_val_ttbar_2groups_improved = mass_ttbar_manual_2groups_improved->GetMaximum();
        // mass_ttbar_manual_2groups_improved->SetMaximum(1.2*max_val_ttbar_2groups_improved);
        // mass_2groups_improved_ttbar_compare->Update();
      }
    }

    mass_ttbar_antikt->SetMaximum(1.2*max_val_ttbar_2groups_improved);
    leg_mass_2groups_improved_ttbar->Draw();
    mass_2groups_improved_ttbar_compare->Write();
    mass_2groups_improved_ttbar_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_mass_2groups_improved_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");

    TCanvas *mass_2groups_improved_dijets_compare = new TCanvas("mass_2groups_improved_dijets_compare", "mass_2groups_improved Comparison", 600, 600);
    mass_2groups_improved_dijets_compare->cd();
    mass_dijets_antikt->SetLineColor(kBlack);
    mass_dijets_antikt->SetTitle("Mass of QCD Jets");
    mass_dijets_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_dijets_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_dijets_antikt->Draw();
    mass_dijets_antikt_6jets->Draw("SAMES");

    TLegend *leg_mass_2groups_improved_dijets = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_mass_2groups_improved_dijets->SetLineColor(0);
    leg_mass_2groups_improved_dijets->SetFillColor(0);
    leg_mass_2groups_improved_dijets->AddEntry(mass_dijets_antikt, "ak_{T} (2-jet)", "L");
    leg_mass_2groups_improved_dijets->AddEntry(mass_dijets_antikt_6jets, "ak_{T} (6-jet)", "L");
    // double max_val_dijets_2groups_improved = (mass_dijets_antikt->GetMaximum() > mass_dijets_antikt_6jets->GetMaximum()) ? mass_dijets_antikt->GetMaximum() : mass_dijets_antikt_6jets->GetMaximum();
    double max_val_dijets_2groups_improved = max_val_dijets_akt;

    qcd_efficiency_mass_2groups_improved[0][i_pt] = mass_dijets_antikt->Integral((double)160/500*mass_dijets_antikt->GetNbinsX(), (double)240/500*mass_dijets_antikt->GetNbinsX());
    qcd_efficiency_mass_2groups_improved[1][i_pt] = mass_dijets_antikt_6jets->Integral((double)160/500*mass_dijets_antikt_6jets->GetNbinsX(), (double)240/500*mass_dijets_antikt_6jets->GetNbinsX());

    double mass_dijets_manual_2groups_improved_scale;
    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_dijets_manual_2groups_improved = (TH1*)mass_dijets_manual_2groups_improved_hists.At(B);

      mass_dijets_manual_2groups_improved_scale = 1/mass_dijets_manual_2groups_improved->Integral(0, mass_dijets_manual_2groups_improved->GetNbinsX() + 1);
      mass_dijets_manual_2groups_improved->Scale(mass_dijets_manual_2groups_improved_scale);

      mass_dijets_manual_2groups_improved->SetStats(0);

      qcd_efficiency_mass_2groups_improved[B+2][i_pt] = mass_dijets_manual_2groups_improved->Integral((double)160/500*mass_dijets_manual_2groups_improved->GetNbinsX(), (double)240/500*mass_dijets_manual_2groups_improved->GetNbinsX());

      mass_dijets_manual_2groups_improved->SetLineColor(colorlist[B]);
      // if (B == 0) mass_dijets_manual_2groups_improved->SetLineColor(kRed);
      // if (B == 1) mass_dijets_manual_2groups_improved->SetLineColor(kBlue);
      // if (B == 2) mass_dijets_manual_2groups_improved->SetLineColor(kRed);
      // if (B == 3) mass_dijets_manual_2groups_improved->SetLineColor(kBlue);

      mass_dijets_manual_2groups_improved->Draw("SAMES");
      leg_mass_2groups_improved_dijets->AddEntry(mass_dijets_manual_2groups_improved, "#beta = " + (TString)ss.str());
      mass_dijets_manual_2groups_improved->Write();

      if (mass_dijets_manual_2groups_improved->GetMaximum() > max_val_dijets_2groups_improved) {
        max_val_dijets_2groups_improved = mass_dijets_manual_2groups_improved->GetMaximum();
        // mass_dijets_manual_2groups_improved->SetMaximum(1.2*max_val_dijets_2groups_improved);
        // mass_2groups_improved_dijets_compare->Update();
      }
    }

    mass_dijets_antikt->SetMaximum(1.2*max_val_dijets_2groups_improved);
    leg_mass_2groups_improved_dijets->Draw();
    mass_2groups_improved_dijets_compare->Write();
    mass_2groups_improved_dijets_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_mass_2groups_improved_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");

    TCanvas *minmass_2groups_ttbar_compare = new TCanvas("minmass_2groups_ttbar_compare", "minmass_2groups Comparison", 600, 600);
    minmass_2groups_ttbar_compare->cd();

    // double minmass_ttbar_antikt_2groups_scale = 1/total_ttbar_jets;
    double minmass_ttbar_antikt_2groups_scale = 1/minmass_ttbar_antikt_2groups->Integral(0, minmass_ttbar_antikt_2groups->GetNbinsX() + 1);
    minmass_ttbar_antikt_2groups->Scale(minmass_ttbar_antikt_2groups_scale);
    minmass_ttbar_antikt_2groups->SetTitle("minmass of ttbar jets");
    minmass_ttbar_antikt_2groups->SetLineColor(8);
    minmass_ttbar_antikt_2groups->Draw();

    TLegend *leg_minmass_2groups_ttbar = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg_minmass_2groups_ttbar->SetLineColor(0);
    leg_minmass_2groups_ttbar->SetFillColor(0);
    leg_minmass_2groups_ttbar->AddEntry(minmass_ttbar_antikt_2groups, "ak_{T} (2 groups)", "L");
    double max_val_ttbar_2groups_minmass = minmass_ttbar_antikt_2groups->GetMaximum();

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* minmass_ttbar_manual_2groups = (TH1*)minmass_ttbar_manual_2groups_hists.At(B);
      TH1* minmass_dijets_manual_2groups = (TH1*)minmass_dijets_manual_2groups_hists.At(B);

      double minmass_ttbar_manual_2groups_scale = 1/minmass_ttbar_manual_2groups->Integral(0, minmass_ttbar_manual_2groups->GetNbinsX() + 1);
      double minmass_dijets_manual_2groups_scale = 1/minmass_dijets_manual_2groups->Integral(0, minmass_dijets_manual_2groups->GetNbinsX() + 1);

      // double minmass_ttbar_manual_2groups_scale = 1/total_ttbar_jets;
      // double minmass_dijets_manual_2groups_scale = 1/total_dijets_jets;
      minmass_ttbar_manual_2groups->Scale(minmass_ttbar_manual_2groups_scale);
      minmass_dijets_manual_2groups->Scale(minmass_dijets_manual_2groups_scale);

      minmass_ttbar_manual_2groups->SetStats(0);
      minmass_dijets_manual_2groups->SetStats(0);

      minmass_ttbar_manual_2groups->SetLineColor(colorlist[B]);
      // if (B == 0) minmass_ttbar_manual_2groups->SetLineColor(kRed);
      // if (B == 1) minmass_ttbar_manual_2groups->SetLineColor(kBlue);
      // if (B == 2) minmass_ttbar_manual_2groups->SetLineColor(kRed);
      // if (B == 3) minmass_ttbar_manual_2groups->SetLineColor(kBlue);

      minmass_ttbar_manual_2groups->Draw("SAMES");
      leg_minmass_2groups_ttbar->AddEntry(minmass_ttbar_manual_2groups, "#beta = " + (TString)ss.str());
      minmass_ttbar_manual_2groups->Write();

      if (minmass_ttbar_manual_2groups->GetMaximum() > max_val_ttbar_2groups_minmass) {
        max_val_ttbar_2groups_minmass = minmass_ttbar_manual_2groups->GetMaximum();
        // minmass_ttbar_manual_2groups->SetMaximum(1.2*max_val_ttbar_2groups);
        // minmass_2groups_ttbar_compare->Update();
      }
    }

    minmass_ttbar_antikt_2groups->SetMaximum(1.2*max_val_ttbar_2groups_minmass);
    leg_minmass_2groups_ttbar->Draw();
    minmass_2groups_ttbar_compare->Write();
    minmass_2groups_ttbar_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_minmass_2groups_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *minmass_2groups_dijets_compare = new TCanvas("minmass_2groups_dijets_compare", "minmass_2groups Comparison", 600, 600);
    minmass_2groups_dijets_compare->cd();

    double minmass_dijets_antikt_2groups_scale = 1/minmass_dijets_antikt_2groups->Integral(0, minmass_dijets_antikt_2groups->GetNbinsX() + 1);
    minmass_dijets_antikt_2groups->Scale(minmass_dijets_antikt_2groups_scale);
    minmass_dijets_antikt_2groups->SetTitle("minmass of QCD Jets");
    minmass_dijets_antikt_2groups->SetLineColor(8);
    minmass_dijets_antikt_2groups->Draw();

    TLegend *leg_minmass_2groups_dijets = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg_minmass_2groups_dijets->SetLineColor(0);
    leg_minmass_2groups_dijets->SetFillColor(0);
    leg_minmass_2groups_dijets->AddEntry(minmass_dijets_antikt_2groups, "ak_{T} (2 groups)", "L");
    // double max_val_dijets_2groups = (minmass_dijets_antikt->GetMaximum() > minmass_dijets_antikt_6jets->GetMaximum()) ? minmass_dijets_antikt->GetMaximum() : minmass_dijets_antikt_6jets->GetMaximum();
    double max_val_dijets_2groups_minmass = minmass_ttbar_antikt_2groups->GetMaximum();

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* minmass_dijets_manual_2groups = (TH1*)minmass_dijets_manual_2groups_hists.At(B);

      double minmass_dijets_manual_2groups_scale = 1/minmass_dijets_manual_2groups->Integral(0, minmass_dijets_manual_2groups->GetNbinsX() + 1);
      minmass_dijets_manual_2groups->Scale(minmass_dijets_manual_2groups_scale);

      minmass_dijets_manual_2groups->SetStats(0);

      minmass_dijets_manual_2groups->SetLineColor(colorlist[B]);
      // if (B == 0) minmass_dijets_manual_2groups->SetLineColor(kRed);
      // if (B == 1) minmass_dijets_manual_2groups->SetLineColor(kBlue);
      // if (B == 2) minmass_dijets_manual_2groups->SetLineColor(kRed);
      // if (B == 3) minmass_dijets_manual_2groups->SetLineColor(kBlue);

      minmass_dijets_manual_2groups->Draw("SAMES");
      leg_minmass_2groups_dijets->AddEntry(minmass_dijets_manual_2groups, "#beta = " + (TString)ss.str());
      minmass_dijets_manual_2groups->Write();


      if (minmass_dijets_manual_2groups->GetMaximum() > max_val_dijets_2groups_minmass) {
        max_val_dijets_2groups_minmass = minmass_dijets_manual_2groups->GetMaximum();
        // minmass_dijets_manual_2groups->SetMaximum(1.2*max_val_dijets_2groups);
        // minmass_2groups_dijets_compare->Update();
      }
    }

    minmass_dijets_antikt_2groups->SetMaximum(1.2*max_val_dijets_2groups_minmass);
    leg_minmass_2groups_dijets->Draw();
    minmass_2groups_dijets_compare->Write();
    minmass_2groups_dijets_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_minmass_2groups_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *minmass_2groups_improved_ttbar_compare = new TCanvas("minmass_2groups_improved_ttbar_compare", "minmass_2groups_improved Comparison", 600, 600);
    minmass_2groups_improved_ttbar_compare->cd();

    // double minmass_ttbar_antikt_2groups_improved_scale = 1/total_ttbar_jets;
    minmass_ttbar_antikt_2groups->SetTitle("minmass of ttbar jets");
    minmass_ttbar_antikt_2groups->SetLineColor(8);
    minmass_ttbar_antikt_2groups->Draw();

    TLegend *leg_minmass_2groups_improved_ttbar = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg_minmass_2groups_improved_ttbar->SetLineColor(0);
    leg_minmass_2groups_improved_ttbar->SetFillColor(0);
    leg_minmass_2groups_improved_ttbar->AddEntry(minmass_ttbar_antikt_2groups, "ak_{T} (2 groups)", "L");
    double max_val_ttbar_2groups_improved_minmass = minmass_ttbar_antikt_2groups->GetMaximum();

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* minmass_ttbar_manual_2groups_improved = (TH1*)minmass_ttbar_manual_2groups_improved_hists.At(B);
      TH1* minmass_dijets_manual_2groups_improved = (TH1*)minmass_dijets_manual_2groups_improved_hists.At(B);

      double minmass_ttbar_manual_2groups_improved_scale = 1/minmass_ttbar_manual_2groups_improved->Integral(0, minmass_ttbar_manual_2groups_improved->GetNbinsX() + 1);
      double minmass_dijets_manual_2groups_improved_scale = 1/minmass_dijets_manual_2groups_improved->Integral(0, minmass_dijets_manual_2groups_improved->GetNbinsX() + 1);

      // double minmass_ttbar_manual_2groups_improved_scale = 1/total_ttbar_jets;
      // double minmass_dijets_manual_2groups_improved_scale = 1/total_dijets_jets;
      minmass_ttbar_manual_2groups_improved->Scale(minmass_ttbar_manual_2groups_improved_scale);
      minmass_dijets_manual_2groups_improved->Scale(minmass_dijets_manual_2groups_improved_scale);

      minmass_ttbar_manual_2groups_improved->SetStats(0);
      minmass_dijets_manual_2groups_improved->SetStats(0);

      minmass_ttbar_manual_2groups_improved->SetLineColor(colorlist[B]);
      // if (B == 0) minmass_ttbar_manual_2groups_improved->SetLineColor(kRed);
      // if (B == 1) minmass_ttbar_manual_2groups_improved->SetLineColor(kBlue);
      // if (B == 2) minmass_ttbar_manual_2groups_improved->SetLineColor(kRed);
      // if (B == 3) minmass_ttbar_manual_2groups_improved->SetLineColor(kBlue);

      minmass_ttbar_manual_2groups_improved->Draw("SAMES");
      leg_minmass_2groups_improved_ttbar->AddEntry(minmass_ttbar_manual_2groups_improved, "#beta = " + (TString)ss.str());
      minmass_ttbar_manual_2groups_improved->Write();

      if (minmass_ttbar_manual_2groups_improved->GetMaximum() > max_val_ttbar_2groups_improved_minmass) {
        max_val_ttbar_2groups_improved_minmass = minmass_ttbar_manual_2groups_improved->GetMaximum();
        // minmass_ttbar_manual_2groups_improved->SetMaximum(1.2*max_val_ttbar_2groups_improved);
        // minmass_2groups_improved_ttbar_compare->Update();
      }
    }

    minmass_ttbar_antikt_2groups->SetMaximum(1.2*max_val_ttbar_2groups_improved_minmass);
    leg_minmass_2groups_improved_ttbar->Draw();
    minmass_2groups_improved_ttbar_compare->Write();
    minmass_2groups_improved_ttbar_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_minmass_2groups_improved_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *minmass_2groups_improved_dijets_compare = new TCanvas("minmass_2groups_improved_dijets_compare", "minmass_2groups_improved Comparison", 600, 600);
    minmass_2groups_improved_dijets_compare->cd();

    minmass_dijets_antikt_2groups->SetTitle("minmass of QCD Jets");
    minmass_dijets_antikt_2groups->SetLineColor(8);
    minmass_dijets_antikt_2groups->Draw();

    TLegend *leg_minmass_2groups_improved_dijets = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg_minmass_2groups_improved_dijets->SetLineColor(0);
    leg_minmass_2groups_improved_dijets->SetFillColor(0);
    leg_minmass_2groups_improved_dijets->AddEntry(minmass_dijets_antikt_2groups, "ak_{T} (2 groups)", "L");
    // double max_val_dijets_2groups_improved = (minmass_dijets_antikt->GetMaximum() > minmass_dijets_antikt_6jets->GetMaximum()) ? minmass_dijets_antikt->GetMaximum() : minmass_dijets_antikt_6jets->GetMaximum();
    double max_val_dijets_2groups_improved_minmass = minmass_ttbar_antikt_2groups->GetMaximum();

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* minmass_dijets_manual_2groups_improved = (TH1*)minmass_dijets_manual_2groups_improved_hists.At(B);

      double minmass_dijets_manual_2groups_improved_scale = 1/minmass_dijets_manual_2groups_improved->Integral(0, minmass_dijets_manual_2groups_improved->GetNbinsX() + 1);
      minmass_dijets_manual_2groups_improved->Scale(minmass_dijets_manual_2groups_improved_scale);

      minmass_dijets_manual_2groups_improved->SetStats(0);

      minmass_dijets_manual_2groups_improved->SetLineColor(colorlist[B]);
      // if (B == 0) minmass_dijets_manual_2groups_improved->SetLineColor(kRed);
      // if (B == 1) minmass_dijets_manual_2groups_improved->SetLineColor(kBlue);
      // if (B == 2) minmass_dijets_manual_2groups_improved->SetLineColor(kRed);
      // if (B == 3) minmass_dijets_manual_2groups_improved->SetLineColor(kBlue);

      minmass_dijets_manual_2groups_improved->Draw("SAMES");
      leg_minmass_2groups_improved_dijets->AddEntry(minmass_dijets_manual_2groups_improved, "#beta = " + (TString)ss.str());
      minmass_dijets_manual_2groups_improved->Write();


      if (minmass_dijets_manual_2groups_improved->GetMaximum() > max_val_dijets_2groups_improved_minmass) {
        max_val_dijets_2groups_improved_minmass = minmass_dijets_manual_2groups_improved->GetMaximum();
        // minmass_dijets_manual_2groups_improved->SetMaximum(1.2*max_val_dijets_2groups_improved);
        // minmass_2groups_improved_dijets_compare->Update();
      }
    }

    minmass_dijets_antikt_2groups->SetMaximum(1.2*max_val_dijets_2groups_improved_minmass);
    leg_minmass_2groups_improved_dijets->Draw();
    minmass_2groups_improved_dijets_compare->Write();
    minmass_2groups_improved_dijets_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_minmass_2groups_improved_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *mass_2groups_Wmasscut_ttbar_compare = new TCanvas("mass_2groups_Wmasscut_ttbar_compare", "mass_2groups_Wmasscut Comparison", 600, 600);
    mass_2groups_Wmasscut_ttbar_compare->cd();
    mass_ttbar_antikt->SetLineColor(kBlack);
    mass_ttbar_antikt->SetTitle("Mass of Top Jets");
    mass_ttbar_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_ttbar_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_ttbar_antikt->Draw();

    // double mass_ttbar_antikt_2groups_Wmasscut_scale = 1/mass_ttbar_antikt_2groups->Integral(0, mass_ttbar_antikt_2groups->GetNbinsX() + 1);
    // double mass_ttbar_antikt_2groups_Wmasscut_scale = 1/total_ttbar_jets;
    mass_ttbar_antikt_2groups_Wmasscut->Scale(mass_ttbar_antikt_2groups_scale);
    mass_ttbar_antikt_2groups_Wmasscut->SetLineColor(8);
    mass_ttbar_antikt_2groups_Wmasscut->Draw("SAMES");

    TLegend *leg_mass_2groups_Wmasscut_ttbar = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_mass_2groups_Wmasscut_ttbar->SetLineColor(0);
    leg_mass_2groups_Wmasscut_ttbar->SetFillColor(0);
    leg_mass_2groups_Wmasscut_ttbar->AddEntry(mass_ttbar_antikt, "ak_{T} (2-jet)", "L");
    leg_mass_2groups_Wmasscut_ttbar->AddEntry(mass_ttbar_antikt_2groups_Wmasscut, "ak_{T} (2 groups)", "L");
    double max_val_ttbar_2groups_Wmasscut = (mass_ttbar_antikt->GetMaximum() > mass_ttbar_antikt_2groups_Wmasscut->GetMaximum()) ? mass_ttbar_antikt->GetMaximum() : mass_ttbar_antikt_2groups_Wmasscut->GetMaximum();

    top_efficiency_mass_2groups_Wmasscut[0][i_pt] = mass_ttbar_antikt->Integral((double)160/500*mass_ttbar_antikt->GetNbinsX(), (double)240/500*mass_ttbar_antikt->GetNbinsX());
    top_efficiency_mass_2groups_Wmasscut[1][i_pt] = mass_ttbar_antikt_2groups_Wmasscut->Integral((double)160/500*mass_ttbar_antikt_2groups_Wmasscut->GetNbinsX(), (double)240/500*mass_ttbar_antikt_2groups_Wmasscut->GetNbinsX());

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_ttbar_manual_2groups_Wmasscut = (TH1*)mass_ttbar_manual_2groups_Wmasscut_hists.At(B);
      // TH1* mass_dijets_manual_2groups_Wmasscut = (TH1*)mass_dijets_manual_2groups_Wmasscut_hists.At(B);
      // TH1* mass_ttbar_manual_2groups = (TH1*)mass_ttbar_manual_2groups_hists.At(B);
      // TH1* mass_dijets_manual_2groups = (TH1*)mass_dijets_manual_2groups_hists.At(B);

      // double mass_ttbar_manual_2groups_Wmasscut_scale = 1/mass_ttbar_manual_2groups->Integral(0, mass_ttbar_manual_2groups->GetNbinsX() + 1);
      // double mass_dijets_manual_2groups_Wmasscut_scale = 1/mass_dijets_manual_2groups->Integral(0, mass_dijets_manual_2groups->GetNbinsX() + 1);
      double mass_ttbar_manual_2groups_Wmasscut_scale = mass_ttbar_manual_2groups_scale;

      // double mass_ttbar_manual_2groups_Wmasscut_scale = 1/total_ttbar_jets;
      // double mass_dijets_manual_2groups_Wmasscut_scale = 1/total_dijets_jets;
      mass_ttbar_manual_2groups_Wmasscut->Scale(mass_ttbar_manual_2groups_Wmasscut_scale);
      // mass_dijets_manual_2groups_Wmasscut->Scale(mass_dijets_manual_2groups_Wmasscut_scale);

      mass_ttbar_manual_2groups_Wmasscut->SetStats(0);
      // mass_dijets_manual_2groups_Wmasscut->SetStats(0);

      top_efficiency_mass_2groups_Wmasscut[B+2][i_pt] = mass_ttbar_manual_2groups_Wmasscut->Integral((double)160/500*mass_ttbar_manual_2groups_Wmasscut->GetNbinsX(), (double)240/500*mass_ttbar_manual_2groups_Wmasscut->GetNbinsX());

      mass_ttbar_manual_2groups_Wmasscut->SetLineColor(colorlist[B]);
      // if (B == 0) mass_ttbar_manual_2groups_Wmasscut->SetLineColor(kRed);
      // if (B == 1) mass_ttbar_manual_2groups_Wmasscut->SetLineColor(kBlue);
      // if (B == 2) mass_ttbar_manual_2groups_Wmasscut->SetLineColor(kRed);
      // if (B == 3) mass_ttbar_manual_2groups_Wmasscut->SetLineColor(kBlue);

      mass_ttbar_manual_2groups_Wmasscut->Draw("SAMES");
      leg_mass_2groups_Wmasscut_ttbar->AddEntry(mass_ttbar_manual_2groups_Wmasscut, "#beta = " + (TString)ss.str());
      mass_ttbar_manual_2groups_Wmasscut->Write();

      if (mass_ttbar_manual_2groups_Wmasscut->GetMaximum() > max_val_ttbar_2groups_Wmasscut) {
        max_val_ttbar_2groups_Wmasscut = mass_ttbar_manual_2groups_Wmasscut->GetMaximum();
        // mass_ttbar_manual_2groups_Wmasscut->SetMaximum(1.2*max_val_ttbar_2groups_Wmasscut);
        // mass_2groups_Wmasscut_ttbar_compare->Update();
      }
    }

    mass_ttbar_antikt->SetMaximum(1.2*max_val_ttbar_2groups_Wmasscut);
    leg_mass_2groups_Wmasscut_ttbar->Draw();
    mass_2groups_Wmasscut_ttbar_compare->Write();
    mass_2groups_Wmasscut_ttbar_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_mass_2groups_Wmasscut_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");

    TCanvas *mass_2groups_Wmasscut_dijets_compare = new TCanvas("mass_2groups_Wmasscut_dijets_compare", "mass_2groups_Wmasscut Comparison", 600, 600);
    mass_2groups_Wmasscut_dijets_compare->cd();
    mass_dijets_antikt->SetLineColor(kBlack);
    mass_dijets_antikt->SetTitle("Mass of QCD Jets");
    mass_dijets_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_dijets_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_dijets_antikt->Draw();

    // double mass_dijets_antikt_2groups_Wmasscut_scale = 1/mass_dijets_antikt_2groups->Integral(0, mass_dijets_antikt_2groups->GetNbinsX() + 1);
    // double mass_dijets_antikt_2groups_Wmasscut_scale = 1/total_dijets_jets;
    mass_dijets_antikt_2groups_Wmasscut->Scale(mass_dijets_antikt_2groups_scale);
    mass_dijets_antikt_2groups_Wmasscut->SetLineColor(8);
    mass_dijets_antikt_2groups_Wmasscut->Draw("SAMES");

    TLegend *leg_mass_2groups_Wmasscut_dijets = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_mass_2groups_Wmasscut_dijets->SetLineColor(0);
    leg_mass_2groups_Wmasscut_dijets->SetFillColor(0);
    leg_mass_2groups_Wmasscut_dijets->AddEntry(mass_dijets_antikt, "ak_{T} (2-jet)", "L");
    leg_mass_2groups_Wmasscut_dijets->AddEntry(mass_dijets_antikt_2groups_Wmasscut, "ak_{T} (2 groups)", "L");
    double max_val_dijets_2groups_Wmasscut = (mass_dijets_antikt->GetMaximum() > mass_dijets_antikt_2groups_Wmasscut->GetMaximum()) ? mass_dijets_antikt->GetMaximum() : mass_dijets_antikt_2groups_Wmasscut->GetMaximum();
    // double max_val_dijets_2groups_Wmasscut = max_val_dijets_akt;

    qcd_efficiency_mass_2groups_Wmasscut[0][i_pt] = mass_dijets_antikt->Integral((double)160/500*mass_dijets_antikt->GetNbinsX(), (double)240/500*mass_dijets_antikt->GetNbinsX());
    qcd_efficiency_mass_2groups_Wmasscut[1][i_pt] = mass_dijets_antikt_2groups_Wmasscut->Integral((double)160/500*mass_dijets_antikt_2groups_Wmasscut->GetNbinsX(), (double)240/500*mass_dijets_antikt_2groups_Wmasscut->GetNbinsX());

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_dijets_manual_2groups_Wmasscut = (TH1*)mass_dijets_manual_2groups_Wmasscut_hists.At(B);
      // TH1* mass_dijets_manual_2groups = (TH1*)mass_dijets_manual_2groups_hists.At(B);

      // double mass_dijets_manual_2groups_Wmasscut_scale = 1/mass_dijets_manual_2groups->Integral(0, mass_dijets_manual_2groups->GetNbinsX() + 1);
      // double mass_dijets_manual_2groups_Wmasscut_scale = 1/total_dijets_jets;
      double mass_dijets_manual_2groups_Wmasscut_scale = mass_dijets_manual_2groups_scale;
      mass_dijets_manual_2groups_Wmasscut->Scale(mass_dijets_manual_2groups_Wmasscut_scale);

      mass_dijets_manual_2groups_Wmasscut->SetStats(0);

      qcd_efficiency_mass_2groups_Wmasscut[B+2][i_pt] = mass_dijets_manual_2groups_Wmasscut->Integral((double)160/500*mass_dijets_manual_2groups_Wmasscut->GetNbinsX(), (double)240/500*mass_dijets_manual_2groups_Wmasscut->GetNbinsX());

      mass_dijets_manual_2groups_Wmasscut->SetLineColor(colorlist[B]);
      // if (B == 0) mass_dijets_manual_2groups_Wmasscut->SetLineColor(kRed);
      // if (B == 1) mass_dijets_manual_2groups_Wmasscut->SetLineColor(kBlue);
      // if (B == 2) mass_dijets_manual_2groups_Wmasscut->SetLineColor(kRed);
      // if (B == 3) mass_dijets_manual_2groups_Wmasscut->SetLineColor(kBlue);

      mass_dijets_manual_2groups_Wmasscut->Draw("SAMES");
      leg_mass_2groups_Wmasscut_dijets->AddEntry(mass_dijets_manual_2groups_Wmasscut, "#beta = " + (TString)ss.str());
      mass_dijets_manual_2groups_Wmasscut->Write();

      if (mass_dijets_manual_2groups_Wmasscut->GetMaximum() > max_val_dijets_2groups_Wmasscut) {
        max_val_dijets_2groups_Wmasscut = mass_dijets_manual_2groups_Wmasscut->GetMaximum();
        // mass_dijets_manual_2groups_Wmasscut->SetMaximum(1.2*max_val_dijets_2groups_Wmasscut);
        // mass_2groups_Wmasscut_dijets_compare->Update();
      }
    }

    mass_dijets_antikt->SetMaximum(1.2*max_val_dijets_2groups_Wmasscut);
    leg_mass_2groups_Wmasscut_dijets->Draw();
    mass_2groups_Wmasscut_dijets_compare->Write();
    mass_2groups_Wmasscut_dijets_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_mass_2groups_Wmasscut_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *mass_2groups_improved_Wmasscut_ttbar_compare = new TCanvas("mass_2groups_improved_Wmasscut_ttbar_compare", "mass_2groups_improved_Wmasscut Comparison", 600, 600);
    mass_2groups_improved_Wmasscut_ttbar_compare->cd();
    mass_ttbar_antikt->SetLineColor(kBlack);
    mass_ttbar_antikt->SetTitle("Mass of Top Jets");
    mass_ttbar_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_ttbar_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_ttbar_antikt->Draw();

    // double mass_ttbar_antikt_2groups_improved_Wmasscut_scale = 1/mass_ttbar_antikt_2groups_improved->Integral(0, mass_ttbar_antikt_2groups_improved->GetNbinsX() + 1);
    // double mass_ttbar_antikt_2groups_improved_Wmasscut_scale = 1/total_ttbar_jets;
    mass_ttbar_antikt_2groups_Wmasscut->Scale(mass_ttbar_antikt_2groups_scale);
    mass_ttbar_antikt_2groups_Wmasscut->SetLineColor(8);
    mass_ttbar_antikt_2groups_Wmasscut->Draw("SAMES");

    TLegend *leg_mass_2groups_improved_Wmasscut_ttbar = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_mass_2groups_improved_Wmasscut_ttbar->SetLineColor(0);
    leg_mass_2groups_improved_Wmasscut_ttbar->SetFillColor(0);
    leg_mass_2groups_improved_Wmasscut_ttbar->AddEntry(mass_ttbar_antikt, "ak_{T} (2-jet)", "L");
    leg_mass_2groups_improved_Wmasscut_ttbar->AddEntry(mass_ttbar_antikt_2groups_Wmasscut, "ak_{T} (2 groups)", "L");
    double max_val_ttbar_2groups_improved_Wmasscut = (mass_ttbar_antikt->GetMaximum() > mass_ttbar_antikt_2groups_Wmasscut->GetMaximum()) ? mass_ttbar_antikt->GetMaximum() : mass_ttbar_antikt_2groups_Wmasscut->GetMaximum();

    top_efficiency_mass_2groups_improved_Wmasscut[0][i_pt] = mass_ttbar_antikt->Integral((double)160/500*mass_ttbar_antikt->GetNbinsX(), (double)240/500*mass_ttbar_antikt->GetNbinsX());
    top_efficiency_mass_2groups_improved_Wmasscut[1][i_pt] = mass_ttbar_antikt_2groups_Wmasscut->Integral((double)160/500*mass_ttbar_antikt_2groups_Wmasscut->GetNbinsX(), (double)240/500*mass_ttbar_antikt_2groups_Wmasscut->GetNbinsX());

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_ttbar_manual_2groups_improved_Wmasscut = (TH1*)mass_ttbar_manual_2groups_improved_Wmasscut_hists.At(B);
      // TH1* mass_dijets_manual_2groups_improved_Wmasscut = (TH1*)mass_dijets_manual_2groups_improved_Wmasscut_hists.At(B);
      // TH1* mass_ttbar_manual_2groups_improved = (TH1*)mass_ttbar_manual_2groups_improved_hists.At(B);
      // TH1* mass_dijets_manual_2groups_improved = (TH1*)mass_dijets_manual_2groups_improved_hists.At(B);

      // double mass_ttbar_manual_2groups_improved_Wmasscut_scale = 1/mass_ttbar_manual_2groups_improved->Integral(0, mass_ttbar_manual_2groups_improved->GetNbinsX() + 1);
      // double mass_dijets_manual_2groups_improved_Wmasscut_scale = 1/mass_dijets_manual_2groups_improved->Integral(0, mass_dijets_manual_2groups_improved->GetNbinsX() + 1);
      double mass_ttbar_manual_2groups_improved_Wmasscut_scale = mass_ttbar_manual_2groups_improved_scale;

      // double mass_ttbar_manual_2groups_improved_Wmasscut_scale = 1/total_ttbar_jets;
      // double mass_dijets_manual_2groups_improved_Wmasscut_scale = 1/total_dijets_jets;
      mass_ttbar_manual_2groups_improved_Wmasscut->Scale(mass_ttbar_manual_2groups_improved_Wmasscut_scale);
      // mass_dijets_manual_2groups_improved_Wmasscut->Scale(mass_dijets_manual_2groups_improved_Wmasscut_scale);

      mass_ttbar_manual_2groups_improved_Wmasscut->SetStats(0);
      // mass_dijets_manual_2groups_improved_Wmasscut->SetStats(0);

      top_efficiency_mass_2groups_improved_Wmasscut[B+2][i_pt] = mass_ttbar_manual_2groups_improved_Wmasscut->Integral((double)160/500*mass_ttbar_manual_2groups_improved_Wmasscut->GetNbinsX(), (double)240/500*mass_ttbar_manual_2groups_improved_Wmasscut->GetNbinsX());

      mass_ttbar_manual_2groups_improved_Wmasscut->SetLineColor(colorlist[B]);
      // if (B == 0) mass_ttbar_manual_2groups_improved_Wmasscut->SetLineColor(kRed);
      // if (B == 1) mass_ttbar_manual_2groups_improved_Wmasscut->SetLineColor(kBlue);
      // if (B == 2) mass_ttbar_manual_2groups_improved_Wmasscut->SetLineColor(kRed);
      // if (B == 3) mass_ttbar_manual_2groups_improved_Wmasscut->SetLineColor(kBlue);

      mass_ttbar_manual_2groups_improved_Wmasscut->Draw("SAMES");
      leg_mass_2groups_improved_Wmasscut_ttbar->AddEntry(mass_ttbar_manual_2groups_improved_Wmasscut, "#beta = " + (TString)ss.str());
      mass_ttbar_manual_2groups_improved_Wmasscut->Write();

      if (mass_ttbar_manual_2groups_improved_Wmasscut->GetMaximum() > max_val_ttbar_2groups_improved_Wmasscut) {
        max_val_ttbar_2groups_improved_Wmasscut = mass_ttbar_manual_2groups_improved_Wmasscut->GetMaximum();
        // mass_ttbar_manual_2groups_improved_Wmasscut->SetMaximum(1.2*max_val_ttbar_2groups_improved_Wmasscut);
        // mass_2groups_improved_Wmasscut_ttbar_compare->Update();
      }
    }

    mass_ttbar_antikt->SetMaximum(1.2*max_val_ttbar_2groups_improved_Wmasscut);
    leg_mass_2groups_improved_Wmasscut_ttbar->Draw();
    mass_2groups_improved_Wmasscut_ttbar_compare->Write();
    mass_2groups_improved_Wmasscut_ttbar_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_mass_2groups_improved_Wmasscut_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");

    TCanvas *mass_2groups_improved_Wmasscut_dijets_compare = new TCanvas("mass_2groups_improved_Wmasscut_dijets_compare", "mass_2groups_improved_Wmasscut Comparison", 600, 600);
    mass_2groups_improved_Wmasscut_dijets_compare->cd();
    mass_dijets_antikt->SetLineColor(kBlack);
    mass_dijets_antikt->SetTitle("Mass of QCD Jets");
    mass_dijets_antikt->GetXaxis()->SetTitle("m (GeV)");
    mass_dijets_antikt->GetYaxis()->SetTitle("Relative Occurrence");
    mass_dijets_antikt->Draw();

    // double mass_dijets_antikt_2groups_improved_Wmasscut_scale = 1/mass_dijets_antikt_2groups_improved->Integral(0, mass_dijets_antikt_2groups_improved->GetNbinsX() + 1);
    // double mass_dijets_antikt_2groups_improved_Wmasscut_scale = 1/total_dijets_jets;
    mass_dijets_antikt_2groups_Wmasscut->Scale(mass_dijets_antikt_2groups_scale);
    mass_dijets_antikt_2groups_Wmasscut->SetLineColor(8);
    mass_dijets_antikt_2groups_Wmasscut->Draw("SAMES");

    TLegend *leg_mass_2groups_improved_Wmasscut_dijets = new TLegend(0.6, 0.6, 0.88, 0.88);
    leg_mass_2groups_improved_Wmasscut_dijets->SetLineColor(0);
    leg_mass_2groups_improved_Wmasscut_dijets->SetFillColor(0);
    leg_mass_2groups_improved_Wmasscut_dijets->AddEntry(mass_dijets_antikt, "ak_{T} (2-jet)", "L");
    leg_mass_2groups_improved_Wmasscut_dijets->AddEntry(mass_dijets_antikt_2groups_Wmasscut, "ak_{T} (2 groups)", "L");
    double max_val_dijets_2groups_improved_Wmasscut = (mass_dijets_antikt->GetMaximum() > mass_dijets_antikt_2groups_Wmasscut->GetMaximum()) ? mass_dijets_antikt->GetMaximum() : mass_dijets_antikt_2groups_Wmasscut->GetMaximum();
    // double max_val_dijets_2groups_improved_Wmasscut = max_val_dijets_akt;

    qcd_efficiency_mass_2groups_improved_Wmasscut[0][i_pt] = mass_dijets_antikt->Integral((double)160/500*mass_dijets_antikt->GetNbinsX(), (double)240/500*mass_dijets_antikt->GetNbinsX());
    qcd_efficiency_mass_2groups_improved_Wmasscut[1][i_pt] = mass_dijets_antikt_2groups_Wmasscut->Integral((double)160/500*mass_dijets_antikt_2groups_Wmasscut->GetNbinsX(), (double)240/500*mass_dijets_antikt_2groups_Wmasscut->GetNbinsX());

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_dijets_manual_2groups_improved_Wmasscut = (TH1*)mass_dijets_manual_2groups_improved_Wmasscut_hists.At(B);
      // TH1* mass_dijets_manual_2groups_improved = (TH1*)mass_dijets_manual_2groups_improved_hists.At(B);

      // double mass_dijets_manual_2groups_improved_Wmasscut_scale = 1/mass_dijets_manual_2groups_improved->Integral(0, mass_dijets_manual_2groups_improved->GetNbinsX() + 1);
      // double mass_dijets_manual_2groups_improved_Wmasscut_scale = 1/total_dijets_jets;
      double mass_dijets_manual_2groups_improved_Wmasscut_scale = mass_dijets_manual_2groups_improved_scale;
      mass_dijets_manual_2groups_improved_Wmasscut->Scale(mass_dijets_manual_2groups_improved_Wmasscut_scale);

      mass_dijets_manual_2groups_improved_Wmasscut->SetStats(0);

      qcd_efficiency_mass_2groups_improved_Wmasscut[B+2][i_pt] = mass_dijets_manual_2groups_improved_Wmasscut->Integral((double)160/500*mass_dijets_manual_2groups_improved_Wmasscut->GetNbinsX(), (double)240/500*mass_dijets_manual_2groups_improved_Wmasscut->GetNbinsX());

      mass_dijets_manual_2groups_improved_Wmasscut->SetLineColor(colorlist[B]);
      // if (B == 0) mass_dijets_manual_2groups_improved_Wmasscut->SetLineColor(kRed);
      // if (B == 1) mass_dijets_manual_2groups_improved_Wmasscut->SetLineColor(kBlue);
      // if (B == 2) mass_dijets_manual_2groups_improved_Wmasscut->SetLineColor(kRed);
      // if (B == 3) mass_dijets_manual_2groups_improved_Wmasscut->SetLineColor(kBlue);

      mass_dijets_manual_2groups_improved_Wmasscut->Draw("SAMES");
      leg_mass_2groups_improved_Wmasscut_dijets->AddEntry(mass_dijets_manual_2groups_improved_Wmasscut, "#beta = " + (TString)ss.str());
      mass_dijets_manual_2groups_improved_Wmasscut->Write();

      if (mass_dijets_manual_2groups_improved_Wmasscut->GetMaximum() > max_val_dijets_2groups_improved_Wmasscut) {
        max_val_dijets_2groups_improved_Wmasscut = mass_dijets_manual_2groups_improved_Wmasscut->GetMaximum();
        // mass_dijets_manual_2groups_improved_Wmasscut->SetMaximum(1.2*max_val_dijets_2groups_improved_Wmasscut);
        // mass_2groups_improved_Wmasscut_dijets_compare->Update();
      }
    }

    mass_dijets_antikt->SetMaximum(1.2*max_val_dijets_2groups_improved_Wmasscut);
    leg_mass_2groups_improved_Wmasscut_dijets->Draw();
    mass_2groups_improved_Wmasscut_dijets_compare->Write();
    mass_2groups_improved_Wmasscut_dijets_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_mass_2groups_improved_Wmasscut_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *mass_ttbar_3jet4jet_2groups_compare = new TCanvas("mass_ttbar_3jet4jet_2groups_compare", "mass_ttbar_3jet4jet_2groups_compare", 600, 600);
    mass_ttbar_3jet4jet_2groups_compare->cd();

    double mass_ttbar_antikt_3jet4jet_2groups_compare_scale = (double)1/nEvent;
    mass_ttbar_antikt_3jet4jet_2groups_compare->Scale(mass_ttbar_antikt_3jet4jet_2groups_compare_scale);
    mass_ttbar_antikt_3jet4jet_2groups_compare->SetMarkerColor(kBlack);
    mass_ttbar_antikt_3jet4jet_2groups_compare->SetTitle("3-jet/4-jet min mass comparison (p_{T} > " + (TString)ss_pt.str() + ")");
    mass_ttbar_antikt_3jet4jet_2groups_compare->GetXaxis()->SetTitle("3-jet mass");
    mass_ttbar_antikt_3jet4jet_2groups_compare->GetYaxis()->SetTitle("4-jet min mass");
    mass_ttbar_antikt_3jet4jet_2groups_compare->SetMinimum(0.0);
    mass_ttbar_antikt_3jet4jet_2groups_compare->Draw("box");
    double max_val_mass_ttbar_3jet4jet_2groups = mass_ttbar_antikt_3jet4jet_2groups_compare->GetMaximum();

    TLegend *leg_mass_ttbar_3jet4jet_2groups = new TLegend(0.14, 0.65, 0.35, 0.88);
    leg_mass_ttbar_3jet4jet_2groups->SetFillColor(kWhite);
    leg_mass_ttbar_3jet4jet_2groups->SetLineColor(kWhite);
    leg_mass_ttbar_3jet4jet_2groups->AddEntry(mass_ttbar_antikt_3jet4jet_2groups_compare, "ak_{T}", "L");

    top_efficiency_mass_3jet4jet_2groups[0][i_pt] = mass_ttbar_antikt->Integral((double)160/500*mass_ttbar_antikt->GetNbinsX(), (double)240/500*mass_ttbar_antikt->GetNbinsX());
    top_efficiency_mass_3jet4jet_2groups[1][i_pt] = (double)total_ttbar_events[1]*mass_ttbar_antikt_2groups_scale;

    for (int B = 0; B < n_betas; B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_ttbar_manual_3jet4jet_2groups_compare = (TH1*)mass_ttbar_manual_3jet4jet_2groups_compare_hists.At(B);
      double mass_ttbar_manual_3jet4jet_2groups_compare_scale = (double)1/nEvent;
      mass_ttbar_manual_3jet4jet_2groups_compare->Scale(mass_ttbar_manual_3jet4jet_2groups_compare_scale);

      mass_ttbar_manual_3jet4jet_2groups_compare->SetStats(0);
      mass_ttbar_manual_3jet4jet_2groups_compare->Write();

      mass_ttbar_manual_3jet4jet_2groups_compare->SetMarkerColor(colorlist[B]);
      mass_ttbar_manual_3jet4jet_2groups_compare->SetLineColor(colorlist[B]);

      top_efficiency_mass_3jet4jet_2groups[B+2][i_pt] = (double)total_ttbar_events[B+2]*mass_ttbar_manual_2groups_scale;

      mass_ttbar_manual_3jet4jet_2groups_compare->Draw("box SAMES");
      leg_mass_ttbar_3jet4jet_2groups->AddEntry(mass_ttbar_manual_3jet4jet_2groups_compare, "#beta = " + (TString)ss.str(), "L");

      if (mass_ttbar_manual_3jet4jet_2groups_compare->GetMaximum() > max_val_mass_ttbar_3jet4jet_2groups) max_val_mass_ttbar_3jet4jet_2groups = mass_ttbar_manual_3jet4jet_2groups_compare->GetMaximum();
    }

    TLine *line1 = new TLine();
    line1->SetX1(160);
    line1->SetVertical(true);
    line1->SetY1(0);
    line1->SetY2(500);
    line1->SetLineColor(kRed);
    line1->SetLineWidth(3);
    line1->SetLineStyle(7);
    line1->Draw("SAMES");

    TLine *line2 = new TLine();
    line2->SetX1(240);
    line2->SetVertical(true);
    line2->SetY1(0);
    line2->SetY2(500);
    line2->SetLineColor(kRed);
    line2->SetLineWidth(3);
    line2->SetLineStyle(7);
    line2->Draw("SAMES");

    TLine *line3 = new TLine();
    line3->SetY1(160);
    line3->SetHorizontal(true);
    line3->SetX1(0);
    line3->SetX2(500);
    line3->SetLineColor(kRed);
    line3->SetLineWidth(3);
    line3->SetLineStyle(7);
    line3->Draw("SAMES");

    TLine *line4 = new TLine();
    line4->SetY1(240);
    line4->SetHorizontal(true);
    line4->SetX1(0);
    line4->SetX2(500);
    line4->SetLineColor(kRed);
    line4->SetLineWidth(3);
    line4->SetLineStyle(7);
    line4->Draw("SAMES");

    mass_ttbar_antikt_3jet4jet_2groups_compare->SetMaximum(1.2*max_val_mass_ttbar_3jet4jet_2groups);
    leg_mass_ttbar_3jet4jet_2groups->Draw("SAMES");
    mass_ttbar_3jet4jet_2groups_compare->Write();
    mass_ttbar_3jet4jet_2groups_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_mass_3jet4jet_2groups_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *mass_dijets_3jet4jet_2groups_compare = new TCanvas("mass_dijets_3jet4jet_2groups_compare", "mass_dijets_3jet4jet_2groups_compare", 600, 600);
    mass_dijets_3jet4jet_2groups_compare->cd();

    double mass_dijets_antikt_3jet4jet_2groups_compare_scale = (double)1/nEvent;
    mass_dijets_antikt_3jet4jet_2groups_compare->Scale(mass_dijets_antikt_3jet4jet_2groups_compare_scale);
    mass_dijets_antikt_3jet4jet_2groups_compare->SetMarkerColor(kBlack);
    mass_dijets_antikt_3jet4jet_2groups_compare->SetTitle("3-jet/4-jet min mass comparison (p_{T} > " + (TString)ss_pt.str() + ")");
    mass_dijets_antikt_3jet4jet_2groups_compare->GetXaxis()->SetTitle("3-jet mass");
    mass_dijets_antikt_3jet4jet_2groups_compare->GetYaxis()->SetTitle("4-jet min mass");
    mass_dijets_antikt_3jet4jet_2groups_compare->SetMinimum(0.0);
    mass_dijets_antikt_3jet4jet_2groups_compare->Draw("box");
    double max_val_mass_dijets_3jet4jet_2groups = mass_dijets_antikt_3jet4jet_2groups_compare->GetMaximum();

    TLegend *leg_mass_dijets_3jet4jet_2groups = new TLegend(0.14, 0.65, 0.35, 0.88);
    leg_mass_dijets_3jet4jet_2groups->SetFillColor(kWhite);
    leg_mass_dijets_3jet4jet_2groups->SetLineColor(kWhite);
    leg_mass_dijets_3jet4jet_2groups->AddEntry(mass_dijets_antikt_3jet4jet_2groups_compare, "ak_{T}", "L");

    qcd_efficiency_mass_3jet4jet_2groups[0][i_pt] = mass_dijets_antikt->Integral((double)160/500*mass_dijets_antikt->GetNbinsX(), (double)240/500*mass_dijets_antikt->GetNbinsX());
    qcd_efficiency_mass_3jet4jet_2groups[1][i_pt] = (double)total_dijets_events[1]*mass_dijets_antikt_2groups_scale;

    for (int B = 0; B < n_betas; B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_dijets_manual_3jet4jet_2groups_compare = (TH1*)mass_dijets_manual_3jet4jet_2groups_compare_hists.At(B);
      double mass_dijets_manual_3jet4jet_2groups_compare_scale = (double)1/nEvent;
      mass_dijets_manual_3jet4jet_2groups_compare->Scale(mass_dijets_manual_3jet4jet_2groups_compare_scale);

      mass_dijets_manual_3jet4jet_2groups_compare->SetStats(0);
      mass_dijets_manual_3jet4jet_2groups_compare->Write();

      mass_dijets_manual_3jet4jet_2groups_compare->SetMarkerColor(colorlist[B]);
      mass_dijets_manual_3jet4jet_2groups_compare->SetLineColor(colorlist[B]);

      qcd_efficiency_mass_3jet4jet_2groups[B+2][i_pt] = (double)total_dijets_events[B+2]*mass_dijets_manual_2groups_scale;

      mass_dijets_manual_3jet4jet_2groups_compare->Draw("box SAMES");
      leg_mass_dijets_3jet4jet_2groups->AddEntry(mass_dijets_manual_3jet4jet_2groups_compare, "#beta = " + (TString)ss.str(), "L");

      if (mass_dijets_manual_3jet4jet_2groups_compare->GetMaximum() > max_val_mass_dijets_3jet4jet_2groups) max_val_mass_dijets_3jet4jet_2groups = mass_dijets_manual_3jet4jet_2groups_compare->GetMaximum();
    }

    line1->Draw("SAMES");
    line2->Draw("SAMES");
    line3->Draw("SAMES");
    line4->Draw("SAMES");

    mass_dijets_antikt_3jet4jet_2groups_compare->SetMaximum(1.2*max_val_mass_dijets_3jet4jet_2groups);
    leg_mass_dijets_3jet4jet_2groups->Draw("SAMES");
    mass_dijets_3jet4jet_2groups_compare->Write();
    mass_dijets_3jet4jet_2groups_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_mass_3jet4jet_2groups_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *minmass_ttbar_2jet3jet_2groups_compare = new TCanvas("minmass_ttbar_2jet3jet_2groups_compare", "minmass_ttbar_2jet3jet_2groups_compare", 600, 600);
    minmass_ttbar_2jet3jet_2groups_compare->cd();

    double minmass_ttbar_antikt_2jet3jet_2groups_compare_scale = (double)1/nEvent;
    minmass_ttbar_antikt_2jet3jet_2groups_compare->Scale(minmass_ttbar_antikt_2jet3jet_2groups_compare_scale);
    minmass_ttbar_antikt_2jet3jet_2groups_compare->SetMarkerColor(kBlack);
    minmass_ttbar_antikt_2jet3jet_2groups_compare->SetTitle("3-jet/4-jet min minmass comparison (p_{T} > " + (TString)ss_pt.str() + ")");
    minmass_ttbar_antikt_2jet3jet_2groups_compare->GetXaxis()->SetTitle("3-jet minmass");
    minmass_ttbar_antikt_2jet3jet_2groups_compare->GetYaxis()->SetTitle("4-jet min minmass");
    minmass_ttbar_antikt_2jet3jet_2groups_compare->SetMinimum(0.0);
    minmass_ttbar_antikt_2jet3jet_2groups_compare->Draw("box");
    double max_val_minmass_ttbar_2jet3jet_2groups = minmass_ttbar_antikt_2jet3jet_2groups_compare->GetMaximum();

    TLegend *leg_minmass_ttbar_2jet3jet_2groups = new TLegend(0.14, 0.65, 0.35, 0.88);
    leg_minmass_ttbar_2jet3jet_2groups->SetFillColor(kWhite);
    leg_minmass_ttbar_2jet3jet_2groups->SetLineColor(kWhite);
    leg_minmass_ttbar_2jet3jet_2groups->AddEntry(minmass_ttbar_antikt_2jet3jet_2groups_compare, "ak_{T}", "L");

    for (int B = 0; B < n_betas; B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* minmass_ttbar_manual_2jet3jet_2groups_compare = (TH1*)minmass_ttbar_manual_2jet3jet_2groups_compare_hists.At(B);
      double minmass_ttbar_manual_2jet3jet_2groups_compare_scale = (double)1/nEvent;
      minmass_ttbar_manual_2jet3jet_2groups_compare->Scale(minmass_ttbar_manual_2jet3jet_2groups_compare_scale);

      minmass_ttbar_manual_2jet3jet_2groups_compare->SetStats(0);
      minmass_ttbar_manual_2jet3jet_2groups_compare->Write();

      minmass_ttbar_manual_2jet3jet_2groups_compare->SetMarkerColor(colorlist[B]);
      minmass_ttbar_manual_2jet3jet_2groups_compare->SetLineColor(colorlist[B]);

      minmass_ttbar_manual_2jet3jet_2groups_compare->Draw("box SAMES");
      leg_minmass_ttbar_2jet3jet_2groups->AddEntry(minmass_ttbar_manual_2jet3jet_2groups_compare, "#beta = " + (TString)ss.str(), "L");

      if (minmass_ttbar_manual_2jet3jet_2groups_compare->GetMaximum() > max_val_minmass_ttbar_2jet3jet_2groups) max_val_minmass_ttbar_2jet3jet_2groups = minmass_ttbar_manual_2jet3jet_2groups_compare->GetMaximum();
    }

    TLine *line5 = new TLine();
    line5->SetX1(50);
    line5->SetVertical(true);
    line5->SetY1(0);
    line5->SetY2(200);
    line5->SetLineColor(kRed);
    line5->SetLineWidth(3);
    line5->SetLineStyle(7);
    line5->Draw("SAMES");

    TLine *line6 = new TLine();
    line6->SetY1(50);
    line6->SetHorizontal(true);
    line6->SetX1(0);
    line6->SetX2(200);
    line6->SetLineColor(kRed);
    line6->SetLineWidth(3);
    line6->SetLineStyle(7);
    line6->Draw("SAMES");

    minmass_ttbar_antikt_2jet3jet_2groups_compare->SetMaximum(1.2*max_val_minmass_ttbar_2jet3jet_2groups);
    leg_minmass_ttbar_2jet3jet_2groups->Draw("SAMES");
    minmass_ttbar_2jet3jet_2groups_compare->Write();
    minmass_ttbar_2jet3jet_2groups_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_minmass_2jet3jet_2groups_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *minmass_dijets_2jet3jet_2groups_compare = new TCanvas("minmass_dijets_2jet3jet_2groups_compare", "minmass_dijets_2jet3jet_2groups_compare", 600, 600);
    minmass_dijets_2jet3jet_2groups_compare->cd();

    double minmass_dijets_antikt_2jet3jet_2groups_compare_scale = (double)1/nEvent;
    minmass_dijets_antikt_2jet3jet_2groups_compare->Scale(minmass_dijets_antikt_2jet3jet_2groups_compare_scale);
    minmass_dijets_antikt_2jet3jet_2groups_compare->SetMarkerColor(kBlack);
    minmass_dijets_antikt_2jet3jet_2groups_compare->SetTitle("3-jet/4-jet min minmass comparison (p_{T} > " + (TString)ss_pt.str() + ")");
    minmass_dijets_antikt_2jet3jet_2groups_compare->GetXaxis()->SetTitle("3-jet minmass");
    minmass_dijets_antikt_2jet3jet_2groups_compare->GetYaxis()->SetTitle("4-jet min minmass");
    minmass_dijets_antikt_2jet3jet_2groups_compare->SetMinimum(0.0);
    minmass_dijets_antikt_2jet3jet_2groups_compare->Draw("box");
    double max_val_minmass_dijets_2jet3jet_2groups = minmass_dijets_antikt_2jet3jet_2groups_compare->GetMaximum();

    TLegend *leg_minmass_dijets_2jet3jet_2groups = new TLegend(0.14, 0.65, 0.35, 0.88);
    leg_minmass_dijets_2jet3jet_2groups->SetFillColor(kWhite);
    leg_minmass_dijets_2jet3jet_2groups->SetLineColor(kWhite);
    leg_minmass_dijets_2jet3jet_2groups->AddEntry(minmass_dijets_antikt_2jet3jet_2groups_compare, "ak_{T}", "L");

    for (int B = 0; B < n_betas; B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* minmass_dijets_manual_2jet3jet_2groups_compare = (TH1*)minmass_dijets_manual_2jet3jet_2groups_compare_hists.At(B);
      double minmass_dijets_manual_2jet3jet_2groups_compare_scale = (double)1/nEvent;
      minmass_dijets_manual_2jet3jet_2groups_compare->Scale(minmass_dijets_manual_2jet3jet_2groups_compare_scale);

      minmass_dijets_manual_2jet3jet_2groups_compare->SetStats(0);
      minmass_dijets_manual_2jet3jet_2groups_compare->Write();

      minmass_dijets_manual_2jet3jet_2groups_compare->SetMarkerColor(colorlist[B]);
      minmass_dijets_manual_2jet3jet_2groups_compare->SetLineColor(colorlist[B]);

      minmass_dijets_manual_2jet3jet_2groups_compare->Draw("box SAMES");
      leg_minmass_dijets_2jet3jet_2groups->AddEntry(minmass_dijets_manual_2jet3jet_2groups_compare, "#beta = " + (TString)ss.str(), "L");

      if (minmass_dijets_manual_2jet3jet_2groups_compare->GetMaximum() > max_val_minmass_dijets_2jet3jet_2groups) max_val_minmass_dijets_2jet3jet_2groups = minmass_dijets_manual_2jet3jet_2groups_compare->GetMaximum();
    }

    line5->Draw("SAMES");
    line6->Draw("SAMES");

    minmass_dijets_antikt_2jet3jet_2groups_compare->SetMaximum(1.2*max_val_minmass_dijets_2jet3jet_2groups);
    leg_minmass_dijets_2jet3jet_2groups->Draw("SAMES");
    minmass_dijets_2jet3jet_2groups_compare->Write();
    minmass_dijets_2jet3jet_2groups_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_minmass_2jet3jet_2groups_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *mass_ttbar_3jet4jet_2groups_Wmasscut_compare = new TCanvas("mass_ttbar_3jet4jet_2groups_Wmasscut_compare", "mass_ttbar_3jet4jet_2groups_Wmasscut_compare", 600, 600);
    mass_ttbar_3jet4jet_2groups_Wmasscut_compare->cd();

    double mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare_scale = (double)1/nEvent;
    mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare->Scale(mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare_scale);
    mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare->SetMarkerColor(kBlack);
    mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare->SetTitle("3-jet/4-jet min mass comparison (p_{T} > " + (TString)ss_pt.str() + ")");
    mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare->GetXaxis()->SetTitle("3-jet mass");
    mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare->GetYaxis()->SetTitle("4-jet min mass");
    mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare->SetMinimum(0.0);
    mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare->Draw("box");
    double max_val_mass_ttbar_3jet4jet_2groups_Wmasscut = mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare->GetMaximum();

    TLegend *leg_mass_ttbar_3jet4jet_2groups_Wmasscut = new TLegend(0.14, 0.65, 0.35, 0.88);
    leg_mass_ttbar_3jet4jet_2groups_Wmasscut->SetFillColor(kWhite);
    leg_mass_ttbar_3jet4jet_2groups_Wmasscut->SetLineColor(kWhite);
    leg_mass_ttbar_3jet4jet_2groups_Wmasscut->AddEntry(mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare, "ak_{T}", "L");

    top_efficiency_mass_3jet4jet_2groups_Wmasscut[0][i_pt] = mass_ttbar_antikt->Integral((double)160/500*mass_ttbar_antikt->GetNbinsX(), (double)240/500*mass_ttbar_antikt->GetNbinsX());
    top_efficiency_mass_3jet4jet_2groups_Wmasscut[1][i_pt] = (double)total_ttbar_Wmasscut_events[1]*mass_ttbar_antikt_2groups_scale;

    for (int B = 0; B < n_betas; B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare = (TH1*)mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare_hists.At(B);
      double mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare_scale = (double)1/nEvent;
      mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare->Scale(mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare_scale);

      mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare->SetStats(0);
      mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare->Write();

      mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare->SetMarkerColor(colorlist[B]);
      mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare->SetLineColor(colorlist[B]);

      top_efficiency_mass_3jet4jet_2groups_Wmasscut[B+2][i_pt] = (double)total_ttbar_Wmasscut_events[B+2]*mass_ttbar_manual_2groups_scale;

      mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare->Draw("box SAMES");
      leg_mass_ttbar_3jet4jet_2groups_Wmasscut->AddEntry(mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare, "#beta = " + (TString)ss.str(), "L");

      if (mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare->GetMaximum() > max_val_mass_ttbar_3jet4jet_2groups_Wmasscut) max_val_mass_ttbar_3jet4jet_2groups_Wmasscut = mass_ttbar_manual_3jet4jet_2groups_Wmasscut_compare->GetMaximum();
    }

    line1->Draw("SAMES");
    line2->Draw("SAMES");
    line3->Draw("SAMES");
    line4->Draw("SAMES");

    mass_ttbar_antikt_3jet4jet_2groups_Wmasscut_compare->SetMaximum(1.2*max_val_mass_ttbar_3jet4jet_2groups_Wmasscut);
    leg_mass_ttbar_3jet4jet_2groups_Wmasscut->Draw("SAMES");
    mass_ttbar_3jet4jet_2groups_Wmasscut_compare->Write();
    mass_ttbar_3jet4jet_2groups_Wmasscut_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_mass_3jet4jet_2groups_Wmasscut_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *mass_dijets_3jet4jet_2groups_Wmasscut_compare = new TCanvas("mass_dijets_3jet4jet_2groups_Wmasscut_compare", "mass_dijets_3jet4jet_2groups_Wmasscut_compare", 600, 600);
    mass_dijets_3jet4jet_2groups_Wmasscut_compare->cd();

    double mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare_scale = (double)1/nEvent;
    mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare->Scale(mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare_scale);
    mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare->SetMarkerColor(kBlack);
    mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare->SetTitle("3-jet/4-jet min mass comparison (p_{T} > " + (TString)ss_pt.str() + ")");
    mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare->GetXaxis()->SetTitle("3-jet mass");
    mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare->GetYaxis()->SetTitle("4-jet min mass");
    mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare->SetMinimum(0.0);
    mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare->Draw("box");
    double max_val_mass_dijets_3jet4jet_2groups_Wmasscut = mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare->GetMaximum();

    TLegend *leg_mass_dijets_3jet4jet_2groups_Wmasscut = new TLegend(0.14, 0.65, 0.35, 0.88);
    leg_mass_dijets_3jet4jet_2groups_Wmasscut->SetFillColor(kWhite);
    leg_mass_dijets_3jet4jet_2groups_Wmasscut->SetLineColor(kWhite);
    leg_mass_dijets_3jet4jet_2groups_Wmasscut->AddEntry(mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare, "ak_{T}", "L");

    qcd_efficiency_mass_3jet4jet_2groups_Wmasscut[0][i_pt] = mass_dijets_antikt->Integral((double)160/500*mass_dijets_antikt->GetNbinsX(), (double)240/500*mass_dijets_antikt->GetNbinsX());
    qcd_efficiency_mass_3jet4jet_2groups_Wmasscut[1][i_pt] = (double)total_dijets_Wmasscut_events[1]*mass_dijets_antikt_2groups_scale;

    for (int B = 0; B < n_betas; B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare = (TH1*)mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare_hists.At(B);
      double mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare_scale = (double)1/nEvent;
      mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare->Scale(mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare_scale);

      mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare->SetStats(0);
      mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare->Write();

      mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare->SetMarkerColor(colorlist[B]);
      mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare->SetLineColor(colorlist[B]);

      qcd_efficiency_mass_3jet4jet_2groups_Wmasscut[B+2][i_pt] = (double)total_dijets_Wmasscut_events[B+2]*mass_dijets_manual_2groups_scale;

      mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare->Draw("box SAMES");
      leg_mass_dijets_3jet4jet_2groups_Wmasscut->AddEntry(mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare, "#beta = " + (TString)ss.str(), "L");

      if (mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare->GetMaximum() > max_val_mass_dijets_3jet4jet_2groups_Wmasscut) max_val_mass_dijets_3jet4jet_2groups_Wmasscut = mass_dijets_manual_3jet4jet_2groups_Wmasscut_compare->GetMaximum();
    }

    line1->Draw("SAMES");
    line2->Draw("SAMES");
    line3->Draw("SAMES");
    line4->Draw("SAMES");

    mass_dijets_antikt_3jet4jet_2groups_Wmasscut_compare->SetMaximum(1.2*max_val_mass_dijets_3jet4jet_2groups_Wmasscut);
    leg_mass_dijets_3jet4jet_2groups_Wmasscut->Draw("SAMES");
    mass_dijets_3jet4jet_2groups_Wmasscut_compare->Write();
    mass_dijets_3jet4jet_2groups_Wmasscut_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_mass_3jet4jet_2groups_Wmasscut_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *volatility_ttbar_compare = new TCanvas("volatility_ttbar_compare", "volatility Comparison", 600, 600);
    volatility_ttbar_compare->cd();

    TLegend *leg_volatility_ttbar = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg_volatility_ttbar->SetLineColor(0);
    leg_volatility_ttbar->SetFillColor(0);

    double max_val_ttbar_vol = 0;

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* volatility_ttbar_manual = (TH1*)volatility_ttbar_manual_hists.At(B);
      double volatility_ttbar_manual_scale = 1/volatility_ttbar_manual->Integral(0, volatility_ttbar_manual->GetNbinsX() + 1);
      volatility_ttbar_manual->Scale(volatility_ttbar_manual_scale);

      volatility_ttbar_manual->SetStats(0);

      volatility_ttbar_manual->SetLineColor(colorlist[B]);
        // if (B == 0) rawvolatility_ttbar_manual_2jets->SetLineColor(kRed);
        // if (B == 1) rawvolatility_ttbar_manual_2jets->SetLineColor(kBlue);
        // if (B == 2) rawvolatility_ttbar_manual_2jets->SetLineColor(kRed);
        // if (B == 3) rawvolatility_ttbar_manual_2jets->SetLineColor(kBlue);

      if (B == 0) volatility_ttbar_manual->Draw();
      else volatility_ttbar_manual->Draw("SAMES");
      leg_volatility_ttbar->AddEntry(volatility_ttbar_manual, "#beta = " + (TString)ss.str());
      volatility_ttbar_manual->Write();

      if (volatility_ttbar_manual->GetMaximum() > max_val_ttbar_vol) {
        max_val_ttbar_vol = volatility_ttbar_manual->GetMaximum();
          // volatility_ttbar_manual->SetMaximum(1.2*max_val_ttbar_vol);
          // volatility_ttbar_compare->Update();
      }
      if (B == 1) volatility_ttbar_manual->SetMaximum(1.2*max_val_ttbar_vol);
    }

    leg_volatility_ttbar->Draw();
    volatility_ttbar_compare->Write();
    volatility_ttbar_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_volatility_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *volatility_dijets_compare = new TCanvas("volatility_dijets_compare", "volatility Comparison", 600, 600);
    volatility_dijets_compare->cd();

    TLegend *leg_volatility_dijets = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg_volatility_dijets->SetLineColor(0);
    leg_volatility_dijets->SetFillColor(0);

    double max_val_dijets_vol = 0;

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* volatility_dijets_manual = (TH1*)volatility_dijets_manual_hists.At(B);

      double volatility_dijets_manual_scale = 1/volatility_dijets_manual->Integral(0, volatility_dijets_manual->GetNbinsX() + 1);
      volatility_dijets_manual->Scale(volatility_dijets_manual_scale);

      volatility_dijets_manual->SetStats(0);

      volatility_dijets_manual->SetLineColor(colorlist[B]);
        // if (B == 0) volatility_dijets_manual->SetLineColor(kRed);
        // if (B == 1) volatility_dijets_manual->SetLineColor(kBlue);
        // if (B == 2) volatility_dijets_manual->SetLineColor(kRed);
        // if (B == 3) volatility_dijets_manual->SetLineColor(kBlue);

      if (B == 0) volatility_dijets_manual->Draw();
      else volatility_dijets_manual->Draw("SAMES");
      leg_volatility_dijets->AddEntry(volatility_dijets_manual, "#beta = " + (TString)ss.str());
      volatility_dijets_manual->Write();

      if (volatility_dijets_manual->GetMaximum() > max_val_dijets_vol) {
        max_val_dijets_vol = volatility_dijets_manual->GetMaximum();
          // volatility_dijets_manual->SetMaximum(1.2*max_val_dijets);
          // volatility_dijets_compare->Update();
      }
      if (B == 1) volatility_dijets_manual->SetMaximum(1.2*max_val_dijets_vol);
    }

    leg_volatility_dijets->Draw();
    volatility_dijets_compare->Write();
    volatility_dijets_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_volatility_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *max_njet_ttbar_compare = new TCanvas("max_njet_ttbar_compare", "max_njet Comparison", 600, 600);
    max_njet_ttbar_compare->cd();
    TLegend *leg_max_njet_ttbar = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg_max_njet_ttbar->SetLineColor(0);
    leg_max_njet_ttbar->SetFillColor(0);
    double max_val_ttbar_njet = 0;

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* max_njet_ttbar_manual = (TH1*)max_njet_ttbar_manual_hists.At(B);
      TH1* max_njet_dijets_manual = (TH1*)max_njet_dijets_manual_hists.At(B);

      double max_njet_ttbar_manual_scale = 1/max_njet_ttbar_manual->Integral(0, max_njet_ttbar_manual->GetNbinsX() + 1);
      double max_njet_dijets_manual_scale = 1/max_njet_dijets_manual->Integral(0, max_njet_dijets_manual->GetNbinsX() + 1);
      max_njet_ttbar_manual->Scale(max_njet_ttbar_manual_scale);
      max_njet_dijets_manual->Scale(max_njet_dijets_manual_scale);

      max_njet_ttbar_manual->SetStats(0);
      max_njet_dijets_manual->SetStats(0);

      max_njet_ttbar_manual->SetTitle("Maximum N-jettiness of Top Jets");
      max_njet_ttbar_manual->GetXaxis()->SetTitle("N-jettiness");
      max_njet_ttbar_manual->GetYaxis()->SetTitle("Relative Occurrence");

      max_njet_ttbar_manual->SetLineColor(colorlist[B]);
        // if (B == 0) max_njet_ttbar_manual->SetLineColor(kRed);
        // if (B == 1) max_njet_ttbar_manual->SetLineColor(kBlue);
        // if (B == 2) max_njet_ttbar_manual->SetLineColor(kRed);
        // if (B == 3) max_njet_ttbar_manual->SetLineColor(kBlue);

      if (B == 0) max_njet_ttbar_manual->Draw();
      else max_njet_ttbar_manual->Draw("SAMES");
      leg_max_njet_ttbar->AddEntry(max_njet_ttbar_manual, "#beta = " + (TString)ss.str());
      max_njet_ttbar_manual->Write();

      if (max_njet_ttbar_manual->GetMaximum() > max_val_ttbar_njet) {
        max_val_ttbar_njet = max_njet_ttbar_manual->GetMaximum();
        max_njet_ttbar_manual->SetMaximum(1.2*max_val_ttbar_njet);
      }
    }

    leg_max_njet_ttbar->Draw();
    max_njet_ttbar_compare->Write();
    max_njet_ttbar_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_max_njet_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");

    TCanvas *max_njet_dijets_compare = new TCanvas("max_njet_dijets_compare", "max_njet Comparison", 600, 600);
    max_njet_dijets_compare->cd();
    TLegend *leg_max_njet_dijets = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg_max_njet_dijets->SetLineColor(0);
    leg_max_njet_dijets->SetFillColor(0);
    double max_val_dijets_njet = 0;

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* max_njet_dijets_manual = (TH1*)max_njet_dijets_manual_hists.At(B);

      double max_njet_dijets_manual_scale = 1/max_njet_dijets_manual->Integral(0, max_njet_dijets_manual->GetNbinsX() + 1);
      max_njet_dijets_manual->Scale(max_njet_dijets_manual_scale);

      max_njet_dijets_manual->SetStats(0);
      max_njet_dijets_manual->SetTitle("Maximum N-jettiness of QCD Jets");
      max_njet_dijets_manual->GetXaxis()->SetTitle("N-jettiness");
      max_njet_dijets_manual->GetYaxis()->SetTitle("Relative Occurrence");

      max_njet_dijets_manual->SetLineColor(colorlist[B]);
        // if (B == 0) max_njet_dijets_manual->SetLineColor(kRed);
        // if (B == 1) max_njet_dijets_manual->SetLineColor(kBlue);
        // if (B == 2) max_njet_dijets_manual->SetLineColor(kRed);
        // if (B == 3) max_njet_dijets_manual->SetLineColor(kBlue);

      if (B == 0) max_njet_dijets_manual->Draw();
      else max_njet_dijets_manual->Draw("SAMES");
      leg_max_njet_dijets->AddEntry(max_njet_dijets_manual, "#beta = " + (TString)ss.str());
      max_njet_dijets_manual->Write();

      if (max_njet_dijets_manual->GetMaximum() > max_val_dijets_njet) {
        max_val_dijets_njet = max_njet_dijets_manual->GetMaximum();
        max_njet_dijets_manual->SetMaximum(1.2*max_val_dijets_njet);
      }

    }

    leg_max_njet_dijets->Draw();
    max_njet_dijets_compare->Write();
    max_njet_dijets_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_max_njet_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


    TCanvas *max_njet_ttbar_2groups_improved_compare = new TCanvas("max_njet_ttbar_2groups_improved_compare", "max_njet Comparison", 600, 600);
    max_njet_ttbar_2groups_improved_compare->cd();
    TLegend *leg_max_njet_2groups_improved_ttbar = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg_max_njet_2groups_improved_ttbar->SetLineColor(0);
    leg_max_njet_2groups_improved_ttbar->SetFillColor(0);
    double max_val_ttbar_njet_2groups_improved = 0;

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* max_njet_ttbar_manual_2groups_improved = (TH1*)max_njet_ttbar_manual_2groups_improved_hists.At(B);
      TH1* max_njet_dijets_manual_2groups_improved = (TH1*)max_njet_dijets_manual_2groups_improved_hists.At(B);

      double max_njet_ttbar_manual_2groups_improved_scale = 1/max_njet_ttbar_manual_2groups_improved->Integral(0, max_njet_ttbar_manual_2groups_improved->GetNbinsX() + 1);
      double max_njet_dijets_manual_2groups_improved_scale = 1/max_njet_dijets_manual_2groups_improved->Integral(0, max_njet_dijets_manual_2groups_improved->GetNbinsX() + 1);
      max_njet_ttbar_manual_2groups_improved->Scale(max_njet_ttbar_manual_2groups_improved_scale);
      max_njet_dijets_manual_2groups_improved->Scale(max_njet_dijets_manual_2groups_improved_scale);

      max_njet_ttbar_manual_2groups_improved->SetStats(0);
      max_njet_dijets_manual_2groups_improved->SetStats(0);

      max_njet_ttbar_manual_2groups_improved->SetLineColor(colorlist[B]);
        // if (B == 0) max_njet_ttbar_manual_2groups_improved->SetLineColor(kRed);
        // if (B == 1) max_njet_ttbar_manual_2groups_improved->SetLineColor(kBlue);
        // if (B == 2) max_njet_ttbar_manual_2groups_improved->SetLineColor(kRed);
        // if (B == 3) max_njet_ttbar_manual_2groups_improved->SetLineColor(kBlue);

      if (B == 0) max_njet_ttbar_manual_2groups_improved->Draw();
      else max_njet_ttbar_manual_2groups_improved->Draw("SAMES");
      leg_max_njet_2groups_improved_ttbar->AddEntry(max_njet_ttbar_manual_2groups_improved, "#beta = " + (TString)ss.str());
      max_njet_ttbar_manual_2groups_improved->Write();

      if (max_njet_ttbar_manual_2groups_improved->GetMaximum() > max_val_ttbar_njet_2groups_improved) {
        max_val_ttbar_njet_2groups_improved = max_njet_ttbar_manual_2groups_improved->GetMaximum();
        max_njet_ttbar_manual_2groups_improved->SetMaximum(1.2*max_val_ttbar_njet_2groups_improved);
      }
    }

    leg_max_njet_2groups_improved_ttbar->Draw();
    max_njet_ttbar_2groups_improved_compare->Write();
    max_njet_ttbar_2groups_improved_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_ttbar_max_njet_2groups_improved_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");

    TCanvas *max_njet_dijets_2groups_improved_compare = new TCanvas("max_njet_dijets_2groups_improved_compare", "max_njet Comparison", 600, 600);
    max_njet_dijets_2groups_improved_compare->cd();
    TLegend *leg_max_njet_2groups_improved_dijets = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg_max_njet_2groups_improved_dijets->SetLineColor(0);
    leg_max_njet_2groups_improved_dijets->SetFillColor(0);
    double max_val_dijets_njet_2groups_improved = 0;

    for (int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];

      ostringstream ss;
      ss << beta;

      TH1* max_njet_dijets_manual_2groups_improved = (TH1*)max_njet_dijets_manual_2groups_improved_hists.At(B);

      double max_njet_dijets_manual_2groups_improved_scale = 1/max_njet_dijets_manual_2groups_improved->Integral(0, max_njet_dijets_manual_2groups_improved->GetNbinsX() + 1);
      max_njet_dijets_manual_2groups_improved->Scale(max_njet_dijets_manual_2groups_improved_scale);

      max_njet_dijets_manual_2groups_improved->SetStats(0);

      max_njet_dijets_manual_2groups_improved->SetLineColor(colorlist[B]);
        // if (B == 0) max_njet_dijets_manual_2groups_improved->SetLineColor(kRed);
        // if (B == 1) max_njet_dijets_manual_2groups_improved->SetLineColor(kBlue);
        // if (B == 2) max_njet_dijets_manual_2groups_improved->SetLineColor(kRed);
        // if (B == 3) max_njet_dijets_manual_2groups_improved->SetLineColor(kBlue);

      if (B == 0) max_njet_dijets_manual_2groups_improved->Draw();
      else max_njet_dijets_manual_2groups_improved->Draw("SAMES");
      leg_max_njet_2groups_improved_dijets->AddEntry(max_njet_dijets_manual_2groups_improved, "#beta = " + (TString)ss.str());
      max_njet_dijets_manual_2groups_improved->Write();

      if (max_njet_dijets_manual_2groups_improved->GetMaximum() > max_val_dijets_njet_2groups_improved) {
        max_val_dijets_njet_2groups_improved = max_njet_dijets_manual_2groups_improved->GetMaximum();
        max_njet_dijets_manual_2groups_improved->SetMaximum(1.2*max_val_dijets_njet_2groups_improved);
      }

    }

    leg_max_njet_2groups_improved_dijets->Draw();
    max_njet_dijets_2groups_improved_compare->Write();
    max_njet_dijets_2groups_improved_compare->Print("ttbarstudy_2x3_plots/ttbarstudy_dijets_max_njet_2groups_improved_compare_pt" + (TString)ss_pt.str() + ".eps", "eps");


  // cout << endl;
  // cout << "mean manual tau3: " << tau3_ttbar_manual_antikt->GetMean() << endl;
  // cout << "sqrt(mean^2) manual tau3: " << TMath::Sqrt(TMath::Power(tau3_ttbar_manual_antikt->GetRMS(),2) + TMath::Power(tau3_ttbar_manual_antikt->GetMean(),2)) << endl;
  // cout << endl;

  // mass_ratios.push_back(mass_ttbar_manual->Integral(75,110)/total_ttbar_jets);
  // mean_values.push_back(tau3_ttbar_manual_antikt->GetMean());
  // rms_values.push_back(TMath::Sqrt(TMath::Power(tau3_ttbar_manual_antikt->GetRMS(),2) + TMath::Power(tau3_ttbar_manual_antikt->GetMean(),2)));

  // double tau6_ttbar_manual_scale = 1/tau6_ttbar_manual->Integral();
  // double tau6_dijets_manual_scale = 1/tau6_dijets_manual->Integral();
  // double mass_ttbar_manual_scale = 1/mass_ttbar_manual->Integral();
  // double mass_dijets_manual_scale = 1/mass_dijets_manual->Integral();
  // double mass_ttbar_antikt_scale = 1/mass_ttbar_antikt->Integral();

  // tau6_ttbar_manual->Scale(tau6_ttbar_manual_scale);
  // tau6_dijets_manual->Scale(tau6_dijets_manual_scale);
  // mass_ttbar_manual->Scale(mass_ttbar_manual_scale);
  // mass_dijets_manual->Scale(mass_dijets_manual_scale);
  // mass_ttbar_antikt->Scale(mass_ttbar_antikt_scale);


  // delete tau6_ttbar_manual;
  // delete tau6_dijets_manual;
  // delete tau32_ttbar_manual;
  // delete tau32_dijets_manual;
  // delete mass_ttbar_manual;
  // delete mass_dijets_manual;
    delete mass_ttbar_antikt;
    delete mass_dijets_antikt;
      // delete ROC_compare;
  // delete tau3_ttbar_manual_antikt;
  // delete tau3_dijets_manual_antikt;
  // delete mass_ttbar_compare;

  }


  TCanvas *ROC_compare_total_big = new TCanvas("ROC_compare_total_big", "ROC Curve Comparison", 600, 600);
  ROC_compare_total_big->cd();
  ROC_compare_total_big->SetLogy();
  TMultiGraph *ROC_multigraph_total_big = new TMultiGraph("ROC_multigraph_total_big", "ROC Comparison for #tau_{3}/#tau_{2} (240 < p_{T} < 800)");
  TLegend *leg_ROC_total_big = new TLegend(0.14, 0.6, 0.4, 0.88);
  leg_ROC_total_big->SetLineColor(0);
  leg_ROC_total_big->SetFillColor(0);

  TObjArray ROC_tau32_antikt_total_allbeta;
  TObjArray ROC_tau32_manual_6jet_total_allbeta;
  TObjArray ROC_tau32_manual_total_allbeta;
  TObjArray ROC_tau32_manual_2groups_total_allbeta;
  TObjArray ROC_tau32_manual_2groups_improved_total_allbeta;
  TObjArray ROC_tau32_manual_2groups_bigcone_total_allbeta;

  for (int B = 0; B < betalist.size(); B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

      // TCanvas *ROC_compare_total = new TCanvas("ROC_compare_total", "ROC Curve Comparison", 600, 600);
      // ROC_compare_total->cd();
      // ROC_compare_total->SetLogy();
      // TMultiGraph *ROC_multigraph_total = new TMultiGraph("ROC_multigraph_total", "ROC Comparison for #tau_{3}/#tau_{2} (240 < p_{T} < 800, #beta = " + (TString)ss.str() + ")");
      // TLegend *leg_ROC_total = new TLegend(0.14, 0.6, 0.4, 0.88);
      // leg_ROC_total->SetLineColor(0);
      // leg_ROC_total->SetFillColor(0);

    TH1F *tau32_ttbar_antikt_total = (TH1F*)tau32_ttbar_antikt_total_hists.At(B);
    TH1F *tau32_dijets_antikt_total = (TH1F*)tau32_dijets_antikt_total_hists.At(B);

    tau32_ttbar_antikt_total->SetStats(0);
    tau32_dijets_antikt_total->SetStats(0);
    tau32_ttbar_antikt_total->Write();
    tau32_dijets_antikt_total->Write();

    int n_antikt_total = tau32_ttbar_antikt_total->GetSize() - 2;
    double integral_dijets_antikt_total[n_antikt_total], integral_ttbar_antikt_total[n_antikt_total];
    for (int i = 0; i < n_antikt_total; i++) {
      integral_ttbar_antikt_total[i] = tau32_ttbar_antikt_total->Integral(0,i)/total_total_ttbar_jets;
      integral_dijets_antikt_total[i] = tau32_dijets_antikt_total->Integral(0,i)/total_total_dijets_jets;
    }
    TGraph* ROC_tau32_antikt_total = new TGraph(n_antikt_total, integral_ttbar_antikt_total, integral_dijets_antikt_total);
    ROC_tau32_antikt_total->GetXaxis()->SetLimits(0, 1);
    ROC_tau32_antikt_total->GetYaxis()->SetLimits(0, 1);
    if (B == 0) ROC_tau32_antikt_total->SetLineColor(8);
    else ROC_tau32_antikt_total->SetLineColor(kOrange);
    ROC_tau32_antikt_total->SetLineWidth(2);
    ROC_tau32_antikt_total->SetMarkerStyle(5);
    ROC_tau32_antikt_total->SetMarkerSize(2);
    ROC_tau32_antikt_total->Write();
      // ROC_multigraph_total->Add(ROC_tau32_antikt_total);
      // leg_ROC_total->AddEntry(ROC_tau32_antikt_total, "Anti-k_{T} Jets", "L");
    ROC_multigraph_total_big->Add(ROC_tau32_antikt_total);
    leg_ROC_total_big->AddEntry(ROC_tau32_antikt_total, "Anti-k_{T} Jets, #beta = " + (TString)ss.str(), "L");

    ROC_tau32_antikt_total_allbeta.Add(ROC_tau32_antikt_total);


    TH1F *tau32_ttbar_manual_6jet_total = (TH1F*)tau32_ttbar_manual_6jet_total_hists.At(B);
    TH1F *tau32_dijets_manual_6jet_total = (TH1F*)tau32_dijets_manual_6jet_total_hists.At(B);

    tau32_ttbar_manual_6jet_total->SetStats(0);
    tau32_dijets_manual_6jet_total->SetStats(0);
    tau32_ttbar_manual_6jet_total->Write();
    tau32_dijets_manual_6jet_total->Write();

    int n_size = tau32_ttbar_manual_6jet_total->GetSize() - 2;
    double integral_dijets[n_size], integral_ttbar[n_size];
    for (int i = 0; i < n_size; i++) {
      integral_ttbar[i] = tau32_ttbar_manual_6jet_total->Integral(0,i)/total_total_ttbar_jets;
      integral_dijets[i] = tau32_dijets_manual_6jet_total->Integral(0,i)/total_total_dijets_jets;
    }
    TGraph* ROC_tau32_manual_6jet_total = new TGraph(n_size, integral_ttbar, integral_dijets);
    ROC_tau32_manual_6jet_total->GetXaxis()->SetLimits(0, 1);
    ROC_tau32_manual_6jet_total->GetYaxis()->SetLimits(0, 1);
    if (B == 0) ROC_tau32_manual_6jet_total->SetLineColor(kOrange);
    else ROC_tau32_manual_6jet_total->SetLineColor(8);
    ROC_tau32_manual_6jet_total->SetLineWidth(2);
    ROC_tau32_manual_6jet_total->SetMarkerStyle(5);
    ROC_tau32_manual_6jet_total->SetMarkerSize(2);
      // ROC_multigraph_total->Add(ROC_tau32_manual_6jet_total);
      // leg_ROC_total->AddEntry(ROC_tau32_manual_6jet_total, "6-jettiness Jets", "L");
    ROC_multigraph_total_big->Add(ROC_tau32_manual_6jet_total);
    leg_ROC_total_big->AddEntry(ROC_tau32_manual_6jet_total, "6-jettiness Jets, #beta = " + (TString)ss.str(), "L");

    ROC_tau32_manual_6jet_total_allbeta.Add(ROC_tau32_manual_6jet_total);


    TH1F *tau32_ttbar_manual_total = (TH1F*)tau32_ttbar_manual_total_hists.At(B);
    TH1F *tau32_dijets_manual_total = (TH1F*)tau32_dijets_manual_total_hists.At(B);

    tau32_ttbar_manual_total->SetStats(0);
    tau32_dijets_manual_total->SetStats(0);
    tau32_ttbar_manual_total->Write();
    tau32_dijets_manual_total->Write();

    int n_size_6jet = tau32_ttbar_manual_total->GetSize() - 2;
    double integral_dijets_6jet[n_size_6jet], integral_ttbar_6jet[n_size_6jet];
    for (int i = 0; i < n_size_6jet; i++) {
      integral_ttbar_6jet[i] = tau32_ttbar_manual_total->Integral(0,i)/total_total_ttbar_jets;
      integral_dijets_6jet[i] = tau32_dijets_manual_total->Integral(0,i)/total_total_dijets_jets;
    }
    TGraph* ROC_tau32_manual_total = new TGraph(n_size_6jet, integral_ttbar_6jet, integral_dijets_6jet);
    ROC_tau32_manual_total->GetXaxis()->SetLimits(0, 1);
    ROC_tau32_manual_total->GetYaxis()->SetLimits(0, 1);
    if (B == 0) ROC_tau32_manual_total->SetLineColor(kBlue);
    else ROC_tau32_manual_total->SetLineColor(8);
    ROC_tau32_manual_total->SetLineWidth(2);
    ROC_tau32_manual_total->SetMarkerStyle(5);
    ROC_tau32_manual_total->SetMarkerSize(2);
      // ROC_multigraph_total->Add(ROC_tau32_manual_total);
      // leg_ROC_total->AddEntry(ROC_tau32_manual_total, "6-jettiness Jets", "L");
    ROC_multigraph_total_big->Add(ROC_tau32_manual_total);
    leg_ROC_total_big->AddEntry(ROC_tau32_manual_total, "6-jettiness Jets, #beta = " + (TString)ss.str(), "L");

    ROC_tau32_manual_total_allbeta.Add(ROC_tau32_manual_total);


    TH1F *tau32_ttbar_manual_2groups_total = (TH1F*)tau32_ttbar_manual_2groups_total_hists.At(B);
    TH1F *tau32_dijets_manual_2groups_total = (TH1F*)tau32_dijets_manual_2groups_total_hists.At(B);

    tau32_ttbar_manual_2groups_total->SetStats(0);
    tau32_dijets_manual_2groups_total->SetStats(0);
    tau32_ttbar_manual_2groups_total->Write();
    tau32_dijets_manual_2groups_total->Write();

    int n_size_manual_2groups_total = tau32_ttbar_manual_2groups_total->GetSize() - 2;
    double integral_dijets_manual_2groups_total[n_size_manual_2groups_total], integral_ttbar_manual_2groups_total[n_size_manual_2groups_total];
    for (int i = 0; i < n_size_manual_2groups_total; i++) {
      integral_ttbar_manual_2groups_total[i] = tau32_ttbar_manual_2groups_total->Integral(0,i)/total_total_ttbar_jets;
      integral_dijets_manual_2groups_total[i] = tau32_dijets_manual_2groups_total->Integral(0,i)/total_total_dijets_jets;
    }
    TGraph* ROC_tau32_manual_2groups_total = new TGraph(n_size_manual_2groups_total, integral_ttbar_manual_2groups_total, integral_dijets_manual_2groups_total);
    ROC_tau32_manual_2groups_total->GetXaxis()->SetLimits(0, 1);
    ROC_tau32_manual_2groups_total->GetYaxis()->SetLimits(0, 1);
    if (B == 0) ROC_tau32_manual_2groups_total->SetLineColor(kRed);
    else ROC_tau32_manual_2groups_total->SetLineColor(6);
    ROC_tau32_manual_2groups_total->SetLineWidth(2);
    ROC_tau32_manual_2groups_total->SetMarkerStyle(5);
    ROC_tau32_manual_2groups_total->SetMarkerSize(2);
      // ROC_multigraph_total->Add(ROC_tau32_manual_2groups_total);
      // leg_ROC_total->AddEntry(ROC_tau32_manual_2groups_total, "2X3-jettiness Jets", "L");
    ROC_multigraph_total_big->Add(ROC_tau32_manual_2groups_total);
    leg_ROC_total_big->AddEntry(ROC_tau32_manual_2groups_total, "2X3-jettiness Jets, #beta = " + (TString)ss.str(), "L");

    ROC_tau32_manual_2groups_total_allbeta.Add(ROC_tau32_manual_2groups_total);


    TH1F *tau32_ttbar_manual_2groups_improved_total = (TH1F*)tau32_ttbar_manual_2groups_improved_total_hists.At(B);
    TH1F *tau32_dijets_manual_2groups_improved_total = (TH1F*)tau32_dijets_manual_2groups_improved_total_hists.At(B);

    tau32_ttbar_manual_2groups_improved_total->SetStats(0);
    tau32_dijets_manual_2groups_improved_total->SetStats(0);
    tau32_ttbar_manual_2groups_improved_total->Write();
    tau32_dijets_manual_2groups_improved_total->Write();

    int n_size_manual_2groups_improved_total = tau32_ttbar_manual_2groups_improved_total->GetSize() - 2;
    double integral_dijets_manual_2groups_improved_total[n_size_manual_2groups_improved_total], integral_ttbar_manual_2groups_improved_total[n_size_manual_2groups_improved_total];
    for (int i = 0; i < n_size_manual_2groups_improved_total; i++) {
      integral_ttbar_manual_2groups_improved_total[i] = tau32_ttbar_manual_2groups_improved_total->Integral(0,i)/total_total_ttbar_jets;
      integral_dijets_manual_2groups_improved_total[i] = tau32_dijets_manual_2groups_improved_total->Integral(0,i)/total_total_dijets_jets;
    }
    TGraph* ROC_tau32_manual_2groups_improved_total = new TGraph(n_size_manual_2groups_improved_total, integral_ttbar_manual_2groups_improved_total, integral_dijets_manual_2groups_improved_total);
    ROC_tau32_manual_2groups_improved_total->GetXaxis()->SetLimits(0, 1);
    ROC_tau32_manual_2groups_improved_total->GetYaxis()->SetLimits(0, 1);
    if (B == 0) ROC_tau32_manual_2groups_improved_total->SetLineColor(kGreen);
    else ROC_tau32_manual_2groups_improved_total->SetLineColor(6);
    ROC_tau32_manual_2groups_improved_total->SetLineWidth(2);
    ROC_tau32_manual_2groups_improved_total->SetMarkerStyle(5);
    ROC_tau32_manual_2groups_improved_total->SetMarkerSize(2);
      // ROC_multigraph_total->Add(ROC_tau32_manual_2groups_improved_total);
      // leg_ROC_total->AddEntry(ROC_tau32_manual_2groups_improved_total, "2X3-jettiness Jets", "L");
    ROC_multigraph_total_big->Add(ROC_tau32_manual_2groups_improved_total);
    leg_ROC_total_big->AddEntry(ROC_tau32_manual_2groups_improved_total, "2X3-jettiness Imp Jets, #beta = " + (TString)ss.str(), "L");

    ROC_tau32_manual_2groups_improved_total_allbeta.Add(ROC_tau32_manual_2groups_improved_total);


    TH1F *tau32_ttbar_manual_2groups_bigcone_total = (TH1F*)tau32_ttbar_manual_2groups_bigcone_total_hists.At(B);
    TH1F *tau32_dijets_manual_2groups_bigcone_total = (TH1F*)tau32_dijets_manual_2groups_bigcone_total_hists.At(B);

    tau32_ttbar_manual_2groups_bigcone_total->SetStats(0);
    tau32_dijets_manual_2groups_bigcone_total->SetStats(0);
    tau32_ttbar_manual_2groups_bigcone_total->Write();
    tau32_dijets_manual_2groups_bigcone_total->Write();

    int n_size_manual_2groups_bigcone_total = tau32_ttbar_manual_2groups_bigcone_total->GetSize() - 2;
    double integral_dijets_manual_2groups_bigcone_total[n_size_manual_2groups_bigcone_total], integral_ttbar_manual_2groups_bigcone_total[n_size_manual_2groups_bigcone_total];
    for (int i = 0; i < n_size_manual_2groups_bigcone_total; i++) {
      integral_ttbar_manual_2groups_bigcone_total[i] = tau32_ttbar_manual_2groups_bigcone_total->Integral(0,i)/total_total_ttbar_jets;
      integral_dijets_manual_2groups_bigcone_total[i] = tau32_dijets_manual_2groups_bigcone_total->Integral(0,i)/total_total_dijets_jets;
    }
    TGraph* ROC_tau32_manual_2groups_bigcone_total = new TGraph(n_size_manual_2groups_bigcone_total, integral_ttbar_manual_2groups_bigcone_total, integral_dijets_manual_2groups_bigcone_total);
    ROC_tau32_manual_2groups_bigcone_total->GetXaxis()->SetLimits(0, 1);
    ROC_tau32_manual_2groups_bigcone_total->GetYaxis()->SetLimits(0, 1);
    if (B == 0) ROC_tau32_manual_2groups_bigcone_total->SetLineColor(kRed);
    else ROC_tau32_manual_2groups_bigcone_total->SetLineColor(6);
    ROC_tau32_manual_2groups_bigcone_total->SetLineWidth(2);
    ROC_tau32_manual_2groups_bigcone_total->SetMarkerStyle(5);
    ROC_tau32_manual_2groups_bigcone_total->SetMarkerSize(2);
      // ROC_multigraph_total->Add(ROC_tau32_manual_2groups_bigcone_total);
      // leg_ROC_total->AddEntry(ROC_tau32_manual_2groups_bigcone_total, "2X3-jettiness Jets", "L");
    ROC_multigraph_total_big->Add(ROC_tau32_manual_2groups_bigcone_total);
    leg_ROC_total_big->AddEntry(ROC_tau32_manual_2groups_bigcone_total, "2X3-jettiness Jets, #beta = " + (TString)ss.str(), "L");

    ROC_tau32_manual_2groups_bigcone_total_allbeta.Add(ROC_tau32_manual_2groups_bigcone_total);


      // ROC_multigraph_total->Draw("AL");
      // ROC_multigraph_total->GetXaxis()->SetTitle("Tagging Efficiency");
      // ROC_multigraph_total->GetYaxis()->SetTitle("Mistag Rate");
      // ROC_multigraph_total->GetYaxis()->SetTitleOffset(0.85);
      // leg_ROC_total->Draw();
      // ROC_compare_total->Write();
      // ROC_compare_total->Print("ttbarstudy_2x3_plots/ttbarstudy_ROC_compare_totalpt_beta" + (TString)ss.str() + ".eps", "eps");
  }

  ROC_multigraph_total_big->Draw("AL");
  ROC_multigraph_total_big->GetXaxis()->SetTitle("Tagging Efficiency");
  ROC_multigraph_total_big->GetYaxis()->SetTitle("Mistag Rate");
  ROC_multigraph_total_big->GetYaxis()->SetTitleOffset(0.85);
  leg_ROC_total_big->Draw();
  ROC_compare_total_big->Write();
  ROC_compare_total_big->Print("ttbarstudy_2x3_plots/ttbarstudy_ROC_compare_totalpt_big.eps", "eps");

  for (int B = 0; B < betalist.size(); B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TCanvas *ROC_compare_total = new TCanvas("ROC_compare_total", "ROC Curve Comparison", 600, 600);
    ROC_compare_total->cd();
    ROC_compare_total->SetLogy();
    TMultiGraph *ROC_multigraph_total = new TMultiGraph("ROC_multigraph_total", "ROC Comparison for #tau_{3}/#tau_{2} (240 < p_{T} < 800, 160 < m < 240, #beta = " + (TString)ss.str() + ")");
    TLegend *leg_ROC_total = new TLegend(0.14, 0.6, 0.4, 0.88);
    leg_ROC_total->SetLineColor(0);
    leg_ROC_total->SetFillColor(0);

    TH1F *tau32_ttbar_antikt_wta_total = (TH1F*)tau32_ttbar_antikt_wta_total_hists.At(B);
    TH1F *tau32_dijets_antikt_wta_total = (TH1F*)tau32_dijets_antikt_wta_total_hists.At(B);
    int n_antikt_wta_total = tau32_ttbar_antikt_wta_total->GetSize() - 2;
    double integral_dijets_antikt_wta_total[n_antikt_wta_total], integral_ttbar_antikt_wta_total[n_antikt_wta_total];
    for (int i = 0; i < n_antikt_wta_total; i++) {
      integral_ttbar_antikt_wta_total[i] = tau32_ttbar_antikt_wta_total->Integral(0,i)/(total_total_ttbar_jets*n_betas);
      integral_dijets_antikt_wta_total[i] = tau32_dijets_antikt_wta_total->Integral(0,i)/(total_total_dijets_jets*n_betas);
    }
    TGraph* ROC_tau32_antikt_wta_total = new TGraph(n_antikt_wta_total, integral_ttbar_antikt_wta_total, integral_dijets_antikt_wta_total);
    ROC_tau32_antikt_wta_total->GetXaxis()->SetLimits(0, 1);
    ROC_tau32_antikt_wta_total->GetYaxis()->SetLimits(0, 1);
    ROC_tau32_antikt_wta_total->SetLineColor(kBlack);
    ROC_tau32_antikt_wta_total->SetLineWidth(2);
    ROC_tau32_antikt_wta_total->SetMarkerStyle(5);
    ROC_tau32_antikt_wta_total->SetMarkerSize(2);
    ROC_tau32_antikt_wta_total->Write();

    ROC_multigraph_total->Add(ROC_tau32_antikt_wta_total);
    leg_ROC_total->AddEntry(ROC_tau32_antikt_wta_total, "ak_{T}", "L");

    TGraph* ROC_tau32_antikt_total = (TGraph*)ROC_tau32_antikt_total_allbeta.At(B);
    TGraph* ROC_tau32_manual_6jet_total = (TGraph*)ROC_tau32_manual_6jet_total_allbeta.At(B);
    TGraph* ROC_tau32_manual_total = (TGraph*)ROC_tau32_manual_total_allbeta.At(B);
    TGraph* ROC_tau32_manual_2groups_total = (TGraph*)ROC_tau32_manual_2groups_total_allbeta.At(B);
    TGraph* ROC_tau32_manual_2groups_improved_total = (TGraph*)ROC_tau32_manual_2groups_improved_total_allbeta.At(B);
    TGraph* ROC_tau32_manual_2groups_bigcone_total = (TGraph*)ROC_tau32_manual_2groups_bigcone_total_allbeta.At(B);

      // ROC_multigraph_total->Add(ROC_tau32_antikt_total);
      // leg_ROC_total->AddEntry(ROC_tau32_antikt_total, "Genk_{T}", "L");
      // ROC_multigraph_total->Add(ROC_tau32_manual_6jet_total);
      // leg_ROC_total->AddEntry(ROC_tau32_manual_6jet_total, "6-jet", "L");

    ROC_multigraph_total->Add(ROC_tau32_manual_total);
    leg_ROC_total->AddEntry(ROC_tau32_manual_total, "6-jet imp", "L");

    ROC_multigraph_total->Add(ROC_tau32_manual_2groups_total);
    leg_ROC_total->AddEntry(ROC_tau32_manual_2groups_total, "2X3-jet", "L");

    ROC_multigraph_total->Add(ROC_tau32_manual_2groups_improved_total);
    leg_ROC_total->AddEntry(ROC_tau32_manual_2groups_improved_total, "2X3-jet imp", "L");

    ROC_multigraph_total->Add(ROC_tau32_manual_2groups_bigcone_total);
    leg_ROC_total->AddEntry(ROC_tau32_manual_2groups_bigcone_total, "2X3-jet bigcone", "L");

    ROC_tau32_antikt_total->SetLineColor(8);
    ROC_tau32_manual_6jet_total->SetLineColor(kOrange);
    ROC_tau32_manual_total->SetLineColor(kBlue);
    ROC_tau32_manual_2groups_total->SetLineColor(kRed);
    ROC_tau32_manual_2groups_improved_total->SetLineColor(kGreen);
    ROC_tau32_manual_2groups_bigcone_total->SetLineColor(kYellow);

    ROC_multigraph_total->Draw("AL");
    ROC_multigraph_total->GetXaxis()->SetTitle("Tagging Efficiency");
    ROC_multigraph_total->GetYaxis()->SetTitle("Mistag Rate");
    ROC_multigraph_total->GetYaxis()->SetTitleOffset(0.85);
    leg_ROC_total->Draw();
    ROC_compare_total->Write();
    ROC_compare_total->Print("ttbarstudy_2x3_plots/ttbarstudy_ROC_compare_totalpt_beta" + (TString)ss.str() + ".eps", "eps");

  }

  TMultiGraph *top_efficiency_rawmass_biggraph = new TMultiGraph();
  TMultiGraph *top_efficiency_mass_biggraph = new TMultiGraph();
  TMultiGraph *top_efficiency_mass_2groups_biggraph = new TMultiGraph();
  TMultiGraph *top_efficiency_mass_2groups_improved_biggraph = new TMultiGraph();
  TMultiGraph *top_efficiency_mass_2groups_Wmasscut_biggraph = new TMultiGraph();
  TMultiGraph *top_efficiency_mass_2groups_improved_Wmasscut_biggraph = new TMultiGraph();
  TMultiGraph *top_efficiency_mass_3jet4jet_2groups_biggraph = new TMultiGraph();
  TMultiGraph *top_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph = new TMultiGraph();
  TMultiGraph *top_efficiency_ratio_biggraph = new TMultiGraph();

  TLegend *leg_top_rawmass = new TLegend(0.14, 0.6, 0.42, 0.88);
  TLegend *leg_top_mass = new TLegend(0.14, 0.6, 0.42, 0.88);
  TLegend *leg_top_mass_2groups = new TLegend(0.14, 0.6, 0.42, 0.88);
  TLegend *leg_top_mass_2groups_improved = new TLegend(0.14, 0.6, 0.42, 0.88);
  TLegend *leg_top_mass_2groups_Wmasscut = new TLegend(0.14, 0.6, 0.42, 0.88);
  TLegend *leg_top_mass_2groups_improved_Wmasscut = new TLegend(0.14, 0.6, 0.42, 0.88);
  TLegend *leg_top_mass_3jet4jet_2groups = new TLegend(0.14, 0.6, 0.42, 0.88);
  TLegend *leg_top_mass_3jet4jet_2groups_Wmasscut = new TLegend(0.14, 0.6, 0.42, 0.88);
  TLegend *leg_top_ratio = new TLegend(0.14, 0.6, 0.42, 0.88);

  leg_top_rawmass->SetFillColor(kWhite);
  leg_top_mass->SetFillColor(kWhite);
  leg_top_mass_2groups->SetFillColor(kWhite);
  leg_top_mass_2groups_improved->SetFillColor(kWhite);
  leg_top_mass_2groups_Wmasscut->SetFillColor(kWhite);
  leg_top_mass_2groups_improved_Wmasscut->SetFillColor(kWhite);
  leg_top_mass_3jet4jet_2groups->SetFillColor(kWhite);
  leg_top_mass_3jet4jet_2groups_Wmasscut->SetFillColor(kWhite);
  leg_top_ratio->SetFillColor(kWhite);
  leg_top_rawmass->SetLineColor(kWhite);
  leg_top_mass->SetLineColor(kWhite);
  leg_top_mass_2groups->SetLineColor(kWhite);
  leg_top_mass_2groups_improved->SetLineColor(kWhite);
  leg_top_mass_2groups_Wmasscut->SetLineColor(kWhite);
  leg_top_mass_2groups_improved_Wmasscut->SetLineColor(kWhite);
  leg_top_mass_3jet4jet_2groups->SetLineColor(kWhite);
  leg_top_mass_3jet4jet_2groups_Wmasscut->SetLineColor(kWhite);
  leg_top_ratio->SetLineColor(kWhite);

  for (int i = 0; i < n_betas + 2; i++) {

    double beta;
    if (i != 0 && i != 1) beta = betalist[i - 2];

    ostringstream ss;
    ss << beta;

    TGraph *top_efficiency_rawmass_graph = new TGraph(n_perps, ptlist, top_efficiency_rawmass[i]);
    TGraph *top_efficiency_mass_graph = new TGraph(n_perps, ptlist, top_efficiency_mass[i]);
    TGraph *top_efficiency_mass_2groups_graph = new TGraph(n_perps, ptlist, top_efficiency_mass_2groups[i]);
    TGraph *top_efficiency_mass_2groups_improved_graph = new TGraph(n_perps, ptlist, top_efficiency_mass_2groups_improved[i]);
    TGraph *top_efficiency_mass_2groups_Wmasscut_graph = new TGraph(n_perps, ptlist, top_efficiency_mass_2groups_Wmasscut[i]);
    TGraph *top_efficiency_mass_2groups_improved_Wmasscut_graph = new TGraph(n_perps, ptlist, top_efficiency_mass_2groups_improved_Wmasscut[i]);
    TGraph *top_efficiency_mass_3jet4jet_2groups_graph = new TGraph(n_perps, ptlist, top_efficiency_mass_3jet4jet_2groups[i]);
    TGraph *top_efficiency_mass_3jet4jet_2groups_Wmasscut_graph = new TGraph(n_perps, ptlist, top_efficiency_mass_3jet4jet_2groups_Wmasscut[i]);

    top_efficiency_rawmass_graph->SetLineWidth(3);
    top_efficiency_mass_graph->SetLineWidth(3);
    top_efficiency_mass_2groups_graph->SetLineWidth(3);
    top_efficiency_mass_2groups_improved_graph->SetLineWidth(3);
    top_efficiency_mass_2groups_Wmasscut_graph->SetLineWidth(3);
    top_efficiency_mass_2groups_improved_Wmasscut_graph->SetLineWidth(3);
    top_efficiency_mass_3jet4jet_2groups_graph->SetLineWidth(3);
    top_efficiency_mass_3jet4jet_2groups_Wmasscut_graph->SetLineWidth(3);

  // higgs_efficiencies_1jet_graph->SetMarkerStyle(3);
  // higgs_efficiencies_1jet_graph->SetMarkerSize(2);
  // top_efficiency_rawmass_graph->SetMarkerStyle(3);
  // top_efficiency_rawmass_graph->SetMarkerSize(2);
  // top_efficiency_mass_graph->SetMarkerStyle(3);
  // top_efficiency_mass_graph->SetMarkerSize(2);
  // top_efficiency_ratio_graph->SetMarkerStyle(3);
  // top_efficiency_ratio_graph->SetMarkerSize(2);
  // angulardiff_1jet_graph->SetMarkerStyle(3);
  // angulardiff_1jet_graph->SetMarkerSize(2);
  // angulardiff_prop_1jet_graph->SetMarkerStyle(3);
  // angulardiff_prop_1jet_graph->SetMarkerSize(2);
  // angulardiff_2jets_graph->SetMarkerStyle(3);
  // angulardiff_2jets_graph->SetMarkerSize(2);
  // jet_distance_diff_1jet_graph->SetMarkerStyle(3);
  // jet_distance_diff_1jet_graph->SetMarkerSize(2);

    if (i == 0) {
      leg_top_rawmass->AddEntry(top_efficiency_rawmass_graph, "ak_{T} (2-jet)", "l");
      leg_top_mass->AddEntry(top_efficiency_mass_graph, "ak_{T} (2-jet)", "l");
      leg_top_mass_2groups->AddEntry(top_efficiency_mass_2groups_graph, "ak_{T} (2-jet)", "l");
      leg_top_mass_2groups_improved->AddEntry(top_efficiency_mass_2groups_improved_graph, "ak_{T} (2-jet)", "l");
      leg_top_mass_2groups_Wmasscut->AddEntry(top_efficiency_mass_2groups_Wmasscut_graph, "ak_{T} (2-jet)", "l");
      leg_top_mass_2groups_improved_Wmasscut->AddEntry(top_efficiency_mass_2groups_improved_Wmasscut_graph, "ak_{T} (2-jet)", "l");
      leg_top_mass_3jet4jet_2groups->AddEntry(top_efficiency_mass_3jet4jet_2groups_graph, "ak_{T} (2-jet)", "l");
      leg_top_mass_3jet4jet_2groups_Wmasscut->AddEntry(top_efficiency_mass_3jet4jet_2groups_Wmasscut_graph, "ak_{T} (2-jet)", "l");
    }

    if (i == 1) {

      top_efficiency_rawmass_graph->SetLineColor(8);
      top_efficiency_mass_graph->SetLineColor(8);
      top_efficiency_mass_2groups_graph->SetLineColor(8);
      top_efficiency_mass_2groups_improved_graph->SetLineColor(8);
      top_efficiency_mass_2groups_Wmasscut_graph->SetLineColor(8);
      top_efficiency_mass_2groups_improved_Wmasscut_graph->SetLineColor(8);
      top_efficiency_mass_3jet4jet_2groups_graph->SetLineColor(8);
      top_efficiency_mass_3jet4jet_2groups_Wmasscut_graph->SetLineColor(8);

      leg_top_rawmass->AddEntry(top_efficiency_rawmass_graph, "ak_{T} (6-jet)", "l");
      leg_top_mass->AddEntry(top_efficiency_mass_graph, "ak_{T} (6-jet)", "l");
      leg_top_mass_2groups->AddEntry(top_efficiency_mass_2groups_graph, "ak_{T} (2 groups)", "l");
      leg_top_mass_2groups_improved->AddEntry(top_efficiency_mass_2groups_improved_graph, "ak_{T} (6-jet)", "l");
      leg_top_mass_2groups_Wmasscut->AddEntry(top_efficiency_mass_2groups_Wmasscut_graph, "ak_{T} (6-jet)", "l");
      leg_top_mass_2groups_Wmasscut->AddEntry(top_efficiency_mass_2groups_improved_Wmasscut_graph, "ak_{T} (6-jet)", "l");
      leg_top_mass_3jet4jet_2groups->AddEntry(top_efficiency_mass_3jet4jet_2groups_graph, "ak_{T} (6-jet)", "l");
      leg_top_mass_3jet4jet_2groups_Wmasscut->AddEntry(top_efficiency_mass_3jet4jet_2groups_Wmasscut_graph, "ak_{T} (6-jet)", "l");
    }

    if (i != 0 && i != 1) {
      TGraph *top_efficiency_ratio_graph = new TGraph(n_perps, ptlist, top_efficiency_ratio[i-2]);

      top_efficiency_rawmass_graph->SetMarkerColor(colorlist[i-2]);
      top_efficiency_mass_graph->SetMarkerColor(colorlist[i-2]);
      top_efficiency_mass_2groups_graph->SetMarkerColor(colorlist[i-2]);
      top_efficiency_mass_2groups_improved_graph->SetMarkerColor(colorlist[i-2]);
      top_efficiency_ratio_graph->SetMarkerColor(colorlist[i-2]);
      top_efficiency_rawmass_graph->SetLineColor(colorlist[i-2]);
      top_efficiency_mass_graph->SetLineColor(colorlist[i-2]);
      top_efficiency_mass_2groups_graph->SetLineColor(colorlist[i-2]);
      top_efficiency_mass_2groups_improved_graph->SetLineColor(colorlist[i-2]);
      top_efficiency_ratio_graph->SetLineColor(colorlist[i-2]);
      top_efficiency_mass_2groups_Wmasscut_graph->SetLineColor(colorlist[i-2]);
      top_efficiency_mass_2groups_improved_Wmasscut_graph->SetLineColor(colorlist[i-2]);
      top_efficiency_mass_3jet4jet_2groups_graph->SetLineColor(colorlist[i-2]);
      top_efficiency_mass_3jet4jet_2groups_Wmasscut_graph->SetLineColor(colorlist[i-2]);

      leg_top_rawmass->AddEntry(top_efficiency_rawmass_graph, "#beta = " + (TString)ss.str(), "l");
      leg_top_mass->AddEntry(top_efficiency_mass_graph, "#beta = " + (TString)ss.str(), "l");
      leg_top_mass_2groups->AddEntry(top_efficiency_mass_2groups_graph, "#beta = " + (TString)ss.str(), "l");
      leg_top_mass_2groups_improved->AddEntry(top_efficiency_mass_2groups_improved_graph, "#beta = " + (TString)ss.str(), "l");
      leg_top_ratio->AddEntry(top_efficiency_ratio_graph, "#beta = " + (TString)ss.str(), "l");
      leg_top_mass_2groups_Wmasscut->AddEntry(top_efficiency_mass_2groups_Wmasscut_graph, "#beta = " + (TString)ss.str(), "l");
      leg_top_mass_2groups_improved_Wmasscut->AddEntry(top_efficiency_mass_2groups_improved_Wmasscut_graph, "#beta = " + (TString)ss.str(), "l");
      leg_top_mass_3jet4jet_2groups->AddEntry(top_efficiency_mass_3jet4jet_2groups_graph, "#beta = " + (TString)ss.str(), "l");
      leg_top_mass_3jet4jet_2groups_Wmasscut->AddEntry(top_efficiency_mass_3jet4jet_2groups_Wmasscut_graph, "#beta = " + (TString)ss.str(), "l");

      top_efficiency_ratio_biggraph->Add(top_efficiency_ratio_graph);
    }

    top_efficiency_rawmass_biggraph->Add(top_efficiency_rawmass_graph);
    top_efficiency_mass_biggraph->Add(top_efficiency_mass_graph);
    top_efficiency_mass_2groups_biggraph->Add(top_efficiency_mass_2groups_graph);
    top_efficiency_mass_2groups_improved_biggraph->Add(top_efficiency_mass_2groups_improved_graph);
    top_efficiency_mass_2groups_Wmasscut_biggraph->Add(top_efficiency_mass_2groups_Wmasscut_graph);
    top_efficiency_mass_2groups_improved_Wmasscut_biggraph->Add(top_efficiency_mass_2groups_improved_Wmasscut_graph);
    top_efficiency_mass_3jet4jet_2groups_biggraph->Add(top_efficiency_mass_3jet4jet_2groups_graph);
    top_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->Add(top_efficiency_mass_3jet4jet_2groups_Wmasscut_graph);
  }

  TCanvas *top_rawmass_can = new TCanvas("top_rawmass_can", "top_rawmass_can", 600, 600);
  top_rawmass_can->cd();
  top_efficiency_rawmass_biggraph->Draw("ALP");
  leg_top_rawmass->Draw();
  top_efficiency_rawmass_biggraph->SetMinimum(0.);     
  top_efficiency_rawmass_biggraph->SetMaximum(1.2);
  top_efficiency_rawmass_biggraph->SetTitle("Top Efficiencies for Raw 6-jets (160 < m < 240)");
  top_efficiency_rawmass_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  top_efficiency_rawmass_biggraph->GetYaxis()->SetTitle("% mass");
// top_efficiency_rawmass_biggraph->GetYaxis()->SetTitleOffset(1.3);
  top_rawmass_can->Write();
  top_rawmass_can->Print("ttbarstudy_2x3_plots/top_efficiency_rawmass.eps", "eps");

  TCanvas *top_mass_can = new TCanvas("top_mass_can", "top_mass_can", 600, 600);
  top_mass_can->cd();
  top_efficiency_mass_biggraph->Draw("ALP");
  leg_top_mass->Draw();
  top_efficiency_mass_biggraph->SetMinimum(0.);     
  top_efficiency_mass_biggraph->SetMaximum(1.2);
  top_efficiency_mass_biggraph->SetTitle("Top Efficiencies for Improved 6-jets (160 < m < 240)");
  top_efficiency_mass_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  top_efficiency_mass_biggraph->GetYaxis()->SetTitle("% mass");
// top_efficiency_mass_biggraph->GetYaxis()->SetTitleOffset(1.3);
  top_mass_can->Write();
  top_mass_can->Print("ttbarstudy_2x3_plots/top_efficiency_mass.eps", "eps");

  TCanvas *top_mass_2groups_can = new TCanvas("top_mass_2groups_can", "top_mass_2groups_can", 600, 600);
  top_mass_2groups_can->cd();
  top_efficiency_mass_2groups_biggraph->Draw("ALP");
  leg_top_mass_2groups->Draw();
  top_efficiency_mass_2groups_biggraph->SetMinimum(0.);     
  top_efficiency_mass_2groups_biggraph->SetMaximum(1.2);
  top_efficiency_mass_2groups_biggraph->SetTitle("Top Efficiencies for 2#times3 Method (160 < m < 240)");
  top_efficiency_mass_2groups_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  top_efficiency_mass_2groups_biggraph->GetYaxis()->SetTitle("% mass");
// top_efficiency_mass_2groups_biggraph->GetYaxis()->SetTitleOffset(1.3);
  top_mass_2groups_can->Write();
  top_mass_2groups_can->Print("ttbarstudy_2x3_plots/top_efficiency_mass_2groups.eps", "eps");

  TCanvas *top_mass_2groups_improved_can = new TCanvas("top_mass_2groups_improved_can", "top_mass_2groups_improved_can", 600, 600);
  top_mass_2groups_improved_can->cd();
  top_efficiency_mass_2groups_improved_biggraph->Draw("ALP");
  leg_top_mass_2groups_improved->Draw();
  top_efficiency_mass_2groups_improved_biggraph->SetMinimum(0.);     
  top_efficiency_mass_2groups_improved_biggraph->SetMaximum(1.2);
  top_efficiency_mass_2groups_improved_biggraph->SetTitle("Top Efficiencies for Improved 2#times3 Method (160 < m < 240)");
  top_efficiency_mass_2groups_improved_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  top_efficiency_mass_2groups_improved_biggraph->GetYaxis()->SetTitle("% mass");
// top_efficiency_mass_2groups_improved_biggraph->GetYaxis()->SetTitleOffset(1.3);
  top_mass_2groups_improved_can->Write();
  top_mass_2groups_improved_can->Print("ttbarstudy_2x3_plots/top_efficiency_mass_2groups_improved.eps", "eps");

  TCanvas *top_mass_2groups_Wmasscut_can = new TCanvas("top_mass_2groups_Wmasscut_can", "top_mass_2groups_Wmasscut_can", 600, 600);
  top_mass_2groups_Wmasscut_can->cd();
  top_efficiency_mass_2groups_Wmasscut_biggraph->Draw("ALP");
  leg_top_mass_2groups_Wmasscut->Draw();
  top_efficiency_mass_2groups_Wmasscut_biggraph->SetMinimum(0.);     
  top_efficiency_mass_2groups_Wmasscut_biggraph->SetMaximum(1.2);
  top_efficiency_mass_2groups_Wmasscut_biggraph->SetTitle("Top Efficiencies for Wmasscut 2#times3 Method (160 < m < 240)");
  top_efficiency_mass_2groups_Wmasscut_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  top_efficiency_mass_2groups_Wmasscut_biggraph->GetYaxis()->SetTitle("% mass");
// top_efficiency_mass_2groups_Wmasscut_biggraph->GetYaxis()->SetTitleOffset(1.3);
  top_mass_2groups_Wmasscut_can->Write();
  top_mass_2groups_Wmasscut_can->Print("ttbarstudy_2x3_plots/top_efficiency_mass_2groups_Wmasscut.eps", "eps");

  TCanvas *top_mass_2groups_improved_Wmasscut_can = new TCanvas("top_mass_2groups_improved_Wmasscut_can", "top_mass_2groups_improved_Wmasscut_can", 600, 600);
  top_mass_2groups_improved_Wmasscut_can->cd();
  top_efficiency_mass_2groups_improved_Wmasscut_biggraph->Draw("ALP");
  leg_top_mass_2groups_improved_Wmasscut->Draw();
  top_efficiency_mass_2groups_improved_Wmasscut_biggraph->SetMinimum(0.);     
  top_efficiency_mass_2groups_improved_Wmasscut_biggraph->SetMaximum(1.2);
  top_efficiency_mass_2groups_improved_Wmasscut_biggraph->SetTitle("Top Efficiencies for Wmasscut 2#times3 Method (160 < m < 240)");
  top_efficiency_mass_2groups_improved_Wmasscut_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  top_efficiency_mass_2groups_improved_Wmasscut_biggraph->GetYaxis()->SetTitle("% mass");
// top_efficiency_mass_2groups_improved_Wmasscut_biggraph->GetYaxis()->SetTitleOffset(1.3);
  top_mass_2groups_improved_Wmasscut_can->Write();
  top_mass_2groups_improved_Wmasscut_can->Print("ttbarstudy_2x3_plots/top_efficiency_mass_2groups_improved_Wmasscut.eps", "eps");


  TCanvas *top_mass_3jet4jet_2groups_can = new TCanvas("top_mass_3jet4jet_2groups_can", "top_mass_3jet4jet_2groups_can", 600, 600);
  top_mass_3jet4jet_2groups_can->cd();
  top_efficiency_mass_3jet4jet_2groups_biggraph->Draw("ALP");
  leg_top_mass_3jet4jet_2groups->Draw();
  top_efficiency_mass_3jet4jet_2groups_biggraph->SetMinimum(0.);     
  top_efficiency_mass_3jet4jet_2groups_biggraph->SetMaximum(1.2);
  top_efficiency_mass_3jet4jet_2groups_biggraph->SetTitle("Top Efficiencies for 2#times3 OR 4 Method (160 < m < 240)");
  top_efficiency_mass_3jet4jet_2groups_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  top_efficiency_mass_3jet4jet_2groups_biggraph->GetYaxis()->SetTitle("% mass");
// top_efficiency_mass_3jet4jet_2groups_biggraph->GetYaxis()->SetTitleOffset(1.3);
  top_mass_3jet4jet_2groups_can->Write();
  top_mass_3jet4jet_2groups_can->Print("ttbarstudy_2x3_plots/top_efficiency_mass_3jet4jet_2groups.eps", "eps");


  TCanvas *top_mass_3jet4jet_2groups_Wmasscut_can = new TCanvas("top_mass_3jet4jet_2groups_Wmasscut_can", "top_mass_3jet4jet_2groups_Wmasscut_can", 600, 600);
  top_mass_3jet4jet_2groups_Wmasscut_can->cd();
  top_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->Draw("ALP");
  leg_top_mass_3jet4jet_2groups_Wmasscut->Draw();
  top_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->SetMinimum(0.);     
  top_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->SetMaximum(1.2);
  top_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->SetTitle("Top Efficiencies for Wmasscut 2#times3 OR 4 Method (160 < m < 240)");
  top_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  top_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->GetYaxis()->SetTitle("% mass");
// top_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->GetYaxis()->SetTitleOffset(1.3);
  top_mass_3jet4jet_2groups_Wmasscut_can->Write();
  top_mass_3jet4jet_2groups_Wmasscut_can->Print("ttbarstudy_2x3_plots/top_efficiency_mass_3jet4jet_2groups_Wmasscut.eps", "eps");


  TCanvas *top_ratio_can = new TCanvas("top_ratio_can", "top_ratio_can", 800, 600);
  top_ratio_can->cd();
  top_efficiency_ratio_biggraph->Draw("ALP");
  leg_top_ratio->Draw();
  top_efficiency_ratio_biggraph->SetMinimum(0.0);     
  top_efficiency_ratio_biggraph->SetMaximum(2.0);
  top_efficiency_ratio_biggraph->SetTitle("Ratio of 3-jettiness/2-jettiness Top Efficiencies");
  top_efficiency_ratio_biggraph->GetXaxis()->SetTitle("pT_{min} (GeV)");
  top_efficiency_ratio_biggraph->GetYaxis()->SetTitle("3-jet eff/2-jet eff");
  top_efficiency_ratio_biggraph->GetYaxis()->SetTitleOffset(1.3);
  top_ratio_can->Write();
  top_ratio_can->Print("ttbarstudy_2x3_plots/top_efficiency_ratio.eps", "eps");

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


  TMultiGraph *qcd_efficiency_rawmass_biggraph = new TMultiGraph();
  TMultiGraph *qcd_efficiency_mass_biggraph = new TMultiGraph();
  TMultiGraph *qcd_efficiency_mass_2groups_biggraph = new TMultiGraph();
  TMultiGraph *qcd_efficiency_mass_2groups_improved_biggraph = new TMultiGraph();
  TMultiGraph *qcd_efficiency_mass_2groups_Wmasscut_biggraph = new TMultiGraph();
  TMultiGraph *qcd_efficiency_mass_2groups_improved_Wmasscut_biggraph = new TMultiGraph();
  TMultiGraph *qcd_efficiency_mass_3jet4jet_2groups_biggraph = new TMultiGraph();
  TMultiGraph *qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph = new TMultiGraph();

  TMultiGraph *qcd_efficiency_ratio_biggraph = new TMultiGraph();

  TLegend *leg_qcd_rawmass = new TLegend(0.14, 0.6, 0.4, 0.88);
  TLegend *leg_qcd_mass = new TLegend(0.14, 0.6, 0.4, 0.88);
  TLegend *leg_qcd_mass_2groups = new TLegend(0.14, 0.6, 0.4, 0.88);
  TLegend *leg_qcd_mass_2groups_improved = new TLegend(0.14, 0.6, 0.4, 0.88);
  TLegend *leg_qcd_mass_2groups_Wmasscut = new TLegend(0.14, 0.6, 0.4, 0.88);
  TLegend *leg_qcd_mass_2groups_improved_Wmasscut = new TLegend(0.14, 0.6, 0.4, 0.88);
  TLegend *leg_qcd_mass_3jet4jet_2groups = new TLegend(0.14, 0.6, 0.4, 0.88);
  TLegend *leg_qcd_mass_3jet4jet_2groups_Wmasscut = new TLegend(0.14, 0.6, 0.4, 0.88);
  TLegend *leg_qcd_ratio = new TLegend(0.14, 0.6, 0.4, 0.88);

  leg_qcd_rawmass->SetFillColor(kWhite);
  leg_qcd_mass->SetFillColor(kWhite);
  leg_qcd_mass_2groups->SetFillColor(kWhite);
  leg_qcd_mass_2groups_improved->SetFillColor(kWhite);
  leg_qcd_mass_2groups_Wmasscut->SetFillColor(kWhite);
  leg_qcd_mass_2groups_improved_Wmasscut->SetFillColor(kWhite);
  leg_qcd_mass_3jet4jet_2groups->SetFillColor(kWhite);
  leg_qcd_mass_3jet4jet_2groups_Wmasscut->SetFillColor(kWhite);
  leg_qcd_ratio->SetFillColor(kWhite);
  leg_qcd_rawmass->SetLineColor(kWhite);
  leg_qcd_mass->SetLineColor(kWhite);
  leg_qcd_mass_2groups->SetLineColor(kWhite);
  leg_qcd_mass_2groups_improved->SetLineColor(kWhite);
  leg_qcd_mass_2groups_Wmasscut->SetLineColor(kWhite);
  leg_qcd_mass_2groups_improved_Wmasscut->SetLineColor(kWhite);
  leg_qcd_mass_3jet4jet_2groups->SetLineColor(kWhite);
  leg_qcd_mass_3jet4jet_2groups_Wmasscut->SetLineColor(kWhite);
  leg_qcd_ratio->SetLineColor(kWhite);

  for (int i = 0; i < n_betas + 2; i++) {

    double beta;
    if (i != 0 && i != 1) beta = betalist[i - 2];

    ostringstream ss;
    ss << beta;

    TGraph *qcd_efficiency_rawmass_graph = new TGraph(n_perps, ptlist, qcd_efficiency_rawmass[i]);
    TGraph *qcd_efficiency_mass_graph = new TGraph(n_perps, ptlist, qcd_efficiency_mass[i]);
    TGraph *qcd_efficiency_mass_2groups_graph = new TGraph(n_perps, ptlist, qcd_efficiency_mass_2groups[i]);
    TGraph *qcd_efficiency_mass_2groups_improved_graph = new TGraph(n_perps, ptlist, qcd_efficiency_mass_2groups_improved[i]);
    TGraph *qcd_efficiency_mass_2groups_Wmasscut_graph = new TGraph(n_perps, ptlist, qcd_efficiency_mass_2groups_Wmasscut[i]);
    TGraph *qcd_efficiency_mass_2groups_improved_Wmasscut_graph = new TGraph(n_perps, ptlist, qcd_efficiency_mass_2groups_improved_Wmasscut[i]);
    TGraph *qcd_efficiency_mass_3jet4jet_2groups_graph = new TGraph(n_perps, ptlist, qcd_efficiency_mass_3jet4jet_2groups[i]);
    TGraph *qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_graph = new TGraph(n_perps, ptlist, qcd_efficiency_mass_3jet4jet_2groups_Wmasscut[i]);

    qcd_efficiency_rawmass_graph->SetLineWidth(3);
    qcd_efficiency_mass_graph->SetLineWidth(3);
    qcd_efficiency_mass_2groups_graph->SetLineWidth(3);
    qcd_efficiency_mass_2groups_improved_graph->SetLineWidth(3);
    qcd_efficiency_mass_2groups_Wmasscut_graph->SetLineWidth(3);
    qcd_efficiency_mass_2groups_improved_Wmasscut_graph->SetLineWidth(3);
    qcd_efficiency_mass_3jet4jet_2groups_graph->SetLineWidth(3);
    qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_graph->SetLineWidth(3);
  // higgs_efficiencies_1jet_graph->SetMarkerStyle(3);
  // higgs_efficiencies_1jet_graph->SetMarkerSize(2);
  // qcd_efficiency_rawmass_graph->SetMarkerStyle(3);
  // qcd_efficiency_rawmass_graph->SetMarkerSize(2);
  // qcd_efficiency_mass_graph->SetMarkerStyle(3);
  // qcd_efficiency_mass_graph->SetMarkerSize(2);
  // qcd_efficiency_ratio_graph->SetMarkerStyle(3);
  // qcd_efficiency_ratio_graph->SetMarkerSize(2);
  // angulardiff_1jet_graph->SetMarkerStyle(3);
  // angulardiff_1jet_graph->SetMarkerSize(2);
  // angulardiff_prop_1jet_graph->SetMarkerStyle(3);
  // angulardiff_prop_1jet_graph->SetMarkerSize(2);
  // angulardiff_2jets_graph->SetMarkerStyle(3);
  // angulardiff_2jets_graph->SetMarkerSize(2);
  // jet_distance_diff_1jet_graph->SetMarkerStyle(3);
  // jet_distance_diff_1jet_graph->SetMarkerSize(2);

    if (i == 0) {
      leg_qcd_rawmass->AddEntry(qcd_efficiency_rawmass_graph, "ak_{T} (2-jet)", "l");
      leg_qcd_mass->AddEntry(qcd_efficiency_mass_graph, "ak_{T} (2-jet)", "l");
      leg_qcd_mass_2groups->AddEntry(qcd_efficiency_mass_2groups_graph, "ak_{T} (2-jet)", "l");
      leg_qcd_mass_2groups_improved->AddEntry(qcd_efficiency_mass_2groups_improved_graph, "ak_{T} (2-jet)", "l");
      leg_qcd_mass_2groups_Wmasscut->AddEntry(qcd_efficiency_mass_2groups_Wmasscut_graph, "ak_{T} (2-jet)", "l");
      leg_qcd_mass_2groups_improved_Wmasscut->AddEntry(qcd_efficiency_mass_2groups_improved_Wmasscut_graph, "ak_{T} (2-jet)", "l");
      leg_qcd_mass_3jet4jet_2groups->AddEntry(qcd_efficiency_mass_3jet4jet_2groups_graph, "ak_{T} (2-jet)", "l");
      leg_qcd_mass_3jet4jet_2groups_Wmasscut->AddEntry(qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_graph, "ak_{T} (2-jet)", "l");
    }

    if (i == 1) {
      qcd_efficiency_rawmass_graph->SetLineColor(8);
      qcd_efficiency_mass_graph->SetLineColor(8);
      qcd_efficiency_mass_2groups_graph->SetLineColor(8);
      qcd_efficiency_mass_2groups_improved_graph->SetLineColor(8);
      qcd_efficiency_mass_2groups_Wmasscut_graph->SetLineColor(8);
      qcd_efficiency_mass_2groups_improved_Wmasscut_graph->SetLineColor(8);
      qcd_efficiency_mass_3jet4jet_2groups_graph->SetLineColor(8);
      qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_graph->SetLineColor(8);

      leg_qcd_rawmass->AddEntry(qcd_efficiency_rawmass_graph, "ak_{T} (6-jet)", "l");
      leg_qcd_mass->AddEntry(qcd_efficiency_mass_graph, "ak_{T} (6-jet)", "l");
      leg_qcd_mass_2groups->AddEntry(qcd_efficiency_mass_2groups_graph, "ak_{T} (6-jet)", "l");
      leg_qcd_mass_2groups_improved->AddEntry(qcd_efficiency_mass_2groups_improved_graph, "ak_{T} (6-jet)", "l");
      leg_qcd_mass_2groups_Wmasscut->AddEntry(qcd_efficiency_mass_2groups_Wmasscut_graph, "ak_{T} (6-jet)", "l");
      leg_qcd_mass_2groups_improved_Wmasscut->AddEntry(qcd_efficiency_mass_2groups_improved_Wmasscut_graph, "ak_{T} (6-jet)", "l");
      leg_qcd_mass_3jet4jet_2groups->AddEntry(qcd_efficiency_mass_3jet4jet_2groups_graph, "ak_{T} (6-jet)", "l");
      leg_qcd_mass_3jet4jet_2groups_Wmasscut->AddEntry(qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_graph, "ak_{T} (6-jet)", "l");
    }

    if (i != 0 && i != 1) {
      TGraph *qcd_efficiency_ratio_graph = new TGraph(n_perps, ptlist, qcd_efficiency_ratio[i-2]);

      qcd_efficiency_rawmass_graph->SetMarkerColor(colorlist[i-2]);
      qcd_efficiency_mass_graph->SetMarkerColor(colorlist[i-2]);
      qcd_efficiency_mass_2groups_graph->SetMarkerColor(colorlist[i-2]);
      qcd_efficiency_mass_2groups_improved_graph->SetMarkerColor(colorlist[i-2]);
      qcd_efficiency_ratio_graph->SetMarkerColor(colorlist[i-2]);
      qcd_efficiency_rawmass_graph->SetLineColor(colorlist[i-2]);
      qcd_efficiency_mass_graph->SetLineColor(colorlist[i-2]);
      qcd_efficiency_mass_2groups_graph->SetLineColor(colorlist[i-2]);
      qcd_efficiency_mass_2groups_improved_graph->SetLineColor(colorlist[i-2]);
      qcd_efficiency_mass_2groups_Wmasscut_graph->SetLineColor(colorlist[i-2]);
      qcd_efficiency_mass_2groups_improved_Wmasscut_graph->SetLineColor(colorlist[i-2]);
      qcd_efficiency_ratio_graph->SetLineColor(colorlist[i-2]);
      qcd_efficiency_mass_3jet4jet_2groups_graph->SetLineColor(colorlist[i-2]);
      qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_graph->SetLineColor(colorlist[i-2]);

      leg_qcd_rawmass->AddEntry(qcd_efficiency_rawmass_graph, "#beta = " + (TString)ss.str(), "l");
      leg_qcd_mass->AddEntry(qcd_efficiency_mass_graph, "#beta = " + (TString)ss.str(), "l");
      leg_qcd_mass_2groups->AddEntry(qcd_efficiency_mass_2groups_graph, "#beta = " + (TString)ss.str(), "l");
      leg_qcd_mass_2groups_improved->AddEntry(qcd_efficiency_mass_2groups_improved_graph, "#beta = " + (TString)ss.str(), "l");
      leg_qcd_mass_2groups_Wmasscut->AddEntry(qcd_efficiency_mass_2groups_Wmasscut_graph, "#beta = " + (TString)ss.str(), "l");
      leg_qcd_mass_2groups_improved_Wmasscut->AddEntry(qcd_efficiency_mass_2groups_improved_Wmasscut_graph, "#beta = " + (TString)ss.str(), "l");
      leg_qcd_mass_3jet4jet_2groups->AddEntry(qcd_efficiency_mass_3jet4jet_2groups_graph, "#beta = " + (TString)ss.str(), "l");
      leg_qcd_mass_3jet4jet_2groups_Wmasscut->AddEntry(qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_graph, "#beta = " + (TString)ss.str(), "l");
      leg_qcd_ratio->AddEntry(qcd_efficiency_ratio_graph, "#beta = " + (TString)ss.str(), "l");

      qcd_efficiency_ratio_biggraph->Add(qcd_efficiency_ratio_graph);
    }

    qcd_efficiency_rawmass_biggraph->Add(qcd_efficiency_rawmass_graph);
    qcd_efficiency_mass_biggraph->Add(qcd_efficiency_mass_graph);
    qcd_efficiency_mass_2groups_biggraph->Add(qcd_efficiency_mass_2groups_graph);
    qcd_efficiency_mass_2groups_improved_biggraph->Add(qcd_efficiency_mass_2groups_improved_graph);
    qcd_efficiency_mass_2groups_Wmasscut_biggraph->Add(qcd_efficiency_mass_2groups_Wmasscut_graph);
    qcd_efficiency_mass_2groups_improved_Wmasscut_biggraph->Add(qcd_efficiency_mass_2groups_improved_Wmasscut_graph);
    qcd_efficiency_mass_3jet4jet_2groups_biggraph->Add(qcd_efficiency_mass_3jet4jet_2groups_graph);
    qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->Add(qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_graph);
  }

  TCanvas *qcd_rawmass_can = new TCanvas("qcd_rawmass_can", "qcd_rawmass_can", 600, 600);
  qcd_rawmass_can->cd();
  qcd_efficiency_rawmass_biggraph->Draw("ALP");
  leg_qcd_rawmass->Draw();
  qcd_efficiency_rawmass_biggraph->SetMinimum(0.);     
  qcd_efficiency_rawmass_biggraph->SetMaximum(1.2);
  qcd_efficiency_rawmass_biggraph->SetTitle("QCD Efficiencies for Raw 6-jets (160 < m < 240)");
  qcd_efficiency_rawmass_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  qcd_efficiency_rawmass_biggraph->GetYaxis()->SetTitle("% mass");
// qcd_efficiency_rawmass_biggraph->GetYaxis()->SetTitleOffset(1.3);
  qcd_rawmass_can->Write();
  qcd_rawmass_can->Print("ttbarstudy_2x3_plots/qcd_efficiency_rawmass.eps", "eps");

  TCanvas *qcd_mass_can = new TCanvas("qcd_mass_can", "qcd_mass_can", 600, 600);
  qcd_mass_can->cd();
  qcd_efficiency_mass_biggraph->Draw("ALP");
  leg_qcd_mass->Draw();
  qcd_efficiency_mass_biggraph->SetMinimum(0.);     
  qcd_efficiency_mass_biggraph->SetMaximum(1.2);
  qcd_efficiency_mass_biggraph->SetTitle("QCD Efficiencies for Improved 6-jets");
  qcd_efficiency_mass_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  qcd_efficiency_mass_biggraph->GetYaxis()->SetTitle("% mass 160-240");
// qcd_efficiency_mass_biggraph->GetYaxis()->SetTitleOffset(1.3);
  qcd_mass_can->Write();
  qcd_mass_can->Print("ttbarstudy_2x3_plots/qcd_efficiency_mass.eps", "eps");

  TCanvas *qcd_mass_2groups_can = new TCanvas("qcd_mass_2groups_can", "qcd_mass_2groups_can", 600, 600);
  qcd_mass_2groups_can->cd();
  qcd_efficiency_mass_2groups_biggraph->Draw("ALP");
  leg_qcd_mass_2groups->Draw();
  qcd_efficiency_mass_2groups_biggraph->SetMinimum(0.);     
  qcd_efficiency_mass_2groups_biggraph->SetMaximum(1.2);
  qcd_efficiency_mass_2groups_biggraph->SetTitle("QCD Efficiencies for 2#times3 Method (160 < m < 240)");
  qcd_efficiency_mass_2groups_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  qcd_efficiency_mass_2groups_biggraph->GetYaxis()->SetTitle("% mass ");
// qcd_efficiency_mass_2groups_biggraph->GetYaxis()->SetTitleOffset(1.3);
  qcd_mass_2groups_can->Write();
  qcd_mass_2groups_can->Print("ttbarstudy_2x3_plots/qcd_efficiency_mass_2groups.eps", "eps");


  TCanvas *qcd_mass_2groups_improved_can = new TCanvas("qcd_mass_2groups_improved_can", "qcd_mass_2groups_improved_can", 600, 600);
  qcd_mass_2groups_improved_can->cd();
  qcd_efficiency_mass_2groups_improved_biggraph->Draw("ALP");
  leg_qcd_mass_2groups_improved->Draw();
  qcd_efficiency_mass_2groups_improved_biggraph->SetMinimum(0.);     
  qcd_efficiency_mass_2groups_improved_biggraph->SetMaximum(1.2);
  qcd_efficiency_mass_2groups_improved_biggraph->SetTitle("QCD Efficiencies for Improved 2#times3 Method (160 < m < 240)");
  qcd_efficiency_mass_2groups_improved_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  qcd_efficiency_mass_2groups_improved_biggraph->GetYaxis()->SetTitle("% mass");
// qcd_efficiency_mass_2groups_improved_biggraph->GetYaxis()->SetTitleOffset(1.3);
  qcd_mass_2groups_improved_can->Write();
  qcd_mass_2groups_improved_can->Print("ttbarstudy_2x3_plots/qcd_efficiency_mass_2groups_improved.eps", "eps");


  TCanvas *qcd_mass_2groups_Wmasscut_can = new TCanvas("qcd_mass_2groups_Wmasscut_can", "qcd_mass_2groups_Wmasscut_can", 600, 600);
  qcd_mass_2groups_Wmasscut_can->cd();
  qcd_efficiency_mass_2groups_Wmasscut_biggraph->Draw("ALP");
  leg_qcd_mass_2groups_Wmasscut->Draw();
  qcd_efficiency_mass_2groups_Wmasscut_biggraph->SetMinimum(0.);     
  qcd_efficiency_mass_2groups_Wmasscut_biggraph->SetMaximum(1.2);
  qcd_efficiency_mass_2groups_Wmasscut_biggraph->SetTitle("QCD Efficiencies for Wmasscut 2#times3 Method (160 < m < 240)");
  qcd_efficiency_mass_2groups_Wmasscut_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  qcd_efficiency_mass_2groups_Wmasscut_biggraph->GetYaxis()->SetTitle("% mass");
// qcd_efficiency_mass_2groups_Wmasscut_biggraph->GetYaxis()->SetTitleOffset(1.3);
  qcd_mass_2groups_Wmasscut_can->Write();
  qcd_mass_2groups_Wmasscut_can->Print("ttbarstudy_2x3_plots/qcd_efficiency_mass_2groups_Wmasscut.eps", "eps");

  TCanvas *qcd_mass_2groups_improved_Wmasscut_can = new TCanvas("qcd_mass_2groups_improved_Wmasscut_can", "qcd_mass_2groups_improved_Wmasscut_can", 600, 600);
  qcd_mass_2groups_improved_Wmasscut_can->cd();
  qcd_efficiency_mass_2groups_improved_Wmasscut_biggraph->Draw("ALP");
  leg_qcd_mass_2groups_improved_Wmasscut->Draw();
  qcd_efficiency_mass_2groups_improved_Wmasscut_biggraph->SetMinimum(0.);     
  qcd_efficiency_mass_2groups_improved_Wmasscut_biggraph->SetMaximum(1.2);
  qcd_efficiency_mass_2groups_improved_Wmasscut_biggraph->SetTitle("QCD Efficiencies for Wmasscut 2#times3 Method (160 < m < 240)");
  qcd_efficiency_mass_2groups_improved_Wmasscut_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  qcd_efficiency_mass_2groups_improved_Wmasscut_biggraph->GetYaxis()->SetTitle("% mass");
// qcd_efficiency_mass_2groups_improved_Wmasscut_biggraph->GetYaxis()->SetTitleOffset(1.3);
  qcd_mass_2groups_improved_Wmasscut_can->Write();
  qcd_mass_2groups_improved_Wmasscut_can->Print("ttbarstudy_2x3_plots/qcd_efficiency_mass_2groups_improved_Wmasscut.eps", "eps");


  TCanvas *qcd_mass_3jet4jet_2groups_can = new TCanvas("qcd_mass_3jet4jet_2groups_can", "qcd_mass_3jet4jet_2groups_can", 600, 600);
  qcd_mass_3jet4jet_2groups_can->cd();
  qcd_efficiency_mass_3jet4jet_2groups_biggraph->Draw("ALP");
  leg_qcd_mass_3jet4jet_2groups->Draw();
  qcd_efficiency_mass_3jet4jet_2groups_biggraph->SetMinimum(0.);     
  qcd_efficiency_mass_3jet4jet_2groups_biggraph->SetMaximum(1.2);
  qcd_efficiency_mass_3jet4jet_2groups_biggraph->SetTitle("QCD Efficiencies for Wmasscut 2#times3 Method (160 < m < 240)");
  qcd_efficiency_mass_3jet4jet_2groups_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  qcd_efficiency_mass_3jet4jet_2groups_biggraph->GetYaxis()->SetTitle("% mass");
// qcd_efficiency_mass_3jet4jet_2groups_biggraph->GetYaxis()->SetTitleOffset(1.3);
  qcd_mass_3jet4jet_2groups_can->Write();
  qcd_mass_3jet4jet_2groups_can->Print("ttbarstudy_2x3_plots/qcd_efficiency_mass_3jet4jet_2groups.eps", "eps");


  TCanvas *qcd_mass_3jet4jet_2groups_Wmasscut_can = new TCanvas("qcd_mass_3jet4jet_2groups_Wmasscut_can", "qcd_mass_3jet4jet_2groups_Wmasscut_can", 600, 600);
  qcd_mass_3jet4jet_2groups_Wmasscut_can->cd();
  qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->Draw("ALP");
  leg_qcd_mass_3jet4jet_2groups_Wmasscut->Draw();
  qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->SetMinimum(0.);     
  qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->SetMaximum(1.2);
  qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->SetTitle("QCD Efficiencies for Wmasscut 2#times3 Method (160 < m < 240)");
  qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->GetYaxis()->SetTitle("% mass");
// qcd_efficiency_mass_3jet4jet_2groups_Wmasscut_biggraph->GetYaxis()->SetTitleOffset(1.3);
  qcd_mass_3jet4jet_2groups_Wmasscut_can->Write();
  qcd_mass_3jet4jet_2groups_Wmasscut_can->Print("ttbarstudy_2x3_plots/qcd_efficiency_mass_3jet4jet_2groups_Wmasscut.eps", "eps");


  TCanvas *qcd_ratio_can = new TCanvas("qcd_ratio_can", "qcd_ratio_can", 800, 600);
  qcd_ratio_can->cd();
  qcd_efficiency_ratio_biggraph->Draw("ALP");
  leg_qcd_ratio->Draw();
  qcd_efficiency_ratio_biggraph->SetMinimum(0.0);     
  qcd_efficiency_ratio_biggraph->SetMaximum(2.0);
  qcd_efficiency_ratio_biggraph->SetTitle("Ratio of 3-jettiness/2-jettiness qcd Efficiencies");
  qcd_efficiency_ratio_biggraph->GetXaxis()->SetTitle("p_{T, min} (GeV)");
  qcd_efficiency_ratio_biggraph->GetYaxis()->SetTitle("3-jet eff/2-jet eff");
  qcd_efficiency_ratio_biggraph->GetYaxis()->SetTitleOffset(1.3);
  qcd_ratio_can->Write();
  qcd_ratio_can->Print("ttbarstudy_2x3_plots/qcd_efficiency_ratio.eps", "eps");

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