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

class GeneralERecombiner : public fastjet::JetDefinition::Recombiner {
public:
// Constructor to choose value of alpha (defaulted to 1 for normal pT sum)
  GeneralERecombiner(double delta) : _delta(delta) {}
  
  std::string description() const {
   return "General E-scheme recombination";
 }

 void recombine(const fastjet::PseudoJet & pa, const fastjet::PseudoJet & pb, fastjet::PseudoJet & pab) const {

  double weighta = pow(pa.perp(), _delta);
  double weightb = pow(pb.perp(), _delta);

  double perp_ab = pa.perp() + pb.perp();
    if (perp_ab != 0.0) { // weights also non-zero...
      double y_ab = (weighta * pa.rap() + weightb * pb.rap())/(weighta+weightb);
      
      // take care with periodicity in phi...
      double phi_a = pa.phi(), phi_b = pb.phi();
      if (phi_a - phi_b > pi)  phi_b += twopi;
      if (phi_a - phi_b < -pi) phi_b -= twopi;
      double phi_ab = (weighta * phi_a + weightb * phi_b)/(weighta+weightb);
      
      pab.reset_PtYPhiM(perp_ab, y_ab, phi_ab);

    }
    else { // weights are zero
      pab.reset(0.0,0.0,0.0,0.0);
    }

  }

private:
  double _delta;
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
}

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

    double manual_tau = 10000;
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

int main(int argc, char* argv[]){

  TFile out("ttbarstudy.root", "RECREATE");
  // TFile out("njettiness_testing_garbage.root", "RECREATE");

  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW"); // Input File Name
  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW"); // Input File Name

  double epsilon = 0.0001;
  double Rparam = 0.4;

  int initial_axes_num = 6;
  int antikt_initial_axes_num = 3;
  double ghost_maxrap = 5.0; // e.g. if particles go up to y=5

  vector<double> mass_ratios;
  vector<double> mean_values;
  vector<double> rms_values;

  //create list of various values of beta
  vector<double> betalist;
  betalist.push_back(0.25);
  betalist.push_back(0.5);
  betalist.push_back(1.0);
  betalist.push_back(2.0);
  // betalist.push_back(3.0);
  int n_betas = betalist.size();

  int total_ttbar_jets = 0;
  int total_dijets_jets = 0;

  TObjArray rawmass_ttbar_manual_hists;
  TObjArray rawmass_dijets_manual_hists;

  TObjArray area_ttbar_manual_hists;
  TObjArray area_dijets_manual_hists;

  TObjArray mass_ttbar_manual_2groups_hists;
  TObjArray mass_dijets_manual_2groups_hists;

  TObjArray mass_ttbar_manual_hists;
  TObjArray mass_dijets_manual_hists;

  TObjArray tau32_ttbar_manual_hists;
  TObjArray tau32_dijets_manual_hists;

  TObjArray max_njet_ttbar_manual_2groups_hists;
  TObjArray max_njet_dijets_manual_2groups_hists;

  TObjArray max_njet_ttbar_manual_hists;
  TObjArray max_njet_dijets_manual_hists;

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1* rawmass_ttbar_manual = new TH1F("rawmass_ttbar_manual", "Raw Mass of reconstructed tops", 100, 0, 500);
    TH1* rawmass_dijets_manual = new TH1F("rawmass_dijets_manual", "Mass of QCD jet", 100, 0, 500);
    TH1* area_ttbar_manual = new TH1F("area_ttbar_manual", "Area of top jets", 50, 0, 1);
    TH1* area_dijets_manual = new TH1F("area_dijets_manual", "Area of QCD jet", 50, 0, 1);
    TH1* mass_ttbar_manual_2groups = new TH1F("mass_ttbar_manual_2groups", "Mass of reconstructed tops from 2 groups", 100, 0, 500);
    TH1* mass_dijets_manual_2groups = new TH1F("mass_dijets_manual_2groups", "Mass of QCD jet from 2 groups", 100, 0, 500);
    TH1* mass_ttbar_manual = new TH1F("mass_ttbar_manual", "Mass of reconstructed tops from 3+3 grouping", 100, 0, 500);
    TH1* mass_dijets_manual = new TH1F("mass_dijets_manual", "Mass of QCD jet from 3+3 grouping", 100, 0, 500);
    TH1* tau32_ttbar_manual = new TH1F("tau32_ttbar_manual", "#tau_{32} of Tops", 50, 0, 1);
    TH1* tau32_dijets_manual = new TH1F("tau32_dijets_manual", "#tau_{32} of QCD", 50, 0, 1);
    TH1* max_njet_ttbar_manual_2groups = new TH1F("max_njet_ttbar_manual_2groups", "Maximum N-jettiness of Tops from 2 groups", 20, 0, 20);
    TH1* max_njet_dijets_manual_2groups = new TH1F("max_njet_dijets_manual_2groups", "Maximum N-jettiness of QCD from 2 groups", 20, 0, 20);
    TH1* max_njet_ttbar_manual = new TH1F("max_njet_ttbar_manual", "Maximum N-jettiness of Tops", 20, 0, 20);
    TH1* max_njet_dijets_manual = new TH1F("max_njet_dijets_manual", "Maximum N-jettiness of QCD", 20, 0, 20);

    rawmass_ttbar_manual_hists.Add(rawmass_ttbar_manual);
    rawmass_dijets_manual_hists.Add(rawmass_dijets_manual);
    area_ttbar_manual_hists.Add(area_ttbar_manual);
    area_dijets_manual_hists.Add(area_dijets_manual);
    mass_ttbar_manual_2groups_hists.Add(mass_ttbar_manual_2groups);
    mass_dijets_manual_2groups_hists.Add(mass_dijets_manual_2groups);
    mass_ttbar_manual_hists.Add(mass_ttbar_manual);
    mass_dijets_manual_hists.Add(mass_dijets_manual);
    tau32_ttbar_manual_hists.Add(tau32_ttbar_manual);
    tau32_dijets_manual_hists.Add(tau32_dijets_manual);
    max_njet_ttbar_manual_2groups_hists.Add(max_njet_ttbar_manual_2groups);
    max_njet_dijets_manual_2groups_hists.Add(max_njet_dijets_manual_2groups);
    max_njet_ttbar_manual_hists.Add(max_njet_ttbar_manual);
    max_njet_dijets_manual_hists.Add(max_njet_dijets_manual);
  } 

  // TH1F *tau6_ttbar_manual = new TH1F("tau6_ttbar_manual", "#tau_{6} for ttbar (Manual)", 200, 0, 1000);
  // TH1F *tau6_dijets_manual = new TH1F("tau6_dijets_manual", "#tau_{6} for dijets (Manual)", 200, 0, 1000);

  // TH1F *tau3_ttbar_manual_antikt = new TH1F("tau3_ttbar_manual_antikt", "#tau_{3} for ttbar (Manual)", 200, 0, 1000);
  // TH1F *tau3_dijets_manual_antikt = new TH1F("tau3_dijets_manual_antikt", "#tau_{3} for dijets (Manual)", 200, 0, 1000);

  TH1F *tau32_ttbar_wta = new TH1F("tau32_ttbar_wta", "#tau_{32} for ttbar (WTA)", 50, 0, 1);
  TH1F *tau32_dijets_wta = new TH1F("tau32_dijets_wta", "#tau_{32} for dijets (WTA)", 50, 0, 1);

  TH1F *mass_ttbar_antikt = new TH1F("mass_ttbar_antikt", "mass ttbar", 100, 0, 500);
  TH1F *mass_dijets_antikt = new TH1F("mass_dijets_antikt", "mass dijets", 100, 0, 500);

  std::string samples[2] = {"/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW", "/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW"};
  // std::string samples[2] = {"/home/tjwilk/Physics/thaler_urop/BOOST_samples/pythia64-tuneDW-lhc7-ttbar2hadrons-pt0500-0600.UW", "/home/tjwilk/Physics/thaler_urop/BOOST_samples/pythia64-tuneDW-lhc7-dijets-pt0500-0600.UW"};

  for (int i_sample = 0; i_sample < 2; i_sample++) {

    const char* current_data = samples[i_sample].c_str();
    ifstream inputStream(current_data); // Input File Name

    const int n_event = 1000; // # of events to be analyzed   
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

        int max_njettiness;

        TH2F *event_display;
        if (i_sample == 0) event_display = new TH2F("event_display", "TTbar event (500 < pT < 600)", 60, -5, 5, 60, 0, 6.4);
        if (i_sample == 1) event_display = new TH2F("event_display", "QCD event (500 < pT < 600)", 60, -5, 5, 60, 0, 6.4);
        TH2F *axes_manual_display = new TH2F("axes_manual_display", "Axes Plot", 300, -5, 5, 300, 0, 6.4);
        TH2F *subjet_display = new TH2F("subjet_display", "Subjet display", 1000, -5, 5, 1000, 0, 6.4);

        vector<PseudoJet> jet_constituents = input_particles;

        for (int j = 0; j < jet_constituents.size(); j++) {
          event_display->Fill(jet_constituents[j].eta(), jet_constituents[j].phi(), jet_constituents[j].perp()*100);
        }

        vector<PseudoJet> final_manual_axes;
    // vector<PseudoJet> final_manual_jets;

        for (unsigned int B = 0; B < betalist.size(); B++) {

          double beta = betalist[B];
          double power = (double)1/beta;
      // double delta;

          const JetDefinition::Recombiner *recombScheme;
          if (beta > 1) recombScheme = new GeneralERecombiner((double)1/(beta - 1));
          else recombScheme = new WinnerTakeAllRecombiner();

          UnnormalizedCutoffMeasure measure_function = UnnormalizedCutoffMeasure(beta, Rparam);
          OnePass_Manual_Axes axes_finder = OnePass_Manual_Axes();

          vector<PseudoJet> final_bigjets;
          vector<PseudoJet> initial_jets;

          int njettiness = 6;

      // Manual Axis Finding for general KT and general Recombiner
      // double Rparam = JetDefinition::max_allowable_R;
      // double Rparam = Rcutoff;
          Strategy strategy = Best;
          JetDefinition *jetDef = new JetDefinition(genkt_algorithm, Rparam, power, recombScheme, strategy);

          while (njettiness <= 20) {

            int initial_axes_num;
            if (njettiness < 20) initial_axes_num = njettiness + 1;
            else initial_axes_num = njettiness;

            ClusterSequence clustSeq(jet_constituents, *jetDef);
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

            AreaDefinition area_def(active_area, GhostedAreaSpec(ghost_maxrap));
            ClusterSequenceArea njet_cluster_min_manual_area(jet_constituents, njet_def_min_manual, area_def);
            ClusterSequence njet_cluster_min_manual(jet_constituents, njet_def_min_manual);
            const NjettinessExtras *extras_min_manual = njettiness_extras(njet_cluster_min_manual);
            vector<PseudoJet> min_manual_jets = njet_cluster_min_manual_area.inclusive_jets();
        // vector<PseudoJet> min_manual_jets = extras_min_manual->jets();
            double min_manual_tau = extras_min_manual->totalTau();

            vector<PseudoJet> temp_bigjets(min_manual_jets.size());
            int index_references[min_manual_jets.size()];

        // initialize the index references to point to their own index
            for (int i_jets = 0; i_jets < min_manual_jets.size(); i_jets++) {
              index_references[i_jets] = i_jets;
            }

            double grouping_radius = 1.5;

        // run through every combination of jets and reset index reference to the index of the closest jet with the lowest index
            for (int j_jets = 0; j_jets < min_manual_jets.size(); j_jets++) {
              for (int k_jets = j_jets + 1; k_jets < min_manual_jets.size(); k_jets++) {
                if (min_manual_axes[j_jets].delta_R(min_manual_axes[k_jets]) < grouping_radius) {
                  if (index_references[k_jets] < index_references[j_jets]) index_references[j_jets] = index_references[k_jets];
                  else if (index_references[j_jets] < index_references[k_jets]) index_references[k_jets] = index_references[j_jets];
                }
              }
            }

        // group the N subjets according to the index references from above
            vector<int> subjet_counter(min_manual_jets.size());
            for (int i = 0; i < min_manual_jets.size(); i++) {
              if (!temp_bigjets[index_references[i]].has_structure()) {
                temp_bigjets[index_references[i]] = min_manual_jets[i];
                subjet_counter[index_references[i]] = 1;
              }
              else {
                temp_bigjets[index_references[i]] = join(temp_bigjets[index_references[i]], min_manual_jets[i]);
                subjet_counter[index_references[i]]++;
              }
            }

        // pick out only the combined jets that have mass (i.e. remove empty jets from the vector)
            vector<PseudoJet> massive_jets;
            vector<int> massive_subjet_counter;
            for (int i = 0; i < temp_bigjets.size(); i++) {
              if (temp_bigjets[i].m() > 0) massive_jets.push_back(temp_bigjets[i]);
              if (subjet_counter[i] > 0) massive_subjet_counter.push_back(subjet_counter[i]);
            }

        //count how many of the jets have 3 subjets
            int triplet_counter = 0;
            if ((beta - 1) < epsilon && njettiness >= 15) cout << "event: " << i_event << " " << endl;
            for (int i = 0; i < massive_subjet_counter.size(); i++) {
              if ((beta - 1) < epsilon && njettiness >= 15) {
                cout << massive_subjet_counter[i] << " ";
              }
              if (massive_subjet_counter[i] >= 3) triplet_counter++;
            }
            if ((beta - 1) < epsilon && njettiness >= 15) cout << endl;

            if (njettiness == 6) {

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
                if (i_sample == 0) {
                  TH1* rawmass_ttbar_manual = (TH1*)rawmass_ttbar_manual_hists.At(B);
                  rawmass_ttbar_manual->Fill(final_bigjets[i_jet].m());
                }
                if (i_sample == 1) {
                  TH1* rawmass_dijets_manual = (TH1*)rawmass_dijets_manual_hists.At(B);
                  rawmass_dijets_manual->Fill(final_bigjets[i_jet].m());
                }
              }
            }

        // if there are at least 2 jets with 3 subjets, fill the histogram with the two hardest jets
            if (triplet_counter >= 2 || njettiness == 20) {
              Selector hard_selector = SelectorNHardest(2);
              final_bigjets = hard_selector(massive_jets);
              for (int i_jet = 0; i_jet < final_bigjets.size(); i_jet++) {
                if (i_sample == 0) {
                  TH1* mass_ttbar_manual = (TH1*)mass_ttbar_manual_hists.At(B);
                  mass_ttbar_manual->Fill(final_bigjets[i_jet].m());
                }
                if (i_sample == 1) {
                  TH1* mass_dijets_manual = (TH1*)mass_dijets_manual_hists.At(B);
                  mass_dijets_manual->Fill(final_bigjets[i_jet].m());
                }
              }
              if (i_sample == 0) {
                TH1* max_njet_ttbar_manual = (TH1*)max_njet_ttbar_manual_hists.At(B);
                max_njet_ttbar_manual->Fill(njettiness);
              }
              if (i_sample == 1) {
                TH1* max_njet_dijets_manual = (TH1*)max_njet_dijets_manual_hists.At(B);
                max_njet_dijets_manual->Fill(njettiness);
              }
              final_manual_axes = min_manual_axes;
              max_njettiness = njettiness;

              break;
            }
        // final_manual_jets = min_manual_jets;

            njettiness++;
          }

          if(abs(beta - 1) < epsilon && max_njettiness >= 15) {

            Double_t ghost_perp = 0.001;
            Double_t n_ghosts = 50;
            for (int i_ghosts = 0; i_ghosts < n_ghosts; i_ghosts++) {
              for (int j_ghosts = 0; j_ghosts < n_ghosts; j_ghosts++) {
                Double_t ghost_eta = -5 + (double)(10/n_ghosts)*i_ghosts;
                Double_t ghost_phi = (double)(2*TMath::Pi()/n_ghosts)*j_ghosts;
                PseudoJet ghost(ghost_perp*TMath::Cos(ghost_phi),ghost_perp*TMath::Sin(ghost_phi),
                  ghost_perp*TMath::SinH(ghost_eta),ghost_perp*TMath::CosH(ghost_eta));
                jet_constituents.push_back(ghost);
              }
            } 

            NjettinessPlugin njet_plugin_min_manual(max_njettiness, axes_finder, measure_function);
            JetDefinition njet_def_min_manual(&njet_plugin_min_manual);
            njet_plugin_min_manual.setAxes(final_manual_axes);
            ClusterSequence njet_cluster_min_manual(jet_constituents, njet_def_min_manual);
            const NjettinessExtras *extras_min_manual = njettiness_extras(njet_cluster_min_manual);
            vector<PseudoJet> final_manual_jets = extras_min_manual->jets();
            final_manual_axes = extras_min_manual->axes();

            for (int a = 0; a < final_manual_axes.size(); a++) {
              axes_manual_display->Fill(final_manual_axes[a].eta(), final_manual_axes[a].phi());
            }

            for (int i = 0; i < final_manual_jets.size(); i++) {
              vector<PseudoJet> constituents = final_manual_jets[i].constituents();

              for (int i_const = 0; i_const < constituents.size(); i_const++) {
                subjet_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
              }
            }
          }

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

            ClusterSequence clustSeq_subjet(hardestJets_std[i].constituents(), *jetDef);
            vector<PseudoJet> subjet_axes3_start = clustSeq_subjet.exclusive_jets(4);
            vector<PseudoJet> subjet_axes2_start = clustSeq_subjet.exclusive_jets(3);
            vector<PseudoJet> subjet_axes3 = findMinAxes(hardestJets_std[i].constituents(), subjet_axes3_start, 3, beta, Rparam);
            vector<PseudoJet> subjet_axes2 = findMinAxes(hardestJets_std[i].constituents(), subjet_axes2_start, 2, beta, Rparam);

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

            NsubjettinessRatio nsub_32_wta(3, 2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta));
            if (hardestJets_std[i].m() > 160 && hardestJets_std[i].m() < 240) {
              if (i_sample == 0) {
                TH1* tau32_ttbar_manual = (TH1*)tau32_ttbar_manual_hists.At(B);
                tau32_ttbar_manual->Fill(nsub3_manual.result(hardestJets_std[i])/nsub2_manual.result(hardestJets_std[i]));
                if (abs(beta - 1) < epsilon) tau32_ttbar_wta->Fill(nsub_32_wta.result(hardestJets_std[i]));
              }
              if (i_sample == 1) {
                TH1* tau32_dijets_manual = (TH1*)tau32_dijets_manual_hists.At(B);
                tau32_dijets_manual->Fill(nsub3_manual.result(hardestJets_std[i])/nsub2_manual.result(hardestJets_std[i]));
                if (abs(beta - 1) < epsilon) tau32_dijets_wta->Fill(nsub_32_wta.result(hardestJets_std[i]));
              }
            }
          }
        }

        event_display->SetStats(0);
        axes_manual_display->SetStats(0);
        subjet_display->SetStats(0);
        // subjet2_display->SetStats(0);
        // subjet3_display->SetStats(0);
        // subjet4_display->SetStats(0);
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
        event_display->SetLineWidth(1);
        subjet_display->SetMarkerColor(kRed);
        // subjet2_display->SetMarkerColor(kBlue);
        // subjet3_display->SetMarkerColor(kGreen);
        // subjet4_display->SetMarkerColor(kYellow);
        // subjet5_display->SetMarkerColor(kOrange);
        // subjet6_display->SetMarkerColor(6);
        // subjet7_display->SetMarkerColor(20);
        subjet_display->SetMarkerStyle(21);
        // subjet2_display->SetMarkerStyle(21);
        // subjet3_display->SetMarkerStyle(21);
        // subjet4_display->SetMarkerStyle(21);
        // subjet5_display->SetMarkerStyle(21);
        // subjet6_display->SetMarkerStyle(21);
        // subjet7_display->SetMarkerStyle(21);
        subjet_display->SetMarkerSize(0.5);
        // subjet2_display->SetMarkerSize(0.5);
        // subjet3_display->SetMarkerSize(0.5);
        // subjet4_display->SetMarkerSize(0.5);
        // subjet5_display->SetMarkerSize(0.5);
        // subjet6_display->SetMarkerSize(0.5);
        // subjet7_display->SetMarkerSize(0.5);
        event_display->Draw("box");

        axes_manual_display->SetMarkerStyle(3);
        axes_manual_display->SetMarkerSize(3);
        axes_manual_display->SetMarkerColor(kGreen);
        axes_manual_display->Draw("SAMES");

        subjet_display->Draw("SAMES");
        // subjet2_display->Draw("SAMES");
        // subjet3_display->Draw("SAMES");
        // subjet4_display->Draw("SAMES");
        // subjet5_display->Draw("SAMES");
        // subjet6_display->Draw("SAMES");
        // subjet7_display->Draw("SAMES");

        if (final_manual_axes.size() >= 15) display->Write();

        // ostringstream ss;
        // ss << i_event;
        // TString title;
        // if (i_sample == 0) title = "ttbar_display" + ss.str() + "_6jets_test_nplus1.pdf";
        // if (i_sample == 1) title = "dijets_display" + ss.str() + "_6jets_test_nplus1.pdf";

        // display->Print(title, "pdf");

        delete display;
        delete event_display;
        delete axes_manual_display;
        delete subjet_display;
        // delete subjet2_display;
        // delete subjet3_display;
        // delete subjet4_display;
        // delete subjet5_display;
        // delete subjet6_display;
        // delete subjet7_display;

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

  mass_ttbar_antikt->SetStats(0);
  mass_dijets_antikt->SetStats(0);
  tau32_ttbar_wta->SetStats(0);
  tau32_dijets_wta->SetStats(0);

  mass_ttbar_antikt->Write();
  mass_dijets_antikt->Write();
  tau32_ttbar_wta->Write();
  tau32_dijets_wta->Write();

  int n_wtakt = tau32_ttbar_wta->GetSize() - 2;
  double integral_dijets_wtakt[n_wtakt], integral_ttbar_wtakt[n_wtakt];
  for (int i = 0; i < n_wtakt; i++) {
    integral_ttbar_wtakt[i] = tau32_ttbar_wta->Integral(0,i)/total_ttbar_jets;
    integral_dijets_wtakt[i] = tau32_dijets_wta->Integral(0,i)/total_dijets_jets;
  }
  TGraph* ROC_tau32_std = new TGraph(n_wtakt, integral_ttbar_wtakt, integral_dijets_wtakt);
  ROC_tau32_std->GetXaxis()->SetLimits(0, 1);
  ROC_tau32_std->GetYaxis()->SetLimits(0, 1);
  ROC_tau32_std->SetLineColor(kBlack);
  ROC_tau32_std->SetLineWidth(2);
  ROC_tau32_std->SetMarkerStyle(5);
  ROC_tau32_std->SetMarkerSize(2);
  ROC_tau32_std->Write();

  TCanvas *ROC_compare = new TCanvas("ROC_compare", "ROC Curve Comparison", 600, 600);
  ROC_compare->cd();
  ROC_compare->SetLogy();
  TMultiGraph *ROC_multigraph = new TMultiGraph("ROC_multigraph", "ROC Comparison for #tau_{3}/#tau_{2} (500 < p_{T} < 600)");
  ROC_multigraph->Add(ROC_tau32_std);
  TLegend *leg_ROC = new TLegend(0.1, 0.6, 0.4, 0.9);
  leg_ROC->AddEntry(ROC_tau32_std, "WTA kT", "L");

  for (int B = 0; B < betalist.size(); B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1F* tau32_ttbar_manual = (TH1F*)tau32_ttbar_manual_hists.At(B);
    TH1F* tau32_dijets_manual = (TH1F*)tau32_dijets_manual_hists.At(B);

    tau32_ttbar_manual->SetStats(0);
    tau32_dijets_manual->SetStats(0);

    tau32_ttbar_manual->Write();
    tau32_dijets_manual->Write();

    int n_size = tau32_ttbar_manual->GetSize() - 2;
    double integral_dijets[n_size], integral_ttbar[n_size];
    for (int i = 0; i < n_size; i++) {
      integral_ttbar[i] = tau32_ttbar_manual->Integral(0,i)/total_ttbar_jets;
      integral_dijets[i] = tau32_dijets_manual->Integral(0,i)/total_dijets_jets;
    }
    TGraph* ROC_tau32 = new TGraph(n_size, integral_ttbar, integral_dijets);
    ROC_tau32->GetXaxis()->SetLimits(0, 1);
    ROC_tau32->GetYaxis()->SetLimits(0, 1);
    if (B == 0) ROC_tau32->SetLineColor(kGreen);
    if (B == 1) ROC_tau32->SetLineColor(kYellow);
    if (B == 2) ROC_tau32->SetLineColor(kBlue);
    if (B == 3) ROC_tau32->SetLineColor(kRed);
    ROC_tau32->SetLineWidth(2);
    ROC_tau32->SetMarkerStyle(5);
    ROC_tau32->SetMarkerSize(2);
    ROC_multigraph->Add(ROC_tau32);
    leg_ROC->AddEntry(ROC_tau32, "Manual Axes (beta = " + (TString)ss.str() + ")", "L");
    // ROC_tau32->Write();

  }

  ROC_multigraph->Draw("AL");
  leg_ROC->Draw();
  ROC_compare->Write();
  ROC_compare->Print("ttbarstudy_ROC_compare_nplus1.pdf", "pdf");

  TCanvas *mass_ttbar_compare = new TCanvas("mass_ttbar_compare", "Mass Comparison", 600, 600);
  mass_ttbar_compare->cd();
  double mass_ttbar_antikt_scale = 1/mass_ttbar_antikt->Integral();
  mass_ttbar_antikt->Scale(mass_ttbar_antikt_scale);
  mass_ttbar_antikt->SetLineColor(kBlack);
  mass_ttbar_antikt->Draw();
  TLegend *leg_mass_ttbar = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg_mass_ttbar->AddEntry(mass_ttbar_antikt, "Anti_KT", "L");
  double max_val_ttbar = mass_ttbar_antikt->GetMaximum();


  for (int B = 0; B < betalist.size(); B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;
    
    TH1* mass_ttbar_manual = (TH1*)mass_ttbar_manual_hists.At(B);
    TH1* mass_dijets_manual = (TH1*)mass_dijets_manual_hists.At(B);

    double mass_ttbar_manual_scale = 1/mass_ttbar_manual->Integral();
    double mass_dijets_manual_scale = 1/mass_dijets_manual->Integral();
    mass_ttbar_manual->Scale(mass_ttbar_manual_scale);
    mass_dijets_manual->Scale(mass_dijets_manual_scale);

    mass_ttbar_manual->SetStats(0);
    mass_dijets_manual->SetStats(0);

    if (B == 0) mass_ttbar_manual->SetLineColor(kGreen);
    if (B == 1) mass_ttbar_manual->SetLineColor(kYellow);
    if (B == 2) mass_ttbar_manual->SetLineColor(kBlue);
    if (B == 3) mass_ttbar_manual->SetLineColor(kRed);

    mass_ttbar_manual->Draw("SAMES");
    leg_mass_ttbar->AddEntry(mass_ttbar_manual, "beta = " + (TString)ss.str());
    mass_ttbar_manual->Write();

    if (mass_ttbar_manual->GetMaximum() > max_val_ttbar) {
      max_val_ttbar = mass_ttbar_manual->GetMaximum();
      mass_ttbar_manual->SetMaximum(max_val_ttbar);
    }

    cout << "manual ttbar mass ratio: " << mass_ttbar_manual->Integral(75, 110) << "(beta = " << beta << ")" << endl;
  }

  leg_mass_ttbar->Draw();
  mass_ttbar_compare->Write();
  mass_ttbar_compare->Print("ttbarstudy_ttbar_mass_compare_nplus1.pdf", "pdf");

  TCanvas *mass_dijets_compare = new TCanvas("mass_dijets_compare", "Mass Comparison", 600, 600);
  mass_dijets_compare->cd();
  double mass_dijets_antikt_scale = 1/mass_dijets_antikt->Integral();
  mass_dijets_antikt->Scale(mass_dijets_antikt_scale);
  mass_dijets_antikt->SetLineColor(kBlack);
  mass_dijets_antikt->Draw();
  TLegend *leg_mass_dijets = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg_mass_dijets->AddEntry(mass_dijets_antikt, "Anti_KT", "L");
  double max_val_dijets = mass_dijets_antikt->GetMaximum();

  for (int B = 0; B < betalist.size(); B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;
    
    TH1* mass_dijets_manual = (TH1*)mass_dijets_manual_hists.At(B);

    double mass_dijets_manual_scale = 1/mass_dijets_manual->Integral();
    mass_dijets_manual->Scale(mass_dijets_manual_scale);

    mass_dijets_manual->SetStats(0);

    if (B == 0) mass_dijets_manual->SetLineColor(kGreen);
    if (B == 1) mass_dijets_manual->SetLineColor(kYellow);
    if (B == 2) mass_dijets_manual->SetLineColor(kBlue);
    if (B == 3) mass_dijets_manual->SetLineColor(kRed);

    mass_dijets_manual->Draw("SAMES");
    leg_mass_dijets->AddEntry(mass_dijets_manual, "beta = " + (TString)ss.str());
    mass_dijets_manual->Write();

    if (mass_dijets_manual->GetMaximum() > max_val_dijets) {
      max_val_dijets = mass_dijets_manual->GetMaximum();
      mass_dijets_manual->SetMaximum(max_val_dijets);
    }

    cout << "manual QCD mass ratio: " << mass_dijets_manual->Integral(75, 110) << "(beta = " << beta << ")" << endl;

  }

  leg_mass_dijets->Draw();
  mass_dijets_compare->Write();
  mass_dijets_compare->Print("ttbarstudy_dijets_mass_compare_nplus1.pdf", "pdf");


  TCanvas *area_ttbar_compare = new TCanvas("area_ttbar_compare", "area Comparison", 600, 600);
  area_ttbar_compare->cd();
  TLegend *leg_area_ttbar = new TLegend(0.1, 0.6, 0.4, 0.9);
  double max_val_ttbar_area = 0;

  for (int B = 0; B < betalist.size(); B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;
    
    TH1* area_ttbar_manual = (TH1*)area_ttbar_manual_hists.At(B);
    TH1* area_dijets_manual = (TH1*)area_dijets_manual_hists.At(B);

    double area_ttbar_manual_scale = 1/area_ttbar_manual->Integral();
    double area_dijets_manual_scale = 1/area_dijets_manual->Integral();
    area_ttbar_manual->Scale(area_ttbar_manual_scale);
    area_dijets_manual->Scale(area_dijets_manual_scale);

    area_ttbar_manual->SetStats(0);
    area_dijets_manual->SetStats(0);

    if (B == 0) area_ttbar_manual->SetLineColor(kGreen);
    if (B == 1) area_ttbar_manual->SetLineColor(kYellow);
    if (B == 2) area_ttbar_manual->SetLineColor(kBlue);
    if (B == 3) area_ttbar_manual->SetLineColor(kRed);

    area_ttbar_manual->Draw("SAMES");
    leg_area_ttbar->AddEntry(area_ttbar_manual, "beta = " + (TString)ss.str());
    area_ttbar_manual->Write();

    if (area_ttbar_manual->GetMaximum() > max_val_ttbar_area) {
      max_val_ttbar_area = area_ttbar_manual->GetMaximum();
      area_ttbar_manual->SetMaximum(max_val_ttbar_area);
    }
  }

  leg_area_ttbar->Draw();
  area_ttbar_compare->Write();
  area_ttbar_compare->Print("ttbarstudy_ttbar_area_compare_nplus1.pdf", "pdf");

  TCanvas *area_dijets_compare = new TCanvas("area_dijets_compare", "area Comparison", 600, 600);
  area_dijets_compare->cd();
  TLegend *leg_area_dijets = new TLegend(0.1, 0.6, 0.4, 0.9);
  double max_val_dijets_area = 0;

  for (int B = 0; B < betalist.size(); B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;
    
    TH1* area_dijets_manual = (TH1*)area_dijets_manual_hists.At(B);

    double area_dijets_manual_scale = 1/area_dijets_manual->Integral();
    area_dijets_manual->Scale(area_dijets_manual_scale);

    area_dijets_manual->SetStats(0);

    if (B == 0) area_dijets_manual->SetLineColor(kGreen);
    if (B == 1) area_dijets_manual->SetLineColor(kYellow);
    if (B == 2) area_dijets_manual->SetLineColor(kBlue);
    if (B == 3) area_dijets_manual->SetLineColor(kRed);

    area_dijets_manual->Draw("SAMES");
    leg_area_dijets->AddEntry(area_dijets_manual, "beta = " + (TString)ss.str());
    area_dijets_manual->Write();

    if (area_dijets_manual->GetMaximum() > max_val_dijets_area) {
      max_val_dijets_area = area_dijets_manual->GetMaximum();
      area_dijets_manual->SetMaximum(max_val_dijets_area);
    }

  }

  leg_area_dijets->Draw();
  area_dijets_compare->Write();
  area_dijets_compare->Print("ttbarstudy_dijets_area_compare_nplus1.pdf", "pdf");


  TCanvas *rawmass_ttbar_compare = new TCanvas("rawmass_ttbar_compare", "rawmass Comparison", 600, 600);
  rawmass_ttbar_compare->cd();
  mass_ttbar_antikt->SetLineColor(kBlack);
  mass_ttbar_antikt->Draw();
  TLegend *leg_rawmass_ttbar = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg_rawmass_ttbar->AddEntry(mass_ttbar_antikt, "Anti_KT", "L");
  double max_val_ttbar_raw = mass_ttbar_antikt->GetMaximum();

  for (int B = 0; B < betalist.size(); B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;
    
    TH1* rawmass_ttbar_manual = (TH1*)rawmass_ttbar_manual_hists.At(B);
    TH1* rawmass_dijets_manual = (TH1*)rawmass_dijets_manual_hists.At(B);

    double rawmass_ttbar_manual_scale = 1/rawmass_ttbar_manual->Integral();
    double rawmass_dijets_manual_scale = 1/rawmass_dijets_manual->Integral();
    rawmass_ttbar_manual->Scale(rawmass_ttbar_manual_scale);
    rawmass_dijets_manual->Scale(rawmass_dijets_manual_scale);

    rawmass_ttbar_manual->SetStats(0);
    rawmass_dijets_manual->SetStats(0);

    if (B == 0) rawmass_ttbar_manual->SetLineColor(kGreen);
    if (B == 1) rawmass_ttbar_manual->SetLineColor(kYellow);
    if (B == 2) rawmass_ttbar_manual->SetLineColor(kBlue);
    if (B == 3) rawmass_ttbar_manual->SetLineColor(kRed);

    rawmass_ttbar_manual->Draw("SAMES");
    leg_rawmass_ttbar->AddEntry(rawmass_ttbar_manual, "beta = " + (TString)ss.str());
    rawmass_ttbar_manual->Write();

    if (rawmass_ttbar_manual->GetMaximum() > max_val_ttbar_raw) {
      max_val_ttbar_raw = rawmass_ttbar_manual->GetMaximum();
      rawmass_ttbar_manual->SetMaximum(max_val_ttbar_raw);
    }

    cout << "manual ttbar rawmass ratio: " << rawmass_ttbar_manual->Integral(75, 110) << "(beta = " << beta << ")" << endl;
  }

  leg_rawmass_ttbar->Draw();
  rawmass_ttbar_compare->Write();
  rawmass_ttbar_compare->Print("ttbarstudy_ttbar_rawmass_compare_nplus1.pdf", "pdf");

  TCanvas *rawmass_dijets_compare = new TCanvas("rawmass_dijets_compare", "rawmass Comparison", 600, 600);
  rawmass_dijets_compare->cd();
  mass_dijets_antikt->SetLineColor(kBlack);
  mass_dijets_antikt->Draw();
  TLegend *leg_rawmass_dijets = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg_rawmass_dijets->AddEntry(mass_dijets_antikt, "Anti_KT", "L");
  double max_val_dijets_raw = mass_dijets_antikt->GetMaximum();

  for (int B = 0; B < betalist.size(); B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;
    
    TH1* rawmass_dijets_manual = (TH1*)rawmass_dijets_manual_hists.At(B);

    double rawmass_dijets_manual_scale = 1/rawmass_dijets_manual->Integral();
    rawmass_dijets_manual->Scale(rawmass_dijets_manual_scale);

    rawmass_dijets_manual->SetStats(0);

    if (B == 0) rawmass_dijets_manual->SetLineColor(kGreen);
    if (B == 1) rawmass_dijets_manual->SetLineColor(kYellow);
    if (B == 2) rawmass_dijets_manual->SetLineColor(kBlue);
    if (B == 3) rawmass_dijets_manual->SetLineColor(kRed);

    rawmass_dijets_manual->Draw("SAMES");
    leg_rawmass_dijets->AddEntry(rawmass_dijets_manual, "beta = " + (TString)ss.str());
    rawmass_dijets_manual->Write();

    if (rawmass_dijets_manual->GetMaximum() > max_val_dijets_raw) {
      max_val_dijets_raw = rawmass_dijets_manual->GetMaximum();
      rawmass_dijets_manual->SetMaximum(max_val_dijets_raw);
    }

    cout << "manual QCD rawmass ratio: " << rawmass_dijets_manual->Integral(75, 110) << "(beta = " << beta << ")" << endl;

  }

  leg_rawmass_dijets->Draw();
  rawmass_dijets_compare->Write();
  rawmass_dijets_compare->Print("ttbarstudy_dijets_rawmass_compare_nplus1.pdf", "pdf");

  TCanvas *max_njet_ttbar_compare = new TCanvas("max_njet_ttbar_compare", "max_njet Comparison", 600, 600);
  max_njet_ttbar_compare->cd();
  TLegend *leg_max_njet_ttbar = new TLegend(0.6, 0.6, 0.9, 0.9);
  double max_val_ttbar_njet = 0;

  for (int B = 0; B < betalist.size(); B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;
    
    TH1* max_njet_ttbar_manual = (TH1*)max_njet_ttbar_manual_hists.At(B);
    TH1* max_njet_dijets_manual = (TH1*)max_njet_dijets_manual_hists.At(B);

    double max_njet_ttbar_manual_scale = 1/max_njet_ttbar_manual->Integral();
    double max_njet_dijets_manual_scale = 1/max_njet_dijets_manual->Integral();
    max_njet_ttbar_manual->Scale(max_njet_ttbar_manual_scale);
    max_njet_dijets_manual->Scale(max_njet_dijets_manual_scale);

    max_njet_ttbar_manual->SetStats(0);
    max_njet_dijets_manual->SetStats(0);

    if (B == 0) max_njet_ttbar_manual->SetLineColor(kGreen);
    if (B == 1) max_njet_ttbar_manual->SetLineColor(kYellow);
    if (B == 2) max_njet_ttbar_manual->SetLineColor(kBlue);
    if (B == 3) max_njet_ttbar_manual->SetLineColor(kRed);

    max_njet_ttbar_manual->Draw("SAMES");
    leg_max_njet_ttbar->AddEntry(max_njet_ttbar_manual, "beta = " + (TString)ss.str());
    max_njet_ttbar_manual->Write();

    if (max_njet_ttbar_manual->GetMaximum() > max_val_ttbar_njet) {
      max_val_ttbar_njet = max_njet_ttbar_manual->GetMaximum();
      max_njet_ttbar_manual->SetMaximum(max_val_ttbar_njet);
    }

    cout << "manual ttbar max_njet ratio: " << max_njet_ttbar_manual->Integral(75, 110) << "(beta = " << beta << ")" << endl;
  }

  leg_max_njet_ttbar->Draw();
  max_njet_ttbar_compare->Write();
  max_njet_ttbar_compare->Print("ttbarstudy_ttbar_max_njet_compare_nplus1.pdf", "pdf");

  TCanvas *max_njet_dijets_compare = new TCanvas("max_njet_dijets_compare", "max_njet Comparison", 600, 600);
  max_njet_dijets_compare->cd();
  TLegend *leg_max_njet_dijets = new TLegend(0.6, 0.6, 0.9, 0.9);
  double max_val_dijets_njet = 0;

  for (int B = 0; B < betalist.size(); B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;
    
    TH1* max_njet_dijets_manual = (TH1*)max_njet_dijets_manual_hists.At(B);

    double max_njet_dijets_manual_scale = 1/max_njet_dijets_manual->Integral();
    max_njet_dijets_manual->Scale(max_njet_dijets_manual_scale);

    max_njet_dijets_manual->SetStats(0);

    if (B == 0) max_njet_dijets_manual->SetLineColor(kGreen);
    if (B == 1) max_njet_dijets_manual->SetLineColor(kYellow);
    if (B == 2) max_njet_dijets_manual->SetLineColor(kBlue);
    if (B == 3) max_njet_dijets_manual->SetLineColor(kRed);

    max_njet_dijets_manual->Draw("SAMES");
    leg_max_njet_dijets->AddEntry(max_njet_dijets_manual, "beta = " + (TString)ss.str());
    max_njet_dijets_manual->Write();

    if (max_njet_dijets_manual->GetMaximum() > max_val_dijets_raw) {
      max_val_dijets_raw = max_njet_dijets_manual->GetMaximum();
      max_njet_dijets_manual->SetMaximum(max_val_dijets_raw);
    }

    cout << "manual QCD max_njet ratio: " << max_njet_dijets_manual->Integral(75, 110) << "(beta = " << beta << ")" << endl;

  }

  leg_max_njet_dijets->Draw();
  max_njet_dijets_compare->Write();
  max_njet_dijets_compare->Print("ttbarstudy_dijets_max_njet_compare_nplus1.pdf", "pdf");


  cout << endl;
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
  delete tau32_ttbar_wta;
  delete tau32_dijets_wta;
  // delete tau32_ttbar_manual;
  // delete tau32_dijets_manual;
  // delete mass_ttbar_manual;
  // delete mass_dijets_manual;
  delete mass_ttbar_antikt;
  delete mass_dijets_antikt;
  delete ROC_compare;
  // delete tau3_ttbar_manual_antikt;
  // delete tau3_dijets_manual_antikt;
  // delete mass_ttbar_compare;

// }

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