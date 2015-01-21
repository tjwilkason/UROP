//Misc. Headers
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

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
#include "fastjet/contrib/WinnerTakeAllRecombiner.hh"

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
#include "TRandom.h"

#include <queue>
#include <unordered_map>

using namespace std;
using namespace Pythia8;
using namespace fastjet;
using namespace fastjet::contrib;

struct BranchFitness {
  PseudoJet branch;
  double    fitness;
  double    zfrac;
};

// BranchFitness calcFitness(PseudoJet particle, double alpha) {
// double calcFitness(PseudoJet particle, double alpha) {
//   PseudoJet parent0, parent1;
//   double fitness;
//   // double zfrac;
  
//   fitness = particle.exclusive_subdmerge(1);
//   return fitness;
// }


vector<BranchFitness> splitBranch(BranchFitness blob, double alpha) {
  PseudoJet branch = blob.branch;
  PseudoJet parent0, parent1;
  branch.has_parents(parent0, parent1);

  BranchFitness parent0Blob, parent1Blob;

  if (parent0.has_structure() && parent1.has_structure()) {

    double angle_distance = TMath::Sqrt(TMath::Power(parent0.eta() - parent1.eta(), 2) + TMath::Power(parent0.phi() - parent1.phi(), 2));
    double parent0frac = (parent0.perp() / (parent0.perp() + parent1.perp()))*TMath::Power(angle_distance, alpha);
    double parent1frac = (parent1.perp() / (parent0.perp() + parent1.perp()))*TMath::Power(angle_distance, alpha);

    // double parent0Fitness = calcFitness(parent0, alpha);
    // double parent1Fitness = calcFitness(parent1, alpha);

    double parent0Fitness = parent0.exclusive_subdmerge(1);
    double parent1Fitness = parent1.exclusive_subdmerge(1);

    parent0Blob = {parent0, parent0Fitness, parent0frac};
    parent1Blob = {parent1, parent1Fitness, parent1frac};
  }

  vector<BranchFitness> parentBlobs;
  parentBlobs.push_back(parent0Blob);
  parentBlobs.push_back(parent1Blob);

  return parentBlobs;
}

// This will allow the priority queue to sort the branches for you.
class CompareJetBranch {
  public:
  bool operator()(BranchFitness b1, BranchFitness b2) {
    // Returns true if b1 is lower in the queue than b2
    return (b1.fitness < b2.fitness);
    // return (b1.fitness > b2.fitness);
  }
};

// NEED HELPER FUNCTION TO CREATE VECTOR FROM QUEUE

vector<PseudoJet> queue2Vector(priority_queue<BranchFitness, vector<BranchFitness>, CompareJetBranch> myQueue) {
  vector<PseudoJet> myVector;
  int size = myQueue.size();
  for (int i = 0; i < size; i++) {
    BranchFitness tempBlob = myQueue.top();
    myVector.push_back(tempBlob.branch);
    myQueue.pop();
  }
  for (int i = 0; i < size; i++) {
    BranchFitness tempBlob = {myVector[i], myVector[i].exclusive_subdmerge(1)};
    myQueue.push(tempBlob);
  }
  return myVector;
}

vector<PseudoJet> runBranchQueue(PseudoJet firstBranch, int n_axes, double alpha, double zcut, double Rmax, double mu = 1.0) {
  priority_queue<BranchFitness, vector<BranchFitness>, CompareJetBranch> branchQueue;

  // double firstBlobfitness = calcFitness(firstBranch, alpha);
  double firstBlobfitness = firstBranch.exclusive_subdmerge(1);
  BranchFitness firstBlob = {firstBranch, firstBlobfitness, 0.0};
  branchQueue.push(firstBlob);

  vector<PseudoJet> branchAxes;
  BranchFitness currentBlob;

  vector<PseudoJet> tempAxes;
  int max_axes_size = 0;

  while (branchQueue.size() < n_axes) {

    if (branchQueue.size() == 0) {
      // cerr << "Error: branch queue empty." << endl;
      break;
    }

    if (branchQueue.size() >= max_axes_size) {
      max_axes_size = branchQueue.size();
      tempAxes = queue2Vector(branchQueue);
      // branchQueue = vector2Queue(tempAxes);
    }

    currentBlob = branchQueue.top();
    branchQueue.pop();

    vector<BranchFitness> parentBlobs = splitBranch(currentBlob, alpha);
    if (parentBlobs[0].branch.has_structure() && parentBlobs[1].branch.has_structure()) {
      double parent0splitting = parentBlobs[0].branch.delta_R(firstBranch);
      double parent1splitting = parentBlobs[1].branch.delta_R(firstBranch);
      // else if (parentBlobs[0].branch.delta_R(parentBlobs[1].branch) < 0.1) {
      //   BranchFitness harderBranch = (parentBlobs[0].branch.perp() > parentBlobs[1].branch.perp()) ? parentBlobs[0] : parentBlobs[1];
      //   branchQueue.push(harderBranch);
      // }
      if (parentBlobs[0].zfrac < zcut || parent0splitting > Rmax) {
        branchQueue.push(parentBlobs[1]);
      }
      else if (parentBlobs[1].zfrac < zcut || parent1splitting > Rmax) {
        branchQueue.push(parentBlobs[0]);
      }
      // check the mass drop (newly added -- 06/02)
      // else if (mu == 1.0 || (max(parentBlobs[0].branch.m2(), parentBlobs[1].branch.m2()) > mu*mu*currentBlob.branch.m2())) {
      else {
        branchQueue.push(parentBlobs[0]);
        branchQueue.push(parentBlobs[1]);
      }
    }
  } 

  int n_found = branchQueue.size();
  if (n_found == 0) branchAxes = tempAxes;

  for (int i = 0; i < n_found; i++) {
    BranchFitness final_blob = branchQueue.top();    
    branchAxes.push_back(final_blob.branch);
    branchQueue.pop();
  }


  return branchAxes;
}

double calculateRadius(vector<PseudoJet> axes) {
// vector<Double_t> nearest_neighbors(vector<PseudoJet> axes, int N = 3, int K = 2) {
  // vector<Double_t> neighbors;

  std::string bitmask(2, 1); // K leading 1's
  bitmask.resize(axes.size(), 0); // N-K trailing 0's

  double avg_distance = 0;
  double max_distance = 0;
  int axis_counter = 0;

  do {
    vector<int> axis_indices;
    for (int i = 0; i < axes.size(); ++i) // [0..N-1] integers
    {
      if (bitmask[i]) axis_indices.push_back(i);
    }
    if (axis_indices.size() > 1) {
      double distance = axes[axis_indices[0]].delta_R(axes[axis_indices[1]]);
      avg_distance += distance;
      if (distance > max_distance) max_distance = distance;
      axis_counter++;
    }
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  avg_distance = avg_distance/(2*axis_counter); 
  // double radius = max(1.0, avg_distance);
  // double radius = avg_distance;
  // double radius = 0.6;
  double radius = max_distance/2;
  return radius;
}

vector<PseudoJet> clusterJets(vector<PseudoJet> my_particles, PseudoJet firstBranch, int n_jets, double alpha, double zcut, double Rmax, double mu = 1.0) {

  JetDefinition::Recombiner *recombScheme = new WinnerTakeAllRecombiner();

  // FIND AXES FOR JETS
  vector<PseudoJet> my_axes = runBranchQueue(firstBranch, n_jets, alpha, zcut, Rmax);
  double Rparam = calculateRadius(my_axes);

  vector<PseudoJet> final_jets;

  // INITIALIZE JETS
  for (int i = 0; i < my_axes.size(); i++) {
    PseudoJet empty_jet(0,0,0,0);
    final_jets.push_back(empty_jet);
  }

  // TAKE INITIAL PARTICLES AND CLUSTER THEM ACCORDING TO NEAREST AXIS
  for (int i_particles = 0; i_particles < my_particles.size(); i_particles++) {
    Double_t min_axis_distance = 10000;
    int min_axis_index = -1;
    for (int i_axis = 0; i_axis < my_axes.size(); i_axis++) {
      Double_t axis_distance = my_particles[i_particles].delta_R(my_axes[i_axis]);
      if ((axis_distance < min_axis_distance) && (axis_distance < Rparam)) {
        min_axis_distance = axis_distance;
        min_axis_index = i_axis;
      }
    }
    if (min_axis_index != -1) {
      final_jets[min_axis_index] = join(final_jets[min_axis_index], my_particles[i_particles], *recombScheme);
      // final_jets[min_axis_index] = join(final_jets[min_axis_index], my_particles[i_particles]);
    }
  }

  return final_jets;
  delete recombScheme;
}

vector<PseudoJet> clusterJets(vector<PseudoJet> my_particles, vector<PseudoJet> my_axes) {

  JetDefinition::Recombiner *recombScheme = new WinnerTakeAllRecombiner();

  // FIND AXES FOR JETS
  double Rparam = calculateRadius(my_axes);

  vector<PseudoJet> final_jets;

  // INITIALIZE JETS
  for (int i = 0; i < my_axes.size(); i++) {
    PseudoJet empty_jet(0,0,0,0);
    final_jets.push_back(empty_jet);
  }

  // TAKE INITIAL PARTICLES AND CLUSTER THEM ACCORDING TO NEAREST AXIS
  for (int i_particles = 0; i_particles < my_particles.size(); i_particles++) {
    Double_t min_axis_distance = 10000;
    int min_axis_index = -1;
    for (int i_axis = 0; i_axis < my_axes.size(); i_axis++) {
      Double_t axis_distance = my_particles[i_particles].delta_R(my_axes[i_axis]);
      if ((axis_distance < min_axis_distance) && (axis_distance < Rparam)) {
        min_axis_distance = axis_distance;
        min_axis_index = i_axis;
      }
    }
    if (min_axis_index != -1) {
      final_jets[min_axis_index] = join(final_jets[min_axis_index], my_particles[i_particles], *recombScheme);
      // final_jets[min_axis_index] = join(final_jets[min_axis_index], my_particles[i_particles]);
    }
  }

  return final_jets;
  delete recombScheme;
}

double calculateMass(vector<PseudoJet> my_jets, int max_perp_location = -1) {
  PseudoJet big_jet(0,0,0,0);
  for (int i_jets = 0; i_jets < my_jets.size(); i_jets++) {
    if (i_jets != max_perp_location) big_jet = join(big_jet, my_jets[i_jets]);
  }
  return big_jet.m();
}

PseudoJet sumJets(vector<PseudoJet> my_jets) {
  PseudoJet big_jet(0,0,0,0);
  for (int i_jets = 0; i_jets < my_jets.size(); i_jets++) {
    big_jet = join(big_jet, my_jets[i_jets]);
  }
  return big_jet;
}

double findMaximumMass(vector<PseudoJet> my_particles, PseudoJet firstBranch, double alpha, double zcut, double Rmax, double mu = 1.0) {
  double max_mass = 0;
  int npronginess;
  for (int i_prong = 2; i_prong < 11; i_prong++) {
    vector<PseudoJet> nprong_jets = clusterJets(my_particles, firstBranch, i_prong, alpha, zcut, Rmax);
    double nprong_mass = calculateMass(nprong_jets);
    if (nprong_mass > max_mass) {
      max_mass = nprong_mass;
      npronginess = i_prong;
    }
  }
  return max_mass;
}

int findMaxPerp(vector<PseudoJet> myJets) {
  double max_perp = 0;
  int max_perp_index = 0;
  for (int i = 0; i < myJets.size(); i++) {
    if (myJets[i].perp() > max_perp) {
      max_perp = myJets[i].perp();
      max_perp_index = i;
    }
  }
  return max_perp_index;
}

double findWMass(vector<PseudoJet> my_particles, PseudoJet firstBranch, double alpha, double zcut, double Rmax) {
  double max_mass = 0;
  int npronginess;
  for (int i_prong = 2; i_prong < 11; i_prong++) {
    vector<PseudoJet> nprong_jets = clusterJets(my_particles, firstBranch, i_prong, alpha, zcut, Rmax);
    int max_perp_location = findMaxPerp(nprong_jets);
    double nprong_mass = calculateMass(nprong_jets, max_perp_location);
    if (nprong_mass > max_mass) {
      max_mass = nprong_mass;
      npronginess = i_prong;
    }
  }
  return max_mass;
}

int main(int argc, char* argv[]){

  TFile out("nprong_invariant_mass.root", "RECREATE");

  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW"); // Input File Name n
  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW"); // Input File Name

  std::string samples[2] = {"/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW", "/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW"};

  TH1F *invariant_mass_ttbar = new TH1F("invariant_mass_ttbar", "Invariant Mass for ttbar", 300, 0, 300);
  TH1F *invariant_mass_dijets = new TH1F("invariant_mass_dijets", "Invariant Mass for dijets", 300, 0, 300);

  TH1F *invariant_mass_ttbar_test = new TH1F("invariant_mass_ttbar_test", "Invariant Mass for ttbar", 300, 0, 300);
  TH1F *invariant_mass_dijets_test = new TH1F("invariant_mass_dijets_test", "Invariant Mass for dijets", 300, 0, 300);

  // TH1F *newvar_ttbar = new TH1F("newvar_ttbar", "New Variable for ttbar", 300, 0, 300);
  // TH1F *newvar_dijets = new TH1F("newvar_dijets", "New for Variable dijets", 300, 0, 300);

  TH1F *newvar_ttbar = new TH1F("newvar_ttbar", "New Variable for ttbar", 100, 0, 1);
  TH1F *newvar_dijets = new TH1F("newvar_dijets", "New for Variable dijets", 100, 0, 1);

  TH1F *invariant_Wmass_ttbar = new TH1F("invariant_Wmass_ttbar", "Invariant WMass for ttbar", 300, 0, 300);
  TH1F *invariant_Wmass_dijets = new TH1F("invariant_Wmass_dijets", "Invariant WMass for dijets", 300, 0, 300);

  TH1F *mass_ratio_ttbar_nocut = new TH1F("mass_ratio_ttbar_nocut", "Wmass/Tmass ratio ttbar", 100, 0, 1);
  TH1F *mass_ratio_dijets_nocut = new TH1F("mass_ratio_dijets_nocut", "Wmass/Tmass ratio dijets", 100, 0, 1);

  TH1F *mass_ratio_ttbar = new TH1F("mass_ratio_ttbar", "Wmass/Tmass ratio ttbar", 100, 0, 1);
  TH1F *mass_ratio_dijets = new TH1F("mass_ratio_dijets", "Wmass/Tmass ratio dijets", 100, 0, 1);

  TH1F *npronginess_ttbar = new TH1F("npronginess_ttbar", "NPronginess ttbar", 10, 0, 10);
  TH1F *npronginess_dijets = new TH1F("npronginess_dijets", "NPronginess dijets", 10, 0, 10);

  TH1F *std_invariant_mass_ttbar = new TH1F("std_invariant_mass_ttbar", "Invariant Mass for ttbar", 300, 0, 300);
  TH1F *std_invariant_mass_dijets = new TH1F("std_invariant_mass_dijets", "Invariant Mass for dijets", 300, 0, 300);

  TH1F *tau32_ttbar = new TH1F("tau32_ttbar", "tau32 for ttbar", 50, 0, 1);
  TH1F *tau32_dijets = new TH1F("tau32_dijets", "tau32 for dijets", 50, 0, 1);

  TH1F *tau32_ttbar_test = new TH1F("tau32_ttbar_test", "tau32 for ttbar", 50, 0, 1);
  TH1F *tau32_dijets_test = new TH1F("tau32_dijets_test", "tau32 for dijets", 50, 0, 1);

  TH1F *min_tau32_ttbar = new TH1F("min_tau32_ttbar", "tau32 for ttbar", 50, 0, 25);
  TH1F *min_tau32_dijets = new TH1F("min_tau32_dijets", "tau32 for dijets", 50, 0, 25);

  TH1F *std_tau32_ttbar = new TH1F("std_tau32_ttbar", "tau32 for ttbar", 50, 0, 1);
  TH1F *std_tau32_dijets = new TH1F("std_tau32_dijets", "tau32 for dijets", 50, 0, 1);

  double alpha = 0.0;
  // double alpha = 0.0;
  double zcut = 0.1;
  // double Rmax = JetDefinition::max_allowable_R;
  double Rmax = 1.2; 
  double mu = 1.0;

  int total_ttbar_jets = 0, total_dijet_jets = 0;
  int total_std_ttbar_jets = 0, total_std_dijet_jets = 0;

  int total_masscut_ttbar_jets = 0, total_masscut_dijets_jets = 0;

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

      TH2F *event_display;
      if (i_sample == 0) event_display = new TH2F("event_display", "ttbar event (500 < pT < 600)", 60, -3.2, 3.2, 60, 0, 6.4);
      if (i_sample == 1) event_display = new TH2F("event_display", "dijet event (500 < pT < 600)", 60, -3.2, 3.2, 60, 0, 6.4);
      TH2F *axes_display = new TH2F("axes_display", "Axes", 300, -3.2, 3.2, 300, 0, 6.4);
      TH2F *ghost_display = new TH2F("ghost_display", "Axes", 300, -3.2, 3.2, 300, 0, 6.4);
      TH2F *subjet1_display = new TH2F("subjet1_display", "Axes", 300, -3.2, 3.2, 300, 0, 6.4);
      TH2F *subjet2_display = new TH2F("subjet2_display", "Axes", 300, -3.2, 3.2, 300, 0, 6.4);
      TH2F *subjet3_display = new TH2F("subjet3_display", "Axes", 300, -3.2, 3.2, 300, 0, 6.4);
      TH2F *subjet4_display = new TH2F("subjet4_display", "Axes", 300, -3.2, 3.2, 300, 0, 6.4);
      TH2F *subjet5_display = new TH2F("subjet5_display", "Axes", 300, -3.2, 3.2, 300, 0, 6.4);
      TH2F *subjet6_display = new TH2F("subjet6_display", "Axes", 300, -3.2, 3.2, 300, 0, 6.4);

      // this should cluster all the particles into the event into a single "jet" with pT of the total sum

      for (int i = 0; i < input_particles.size(); i++) {
        event_display->Fill(input_particles[i].eta(), input_particles[i].phi(), input_particles[i].perp());
      }

      double Rparam = JetDefinition::max_allowable_R;
      Strategy strategy = Best;
      const JetDefinition::Recombiner *recombScheme = new WinnerTakeAllRecombiner();
      // RecombinationScheme recombScheme = E_scheme;
      JetDefinition *jetDef = NULL;
      jetDef = new JetDefinition(genkt_algorithm, Rparam, 0.2, recombScheme, strategy);
      // jetDef = new JetDefinition(genkt_algorithm, Rparam, 1.0);
      // jetDef->set_recombiner(recombScheme);
      ClusterSequence clustSeq(input_particles, *jetDef);
      vector <PseudoJet> inclusiveJets = clustSeq.inclusive_jets(0);
      PseudoJet total_jet = inclusiveJets[0];

      vector<PseudoJet> seedAxes = runBranchQueue(total_jet, 2, 0, 0.1, JetDefinition::max_allowable_R);

      for (int i_seed = 0; i_seed < seedAxes.size(); i_seed++) {

        double maxmass_test = 0;

        for (int i_nsub = 2; i_nsub < 4; i_nsub++) {
          double min_mass_diff = 100000;
          double min_tau3 = 100000;
          PseudoJet min_test_bigjet;
          double min_test_zcut = 0;
          vector<PseudoJet> min_test_axes;
          double max_ptdr = 0;
          double max_mass_nsub = 0;

          for (int i_zcut = 0; i_zcut < 10; i_zcut++) {
            double zcut_test = 0.01*i_zcut + 0.05;
            vector<PseudoJet> test_axes = runBranchQueue(seedAxes[i_seed], i_nsub, alpha, zcut_test, Rmax);
            vector<PseudoJet> test_jets = clusterJets(input_particles, test_axes);
            PseudoJet test_bigjet = sumJets(test_jets);
            Nsubjettiness tau3_test(i_nsub, Njettiness::manual_axes, Njettiness::unnormalized_measure, 1.0);
            if (test_axes.size() == i_nsub) {
              tau3_test.setAxes(test_axes);
              double test_ptdr = 0;
              for (int i = 0; i < test_jets.size(); i++) {
                test_ptdr += TMath::Power(test_jets[i].perp(),1.0)*TMath::Power(test_jets[i].delta_R(seedAxes[i_seed]),1.0);
              }
              double tau3_test_result = tau3_test.result(test_bigjet);

              // if (tau3_test_result < min_tau3) {
              // if (abs(test_bigjet.m() - 175) < min_mass_diff) {
              if (test_ptdr > max_ptdr) {
              // if (test_bigjet.m() > max_mass_nsub) {
                min_tau3 = tau3_test_result;
                min_mass_diff = abs(test_bigjet.m() - 175);
                min_test_bigjet = test_bigjet;
                min_test_zcut = zcut_test;
                max_ptdr = test_ptdr;
                // max_mass_nsub = test_bigjet.m();

                if (i_nsub == 3) min_test_axes = test_axes;
              }
            }
          }
          // cout << min_test_zcut << endl;
          if (i_nsub == 3) {
            for (int i = 0; i < min_test_axes.size(); i++) {
              axes_display->Fill(min_test_axes[i].eta(), min_test_axes[i].phi());
            }
          }
          if (min_test_bigjet.m() > maxmass_test) maxmass_test = min_test_bigjet.m();
        }
        


        vector<PseudoJet> jets = clusterJets(input_particles, seedAxes[i_seed], 3, alpha, zcut, Rmax, mu);
        PseudoJet big_jet = sumJets(jets);

        double maxmass = findMaximumMass(input_particles, seedAxes[i_seed], alpha, zcut, Rmax, mu);
        // double maxWmass = findWMass(input_particles, seedAxes[i_seed], alpha, zcut, Rmax);

        if (big_jet.perp() > 100) {
          if (i_sample == 0) {
            invariant_mass_ttbar->Fill(maxmass);
            invariant_mass_ttbar_test->Fill(maxmass_test);
            // newvar_ttbar->Fill(test_var);
            // invariant_Wmass_ttbar->Fill(maxWmass);
            // mass_ratio_ttbar_nocut->Fill(maxWmass/maxmass);
            total_ttbar_jets++;
          }
          if (i_sample == 1) {
            invariant_mass_dijets->Fill(maxmass);
            invariant_mass_dijets_test->Fill(maxmass_test);
            // newvar_dijets->Fill(test_var);
            // invariant_Wmass_dijets->Fill(maxWmass);
            // mass_ratio_dijets_nocut->Fill(maxWmass/maxmass);
            total_dijet_jets++;
          }
        }

        Nsubjettiness tau3(3, Njettiness::manual_axes, Njettiness::unnormalized_measure, 1.0);
        Nsubjettiness tau2(2, Njettiness::manual_axes, Njettiness::unnormalized_measure, 1.0);

        vector<PseudoJet> axes_3 = runBranchQueue(seedAxes[i_seed], 3, alpha, zcut, Rmax, mu);
        vector<PseudoJet> axes_2 = runBranchQueue(seedAxes[i_seed], 2, alpha, zcut, Rmax, mu);

        // for (int i = 0; i < axes_3.size(); i++) {
        //   axes_display->Fill(axes_3[i].eta(), axes_3[i].phi());
        // }
        
        if (axes_3.size() == 3 && axes_2.size() == 2) {
          tau3.setAxes(axes_3);
          tau2.setAxes(axes_2);
        
          if (maxmass > 160 && maxmass < 240 && big_jet.perp() > 100) {
            // if (tau3.result(big_jet)/tau2.result(big_jet) > 0.8) {
            // }
            if (i_sample == 0) {
              tau32_ttbar->Fill(tau3.result(big_jet)/tau2.result(big_jet));
              // tau32_ttbar_test->Fill((double)min_test_tauN[1]/min_test_tauN[0]);

              // mass_ratio_ttbar->Fill(maxWmass/maxmass);
            }
            if (i_sample == 1) {
              tau32_dijets->Fill(tau3.result(big_jet)/tau2.result(big_jet));
              // tau32_dijets_test->Fill((double)min_test_tauN[1]/min_test_tauN[0]);

              // mass_ratio_dijets->Fill(maxWmass/maxmass);
            }
          }
        }

        // double max_mass = 0;
        // int npronginess;
        // for (int i_prong = 2; i_prong < 11; i_prong++) {
        //   vector<PseudoJet> nprong_jets = clusterJets(input_particles, seedAxes[i_seed], i_prong, alpha, zcut, Rmax);
        //   double nprong_mass = calculateMass(nprong_jets);
        //   nprong_masses[i_prong - 2] = nprong_mass;
        //   if (nprong_mass > max_mass) {
        //     max_mass = nprong_mass;
        //     npronginess = i_prong;
        //   }
        // }

      }

      /******************** CALCULATE NORMAL INVARIANT MASS THROUGH STANDARD PROCEDURE HERE ***************************/
      double Rparam_2 = 1.0;
      Strategy strategy_2 = Best;
      RecombinationScheme recombScheme_2 = E_scheme;
      JetDefinition *jetDef_2 = NULL;
      jetDef_2 = new JetDefinition(antikt_algorithm, Rparam_2, recombScheme_2, strategy_2);

      vector <PseudoJet> inclusiveJets_quark, sortedJets_quark, centralJets_quark;
      ClusterSequence clustSeq_quark(input_particles, *jetDef_2);

      inclusiveJets_quark = clustSeq_quark.inclusive_jets(100.0);
      sortedJets_quark = sorted_by_pt(inclusiveJets_quark);

      Selector eta_selector = SelectorAbsEtaMax(2.0);
      centralJets_quark = eta_selector(sortedJets_quark);

      Selector jet_selector = SelectorNHardest(2);
      vector<PseudoJet> hardest_jet_quark = jet_selector(centralJets_quark);

      for (int i = 0; i < hardest_jet_quark.size(); i++) {
        if (i_sample == 0) {
          std_invariant_mass_ttbar->Fill(hardest_jet_quark[i].m());
          total_std_ttbar_jets++;
        }
        if (i_sample == 1) {
          std_invariant_mass_dijets->Fill(hardest_jet_quark[i].m());
          total_std_dijet_jets++;
        }

        Nsubjettiness std_tau3(3, Njettiness::wta_kt_axes, Njettiness::unnormalized_measure, 1.0);
        Nsubjettiness std_tau2(2, Njettiness::wta_kt_axes, Njettiness::unnormalized_measure, 1.0);

        if (hardest_jet_quark[i].m() > 160 && hardest_jet_quark[i].m() < 240) {
          if (i_sample == 0) std_tau32_ttbar->Fill(std_tau3.result(hardest_jet_quark[i])/std_tau2.result(hardest_jet_quark[i]));
          if (i_sample == 1) std_tau32_dijets->Fill(std_tau3.result(hardest_jet_quark[i])/std_tau2.result(hardest_jet_quark[i]));
        }

      }

      // // FILL EVENT WITH GHOSTS TO MAKE FIREWORKS
      // Double_t ghost_perp = 0.001;
      // Double_t n_ghosts = 200;
      // for (int i_ghosts = 0; i_ghosts < n_ghosts; i_ghosts++) {
      //   for (int j_ghosts = 0; j_ghosts < n_ghosts; j_ghosts++) {
      //     Double_t ghost_eta = -5 + (double)(10/n_ghosts)*i_ghosts;
      //     Double_t ghost_phi = (double)(2*TMath::Pi()/n_ghosts)*j_ghosts;
      //     PseudoJet ghost(ghost_perp*TMath::Cos(ghost_phi),ghost_perp*TMath::Sin(ghost_phi),
      //        ghost_perp*TMath::SinH(ghost_eta),ghost_perp*TMath::CosH(ghost_eta));
      //     Double_t min_axis_distance = 10000;
      //     int min_axis_index = -1;
      //     for (int i_axis = 0; i_axis < axes.size(); i_axis++) {
      //       Double_t Rparam_subjet = axis_neighbors[i_axis];
      //       Double_t axis_distance = ghost.delta_R(axes[i_axis]);
      //       if (axis_distance < min_axis_distance && axis_distance < Rparam_subjet) {
      //         min_axis_distance = axis_distance;
      //         min_axis_index = i_axis;
      //       }
      //     }
      //     // if (min_axis_index != -1) final_jets[min_axis_index] = join(final_jets[min_axis_index], ghost, *recombScheme);
      //     if (min_axis_index == 0) subjet1_display->Fill(ghost_eta, ghost_phi);
      //     if (min_axis_index == 1) subjet2_display->Fill(ghost_eta, ghost_phi);
      //     if (min_axis_index == 2) subjet3_display->Fill(ghost_eta, ghost_phi);
      //     if (min_axis_index == 3) subjet4_display->Fill(ghost_eta, ghost_phi);
      //     if (min_axis_index == 4) subjet5_display->Fill(ghost_eta, ghost_phi);
      //     if (min_axis_index == 5) subjet6_display->Fill(ghost_eta, ghost_phi);
      //   }
      // }

      event_display->SetStats(0);
      axes_display->SetStats(0);
      ghost_display->SetStats(0);
      subjet1_display->SetStats(0);
      subjet2_display->SetStats(0);
      subjet3_display->SetStats(0);
      subjet4_display->SetStats(0);
      subjet5_display->SetStats(0);
      subjet6_display->SetStats(0);
      TCanvas *display = new TCanvas("display", "Event Display", 700, 700);
      display->cd();
      event_display->GetXaxis()->SetTitle("#eta");
      event_display->GetYaxis()->SetTitle("#phi");
      event_display->SetFillColor(kBlack);
      axes_display->SetMarkerStyle(3);
      axes_display->SetMarkerSize(3);
      axes_display->SetMarkerColor(kRed);
      ghost_display->SetFillColor(kBlack);
      subjet1_display->SetMarkerColor(kRed);
      subjet2_display->SetMarkerColor(kBlue);
      subjet3_display->SetMarkerColor(kGreen);
      subjet4_display->SetMarkerColor(kYellow);
      subjet5_display->SetMarkerColor(kOrange);
      subjet6_display->SetMarkerColor(6);
      subjet1_display->SetMarkerStyle(21);
      subjet2_display->SetMarkerStyle(21);
      subjet3_display->SetMarkerStyle(21);
      subjet4_display->SetMarkerStyle(21);
      subjet5_display->SetMarkerStyle(21);
      subjet6_display->SetMarkerStyle(21);
      subjet1_display->SetMarkerSize(0.2);
      subjet2_display->SetMarkerSize(0.2);
      subjet3_display->SetMarkerSize(0.2);
      subjet4_display->SetMarkerSize(0.2);
      subjet5_display->SetMarkerSize(0.2);
      subjet6_display->SetMarkerSize(0.2);
      // // axes_display->SetFillColor(kRed);
      event_display->Draw("box");
      axes_display->Draw("SAMES");

      // ghost_display->Draw("boxSAMES");

      // subjet1_display->Draw("SAMES");
      // subjet2_display->Draw("SAMES");
      // subjet3_display->Draw("SAMES");
      // subjet4_display->Draw("SAMES");
      // subjet5_display->Draw("SAMES");
      // subjet6_display->Draw("SAMES");

      // if (axes_display->Integral() != 0) display->Write();

      // display->Write();

      // ostringstream ss;
      // ss << i_event;
      // TString title;
      // if (i_sample == 0) title = "ttbar_display" + ss.str() + ".pdf";
      // if (i_sample == 1) title = "dijets_display" + ss.str() + ".pdf";

      // display->Print(title, "pdf");

      delete display;
      delete event_display;
      delete axes_display;
      delete ghost_display;
      delete subjet1_display;
      delete subjet2_display;
      delete subjet3_display;
      delete subjet4_display;
      delete subjet5_display;
      delete subjet6_display;

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
  // out.Close();
  }

  npronginess_ttbar->SetStats(0);
  npronginess_dijets->SetStats(0);
  npronginess_ttbar->Write();
  npronginess_dijets->Write();

  newvar_ttbar->SetStats(0);
  newvar_dijets->SetStats(0);
  newvar_ttbar->Write();
  newvar_dijets->Write();

  invariant_mass_ttbar_test->SetStats(0);
  invariant_mass_ttbar_test->Write();
  invariant_mass_dijets_test->SetStats(0);
  invariant_mass_dijets_test->Write();
  tau32_ttbar_test->SetStats(0);
  tau32_ttbar_test->Write();
  tau32_dijets_test->SetStats(0);
  tau32_dijets_test->Write();

  invariant_mass_ttbar->SetStats(0);
  invariant_mass_dijets->SetStats(0);
  invariant_Wmass_ttbar->SetStats(0);
  invariant_Wmass_dijets->SetStats(0);
  std_invariant_mass_ttbar->SetStats(0);
  std_invariant_mass_dijets->SetStats(0);
  mass_ratio_ttbar->SetStats(0);
  mass_ratio_dijets->SetStats(0);
  mass_ratio_ttbar_nocut->SetStats(0);
  mass_ratio_dijets_nocut->SetStats(0);
  tau32_ttbar->SetStats(0);
  tau32_dijets->SetStats(0);
  min_tau32_ttbar->SetStats(0);
  min_tau32_dijets->SetStats(0);

  invariant_mass_ttbar->Write();
  invariant_mass_dijets->Write();
  invariant_Wmass_ttbar->Write();
  invariant_Wmass_dijets->Write();
  std_invariant_mass_ttbar->Write();
  std_invariant_mass_dijets->Write();
  mass_ratio_ttbar->Write();
  mass_ratio_dijets->Write();
  mass_ratio_ttbar_nocut->Write();
  mass_ratio_dijets_nocut->Write();

  tau32_ttbar->Write();
  tau32_dijets->Write();
  std_tau32_ttbar->Write();
  std_tau32_dijets->Write();
  min_tau32_ttbar->Write();
  min_tau32_dijets->Write();

  Double_t invariant_mass_ttbar_scale = 1/invariant_mass_ttbar->Integral();
  Double_t invariant_mass_dijets_scale = 1/invariant_mass_dijets->Integral();
  Double_t std_invariant_mass_ttbar_scale = 1/std_invariant_mass_ttbar->Integral();
  Double_t std_invariant_mass_dijets_scale = 1/std_invariant_mass_dijets->Integral();

  Double_t invariant_mass_ttbar_test_scale = 1/invariant_mass_ttbar_test->Integral();
  Double_t invariant_mass_dijets_test_scale = 1/invariant_mass_dijets_test->Integral();

  // Double_t tau32_ttbar_scale = 1/tau32_ttbar->Integral();
  // Double_t tau32_dijets_scale = 1/tau32_dijets->Integral();
  // Double_t std_tau32_ttbar_scale = 1/std_tau32_ttbar->Integral();
  // Double_t std_tau32_dijets_scale = 1/std_tau32_dijets->Integral();

  invariant_mass_ttbar->Scale(invariant_mass_ttbar_scale);
  invariant_mass_dijets->Scale(invariant_mass_dijets_scale);
  std_invariant_mass_ttbar->Scale(std_invariant_mass_ttbar_scale);
  std_invariant_mass_dijets->Scale(std_invariant_mass_dijets_scale);

  invariant_mass_ttbar_test->Scale(invariant_mass_ttbar_test_scale);
  invariant_mass_dijets_test->Scale(invariant_mass_dijets_test_scale);

  // tau32_ttbar->Scale(tau32_ttbar_scale);
  // tau32_dijets->Scale(tau32_dijets_scale);
  // std_tau32_ttbar->Scale(std_tau32_ttbar_scale);
  // std_tau32_dijets->Scale(std_tau32_dijets_scale);

  TCanvas *mass_compare = new TCanvas("mass_compare", "Mass Comparison", 600, 600);
  mass_compare->cd();
  invariant_mass_ttbar->SetLineColor(kRed);
  std_invariant_mass_ttbar->SetLineColor(kBlue);
  invariant_mass_ttbar->Draw();
  std_invariant_mass_ttbar->Draw("same");
  TLegend *leg_mass = new TLegend(0.1, 0.7, 0.3, 0.9);
  leg_mass->AddEntry(std_invariant_mass_ttbar, "WTA kT", "L");
  leg_mass->AddEntry(invariant_mass_ttbar, "Declustering", "L");
  leg_mass->Draw();

  mass_compare->Write();

  TCanvas *mass_compare_2 = new TCanvas("mass_compare_2", "Mass Comparison", 600, 600);
  mass_compare_2->cd();
  invariant_mass_dijets->SetLineColor(kRed);
  std_invariant_mass_dijets->SetLineColor(kBlue);
  invariant_mass_dijets->Draw();
  std_invariant_mass_dijets->Draw("same");
  TLegend *leg_mass_2 = new TLegend(0.1, 0.7, 0.3, 0.9);
  leg_mass_2->AddEntry(std_invariant_mass_dijets, "WTA kT", "L");
  leg_mass_2->AddEntry(invariant_mass_dijets, "Declustering", "L");
  leg_mass_2->Draw();

  mass_compare_2->Write();

  TCanvas *tau32_compare = new TCanvas("tau32_compare", "Mass Comparison", 600, 600);
  tau32_compare->cd();
  tau32_ttbar->SetLineColor(kRed);
  tau32_dijets->SetLineColor(kBlue);
  tau32_ttbar->Draw();
  tau32_dijets->Draw("same");
  TLegend *leg_tau32 = new TLegend(0.1, 0.7, 0.3, 0.9);
  leg_tau32->AddEntry(tau32_ttbar, "TTbar", "L");
  leg_tau32->AddEntry(tau32_dijets, "QCD", "L");
  leg_tau32->Draw();

  tau32_compare->Write();

  double ttbar_efficiencies[4];
  double dijets_efficiencies[4];

  ttbar_efficiencies[0] = (double)invariant_mass_ttbar->Integral(120, 240);
  dijets_efficiencies[0] = (double)invariant_mass_dijets->Integral(120, 240);

  ttbar_efficiencies[1] = (double)invariant_mass_ttbar->Integral(160, 240);
  dijets_efficiencies[1] = (double)invariant_mass_dijets->Integral(160, 240);

  ttbar_efficiencies[2] = (double)invariant_mass_ttbar->Integral(160, 200);
  dijets_efficiencies[2] = (double)invariant_mass_dijets->Integral(160, 200);

  ttbar_efficiencies[3] = (double)invariant_mass_ttbar->Integral(160, 180);
  dijets_efficiencies[3] = (double)invariant_mass_dijets->Integral(160, 180);

  TGraph *mass_efficiencies = new TGraph(4, ttbar_efficiencies, dijets_efficiencies);
  mass_efficiencies->SetMarkerColor(kBlue);
  mass_efficiencies->SetLineColor(kBlue);
  mass_efficiencies->SetLineWidth(2);
  mass_efficiencies->SetMarkerStyle(3);
  mass_efficiencies->Write();

  double ttbar_efficiencies_test[4];
  double dijets_efficiencies_test[4];

  ttbar_efficiencies_test[0] = (double)invariant_mass_ttbar_test->Integral(120, 240);
  dijets_efficiencies_test[0] = (double)invariant_mass_dijets_test->Integral(120, 240);

  ttbar_efficiencies_test[1] = (double)invariant_mass_ttbar_test->Integral(160, 240);
  dijets_efficiencies_test[1] = (double)invariant_mass_dijets_test->Integral(160, 240);

  ttbar_efficiencies_test[2] = (double)invariant_mass_ttbar_test->Integral(160, 200);
  dijets_efficiencies_test[2] = (double)invariant_mass_dijets_test->Integral(160, 200);

  ttbar_efficiencies_test[3] = (double)invariant_mass_ttbar_test->Integral(160, 180);
  dijets_efficiencies_test[3] = (double)invariant_mass_dijets_test->Integral(160, 180);

  TGraph *mass_efficiencies_test = new TGraph(4, ttbar_efficiencies_test, dijets_efficiencies_test);
  mass_efficiencies_test->SetMarkerColor(kGreen);
  mass_efficiencies_test->SetLineColor(kGreen);
  mass_efficiencies_test->SetLineWidth(2);
  mass_efficiencies_test->SetMarkerStyle(3);
  mass_efficiencies_test->Write();

  int n_wtakt = std_tau32_ttbar->GetSize() - 2;
  double integral_dijets_wtakt[n_wtakt], integral_ttbar_wtakt[n_wtakt];
  for (int i = 0; i < n_wtakt; i++) {
    // integral_ttbar_wtakt[i] = std_tau32_ttbar->Integral(0,i)/total_std_ttbar_jets;
    // integral_dijets_wtakt[i] = std_tau32_dijets->Integral(0,i)/total_std_dijet_jets;
    integral_ttbar_wtakt[i] = std_tau32_ttbar->Integral(0,i)*std_invariant_mass_ttbar_scale;
    integral_dijets_wtakt[i] = std_tau32_dijets->Integral(0,i)*std_invariant_mass_dijets_scale;
  }
  TGraph* ROC_tau32_std = new TGraph(n_wtakt, integral_ttbar_wtakt, integral_dijets_wtakt);
  ROC_tau32_std->GetXaxis()->SetLimits(0, 1);
  ROC_tau32_std->GetYaxis()->SetLimits(0, 1);
  ROC_tau32_std->SetLineColor(kBlack);
  ROC_tau32_std->SetLineWidth(2);
  ROC_tau32_std->SetMarkerStyle(5);
  ROC_tau32_std->SetMarkerSize(2);
  ROC_tau32_std->Write();

  int n_size = tau32_ttbar->GetSize() - 2;
  double integral_dijets[n_size], integral_ttbar[n_size];
  for (int i = 0; i < n_size; i++) {
    integral_ttbar[i] = tau32_ttbar->Integral(0,i)/total_ttbar_jets;
    integral_dijets[i] = tau32_dijets->Integral(0,i)/total_dijet_jets;
  }
  TGraph* ROC_tau32 = new TGraph(n_size, integral_ttbar, integral_dijets);
  ROC_tau32->GetXaxis()->SetLimits(0, 1);
  ROC_tau32->GetYaxis()->SetLimits(0, 1);
  ROC_tau32->SetLineColor(kRed);
  ROC_tau32->SetLineWidth(2);
  ROC_tau32->SetMarkerStyle(5);
  ROC_tau32->SetMarkerSize(2);
  ROC_tau32->Write();

  int n_size_newvar = newvar_ttbar->GetSize() - 2;
  double integral_dijets_newvar[n_size_newvar], integral_ttbar_newvar[n_size_newvar];
  for (int i = 0; i < n_size_newvar; i++) {
    integral_ttbar_newvar[i] = newvar_ttbar->Integral(n_size_newvar - i,n_size_newvar)/total_ttbar_jets;
    integral_dijets_newvar[i] = newvar_dijets->Integral(n_size_newvar - i,n_size_newvar)/total_dijet_jets;
  }
  TGraph* ROC_newvar = new TGraph(n_size_newvar, integral_ttbar_newvar, integral_dijets_newvar);
  ROC_newvar->GetXaxis()->SetLimits(0, 1);
  ROC_newvar->GetYaxis()->SetLimits(0, 1);
  ROC_newvar->SetLineColor(kYellow);
  ROC_newvar->SetLineWidth(2);
  ROC_newvar->SetMarkerStyle(5);
  ROC_newvar->SetMarkerSize(2);
  ROC_newvar->Write();

  int n_ratio = mass_ratio_ttbar->GetSize() - 2;
  double integral_dijets_ratio[n_ratio], integral_ttbar_ratio[n_ratio];
  for (int i = 0; i < n_ratio; i++) {
    integral_ttbar_ratio[i] = mass_ratio_ttbar->Integral(n_ratio - i, n_ratio)/total_ttbar_jets;
    integral_dijets_ratio[i] = mass_ratio_dijets->Integral(n_ratio - i, n_ratio)/total_dijet_jets;
  }
  TGraph* ROC_ratio = new TGraph(n_ratio, integral_ttbar_ratio, integral_dijets_ratio);
  ROC_ratio->GetXaxis()->SetLimits(0, 1);
  ROC_ratio->GetYaxis()->SetLimits(0, 1);
  ROC_ratio->SetLineColor(kGreen);
  ROC_ratio->SetLineWidth(2);
  ROC_ratio->SetMarkerStyle(5);
  ROC_ratio->SetMarkerSize(2);
  ROC_ratio->Write();

  TCanvas *ROC_compare = new TCanvas("ROC_compare", "ROC Curve Comparison", 600, 600);
  ROC_compare->cd();
  ROC_compare->SetLogy();
  TMultiGraph *ROC_multigraph = new TMultiGraph("ROC_multigraph", "ROC Comparison for #tau_{3}/#tau_{2} (500 < p_{T} < 600)");
  ROC_multigraph->Add(ROC_tau32_std);
  ROC_multigraph->Add(ROC_tau32);
  ROC_multigraph->Add(mass_efficiencies);
  ROC_multigraph->Add(mass_efficiencies_test);
  ROC_multigraph->Add(ROC_ratio);
  ROC_multigraph->Add(ROC_newvar);
  ROC_multigraph->Draw("AL");
  TLegend *leg_ROC = new TLegend(0.1, 0.7, 0.3, 0.9);
  leg_ROC->AddEntry(ROC_tau32_std, "WTA kT", "L");
  leg_ROC->AddEntry(ROC_tau32, "Declustering", "L");
  leg_ROC->AddEntry(mass_efficiencies, "Max Mass Alone", "L");
  leg_ROC->AddEntry(mass_efficiencies_test, "Max Mass Test", "L");
  // leg_ROC->AddEntry(ROC_newvar, "New Variable", "L");
  // leg_ROC->AddEntry(ROC_ratio, "Mass Ratio", "L");
  leg_ROC->Draw();
  ROC_compare->Write();

}