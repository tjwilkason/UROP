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
double calcFitness(PseudoJet particle, double alpha) {
  PseudoJet parent0, parent1;
  double fitness;
  // double zfrac;
  
  fitness = particle.exclusive_subdmerge(1);

  // bool has_parents = particle.has_parents(parent0, parent1);
  // if (particle.has_parents(parent0, parent1)) {

    // if (parent0.has_structure() && parent1.has_structure()) {
      // double zfrac_temp = min(parent0.perp(), parent1.perp())/(parent0.perp() + parent1.perp());
      // double zfrac_temp = 4*parent0.perp()*parent1.perp()/(TMath::Power(parent0.perp() + parent1.perp(), 2));
      // zfrac = parent0.perp()*parent1.perp()/(parent0.perp() + parent1.perp());
      // zfrac = TMath::Power((parent0.perp() + parent1.perp()),2)/(parent0.perp()*parent1.perp());
      // zfrac = parent0.perp()*parent1.perp();
      // zfrac = min(parent0.perp(),parent1.perp());
      // zfrac = max(parent0.perp(),parent1.perp());
      // zfrac = parent0.perp() + parent1.perp();
      // zfrac = max(TMath::Power(parent0.perp(),2)/parent0.E(),TMath::Power(parent1.perp(),2)/parent1.E()); //WORKS FAIRLY DECENTLY
      // zfrac = (TMath::Power(parent0.perp(),2) + TMath::Power(parent0.perp(),2))/max(parent0.E(), parent1.E()); 
      // zfrac = TMath::Power(parent0.perp(),2) + TMath::Power(parent0.perp(),2); 
      // double angle_distance = TMath::Sqrt(TMath::Power(parent0.eta() - parent1.eta(), 2) + TMath::Power(parent0.phi() - parent1.phi(), 2));

      // zfrac = zfrac_temp;
      // zfrac = zfrac_temp*TMath::Power(angle_distance, alpha);
      // fitness = zfrac*TMath::Power(angle_distance, 0.5);
      // fitness = zfrac;
      // fitness = particle.perp();
    // }
    // else {
    //   fitness = 0;
    //   zfrac = 0;
    // }

  // }
  
  // BranchFitness branchfit = {particle, fitness, zfrac};
  // return branchfit;
  return fitness;
}

vector<BranchFitness> splitBranch(BranchFitness blob, double alpha) {
  PseudoJet branch = blob.branch;
  PseudoJet parent0, parent1;
  branch.has_parents(parent0, parent1);

  BranchFitness parent0Blob, parent1Blob;

  if (parent0.has_structure() && parent1.has_structure()) {

    double angle_distance = TMath::Sqrt(TMath::Power(parent0.eta() - parent1.eta(), 2) + TMath::Power(parent0.phi() - parent1.phi(), 2));
    double parent0frac = (parent0.perp() / (parent0.perp() + parent1.perp()))*TMath::Power(angle_distance, alpha);
    double parent1frac = (parent1.perp() / (parent0.perp() + parent1.perp()))*TMath::Power(angle_distance, alpha);

    double parent0Fitness = calcFitness(parent0, alpha);
    double parent1Fitness = calcFitness(parent1, alpha);

    parent0Blob = {parent0, parent0Fitness, parent0frac};
    parent1Blob = {parent1, parent1Fitness, parent1frac};

    // parent0Blob = calcFitness(parent0, alpha);
    // parent1Blob = calcFitness(parent1, alpha);
  }

  // else {
  //   parent0Blob = {parent0, 0, 0};
  //   parent1Blob = {parent1, 0, 0};
  // }

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

vector<PseudoJet> runBranchQueue(PseudoJet firstBranch, int n_axes, double alpha, double zcut) {
  priority_queue<BranchFitness, vector<BranchFitness>, CompareJetBranch> branchQueue;

  double firstBlobfitness = calcFitness(firstBranch, alpha);
  BranchFitness firstBlob = {firstBranch, firstBlobfitness, 0.0};
  branchQueue.push(firstBlob);

  vector<PseudoJet> branchAxes;
  BranchFitness currentBlob;

  vector<PseudoJet> tempAxes;
  int max_axes_size = 0;

  while (branchQueue.size() < n_axes) {

    // if (branchQueue.size() == 0) {
    //   cerr << "Error: branch queue empty." << endl;
    //   break;
    // }

    if (branchQueue.size() < max_axes_size) {
      cerr << "Warning: branch queue ending early." << endl;
      branchAxes = tempAxes;
      break;
    }

    currentBlob = branchQueue.top();
    if (branchQueue.size() > max_axes_size) {
      tempAxes.push_back(currentBlob.branch);
      max_axes_size = branchQueue.size();
    }
    branchQueue.pop();

    vector<BranchFitness> parentBlobs = splitBranch(currentBlob, alpha);
    if (parentBlobs[0].branch.has_structure() && parentBlobs[1].branch.has_structure()) {
      // cout << currentBlob.zfrac << endl;

      if (parentBlobs[0].zfrac < zcut/* || TMath::Abs(parentBlobs[0].branch.eta()) > 2*/) {
        branchQueue.push(parentBlobs[1]);
      }
      else if (parentBlobs[1].zfrac < zcut/* || TMath::Abs(parentBlobs[0].branch.eta()) > 2*/) {
        branchQueue.push(parentBlobs[0]);
      }
      else{
        branchQueue.push(parentBlobs[0]);
        branchQueue.push(parentBlobs[1]);
      }
    }
  } 

  int n_found = branchQueue.size();

  for (int i = 0; i < n_found; i++) {
    BranchFitness final_blob = branchQueue.top();    
    branchAxes.push_back(final_blob.branch);
    // cout << TMath::Sqrt(final_blob.branch.exclusive_subdmerge(1))*JetDefinition::max_allowable_R << endl;
    branchQueue.pop();
  }

  return branchAxes;

}

vector<PseudoJet> getBranchAxes(PseudoJet firstBranch, int n_axes, double alpha, double zcut) {
  // priority_queue<BranchFitness, vector<BranchFitness>, CompareJetBranch> branchQueue;
  // priority_queue<BranchFitness, vector<BranchFitness>, CompareJetBranch> branchQueue_left;
  // priority_queue<BranchFitness, vector<BranchFitness>, CompareJetBranch> branchQueue_right;

  // unordered_map<double, int> hashtable;
  // Make a blob from the first branch. We're not using this in the queue so
  // 0.0 is a dummy fitness value.

  // double firstBlobfitness = calcFitness(firstBranch, alpha);
  // BranchFitness firstBlob = {firstBranch, firstBlobfitness, 0.0};

  // vector<BranchFitness> parentBlobs = splitBranch(firstBlob, alpha);

  // branchQueue.push(firstBlob);
  // branchQueue.push(parentBlobs[0]);
  // branchQueue_right.push(parentBlobs[1]);
  
  // cout << parentBlobs[0].branch.eta() << " " << parentBlobs[0].branch.phi() << endl;
  // cout << parentBlobs[1].branch.eta() << " " << parentBlobs[1].branch.phi() << endl;

  vector<PseudoJet> branchAxes;

  vector<PseudoJet> initialAxes = runBranchQueue(firstBranch, 2, alpha, zcut);

  branchAxes = runBranchQueue(initialAxes[0], n_axes/2, alpha, zcut);
  vector<PseudoJet> tempBranchAxes = runBranchQueue(initialAxes[1], n_axes/2, alpha, zcut);
  branchAxes.insert(branchAxes.end(), tempBranchAxes.begin(), tempBranchAxes.end());
  
  return branchAxes;  
}

double nearest_neighbors(vector<PseudoJet> axes, int N = 3, int K = 2) {
// vector<Double_t> nearest_neighbors(vector<PseudoJet> axes, int N = 3, int K = 2) {
  // vector<Double_t> neighbors;

  std::string bitmask(K, 1); // K leading 1's
  bitmask.resize(N, 0); // N-K trailing 0's
 
  double avg_distance_new = 0;
  int axis_counter = 0;

  // print integers and permute bitmask
  do {
    vector<int> axis_indices;
    for (int i = 0; i < N; ++i) // [0..N-1] integers
    {
      if (bitmask[i]) axis_indices.push_back(i);
    }
    avg_distance_new += axes[axis_indices[0] + 3].delta_R(axes[axis_indices[1] + 3]);
    axis_counter++;
    // cout << axis_indices.size() << endl;
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
 
  avg_distance_new = avg_distance_new/(2*axis_counter);

  // for (int i = 0; i < axes.size()/3; i++) {
  //   Double_t distance_12 = axes[3*i].delta_R(axes[3*i+1]); 
  //   Double_t distance_13 = axes[3*i].delta_R(axes[3*i+2]); 
  //   Double_t distance_23 = axes[3*i + 1].delta_R(axes[3*i+2]); 
  //   Double_t neighbor1 = (distance_12 < distance_13) ? distance_12/2 : distance_13/2;
  //   Double_t neighbor2 = (distance_12 < distance_23) ? distance_12/2 : distance_23/2;
  //   Double_t neighbor3 = (distance_13 < distance_23) ? distance_13/2 : distance_23/2;
  //   Double_t avg_distance = (distance_12 + distance_13 + distance_23)/6;
  //   // neighbors.push_back(neighbor1);
  //   // neighbors.push_back(neighbor2);
  //   // neighbors.push_back(neighbor3);
  //   neighbors.push_back(avg_distance);
  //   neighbors.push_back(avg_distance);
  //   neighbors.push_back(avg_distance);

  //   cout << avg_distance << " ";
  // }
  // cout << endl;

  // return neighbors;
  return avg_distance_new
}

vector<PseudoJet> cluster_jets(vector<PseudoJet> my_particles, vector<PseudoJet> my_axes, vector<double> radii) {

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
      Double_t Rparam_subjet = radii[i_axis];
      Double_t axis_distance = my_particles[i_particles].delta_R(my_axes[i_axis]);
      if ((axis_distance < min_axis_distance) && (axis_distance < Rparam_subjet)) {
        min_axis_distance = axis_distance;
        min_axis_index = i_axis;
      }
    }
    if (min_axis_index != -1) {
      final_jets[min_axis_index] = join(final_jets[min_axis_index], my_particles[i_particles], *recombScheme);
    }
  }


}

int main(int argc, char* argv[]){

  TFile out("ttbar_recursive_test.root", "RECREATE");

  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW"); // Input File Name n
  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW"); // Input File Name

  std::string samples[2] = {"/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW", "/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW"};

  TH1F *invariant_mass_ttbar = new TH1F("invariant_mass_ttbar", "Invariant Mass for ttbar", 1000, 0, 1000);
  TH1F *invariant_mass_dijets = new TH1F("invariant_mass_dijets", "Invariant Mass for dijets", 1000, 0, 1000);

  TH1F *tau32_ttbar = new TH1F("tau32_ttbar", "tau32 for ttbar", 100, 0, 1);
  TH1F *tau32_dijets = new TH1F("tau32_dijets", "tau32 for dijets", 100, 0, 1);

  double alpha = 0.0;
  double zcut = 0.1;

  for (int i_sample = 0; i_sample < 2; i_sample++) {

  const char* current_data = samples[i_sample].c_str();
  ifstream inputStream(current_data); // Input File Name

  const int n_event = 10; // # of events to be analyzed   

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

      double Rparam = JetDefinition::max_allowable_R;
      Strategy strategy = Best;
      const JetDefinition::Recombiner *recombScheme = new WinnerTakeAllRecombiner();
      JetDefinition *jetDef = NULL;
      jetDef = new JetDefinition(cambridge_algorithm, Rparam, recombScheme, strategy);
      ClusterSequence clustSeq(input_particles, *jetDef);
      vector <PseudoJet> inclusiveJets = clustSeq.inclusive_jets(0);
      PseudoJet total_jet = inclusiveJets[0];

      vector<PseudoJet> axes = getBranchAxes(total_jet, 6, alpha, zcut);
      vector<Double_t> axis_neighbors = nearest_neighbors(axes);
      vector<PseudoJet> final_jets;

      // INITIALIZE JETS
      for (int i = 0; i < axes.size(); i++) {
        axes_display->Fill(axes[i].eta(), axes[i].phi());
        PseudoJet empty_jet(0,0,0,0);
        final_jets.push_back(empty_jet);
      }

      // TAKE INITIAL PARTICLES AND CLUSTER THEM ACCORDING TO NEAREST AXIS
      for (int i_particles = 0; i_particles < input_particles.size(); i_particles++) {
        event_display->Fill(input_particles[i_particles].eta(), input_particles[i_particles].phi(), input_particles[i_particles].perp());
        Double_t min_axis_distance = 10000;
        int min_axis_index = -1;
        for (int i_axis = 0; i_axis < axes.size(); i_axis++) {
          Double_t Rparam_subjet = axis_neighbors[i_axis];
          Double_t axis_distance = input_particles[i_particles].delta_R(axes[i_axis]);
          if ((axis_distance < min_axis_distance) && (axis_distance < Rparam_subjet)) {
            min_axis_distance = axis_distance;
            min_axis_index = i_axis;
          }
        }
        if (min_axis_index != -1) {
          final_jets[min_axis_index] = join(final_jets[min_axis_index], input_particles[i_particles], *recombScheme);
        }
      }

      // FILL EVENT WITH GHOSTS TO MAKE FIREWORKS
      Double_t ghost_perp = 0.001;
      Double_t n_ghosts = 200;
      for (int i_ghosts = 0; i_ghosts < n_ghosts; i_ghosts++) {
        for (int j_ghosts = 0; j_ghosts < n_ghosts; j_ghosts++) {
          Double_t ghost_eta = -5 + (double)(10/n_ghosts)*i_ghosts;
          Double_t ghost_phi = (double)(2*TMath::Pi()/n_ghosts)*j_ghosts;
          PseudoJet ghost(ghost_perp*TMath::Cos(ghost_phi),ghost_perp*TMath::Sin(ghost_phi),
             ghost_perp*TMath::SinH(ghost_eta),ghost_perp*TMath::CosH(ghost_eta));
          Double_t min_axis_distance = 10000;
          int min_axis_index = -1;
          for (int i_axis = 0; i_axis < axes.size(); i_axis++) {
            Double_t Rparam_subjet = axis_neighbors[i_axis];
            Double_t axis_distance = ghost.delta_R(axes[i_axis]);
            if (axis_distance < min_axis_distance && axis_distance < Rparam_subjet) {
              min_axis_distance = axis_distance;
              min_axis_index = i_axis;
            }
          }
          // if (min_axis_index != -1) final_jets[min_axis_index] = join(final_jets[min_axis_index], ghost, *recombScheme);
          if (min_axis_index == 0) subjet1_display->Fill(ghost_eta, ghost_phi);
          if (min_axis_index == 1) subjet2_display->Fill(ghost_eta, ghost_phi);
          if (min_axis_index == 2) subjet3_display->Fill(ghost_eta, ghost_phi);
          if (min_axis_index == 3) subjet4_display->Fill(ghost_eta, ghost_phi);
          if (min_axis_index == 4) subjet5_display->Fill(ghost_eta, ghost_phi);
          if (min_axis_index == 5) subjet6_display->Fill(ghost_eta, ghost_phi);
        }
      }

      // SUM EACH JET TRIPLET TO FIND INVARIANT MASS AND NSUBJETTINESS
      vector<PseudoJet> big_jets;
      for (int i = 0; i < final_jets.size()/3; i++) {
        PseudoJet summed_jet = join(final_jets[3*i],final_jets[3*i + 1],final_jets[3*i + 2]);
        big_jets.push_back(summed_jet);
        if (i_sample == 0) invariant_mass_ttbar->Fill(summed_jet.m());
        if (i_sample == 1) invariant_mass_dijets->Fill(summed_jet.m());
      }

      Nsubjettiness tau3(3, Njettiness::manual_axes, Njettiness::unnormalized_measure, 1.0);
      Nsubjettiness tau2(2, Njettiness::manual_axes, Njettiness::unnormalized_measure, 1.0);

      vector<PseudoJet> axes_4 = getBranchAxes(total_jet, 4, alpha, zcut);

      for (int i_jets = 0; i_jets < big_jets.size(); i_jets++) {
        vector<PseudoJet> threeprong_axes;
        threeprong_axes.push_back(axes[3*i_jets]);
        threeprong_axes.push_back(axes[3*i_jets + 1]);
        threeprong_axes.push_back(axes[3*i_jets + 2]);
        tau3.setAxes(threeprong_axes);

        vector<PseudoJet> twoprong_axes;
        twoprong_axes.push_back(axes[2*i_jets]);
        twoprong_axes.push_back(axes[2*i_jets + 1]);
        tau2.setAxes(twoprong_axes);

        if (big_jets[i_jets].m() > 160 && big_jets[i_jets].m() < 240) {
          if (i_sample == 0) tau32_ttbar->Fill(tau3.result(big_jets[i_jets])/tau2.result(big_jets[i_jets]));
          if (i_sample == 1) tau32_dijets->Fill(tau3.result(big_jets[i_jets])/tau2.result(big_jets[i_jets]));
        }
      }

      // for (int i_jets = 0; i_jets < final_jets.size(); i_jets++) {
      //   PseudoJet temp_jet = final_jets[i_jets];
      //   if (temp_jet.has_constituents()) {
      //     for (int i_const = 0; i_const < temp_jet.constituents().size(); i_const++) {
      //       if (i_jets == 0) subjet1_display->Fill(temp_jet.constituents()[i_const].eta(), temp_jet.constituents()[i_const].phi());
      //       if (i_jets == 1) subjet2_display->Fill(temp_jet.constituents()[i_const].eta(), temp_jet.constituents()[i_const].phi());
      //       if (i_jets == 2) subjet3_display->Fill(temp_jet.constituents()[i_const].eta(), temp_jet.constituents()[i_const].phi());
      //       if (i_jets == 3) subjet4_display->Fill(temp_jet.constituents()[i_const].eta(), temp_jet.constituents()[i_const].phi());
      //       if (i_jets == 4) subjet5_display->Fill(temp_jet.constituents()[i_const].eta(), temp_jet.constituents()[i_const].phi());
      //       if (i_jets == 5) subjet6_display->Fill(temp_jet.constituents()[i_const].eta(), temp_jet.constituents()[i_const].phi());
      //     }
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
      // event_display->Draw("box");
      // axes_display->Draw("SAMES");

      // ghost_display->Draw("boxSAMES");

      // subjet1_display->Draw("SAMES");
      // subjet2_display->Draw("SAMES");
      // subjet3_display->Draw("SAMES");
      // subjet4_display->Draw("SAMES");
      // subjet5_display->Draw("SAMES");
      // subjet6_display->Draw("SAMES");

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

  invariant_mass_ttbar->Write();
  invariant_mass_dijets->Write();

  tau32_ttbar->Write();
  tau32_dijets->Write();

//   ttbar_tau3_tau2_wtakt_beta1->Write();
//   ttbar_tau3_tau2_declustering_beta1->Write();
//   dijets_tau3_tau2_wtakt_beta1->Write();
//   dijets_tau3_tau2_declustering_beta1->Write();

//   Double_t ttbar_tau3_tau2_wtakt_beta1_scale = 1/ttbar_tau3_tau2_wtakt_beta1->Integral();
//   Double_t ttbar_tau3_tau2_declustering_beta1_scale = 1/ttbar_tau3_tau2_declustering_beta1->Integral();
//   Double_t dijets_tau3_tau2_wtakt_beta1_scale = 1/dijets_tau3_tau2_wtakt_beta1->Integral();
//   Double_t dijets_tau3_tau2_declustering_beta1_scale = 1/dijets_tau3_tau2_declustering_beta1->Integral();

//   ttbar_tau3_tau2_wtakt_beta1->Scale(ttbar_tau3_tau2_wtakt_beta1_scale);
//   ttbar_tau3_tau2_declustering_beta1->Scale(ttbar_tau3_tau2_declustering_beta1_scale);
//   dijets_tau3_tau2_wtakt_beta1->Scale(dijets_tau3_tau2_wtakt_beta1_scale);
//   dijets_tau3_tau2_declustering_beta1->Scale(dijets_tau3_tau2_declustering_beta1_scale);

//   int n_wtakt = ttbar_tau3_tau2_wtakt_beta1->GetSize() - 2;
//   double integral_dijets_wtakt[n_wtakt], integral_ttbar_wtakt[n_wtakt];
//   for (int i = 0; i < n_wtakt; i++) {
//     integral_ttbar_wtakt[i] = ttbar_tau3_tau2_wtakt_beta1->Integral(0,i);
//     integral_dijets_wtakt[i] = dijets_tau3_tau2_declustering_beta1->Integral(0,i);
//   }
//   TGraph* ROC_tau3_tau2_wtakt = new TGraph(n_wtakt, integral_ttbar_wtakt, integral_dijets_wtakt);
//   ROC_tau3_tau2_wtakt->GetXaxis()->SetLimits(0, 1);
//   ROC_tau3_tau2_wtakt->GetYaxis()->SetLimits(0, 1);
//   ROC_tau3_tau2_wtakt->SetLineColor(kBlack);
//   ROC_tau3_tau2_wtakt->SetLineWidth(2);
//   ROC_tau3_tau2_wtakt->SetMarkerStyle(5);
//   ROC_tau3_tau2_wtakt->SetMarkerSize(2);
//   ROC_tau3_tau2_wtakt->Write();

//   int n_declustering = ttbar_tau3_tau2_declustering_beta1->GetSize() - 2;
//   double integral_dijets_declustering[n_declustering], integral_ttbar_declustering[n_declustering];
//   for (int i = 0; i < n_declustering; i++) {
//     integral_ttbar_declustering[i] = ttbar_tau3_tau2_declustering_beta1->Integral(0,i);
//     integral_dijets_declustering[i] = dijets_tau3_tau2_declustering_beta1->Integral(0,i);
//   }
//   TGraph* ROC_tau3_tau2_declustering = new TGraph(n_declustering, integral_ttbar_declustering, integral_dijets_declustering);
//   ROC_tau3_tau2_declustering->GetXaxis()->SetLimits(0, 1);
//   ROC_tau3_tau2_declustering->GetYaxis()->SetLimits(0, 1);
//   ROC_tau3_tau2_declustering->SetLineColor(kRed);
//   ROC_tau3_tau2_declustering->SetLineWidth(2);
//   ROC_tau3_tau2_declustering->SetMarkerStyle(5);
//   ROC_tau3_tau2_declustering->SetMarkerSize(2);
//   ROC_tau3_tau2_declustering->Write();

//   TCanvas *ROC_compare = new TCanvas("ROC_compare", "ROC Curve Comparison", 600, 600);
//   ROC_compare->cd();
//   ROC_compare->SetLogy();
//   TMultiGraph *ROC_multigraph = new TMultiGraph("ROC_multigraph", "ROC Comparison for #tau_{3}/#tau_{2} (500 < p_{T} < 600)");
//   ROC_multigraph->Add(ROC_tau3_tau2_wtakt);
//   ROC_multigraph->Add(ROC_tau3_tau2_declustering);
//   ROC_multigraph->Draw("AL");
//   TLegend *leg_ROC = new TLegend(0.1, 0.7, 0.3, 0.9);
//   leg_ROC->AddEntry(ROC_tau3_tau2_wtakt, "WTA kT", "L");
//   leg_ROC->AddEntry(ROC_tau3_tau2_declustering, "Declustering", "L");
//   leg_ROC->Draw();
//   ROC_compare->Write();

}