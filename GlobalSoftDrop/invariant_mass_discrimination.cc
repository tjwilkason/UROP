//Misc. Headers
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>

//Fastjet headers
#include "FastJet3.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/CompositeJetStructure.hh"
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
#include "TLine.h"
#include "TString.h"

using namespace std;
using namespace Pythia8;
using namespace fastjet;
using namespace fastjet::contrib;
 
bool inRange(double val, double x, double y) {
  return val > x && val < y;
}

int main(int argc, char* argv[]){

  // TFile out("volatility_test_nomasscut_p2p5.root", "RECREATE");
  // TFile out("volatility_test_masscuts_p2p5_testing_more_2.root", "RECREATE");
  TFile out("phase_space_testing_singlemasscut.root", "RECREATE");
  // TFile out("volatility_test.root", "RECREATE");

  TH1F *ttbar_wmass_distribution_rp2 = new TH1F("ttbar_wmass_distribution_rp2", "Histogram of W mass", 60, 0, 300);
  TH1F *dijets_wmass_distribution_rp2 = new TH1F("dijets_wmass_distribution_rp2", "Histogram of W mass", 60, 0, 300);
  TH1F *ttbar_wmass_distribution_rp5 = new TH1F("ttbar_wmass_distribution_rp5", "Histogram of W mass", 60, 0, 300);
  TH1F *dijets_wmass_distribution_rp5 = new TH1F("dijets_wmass_distribution_rp5", "Histogram of W mass", 60, 0, 300);

  // TString hist_title = "Histogram of volatility" + event_string.str();
  TH1F *volatility_hist_ttbar_nocut = new TH1F("volatility_hist_ttbar_nocut", "Histogram of volatility (ttbar)", 100, 0, 1);
  TH1F *volatility_hist_dijets_nocut = new TH1F("volatility_hist_dijets_nocut", "Histogram of volatility (dijets)", 100, 0, 1);
  TH1F *volatility_hist_ttbar_masscut = new TH1F("volatility_hist_ttbar_masscut", "Histogram of volatility (ttbar)", 100, 0, 1);
  TH1F *volatility_hist_dijets_masscut = new TH1F("volatility_hist_dijets_masscut", "Histogram of volatility (dijets)", 100, 0, 1);
  TH1F *volatility_hist_ttbar_nsubcut = new TH1F("volatility_hist_ttbar_nsubcut", "Histogram of volatility (ttbar)", 100, 0, 1);
  TH1F *volatility_hist_dijets_nsubcut = new TH1F("volatility_hist_dijets_nsubcut", "Histogram of volatility (dijets)", 100, 0, 1);

  TH1F *tau32_wta_ttbar_nocut = new TH1F("tau32_wta_ttbar_nocut", "Histogram of tau32_wta (ttbar)", 100, 0, 1);
  TH1F *tau32_wta_dijets_nocut = new TH1F("tau32_wta_dijets_nocut", "Histogram of tau32_wta (dijets)", 100, 0, 1);
  TH1F *tau32_wta_ttbar_masscut = new TH1F("tau32_wta_ttbar_masscut", "Histogram of tau32_wta (ttbar)", 100, 0, 1);
  TH1F *tau32_wta_dijets_masscut = new TH1F("tau32_wta_dijets_masscut", "Histogram of tau32_wta (dijets)", 100, 0, 1);
  TH1F *tau32_wta_ttbar_vol2cut = new TH1F("tau32_wta_ttbar_vol2cut", "Histogram of tau32_wta (ttbar)", 100, 0, 1);
  TH1F *tau32_wta_dijets_vol2cut = new TH1F("tau32_wta_dijets_vol2cut", "Histogram of tau32_wta (dijets)", 100, 0, 1);
  TH1F *tau32_wta_ttbar_vol4cut = new TH1F("tau32_wta_ttbar_vol4cut", "Histogram of tau32_wta (ttbar)", 100, 0, 1);
  TH1F *tau32_wta_dijets_vol4cut = new TH1F("tau32_wta_dijets_vol4cut", "Histogram of tau32_wta (dijets)", 100, 0, 1);
  TH1F *tau32_wta_ttbar_vol6cut = new TH1F("tau32_wta_ttbar_vol6cut", "Histogram of tau32_wta (ttbar)", 100, 0, 1);
  TH1F *tau32_wta_dijets_vol6cut = new TH1F("tau32_wta_dijets_vol6cut", "Histogram of tau32_wta (dijets)", 100, 0, 1);
  TH1F *tau32_wta_ttbar_vol8cut = new TH1F("tau32_wta_ttbar_vol8cut", "Histogram of tau32_wta (ttbar)", 100, 0, 1);
  TH1F *tau32_wta_dijets_vol8cut = new TH1F("tau32_wta_dijets_vol8cut", "Histogram of tau32_wta (dijets)", 100, 0, 1);
  TH1F *tau32_wta_ttbar_vol10cut = new TH1F("tau32_wta_ttbar_vol10cut", "Histogram of tau32_wta (ttbar)", 100, 0, 1);
  TH1F *tau32_wta_dijets_vol10cut = new TH1F("tau32_wta_dijets_vol10cut", "Histogram of tau32_wta (dijets)", 100, 0, 1);

  TH2F *phase_space_ttbar = new TH2F("phase_space", "Phase Space for Volatility vs. Tau32", 100, 0, 1, 100, 0, 1);
  TH2F *phase_space_dijets = new TH2F("phase_space", "Phase Space for Volatility vs. Tau32", 100, 0, 1, 100, 0, 1);

  int total_ttbar_jets = 0;
  int total_dijets_jets = 0;

  const int n_samples = 2;
  std::string samples[n_samples] = {"/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW", "/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW"};
  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-dijets-pt0500-0600.UW"); // Input File Name
  // ifstream inputStream("/home/tjwilk/Physics/thaler_urop/BOOST_samples/herwig65-lhc7-ttbar2hadrons-pt0500-0600.UW"); // Input File Name
  int n_R = 2;

  // int ttbar_jets_nsub_cut = 0;
  // int dijets_jets_nsub_cut = 0;

  // int ttbar_jets_rad_cut = 0;
  // int dijets_jets_rad_cut = 0;

  // int ttbar_jets_vol_cut = 0;
  // int dijets_jets_vol_cut = 0;

  // int ttbar_double_cut = 0;
  // int dijets_double_cut = 0;

  for (int i_sample = 0; i_sample < n_samples; i_sample++) {

    const char* current_data = samples[i_sample].c_str();
    ifstream inputStream(current_data); // Input File Name

    const int n_event = 10000; // # of events to be analyzed   

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

        double Rparam = 1.0;
        Strategy strategy = Best;
        RecombinationScheme recombScheme = E_scheme;
        JetDefinition *jetDef = NULL;
        jetDef = new JetDefinition(antikt_algorithm, Rparam, recombScheme, strategy);

        //Run Fastjet algorithm
        vector <PseudoJet> inclusiveJets, sortedJets, centralJets, massiveJets;
        ClusterSequence clustSeq(input_particles, *jetDef);

        //Extract inclusive jets sorted by pT
        inclusiveJets = clustSeq.inclusive_jets(200.0);

        // Sort central jets by pT
        sortedJets = sorted_by_pt(inclusiveJets);

        //Create selectors for hardest jets in the event
        Selector jet_selector = SelectorNHardest(2);
        vector<PseudoJet> hardest_jet = jet_selector(sortedJets);

        for (int i = 0; i < hardest_jet.size(); i++) {

          TH1F *event_mass = new TH1F("event_mass", "Mass distribution for individual event", 300, 0, 300);

          NsubjettinessRatio nsub_3_2_wta_kt(3, 2, Njettiness::wta_kt_axes, Njettiness::normalized_cutoff_measure, 1.0, 1.0, 1.0);

          int mass1, mass2, Wmass1, Wmass2;

          for (int i_R = 1; i_R <= n_R; i_R++) {

            double R_value;
            if (i_R == 1) R_value = 0.2;
            if (i_R == 2) R_value = 0.5;

            NjettinessPlugin njet_plugin_wta(3, Njettiness::wta_kt_axes, Njettiness::unnormalized_cutoff_measure, 1.0, R_value);
            JetDefinition njet_def_wta(&njet_plugin_wta);

            vector<PseudoJet> temp_constituents = hardest_jet[i].constituents();

            ClusterSequence njet_cluster_wta(temp_constituents, njet_def_wta);
            vector<PseudoJet> jets_wta, axes_wta;
            jets_wta = njet_cluster_wta.inclusive_jets();

            PseudoJet summed_jets;
            for (int n_jets = 0; n_jets < jets_wta.size(); n_jets++) {
              summed_jets += jets_wta[n_jets];
            }

            PseudoJet W_jet;
            if ((jets_wta[0].perp() < jets_wta[2].perp()) && jets_wta[1].perp() < jets_wta[2].perp()) W_jet = jets_wta[0] + jets_wta[1];
            else if ((jets_wta[0].perp() < jets_wta[1].perp()) && jets_wta[2].perp() < jets_wta[1].perp()) W_jet = jets_wta[0] + jets_wta[2];
            else if ((jets_wta[1].perp() < jets_wta[0].perp()) && jets_wta[2].perp() < jets_wta[0].perp()) W_jet = jets_wta[1] + jets_wta[2];

            if (i_sample == 0) {
              if (i_R == 1) dijets_wmass_distribution_rp2->Fill(W_jet.m());
              if (i_R == 2) dijets_wmass_distribution_rp5->Fill(W_jet.m());
            }
            if (i_sample == 1) {
              if (i_R == 1) ttbar_wmass_distribution_rp2->Fill(W_jet.m());
              if (i_R == 2) ttbar_wmass_distribution_rp5->Fill(W_jet.m());
            }

            event_mass->Fill(summed_jets.m());
            if (i_R == 1) {
              mass1 = summed_jets.m();
              Wmass1 = W_jet.m();
            }
            if (i_R == 2) {
              mass2 = summed_jets.m();
              Wmass2 = W_jet.m();
            }

            // in case axes information is also desired
            const NjettinessExtras *extras_wta = njettiness_extras(njet_cluster_wta);
            axes_wta = extras_wta->axes();
          }

          double mean_mass = (mass1 + mass2)/2;
          double mean_Wmass = (Wmass1 + Wmass2)/2;
          // double volatility = event_mass->GetRMS()/event_mass->GetMean();
          double volatility = (TMath::Abs(mass1 - mass2))/mean_mass;
          if (i_sample == 0) {
            total_dijets_jets++;

            volatility_hist_dijets_nocut->Fill(volatility);
            if (inRange(hardest_jet[i].m(), 160, 240) && inRange(mean_mass, 160, 240) && inRange(mean_Wmass, 60, 100)) {
              volatility_hist_dijets_masscut->Fill(volatility);
              if (nsub_3_2_wta_kt.result(hardest_jet[i]) < 0.7) {
                volatility_hist_dijets_nsubcut->Fill(volatility);
              }
            }

            tau32_wta_dijets_nocut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
            if (inRange(mean_mass, 160, 240)/* && inRange(mean_Wmass, 60, 100)*/) {
              tau32_wta_dijets_masscut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
              if (volatility < 0.2) {
                tau32_wta_dijets_vol2cut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
              }
              if (volatility < 0.4) {
                tau32_wta_dijets_vol4cut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
              }
              if (volatility < 0.6) {
                tau32_wta_dijets_vol6cut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
              }
              if (volatility < 0.8) {
                tau32_wta_dijets_vol8cut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
              }
              if (volatility < 1.0) {
                tau32_wta_dijets_vol10cut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
              }
            }

            // if (inRange(hardest_jet[i].m(), 160, 240) && inRange(mean_mass, 160, 240) && inRange(mean_Wmass, 60, 100)) {
              phase_space_dijets->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]), volatility);
            // }
          }
          else if (i_sample == 1) {
            total_ttbar_jets++;

            volatility_hist_ttbar_nocut->Fill(volatility);
            if (inRange(hardest_jet[i].m(), 160, 240) && inRange(mean_mass, 160, 240) && inRange(mean_Wmass, 60, 100)) {
              volatility_hist_ttbar_masscut->Fill(volatility);
              if (nsub_3_2_wta_kt.result(hardest_jet[i]) < 0.7) {
                volatility_hist_ttbar_nsubcut->Fill(volatility);
              }          
            }

            tau32_wta_ttbar_nocut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
            if (inRange(mean_mass, 160, 240)/* && inRange(mean_Wmass, 60, 100)*/) {
              tau32_wta_ttbar_masscut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
              if (volatility < 0.2) {
                tau32_wta_ttbar_vol2cut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
              }
              if (volatility < 0.4) {
                tau32_wta_ttbar_vol4cut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
              }
              if (volatility < 0.6) {
                tau32_wta_ttbar_vol6cut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
              }
              if (volatility < 0.8) {
                tau32_wta_ttbar_vol8cut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
              }
              if (volatility < 1.0) {
                tau32_wta_ttbar_vol10cut->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]));
              }
            }

            // if (inRange(hardest_jet[i].m(), 160, 240) && inRange(mean_mass, 160, 240) && inRange(mean_Wmass, 60, 100)) {
              phase_space_ttbar->Fill(nsub_3_2_wta_kt.result(hardest_jet[i]), volatility);
            // }
          }
          delete event_mass;
        }

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

  volatility_hist_ttbar_nocut->SetStats(0);
  volatility_hist_ttbar_nocut->Write();
  volatility_hist_dijets_nocut->SetStats(0);
  volatility_hist_dijets_nocut->Write();
  volatility_hist_ttbar_masscut->SetStats(0);
  volatility_hist_ttbar_masscut->Write();
  volatility_hist_dijets_masscut->SetStats(0);
  volatility_hist_dijets_masscut->Write();
  volatility_hist_ttbar_nsubcut->SetStats(0);
  volatility_hist_ttbar_nsubcut->Write();
  volatility_hist_dijets_nsubcut->SetStats(0);
  volatility_hist_dijets_nsubcut->Write();

  tau32_wta_ttbar_nocut->SetStats(0);
  tau32_wta_ttbar_nocut->Write();
  tau32_wta_dijets_nocut->SetStats(0);
  tau32_wta_dijets_nocut->Write();
  tau32_wta_ttbar_masscut->SetStats(0);
  tau32_wta_ttbar_masscut->Write();
  tau32_wta_dijets_masscut->SetStats(0);
  tau32_wta_dijets_masscut->Write();
  tau32_wta_ttbar_vol2cut->SetStats(0);
  tau32_wta_ttbar_vol2cut->Write();
  tau32_wta_dijets_vol2cut->SetStats(0);
  tau32_wta_dijets_vol2cut->Write();
  tau32_wta_ttbar_vol4cut->SetStats(0);
  tau32_wta_ttbar_vol4cut->Write();
  tau32_wta_dijets_vol4cut->SetStats(0);
  tau32_wta_dijets_vol4cut->Write();
  tau32_wta_ttbar_vol6cut->SetStats(0);
  tau32_wta_ttbar_vol6cut->Write();
  tau32_wta_dijets_vol6cut->SetStats(0);
  tau32_wta_dijets_vol6cut->Write();
  tau32_wta_ttbar_vol8cut->SetStats(0);
  tau32_wta_ttbar_vol8cut->Write();
  tau32_wta_dijets_vol8cut->SetStats(0);
  tau32_wta_dijets_vol8cut->Write();
  tau32_wta_ttbar_vol10cut->SetStats(0);
  tau32_wta_ttbar_vol10cut->Write();
  tau32_wta_dijets_vol10cut->SetStats(0);
  tau32_wta_dijets_vol10cut->Write();
  phase_space_ttbar->SetStats(0);
  phase_space_dijets->SetStats(0);

  ttbar_wmass_distribution_rp2->SetStats(0);
  ttbar_wmass_distribution_rp2->Write();
  ttbar_wmass_distribution_rp5->SetStats(0);
  ttbar_wmass_distribution_rp5->Write();
  dijets_wmass_distribution_rp2->SetStats(0);
  dijets_wmass_distribution_rp2->Write();
  dijets_wmass_distribution_rp5->SetStats(0);
  dijets_wmass_distribution_rp5->Write();

  // cout << "ttbar nsub 0.6 efficiency: " << (double)ttbar_jets_nsub_cut/total_ttbar_jets << endl;
  // cout << "dijet nsub 0.6 efficiency: " << (double)dijets_jets_nsub_cut/total_dijets_jets << endl;

  // cout << "ttbar vol cut 0.2 efficiency: " << (double)ttbar_jets_vol_cut/total_ttbar_jets << endl;
  // cout << "dijet vol cut 0.2 efficiency: " << (double)dijets_jets_vol_cut/total_dijets_jets << endl;

  // cout << "ttbar double cut efficiency: " << (double)ttbar_double_cut/total_ttbar_jets << endl;
  // cout << "dijet double cut efficiency: " << (double)dijets_double_cut/total_dijets_jets << endl;

  TCanvas *phase_space_can = new TCanvas("phase_space_can", "Phase Space", 600, 600);
  phase_space_dijets->SetMarkerColor(kRed);
  phase_space_ttbar->SetMarkerColor(kBlue);
  phase_space_dijets->Draw();
  phase_space_ttbar->Draw("SAME");
  phase_space_can->Write();  

  int nvol_bins_nocut = volatility_hist_dijets_nocut->GetSize() - 2;
  double integral_dijets_vol_nocut[nvol_bins_nocut], integral_ttbar_vol_nocut[nvol_bins_nocut];
  for (int i = 0; i < nvol_bins_nocut; i++) {
    integral_dijets_vol_nocut[i] = volatility_hist_dijets_nocut->Integral(0, i)/total_dijets_jets;
    integral_ttbar_vol_nocut[i] = volatility_hist_ttbar_nocut->Integral(0, i)/total_ttbar_jets;
  }
  TGraph* ROC_vol_nocut = new TGraph(nvol_bins_nocut, integral_ttbar_vol_nocut, integral_dijets_vol_nocut);
  ROC_vol_nocut->GetXaxis()->SetLimits(0, 1);
  ROC_vol_nocut->GetYaxis()->SetLimits(0, 1);
  ROC_vol_nocut->SetLineColor(kBlack);
  // ROC_vol_nocut->SetLineWidth(2);
  ROC_vol_nocut->SetMarkerStyle(5);
  ROC_vol_nocut->SetMarkerSize(2);
  ROC_vol_nocut->Write();

  int nsub_bins_nocut = tau32_wta_dijets_nocut->GetSize() - 2;
  double integral_dijets_nsub_nocut[nsub_bins_nocut], integral_ttbar_nsub_nocut[nsub_bins_nocut];
  for (int i = 0; i < nsub_bins_nocut; i++) {
    integral_dijets_nsub_nocut[i] = tau32_wta_dijets_nocut->Integral(0, i)/total_dijets_jets;
    integral_ttbar_nsub_nocut[i] = tau32_wta_ttbar_nocut->Integral(0, i)/total_ttbar_jets;
  }
  TGraph* ROC_nsub_nocut = new TGraph(nsub_bins_nocut, integral_ttbar_nsub_nocut, integral_dijets_nsub_nocut);
  ROC_nsub_nocut->GetXaxis()->SetLimits(0, 1);
  ROC_nsub_nocut->GetYaxis()->SetLimits(0, 1);
  ROC_nsub_nocut->SetLineColor(kBlack);
  // ROC_nsub_nocut->SetLineWidth(2);
  ROC_nsub_nocut->SetMarkerStyle(5);
  ROC_nsub_nocut->SetMarkerSize(2);
  ROC_nsub_nocut->Write();

  int nvol_bins_140_200 = volatility_hist_dijets_masscut->GetSize() - 2;
  double integral_dijets_vol_140_200[nvol_bins_140_200], integral_ttbar_vol_140_200[nvol_bins_140_200];
  for (int i = 0; i < nvol_bins_140_200; i++) {
    integral_dijets_vol_140_200[i] = volatility_hist_dijets_masscut->Integral(0, i)/total_dijets_jets;
    integral_ttbar_vol_140_200[i] = volatility_hist_ttbar_masscut->Integral(0, i)/total_ttbar_jets;
  }
  TGraph* ROC_vol_masscut = new TGraph(nvol_bins_140_200, integral_ttbar_vol_140_200, integral_dijets_vol_140_200);
  ROC_vol_masscut->GetXaxis()->SetLimits(0, 1);
  ROC_vol_masscut->GetYaxis()->SetLimits(0, 1);
  ROC_vol_masscut->SetLineColor(kBlack);
  // ROC_vol_masscut->SetLineWidth(2);
  ROC_vol_masscut->SetMarkerStyle(5);
  ROC_vol_masscut->SetMarkerSize(2);
  ROC_vol_masscut->Write();

  int nsub_bins_140_200 = tau32_wta_dijets_masscut->GetSize() - 2;
  double integral_dijets_nsub_140_200[nsub_bins_140_200], integral_ttbar_nsub_140_200[nsub_bins_140_200];
  for (int i = 0; i < nsub_bins_140_200; i++) {
    integral_dijets_nsub_140_200[i] = tau32_wta_dijets_masscut->Integral(0, i)/total_dijets_jets;
    integral_ttbar_nsub_140_200[i] = tau32_wta_ttbar_masscut->Integral(0, i)/total_ttbar_jets;
  }
  TGraph* ROC_nsub_masscut = new TGraph(nsub_bins_140_200, integral_ttbar_nsub_140_200, integral_dijets_nsub_140_200);
  ROC_nsub_masscut->GetXaxis()->SetLimits(0, 1);
  ROC_nsub_masscut->GetYaxis()->SetLimits(0, 1);
  ROC_nsub_masscut->SetLineColor(kBlack);
  // ROC_nsub_masscut->SetLineWidth(2);
  ROC_nsub_masscut->SetMarkerStyle(5);
  ROC_nsub_masscut->SetMarkerSize(2);
  ROC_nsub_masscut->Write();

  int nvol_bins_160_180 = volatility_hist_dijets_nsubcut->GetSize() - 2;
  double integral_dijets_vol_160_180[nvol_bins_160_180], integral_ttbar_vol_160_180[nvol_bins_160_180];
  for (int i = 0; i < nvol_bins_160_180; i++) {
    integral_dijets_vol_160_180[i] = volatility_hist_dijets_nsubcut->Integral(0, i)/total_dijets_jets;
    integral_ttbar_vol_160_180[i] = volatility_hist_ttbar_nsubcut->Integral(0, i)/total_ttbar_jets;
  }
  TGraph* ROC_vol_nsubcut = new TGraph(nvol_bins_160_180, integral_ttbar_vol_160_180, integral_dijets_vol_160_180);
  ROC_vol_nsubcut->GetXaxis()->SetLimits(0, 1);
  ROC_vol_nsubcut->GetYaxis()->SetLimits(0, 1);
  ROC_vol_nsubcut->SetLineColor(kBlack);
  // ROC_vol_nsubcut->SetLineWidth(2);
  ROC_vol_nsubcut->SetMarkerStyle(5);
  ROC_vol_nsubcut->SetMarkerSize(2);
  ROC_vol_nsubcut->Write();

  int nsub_bins_vol2cut = tau32_wta_dijets_vol2cut->GetSize() - 2;
  double integral_dijets_nsub_vol2cut[nsub_bins_vol2cut], integral_ttbar_nsub_vol2cut[nsub_bins_vol2cut];
  for (int i = 0; i < nsub_bins_vol2cut; i++) {
    integral_dijets_nsub_vol2cut[i] = tau32_wta_dijets_vol2cut->Integral(0, i)/total_dijets_jets;
    integral_ttbar_nsub_vol2cut[i] = tau32_wta_ttbar_vol2cut->Integral(0, i)/total_ttbar_jets;
  }
  TGraph* ROC_nsub_vol2cut = new TGraph(nsub_bins_vol2cut, integral_ttbar_nsub_vol2cut, integral_dijets_nsub_vol2cut);
  ROC_nsub_vol2cut->GetXaxis()->SetLimits(0, 1);
  ROC_nsub_vol2cut->GetYaxis()->SetLimits(0, 1);
  ROC_nsub_vol2cut->SetLineColor(kBlack);
  // ROC_nsub_vol2cut->SetLineWidth(2);
  ROC_nsub_vol2cut->SetMarkerStyle(5);
  ROC_nsub_vol2cut->SetMarkerSize(2);
  ROC_nsub_vol2cut->Write();

  int nsub_bins_vol4cut = tau32_wta_dijets_vol4cut->GetSize() - 2;
  double integral_dijets_nsub_vol4cut[nsub_bins_vol4cut], integral_ttbar_nsub_vol4cut[nsub_bins_vol4cut];
  for (int i = 0; i < nsub_bins_vol4cut; i++) {
    integral_dijets_nsub_vol4cut[i] = tau32_wta_dijets_vol4cut->Integral(0, i)/total_dijets_jets;
    integral_ttbar_nsub_vol4cut[i] = tau32_wta_ttbar_vol4cut->Integral(0, i)/total_ttbar_jets;
  }
  TGraph* ROC_nsub_vol4cut = new TGraph(nsub_bins_vol4cut, integral_ttbar_nsub_vol4cut, integral_dijets_nsub_vol4cut);
  ROC_nsub_vol4cut->GetXaxis()->SetLimits(0, 1);
  ROC_nsub_vol4cut->GetYaxis()->SetLimits(0, 1);
  ROC_nsub_vol4cut->SetLineColor(kBlack);
  // ROC_nsub_vol4cut->SetLineWidth(4);
  ROC_nsub_vol4cut->SetMarkerStyle(5);
  ROC_nsub_vol4cut->SetMarkerSize(4);
  ROC_nsub_vol4cut->Write();

  int nsub_bins_vol6cut = tau32_wta_dijets_vol6cut->GetSize() - 2;
  double integral_dijets_nsub_vol6cut[nsub_bins_vol6cut], integral_ttbar_nsub_vol6cut[nsub_bins_vol6cut];
  for (int i = 0; i < nsub_bins_vol6cut; i++) {
    integral_dijets_nsub_vol6cut[i] = tau32_wta_dijets_vol6cut->Integral(0, i)/total_dijets_jets;
    integral_ttbar_nsub_vol6cut[i] = tau32_wta_ttbar_vol6cut->Integral(0, i)/total_ttbar_jets;
  }
  TGraph* ROC_nsub_vol6cut = new TGraph(nsub_bins_vol6cut, integral_ttbar_nsub_vol6cut, integral_dijets_nsub_vol6cut);
  ROC_nsub_vol6cut->GetXaxis()->SetLimits(0, 1);
  ROC_nsub_vol6cut->GetYaxis()->SetLimits(0, 1);
  ROC_nsub_vol6cut->SetLineColor(kBlack);
  // ROC_nsub_vol6cut->SetLineWidth(6);
  ROC_nsub_vol6cut->SetMarkerStyle(5);
  ROC_nsub_vol6cut->SetMarkerSize(6);
  ROC_nsub_vol6cut->Write();

  int nsub_bins_vol8cut = tau32_wta_dijets_vol8cut->GetSize() - 2;
  double integral_dijets_nsub_vol8cut[nsub_bins_vol8cut], integral_ttbar_nsub_vol8cut[nsub_bins_vol8cut];
  for (int i = 0; i < nsub_bins_vol8cut; i++) {
    integral_dijets_nsub_vol8cut[i] = tau32_wta_dijets_vol8cut->Integral(0, i)/total_dijets_jets;
    integral_ttbar_nsub_vol8cut[i] = tau32_wta_ttbar_vol8cut->Integral(0, i)/total_ttbar_jets;
  }
  TGraph* ROC_nsub_vol8cut = new TGraph(nsub_bins_vol8cut, integral_ttbar_nsub_vol8cut, integral_dijets_nsub_vol8cut);
  ROC_nsub_vol8cut->GetXaxis()->SetLimits(0, 1);
  ROC_nsub_vol8cut->GetYaxis()->SetLimits(0, 1);
  ROC_nsub_vol8cut->SetLineColor(kBlack);
  // ROC_nsub_vol8cut->SetLineWidth(8);
  ROC_nsub_vol8cut->SetMarkerStyle(5);
  ROC_nsub_vol8cut->SetMarkerSize(8);
  ROC_nsub_vol8cut->Write();

  int nsub_bins_vol10cut = tau32_wta_dijets_vol10cut->GetSize() - 2;
  double integral_dijets_nsub_vol10cut[nsub_bins_vol10cut], integral_ttbar_nsub_vol10cut[nsub_bins_vol10cut];
  for (int i = 0; i < nsub_bins_vol10cut; i++) {
    integral_dijets_nsub_vol10cut[i] = tau32_wta_dijets_vol10cut->Integral(0, i)/total_dijets_jets;
    integral_ttbar_nsub_vol10cut[i] = tau32_wta_ttbar_vol10cut->Integral(0, i)/total_ttbar_jets;
  }
  TGraph* ROC_nsub_vol10cut = new TGraph(nsub_bins_vol10cut, integral_ttbar_nsub_vol10cut, integral_dijets_nsub_vol10cut);
  ROC_nsub_vol10cut->GetXaxis()->SetLimits(0, 1);
  ROC_nsub_vol10cut->GetYaxis()->SetLimits(0, 1);
  ROC_nsub_vol10cut->SetLineColor(kBlack);
  // ROC_nsub_vol10cut->SetLineWidth(10);
  ROC_nsub_vol10cut->SetMarkerStyle(5);
  ROC_nsub_vol10cut->SetMarkerSize(10);
  ROC_nsub_vol10cut->Write();

  TCanvas *ROC_compare = new TCanvas("ROC_compare", "ROC Curve Comparison", 600, 600);
  ROC_compare->cd();
  ROC_compare->SetLogy();
  // ROC_vol_nocut->SetLineColor(1);
  // ROC_nsub_nocut->SetLineColor(1);
  // ROC_nsub_nocut->SetLineStyle(2);
  // ROC_vol_masscut->SetLineColor(2);
  // ROC_nsub_masscut->SetLineColor(kGreen);
  // ROC_nsub_masscut->SetLineStyle(2);
  // ROC_vol_nsubcut->SetLineColor(4);
  // ROC_nsub_vol2cut->SetLineColor(kBlack);
  ROC_nsub_masscut->SetLineStyle(2);
  ROC_nsub_vol2cut->SetLineColor(kRed);
  ROC_nsub_vol4cut->SetLineColor(kBlue);
  ROC_nsub_vol6cut->SetLineColor(kGreen);
  ROC_nsub_vol8cut->SetLineColor(kYellow);
  ROC_nsub_vol10cut->SetLineColor(kRed);
  ROC_nsub_vol10cut->SetLineStyle(2);
  TMultiGraph *ROC_multigraph_beta_compare = new TMultiGraph("ROC_multigraph_beta_compare", "ROC comparison for volatility vs. #tau_{3}/#tau_{2}");
  // ROC_multigraph_beta_compare->Add(ROC_vol_nocut);
  ROC_multigraph_beta_compare->Add(ROC_nsub_nocut);
  // ROC_multigraph_beta_compare->Add(ROC_vol_masscut);
  ROC_multigraph_beta_compare->Add(ROC_nsub_masscut);
  // ROC_multigraph_beta_compare->Add(ROC_vol_nsubcut);
  ROC_multigraph_beta_compare->Add(ROC_nsub_vol2cut);
  ROC_multigraph_beta_compare->Add(ROC_nsub_vol4cut);
  ROC_multigraph_beta_compare->Add(ROC_nsub_vol6cut);
  ROC_multigraph_beta_compare->Add(ROC_nsub_vol8cut);
  // ROC_multigraph_beta_compare->Add(ROC_nsub_vol10cut);
  ROC_multigraph_beta_compare->Draw("AL");

  TLegend *leg_ROC_compare = new TLegend(0.1, 0.7, 0.3, 0.9);
  // leg_ROC_compare->AddEntry(ROC_vol_nocut, "Volatility No cut", "L");
  leg_ROC_compare->AddEntry(ROC_nsub_nocut, "#tau_{3}/#tau_{2} No cut", "L");
  // leg_ROC_compare->AddEntry(ROC_nsub_masscut, "#tau_{3}/#tau_{2} 160 < m < 240, 60 < m_{W} < 100", "L");
  // leg_ROC_compare->AddEntry(ROC_nsub_vol2cut, "#tau_{3}/#tau_{2} 160 < m < 240, 60 < m_{W} < 100 + volatility < 0.2", "L");
  // leg_ROC_compare->AddEntry(ROC_nsub_vol4cut, "#tau_{3}/#tau_{2} 160 < m < 240, 60 < m_{W} < 100 + volatility < 0.4", "L");
  // leg_ROC_compare->AddEntry(ROC_nsub_vol6cut, "#tau_{3}/#tau_{2} 160 < m < 240, 60 < m_{W} < 100 + volatility < 0.6", "L");
  // leg_ROC_compare->AddEntry(ROC_nsub_vol8cut, "#tau_{3}/#tau_{2} 160 < m < 240, 60 < m_{W} < 100 + volatility < 0.8", "L");
  // leg_ROC_compare->AddEntry(ROC_nsub_vol2cut, "#tau_{3}/#tau_{2} 160 < m < 240, 60 < m_{W} < 100 + volatility < 1.0", "L");
  leg_ROC_compare->AddEntry(ROC_nsub_masscut, "#tau_{3}/#tau_{2} 160 < m < 240", "L");
  leg_ROC_compare->AddEntry(ROC_nsub_vol2cut, "#tau_{3}/#tau_{2} 160 < m < 240 + volatility < 0.2", "L");
  leg_ROC_compare->AddEntry(ROC_nsub_vol4cut, "#tau_{3}/#tau_{2} 160 < m < 240 + volatility < 0.4", "L");
  leg_ROC_compare->AddEntry(ROC_nsub_vol6cut, "#tau_{3}/#tau_{2} 160 < m < 240 + volatility < 0.6", "L");
  leg_ROC_compare->AddEntry(ROC_nsub_vol8cut, "#tau_{3}/#tau_{2} 160 < m < 240 + volatility < 0.8", "L");
  // leg_ROC_compare->AddEntry(ROC_nsub_vol2cut, "#tau_{3}/#tau_{2} 160 < m < 240 + volatility < 1.0", "L");
  // leg_ROC_compare->AddEntry(ROC_vol_masscut, "Volatility 160 < m < 180, 70 < m_{W} < 90", "L");
  // leg_ROC_compare->AddEntry(ROC_vol_nsubcut, "Volatility 160 < m 180, 70 < m_{W} < 90 + #tau_{3}/#tau_{2} < 0.6", "L"); 
  leg_ROC_compare->Draw();
  ROC_compare->Write();
  // ROC_compare->Print("ROC_compare_800to800_log.pdf", "pdf");

  out.Close();
  // delete volatility_hist_dijets;
  // delete volatility_hist_ttbar;

  return 0;
}