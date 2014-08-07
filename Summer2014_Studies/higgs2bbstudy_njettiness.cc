//Misc. Headers
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

//Pythia headers
#include "Pythia.h"

//Fastjet headers
#include "FastJet3.h"
#include "fastjet/Selector.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"

//Root headers
#include "TH1.h"
#include "TH2.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"

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


// Simple class to store Axes along with a name for display
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

double calcDphi(PseudoJet jet1, PseudoJet jet2) {
  double dphi;
  dphi = jet1.phi() - jet2.phi();
  if (dphi >  TMath::Pi()) dphi -= 2*TMath::Pi();
  if (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
  return dphi; 
}

Double_t calcMedian1(TH1 *h1) { 
   //compute the median for 1-d histogram h1 
 Int_t nbins = h1->GetXaxis()->GetNbins(); 
 Double_t *x = new Double_t[nbins]; 
 Double_t *y = new Double_t[nbins]; 
 for (Int_t i=0;i<nbins;i++) {
  x[i] = h1->GetXaxis()->GetBinCenter(i+1); 
  y[i] = h1->GetBinContent(i+1); 
} 
Double_t median = TMath::Median(nbins,x,y); 
delete [] x; 
delete [] y; 
return median; 
} 

int main(int argc, char* argv[]) {

  double epsilon = 0.0001;
  int nEvents = 100;
  gStyle->SetOptStat(0);
  // Create the ROOT application environment.
  TApplication theApp("hist", &argc, argv);

  // Create file on which histograms can be saved.  
  TFile* outFile = new TFile("higgs2bbstudy_allpt.root", "RECREATE");

  int n_perps = 11;

  //create list of various values of beta
  vector<double> betalist;
  betalist.push_back(0.25);
  betalist.push_back(0.5);
  betalist.push_back(1.0);
  betalist.push_back(2.0);
  // betalist.push_back(3.0);
  int n_betas = betalist.size();

  double perplist[n_perps];
  double higgs_efficiencies_1jet[n_betas + 1][n_perps];
  double higgs_efficiencies_2jets[n_betas + 1][n_perps];
  double higgs_efficiencies_3jets[n_betas + 1][n_perps];
  double higgs_efficiencies_32ratio[n_betas + 1][n_perps];
  double angulardiff_2jets[n_betas + 1][n_perps];
  double angulardiff_width_2jets[n_betas + 1][n_perps];
  double angulardiff_1jet[n_betas + 1][n_perps];
  double angulardiff_prop_1jet[n_betas + 1][n_perps];
  double jet_distance_1jet_median[n_betas + 1][n_perps];
  double jet_distance_1jet_firstquartile[n_betas + 1][n_perps];
  double jet_distance_1jet_thirdquartile[n_betas + 1][n_perps];

  for (int i_perp = 0; i_perp < n_perps; i_perp++) {

  // double ptcut = atoi(argv[1]);
    double ptcut = i_perp*100.0;
  // perplist.push_back(ptcut);
    perplist[i_perp] = ptcut;
    string pythia_ptcut;
    ostringstream convert;
    convert << ptcut;
    pythia_ptcut = convert.str();

  // Generator. Process selection. LHC upgrade initialization. 
    Pythia pythia;
  pythia.readString("Beams:eCM = 14000.");    //LHC Upgrade COM Energy
  pythia.readString("HiggsSM:ffbar2HZ = on"); //Turn on Higgs + Z
  pythia.readString("23:onMode = off"); 
  pythia.readString("23:onIfAny = 12 14 16"); //Only turn on Z -> neutrino decays
  pythia.readString("25:onMode = off");  
  pythia.readString("25:onIfAny = 5 -5"); // only turn on Higgs -> bbbar
  pythia.readString("PhaseSpace:pTHatMin = " + pythia_ptcut);
  pythia.init();

  TH1* particle_mass = new TH1F("particle_mass", "Mass of Higgs particle", 100, 0, 500);

  TH1* akt_jet_mass = new TH1F("akt_jet_mass", "Mass of individual akt jet", 100, 0, 500);
  TH1* akt_jet_perp = new TH1F("akt_jet_perp", "Perp of individual akt jet", 100, 0, 500);
  TH1* akt_jet_phi = new TH1F("akt_jet_phi", "Phi of individual akt jet", 32, 0, 6.4);
  TH1* akt_jet_eta = new TH1F("akt_jet_eta", "Eta of individual akt jet", 50, -5, 5);

  TH1* higgs_invmass_akt = new TH1F("higgs_invmass_akt", "Invariant Mass of Higgs (antikt)", 100, 0, 500);
  TH1* jet_ptspectrum_akt = new TH1F("jet_ptspectrum_akt", "pT spectrum of jets (antikt)", 100, 0, 500);
  TH1* jet_phispectrum_akt = new TH1F("jet_phispectrum_akt", "#phi spectrum of jets (antikt)", 32, 0, 3.2);
  TH1* jet_distance_akt = new TH1F("jet_distance_akt", "Distance between two jets (antikt)", 100, 0, 10.0);
  TH1* jet_distance_diff_akt = new TH1F("jet_distance_diff_akt", "Distance between individual jet and closest of two jets", 60, 0, 3.0);

  TH1* threejet_invmass_akt = new TH1F("threejet_invmass_akt", "3-jet invariant mass (antikt)", 100, 0, 500);
  TH1* threejet_thirdjet_mass_akt = new TH1F("threejet_thirdjet_mass_akt", "Third jet mass (antikt)", 100, 0, 500);
  TH1* threejet_thirdjet_perp_akt = new TH1F("threejet_thirdjet_perp_akt", "Third jet perp (antikt)", 100, 0, 500);

  TObjArray njet_jet_mass_hists(n_betas);  
  TObjArray njet_jet_perp_hists(n_betas);
  TObjArray njet_jet_phi_hists(n_betas);
  TObjArray njet_jet_eta_hists(n_betas);
  TObjArray jetmass_diff_hists(n_betas);
  TObjArray jetperp_diff_hists(n_betas);
  TObjArray jetphi_diff_hists(n_betas);
  TObjArray tau2_diff_hists(n_betas);
  TObjArray tau1_diff_hists(n_betas);

  TObjArray higgs_invmass_njet_hists(n_betas);
  TObjArray jet_ptspectrum_njet_hists(n_betas);
  TObjArray jet_phispectrum_njet_hists(n_betas);
  TObjArray jet_distance_njet_hists(n_betas);
  TObjArray jet_distance_diff_njet_hists(n_betas);

  TObjArray threejet_invmass_njet_hists(n_betas);
  TObjArray threejet_thirdjet_mass_njet_hists(n_betas);
  TObjArray threejet_thirdjet_perp_njet_hists(n_betas);

  for (int Bs = 0; Bs < n_betas; Bs++) {

    double beta = betalist[Bs];

    ostringstream ss;
    ss << beta;
    TString title;

    TH1* njet_jet_mass = new TH1F("njet_jet_mass", "Mass of individual njet jet (#beta = " + (TString)ss.str() + ")", 100, 0, 500);
    TH1* njet_jet_perp = new TH1F("njet_jet_perp", "Perp of individual njet jet (#beta = " + (TString)ss.str() + ")", 100, 0, 500);
    TH1* njet_jet_phi = new TH1F("njet_jet_phi", "Phi of individual njet jet (#beta = " + (TString)ss.str() + ")", 32, 0, 6.4);
    TH1* njet_jet_eta = new TH1F("njet_jet_eta", "Eta of individual njet jet (#beta = " + (TString)ss.str() + ")", 50, -5, 5);
    TH1* jetmass_diff = new TH1F("jetmass_diff", "difference in mass of akt and njet jets (#beta = " + (TString)ss.str() + ")", 100, 0, 500);
    TH1* jetperp_diff = new TH1F("jetperp_diff", "difference in perp of akt and njet jets (#beta = " + (TString)ss.str() + ")", 100, -200, 200);
    TH1* jetphi_diff = new TH1F("jetphi_diff", "difference in phi of akt and njet jets (#beta = " + (TString)ss.str() + ")", 100, 0, 1.0);
    TH1* tau1_diff = new TH1F("tau1_diff", "difference in #tau_{1} jets (akt - njet) (#beta = " + (TString)ss.str() + ")", 100, -200, 200);
    TH1* tau2_diff = new TH1F("tau2_diff", "difference in #tau_{2} jets (akt - njet) (#beta = " + (TString)ss.str() + ")", 100, -200, 200);

    TH1* higgs_invmass_njet = new TH1F("higgs_invmass_njet", "Invariant Mass of Higgs (njet, beta = " + (TString)ss.str() + ")", 100, 0, 500);
    TH1* jet_ptspectrum_njet = new TH1F("jet_ptspectrum_njet", "pT spectrum of jets (njet, beta = " + (TString)ss.str() + ")", 100, 0, 500);
    TH1* jet_phispectrum_njet = new TH1F("jet_phispectrum_njet", "#phi spectrum of jets (njet, beta = " + (TString)ss.str() + ")", 32, 0, 3.2);
    TH1* jet_distance_njet = new TH1F("jet_distance_njet", "Distance between two jets (njet, beta = " + (TString)ss.str() + ")", 100, 0, 10.0);
    TH1* jet_distance_diff_njet = new TH1F("jet_distance_diff_njet", "Distance between individual jet and closest of two jets (njet, beta = " + (TString)ss.str() + ")", 100, 0, 0.5);

    TH1* threejet_invmass_njet = new TH1F("threejet_invmass_njet", "3-jet invariant mass (njet, beta = " + (TString)ss.str() + ")", 100, 0, 500);
    TH1* threejet_thirdjet_mass_njet = new TH1F("threejet_thirdjet_mass_njet", "Third jet mass (njet, beta = " + (TString)ss.str() + ")", 100, 0, 500);
    TH1* threejet_thirdjet_perp_njet = new TH1F("threejet_thirdjet_perp_njet", "Third jet perp (njet, beta = " + (TString)ss.str() + ")", 100, 0, 500);

    njet_jet_mass_hists.Add(njet_jet_mass);
    njet_jet_perp_hists.Add(njet_jet_perp);
    njet_jet_phi_hists.Add(njet_jet_phi);
    njet_jet_eta_hists.Add(njet_jet_eta);
    jetmass_diff_hists.Add(jetmass_diff);
    jetperp_diff_hists.Add(jetperp_diff);
    jetphi_diff_hists.Add(jetphi_diff);
    tau2_diff_hists.Add(tau2_diff);
    tau1_diff_hists.Add(tau1_diff);

    higgs_invmass_njet_hists.Add(higgs_invmass_njet);
    jet_ptspectrum_njet_hists.Add(jet_ptspectrum_njet);
    jet_phispectrum_njet_hists.Add(jet_phispectrum_njet);
    jet_distance_njet_hists.Add(jet_distance_njet);
    jet_distance_diff_njet_hists.Add(jet_distance_diff_njet);

    threejet_invmass_njet_hists.Add(threejet_invmass_njet);
    threejet_thirdjet_mass_njet_hists.Add(threejet_thirdjet_mass_njet);
    threejet_thirdjet_perp_njet_hists.Add(threejet_thirdjet_perp_njet);
  }

  // Fastjet input (different inputs for quark vs gluon)
  vector <PseudoJet> fjInputs;

  double Rparam = 0.5;
  Strategy strategy = Best;
  RecombinationScheme recombScheme_akt = E_scheme;
  JetDefinition *jetDef = new JetDefinition(antikt_algorithm, Rparam, recombScheme_akt, strategy);

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

    if (!pythia.next()) continue;

    // Reset Fastjet input
    fjInputs.resize(0);

    int iH = 0;

    //store particle as input for Fastjet
    for (int i = 0; i < pythia.event.size(); ++i) {
      // Locate the final copy of the Higgs.
      if (pythia.event[i].id() == 25) iH = i;

      // Final State only
      if (!pythia.event[i].isFinal()) continue;

      //Exclude neutrinos and leptons, which should only appear as final products of bosons, which should not be in the jets
      // if (pythia.event[pythia.event[i].mother1()].id() == 23) continue;
      if (pythia.event[i].idAbs() == 12 || pythia.event[i].idAbs() == 14 || pythia.event[i].idAbs() == 16) continue;
      // if (pythia.event[i].idAbs() == 11 || pythia.event[i].idAbs() == 13 || pythia.event[i].idAbs() == 15) continue;
      PseudoJet fj_particle = pythia.event[i];
      fjInputs.push_back(fj_particle);
    }

    TH2F *event_display = new TH2F("event_display", "Event Display", 60, -5, 5, 60, 0, 6.4);
  
    particle_mass->Fill(pythia.event[iH].m());

    for (int i_part = 0; i_part < fjInputs.size(); i_part++) {
      event_display->Fill(fjInputs[i_part].eta(), fjInputs[i_part].phi(), fjInputs[i_part].perp());
    }

    //Run Fastjet algorithm
    vector <PseudoJet> inclusiveJets, sortedJets, centralJets;
    ClusterSequence clustSeq(fjInputs, *jetDef);

    //Extract inclusive jets sorted by pT
    inclusiveJets = clustSeq.inclusive_jets(0.0);

    // Sort central jets by pT
    sortedJets = sorted_by_pt(inclusiveJets);

    //Cut on jets within eta range
    Selector eta_selector = SelectorAbsEtaMax(5.0);
    centralJets = eta_selector(sortedJets);

    Selector jet_selector = SelectorNHardest(1);
    vector<PseudoJet> hardest_jet = jet_selector(centralJets);

    if (hardest_jet.size() == 1) {
      akt_jet_mass->Fill(hardest_jet[0].m());
      akt_jet_perp->Fill(hardest_jet[0].perp());
      akt_jet_phi->Fill(hardest_jet[0].phi());
      akt_jet_eta->Fill(hardest_jet[0].eta());
    }

    //Create selectors for hardest jet in the event
    Selector twojet_selector = SelectorNHardest(2);
    vector<PseudoJet> hardest_twojets = twojet_selector(centralJets);

    PseudoJet big_jet(0,0,0,0);
    for (int i_jet = 0; i_jet < hardest_twojets.size(); i_jet++) {
      big_jet = join(big_jet, hardest_twojets[i_jet]);
      jet_ptspectrum_akt->Fill(hardest_twojets[i_jet].perp());
    }
    if (hardest_twojets.size() == 2) {
      jet_phispectrum_akt->Fill(abs(calcDphi(hardest_twojets[0], hardest_twojets[1])));
      jet_distance_akt->Fill(hardest_twojets[0].delta_R(hardest_twojets[1]));
      if (hardest_jet.size() == 1) {
        double smaller_distance = (hardest_jet[0].delta_R(hardest_twojets[0]) < hardest_jet[0].delta_R(hardest_twojets[1])) ? hardest_jet[0].delta_R(hardest_twojets[0]) : hardest_jet[0].delta_R(hardest_twojets[1]);
        jet_distance_diff_akt->Fill(smaller_distance);
      }
    }
    higgs_invmass_akt->Fill(big_jet.m());

    Selector threejet_selector = SelectorNHardest(3);
    vector<PseudoJet> hardest_threejets = threejet_selector(centralJets);

    PseudoJet third_jet = hardest_threejets[hardest_threejets.size() - 1];
    threejet_thirdjet_mass_akt->Fill(third_jet.m());
    threejet_thirdjet_perp_akt->Fill(third_jet.perp());

    PseudoJet big_jet_3jets(0,0,0,0);
    // for (int i_jet = 0; i_jet < hardest_threejets.size(); i_jet++) {
    if ((third_jet.delta_R(hardest_twojets[0]) < Rparam || third_jet.delta_R(hardest_twojets[1]) < Rparam) && third_jet.perp() > 50) {
      big_jet_3jets = join(big_jet, third_jet);
    }
    else big_jet_3jets = big_jet;
    // }
    threejet_invmass_akt->Fill(big_jet_3jets.m());


    for (int B = 0; B < n_betas; B++) {

      TH2F *axes_display = new TH2F("axes_display", "Axes Plot", 300, -5, 5, 300, 0, 6.4);
      TH2F *axes_njet_display = new TH2F("axes_njet_display", "Axes Plot (Njettiness)", 300, -5, 5, 300, 0, 6.4);
      TH2F *akt_jets_display = new TH2F("akt_jets_display", "Anti-KT Jets", 300, -5, 5, 300, 0, 6.4);
      TH2F *njet_jets_display = new TH2F("njet_jets_display", "N-jettiness Jets", 300, -5, 5, 300, 0, 6.4);

      double beta = betalist[B];
      double power = (double)1/beta;
      // double delta;

      const JetDefinition::Recombiner *recombScheme;
      if (beta > 1) recombScheme = new GeneralERecombiner((double)1/(beta - 1));
      else recombScheme = new WinnerTakeAllRecombiner();

      UnnormalizedCutoffMeasure measure_function = UnnormalizedCutoffMeasure(beta, Rparam);
      // OnePass_Manual_Axes axes_finder = OnePass_Manual_Axes();

      AxesStruct *axes_finder;
      if (beta < 1 || beta > 3) {
        axes_finder = new AxesStruct(Manual_Axes());
      }
      else {
        axes_finder = new AxesStruct(OnePass_Manual_Axes());
      }

      JetDefinition *jetDef = new JetDefinition(genkt_algorithm, Rparam, power, recombScheme, strategy);
      ClusterSequence clustSeq(fjInputs, *jetDef);
      vector<PseudoJet> exclusive_1jet_start = clustSeq.exclusive_jets(2);
      vector<PseudoJet> exclusive_1jet = findMinAxes(fjInputs, exclusive_1jet_start, 1, beta, Rparam);
      // vector<PseudoJet> exclusive_1jet = findMinAxes(fjInputs, exclusive_1jet_start, 1, beta, Rparam);

      //Use Njettiness algorithm for comparison
      NjettinessPlugin njet_plugin_1jet(1, axes_finder->def(), measure_function);
      JetDefinition njet_def_1jet(&njet_plugin_1jet);
      njet_plugin_1jet.setAxes(exclusive_1jet);
      ClusterSequence njet_cluster_1jet(fjInputs, njet_def_1jet);
      const NjettinessExtras *extras_1jet = njettiness_extras(njet_cluster_1jet);
      vector<PseudoJet> njet_1jet = extras_1jet->jets();

      NjettinessPlugin njet_plugin_manual_1jet(1, OnePass_Manual_Axes(), measure_function);
      JetDefinition njet_def_manual_1jet(&njet_plugin_manual_1jet);

      if (njet_1jet.size() == 1) {

        TH1* njet_jet_mass = (TH1*)njet_jet_mass_hists.At(B);
        TH1* njet_jet_perp = (TH1*)njet_jet_perp_hists.At(B);
        TH1* njet_jet_phi = (TH1*)njet_jet_phi_hists.At(B);
        TH1* njet_jet_eta = (TH1*)njet_jet_eta_hists.At(B);

        njet_jet_mass->Fill(njet_1jet[0].m());
        njet_jet_perp->Fill(njet_1jet[0].perp());
        njet_jet_phi->Fill(njet_1jet[0].phi());
        njet_jet_eta->Fill(njet_1jet[0].eta());

        if (hardest_jet.size() == 1) {
          njet_plugin_manual_1jet.setAxes(hardest_jet);
          ClusterSequence njet_cluster_manual_1jet(fjInputs, njet_def_manual_1jet);
          const NjettinessExtras *extras_manual_1jet = njettiness_extras(njet_cluster_manual_1jet);

          TH1* tau1_diff = (TH1*)tau1_diff_hists.At(B);
          tau1_diff->Fill(extras_manual_1jet->totalTau() - extras_1jet->totalTau());

          TH1* jetmass_diff = (TH1*)jetmass_diff_hists.At(B);
          TH1* jetperp_diff = (TH1*)jetperp_diff_hists.At(B);
          TH1* jetphi_diff = (TH1*)jetphi_diff_hists.At(B);

          jetmass_diff->Fill(abs(hardest_jet[0].m() - njet_1jet[0].m()));
          jetperp_diff->Fill(hardest_jet[0].perp() - njet_1jet[0].perp());
          jetphi_diff->Fill(njet_1jet[0].delta_R(hardest_jet[0]));
        }
      }

      vector<PseudoJet> exclusive_2jets_start = clustSeq.exclusive_jets(2);
      vector<PseudoJet> exclusive_2jets = findMinAxes(fjInputs, exclusive_2jets_start, 2, beta, Rparam);
      // vector<PseudoJet> exclusive_2jets = clustSeq.exclusive_jets(2);
      NjettinessPlugin njet_plugin_twojet(2, axes_finder->def(), measure_function);
      JetDefinition njet_def_twojet(&njet_plugin_twojet);
      njet_plugin_twojet.setAxes(exclusive_2jets);
      ClusterSequence njet_cluster_twojet(fjInputs, njet_def_twojet);
      const NjettinessExtras *extras_twojet = njettiness_extras(njet_cluster_twojet);
      vector<PseudoJet> njet_jets_2jets = extras_twojet->jets();

      PseudoJet big_jet_njet(0,0,0,0);
      TH1* higgs_invmass_njet = (TH1*)higgs_invmass_njet_hists.At(B);
      TH1* jet_ptspectrum_njet = (TH1*)jet_ptspectrum_njet_hists.At(B);
      TH1* jet_phispectrum_njet = (TH1*)jet_phispectrum_njet_hists.At(B);
      TH1* jet_distance_njet = (TH1*)jet_distance_njet_hists.At(B);
      TH1* jet_distance_diff_njet = (TH1*)jet_distance_diff_njet_hists.At(B);

      for (int i_jet = 0; i_jet < njet_jets_2jets.size(); i_jet++) {
        big_jet_njet = join(big_jet_njet, njet_jets_2jets[i_jet]);
        jet_ptspectrum_njet->Fill(njet_jets_2jets[i_jet].perp());
      }
      jet_distance_njet->Fill(njet_jets_2jets[0].delta_R(njet_jets_2jets[1]));
      jet_phispectrum_njet->Fill(abs(calcDphi(njet_jets_2jets[0], njet_jets_2jets[1])));
      double smaller_distance_njet = (njet_1jet[0].delta_R(njet_jets_2jets[0]) < njet_1jet[0].delta_R(njet_jets_2jets[1])) ? njet_1jet[0].delta_R(njet_jets_2jets[0]) : njet_1jet[0].delta_R(njet_jets_2jets[1]);
      if (smaller_distance_njet > epsilon) jet_distance_diff_njet->Fill(smaller_distance_njet);
      higgs_invmass_njet->Fill(big_jet_njet.m());

      vector<PseudoJet> exclusive_threejets_start = clustSeq.exclusive_jets(3);
      vector<PseudoJet> exclusive_threejets = findMinAxes(fjInputs, exclusive_threejets_start, 3, beta, Rparam);
      // vector<PseudoJet> exclusive_threejets = clustSeq.exclusive_jets(2);
      NjettinessPlugin njet_plugin_threejet(3, axes_finder->def(), measure_function);
      JetDefinition njet_def_threejet(&njet_plugin_threejet);
      njet_plugin_threejet.setAxes(exclusive_threejets);
      ClusterSequence njet_cluster_threejet(fjInputs, njet_def_threejet);
      const NjettinessExtras *extras_threejet = njettiness_extras(njet_cluster_threejet);
      vector<PseudoJet> njet_jets_threejets = extras_threejet->jets();

      vector<PseudoJet> sorted_njet_jets_threejets = sorted_by_pt(njet_jets_threejets);
      PseudoJet third_jet_njet = sorted_njet_jets_threejets[sorted_njet_jets_threejets.size() - 1];


      PseudoJet big_jet_njet_3jets(0,0,0,0);
      // for (int i = 0; i < njet_jets_threejets.size(); i++) {
      if ((third_jet_njet.delta_R(njet_jets_2jets[0]) < Rparam || third_jet_njet.delta_R(njet_jets_2jets[1]) < Rparam) && third_jet_njet.perp() > 50) {
        big_jet_njet_3jets = join(big_jet_njet, third_jet_njet);
      }
      else big_jet_njet_3jets = big_jet_njet;
        // if (njet_jets_threejets[0].delta_R(njet_jets_threejets[1]) < 1.0 && njet_jets_threejets[1].delta_R(njet_jets_threejets[2]) < 1.0
        //   && njet_jets_threejets[0].delta_R(njet_jets_threejets[2]) < 1.0)
        //   big_jet_njet_3jets = join(big_jet_njet_3jets, njet_jets_threejets[i]);
      // }

      TH1* threejet_invmass_njet = (TH1*)threejet_invmass_njet_hists.At(B);
      TH1* threejet_thirdjet_mass_njet = (TH1*)threejet_thirdjet_mass_njet_hists.At(B);
      TH1* threejet_thirdjet_perp_njet = (TH1*)threejet_thirdjet_perp_njet_hists.At(B);

      threejet_invmass_njet->Fill(big_jet_njet_3jets.m());
      threejet_thirdjet_mass_njet->Fill(third_jet_njet.m());
      threejet_thirdjet_perp_njet->Fill(third_jet_njet.perp());
    
      if ((B == 0 || B == 1) && big_jet_njet.m() < 75) {
      event_display->SetStats(0);
      axes_display->SetStats(0);
      axes_njet_display->SetStats(0);
      akt_jets_display->SetStats(0);
      njet_jets_display->SetStats(0);

      TCanvas *display = new TCanvas("display", "Event Display", 1093, 700);
      display->cd();
      display->SetFixedAspectRatio();
      event_display->GetXaxis()->SetTitle("#eta");
      event_display->GetYaxis()->SetTitle("#phi");
      event_display->SetFillColor(kBlack);
      event_display->SetLineColor(kBlack);
      event_display->SetLineWidth(1);
      event_display->Draw("box");

      for (int a = 0; a < hardest_twojets.size(); a++) {
        axes_display->Fill(hardest_twojets[a].eta(), hardest_twojets[a].phi());
        vector<PseudoJet> constituents = hardest_twojets[a].constituents();
        for (int i_const = 0; i_const < constituents.size(); i_const++) {
          akt_jets_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
        }          
      }
      axes_display->SetMarkerStyle(3);
      axes_display->SetMarkerSize(3);
      axes_display->SetMarkerColor(kBlue);
      akt_jets_display->SetMarkerStyle(21);
      akt_jets_display->SetMarkerSize(0.5);
      akt_jets_display->SetMarkerColor(kBlue);
      axes_display->Draw("SAMES");
      akt_jets_display->Draw("SAMES");

      for (int a = 0; a < njet_jets_2jets.size(); a++) {
        axes_njet_display->Fill(njet_jets_2jets[a].eta(), njet_jets_2jets[a].phi());
        vector<PseudoJet> constituents = njet_jets_2jets[a].constituents();
        for (int i_const = 0; i_const < constituents.size(); i_const++) {
          njet_jets_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
        }          
      }
      axes_njet_display->SetMarkerStyle(3);
      axes_njet_display->SetMarkerSize(3);
      axes_njet_display->SetMarkerColor(kRed);
      njet_jets_display->SetMarkerStyle(21);
      njet_jets_display->SetMarkerSize(0.5);
      njet_jets_display->SetMarkerColor(kRed);
      axes_njet_display->Draw("SAMES");
      njet_jets_display->Draw("SAMES");

      display->Write();
      }

      // delete display;
      delete axes_display;
      delete axes_njet_display;
      delete njet_jets_display;
      delete akt_jets_display;
    }
    delete event_display;

  }

  // pythia.stat();
  particle_mass->SetStats(0);
  particle_mass->Write();

  higgs_invmass_akt->Write();
  jet_ptspectrum_akt->Write();
  jet_phispectrum_akt->Write();
  jet_distance_akt->Write();
  jet_distance_diff_akt->Write();

  TCanvas *mass_compare_2jets = new TCanvas("mass_compare_2jets", "mass_compare_2jets", 600, 600);
  mass_compare_2jets->cd();

  double higgs_invmass_akt_scale = 1/higgs_invmass_akt->Integral(0, 101);
  higgs_invmass_akt->Scale(higgs_invmass_akt_scale);
  higgs_invmass_akt->SetLineColor(kBlack);
  higgs_invmass_akt->SetTitle("2-jet Mass (pT > " + (TString)pythia_ptcut + ", R_{0} = 0.5)");
  higgs_invmass_akt->GetXaxis()->SetTitle("m_{jj} [GeV]");
  higgs_invmass_akt->SetMinimum(0.0);
  double max_val = higgs_invmass_akt->GetMaximum();

  higgs_invmass_akt->Draw();
  TLegend *leg_mass_2jets = new TLegend(0.48, 0.48, 0.88, 0.88);
  leg_mass_2jets->AddEntry(higgs_invmass_akt, "Anti-KT", "L");
  leg_mass_2jets->SetFillColor(kWhite);
  leg_mass_2jets->SetLineColor(kWhite);
  double akt_efficiency = higgs_invmass_akt->Integral(20, 30);
  higgs_efficiencies_2jets[0][i_perp] = akt_efficiency;

  // double akt_angulardiff = jet_phispectrum_akt->GetMean();
  // double akt_angulardiff = calcMedian1(jet_phispectrum_akt);
  double akt_angulardiff = calcMedian1(jet_distance_akt);
  double akt_angulardiff_width = jet_phispectrum_akt->GetRMS();
  angulardiff_2jets[0][i_perp] = akt_angulardiff;
  angulardiff_width_2jets[0][i_perp] = akt_angulardiff_width;
  angulardiff_1jet[0][i_perp] = jet_distance_diff_akt->GetMean();

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1* njet_jet_perp = (TH1*)njet_jet_perp_hists.At(B);
    TH1* njet_jet_phi = (TH1*)njet_jet_phi_hists.At(B);
    TH1* njet_jet_eta = (TH1*)njet_jet_eta_hists.At(B);

    TH1* jetmass_diff = (TH1*)jetmass_diff_hists.At(B);
    TH1* jetperp_diff = (TH1*)jetperp_diff_hists.At(B);
    TH1* jetphi_diff = (TH1*)jetphi_diff_hists.At(B);
    TH1* tau2_diff = (TH1*)tau2_diff_hists.At(B);
    TH1* tau1_diff = (TH1*)tau1_diff_hists.At(B);

    TH1* higgs_invmass_njet = (TH1*)higgs_invmass_njet_hists.At(B);
    TH1* jet_ptspectrum_njet = (TH1*)jet_ptspectrum_njet_hists.At(B);
    TH1* jet_phispectrum_njet = (TH1*)jet_phispectrum_njet_hists.At(B);
    TH1* jet_distance_njet = (TH1*)jet_distance_njet_hists.At(B);
    TH1* jet_distance_diff_njet = (TH1*)jet_distance_diff_njet_hists.At(B);

    njet_jet_perp->Write();
    njet_jet_phi->Write();
    njet_jet_eta->Write();
    jetmass_diff->Write();
    jetperp_diff->Write();
    jetphi_diff->Write();
    tau2_diff->Write();
    tau1_diff->Write();
    higgs_invmass_njet->Write();
    jet_ptspectrum_njet->Write();
    jet_phispectrum_njet->Write();
    jet_distance_njet->Write();
    jet_distance_diff_njet->Write();

    double higgs_invmass_njet_scale = 1/higgs_invmass_njet->Integral(0, 101);
    higgs_invmass_njet->Scale(higgs_invmass_akt_scale);

    if (B == 0) higgs_invmass_njet->SetLineColor(kGreen);
    if (B == 1) higgs_invmass_njet->SetLineColor(kYellow);
    if (B == 2) higgs_invmass_njet->SetLineColor(kBlue);
    if (B == 3) higgs_invmass_njet->SetLineColor(kRed);

    if (B == 2 || B == 3) {
      higgs_invmass_njet->Draw("SAMES");
      leg_mass_2jets->AddEntry(higgs_invmass_njet, "#beta = " + (TString)ss.str(), "L");
    }
    if (higgs_invmass_njet->GetMaximum() > max_val) max_val = higgs_invmass_njet->GetMaximum();
    double efficiency = higgs_invmass_njet->Integral(20, 30);
    higgs_efficiencies_2jets[B + 1][i_perp] = efficiency;

    // double njet_angulardiff = jet_phispectrum_njet->GetMean();
    // double njet_angulardiff = calcMedian1(jet_phispectrum_njet);
    double njet_angulardiff = calcMedian1(jet_distance_njet);
    double njet_angulardiff_width = jet_phispectrum_njet->GetRMS();
    angulardiff_2jets[B + 1][i_perp] = njet_angulardiff;
    angulardiff_width_2jets[B + 1][i_perp] = njet_angulardiff_width;
    angulardiff_1jet[B + 1][i_perp] = jet_distance_diff_njet->GetMean();
    // angulardiff_1jet[B + 1][i_perp] = calcMedian1(jet_distance_diff_njet);
    angulardiff_prop_1jet[B + 1][i_perp] = (double)jet_distance_diff_njet->Integral(0, (double)jet_distance_diff_njet->GetNbinsX()/10.0)/jet_distance_diff_njet->Integral(0, jet_distance_diff_njet->GetNbinsX() + 1);
  }
  higgs_invmass_akt->SetMaximum(max_val);
  leg_mass_2jets->Draw("SAMES");
  mass_compare_2jets->Write();
  mass_compare_2jets->Print("mass_compare_2jets_pt" + (TString)pythia_ptcut + ".eps", "eps");


  TCanvas *phi_compare_1jet = new TCanvas("phi_compare_1jet", "phi_compare_1jet", 600, 600);
  phi_compare_1jet->cd();
  phi_compare_1jet->SetLogy();

  double max_val_1jet_phi = 0;
  TLegend *leg_phi_1jet = new TLegend(0.58, 0.58, 0.88, 0.88);
  leg_phi_1jet->SetFillColor(kWhite);
  leg_phi_1jet->SetLineColor(kWhite);

  double firstquartiles[n_betas];
  double medians[n_betas];
  double thirdquartiles[n_betas];

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1* jetphi_diff = (TH1*)jetphi_diff_hists.At(B);

    const Int_t nq = 100;
    Double_t xq[nq];  // position where to compute the quantiles in [0,1]
    Double_t yq[nq];  // array to contain the quantiles
    for (Int_t i=0;i<nq;i++) { 
      xq[i] = Float_t(i+1)/nq;
    }

    jetphi_diff->GetQuantiles(nq,yq,xq);
    jet_distance_1jet_firstquartile[B + 1][i_perp] = yq[25];
    // jet_distance_1jet_median[B + 1][i_perp] = calcMedian1(jetphi_diff)  ;
    jet_distance_1jet_median[B + 1][i_perp] = jetphi_diff->GetMean();
    jet_distance_1jet_thirdquartile[B + 1][i_perp] = yq[75];

    double jetphi_diff_scale = 1/jetphi_diff->Integral(0, 101);
    jetphi_diff->Scale(jetphi_diff_scale);
    jetphi_diff->SetTitle("1-jet and hardest AKT Distance(pT > " + (TString)pythia_ptcut + ", R_{0} = 0.5)");
    jetphi_diff->GetXaxis()->SetTitle("#Delta R");

    if (B == 0) jetphi_diff->SetLineColor(kGreen);
    if (B == 1) jetphi_diff->SetLineColor(kYellow);
    if (B == 2) {
      jetphi_diff->SetLineColor(kBlue);
      jetphi_diff->Draw();
    }
    if (B == 3) {
      jetphi_diff->SetLineColor(kRed);
      jetphi_diff->Draw("SAMES");
    }

    if (B == 2 || B == 3) {
      leg_phi_1jet->AddEntry(jetphi_diff, "#beta = " + (TString)ss.str(), "L");
    }
    if (jetphi_diff->GetMaximum() > max_val_1jet_phi) {
      max_val_1jet_phi = jetphi_diff->GetMaximum();
      jetphi_diff->SetMaximum(1.2*max_val_1jet_phi);
    }

  }

  leg_phi_1jet->Draw("SAMES");
  phi_compare_1jet->Write();
  phi_compare_1jet->Print("phi_compare_1jet_pt" + (TString)pythia_ptcut + ".eps", "eps");


  TCanvas *perp_compare_1jet = new TCanvas("perp_compare_1jet", "perp_compare_1jet", 600, 600);
  perp_compare_1jet->cd();
  perp_compare_1jet->SetLogy();

  double akt_jet_perp_scale = 1/akt_jet_perp->Integral(0, 101);
  double max_val_1jet_perp = 0;
  TLegend *leg_perp_1jet = new TLegend(0.12, 0.48, 0.48, 0.88);
  leg_perp_1jet->SetFillColor(kWhite);
  leg_perp_1jet->SetLineColor(kWhite);

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1* jetperp_diff = (TH1*)jetperp_diff_hists.At(B);

    double jetperp_diff_scale = 1/jetperp_diff->Integral(0, 101);
    jetperp_diff->Scale(akt_jet_perp_scale);

    jetperp_diff->SetTitle("1-jet and hardest AKT pT difference (pT > " + (TString)pythia_ptcut + ", R_{0} = 0.5)");
    jetperp_diff->GetXaxis()->SetTitle("p_{T}_{akt} - p_{T}_{1jet}");

    if (B == 0) jetperp_diff->SetLineColor(kGreen);
    if (B == 1) jetperp_diff->SetLineColor(kYellow);
    if (B == 2) jetperp_diff->SetLineColor(kBlue);
    if (B == 3) jetperp_diff->SetLineColor(kRed);

    if (B == 2 || B == 3) {
      jetperp_diff->Draw("SAMES");
      leg_perp_1jet->AddEntry(jetperp_diff, "#beta = " + (TString)ss.str(), "L");
    }
    if (jetperp_diff->GetMaximum() > max_val_1jet_perp) {
      max_val_1jet_perp = jetperp_diff->GetMaximum();
      jetperp_diff->SetMaximum(1.2*max_val_1jet_perp);
    }
  }

  leg_perp_1jet->Draw("SAMES");
  perp_compare_1jet->Write();
  perp_compare_1jet->Print("perp_compare_1jet_pt" + (TString)pythia_ptcut + ".eps", "eps");

  TCanvas *mass_compare_1jet = new TCanvas("mass_compare_1jet", "mass_compare_1jet", 600, 600);
  mass_compare_1jet->cd();

  double akt_jet_mass_scale = 1/akt_jet_mass->Integral(0, 101);
  akt_jet_mass->Scale(akt_jet_mass_scale);
  akt_jet_mass->SetLineColor(kBlack);
  akt_jet_mass->SetTitle("1-jet Mass (pT > " + (TString)pythia_ptcut + ", R_{0} = 0.5)");
  akt_jet_mass->GetXaxis()->SetTitle("m_{j} [GeV]");
  akt_jet_mass->SetMinimum(0.0);
  double max_val_1jet = akt_jet_mass->GetMaximum();

  akt_jet_mass->Draw();
  TLegend *leg_mass_1jet = new TLegend(0.48, 0.48, 0.88, 0.88);
  leg_mass_1jet->SetFillColor(kWhite);
  leg_mass_1jet->SetLineColor(kWhite);
  leg_mass_1jet->AddEntry(akt_jet_mass, "Anti-KT", "L");
  double akt_efficiency_1jet = akt_jet_mass->Integral(20, 30);
  higgs_efficiencies_1jet[0][i_perp] = akt_efficiency_1jet;

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1* njet_jet_mass = (TH1*)njet_jet_mass_hists.At(B);

    double njet_jet_mass_scale = 1/njet_jet_mass->Integral(0, 101);
    njet_jet_mass->Scale(akt_jet_mass_scale);

    if (B == 0) njet_jet_mass->SetLineColor(kGreen);
    if (B == 1) njet_jet_mass->SetLineColor(kYellow);
    if (B == 2) njet_jet_mass->SetLineColor(kBlue);
    if (B == 3) njet_jet_mass->SetLineColor(kRed);

    if (B == 2 || B == 3) {
      njet_jet_mass->Draw("SAMES");
      leg_mass_1jet->AddEntry(njet_jet_mass, "#beta = " + (TString)ss.str(), "L");
    }
    if (njet_jet_mass->GetMaximum() > max_val_1jet) max_val_1jet = njet_jet_mass->GetMaximum();
    double efficiency = njet_jet_mass->Integral(20, 30);
    higgs_efficiencies_1jet[B + 1][i_perp] = efficiency;

  }

  akt_jet_mass->SetMaximum(max_val);
  leg_mass_1jet->Draw("SAMES");
  mass_compare_1jet->Write();
  mass_compare_1jet->Print("mass_compare_1jet_pt" + (TString)pythia_ptcut + ".eps", "eps");

  TCanvas *mass_compare_3jets = new TCanvas("mass_compare_3jets", "mass_compare_3jets", 600, 600);
  mass_compare_3jets->cd();

  double threejet_invmass_akt_scale = 1/threejet_invmass_akt->Integral(0, 101);
  threejet_invmass_akt->Scale(threejet_invmass_akt_scale);
  threejet_invmass_akt->SetLineColor(kBlack);
  threejet_invmass_akt->SetTitle("3-jet Mass (pT > " + (TString)pythia_ptcut + ", R_{0} = 0.5)");
  threejet_invmass_akt->GetXaxis()->SetTitle("m_{jjj} [GeV]");
  threejet_invmass_akt->SetMinimum(0.0);
  double max_val_mass_3jets = threejet_invmass_akt->GetMaximum();

  threejet_invmass_akt->Draw();
  TLegend *leg_mass_3jets = new TLegend(0.48, 0.48, 0.88, 0.88);
  leg_mass_3jets->SetFillColor(kWhite);
  leg_mass_3jets->SetLineColor(kWhite);
  leg_mass_3jets->AddEntry(threejet_invmass_akt, "Anti-KT", "L");

  double akt_efficiency_3jets = threejet_invmass_akt->Integral(20, 30);
  higgs_efficiencies_3jets[0][i_perp] = akt_efficiency_3jets;
  higgs_efficiencies_32ratio[0][i_perp] = (double)akt_efficiency_3jets/akt_efficiency;


  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1* threejet_invmass_njet = (TH1*)threejet_invmass_njet_hists.At(B);
    TH1* higgs_invmass_njet = (TH1*)higgs_invmass_njet_hists.At(B);
    TH1* threejet_thirdjet_mass_njet = (TH1*)threejet_thirdjet_mass_njet_hists.At(B);
    TH1* threejet_thirdjet_perp_njet = (TH1*)threejet_thirdjet_perp_njet_hists.At(B);

    threejet_invmass_njet->Write();
    threejet_thirdjet_mass_njet->Write();
    threejet_thirdjet_perp_njet->Write();

    double threejet_invmass_njet_scale = 1/threejet_invmass_njet->Integral(0, 101);
    threejet_invmass_njet->Scale(threejet_invmass_njet_scale);
    double higgs_invmass_njet_scale = 1/higgs_invmass_njet->Integral(0, 101);
    higgs_invmass_njet->Scale(higgs_invmass_njet_scale);


    if (B == 0) threejet_invmass_njet->SetLineColor(kGreen);
    if (B == 1) threejet_invmass_njet->SetLineColor(kYellow);
    if (B == 2) threejet_invmass_njet->SetLineColor(kBlue);
    if (B == 3) threejet_invmass_njet->SetLineColor(kRed);

    if (B == 2 || B == 3) {
      threejet_invmass_njet->Draw("SAMES");
      leg_mass_3jets->AddEntry(threejet_invmass_njet, "#beta = " + (TString)ss.str(), "L");
    }
    if (threejet_invmass_njet->GetMaximum() > max_val_mass_3jets) max_val_mass_3jets = threejet_invmass_njet->GetMaximum();

    double twojet_efficiency = higgs_invmass_njet->Integral(20, 30);
    double efficiency = threejet_invmass_njet->Integral(20, 30);
    higgs_efficiencies_3jets[B + 1][i_perp] = efficiency;
    higgs_efficiencies_32ratio[B + 1][i_perp] = (double)efficiency/twojet_efficiency;

  }

  threejet_invmass_akt->SetMaximum(max_val);
  leg_mass_3jets->Draw("SAMES");
  mass_compare_3jets->Write();
  mass_compare_3jets->Print("mass_compare_3jets_pt" + (TString)pythia_ptcut + ".eps", "eps");


  //Prevent memory leaks  
  delete jetDef;
  // delete mass_compare;
  // delete leg_mass;

}

TMultiGraph *higgs_efficiencies_1jet_biggraph = new TMultiGraph();
TMultiGraph *higgs_efficiencies_2jets_biggraph = new TMultiGraph();
TMultiGraph *higgs_efficiencies_3jets_biggraph = new TMultiGraph();
TMultiGraph *higgs_efficiencies_32ratio_biggraph = new TMultiGraph();
TMultiGraph *angulardiff_1jet_biggraph = new TMultiGraph();
TMultiGraph *angulardiff_prop_1jet_biggraph = new TMultiGraph();
TMultiGraph *angulardiff_2jets_biggraph = new TMultiGraph();
TMultiGraph *jet_distance_diff_1jet_biggraph = new TMultiGraph();

TLegend *leg_higgs_1jet = new TLegend(0.12, 0.65, 0.35, 0.88);
TLegend *leg_higgs_2jets = new TLegend(0.12, 0.65, 0.35, 0.88);
TLegend *leg_higgs_3jets = new TLegend(0.12, 0.65, 0.35, 0.88);
TLegend *leg_higgs_32ratio = new TLegend(0.12, 0.65, 0.35, 0.88);
TLegend *leg_angulardiff_1jet = new TLegend(0.6, 0.6, 0.88, 0.88);
TLegend *leg_angulardiff_prop_1jet = new TLegend(0.6, 0.6, 0.88, 0.88);
TLegend *leg_angulardiff_2jets = new TLegend(0.12, 0.6, 0.4, 0.88);
TLegend *leg_jet_distance_diff_1jet = new TLegend(0.12, 0.6, 0.4, 0.88);

leg_higgs_1jet->SetFillColor(kWhite);
leg_higgs_1jet->SetLineWidth(0);
leg_higgs_1jet->SetLineColor(kWhite);
leg_higgs_2jets->SetFillColor(kWhite);
leg_higgs_2jets->SetLineWidth(0);
leg_higgs_2jets->SetLineColor(kWhite);
leg_higgs_3jets->SetFillColor(kWhite);
leg_higgs_3jets->SetLineWidth(0);
leg_higgs_3jets->SetLineColor(kWhite);
leg_higgs_32ratio->SetFillColor(kWhite);
leg_higgs_32ratio->SetLineWidth(0);
leg_higgs_32ratio->SetLineColor(kWhite);
leg_angulardiff_1jet->SetFillColor(kWhite);
leg_angulardiff_1jet->SetLineWidth(0);
leg_angulardiff_1jet->SetLineColor(kWhite);
leg_angulardiff_prop_1jet->SetFillColor(kWhite);
leg_angulardiff_prop_1jet->SetLineWidth(0);
leg_angulardiff_prop_1jet->SetLineColor(kWhite);
leg_angulardiff_2jets->SetFillColor(kWhite);
leg_angulardiff_2jets->SetLineWidth(0);
leg_angulardiff_2jets->SetLineColor(kWhite);
leg_jet_distance_diff_1jet->SetFillColor(kWhite);
leg_jet_distance_diff_1jet->SetLineColor(kWhite);

for (int i = 0; i < n_betas + 1; i++) {

  double beta;
  if (i != 0) beta = betalist[i - 1];

  ostringstream ss;
  ss << beta;



  TGraph *higgs_efficiencies_1jet_graph = new TGraph(n_perps, perplist, higgs_efficiencies_1jet[i]);
  TGraph *higgs_efficiencies_2jets_graph = new TGraph(n_perps, perplist, higgs_efficiencies_2jets[i]);
  TGraph *higgs_efficiencies_3jets_graph = new TGraph(n_perps, perplist, higgs_efficiencies_3jets[i]);
  TGraph *higgs_efficiencies_32ratio_graph = new TGraph(n_perps, perplist, higgs_efficiencies_32ratio[i]);
  TGraph *angulardiff_1jet_graph = new TGraph(n_perps, perplist, angulardiff_1jet[i]);
  TGraph *angulardiff_prop_1jet_graph = new TGraph(n_perps, perplist, angulardiff_prop_1jet[i]);
  TGraph *angulardiff_2jets_graph = new TGraph(n_perps, perplist, angulardiff_2jets[i]);
  // TGraphAsymmErrors *jet_distance_diff_1jet_graph = new TGraphAsymmErrors(n_perps, perplist, jet_distance_1jet_median[i], 0, 0, jet_distance_1jet_firstquartile[i], jet_distance_1jet_thirdquartile[i]);
  TGraph *jet_distance_diff_1jet_graph = new TGraph(n_perps, perplist, jet_distance_1jet_median[i]);

  higgs_efficiencies_1jet_graph->SetMarkerStyle(3);
  higgs_efficiencies_1jet_graph->SetMarkerSize(2);
  higgs_efficiencies_2jets_graph->SetMarkerStyle(3);
  higgs_efficiencies_2jets_graph->SetMarkerSize(2);
  higgs_efficiencies_3jets_graph->SetMarkerStyle(3);
  higgs_efficiencies_3jets_graph->SetMarkerSize(2);
  higgs_efficiencies_32ratio_graph->SetMarkerStyle(3);
  higgs_efficiencies_32ratio_graph->SetMarkerSize(2);
  angulardiff_1jet_graph->SetMarkerStyle(3);
  angulardiff_1jet_graph->SetMarkerSize(2);
  angulardiff_prop_1jet_graph->SetMarkerStyle(3);
  angulardiff_prop_1jet_graph->SetMarkerSize(2);
  angulardiff_2jets_graph->SetMarkerStyle(3);
  angulardiff_2jets_graph->SetMarkerSize(2);
  jet_distance_diff_1jet_graph->SetMarkerStyle(3);
  jet_distance_diff_1jet_graph->SetMarkerSize(2);

  if (i == 0) {
    higgs_efficiencies_1jet_graph->SetMarkerColor(kBlack);
    higgs_efficiencies_2jets_graph->SetMarkerColor(kBlack);
    higgs_efficiencies_3jets_graph->SetMarkerColor(kBlack);
    higgs_efficiencies_32ratio_graph->SetMarkerColor(kBlack);
    angulardiff_2jets_graph->SetMarkerColor(kBlack);
    angulardiff_2jets_graph->SetLineColor(kBlack);
    leg_higgs_1jet->AddEntry(higgs_efficiencies_1jet_graph, "Anti-KT", "p");
    leg_higgs_2jets->AddEntry(higgs_efficiencies_2jets_graph, "Anti-KT", "p");
    leg_higgs_3jets->AddEntry(higgs_efficiencies_3jets_graph, "Anti-KT", "p");
    leg_higgs_32ratio->AddEntry(higgs_efficiencies_32ratio_graph, "Anti-KT", "p");
    leg_angulardiff_2jets->AddEntry(angulardiff_2jets_graph, "Anti-KT", "p");
  }
  if (i == 1) {
    higgs_efficiencies_1jet_graph->SetMarkerColor(kGreen);
    higgs_efficiencies_2jets_graph->SetMarkerColor(kGreen);
    higgs_efficiencies_3jets_graph->SetMarkerColor(kGreen);
    higgs_efficiencies_1jet_graph->SetLineColor(kGreen);
    higgs_efficiencies_2jets_graph->SetLineColor(kGreen);
    higgs_efficiencies_3jets_graph->SetLineColor(kGreen);
    higgs_efficiencies_32ratio_graph->SetMarkerColor(kGreen);
    higgs_efficiencies_32ratio_graph->SetLineColor(kGreen);
    angulardiff_1jet_graph->SetMarkerColor(kGreen);
    angulardiff_prop_1jet_graph->SetMarkerColor(kGreen);
    angulardiff_2jets_graph->SetMarkerColor(kGreen);
    angulardiff_2jets_graph->SetLineColor(kGreen);
    jet_distance_diff_1jet_graph->SetMarkerColor(kGreen);
    jet_distance_diff_1jet_graph->SetLineColor(kGreen);
  }
  if (i == 2) {
    higgs_efficiencies_1jet_graph->SetMarkerColor(kYellow);
    higgs_efficiencies_2jets_graph->SetMarkerColor(kYellow);
    higgs_efficiencies_3jets_graph->SetMarkerColor(kYellow);
    higgs_efficiencies_1jet_graph->SetLineColor(kYellow);
    higgs_efficiencies_2jets_graph->SetLineColor(kYellow);
    higgs_efficiencies_3jets_graph->SetLineColor(kYellow);
    higgs_efficiencies_32ratio_graph->SetMarkerColor(kYellow);
    higgs_efficiencies_32ratio_graph->SetLineColor(kYellow);
    angulardiff_1jet_graph->SetMarkerColor(kYellow);
    angulardiff_prop_1jet_graph->SetMarkerColor(kYellow);
    angulardiff_2jets_graph->SetMarkerColor(kYellow);
    angulardiff_2jets_graph->SetLineColor(kYellow);
    jet_distance_diff_1jet_graph->SetMarkerColor(kYellow);
    jet_distance_diff_1jet_graph->SetLineColor(kYellow);
  }
  if (i == 3) {
    higgs_efficiencies_1jet_graph->SetMarkerColor(kBlue);
    higgs_efficiencies_2jets_graph->SetMarkerColor(kBlue);
    higgs_efficiencies_3jets_graph->SetMarkerColor(kBlue);
    higgs_efficiencies_1jet_graph->SetLineColor(kBlue);
    higgs_efficiencies_2jets_graph->SetLineColor(kBlue);
    higgs_efficiencies_3jets_graph->SetLineColor(kBlue);
    higgs_efficiencies_32ratio_graph->SetMarkerColor(kBlue);
    higgs_efficiencies_32ratio_graph->SetLineColor(kBlue);
    angulardiff_1jet_graph->SetMarkerColor(kBlue);
    angulardiff_prop_1jet_graph->SetMarkerColor(kBlue);
    angulardiff_2jets_graph->SetMarkerColor(kBlue);
    angulardiff_2jets_graph->SetLineColor(kBlue);
    jet_distance_diff_1jet_graph->SetMarkerColor(kBlue);
    jet_distance_diff_1jet_graph->SetLineColor(kBlue);

  }
  if (i == 4) {
    higgs_efficiencies_1jet_graph->SetMarkerColor(kRed);
    higgs_efficiencies_2jets_graph->SetMarkerColor(kRed);
    higgs_efficiencies_3jets_graph->SetMarkerColor(kRed);
    higgs_efficiencies_1jet_graph->SetLineColor(kRed);
    higgs_efficiencies_2jets_graph->SetLineColor(kRed);
    higgs_efficiencies_3jets_graph->SetLineColor(kRed);
    higgs_efficiencies_32ratio_graph->SetMarkerColor(kRed);
    higgs_efficiencies_32ratio_graph->SetLineColor(kRed);
    angulardiff_1jet_graph->SetMarkerColor(kRed);
    angulardiff_prop_1jet_graph->SetMarkerColor(kRed);
    angulardiff_2jets_graph->SetMarkerColor(kRed);
    angulardiff_2jets_graph->SetLineColor(kRed);
    jet_distance_diff_1jet_graph->SetMarkerColor(kRed);
    jet_distance_diff_1jet_graph->SetLineColor(kRed);

  }

  if (i == 0 || i == 3 || i == 4) {
    higgs_efficiencies_1jet_biggraph->Add(higgs_efficiencies_1jet_graph);
    higgs_efficiencies_2jets_biggraph->Add(higgs_efficiencies_2jets_graph);
    higgs_efficiencies_3jets_biggraph->Add(higgs_efficiencies_3jets_graph);
    higgs_efficiencies_32ratio_biggraph->Add(higgs_efficiencies_32ratio_graph);
    angulardiff_2jets_biggraph->Add(angulardiff_2jets_graph);
    angulardiff_1jet_biggraph->Add(angulardiff_1jet_graph);
    angulardiff_prop_1jet_biggraph->Add(angulardiff_prop_1jet_graph);
    if (i != 0) {
      jet_distance_diff_1jet_biggraph->Add(jet_distance_diff_1jet_graph);
      leg_higgs_1jet->AddEntry(higgs_efficiencies_1jet_graph, "#beta = " + (TString)ss.str(), "p");
      leg_higgs_2jets->AddEntry(higgs_efficiencies_2jets_graph, "#beta = " + (TString)ss.str(), "p");
      leg_higgs_3jets->AddEntry(higgs_efficiencies_3jets_graph, "#beta = " + (TString)ss.str(), "p");
      leg_higgs_32ratio->AddEntry(higgs_efficiencies_32ratio_graph, "#beta = " + (TString)ss.str(), "p");
      leg_angulardiff_1jet->AddEntry(angulardiff_1jet_graph, "#beta = " + (TString)ss.str(), "p");
      leg_angulardiff_prop_1jet->AddEntry(angulardiff_prop_1jet_graph, "#beta = " + (TString)ss.str(), "p");
      leg_angulardiff_2jets->AddEntry(angulardiff_2jets_graph, "#beta = " + (TString)ss.str(), "p");
      leg_jet_distance_diff_1jet->AddEntry(jet_distance_diff_1jet_graph, "#beta = " + (TString)ss.str(), "p");
    }
  }
}

TCanvas *higgs_1jet_can = new TCanvas("higgs_1jet_can", "higgs_1jet_can", 800, 600);
higgs_1jet_can->cd();
higgs_efficiencies_1jet_biggraph->Draw("ACP");
leg_higgs_1jet->Draw();
higgs_efficiencies_1jet_biggraph->SetMinimum(0.);     
higgs_efficiencies_1jet_biggraph->SetMaximum(1.0);
higgs_efficiencies_1jet_biggraph->SetTitle("Higgs Efficiencies for 1-jettiness with 100 < m_{j} < 150 (R_{0} = 0.5)");
higgs_efficiencies_1jet_biggraph->GetXaxis()->SetTitle("pT_{min} (GeV)");
higgs_efficiencies_1jet_biggraph->GetYaxis()->SetTitle("% mass 100-150");
higgs_efficiencies_1jet_biggraph->GetYaxis()->SetTitleOffset(1.3);
higgs_1jet_can->Write();
higgs_1jet_can->Print("higgs_efficiencies_1jet.eps", "eps");

TCanvas *higgs_2jets_can = new TCanvas("higgs_2jets_can", "higgs_2jets_can", 800, 600);
higgs_2jets_can->cd();
higgs_efficiencies_2jets_biggraph->Draw("ACP");
leg_higgs_2jets->Draw();
higgs_efficiencies_2jets_biggraph->SetMinimum(0.);     
higgs_efficiencies_2jets_biggraph->SetMaximum(1.0);
higgs_efficiencies_2jets_biggraph->SetTitle("Higgs Efficiencies for 2-jettiness with 100 < m_{jj} < 150 (R_{0} = 0.5)");
higgs_efficiencies_2jets_biggraph->GetXaxis()->SetTitle("pT_{min} (GeV)");
higgs_efficiencies_2jets_biggraph->GetYaxis()->SetTitle("% mass 100-150");
higgs_efficiencies_2jets_biggraph->GetYaxis()->SetTitleOffset(1.3);
higgs_2jets_can->Write();
higgs_2jets_can->Print("higgs_efficiencies_2jets.eps", "eps");

TCanvas *higgs_3jets_can = new TCanvas("higgs_3jets_can", "higgs_3jets_can", 800, 600);
higgs_3jets_can->cd();
higgs_efficiencies_3jets_biggraph->Draw("ACP");
leg_higgs_3jets->Draw();
higgs_efficiencies_3jets_biggraph->SetMinimum(0.);     
higgs_efficiencies_3jets_biggraph->SetMaximum(1.0);
higgs_efficiencies_3jets_biggraph->SetTitle("Higgs Efficiencies for 3-jettiness with 100 < m_{jjj} < 150 (R_{0} = 0.5, R_{cut} = 1.0)");
higgs_efficiencies_3jets_biggraph->GetXaxis()->SetTitle("pT_{min} (GeV)");
higgs_efficiencies_3jets_biggraph->GetYaxis()->SetTitle("% mass 100-150");
higgs_efficiencies_3jets_biggraph->GetYaxis()->SetTitleOffset(1.3);
higgs_3jets_can->Write();
higgs_3jets_can->Print("higgs_efficiencies_3jets.eps", "eps");

TCanvas *higgs_32ratio_can = new TCanvas("higgs_32ratio_can", "higgs_32ratio_can", 800, 600);
higgs_32ratio_can->cd();
higgs_efficiencies_32ratio_biggraph->Draw("ACP");
leg_higgs_32ratio->Draw();
higgs_efficiencies_32ratio_biggraph->SetMinimum(0.5);     
higgs_efficiencies_32ratio_biggraph->SetMaximum(1.5);
higgs_efficiencies_32ratio_biggraph->SetTitle("Ratio of 3-jettiness/2-jettiness Higgs Efficiencies (R_{0} = 0.5, R_{cut} = 1.0)");
higgs_efficiencies_32ratio_biggraph->GetXaxis()->SetTitle("pT_{min} (GeV)");
higgs_efficiencies_32ratio_biggraph->GetYaxis()->SetTitle("3-jet eff/2-jet eff");
higgs_efficiencies_32ratio_biggraph->GetYaxis()->SetTitleOffset(1.3);
higgs_32ratio_can->Write();
higgs_32ratio_can->Print("higgs_efficiencies_32ratio.eps", "eps");

TCanvas *angulardiff_1jet_can = new TCanvas("angulardiff_1jet_can", "angulardiff_1jet_can", 800, 600);
angulardiff_1jet_can->cd();
angulardiff_1jet_biggraph->Draw("ACP");
leg_angulardiff_1jet->Draw();
angulardiff_1jet_biggraph->SetMinimum(0.);
angulardiff_1jet_biggraph->SetMaximum(0.25);
angulardiff_1jet_biggraph->SetTitle("Distance between 1-jettiness jet and closest 2-jettiness jet");
angulardiff_1jet_biggraph->GetXaxis()->SetTitle("pT_{min} (GeV)");
angulardiff_1jet_biggraph->GetYaxis()->SetTitle("Angular Distance");
angulardiff_1jet_biggraph->GetYaxis()->SetTitleOffset(1.3);
angulardiff_1jet_can->Write();
angulardiff_1jet_can->Print("angulardiff_1jet.eps", "eps");


TCanvas *angulardiff_prop_1jet_can = new TCanvas("angulardiff_prop_1jet_can", "angulardiff_prop_1jet_can", 800, 600);
angulardiff_prop_1jet_can->cd();
angulardiff_prop_1jet_biggraph->Draw("ACP");
leg_angulardiff_prop_1jet->Draw();
angulardiff_prop_1jet_biggraph->SetMinimum(0.);
angulardiff_prop_1jet_biggraph->SetMaximum(1.0);
angulardiff_prop_1jet_biggraph->SetTitle("Proportion of 1-jettiness jets aligned with 2-jettiness jets");
angulardiff_prop_1jet_biggraph->GetXaxis()->SetTitle("pT_{min} (GeV)");
angulardiff_prop_1jet_biggraph->GetYaxis()->SetTitle("Angular Distance");
angulardiff_prop_1jet_biggraph->GetYaxis()->SetTitleOffset(1.3);
angulardiff_prop_1jet_can->Write();
angulardiff_prop_1jet_can->Print("angulardiff_prop_1jet.eps", "eps");

TCanvas *angulardiff_2jets_can = new TCanvas("angulardiff_2jets_can", "angulardiff_2jets_can", 800, 600);
angulardiff_2jets_can->cd();
angulardiff_2jets_biggraph->Draw("ACP");
leg_angulardiff_2jets->Draw();
angulardiff_2jets_biggraph->SetMinimum(0.);
angulardiff_2jets_biggraph->SetMaximum(4.0);
angulardiff_2jets_biggraph->SetTitle("Angular Distance between two jets");
angulardiff_2jets_biggraph->GetXaxis()->SetTitle("pT_{min} (GeV)");
angulardiff_2jets_biggraph->GetYaxis()->SetTitle("Angular Distance");
angulardiff_2jets_biggraph->GetYaxis()->SetTitleOffset(1.3);
angulardiff_2jets_can->Write();
angulardiff_2jets_can->Print("angulardiff_2jets.eps", "eps");

TCanvas* jet_distance_diff_1jet_can = new TCanvas("jet_distance_diff_1jet_can", "jet_distance_diff_1jet_can", 800, 600);
jet_distance_diff_1jet_can->cd();
jet_distance_diff_1jet_biggraph->Draw("ACP");
leg_jet_distance_diff_1jet->Draw();
jet_distance_diff_1jet_biggraph->SetMinimum(0.);
jet_distance_diff_1jet_biggraph->SetMaximum(0.1);
jet_distance_diff_1jet_biggraph->SetTitle("Avg Angular Distance between AKT and 1-jet");
jet_distance_diff_1jet_biggraph->GetXaxis()->SetTitle("pT_{min} (GeV)");
jet_distance_diff_1jet_biggraph->GetYaxis()->SetTitle("Avg Angular Distance");
jet_distance_diff_1jet_biggraph->GetYaxis()->SetTitleOffset(1.3);
jet_distance_diff_1jet_can->Write();
jet_distance_diff_1jet_can->Print("jet_distance_diff_1jet.eps", "eps");

delete outFile;

return 0;
}