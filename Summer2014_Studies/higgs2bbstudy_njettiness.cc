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
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"

//Root headers
#include "TROOT.h"
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
#include "TPaveText.h"

using namespace std;
using namespace Pythia8;
using namespace fastjet;
using namespace fastjet::contrib;

// class GeneralERecombiner : public fastjet::JetDefinition::Recombiner {
// public:
// // Constructor to choose value of alpha (defaulted to 1 for normal p_{T} sum)
//   GeneralERecombiner(double delta) : _delta(delta) {}
  
//   std::string description() const {
//    return "General E-scheme recombination";
//  }

//  void recombine(const fastjet::PseudoJet & pa, const fastjet::PseudoJet & pb, fastjet::PseudoJet & pab) const {

//   double weighta = pow(pa.perp(), _delta);
//   double weightb = pow(pb.perp(), _delta);

//   double perp_ab = pa.perp() + pb.perp();
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


PseudoJet findMinMass(vector<PseudoJet> initial_jets, double zcut) {

  PseudoJet minmass_jet;

  std::string bitmask(2, 1); // K leading 1's
  bitmask.resize(initial_jets.size(), 0); // N-K trailing 0's

  double minmass = std::numeric_limits<int>::max();
  // double maxperp = 0;
  do {
    vector<int> axis_indices;
    for (int i = 0; i < initial_jets.size(); ++i) {
      if (bitmask[i]) axis_indices.push_back(i);
    }
    PseudoJet temp_jet(0,0,0,0);
    double min_perp = std::numeric_limits<int>::max();
    double sum_perp = 0;
    for (int j = 0; j < axis_indices.size(); j++) {
      temp_jet = join(temp_jet,initial_jets[axis_indices[j]]);
      // if (initial_jets[axis_indices[j]].perp() < min_perp) min_perp = initial_jets[axis_indices[j]].perp();
      // sum_perp += initial_jets[axis_indices[j]].perp();
    }
    // double temp_zfrac = min_perp/sum_perp;

    if (temp_jet.m() < minmass) {
      minmass = temp_jet.m();
      minmass_jet = temp_jet;
    }
    // if (temp_jet.perp() > maxperp && temp_zfrac > zcut) {
    //   maxperp = temp_jet.perp();
    //   minmass_jet = temp_jet;
    // }
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  return minmass_jet;
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
  int nEvents = 1000;
  gStyle->SetOptStat(0);
  // Create the ROOT application environment.
  TApplication theApp("hist", &argc, argv);
  
  TStyle *plain  = new TStyle("Plain","plain");
  plain->SetLegendBorderSize(0);
  plain->SetLegendFillColor(0);
  plain->SetTitleFont(132, "a");
  plain->SetTitleFont(132, "xy");
  plain->SetLegendFont(132);
  // plain->SetTextSize(0.05);
  plain->SetLabelSize(0.05, "xy");
  // plain->SetLabelOffset(0.003, "xy");
  plain->SetTitleSize(0.09, "a");
  plain->SetTitleSize(0.07, "y");
  plain->SetTitleSize(0.08, "x");
  plain->SetPadTopMargin(0.12);
  plain->SetPadLeftMargin(0.15);
  plain->SetPadBottomMargin(0.15);
  // plain->SetTitleOffset(1.4, "a");
  plain->SetTitleOffset(0.8, "x");
  plain->SetTitleOffset(1.1, "y");
  plain->SetHistLineWidth(4);

  // plain->SetLegendSize(12);
  plain->SetTitleBorderSize(0);
  plain->SetTitleX(0.1f);
  plain->SetTitleW(0.8f);
  plain->SetOptStat(0);
  plain->cd();

  gROOT->SetStyle("plain");
  // Create file on which histograms can be saved.  
  TFile* outFile = new TFile("higgs2bbstudy_onepassXcone_test_allpt.root", "RECREATE");

  int n_perps = 9;
  // int n_perps = 1;
  bool willDisplay = false;

  //create list of various values of beta
  vector<double> betalist;
  betalist.push_back(0.25);
  betalist.push_back(0.5);
  betalist.push_back(2.0);
  betalist.push_back(1.0);
  // betalist.push_back(3.0);
  int n_betas = betalist.size();


  vector<Color_t> colorlist;
  colorlist.push_back(kYellow);
  colorlist.push_back(8);
  colorlist.push_back(kBlue);
  colorlist.push_back(kRed);


  double perplist[n_perps];
  double higgs_efficiencies_1jet[n_betas + 1][n_perps];
  double higgs_efficiencies_2jets[n_betas + 1][n_perps];
  double higgs_efficiencies_3jets[n_betas + 1][n_perps];
  double higgs_efficiencies_32ratio[n_betas + 1][n_perps];
  double higgs_efficiencies_2jets_combined[n_betas + 1][n_perps];
  double angulardiff_2jets[n_betas + 1][n_perps];
  double angulardiff_width_2jets[n_betas + 1][n_perps];
  double angulardiff_1jet[n_betas + 1][n_perps];
  double angulardiff_prop_1jet[n_betas + 1][n_perps];
  double jet_distance_1jet_median[n_betas + 1][n_perps];
  double jet_distance_1jet_firstquartile[n_betas + 1][n_perps];
  double jet_distance_1jet_thirdquartile[n_betas + 1][n_perps];

  TObjArray jetphi_diff_allpt_allbeta(n_perps*n_betas);

  for (int i_perp = 0; i_perp < n_perps; i_perp++) {

    // double ptcut = atoi(argv[1]);
    double ptcut = (i_perp + 2)*100.0;
    // double ptcut = 800;
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

  TH1* particle_mass = new TH1F("particle_mass", "Mass of Higgs particle", 50, 0, 500);

  TH1* akt_jet_mass = new TH1F("akt_jet_mass", "Mass of individual akt jet", 50, 0, 250);
  TH1* akt_jet_perp = new TH1F("akt_jet_perp", "Perp of individual akt jet", 50, 0, 250);
  TH1* akt_jet_phi = new TH1F("akt_jet_phi", "Phi of individual akt jet", 32, 0, 6.4);
  TH1* akt_jet_eta = new TH1F("akt_jet_eta", "Eta of individual akt jet", 50, -5, 5);
  TH1* akt_jet_area = new TH1F("akt_jet_area", "Area of individual akt jet", 40, 0.4, 1.2);

  TH1* higgs_invmass_akt = new TH1F("higgs_invmass_akt", "Invariant Mass of Higgs (antikt)", 50, 0, 250);
  TH1* higgs_invmass_akt_combined = new TH1F("higgs_invmass_akt_combined", "", 50, 0, 500);
  TH1* jet_ptspectrum_akt = new TH1F("jet_ptspectrum_akt", "p_{T} spectrum of jets (antikt)", 50, 0, 250);
  TH1* jet_phispectrum_akt = new TH1F("jet_phispectrum_akt", "#phi spectrum of jets (antikt)", 32, 0, 3.2);
  TH1* jet_distance_akt = new TH1F("jet_distance_akt", "Distance between two jets (antikt)", 50, 0, 10.0);
  TH1* jet_distance_diff_akt = new TH1F("jet_distance_diff_akt", "Distance between individual jet and closest of two jets", 60, 0, 3.0);
  TH1* twojet_akt_jet_area = new TH1F("twojet_akt_jet_area", "Area of individual akt jet", 40, 0.4, 1.2);

  TH1* threejet_invmass_akt = new TH1F("threejet_invmass_akt", "3-jet invariant mass (antikt)", 50, 0, 500);
  TH1* threejet_thirdjet_mass_akt = new TH1F("threejet_thirdjet_mass_akt", "Third jet mass (antikt)", 50, 0, 250);
  TH1* threejet_thirdjet_perp_akt = new TH1F("threejet_thirdjet_perp_akt", "Third jet perp (antikt)", 50, 0, 250);
  TH1* threejet_thirdjet_phidiff_akt = new TH1F("threejet_thirdjet_phidiff_akt", "", 16, 0, 3.2);
  TH1* threejet_thirdjet_area_akt = new TH1F("threejet_thirdjet_area_akt", "", 50, 0, 1.5);
  TH2* threejet_mass_2jet3jet_compare_akt = new TH2F("threejet_mass_2jet3jet_compare_akt", "", 50, 0, 250, 50, 0, 250);

  TObjArray njet_jet_mass_hists(n_betas);  
  TObjArray njet_jet_perp_hists(n_betas);
  TObjArray njet_jet_phi_hists(n_betas);
  TObjArray njet_jet_eta_hists(n_betas);
  TObjArray jetmass_diff_hists(n_betas);
  TObjArray jetperp_diff_hists(n_betas);
  TObjArray jetphi_diff_hists(n_betas);
  TObjArray tau2_diff_hists(n_betas);
  TObjArray tau1_diff_hists(n_betas);
  TObjArray njet_jet_area_hists(n_betas);
  TObjArray njet_jet_areadiff_hists(n_betas);
  TObjArray twojet_njet_jet_area_hists(n_betas);

  TObjArray higgs_invmass_njet_hists(n_betas);
  TObjArray higgs_invmass_njet_combined_hists(n_betas);
  TObjArray jet_ptspectrum_njet_hists(n_betas);
  TObjArray jet_phispectrum_njet_hists(n_betas);
  TObjArray jet_distance_njet_hists(n_betas);
  TObjArray jet_distance_diff_njet_hists(n_betas);

  TObjArray threejet_invmass_njet_hists(n_betas);
  TObjArray threejet_thirdjet_mass_njet_hists(n_betas);
  TObjArray threejet_thirdjet_perp_njet_hists(n_betas);
  TObjArray threejet_thirdjet_phidiff_njet_hists(n_betas);
  TObjArray threejet_thirdjet_area_njet_hists(n_betas);
  TObjArray threejet_mass_2jet3jet_compare_njet_hists(n_betas);

  vector<int> higgs_events_2jets_combined;
  higgs_events_2jets_combined.push_back(0);

  for (int Bs = 0; Bs < n_betas; Bs++) {

    double beta = betalist[Bs];

    ostringstream ss;
    ss << beta;
    TString title;

    TH1* njet_jet_mass = new TH1F("njet_jet_mass", "Mass of individual njet jet (#beta = " + (TString)ss.str() + ")", 50, 0, 250);
    TH1* njet_jet_perp = new TH1F("njet_jet_perp", "Perp of individual njet jet (#beta = " + (TString)ss.str() + ")", 50, 0, 500);
    TH1* njet_jet_phi = new TH1F("njet_jet_phi", "Phi of individual njet jet (#beta = " + (TString)ss.str() + ")", 32, 0, 6.4);
    TH1* njet_jet_eta = new TH1F("njet_jet_eta", "Eta of individual njet jet (#beta = " + (TString)ss.str() + ")", 50, -5, 5);
    TH1* jetmass_diff = new TH1F("jetmass_diff", "difference in mass of akt and njet jets (#beta = " + (TString)ss.str() + ")", 50, 0, 250);
    TH1* jetperp_diff = new TH1F("jetperp_diff", "difference in perp of akt and njet jets (#beta = " + (TString)ss.str() + ")", 50, -200, 200);
    TH1* jetphi_diff = new TH1F("jetphi_diff", "difference in phi of akt and njet jets (#beta = " + (TString)ss.str() + ")", 50, 0, 1.0);
    TH1* tau1_diff = new TH1F("tau1_diff", "difference in #tau_{1} jets (akt - njet) (#beta = " + (TString)ss.str() + ")", 50, -200, 200);
    TH1* tau2_diff = new TH1F("tau2_diff", "difference in #tau_{2} jets (akt - njet) (#beta = " + (TString)ss.str() + ")", 50, -200, 200);
    TH1* njet_jet_area = new TH1F("njet_jet_area", "", 40, 0.4, 1.2);
    TH1* njet_jet_areadiff = new TH1F("njet_jet_area", "", 50, -0.5, 0.5);

    TH1* higgs_invmass_njet = new TH1F("higgs_invmass_njet", "Invariant Mass of Higgs (njet, beta = " + (TString)ss.str() + ")", 50, 0, 250);
    TH1* higgs_invmass_njet_combined = new TH1F("higgs_invmass_njet_combined", "Invariant Mass of Higgs (njet, beta = " + (TString)ss.str() + ")", 50, 0, 250);
    TH1* jet_ptspectrum_njet = new TH1F("jet_ptspectrum_njet", "p_{T} spectrum of jets (njet, beta = " + (TString)ss.str() + ")", 50, 0, 500);
    TH1* jet_phispectrum_njet = new TH1F("jet_phispectrum_njet", "#phi spectrum of jets (njet, beta = " + (TString)ss.str() + ")", 32, 0, 3.2);
    TH1* jet_distance_njet = new TH1F("jet_distance_njet", "Distance between two jets (njet, beta = " + (TString)ss.str() + ")", 50, 0, 10.0);
    TH1* jet_distance_diff_njet = new TH1F("jet_distance_diff_njet", "Distance between individual jet and closest of two jets (njet, beta = " + (TString)ss.str() + ")", 50, 0, 0.5);
    TH1* twojet_njet_jet_area = new TH1F("twojet_njet_jet_area", "", 40, 0.4, 1.2);

    TH1* threejet_invmass_njet = new TH1F("threejet_invmass_njet", "3-jet invariant mass (njet, beta = " + (TString)ss.str() + ")", 50, 0, 250);
    TH1* threejet_thirdjet_mass_njet = new TH1F("threejet_thirdjet_mass_njet", "Third jet mass (njet, beta = " + (TString)ss.str() + ")", 50, 0, 250);
    TH1* threejet_thirdjet_perp_njet = new TH1F("threejet_thirdjet_perp_njet", "Third jet perp (njet, beta = " + (TString)ss.str() + ")", 50, 0, 500);
    TH1* threejet_thirdjet_phidiff_njet = new TH1F("threejet_thirdjet_phidiff_njet", "", 16, 0, 3.2);
    TH1* threejet_thirdjet_area_njet = new TH1F("threejet_thirdjet_area_njet", "", 50, 0, 1.5);
    TH2* threejet_mass_2jet3jet_compare_njet = new TH2F("threejet_mass_2jet3jet_compare_njet", "", 50, 0, 500, 100, 0, 500);

    njet_jet_mass_hists.Add(njet_jet_mass);
    njet_jet_perp_hists.Add(njet_jet_perp);
    njet_jet_phi_hists.Add(njet_jet_phi);
    njet_jet_eta_hists.Add(njet_jet_eta);
    jetmass_diff_hists.Add(jetmass_diff);
    jetperp_diff_hists.Add(jetperp_diff);
    jetphi_diff_hists.Add(jetphi_diff);
    tau2_diff_hists.Add(tau2_diff);
    tau1_diff_hists.Add(tau1_diff);
    njet_jet_area_hists.Add(njet_jet_area);
    njet_jet_areadiff_hists.Add(njet_jet_areadiff);
    twojet_njet_jet_area_hists.Add(twojet_njet_jet_area);

    higgs_invmass_njet_hists.Add(higgs_invmass_njet);
    higgs_invmass_njet_combined_hists.Add(higgs_invmass_njet_combined);
    jet_ptspectrum_njet_hists.Add(jet_ptspectrum_njet);
    jet_phispectrum_njet_hists.Add(jet_phispectrum_njet);
    jet_distance_njet_hists.Add(jet_distance_njet);
    jet_distance_diff_njet_hists.Add(jet_distance_diff_njet);

    threejet_invmass_njet_hists.Add(threejet_invmass_njet);
    threejet_thirdjet_mass_njet_hists.Add(threejet_thirdjet_mass_njet);
    threejet_thirdjet_perp_njet_hists.Add(threejet_thirdjet_perp_njet);
    threejet_thirdjet_phidiff_njet_hists.Add(threejet_thirdjet_phidiff_njet);
    threejet_thirdjet_area_njet_hists.Add(threejet_thirdjet_area_njet);
    threejet_mass_2jet3jet_compare_njet_hists.Add(threejet_mass_2jet3jet_compare_njet);

    higgs_events_2jets_combined.push_back(0);
  }

  // Fastjet input (different inputs for quark vs gluon)
  vector <PseudoJet> fjInputs;

  double Rparam = 0.5;
  // double Rparam = min(0.3, (double)125/ptcut);
  double zfrac = 0.0;
  double merging_point = 2*125/Rparam;
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
      // if (abs(fj_particle.eta()) > 3.0) continue;

      fjInputs.push_back(fj_particle);
    }

    TH2F *event_display = new TH2F("event_display", "Event Display", 60, -5, 5, 60, 0, 6.4);

    particle_mass->Fill(pythia.event[iH].m());

    Selector eta_selector = SelectorAbsEtaMax(3.0);
    fjInputs = eta_selector(fjInputs);

    for (int i_part = 0; i_part < fjInputs.size(); i_part++) {
      event_display->Fill(fjInputs[i_part].eta(), fjInputs[i_part].phi(), fjInputs[i_part].perp());
    }
    
    if (iEvent % 10 == 0) cout << iEvent << endl;
    // cout << iEvent << endl << endl;

    //Run Fastjet algorithm
    vector <PseudoJet> inclusiveJets, sortedJets, centralJets;
    // ClusterSequence clustSeq(fjInputs, *jetDef);
    AreaDefinition area_def(active_area, GhostedAreaSpec(3.0));
    ClusterSequenceArea clustSeq(fjInputs, *jetDef, area_def);

    //Extract inclusive jets sorted by p_{T}
    inclusiveJets = clustSeq.inclusive_jets(0.0);

    // Sort central jets by p_{T}
    sortedJets = sorted_by_pt(inclusiveJets);

    //Cut on jets within eta range
    // Selector eta_selector = SelectorAbsEtaMax(3.0);
    centralJets = eta_selector(sortedJets);

    Selector jet_selector = SelectorNHardest(1);
    vector<PseudoJet> hardest_jet = jet_selector(centralJets);

    if (hardest_jet.size() == 1) {
      akt_jet_mass->Fill(hardest_jet[0].m());
      akt_jet_perp->Fill(hardest_jet[0].perp());
      akt_jet_phi->Fill(hardest_jet[0].phi());
      akt_jet_eta->Fill(hardest_jet[0].eta());
      akt_jet_area->Fill(hardest_jet[0].area());
    }

    //Create selectors for hardest jet in the event
    Selector twojet_selector = SelectorNHardest(2);
    vector<PseudoJet> hardest_twojets = twojet_selector(centralJets);

    PseudoJet big_jet(0,0,0,0);
    for (int i_jet = 0; i_jet < hardest_twojets.size(); i_jet++) {
      big_jet = join(big_jet, hardest_twojets[i_jet]);
      jet_ptspectrum_akt->Fill(hardest_twojets[i_jet].perp());
      twojet_akt_jet_area->Fill(hardest_twojets[i_jet].area());
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

    // if (hardest_jet[0].perp() > merging_point && (hardest_twojets[0].delta_R(hardest_twojets[1]) > 2*Rparam)) {
      // higgs_invmass_akt_combined->Fill(hardest_jet[0].m());
    // }
    // else higgs_invmass_akt_combined->Fill(big_jet.m());


    Selector threejet_selector = SelectorNHardest(3);
    vector<PseudoJet> hardest_threejets = threejet_selector(centralJets);

    PseudoJet threejet_minmassjet_akt = findMinMass(hardest_threejets, zfrac);
    PseudoJet big_jet_3jets(0,0,0,0);
    // for (int i_jet = 0; i_jet < hardest_threejets.size(); i_jet++) {
    //   big_jet_3jets = join(big_jet_3jets, hardest_threejets[i_jet]);
    // }

    // double mass_ratio_akt = (double)big_jet_3jets.m()/big_jet.m();

    // if (mass_ratio_akt > 0.9 && mass_ratio_akt < 1.1) threejet_invmass_akt->Fill(big_jet_3jets.m());
    // if (mass_ratio_akt > 0.9 && mass_ratio_akt < 1.1) threejet_invmass_akt->Fill(big_jet.m());
    // if (big_jet.m() > 100 && big_jet.m() < 150) threejet_invmass_akt->Fill(big_jet.m());
    // else threejet_invmass_akt->Fill(threejet_minmassjet_akt.m());
    // else threejet_invmass_akt->Fill(threejet_minmassjet_akt.m());

    PseudoJet third_jet = hardest_threejets[hardest_threejets.size() - 1];
    PseudoJet closer_jet = (third_jet.delta_R(hardest_twojets[0]) < third_jet.delta_R(hardest_twojets[1])) ? hardest_twojets[0] : hardest_twojets[1];

    // if (threejet_minmassjet_akt.m() < big_jet.m()) threejet_invmass_akt->Fill(big_jet_3jets.m());
    // else threejet_invmass_akt->Fill(threejet_minmassjet_akt.m());

    threejet_invmass_akt->Fill(threejet_minmassjet_akt.m());

    if ((threejet_minmassjet_akt.m() > 100 && threejet_minmassjet_akt.m() < 150) || (big_jet.m() > 100 && big_jet.m() < 150)) {
      higgs_events_2jets_combined[4]++;
    }
    // if ((third_jet.delta_R(closer_jet) < 2*Rparam) && third_jet.perp() > 50) threejet_invmass_akt->Fill(big_jet_3jets.m());
    // else threejet_invmass_akt->Fill(threejet_minmassjet_akt.m());

    threejet_thirdjet_phidiff_akt->Fill(third_jet.delta_R(closer_jet));
    threejet_thirdjet_area_akt->Fill(third_jet.area());
    // threejet_mass_2jet3jet_compare_akt->Fill(threejet_minmassjet_akt.m(), third_jet.delta_R(closer_jet));
    threejet_mass_2jet3jet_compare_akt->Fill(big_jet.m(), threejet_minmassjet_akt.m());

    // threejet_thirdjet_mass_akt->Fill(third_jet.m());
    // threejet_thirdjet_perp_akt->Fill(third_jet.perp());

    // else big_jet_3jets = big_jet;
    // // }
    // threejet_invmass_akt->Fill(big_jet_3jets.m());


    for (int B = 0; B < n_betas; B++) {

      double beta = betalist[B];
      double power = (double)1/beta;
      double delta;

      // const JetDefinition::Recombiner *recombScheme;
      // if (beta > 1) recombScheme = new GeneralERecombiner((double)1/(beta - 1));
      // else recombScheme = new WinnerTakeAllRecombiner();

      if (beta > 1) delta = (double)1/(beta - 1);
      else delta = std::numeric_limits<int>::max();

      UnnormalizedCutoffMeasure measure_function = UnnormalizedCutoffMeasure(beta, Rparam);
      // XConeCutoffMeasure measure_function = XConeCutoffMeasure(beta, Rparam);
      XConeCutoffMeasure measure_function_xcone = XConeCutoffMeasure(beta, Rparam);

      AxesStruct *axes_finder = new AxesStruct(GenRecomb_GenKT_Axes(delta, power, Rparam));

      AxesStruct *axes_finder_onepass;
      if (beta >= 1 && beta <= 3) axes_finder_onepass = new AxesStruct(OnePass_GenRecomb_GenKT_Axes(delta, power, Rparam));
      else axes_finder_onepass = axes_finder;

      // JetDefinition *jetDef = new JetDefinition(genkt_algorithm, Rparam, power, recombScheme, strategy);
      // ClusterSequence clustSeq(fjInputs, *jetDef);
      // vector<PseudoJet> exclusive_1jet_start = clustSeq.exclusive_jets(1);
      // vector<PseudoJet> exclusive_1jet = findMinAxes(fjInputs, exclusive_1jet_start, 1, beta, Rparam);

      //Use Njettiness algorithm for comparison
      // NjettinessPlugin njet_plugin_conical_1jet(1, axes_finder->def(), measure_function);
      // JetDefinition njet_def_conical_1jet(&njet_plugin_conical_1jet);
      // ClusterSequence njet_cluster_conical_1jet(fjInputs, njet_def_conical_1jet);
      // const NjettinessExtras *extras_1jet = njettiness_extras(njet_cluster_conical_1jet);
      // vector<PseudoJet> njet_conical_1jet = extras_1jet->jets();
      // vector<PseudoJet> njet_conical_1jet_axes = extras_1jet->axes();

      NjettinessPlugin njet_plugin_xcone_1jet(1, axes_finder_onepass->def(), measure_function_xcone);
      JetDefinition njet_def_xcone_1jet(&njet_plugin_xcone_1jet);
      ClusterSequence njet_cluster_xcone_1jet(fjInputs, njet_def_xcone_1jet);
      const NjettinessExtras *njet_cluster_xcone_1jet_extras = njettiness_extras(njet_cluster_xcone_1jet);
      vector<PseudoJet> njet_1jet = njet_cluster_xcone_1jet.inclusive_jets();
      vector<PseudoJet> njet_1jet_axes = njet_cluster_xcone_1jet_extras->axes();

      NjettinessPlugin njet_plugin_xcone_1jet_area(1, Manual_Axes(), measure_function_xcone);
      JetDefinition njet_def_xcone_1jet_area(&njet_plugin_xcone_1jet_area);
      njet_plugin_xcone_1jet_area.setAxes(njet_1jet_axes);
      ClusterSequenceArea njet_cluster_xcone_1jet_area(fjInputs, njet_def_xcone_1jet_area, area_def);
      vector<PseudoJet> njet_1jet_area = njet_cluster_xcone_1jet_area.inclusive_jets();

      // cout << "tau1" << endl;

      if (njet_1jet.size() == 1) {

        TH1* njet_jet_mass = (TH1*)njet_jet_mass_hists.At(B);
        TH1* njet_jet_perp = (TH1*)njet_jet_perp_hists.At(B);
        TH1* njet_jet_phi = (TH1*)njet_jet_phi_hists.At(B);
        TH1* njet_jet_eta = (TH1*)njet_jet_eta_hists.At(B);
        TH1* njet_jet_area = (TH1*)njet_jet_area_hists.At(B);

        njet_jet_mass->Fill(njet_1jet[0].m());
        njet_jet_perp->Fill(njet_1jet[0].perp());
        njet_jet_phi->Fill(njet_1jet[0].phi());
        njet_jet_eta->Fill(njet_1jet[0].eta());
        njet_jet_area->Fill(njet_1jet_area[0].area());

        if (hardest_jet.size() == 1) {
          // njet_plugin_manual_1jet.setAxes(hardest_jet);
          // ClusterSequence njet_cluster_manual_1jet(fjInputs, njet_def_manual_1jet);
          // const NjettinessExtras *extras_manual_1jet = njettiness_extras(njet_cluster_manual_1jet);

          // TH1* tau1_diff = (TH1*)tau1_diff_hists.At(B);
          // tau1_diff->Fill(extras_manual_1jet->totalTau() - extras_1jet->totalTau());

          TH1* jetmass_diff = (TH1*)jetmass_diff_hists.At(B);
          TH1* jetperp_diff = (TH1*)jetperp_diff_hists.At(B);
          TH1* jetphi_diff = (TH1*)jetphi_diff_hists.At(B);
          TH1* njet_jet_areadiff = (TH1*)njet_jet_areadiff_hists.At(B);

          jetmass_diff->Fill(abs(hardest_jet[0].m() - njet_1jet[0].m()));
          jetperp_diff->Fill(hardest_jet[0].perp() - njet_1jet[0].perp());
          jetphi_diff->Fill(njet_1jet[0].delta_R(hardest_jet[0]));
          njet_jet_areadiff->Fill(hardest_jet[0].area() - njet_1jet_area[0].area());
        }
      }

      // vector<PseudoJet> exclusive_2jets_start = clustSeq.exclusive_jets(2);
      // vector<PseudoJet> exclusive_2jets = findMinAxes(fjInputs, exclusive_2jets_start, 2, beta, Rparam);
      // vector<PseudoJet> exclusive_2jets = clustSeq.exclusive_jets(2);

      // NjettinessPlugin njet_plugin_conical_twojet(2, axes_finder->def(), measure_function);
      // JetDefinition njet_def_conical_twojet(&njet_plugin_conical_twojet);
      // ClusterSequence njet_cluster_conical_twojet(fjInputs, njet_def_conical_twojet);
      // const NjettinessExtras *extras_twojet = njettiness_extras(njet_cluster_conical_twojet);
      // vector<PseudoJet> njet_conical_2jets = extras_twojet->jets();
      // vector<PseudoJet> njet_conical_2jets_axes = extras_twojet->axes();

      NjettinessPlugin njet_plugin_xcone_2jets(2, axes_finder_onepass->def(), measure_function_xcone);
      JetDefinition njet_def_xcone_2jets(&njet_plugin_xcone_2jets);
      // njet_plugin_xcone_2jets.setAxes(njet_conical_2jets_axes);
      ClusterSequence njet_cluster_xcone_2jets(fjInputs, njet_def_xcone_2jets);
      const NjettinessExtras *njet_cluster_xcone_2jets_extras = njettiness_extras(njet_cluster_xcone_2jets);
      vector<PseudoJet> njet_jets_2jets = njet_cluster_xcone_2jets.inclusive_jets();
      vector<PseudoJet> njet_2jets_axes = njet_cluster_xcone_2jets_extras->axes();

      NjettinessPlugin njet_plugin_xcone_2jets_area(2, Manual_Axes(), measure_function_xcone);
      JetDefinition njet_def_xcone_2jets_area(&njet_plugin_xcone_2jets_area);
      njet_plugin_xcone_2jets_area.setAxes(njet_2jets_axes);
      ClusterSequenceArea njet_cluster_xcone_2jets_area(fjInputs, njet_def_xcone_2jets_area, area_def);
      vector<PseudoJet> njet_2jets_area = njet_cluster_xcone_2jets_area.inclusive_jets();

      // cout << "tau2" << endl;

      PseudoJet big_jet_njet(0,0,0,0);
      TH1* higgs_invmass_njet = (TH1*)higgs_invmass_njet_hists.At(B);
      TH1* higgs_invmass_njet_combined = (TH1*)higgs_invmass_njet_combined_hists.At(B);
      TH1* jet_ptspectrum_njet = (TH1*)jet_ptspectrum_njet_hists.At(B);
      TH1* jet_phispectrum_njet = (TH1*)jet_phispectrum_njet_hists.At(B);
      TH1* jet_distance_njet = (TH1*)jet_distance_njet_hists.At(B);
      TH1* jet_distance_diff_njet = (TH1*)jet_distance_diff_njet_hists.At(B);
      TH1* twojet_njet_jet_area = (TH1*)twojet_njet_jet_area_hists.At(B);

      for (int i_jet = 0; i_jet < njet_jets_2jets.size(); i_jet++) {
        big_jet_njet = join(big_jet_njet, njet_jets_2jets[i_jet]);
        jet_ptspectrum_njet->Fill(njet_jets_2jets[i_jet].perp());
        twojet_njet_jet_area->Fill(njet_2jets_area[i_jet].area());
      }
      jet_distance_njet->Fill(njet_jets_2jets[0].delta_R(njet_jets_2jets[1]));
      jet_phispectrum_njet->Fill(abs(calcDphi(njet_jets_2jets[0], njet_jets_2jets[1])));
      double smaller_distance_njet = (njet_1jet[0].delta_R(njet_jets_2jets[0]) < njet_1jet[0].delta_R(njet_jets_2jets[1])) ? njet_1jet[0].delta_R(njet_jets_2jets[0]) : njet_1jet[0].delta_R(njet_jets_2jets[1]);
      if (smaller_distance_njet > epsilon) jet_distance_diff_njet->Fill(smaller_distance_njet);

      higgs_invmass_njet->Fill(big_jet_njet.m());

      // if (njet_1jet[0].perp() > merging_point && (njet_jets_2jets[0].delta_R(njet_jets_2jets[1]) > 2*Rparam)) {
        // higgs_invmass_njet_combined->Fill(njet_1jet[0].m());
      // }
      // else higgs_invmass_njet_combined->Fill(big_jet_njet.m());


      // vector<PseudoJet> exclusive_threejets_start = clustSeq.exclusive_jets(3);
      // vector<PseudoJet> exclusive_threejets = findMinAxes(fjInputs, exclusive_threejets_start, 3, beta, Rparam);
      // vector<PseudoJet> exclusive_threejets = clustSeq.exclusive_jets(2);
      // NjettinessPlugin njet_plugin_conical_3jets(3, axes_finder->def(), measure_function);
      // JetDefinition njet_def_conical_3jets(&njet_plugin_conical_3jets);
      // ClusterSequence njet_cluster_conical_3jets(fjInputs, njet_def_conical_3jets);
      // const NjettinessExtras *extras_3jets = njettiness_extras(njet_cluster_conical_3jets);
      // vector<PseudoJet> njet_conical_jets_3jets = extras_3jets->jets();
      // vector<PseudoJet> njet_conical_3jets_axes = extras_3jets->axes();

      NjettinessPlugin njet_plugin_xcone_3jets(3, axes_finder_onepass->def(), measure_function_xcone);
      JetDefinition njet_def_xcone_3jets(&njet_plugin_xcone_3jets);
      // njet_plugin_xcone_3jets.setAxes(njet_conical_3jets_axes);
      ClusterSequence njet_cluster_xcone_3jets(fjInputs, njet_def_xcone_3jets);
      vector<PseudoJet> njet_jets_3jets = njet_cluster_xcone_3jets.inclusive_jets();

      // cout << "tau3" << endl;

      PseudoJet threejet_minmassjet = findMinMass(njet_jets_3jets, zfrac);
      PseudoJet big_jet_njet_3jets(0,0,0,0);
      // for (int i_jet = 0; i_jet < njet_jets_3jets.size(); i_jet++) {
      //   big_jet_njet_3jets = join(big_jet_njet_3jets, njet_jets_3jets[i_jet]);
      // }

      TH1* threejet_invmass_njet = (TH1*)threejet_invmass_njet_hists.At(B);
      double mass_ratio = (double)big_jet_njet_3jets.m()/big_jet_njet.m();

      // if (mass_ratio > 0.8 && mass_ratio < 1.2) threejet_invmass_njet->Fill(big_jet_njet_3jets.m());
      // if (mass_ratio > 0.9 && mass_ratio < 1.1) threejet_invmass_njet->Fill(big_jet_njet.m());
      // if (big_jet_njet.m() > 100 && big_jet_njet.m() < 150) threejet_invmass_njet->Fill(big_jet_njet.m());
      // else threejet_invmass_njet->Fill(threejet_minmassjet.m());
      // else threejet_invmass_njet->Fill(threejet_bigjet.m());

      vector<PseudoJet> sorted_njet_jets_3jets = sorted_by_pt(njet_jets_3jets);
      PseudoJet third_jet_njet = sorted_njet_jets_3jets[sorted_njet_jets_3jets.size() - 1];
      PseudoJet closer_jet_njet = (third_jet_njet.delta_R(njet_jets_2jets[0]) < third_jet.delta_R(njet_jets_2jets[1])) ? njet_jets_2jets[0] : njet_jets_2jets[1];

      // if (threejet_minmassjet.m() < big_jet_njet.m()) threejet_invmass_njet->Fill(big_jet_njet_3jets.m());
      // else threejet_invmass_njet->Fill(threejet_minmassjet.m());

      threejet_invmass_njet->Fill(threejet_minmassjet.m());

      // if ((third_jet_njet.delta_R(closer_jet_njet) < 2*Rparam) && third_jet_njet.perp() > 50) threejet_invmass_njet->Fill(big_jet_njet_3jets.m());
      // else threejet_invmass_njet->Fill(threejet_minmassjet.m());

      TH1* threejet_thirdjet_phidiff_njet = (TH1*)threejet_thirdjet_phidiff_njet_hists.At(B);
      // TH1* threejet_thirdjet_area_njet = (TH1*)threejet_thirdjet_area_njet_hists.At(B);
      TH1* threejet_mass_2jet3jet_compare_njet = (TH1*)threejet_mass_2jet3jet_compare_njet_hists.At(B);

      threejet_thirdjet_phidiff_njet->Fill(third_jet_njet.delta_R(closer_jet_njet));
      // threejet_thirdjet_area_njet->Fill(third_jet_njet.area());
      // threejet_mass_2jet3jet_compare_njet->Fill(threejet_minmassjet.m(), third_jet_njet.delta_R(closer_jet_njet));
      threejet_mass_2jet3jet_compare_njet->Fill(big_jet_njet.m(), threejet_minmassjet.m());

      if ((threejet_minmassjet.m() > 100 && threejet_minmassjet.m() < 150) || (big_jet_njet.m() > 100 && big_jet_njet.m() < 150)) {
        higgs_events_2jets_combined[B]++;
      }

      // else big_jet_njet_3jets = big_jet_njet;
      //   // if (njet_jets_3jets[0].delta_R(njet_jets_3jets[1]) < 1.0 && njet_jets_3jets[1].delta_R(njet_jets_3jets[2]) < 1.0
      //   //   && njet_jets_3jets[0].delta_R(njet_jets_3jets[2]) < 1.0)
      //   //   big_jet_njet_3jets = join(big_jet_njet_3jets, njet_jets_3jets[i]);
      // // }

      // TH1* threejet_invmass_njet = (TH1*)threejet_invmass_njet_hists.At(B);
      // TH1* threejet_thirdjet_mass_njet = (TH1*)threejet_thirdjet_mass_njet_hists.At(B);
      // TH1* threejet_thirdjet_perp_njet = (TH1*)threejet_thirdjet_perp_njet_hists.At(B);

      // threejet_invmass_njet->Fill(big_jet_njet_3jets.m());
      // threejet_thirdjet_mass_njet->Fill(third_jet_njet.m());
      // threejet_thirdjet_perp_njet->Fill(third_jet_njet.perp());

      if ((B == 3) && threejet_minmassjet.m() < 100 && willDisplay) {

        TH2F *axes_display = new TH2F("axes_display", "Axes Plot", 300, -5, 5, 300, 0, 6.4);
        TH2F *axes_njet_display = new TH2F("axes_njet_display", "Axes Plot (Njettiness)", 300, -5, 5, 300, 0, 6.4);
        TH2F *akt_jets_display = new TH2F("akt_jets_display", "ak_{T} Jets", 300, -5, 5, 300, 0, 6.4);
        TH2F *njet_jets_display = new TH2F("njet_jets_display", "N-jettiness Jets", 300, -5, 5, 300, 0, 6.4);

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

        // for (int a = 0; a < hardest_twojets.size(); a++) {
        //   axes_display->Fill(hardest_twojets[a].eta(), hardest_twojets[a].phi());
        //   vector<PseudoJet> constituents = hardest_twojets[a].constituents();
        //   for (int i_const = 0; i_const < constituents.size(); i_const++) {
        //     akt_jets_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
        //   }          
        // }
        axes_display->SetMarkerStyle(3);
        axes_display->SetMarkerSize(3);
        axes_display->SetMarkerColor(kRed);
        akt_jets_display->SetMarkerStyle(21);
        akt_jets_display->SetMarkerSize(0.5);
        akt_jets_display->SetMarkerColor(kRed);
        // axes_display->Draw("SAMES");
        // akt_jets_display->Draw("SAMES");

        for (int a = 0; a < njet_jets_3jets.size(); a++) {
          axes_njet_display->Fill(njet_jets_3jets[a].eta(), njet_jets_3jets[a].phi());
          vector<PseudoJet> constituents = njet_jets_3jets[a].constituents();
          for (int i_const = 0; i_const < constituents.size(); i_const++) {
            njet_jets_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
          }          
        }
        axes_njet_display->SetMarkerStyle(3);
        axes_njet_display->SetMarkerSize(3);
        axes_njet_display->SetMarkerColor(kBlue);
        njet_jets_display->SetMarkerStyle(21);
        njet_jets_display->SetMarkerSize(0.5);
        njet_jets_display->SetMarkerColor(kBlue);
        axes_njet_display->Draw("SAMES");
        njet_jets_display->Draw("SAMES");

        display->Write();
        delete display;
        delete axes_display;
        delete axes_njet_display;
        delete njet_jets_display;
        delete akt_jets_display;

      }
    }
    delete event_display;

  }

  // pythia.stat();
  particle_mass->SetStats(0);
  particle_mass->Write();

  higgs_invmass_akt->Write();
  higgs_invmass_akt_combined->Write();
  jet_ptspectrum_akt->Write();
  jet_phispectrum_akt->Write();
  jet_distance_akt->Write();
  jet_distance_diff_akt->Write();

  TCanvas *mass_compare_2jets = new TCanvas("mass_compare_2jets", "mass_compare_2jets", 600, 600);
  mass_compare_2jets->cd();


  double higgs_invmass_akt_scale = 1/higgs_invmass_akt->Integral(0, 101);
  higgs_invmass_akt->Scale(higgs_invmass_akt_scale);
  higgs_invmass_akt->SetLineColor(kBlack);
  higgs_invmass_akt->SetLineStyle(7);
  higgs_invmass_akt->SetTitle("Dijet Mass (p_{T} > " + (TString)pythia_ptcut + ")");
  higgs_invmass_akt->GetXaxis()->SetTitle("m_{jj} (GeV)");
  higgs_invmass_akt->GetYaxis()->SetTitle("Relative Occurrence");
  higgs_invmass_akt->SetMinimum(0.0);
  double max_val = higgs_invmass_akt->GetMaximum();

  higgs_invmass_akt->Draw();
  TLegend *leg_mass_2jets = new TLegend(0.6, 0.55, 0.86, 0.86);
  leg_mass_2jets->SetFillColor(kWhite);
  leg_mass_2jets->SetLineColor(kWhite);
  // TPaveText *radius_text = new TPaveText(0.5, 0.82, 0.6, 0.88, "brNDC");
  TPaveText *radius_text = new TPaveText(0.6, 0.4, 0.86, 0.5, "brNDC");
  radius_text->SetTextFont(132);
  radius_text->SetTextSize(0.08);
  radius_text->SetFillColor(kWhite);
  radius_text->SetFillStyle(0);
  radius_text->SetLineColor(kWhite);
  radius_text->AddText("R = 0.5");
  double akt_efficiency = higgs_invmass_akt->Integral((double)100/250*higgs_invmass_akt->GetNbinsX(), (double)150/250*higgs_invmass_akt->GetNbinsX());
  higgs_efficiencies_2jets[4][i_perp] = akt_efficiency;

  // double akt_angulardiff = jet_phispectrum_akt->GetMean();
  // double akt_angulardiff = calcMedian1(jet_phispectrum_akt);
  double akt_angulardiff = calcMedian1(jet_distance_akt);
  double akt_angulardiff_width = jet_phispectrum_akt->GetRMS();
  angulardiff_2jets[4][i_perp] = akt_angulardiff;
  angulardiff_width_2jets[4][i_perp] = akt_angulardiff_width;
  angulardiff_1jet[4][i_perp] = jet_distance_diff_akt->GetMean();


  TH1* final_drawn_histogram;
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
    TH1* higgs_invmass_njet_combined = (TH1*)higgs_invmass_njet_combined_hists.At(B);
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
    higgs_invmass_njet_combined->Write();
    jet_ptspectrum_njet->Write();
    jet_phispectrum_njet->Write();
    jet_distance_njet->Write();
    jet_distance_diff_njet->Write();

    double higgs_invmass_njet_scale = 1/higgs_invmass_njet->Integral(0, 101);
    higgs_invmass_njet->Scale(higgs_invmass_akt_scale);


    higgs_invmass_njet->SetLineColor(colorlist[B]);
    // if (B == 0) higgs_invmass_njet->SetLineColor(kGreen);
    // if (B == 1) higgs_invmass_njet->SetLineColor(kYellow);
    // if (B == 2) higgs_invmass_njet->SetLineColor(kBlue);
    // if (B == 3) higgs_invmass_njet->SetLineColor(kRed);

    if (B == 2 || B == 3) {
      // higgs_invmass_njet->Draw("SAMES"); 
      leg_mass_2jets->AddEntry(higgs_invmass_njet, "#beta = " + (TString)ss.str(), "L");
    }
    if (B == 2) final_drawn_histogram = higgs_invmass_njet;
    if (B == 3) higgs_invmass_njet->Draw("SAMES");
    if (higgs_invmass_njet->GetMaximum() > max_val) max_val = higgs_invmass_njet->GetMaximum();
    double efficiency = higgs_invmass_njet->Integral((double)100/250*higgs_invmass_njet->GetNbinsX(), (double)150/250*higgs_invmass_njet->GetNbinsX());
    higgs_efficiencies_2jets[B][i_perp] = efficiency;

    // double njet_angulardiff = jet_phispectrum_njet->GetMean();
    // double njet_angulardiff = calcMedian1(jet_phispectrum_njet);
    double njet_angulardiff = calcMedian1(jet_distance_njet);
    double njet_angulardiff_width = jet_phispectrum_njet->GetRMS();
    angulardiff_2jets[B][i_perp] = njet_angulardiff;
    angulardiff_width_2jets[B][i_perp] = njet_angulardiff_width;
    angulardiff_1jet[B][i_perp] = jet_distance_diff_njet->GetMean();
    // angulardiff_1jet[B + 1][i_perp] = calcMedian1(jet_distance_diff_njet);
    angulardiff_prop_1jet[B][i_perp] = (double)jet_distance_diff_njet->Integral(0, (double)jet_distance_diff_njet->GetNbinsX()/10.0)/jet_distance_diff_njet->Integral(0, jet_distance_diff_njet->GetNbinsX() + 1);
  }
  // higgs_invmass_akt->SetMaximum(1.2*max_val);
  final_drawn_histogram->Draw("SAMES");
  leg_mass_2jets->AddEntry(higgs_invmass_akt, "ak_{T}", "L");
  higgs_invmass_akt->SetMaximum(0.25);
  leg_mass_2jets->Draw("SAMES");
  radius_text->Draw("SAMES");
  mass_compare_2jets->Write();
  mass_compare_2jets->Print("higgs2bbstudy_onepassXcone_test_plots/mass_compare_2jets_pt" + (TString)pythia_ptcut + ".eps", "eps");


  TCanvas *mass_compare_2jets_combined = new TCanvas("mass_compare_2jets_combined", "mass_compare_2jets_combined", 600, 600);
  mass_compare_2jets_combined->cd();

  // double higgs_invmass_akt_combined_scale = 1/higgs_invmass_akt_combined->Integral(0, 101);
  // higgs_invmass_akt_combined->Scale(higgs_invmass_akt_combined_scale);
  // higgs_invmass_akt_combined->SetLineColor(kBlack);
  // higgs_invmass_akt_combined->SetTitle("Invariant Mass (p_{T} > " + (TString)pythia_ptcut + ")");
  // higgs_invmass_akt_combined->GetXaxis()->SetTitle("m_{jj} (GeV)");
  // higgs_invmass_akt_combined->GetYaxis()->SetTitle("Relative Occurrence");
  // higgs_invmass_akt_combined->SetMinimum(0.0);
  // double max_val_combined = higgs_invmass_akt_combined->GetMaximum();

  // higgs_invmass_akt_combined->Draw();
  // TLegend *leg_mass_2jets_combined = new TLegend(0.65, 0.65, 0.88, 0.88);
  // leg_mass_2jets_combined->AddEntry(higgs_invmass_akt_combined, "ak_{T}", "L");
  // leg_mass_2jets_combined->SetFillColor(kWhite);
  // leg_mass_2jets_combined->SetLineColor(kWhite);
  // double akt_efficiency_combined = higgs_invmass_akt_combined->Integral(20, 30);
  // higgs_efficiencies_2jets_combined[0][i_perp] = akt_efficiency_combined;

  // for (int B = 0; B < n_betas; B++) {

  //   double beta = betalist[B];

  //   ostringstream ss;
  //   ss << beta;

  //   TH1* higgs_invmass_njet_combined = (TH1*)higgs_invmass_njet_combined_hists.At(B);

  //   higgs_invmass_njet_combined->Write();

  //   double higgs_invmass_njet_combined_scale = 1/higgs_invmass_njet_combined->Integral(0, 101);
  //   higgs_invmass_njet_combined->Scale(higgs_invmass_njet_combined_scale);

  //   if (B == 0) higgs_invmass_njet_combined->SetLineColor(kGreen);
  //   if (B == 1) higgs_invmass_njet_combined->SetLineColor(kYellow);
  //   if (B == 2) higgs_invmass_njet_combined->SetLineColor(kRed);
  //   if (B == 3) higgs_invmass_njet_combined->SetLineColor(kBlue);

  //   if (B == 2 || B == 3) {
  //     higgs_invmass_njet_combined->Draw("SAMES");
  //     leg_mass_2jets_combined->AddEntry(higgs_invmass_njet_combined, "#beta = " + (TString)ss.str(), "L");
  //   }
  //   if (higgs_invmass_njet_combined->GetMaximum() > max_val_combined) max_val_combined = higgs_invmass_njet_combined->GetMaximum();
  //   double efficiency_combined = higgs_invmass_njet_combined->Integral(20, 30);
  //   higgs_efficiencies_2jets_combined[B + 1][i_perp] = efficiency_combined;

  // }

  // higgs_invmass_akt_combined->SetMaximum(1.2*max_val_combined);
  // leg_mass_2jets_combined->Draw("SAMES");
  // radius_text->Draw("SAMES");
  // mass_compare_2jets_combined->Write();
  // mass_compare_2jets_combined->Print("higgs2bbstudy_onepassXcone_test_plots/mass_compare_2jets_combined_pt" + (TString)pythia_ptcut + ".eps", "eps");

  for (int i = 0; i < n_betas + 1; i++) {
    higgs_efficiencies_2jets_combined[i][i_perp] = (double)higgs_events_2jets_combined[i]/nEvents;
  }

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

    jetphi_diff_allpt_allbeta.Add(jetphi_diff);

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
    jetphi_diff->SetTitle("1-jet and hardest AKT Distance(p_{T} > " + (TString)pythia_ptcut + ")");
    jetphi_diff->GetXaxis()->SetTitle("#Delta R");
    jetphi_diff->GetYaxis()->SetTitle("Relative Occurrence");

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
  phi_compare_1jet->Print("higgs2bbstudy_onepassXcone_test_plots/phi_compare_1jet_pt" + (TString)pythia_ptcut + ".eps", "eps");


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

    jetperp_diff->SetTitle("1-jet and hardest AKT p_{T} difference (p_{T} > " + (TString)pythia_ptcut + ")");
    jetperp_diff->GetXaxis()->SetTitle("p_{T, akt} - p_{T, njet}");
    jetperp_diff->GetYaxis()->SetTitle("Relative Occurrence");

    jetperp_diff->SetLineColor(colorlist[B]);


    // if (B == 0) jetperp_diff->SetLineColor(kGreen);
    // if (B == 1) jetperp_diff->SetLineColor(kYellow);
    // if (B == 2) jetperp_diff->SetLineColor(kRed);
    // if (B == 3) jetperp_diff->SetLineColor(kBlue);

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
  perp_compare_1jet->Print("higgs2bbstudy_onepassXcone_test_plots/perp_compare_1jet_pt" + (TString)pythia_ptcut + ".eps", "eps");

  TCanvas *mass_compare_1jet = new TCanvas("mass_compare_1jet", "mass_compare_1jet", 600, 600);
  mass_compare_1jet->cd();

  double akt_jet_mass_scale = 1/akt_jet_mass->Integral(0, akt_jet_mass->GetNbinsX() + 1);
  akt_jet_mass->Scale(akt_jet_mass_scale);
  akt_jet_mass->SetLineColor(kBlack);
  akt_jet_mass->SetLineStyle(7);
  akt_jet_mass->SetTitle("Jet Mass (p_{T} > " + (TString)pythia_ptcut + ")");
  akt_jet_mass->GetXaxis()->SetTitle("m_{j} (GeV)");
  akt_jet_mass->GetYaxis()->SetTitle("Relative Occurrence");
  akt_jet_mass->SetMinimum(0.0);
  double max_val_1jet = akt_jet_mass->GetMaximum();

  akt_jet_mass->Draw();
  TLegend *leg_mass_1jet = new TLegend(0.6, 0.55, 0.86, 0.86);
  leg_mass_1jet->SetFillColor(kWhite);
  leg_mass_1jet->SetLineColor(kWhite);
  double akt_efficiency_1jet = akt_jet_mass->Integral((double)100/250*akt_jet_mass->GetNbinsX(), (double)150/250*akt_jet_mass->GetNbinsX());
  higgs_efficiencies_1jet[4][i_perp] = akt_efficiency_1jet;

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1* njet_jet_mass = (TH1*)njet_jet_mass_hists.At(B);

    double njet_jet_mass_scale = 1/njet_jet_mass->Integral(0, njet_jet_mass->GetNbinsX() + 1);
    njet_jet_mass->Scale(akt_jet_mass_scale);

    njet_jet_mass->SetLineColor(colorlist[B]);


    // if (B == 0) njet_jet_mass->SetLineColor(kGreen);
    // if (B == 1) njet_jet_mass->SetLineColor(kYellow);
    // if (B == 3) njet_jet_mass->SetLineColor(kRed);
    // if (B == 2) njet_jet_mass->SetLineColor(kBlue);

    if (B == 2 || B == 3) {
      // njet_jet_mass->Draw("SAMES");
      leg_mass_1jet->AddEntry(njet_jet_mass, "#beta = " + (TString)ss.str(), "L");
    }
    if (B == 2) final_drawn_histogram = njet_jet_mass;
    if (B == 3) njet_jet_mass->Draw("SAMES");
    if (njet_jet_mass->GetMaximum() > max_val_1jet) max_val_1jet = njet_jet_mass->GetMaximum();
    double efficiency = njet_jet_mass->Integral((double)100/250*njet_jet_mass->GetNbinsX(), (double)150/250*njet_jet_mass->GetNbinsX());
    higgs_efficiencies_1jet[B][i_perp] = efficiency;

  }
  final_drawn_histogram->Draw("SAMES");

  leg_mass_1jet->AddEntry(akt_jet_mass, "ak_{T}", "L");

  // akt_jet_mass->SetMaximum(1.2*max_val_1jet);
  akt_jet_mass->SetMaximum(0.3);
  leg_mass_1jet->Draw("SAMES");
  radius_text->Draw("SAMES");
  mass_compare_1jet->Write();
  mass_compare_1jet->Print("higgs2bbstudy_onepassXcone_test_plots/mass_compare_1jet_pt" + (TString)pythia_ptcut + ".eps", "eps");

  TCanvas *mass_compare_3jets = new TCanvas("mass_compare_3jets", "mass_compare_3jets", 600, 600);
  mass_compare_3jets->cd();

  // double threejet_invmass_akt_scale = 1/threejet_invmass_akt->Integral(0, 101);
  double threejet_invmass_akt_scale = (double)1/nEvents;
  threejet_invmass_akt->Scale(threejet_invmass_akt_scale);
  threejet_invmass_akt->SetLineColor(kBlack);
  threejet_invmass_akt->SetTitle("3-jet Min Mass (p_{T} > " + (TString)pythia_ptcut + ")");
  threejet_invmass_akt->GetXaxis()->SetTitle("m_{jjj} (GeV)");
  threejet_invmass_akt->GetYaxis()->SetTitle("Relative Occurrence");
  threejet_invmass_akt->SetMinimum(0.0);
  double max_val_mass_3jets = threejet_invmass_akt->GetMaximum();

  threejet_invmass_akt->Draw();
  TLegend *leg_mass_3jets = new TLegend(0.6, 0.55, 0.86, 0.86);
  leg_mass_3jets->SetFillColor(kWhite);
  leg_mass_3jets->SetLineColor(kWhite);

  double akt_efficiency_3jets = threejet_invmass_akt->Integral((double)100/250*threejet_invmass_akt->GetNbinsX(), (double)150/250*threejet_invmass_akt->GetNbinsX());
  higgs_efficiencies_3jets[4][i_perp] = akt_efficiency_3jets;
  higgs_efficiencies_32ratio[4][i_perp] = (double)akt_efficiency_3jets/akt_efficiency;


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

    // double threejet_invmass_njet_scale = 1/threejet_invmass_njet->Integral(0, 101);
    double threejet_invmass_njet_scale = (double)1/nEvents;
    threejet_invmass_njet->Scale(threejet_invmass_njet_scale);
    double higgs_invmass_njet_scale = (double)1/nEvents;
    // double higgs_invmass_njet_scale = 1/higgs_invmass_njet->Integral(0, 101);
    higgs_invmass_njet->Scale(higgs_invmass_njet_scale);

    threejet_invmass_njet->SetLineColor(colorlist[B]);


    // if (B == 0) threejet_invmass_njet->SetLineColor(kGreen);
    // if (B == 1) threejet_invmass_njet->SetLineColor(kYellow);
    // if (B == 2) threejet_invmass_njet->SetLineColor(kRed);
    // if (B == 3) threejet_invmass_njet->SetLineColor(kBlue);

    if (B == 2 || B == 3) {
      threejet_invmass_njet->Draw("SAMES");
      leg_mass_3jets->AddEntry(threejet_invmass_njet, "#beta = " + (TString)ss.str(), "L");
    }
    if (threejet_invmass_njet->GetMaximum() > max_val_mass_3jets) max_val_mass_3jets = threejet_invmass_njet->GetMaximum();

    double twojet_efficiency = higgs_invmass_njet->Integral((double)100/250*higgs_invmass_njet->GetNbinsX(), (double)150/250*higgs_invmass_njet->GetNbinsX());
    double efficiency = threejet_invmass_njet->Integral((double)100/250*threejet_invmass_njet->GetNbinsX(), (double)150/250*threejet_invmass_njet->GetNbinsX());
    higgs_efficiencies_3jets[B][i_perp] = efficiency;
    higgs_efficiencies_32ratio[B][i_perp] = (double)efficiency/twojet_efficiency;

  }
  leg_mass_3jets->AddEntry(threejet_invmass_akt, "ak_{T}", "L");

  threejet_invmass_akt->SetMaximum(1.2*max_val_mass_3jets);
  leg_mass_3jets->Draw("SAMES");
  mass_compare_3jets->Write();
  mass_compare_3jets->Print("higgs2bbstudy_onepassXcone_test_plots/mass_compare_3jets_pt" + (TString)pythia_ptcut + ".eps", "eps");


  TCanvas *phidiff_compare_3jets = new TCanvas("phidiff_compare_3jets", "phidiff_compare_3jets", 600, 600);
  phidiff_compare_3jets->cd();

  // double threejet_thirdjet_phidiff_akt_scale = 1/threejet_thirdjet_phidiff_akt->Integral(0, 101);
  double threejet_thirdjet_phidiff_akt_scale = (double)1/nEvents;
  threejet_thirdjet_phidiff_akt->Scale(threejet_thirdjet_phidiff_akt_scale);
  threejet_thirdjet_phidiff_akt->SetLineColor(kBlack);
  threejet_thirdjet_phidiff_akt->SetTitle("#Delta #phi");
  threejet_thirdjet_phidiff_akt->GetXaxis()->SetTitle("#phi");
  threejet_thirdjet_phidiff_akt->GetYaxis()->SetTitle("Relative Occurrence");
  threejet_thirdjet_phidiff_akt->SetMinimum(0.0);
  double max_val_phidiff_3jets = threejet_thirdjet_phidiff_akt->GetMaximum();

  threejet_thirdjet_phidiff_akt->Draw();
  TLegend *leg_phidiff_3jets = new TLegend(0.65, 0.65, 0.88, 0.88);
  leg_phidiff_3jets->SetFillColor(kWhite);
  leg_phidiff_3jets->SetLineColor(kWhite);

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1* threejet_thirdjet_phidiff_njet = (TH1*)threejet_thirdjet_phidiff_njet_hists.At(B);
    double threejet_thirdjet_phidiff_njet_scale = (double)1/nEvents;
    threejet_thirdjet_phidiff_njet->Scale(threejet_thirdjet_phidiff_njet_scale);

    threejet_thirdjet_phidiff_njet->Write();

    threejet_thirdjet_phidiff_njet->SetLineColor(colorlist[B]);


    // if (B == 0) threejet_thirdjet_phidiff_njet->SetLineColor(kGreen);
    // if (B == 1) threejet_thirdjet_phidiff_njet->SetLineColor(kYellow);
    // if (B == 3) threejet_thirdjet_phidiff_njet->SetLineColor(kRed);
    // if (B == 2) threejet_thirdjet_phidiff_njet->SetLineColor(kBlue);

    if (B == 2 || B == 3) {
      threejet_thirdjet_phidiff_njet->Draw("SAMES");
      leg_phidiff_3jets->AddEntry(threejet_thirdjet_phidiff_njet, "#beta = " + (TString)ss.str(), "L");
    }
    if (threejet_thirdjet_phidiff_njet->GetMaximum() > max_val_phidiff_3jets) max_val_phidiff_3jets = threejet_thirdjet_phidiff_njet->GetMaximum();
  }
  leg_phidiff_3jets->AddEntry(threejet_thirdjet_phidiff_akt, "ak_{T}", "L");

  threejet_thirdjet_phidiff_akt->SetMaximum(1.2*max_val_phidiff_3jets);
  leg_phidiff_3jets->Draw("SAMES");
  phidiff_compare_3jets->Write();
  phidiff_compare_3jets->Print("higgs2bbstudy_onepassXcone_test_plots/phidiff_compare_3jets_pt" + (TString)pythia_ptcut + ".eps", "eps");


  TCanvas *area_compare_1jet = new TCanvas("area_compare_1jet", "area_compare_1jet", 600, 600);
  area_compare_1jet->cd();

  // double threejet_thirdjet_area_akt_scale = 1/threejet_thirdjet_area_akt->Integral(0, 101);
  double akt_jet_area_scale = (double)1/akt_jet_area->Integral(0, akt_jet_area->GetNbinsX() + 1);
  akt_jet_area->Scale(akt_jet_area_scale);
  akt_jet_area->SetLineColor(kBlack);
  akt_jet_area->SetLineStyle(7);
  akt_jet_area->SetTitle("Area of Jet (p_{T} > " + (TString)pythia_ptcut + ")");
  akt_jet_area->GetXaxis()->SetTitle("Area");
  akt_jet_area->GetYaxis()->SetTitle("Relative Occurrence");
  akt_jet_area->SetMinimum(0.0);
  double max_val_area_1jet = akt_jet_area->GetMaximum();

  akt_jet_area->Draw();
  TLegend *leg_area_1jet = new TLegend(0.6, 0.55, 0.86, 0.86);
  leg_area_1jet->SetFillColor(kWhite);
  leg_area_1jet->SetLineColor(kWhite);

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1* njet_jet_area = (TH1*)njet_jet_area_hists.At(B);
    double njet_jet_area_scale = (double)1/njet_jet_area->Integral(0, njet_jet_area->GetNbinsX() + 1);
    njet_jet_area->Scale(njet_jet_area_scale);

    njet_jet_area->Write();

    njet_jet_area->SetLineColor(colorlist[B]);

    // if (B == 0) njet_jet_area->SetLineColor(kGreen);
    // if (B == 1) njet_jet_area->SetLineColor(kYellow);
    // if (B == 3) njet_jet_area->SetLineColor(kRed);
    // if (B == 2) njet_jet_area->SetLineColor(kBlue);

    if (B == 2 || B == 3) {
      // njet_jet_area->Draw("SAMES");
      leg_area_1jet->AddEntry(njet_jet_area, "#beta = " + (TString)ss.str(), "L");
    }
    if (B == 2) final_drawn_histogram = njet_jet_area;
    if (B == 3) njet_jet_area->Draw("SAMES");
    if (njet_jet_area->GetMaximum() > max_val_area_1jet) max_val_area_1jet = njet_jet_area->GetMaximum();
  }
  leg_area_1jet->AddEntry(akt_jet_area, "ak_{T}", "L");
  final_drawn_histogram->Draw("SAMES");

  akt_jet_area->SetMaximum(1.2*max_val_area_1jet);
  leg_area_1jet->Draw("SAMES");
  radius_text->Draw("SAMES");
  TLine *true_area = new TLine();
  true_area->SetX1(TMath::Pi()*Rparam*Rparam);
  true_area->SetVertical(true);
  true_area->SetY1(0);
  true_area->SetY2(1.2*max_val_area_1jet);
  true_area->SetLineColor(kRed);
  true_area->SetLineWidth(4);
  true_area->SetLineStyle(7);
  true_area->Draw();

  area_compare_1jet->Write();
  area_compare_1jet->Print("higgs2bbstudy_onepassXcone_test_plots/area_compare_1jet_pt" + (TString)pythia_ptcut + ".eps", "eps");


  TCanvas *area_compare_2jets = new TCanvas("area_compare_2jets", "area_compare_2jets", 600, 600);
  area_compare_2jets->cd();

  // double threejet_thirdjet_area_akt_scale = 1/threejet_thirdjet_area_akt->Integral(0, 101);
  double twojet_akt_jet_area_scale = (double)1/twojet_akt_jet_area->Integral(0, twojet_akt_jet_area->GetNbinsX() + 1);
  twojet_akt_jet_area->Scale(twojet_akt_jet_area_scale);
  twojet_akt_jet_area->SetLineColor(kBlack);
  twojet_akt_jet_area->SetLineStyle(7);
  twojet_akt_jet_area->SetTitle("Area of Jet (p_{T} > " + (TString)pythia_ptcut + ")");
  twojet_akt_jet_area->GetXaxis()->SetTitle("Area");
  twojet_akt_jet_area->GetYaxis()->SetTitle("Relative Occurrence");
  twojet_akt_jet_area->SetMinimum(0.0);
  double max_val_area_2jets = twojet_akt_jet_area->GetMaximum();

  twojet_akt_jet_area->Draw();
  TLegend *leg_area_2jets = new TLegend(0.6, 0.55, 0.86, 0.86);
  leg_area_2jets->SetFillColor(kWhite);
  leg_area_2jets->SetLineColor(kWhite);

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1* twojet_njet_jet_area = (TH1*)twojet_njet_jet_area_hists.At(B);
    double twojet_njet_jet_area_scale = (double)1/twojet_njet_jet_area->Integral(0, twojet_njet_jet_area->GetNbinsX() + 1);
    twojet_njet_jet_area->Scale(twojet_njet_jet_area_scale);

    twojet_njet_jet_area->Write();

    twojet_njet_jet_area->SetLineColor(colorlist[B]);

    // if (B == 0) twojet_njet_jet_area->SetLineColor(kGreen);
    // if (B == 1) twojet_njet_jet_area->SetLineColor(kYellow);
    // if (B == 3) twojet_njet_jet_area->SetLineColor(kRed);
    // if (B == 2) twojet_njet_jet_area->SetLineColor(kBlue);

    if (B == 2 || B == 3) {
      // twojet_njet_jet_area->Draw("SAMES");
      leg_area_2jets->AddEntry(twojet_njet_jet_area, "#beta = " + (TString)ss.str(), "L");
    }
    if (B == 2) final_drawn_histogram = twojet_njet_jet_area;
    if (B == 3) twojet_njet_jet_area->Draw("SAMES");
    if (twojet_njet_jet_area->GetMaximum() > max_val_area_2jets) max_val_area_2jets = twojet_njet_jet_area->GetMaximum();
  }
  leg_area_2jets->AddEntry(twojet_akt_jet_area, "ak_{T}", "L");
  final_drawn_histogram->Draw("SAMES");

  twojet_akt_jet_area->SetMaximum(1.2*max_val_area_2jets);
  leg_area_2jets->Draw("SAMES");
  radius_text->Draw("SAMES");
  TLine *true_area_2 = new TLine();
  true_area_2->SetX1(TMath::Pi()*Rparam*Rparam);
  true_area_2->SetVertical(true);
  true_area_2->SetY1(0);
  true_area_2->SetY2(1.2*max_val_area_2jets);
  true_area_2->SetLineColor(kRed);
  true_area_2->SetLineWidth(4);
  true_area_2->SetLineStyle(7);
  true_area_2->Draw();

  area_compare_2jets->Write();
  area_compare_2jets->Print("higgs2bbstudy_onepassXcone_test_plots/area_compare_2jets_pt" + (TString)pythia_ptcut + ".eps", "eps");



  TCanvas *areadiff_compare_1jet = new TCanvas("areadiff_compare_1jet", "areadiff_compare_1jet", 800, 800);
  areadiff_compare_1jet->cd();
  areadiff_compare_1jet->SetLogy();
  double max_val_areadiff_1jet = 0;
  TLegend *leg_areadiff_1jet = new TLegend(0.7, 0.7, 0.88, 0.88);
  leg_areadiff_1jet->SetFillColor(kWhite);
  leg_areadiff_1jet->SetLineColor(kWhite);

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* njet_jet_areadiff = (TH1*)njet_jet_areadiff_hists.At(B);
    double njet_jet_areadiff_scale = 1/njet_jet_areadiff->Integral(0, njet_jet_areadiff->GetNbinsX() + 1);
    njet_jet_areadiff->Scale(njet_jet_areadiff_scale);

    njet_jet_areadiff->SetLineColor(colorlist[B]);

    // if (B == 0) njet_jet_areadiff->SetLineColor(kRed);
    // if (B == 1) njet_jet_areadiff->SetLineColor(kBlue);
    // if (B == 3) njet_jet_areadiff->SetLineColor(kRed);
    // if (B == 2) njet_jet_areadiff->SetLineColor(kBlue);

    if (B == 2) {
      // njet_jet_areadiff->Draw();
      final_drawn_histogram = njet_jet_areadiff;
      leg_areadiff_1jet->AddEntry(njet_jet_areadiff, "#beta = " + (TString)ss.str());
    }
    else if (B == 3) {
      njet_jet_areadiff->Draw();
      leg_areadiff_1jet->AddEntry(njet_jet_areadiff, "#beta = " + (TString)ss.str());
    }
    njet_jet_areadiff->SetTitle("Jet Area Difference (p_{T} > " + (TString)pythia_ptcut + ")");
    njet_jet_areadiff->GetXaxis()->SetTitle("A_{akt} - A_{XCone}");
    njet_jet_areadiff->GetYaxis()->SetTitle("Relative Occurrence");

    if (njet_jet_areadiff->GetMaximum() > max_val_areadiff_1jet) {
      max_val_areadiff_1jet = njet_jet_areadiff->GetMaximum();
      njet_jet_areadiff->SetMaximum(1.5*max_val_areadiff_1jet);
    }
  }
  final_drawn_histogram->Draw("SAMES");
  leg_areadiff_1jet->Draw("SAMES");
  areadiff_compare_1jet->Write();
  areadiff_compare_1jet->Print("higgs2bbstudy_onepassXcone_test_plots/areadiff_compare_1jet_pt" + (TString)pythia_ptcut + ".eps", "eps");

  TCanvas *area_compare_3jets = new TCanvas("area_compare_3jets", "area_compare_3jets", 600, 600);
  area_compare_3jets->cd();

  // double threejet_thirdjet_area_akt_scale = 1/threejet_thirdjet_area_akt->Integral(0, 101);
  double threejet_thirdjet_area_akt_scale = (double)1/nEvents;
  threejet_thirdjet_area_akt->Scale(threejet_thirdjet_area_akt_scale);
  threejet_thirdjet_area_akt->SetLineColor(kBlack);
  threejet_thirdjet_area_akt->SetTitle("Area of Third Jet");
  threejet_thirdjet_area_akt->GetXaxis()->SetTitle("A");
  threejet_thirdjet_area_akt->GetYaxis()->SetTitle("Relative Occurrence");
  threejet_thirdjet_area_akt->SetMinimum(0.0);
  double max_val_area_3jets = threejet_thirdjet_area_akt->GetMaximum();

  threejet_thirdjet_area_akt->Draw();
  TLegend *leg_area_3jets = new TLegend(0.65, 0.65, 0.88, 0.88);
  leg_area_3jets->SetFillColor(kWhite);
  leg_area_3jets->SetLineColor(kWhite);

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1* threejet_thirdjet_area_njet = (TH1*)threejet_thirdjet_area_njet_hists.At(B);
    double threejet_thirdjet_area_njet_scale = (double)1/nEvents;
    threejet_thirdjet_area_njet->Scale(threejet_thirdjet_area_njet_scale);

    threejet_thirdjet_area_njet->Write();

    threejet_thirdjet_area_njet->SetLineColor(colorlist[B]);

    // if (B == 0) threejet_thirdjet_area_njet->SetLineColor(kGreen);
    // if (B == 1) threejet_thirdjet_area_njet->SetLineColor(kYellow);
    // if (B == 2) threejet_thirdjet_area_njet->SetLineColor(kRed);
    // if (B == 3) threejet_thirdjet_area_njet->SetLineColor(kBlue);

    if (B == 2 || B == 3) {
      threejet_thirdjet_area_njet->Draw("SAMES");
      leg_area_3jets->AddEntry(threejet_thirdjet_area_njet, "#beta = " + (TString)ss.str(), "L");
    }
    if (threejet_thirdjet_area_njet->GetMaximum() > max_val_area_3jets) max_val_area_3jets = threejet_thirdjet_area_njet->GetMaximum();
  }
  leg_area_3jets->AddEntry(threejet_thirdjet_area_akt, "ak_{T}", "L");

  threejet_thirdjet_area_akt->SetMaximum(1.2*max_val_area_3jets);
  leg_area_3jets->Draw("SAMES");
  area_compare_3jets->Write();
  area_compare_3jets->Print("higgs2bbstudy_onepassXcone_test_plots/area_compare_3jets_pt" + (TString)pythia_ptcut + ".eps", "eps");


  TCanvas *mass_2jet3jet_compare_3jets = new TCanvas("mass_2jet3jet_compare_3jets", "mass_2jet3jet_compare_3jets", 600, 600);
  mass_2jet3jet_compare_3jets->cd();

  // double threejet_mass_2jet3jet_compare_akt_scale = 1/threejet_mass_2jet3jet_compare_akt->Integral(0, 101);
  double threejet_mass_2jet3jet_compare_akt_scale = (double)1/nEvents;
  threejet_mass_2jet3jet_compare_akt->Scale(threejet_mass_2jet3jet_compare_akt_scale);
  threejet_mass_2jet3jet_compare_akt->SetMarkerColor(kBlack);
  threejet_mass_2jet3jet_compare_akt->SetTitle("2-jet/3-jet min mass comparison (p_{T} > " + (TString)pythia_ptcut + ")");
  threejet_mass_2jet3jet_compare_akt->GetXaxis()->SetTitle("2-jet mass");
  threejet_mass_2jet3jet_compare_akt->GetYaxis()->SetTitle("3-jet min mass");
  threejet_mass_2jet3jet_compare_akt->SetMinimum(0.0);
  double max_val_mass_2jet3jet_3jets = threejet_mass_2jet3jet_compare_akt->GetMaximum();

  // for (int i = 1; i <= threejet_mass_2jet3jet_compare_akt->GetNbinsX(); i++) {
  //   for (int j = 1; j <= threejet_mass_2jet3jet_compare_akt->GetNbinsY(); j++) {
  //     double log_currentcontent;
  //     double currentcontent = threejet_mass_2jet3jet_compare_akt->GetCellContent(i,j);
  //     if (currentcontent == 0) log_currentcontent = 0;
  //     else log_currentcontent = TMath::Log10(currentcontent);
      
  //     threejet_mass_2jet3jet_compare_akt->SetCellContent(i, j, log_currentcontent);
  //   }
  // }

  threejet_mass_2jet3jet_compare_akt->Draw("box");


  TLegend *leg_mass_2jet3jet_3jets = new TLegend(0.40, 0.65, 0.65, 0.88);
  leg_mass_2jet3jet_3jets->SetFillColor(kWhite);
  leg_mass_2jet3jet_3jets->SetLineColor(kWhite);

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];

    ostringstream ss;
    ss << beta;

    TH1* threejet_mass_2jet3jet_compare_njet = (TH1*)threejet_mass_2jet3jet_compare_njet_hists.At(B);
    double threejet_mass_2jet3jet_compare_njet_scale = (double)1/nEvents;
    threejet_mass_2jet3jet_compare_njet->Scale(threejet_mass_2jet3jet_compare_njet_scale);

    threejet_mass_2jet3jet_compare_njet->Write();

    // for (int i = 1; i <= threejet_mass_2jet3jet_compare_njet->GetNbinsX(); i++) {
    //   for (int j = 1; j <= threejet_mass_2jet3jet_compare_njet->GetNbinsY(); j++) {
    //     double log_currentcontent;
    //     double currentcontent = threejet_mass_2jet3jet_compare_njet->GetCellContent(i,j);
    //     if (currentcontent == 0) log_currentcontent = 0;
    //     else log_currentcontent = TMath::Log10(currentcontent);
    //     threejet_mass_2jet3jet_compare_njet->SetCellContent(i, j, log_currentcontent);
    //   }
    // }

    threejet_mass_2jet3jet_compare_njet->SetMarkerColor(colorlist[B]);
    threejet_mass_2jet3jet_compare_njet->SetLineColor(colorlist[B]);

    // if (B == 0) {
    //   threejet_mass_2jet3jet_compare_njet->SetMarkerColor(kGreen);
    //   threejet_mass_2jet3jet_compare_njet->SetLineColor(kGreen);
    // }
    // if (B == 1) {
    //   threejet_mass_2jet3jet_compare_njet->SetMarkerColor(kYellow);
    //   threejet_mass_2jet3jet_compare_njet->SetLineColor(kYellow);
    // }
    // if (B == 3) {
    //   threejet_mass_2jet3jet_compare_njet->SetMarkerColor(kRed);
    //   threejet_mass_2jet3jet_compare_njet->SetLineColor(kRed);
    // }
    // if (B == 2) {
    //   threejet_mass_2jet3jet_compare_njet->SetMarkerColor(kBlue);
    //   threejet_mass_2jet3jet_compare_njet->SetLineColor(kBlue);
    // }

    if (B == 2 || B == 3) {
      threejet_mass_2jet3jet_compare_njet->Draw("box SAMES");
      leg_mass_2jet3jet_3jets->AddEntry(threejet_mass_2jet3jet_compare_njet, "#beta = " + (TString)ss.str(), "L");
    }
    if (threejet_mass_2jet3jet_compare_njet->GetMaximum() > max_val_mass_2jet3jet_3jets) max_val_mass_2jet3jet_3jets = threejet_mass_2jet3jet_compare_njet->GetMaximum();
  }

  TLine *line1 = new TLine();
  line1->SetX1(100);
  line1->SetVertical(true);
  line1->SetY1(0);
  line1->SetY2(500);
  line1->SetLineColor(kRed);
  line1->SetLineWidth(4);
  line1->SetLineStyle(7);
  line1->Draw("SAMES");

  TLine *line2 = new TLine();
  line2->SetX1(150);
  line2->SetVertical(true);
  line2->SetY1(0);
  line2->SetY2(500);
  line2->SetLineColor(kRed);
  line2->SetLineWidth(4);
  line2->SetLineStyle(7);
  line2->Draw("SAMES");

  TLine *line3 = new TLine();
  line3->SetY1(100);
  line3->SetHorizontal(true);
  line3->SetX1(0);
  line3->SetX2(500);
  line3->SetLineColor(kRed);
  line3->SetLineWidth(4);
  line3->SetLineStyle(7);
  line3->Draw("SAMES");

  TLine *line4 = new TLine();
  line4->SetY1(150);
  line4->SetHorizontal(true);
  line4->SetX1(0);
  line4->SetX2(500);
  line4->SetLineColor(kRed);
  line4->SetLineWidth(4);
  line4->SetLineStyle(7);
  line4->Draw("SAMES");
  leg_mass_2jet3jet_3jets->AddEntry(threejet_mass_2jet3jet_compare_akt, "ak_{T}", "L");

  threejet_mass_2jet3jet_compare_akt->SetMaximum(1.2*max_val_mass_2jet3jet_3jets);
  leg_mass_2jet3jet_3jets->Draw("SAMES");
  mass_2jet3jet_compare_3jets->Write();
  mass_2jet3jet_compare_3jets->Print("higgs2bbstudy_onepassXcone_test_plots/mass_2jet3jet_compare_pt" + (TString)pythia_ptcut + ".eps", "eps");


  //Prevent memory leaks  
  delete jetDef;
  // delete mass_compare;
  // delete leg_mass;

}

TMultiGraph *higgs_efficiencies_1jet_biggraph = new TMultiGraph();
TMultiGraph *higgs_efficiencies_2jets_biggraph = new TMultiGraph();
TMultiGraph *higgs_efficiencies_3jets_biggraph = new TMultiGraph();
TMultiGraph *higgs_efficiencies_32ratio_biggraph = new TMultiGraph();
TMultiGraph *higgs_efficiencies_2jets_combined_biggraph = new TMultiGraph();
TMultiGraph *angulardiff_1jet_biggraph = new TMultiGraph();
TMultiGraph *angulardiff_prop_1jet_biggraph = new TMultiGraph();
TMultiGraph *angulardiff_2jets_biggraph = new TMultiGraph();
TMultiGraph *jet_distance_diff_1jet_biggraph = new TMultiGraph();

TLegend *leg_higgs_1jet = new TLegend(0.16, 0.55, 0.45, 0.86);
TLegend *leg_higgs_2jets = new TLegend(0.16, 0.16, 0.45, 0.45);
TLegend *leg_higgs_2jets_combined = new TLegend(0.16, 0.16, 0.45, 0.45);

leg_higgs_1jet->SetFillStyle(0);
leg_higgs_2jets->SetFillStyle(0);
leg_higgs_2jets_combined->SetFillStyle(0);

TLegend *leg_higgs_3jets = new TLegend(0.14, 0.68, 0.35, 0.88);
TLegend *leg_higgs_32ratio = new TLegend(0.14, 0.68, 0.35, 0.88);
TLegend *leg_angulardiff_1jet = new TLegend(0.6, 0.6, 0.88, 0.88);
TLegend *leg_angulardiff_prop_1jet = new TLegend(0.6, 0.6, 0.88, 0.88);
TLegend *leg_angulardiff_2jets = new TLegend(0.12, 0.6, 0.4, 0.88);
TLegend *leg_jet_distance_diff_1jet = new TLegend(0.12, 0.6, 0.4, 0.88);

TPaveText* radius_text = new TPaveText(0.64, 0.76, 0.89, 0.87, "brNDC");
radius_text->SetTextFont(132);
radius_text->SetTextSize(0.08);
radius_text->SetFillColor(kWhite);
radius_text->SetLineColor(kWhite);
radius_text->SetFillStyle(0);
radius_text->SetBorderSize(0);
radius_text->AddText("R = 0.5");

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
leg_higgs_2jets_combined->SetFillColor(kWhite);
leg_higgs_2jets_combined->SetLineWidth(0);
leg_higgs_2jets_combined->SetLineColor(kWhite);
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
  if (i != 4) beta = betalist[i];

  ostringstream ss;
  ss << beta;

  TGraph *higgs_efficiencies_1jet_graph = new TGraph(n_perps, perplist, higgs_efficiencies_1jet[i]);
  TGraph *higgs_efficiencies_2jets_graph = new TGraph(n_perps, perplist, higgs_efficiencies_2jets[i]);
  TGraph *higgs_efficiencies_3jets_graph = new TGraph(n_perps, perplist, higgs_efficiencies_3jets[i]);
  TGraph *higgs_efficiencies_32ratio_graph = new TGraph(n_perps, perplist, higgs_efficiencies_32ratio[i]);
  TGraph *higgs_efficiencies_2jets_combined_graph = new TGraph(n_perps, perplist, higgs_efficiencies_2jets_combined[i]);
  TGraph *angulardiff_1jet_graph = new TGraph(n_perps, perplist, angulardiff_1jet[i]);
  TGraph *angulardiff_prop_1jet_graph = new TGraph(n_perps, perplist, angulardiff_prop_1jet[i]);
  TGraph *angulardiff_2jets_graph = new TGraph(n_perps, perplist, angulardiff_2jets[i]);
  // TGraphAsymmErrors *jet_distance_diff_1jet_graph = new TGraphAsymmErrors(n_perps, perplist, jet_distance_1jet_median[i], 0, 0, jet_distance_1jet_firstquartile[i], jet_distance_1jet_thirdquartile[i]);
  TGraph *jet_distance_diff_1jet_graph = new TGraph(n_perps, perplist, jet_distance_1jet_median[i]);

  higgs_efficiencies_1jet_graph->SetLineWidth(4);
  higgs_efficiencies_2jets_graph->SetLineWidth(4);
  higgs_efficiencies_2jets_combined_graph->SetLineWidth(4);
  higgs_efficiencies_3jets_graph->SetLineWidth(4);

  // higgs_efficiencies_1jet_graph->SetMarkerStyle(3);
  // higgs_efficiencies_1jet_graph->SetMarkerSize(2);
  // higgs_efficiencies_2jets_graph->SetMarkerStyle(3);
  // higgs_efficiencies_2jets_graph->SetMarkerSize(2);
  // higgs_efficiencies_3jets_graph->SetMarkerStyle(3);
  // higgs_efficiencies_3jets_graph->SetMarkerSize(2);
  // higgs_efficiencies_32ratio_graph->SetMarkerStyle(3);
  // higgs_efficiencies_32ratio_graph->SetMarkerSize(2);
  // angulardiff_1jet_graph->SetMarkerStyle(3);
  // angulardiff_1jet_graph->SetMarkerSize(2);
  // angulardiff_prop_1jet_graph->SetMarkerStyle(3);
  // angulardiff_prop_1jet_graph->SetMarkerSize(2);
  // angulardiff_2jets_graph->SetMarkerStyle(3);
  // angulardiff_2jets_graph->SetMarkerSize(2);
  // jet_distance_diff_1jet_graph->SetMarkerStyle(3);
  // jet_distance_diff_1jet_graph->SetMarkerSize(2);

  if (i == 4) {
    higgs_efficiencies_1jet_graph->SetMarkerColor(kBlack);
    higgs_efficiencies_2jets_graph->SetMarkerColor(kBlack);
    higgs_efficiencies_3jets_graph->SetMarkerColor(kBlack);
    higgs_efficiencies_32ratio_graph->SetMarkerColor(kBlack);
    higgs_efficiencies_2jets_combined_graph->SetMarkerColor(kBlack);
    angulardiff_2jets_graph->SetMarkerColor(kBlack);
    angulardiff_2jets_graph->SetLineColor(kBlack);

    higgs_efficiencies_1jet_graph->SetLineStyle(7);
    higgs_efficiencies_2jets_graph->SetLineStyle(7);
    higgs_efficiencies_2jets_combined_graph->SetLineStyle(7);

    leg_higgs_1jet->AddEntry(higgs_efficiencies_1jet_graph, "ak_{T}", "l");
    leg_higgs_2jets->AddEntry(higgs_efficiencies_2jets_graph, "ak_{T}", "l");
    leg_higgs_3jets->AddEntry(higgs_efficiencies_3jets_graph, "ak_{T}", "l");
    leg_higgs_32ratio->AddEntry(higgs_efficiencies_32ratio_graph, "ak_{T}", "l");
    leg_higgs_2jets_combined->AddEntry(higgs_efficiencies_2jets_combined_graph, "ak_{T}", "l");
    leg_angulardiff_2jets->AddEntry(angulardiff_2jets_graph, "ak_{T}", "l");
  }
  if (i == 0) {
    higgs_efficiencies_1jet_graph->SetMarkerColor(kGreen);
    higgs_efficiencies_2jets_graph->SetMarkerColor(kGreen);
    higgs_efficiencies_3jets_graph->SetMarkerColor(kGreen);
    higgs_efficiencies_1jet_graph->SetLineColor(kGreen);
    higgs_efficiencies_2jets_graph->SetLineColor(kGreen);
    higgs_efficiencies_3jets_graph->SetLineColor(kGreen);
    higgs_efficiencies_32ratio_graph->SetMarkerColor(kGreen);
    higgs_efficiencies_32ratio_graph->SetLineColor(kGreen);
    higgs_efficiencies_2jets_combined_graph->SetMarkerColor(kGreen);
    higgs_efficiencies_2jets_combined_graph->SetLineColor(kGreen);
    angulardiff_1jet_graph->SetMarkerColor(kGreen);
    angulardiff_prop_1jet_graph->SetMarkerColor(kGreen);
    angulardiff_2jets_graph->SetMarkerColor(kGreen);
    angulardiff_2jets_graph->SetLineColor(kGreen);
    jet_distance_diff_1jet_graph->SetMarkerColor(kGreen);
    jet_distance_diff_1jet_graph->SetLineColor(kGreen);
  }
  if (i == 1) {
    higgs_efficiencies_1jet_graph->SetMarkerColor(kYellow);
    higgs_efficiencies_2jets_graph->SetMarkerColor(kYellow);
    higgs_efficiencies_3jets_graph->SetMarkerColor(kYellow);
    higgs_efficiencies_1jet_graph->SetLineColor(kYellow);
    higgs_efficiencies_2jets_graph->SetLineColor(kYellow);
    higgs_efficiencies_3jets_graph->SetLineColor(kYellow);
    higgs_efficiencies_32ratio_graph->SetMarkerColor(kYellow);
    higgs_efficiencies_32ratio_graph->SetLineColor(kYellow);
    higgs_efficiencies_2jets_combined_graph->SetMarkerColor(kYellow);
    higgs_efficiencies_2jets_combined_graph->SetLineColor(kYellow);
    angulardiff_1jet_graph->SetMarkerColor(kYellow);
    angulardiff_prop_1jet_graph->SetMarkerColor(kYellow);
    angulardiff_2jets_graph->SetMarkerColor(kYellow);
    angulardiff_2jets_graph->SetLineColor(kYellow);
    jet_distance_diff_1jet_graph->SetMarkerColor(kYellow);
    jet_distance_diff_1jet_graph->SetLineColor(kYellow);
  }
  if (i == 3) {
    higgs_efficiencies_1jet_graph->SetMarkerColor(kRed);
    higgs_efficiencies_2jets_graph->SetMarkerColor(kRed);
    higgs_efficiencies_3jets_graph->SetMarkerColor(kRed);
    higgs_efficiencies_1jet_graph->SetLineColor(kRed);
    higgs_efficiencies_2jets_graph->SetLineColor(kRed);
    higgs_efficiencies_3jets_graph->SetLineColor(kRed);
    higgs_efficiencies_32ratio_graph->SetMarkerColor(kRed);
    higgs_efficiencies_32ratio_graph->SetLineColor(kRed);
    higgs_efficiencies_2jets_combined_graph->SetMarkerColor(kRed);
    higgs_efficiencies_2jets_combined_graph->SetLineColor(kRed);
    angulardiff_1jet_graph->SetMarkerColor(kRed);
    angulardiff_prop_1jet_graph->SetMarkerColor(kRed);
    angulardiff_2jets_graph->SetMarkerColor(kRed);
    angulardiff_2jets_graph->SetLineColor(kRed);
    jet_distance_diff_1jet_graph->SetMarkerColor(kRed);
    jet_distance_diff_1jet_graph->SetLineColor(kRed);

  }
  if (i == 2) {
    higgs_efficiencies_1jet_graph->SetMarkerColor(kBlue);
    higgs_efficiencies_2jets_graph->SetMarkerColor(kBlue);
    higgs_efficiencies_3jets_graph->SetMarkerColor(kBlue);
    higgs_efficiencies_1jet_graph->SetLineColor(kBlue);
    higgs_efficiencies_2jets_graph->SetLineColor(kBlue);
    higgs_efficiencies_3jets_graph->SetLineColor(kBlue);
    higgs_efficiencies_32ratio_graph->SetMarkerColor(kBlue);
    higgs_efficiencies_32ratio_graph->SetLineColor(kBlue);
    higgs_efficiencies_2jets_combined_graph->SetMarkerColor(kBlue);
    higgs_efficiencies_2jets_combined_graph->SetLineColor(kBlue);
    angulardiff_1jet_graph->SetMarkerColor(kBlue);
    angulardiff_prop_1jet_graph->SetMarkerColor(kBlue);
    angulardiff_2jets_graph->SetMarkerColor(kBlue);
    angulardiff_2jets_graph->SetLineColor(kBlue);
    jet_distance_diff_1jet_graph->SetMarkerColor(kBlue);
    jet_distance_diff_1jet_graph->SetLineColor(kBlue);

  }

  if (i == 2 || i == 3 || i == 4) {
    higgs_efficiencies_1jet_biggraph->Add(higgs_efficiencies_1jet_graph);
    higgs_efficiencies_2jets_biggraph->Add(higgs_efficiencies_2jets_graph);
    higgs_efficiencies_3jets_biggraph->Add(higgs_efficiencies_3jets_graph);
    higgs_efficiencies_32ratio_biggraph->Add(higgs_efficiencies_32ratio_graph);
    higgs_efficiencies_2jets_combined_biggraph->Add(higgs_efficiencies_2jets_combined_graph);
    angulardiff_2jets_biggraph->Add(angulardiff_2jets_graph);
    angulardiff_1jet_biggraph->Add(angulardiff_1jet_graph);
    angulardiff_prop_1jet_biggraph->Add(angulardiff_prop_1jet_graph);
    if (i != 4) {
      jet_distance_diff_1jet_biggraph->Add(jet_distance_diff_1jet_graph);
      leg_higgs_1jet->AddEntry(higgs_efficiencies_1jet_graph, "#beta = " + (TString)ss.str(), "l");
      leg_higgs_2jets->AddEntry(higgs_efficiencies_2jets_graph, "#beta = " + (TString)ss.str(), "l");
      leg_higgs_3jets->AddEntry(higgs_efficiencies_3jets_graph, "#beta = " + (TString)ss.str(), "l");
      leg_higgs_32ratio->AddEntry(higgs_efficiencies_32ratio_graph, "#beta = " + (TString)ss.str(), "l");
      leg_higgs_2jets_combined->AddEntry(higgs_efficiencies_2jets_combined_graph, "#beta = " + (TString)ss.str(), "l");
      leg_angulardiff_1jet->AddEntry(angulardiff_1jet_graph, "#beta = " + (TString)ss.str(), "l");
      leg_angulardiff_prop_1jet->AddEntry(angulardiff_prop_1jet_graph, "#beta = " + (TString)ss.str(), "l");
      leg_angulardiff_2jets->AddEntry(angulardiff_2jets_graph, "#beta = " + (TString)ss.str(), "l");
      leg_jet_distance_diff_1jet->AddEntry(jet_distance_diff_1jet_graph, "#beta = " + (TString)ss.str(), "l");
    }
  }
}

TLine *merging_point = new TLine();
merging_point->SetX1(500);
merging_point->SetVertical(true);
merging_point->SetY1(0);
merging_point->SetY2(1.0);
merging_point->SetLineColor(kRed);
merging_point->SetLineWidth(4);
merging_point->SetLineStyle(7);

TCanvas *higgs_1jet_can = new TCanvas("higgs_1jet_can", "higgs_1jet_can", 800, 600);
higgs_1jet_can->cd();
// gStyle->SetTitleOffset(0.8, "y");
higgs_efficiencies_1jet_biggraph->Draw("ALP");
leg_higgs_1jet->Draw();
radius_text->Draw();
merging_point->Draw();
higgs_efficiencies_1jet_biggraph->SetMinimum(0.);     
higgs_efficiencies_1jet_biggraph->SetMaximum(1.0);
higgs_efficiencies_1jet_biggraph->SetTitle("Higgs Eff. for 1-jet");
higgs_efficiencies_1jet_biggraph->GetXaxis()->SetTitle("p_{T,min} (GeV)");
higgs_efficiencies_1jet_biggraph->GetYaxis()->SetTitle("% mass 100-150");
// higgs_efficiencies_1jet_biggraph->GetYaxis()->SetTitleOffset(1.3);
higgs_1jet_can->Write();
higgs_1jet_can->Print("higgs2bbstudy_onepassXcone_test_plots/higgs_efficiencies_1jet.eps", "eps");

TCanvas *higgs_2jets_can = new TCanvas("higgs_2jets_can", "higgs_2jets_can", 800, 600);
higgs_2jets_can->cd();
higgs_efficiencies_2jets_biggraph->Draw("ALP");
leg_higgs_2jets->Draw();
radius_text->Draw();
merging_point->Draw();
higgs_efficiencies_2jets_biggraph->SetMinimum(0.);     
higgs_efficiencies_2jets_biggraph->SetMaximum(1.0);
higgs_efficiencies_2jets_biggraph->SetTitle("Higgs Eff. for 2-jet");
higgs_efficiencies_2jets_biggraph->GetXaxis()->SetTitle("p_{T,min} (GeV)");
higgs_efficiencies_2jets_biggraph->GetYaxis()->SetTitle("% mass 100-150");
// higgs_efficiencies_2jets_biggraph->GetYaxis()->SetTitleOffset(1.3);
higgs_2jets_can->Write();
higgs_2jets_can->Print("higgs2bbstudy_onepassXcone_test_plots/higgs_efficiencies_2jets.eps", "eps");

TCanvas *higgs_3jets_can = new TCanvas("higgs_3jets_can", "higgs_3jets_can", 800, 600);
higgs_3jets_can->cd();
higgs_efficiencies_3jets_biggraph->Draw("ALP");
leg_higgs_3jets->Draw();
radius_text->Draw();
merging_point->Draw();
higgs_efficiencies_3jets_biggraph->SetMinimum(0.);     
higgs_efficiencies_3jets_biggraph->SetMaximum(1.0);
higgs_efficiencies_3jets_biggraph->SetTitle("Higgs Eff. for 3-jettiness Min Mass");
higgs_efficiencies_3jets_biggraph->GetXaxis()->SetTitle("p_{T,min} (GeV)");
higgs_efficiencies_3jets_biggraph->GetYaxis()->SetTitle("% mass 100-150");
higgs_efficiencies_3jets_biggraph->GetYaxis()->SetTitleOffset(1.3);
higgs_3jets_can->Write();
higgs_3jets_can->Print("higgs2bbstudy_onepassXcone_test_plots/higgs_efficiencies_3jets.eps", "eps");

TCanvas *higgs_32ratio_can = new TCanvas("higgs_32ratio_can", "higgs_32ratio_can", 800, 600);
higgs_32ratio_can->cd();
higgs_efficiencies_32ratio_biggraph->Draw("ALP");
leg_higgs_32ratio->Draw();
radius_text->Draw();
higgs_efficiencies_32ratio_biggraph->SetMinimum(0.5);     
higgs_efficiencies_32ratio_biggraph->SetMaximum(1.5);
higgs_efficiencies_32ratio_biggraph->SetTitle("Ratio of 3-jettiness/2-jettiness Higgs higgs_efficiencies_32ratio");
higgs_efficiencies_32ratio_biggraph->GetXaxis()->SetTitle("p_{T,min} (GeV)");
higgs_efficiencies_32ratio_biggraph->GetYaxis()->SetTitle("3-jet eff/2-jet eff");
higgs_efficiencies_32ratio_biggraph->GetYaxis()->SetTitleOffset(1.3);
higgs_32ratio_can->Write();
higgs_32ratio_can->Print("higgs2bbstudy_onepassXcone_test_plots/higgs_efficiencies_32ratio.eps", "eps");

TCanvas *higgs_2jets_combined_can = new TCanvas("higgs_2jets_combined_can", "higgs_2jets_combined_can", 800, 600);
higgs_2jets_combined_can->cd();
higgs_efficiencies_2jets_combined_biggraph->Draw("ALP");
leg_higgs_2jets_combined->Draw();
radius_text->Draw();
merging_point->Draw();
higgs_efficiencies_2jets_combined_biggraph->SetMinimum(0.0);     
higgs_efficiencies_2jets_combined_biggraph->SetMaximum(1.0);
higgs_efficiencies_2jets_combined_biggraph->SetTitle("Higgs Eff. for 2-jet + ISR tag");
higgs_efficiencies_2jets_combined_biggraph->GetXaxis()->SetTitle("p_{T,min} (GeV)");
higgs_efficiencies_2jets_combined_biggraph->GetYaxis()->SetTitle("% mass 100-150");
// higgs_efficiencies_2jets_combined_biggraph->GetYaxis()->SetTitleOffset(1.3);
higgs_2jets_combined_can->Write();
higgs_2jets_combined_can->Print("higgs2bbstudy_onepassXcone_test_plots/higgs_efficiencies_2jets_combined.eps", "eps");


TCanvas *angulardiff_1jet_can = new TCanvas("angulardiff_1jet_can", "angulardiff_1jet_can", 800, 600);
angulardiff_1jet_can->cd();
angulardiff_1jet_biggraph->Draw("ALP");
leg_angulardiff_1jet->Draw();
radius_text->Draw();
angulardiff_1jet_biggraph->SetMinimum(0.);
angulardiff_1jet_biggraph->SetMaximum(0.25);
angulardiff_1jet_biggraph->SetTitle("Distance between 1-jettiness jet and closest 2-jettiness jet");
angulardiff_1jet_biggraph->GetXaxis()->SetTitle("p_{T,min} (GeV)");
angulardiff_1jet_biggraph->GetYaxis()->SetTitle("Angular Distance");
angulardiff_1jet_biggraph->GetYaxis()->SetTitleOffset(1.3);
angulardiff_1jet_can->Write();
angulardiff_1jet_can->Print("higgs2bbstudy_onepassXcone_test_plots/angulardiff_1jet.eps", "eps");


TCanvas *angulardiff_prop_1jet_can = new TCanvas("angulardiff_prop_1jet_can", "angulardiff_prop_1jet_can", 800, 600);
angulardiff_prop_1jet_can->cd();
angulardiff_prop_1jet_biggraph->Draw("ALP");
leg_angulardiff_prop_1jet->Draw();
radius_text->Draw();
angulardiff_prop_1jet_biggraph->SetMinimum(0.);
angulardiff_prop_1jet_biggraph->SetMaximum(1.0);
angulardiff_prop_1jet_biggraph->SetTitle("Proportion of 1-jettiness jets aligned with 2-jettiness jets");
angulardiff_prop_1jet_biggraph->GetXaxis()->SetTitle("p_{T,min} (GeV)");
angulardiff_prop_1jet_biggraph->GetYaxis()->SetTitle("Angular Distance");
angulardiff_prop_1jet_biggraph->GetYaxis()->SetTitleOffset(1.3);
angulardiff_prop_1jet_can->Write();
angulardiff_prop_1jet_can->Print("higgs2bbstudy_onepassXcone_test_plots/angulardiff_prop_1jet.eps", "eps");

TCanvas *angulardiff_2jets_can = new TCanvas("angulardiff_2jets_can", "angulardiff_2jets_can", 800, 600);
angulardiff_2jets_can->cd();
angulardiff_2jets_biggraph->Draw("ALP");
leg_angulardiff_2jets->Draw();
radius_text->Draw();
angulardiff_2jets_biggraph->SetMinimum(0.);
angulardiff_2jets_biggraph->SetMaximum(4.0);
angulardiff_2jets_biggraph->SetTitle("Angular Distance between two jets");
angulardiff_2jets_biggraph->GetXaxis()->SetTitle("p_{T,min} (GeV)");
angulardiff_2jets_biggraph->GetYaxis()->SetTitle("Angular Distance");
angulardiff_2jets_biggraph->GetYaxis()->SetTitleOffset(1.3);
angulardiff_2jets_can->Write();
angulardiff_2jets_can->Print("higgs2bbstudy_onepassXcone_test_plots/angulardiff_2jets.eps", "eps");

TCanvas* jet_distance_diff_1jet_can = new TCanvas("jet_distance_diff_1jet_can", "jet_distance_diff_1jet_can", 800, 600);
jet_distance_diff_1jet_can->cd();

for (int B = 0; B < n_betas; B++) {
  TH2* jetphi_violin = new TH2F("jetphi_violin", "Violin plot", n_perps, 0, 1000, 100, 0, 1.0);

  for (int i_perp = 0; i_perp < n_perps; i_perp++) {
    TH1 *jetphi_diff = (TH1*)jetphi_diff_allpt_allbeta.At(4*i_perp + B);
    for (int n = 0; n < jetphi_diff->GetNbinsX(); n++) {
      jetphi_violin->Fill(perplist[i_perp], jetphi_diff->GetBinContent(n));
    }
  }
  jetphi_violin->SetMarkerStyle(20);
  jetphi_violin->SetMarkerSize(0.5);
  jetphi_violin->SetFillStyle(1001);
  if (B == 2) {
    jetphi_violin->SetFillColor(kRed);
    jetphi_violin->SetMarkerColor(kRed);
    jetphi_violin->Draw("VIOLIN");
  }
  if (B == 3) {
    jetphi_violin->SetFillColor(kBlue);
    jetphi_violin->SetMarkerColor(kBlue);
    jetphi_violin->Draw("VIOLIN SAMES");
  }
}
// jet_distance_diff_1jet_biggraph->Draw("ALP");
// leg_jet_distance_diff_1jet->Draw();
// radius_text->Draw();
// jet_distance_diff_1jet_biggraph->SetMinimum(0.);
// jet_distance_diff_1jet_biggraph->SetMaximum(0.1);
// jet_distance_diff_1jet_biggraph->SetTitle("Avg Angular Distance between AKT and 1-jet");
// jet_distance_diff_1jet_biggraph->GetXaxis()->SetTitle("p_{T,min} (GeV)");
// jet_distance_diff_1jet_biggraph->GetYaxis()->SetTitle("Avg Angular Distance");
// jet_distance_diff_1jet_biggraph->GetYaxis()->SetTitleOffset(1.3);
jet_distance_diff_1jet_can->Write();
jet_distance_diff_1jet_can->Print("higgs2bbstudy_onepassXcone_test_plots/jet_distance_diff_1jet.eps", "eps");

delete outFile;

return 0;
}