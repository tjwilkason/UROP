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
#include "TString.h"

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
  
double calcDphi(PseudoJet jet1, PseudoJet jet2) {
  double dphi;
  dphi = jet1.phi() - jet2.phi();
  if (dphi >  TMath::Pi()) dphi -= 2*TMath::Pi();
  if (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
  return dphi; 
}

double calcDphi(double phi1, double phi2) {
  double dphi;
  dphi = phi1 - phi2;
  if (dphi >  TMath::Pi()) dphi -= 2*TMath::Pi();
  if (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
  return dphi; 
}
vector<PseudoJet> findMinAxes(vector<PseudoJet> input_particles, vector<PseudoJet> starting_axes, int njettiness, double beta, double Rcutoff) {

  vector<PseudoJet> min_manual_axes;
  // JetInformation min_manual_set;

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
        // min_manual_set.tau = manual_tau;
        // min_manual_set.axes = extras_manual->axes();
        // min_manual_set.jets = extras_manual->jets();
      }
    }
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  return min_manual_axes;
}

int main(int argc, char* argv[]) {

  double epsilon = 0.0001;
  int nEvents = 10000;
  gStyle->SetOptStat(0);
  // Create the ROOT application environment.
  TApplication theApp("hist", &argc, argv);

  // Generator. Process selection. LHC upgrade initialization. Create two pythia processes for gluons vs quarks
  Pythia pythia;
  pythia.readString("Beams:eCM = 14000.");    //LHC Upgrade COM Energy
  pythia.readString("NewGaugeBoson:ffbar2gmZZprime = on"); //Turn on Z'
  pythia.readString("Zprime:gmZmode = 3"); // No Z'/Z/gamma interference (only Z')
  pythia.readString("32:onMode = off"); // Turn off all Z decays
  pythia.readString("32:onIfAny = 1 2 3"); //Only turn on Z'->u, d, s decays
  pythia.readString("32:m0 = 1000");
  pythia.readString("PhaseSpace:pTHatMin = 100");
  pythia.init();

  // Create file on which histograms can be saved.
  TFile* outFile = new TFile("zprimestudy.root", "RECREATE");

    //create list of various values of beta
  vector<double> betalist;
  betalist.push_back(0.25);
  betalist.push_back(0.5);
  betalist.push_back(1.0);
  betalist.push_back(2.0);
  // betalist.push_back(3.0);
  int n_betas = betalist.size();

  TObjArray njet_jet_mass_hists(n_betas);
  TObjArray njet_jet_perp_hists(n_betas);
  TObjArray njet_jet_phi_hists(n_betas);
  TObjArray njet_jet_eta_hists(n_betas);
  
  TObjArray jetmass_diff_hists(n_betas);
  TObjArray jetperp_diff_hists(n_betas);
  TObjArray jetphi_diff_hists(n_betas);
  TObjArray tau1_diff_hists(n_betas);
  TObjArray jetmass_diff_wta_hists(n_betas);
  TObjArray jetperp_diff_wta_hists(n_betas);
  TObjArray jetphi_diff_wta_hists(n_betas);
  TObjArray tau1_diff_wta_hists(n_betas);

  TObjArray zprime_invmass_njet_hists(n_betas);
  TObjArray zprime_invperp_njet_hists(n_betas);
  TObjArray jet_ptspectrum_njet_hists(n_betas);
  TObjArray jet_phispectrum_njet_hists(n_betas);
  TObjArray twojet_bothjets_perp_diff_hists(n_betas);
  TObjArray twojet_bothjets_phi_diff_hists(n_betas);

  TObjArray twojet_firstjet_mass_njet_hists(n_betas);
  TObjArray twojet_firstjet_perp_njet_hists(n_betas);
  TObjArray twojet_firstjet_phi_njet_hists(n_betas);

  TObjArray twojet_secondjet_mass_njet_hists(n_betas);
  TObjArray twojet_secondjet_perp_njet_hists(n_betas);
  TObjArray twojet_secondjet_phi_njet_hists(n_betas);

  TObjArray zprime_invmass_diff_hists(n_betas);
  TObjArray jet_phispectrum_diff_hists(n_betas);
  TObjArray tau2_diff_hists(n_betas);

  TObjArray twojet_firstjet_mass_diff_hists(n_betas);
  TObjArray twojet_firstjet_perp_diff_hists(n_betas);
  TObjArray twojet_firstjet_phi_diff_hists(n_betas);

  TObjArray twojet_secondjet_mass_diff_hists(n_betas);
  TObjArray twojet_secondjet_perp_diff_hists(n_betas);
  TObjArray twojet_secondjet_phi_diff_hists(n_betas);

  TObjArray zprime_invmass_diff_wta_hists(n_betas);
  TObjArray jet_phispectrum_diff_wta_hists(n_betas);
  TObjArray tau2_diff_wta_hists(n_betas);  

  TObjArray twojet_firstjet_mass_diff_wta_hists(n_betas);
  TObjArray twojet_firstjet_perp_diff_wta_hists(n_betas);
  TObjArray twojet_firstjet_phi_diff_wta_hists(n_betas);

  TObjArray twojet_secondjet_mass_diff_wta_hists(n_betas);
  TObjArray twojet_secondjet_perp_diff_wta_hists(n_betas);
  TObjArray twojet_secondjet_phi_diff_wta_hists(n_betas);

  TObjArray threejet_invmass_njet_hists(n_betas);

  // TObjArray threejet_firstjet_mass_njet_hists(n_betas);
  // TObjArray threejet_firstjet_perp_njet_hists(n_betas);
  // TObjArray threejet_firstjet_phi_njet_hists(n_betas);

  // TObjArray threejet_secondjet_mass_njet_hists(n_betas);
  // TObjArray threejet_secondjet_perp_njet_hists(n_betas);
  // TObjArray threejet_secondjet_phi_njet_hists(n_betas);

  TObjArray threejet_thirdjet_mass_njet_hists(n_betas);
  TObjArray threejet_thirdjet_perp_njet_hists(n_betas);
  TObjArray threejet_thirdjet_phi_njet_hists(n_betas);

  TObjArray threejet_invmass_diff_hists(n_betas);
  // TObjArray threejet_firstjet_mass_diff_hists(n_betas);
  // TObjArray threejet_firstjet_perp_diff_hists(n_betas);
  // TObjArray threejet_firstjet_phi_diff_hists(n_betas);

  // TObjArray threejet_secondjet_mass_diff_hists(n_betas);
  // TObjArray threejet_secondjet_perp_diff_hists(n_betas);
  // TObjArray threejet_secondjet_phi_diff_hists(n_betas);

  TObjArray threejet_thirdjet_mass_diff_hists(n_betas);
  TObjArray threejet_thirdjet_perp_diff_hists(n_betas);
  TObjArray threejet_thirdjet_phi_diff_hists(n_betas);
  TObjArray tau3_diff_hists(n_betas);

  TObjArray threejet_invmass_diff_wta_hists(n_betas);
  // TObjArray threejet_firstjet_mass_diff_wta_hists(n_betas);
  // TObjArray threejet_firstjet_perp_diff_wta_hists(n_betas);
  // TObjArray threejet_firstjet_phi_diff_wta_hists(n_betas);

  // TObjArray threejet_secondjet_mass_diff_wta_hists(n_betas);
  // TObjArray threejet_secondjet_perp_diff_wta_hists(n_betas);
  // TObjArray threejet_secondjet_phi_diff_wta_hists(n_betas);

  TObjArray threejet_thirdjet_mass_diff_wta_hists(n_betas);
  TObjArray threejet_thirdjet_perp_diff_wta_hists(n_betas);
  TObjArray threejet_thirdjet_phi_diff_wta_hists(n_betas);
  TObjArray tau3_diff_wta_hists(n_betas);


  for (int i_betas = 0; i_betas < n_betas; i_betas++) {

    double beta = betalist[i_betas];

    ostringstream ss;
    ss << beta;
    TString title;
    
    TH1* njet_jet_mass = new TH1F("njet_jet_mass", "Mass of individual njet jet (#beta = " + (TString)ss.str() + ")", 50, 0, 500);
    TH1* njet_jet_perp = new TH1F("njet_jet_perp", "Perp of individual njet jet (#beta = " + (TString)ss.str() + ")", 50, 0, 1000);
    TH1* njet_jet_phi = new TH1F("njet_jet_phi", "Phi of individual njet jet (#beta = " + (TString)ss.str() + ")", 32, 0, 6.4);
    TH1* njet_jet_eta = new TH1F("njet_jet_eta", "Eta of individual njet jet (#beta = " + (TString)ss.str() + ")", 50, -5, 5);

    TH1* jetmass_diff = new TH1F("jetmass_diff", "difference in mass of akt and njet jets (#beta = " + (TString)ss.str() + ")", 100, -100, 200);
    TH1* jetperp_diff = new TH1F("jetperp_diff", "difference in perp of akt and njet jets (#beta = " + (TString)ss.str() + ")", 100, -100, 200);
    TH1* jetphi_diff = new TH1F("jetphi_diff", "difference in phi of akt and njet jets (#beta = " + (TString)ss.str() + ")", 32, 0, 3.2);
    TH1* tau1_diff = new TH1F("tau1_diff", "difference in #tau_{1} jets (akt - njet) (#beta = " + (TString)ss.str() + ")", 100, -200, 200);
    TH1* tau2_diff = new TH1F("tau2_diff", "difference in #tau_{2} jets (akt - njet) (#beta = " + (TString)ss.str() + ")", 100, -200, 200);
    TH1* jetmass_diff_wta = new TH1F("jetmass_diff_wta", "difference in mass of akt and njet jets (#beta = " + (TString)ss.str() + ")", 100, -100, 100);
    TH1* jetperp_diff_wta = new TH1F("jetperp_diff_wta", "difference in perp of akt and njet jets (#beta = " + (TString)ss.str() + ")", 100, -100, 100);
    TH1* jetphi_diff_wta = new TH1F("jetphi_diff_wta", "difference in phi of akt and njet jets (#beta = " + (TString)ss.str() + ")", 32, 0, 3.2);
    TH1* tau1_diff_wta = new TH1F("tau1_diff_wta", "difference in #tau_{1} jets (akt - njet) (#beta = " + (TString)ss.str() + ")", 100, -200, 200);
    TH1* tau2_diff_wta = new TH1F("tau2_diff_wta", "difference in #tau_{2} jets (akt - njet) (#beta = " + (TString)ss.str() + ")", 100, -200, 200);

    njet_jet_mass_hists.Add(njet_jet_mass);
    njet_jet_perp_hists.Add(njet_jet_perp);
    njet_jet_phi_hists.Add(njet_jet_phi);
    njet_jet_eta_hists.Add(njet_jet_eta);

    jetmass_diff_hists.Add(jetmass_diff);
    jetperp_diff_hists.Add(jetperp_diff);
    jetphi_diff_hists.Add(jetphi_diff);
    tau2_diff_hists.Add(tau2_diff);
    tau1_diff_hists.Add(tau1_diff);

    jetmass_diff_wta_hists.Add(jetmass_diff_wta);
    jetperp_diff_wta_hists.Add(jetperp_diff_wta);
    jetphi_diff_wta_hists.Add(jetphi_diff_wta);
    tau2_diff_wta_hists.Add(tau2_diff_wta);
    tau1_diff_wta_hists.Add(tau1_diff_wta);

    TH1* zprime_invmass_njet = new TH1F("zprime_invmass_njet", "Invariant Mass of Z' (#beta = " + (TString)ss.str() + ")", 50, 0, 1500);
    TH1* zprime_invperp_njet = new TH1F("zprime_invperp_njet", "pT of recombined dijet (#beta = " + (TString)ss.str() + ")", 50, 0, 1000);
    TH1* jet_ptspectrum_njet = new TH1F("jet_ptspectrum_njet", "pT spectrum of jets (#beta = " + (TString)ss.str() + ")", 50, 0, 1000);
    TH1* jet_phispectrum_njet = new TH1F("jet_phispectrum_njet", "#phi spectrum of jets (#beta = " + (TString)ss.str() + ")", 32, 0, 3.2);
    TH1* twojet_bothjets_perp_diff = new TH1F("twojet_bothjets_perp_diff", "pT spectrum of jets (#beta = " + (TString)ss.str() + ")", 100, -200, 200);
    TH1* twojet_bothjets_phi_diff = new TH1F("twojet_bothjets_phi_diff", "#phi spectrum of jets (#beta = " + (TString)ss.str() + ")", 32, 0, 3.2);

    TH1* twojet_firstjet_mass_njet = new TH1F("twojet_firstjet_mass_njet", "Mass of twojet_first Hardest Jet (#beta = " + (TString)ss.str() + ")", 50, 0, 500);
    TH1* twojet_firstjet_perp_njet = new TH1F("twojet_firstjet_perp_njet", "pT of twojet_first Hardest Jet (#beta = " + (TString)ss.str() + ")", 50, 0, 1000);
    TH1* twojet_firstjet_phi_njet = new TH1F("twojet_firstjet_phi_njet", "Azimuth of twojet_first Hardest Jet (#beta = " + (TString)ss.str() + ")", 32, 0, 6.4);

    TH1* twojet_secondjet_mass_njet = new TH1F("twojet_secondjet_mass_njet", "Mass of twojet_second Hardest Jet (#beta = " + (TString)ss.str() + ")", 50, 0, 500);
    TH1* twojet_secondjet_perp_njet = new TH1F("twojet_secondjet_perp_njet", "pT of twojet_second Hardest Jet (#beta = " + (TString)ss.str() + ")", 50, 0, 1000);
    TH1* twojet_secondjet_phi_njet = new TH1F("twojet_secondjet_phi_njet", "Azimuth of twojet_second Hardest Jet (#beta = " + (TString)ss.str() + ")", 32, 0, 6.4);

    TH1* threejet_invmass_njet = new TH1F("threejet_invmass_njet", "Invariant Mass of three reconstructed jets (#beta = " + (TString)ss.str() + ")", 50, 0, 1500);

    // TH1* threejet_firstjet_mass_njet = new TH1F("threejet_firstjet_mass_njet", "Mass of threejet_first Hardest Jet (#beta = " + (TString)ss.str() + ")", 50, 0, 500);
    // TH1* threejet_firstjet_perp_njet = new TH1F("threejet_firstjet_perp_njet", "pT of threejet_first Hardest Jet (#beta = " + (TString)ss.str() + ")", 50, 0, 1000);
    // TH1* threejet_firstjet_phi_njet = new TH1F("threejet_firstjet_phi_njet", "Azimuth of threejet_first Hardest Jet (#beta = " + (TString)ss.str() + ")", 32, 0, 6.4);

    // TH1* threejet_secondjet_mass_njet = new TH1F("threejet_secondjet_mass_njet", "Mass of threejet_second Hardest Jet (#beta = " + (TString)ss.str() + ")", 50, 0, 500);
    // TH1* threejet_secondjet_perp_njet = new TH1F("threejet_secondjet_perp_njet", "pT of threejet_second Hardest Jet (#beta = " + (TString)ss.str() + ")", 50, 0, 1000);
    // TH1* threejet_secondjet_phi_njet = new TH1F("threejet_secondjet_phi_njet", "Azimuth of threejet_second Hardest Jet (#beta = " + (TString)ss.str() + ")", 32, 0, 6.4);

    TH1* threejet_thirdjet_mass_njet = new TH1F("threejet_thirdjet_mass_njet", "Mass of threejet_third Hardest Jet (#beta = " + (TString)ss.str() + ")", 50, 0, 200);
    TH1* threejet_thirdjet_perp_njet = new TH1F("threejet_thirdjet_perp_njet", "pT of threejet_third Hardest Jet (#beta = " + (TString)ss.str() + ")", 50, 0, 1000);
    TH1* threejet_thirdjet_phi_njet = new TH1F("threejet_thirdjet_phi_njet", "Azimuth of threejet_third Hardest Jet (#beta = " + (TString)ss.str() + ")", 32, 0, 6.4);

    TH1* threejet_invmass_diff = new TH1F("threejet_invmass_diff", "Invariant Mass of three reconstructed jets (#beta = " + (TString)ss.str() + ")", 50, 0, 1500);
    TH1* threejet_thirdjet_mass_diff = new TH1F("threejet_thirdjet_mass_diff", "Mass of threejet_third Hardest Jet (#beta = " + (TString)ss.str() + ")", 50, -50, 50);
    TH1* threejet_thirdjet_perp_diff = new TH1F("threejet_thirdjet_perp_diff", "pT of threejet_third Hardest Jet (#beta = " + (TString)ss.str() + ")", 100, -200, 200);
    TH1* threejet_thirdjet_phi_diff = new TH1F("threejet_thirdjet_phi_diff", "Azimuth of threejet_third Hardest Jet (#beta = " + (TString)ss.str() + ")", 32, 0, 3.2);

    TH1* threejet_invmass_diff_wta = new TH1F("threejet_invmass_diff_wta", "Invariant Mass of three reconstructed jets (#beta = " + (TString)ss.str() + ")", 50, 0, 1500);
    TH1* threejet_thirdjet_mass_diff_wta = new TH1F("threejet_thirdjet_mass_diff_wta", "Mass of threejet_third Hardest Jet (#beta = " + (TString)ss.str() + ")", 50, -50, 50);
    TH1* threejet_thirdjet_perp_diff_wta = new TH1F("threejet_thirdjet_perp_diff_wta", "pT of threejet_third Hardest Jet (#beta = " + (TString)ss.str() + ")", 100, -200, 200);
    TH1* threejet_thirdjet_phi_diff_wta = new TH1F("threejet_thirdjet_phi_diff_wta", "Azimuth of threejet_third Hardest Jet (#beta = " + (TString)ss.str() + ")", 32, 0, 3.2);

    TH1* tau3_diff = new TH1F("tau3_diff", "Difference in tau3 (#beta = " + (TString)ss.str() + ")", 100, -200, 200);
    TH1* tau3_diff_wta = new TH1F("tau3_diff_wta", "Difference in tau3 (#beta = " + (TString)ss.str() + ")", 100, -200, 200);

    TH1* zprime_invmass_diff = new TH1F("zprime_invmass_diff", "Invariant Mass of Z' (#beta = " + (TString)ss.str() + ")", 100, -500, 500);
    TH1* jet_phispectrum_diff = new TH1F("jet_phispectrum_diff", "#phi spectrum of jets (#beta = " + (TString)ss.str() + ")", 32, 0, 3.2);

    TH1* zprime_invmass_diff_wta = new TH1F("zprime_invmass_diff_wta", "Invariant Mass of Z' (#beta = " + (TString)ss.str() + ")", 100, -500, 500);
    TH1* jet_phispectrum_diff_wta = new TH1F("jet_phispectrum_diff_wta", "#phi spectrum of jets (#beta = " + (TString)ss.str() + ")", 32, 0, 3.2);

    zprime_invmass_njet_hists.Add(zprime_invmass_njet);
    zprime_invperp_njet_hists.Add(zprime_invperp_njet);
    jet_ptspectrum_njet_hists.Add(jet_ptspectrum_njet);
    jet_phispectrum_njet_hists.Add(jet_phispectrum_njet);
    twojet_bothjets_perp_diff_hists.Add(twojet_bothjets_perp_diff);
    twojet_bothjets_phi_diff_hists.Add(twojet_bothjets_phi_diff);

    zprime_invmass_diff_hists.Add(zprime_invmass_diff);
    jet_phispectrum_diff_hists.Add(jet_phispectrum_diff);

    zprime_invmass_diff_wta_hists.Add(zprime_invmass_diff_wta);
    jet_phispectrum_diff_wta_hists.Add(jet_phispectrum_diff_wta);

    twojet_firstjet_mass_njet_hists.Add(twojet_firstjet_mass_njet);
    twojet_firstjet_perp_njet_hists.Add(twojet_firstjet_perp_njet);
    twojet_firstjet_phi_njet_hists.Add(twojet_firstjet_phi_njet);

    twojet_secondjet_mass_njet_hists.Add(twojet_secondjet_mass_njet);
    twojet_secondjet_perp_njet_hists.Add(twojet_secondjet_perp_njet);
    twojet_secondjet_phi_njet_hists.Add(twojet_secondjet_phi_njet);

    threejet_invmass_njet_hists.Add(threejet_invmass_njet);

    // threejet_firstjet_mass_njet_hists.Add(threejet_firstjet_mass_njet);
    // threejet_firstjet_perp_njet_hists.Add(threejet_firstjet_perp_njet);
    // threejet_firstjet_phi_njet_hists.Add(threejet_firstjet_phi_njet);

    // threejet_secondjet_mass_njet_hists.Add(threejet_secondjet_mass_njet);
    // threejet_secondjet_perp_njet_hists.Add(threejet_secondjet_perp_njet);
    // threejet_secondjet_phi_njet_hists.Add(threejet_secondjet_phi_njet);

    threejet_thirdjet_mass_njet_hists.Add(threejet_thirdjet_mass_njet);
    threejet_thirdjet_perp_njet_hists.Add(threejet_thirdjet_perp_njet);
    threejet_thirdjet_phi_njet_hists.Add(threejet_thirdjet_phi_njet);

    threejet_invmass_diff_hists.Add(threejet_invmass_diff);
    threejet_thirdjet_mass_diff_hists.Add(threejet_thirdjet_mass_diff);
    threejet_thirdjet_perp_diff_hists.Add(threejet_thirdjet_perp_diff);
    threejet_thirdjet_phi_diff_hists.Add(threejet_thirdjet_phi_diff);

    threejet_invmass_diff_wta_hists.Add(threejet_invmass_diff_wta);
    threejet_thirdjet_mass_diff_wta_hists.Add(threejet_thirdjet_mass_diff_wta);
    threejet_thirdjet_perp_diff_wta_hists.Add(threejet_thirdjet_perp_diff_wta);
    threejet_thirdjet_phi_diff_wta_hists.Add(threejet_thirdjet_phi_diff_wta);

    tau3_diff_hists.Add(tau3_diff);
    tau3_diff_wta_hists.Add(tau3_diff_wta);

  }

  TH1* akt_jet_mass = new TH1F("akt_jet_mass", "Mass of individual akt jet", 50, 0, 500);
  TH1* akt_jet_perp = new TH1F("akt_jet_perp", "Perp of individual akt jet", 50, 0, 1000);
  TH1* akt_jet_phi = new TH1F("akt_jet_phi", "Phi of individual akt jet", 32, 0, 6.4);
  TH1* akt_jet_eta = new TH1F("akt_jet_eta", "Eta of individual akt jet", 50, -5, 5);

  TH1* akt_wta_jet_mass = new TH1F("akt_wta_jet_mass", "Mass of individual akt_wta jet", 50, 0, 500);
  TH1* akt_wta_jet_perp = new TH1F("akt_wta_jet_perp", "Perp of individual akt_wta jet", 50, 0, 1000);
  TH1* akt_wta_jet_phi = new TH1F("akt_wta_jet_phi", "Phi of individual akt_wta jet", 32, 0, 6.4);
  TH1* akt_wta_jet_eta = new TH1F("akt_wta_jet_eta", "Eta of individual akt_wta jet", 50, -5, 5);

  TH1* zprime_invmass_akt = new TH1F("zprime_invmass_akt", "Invariant Mass of Z'", 50, 0, 1500);
  TH1* zprime_invperp_akt = new TH1F("zprime_invperp_akt", "pT of recombined dijet", 50, 0, 1000);
  TH1* jet_ptspectrum_akt = new TH1F("jet_ptspectrum_akt", "pT spectrum of jets", 50, 0, 1000);
  TH1* jet_phispectrum_akt = new TH1F("jet_phispectrum_akt", "#phi spectrum of jets", 32, 0, 3.2);

  TH1* zprime_invmass_akt_wta = new TH1F("zprime_invmass_akt_wta", "Invariant Mass of Z'", 50, 0, 1500);
  TH1* zprime_invperp_akt_wta = new TH1F("zprime_invperp_akt_wta", "pT of recombined dijet", 50, 0, 1000);
  TH1* jet_ptspectrum_akt_wta = new TH1F("jet_ptspectrum_akt_wta", "pT spectrum of jets", 50, 0, 1000);
  TH1* jet_phispectrum_akt_wta = new TH1F("jet_phispectrum_akt_wta", "#phi spectrum of jets", 32, 0, 3.2);

  TH1* secondjet_mass_akt = new TH1F("secondjet_mass_akt", "Mass of second akt jet", 50, 0, 500);
  TH1* secondjet_perp_akt = new TH1F("secondjet_perp_akt", "pT of second akt jet", 50, 0, 1000);
  TH1* secondjet_phi_akt = new TH1F("secondjet_phi_akt", "Azimuth of second akt jet", 32, 0, 6.4);

  TH1* secondjet_mass_akt_wta = new TH1F("secondjet_mass_akt_wta", "Mass of second akt_wta jet", 50, 0, 500);
  TH1* secondjet_perp_akt_wta = new TH1F("secondjet_perp_akt_wta", "pT of second akt_wta jet", 50, 0, 1000);
  TH1* secondjet_phi_akt_wta = new TH1F("secondjet_phi_akt_wta", "Azimuth of second akt_wta jet", 32, 0, 6.4);

  TH1* threejet_invmass_akt = new TH1F("threejet_invmass_akt", "Invariant mass of three akt jets", 50, 0, 1500);
  TH1* threejet_thirdjet_mass_akt = new TH1F("thirdjet_mass_akt", "Mass of third akt jet", 50, 0, 500);
  TH1* threejet_thirdjet_perp_akt = new TH1F("thirdjet_perp_akt", "pT of third akt jet", 50, 0, 1000);
  TH1* threejet_thirdjet_phi_akt = new TH1F("thirdjet_phi_akt", "Azimuth of third akt jet", 32, 0, 6.4);

  TH1* threejet_invmass_akt_wta = new TH1F("threejet_invmass_akt_wta", "Invariant mass of three akt (WTA) jets", 50, 0, 1500);
  TH1* threejet_thirdjet_mass_akt_wta = new TH1F("thirdjet_mass_akt_wta", "Mass of third akt_wta jet", 50, 0, 500);
  TH1* threejet_thirdjet_perp_akt_wta = new TH1F("thirdjet_perp_akt_wta", "pT of third akt_wta jet", 50, 0, 1000);
  TH1* threejet_thirdjet_phi_akt_wta = new TH1F("thirdjet_phi_akt_wta", "Azimuth of third akt_wta jet", 32, 0, 6.4);

  // Fastjet input (different inputs for quark vs gluon)
  vector <PseudoJet> fjInputs;

  double Rparam = 0.6;
  Strategy strategy = Best;
  const JetDefinition::Recombiner *recombScheme_akt_wta = new WinnerTakeAllRecombiner();
  RecombinationScheme recombScheme_akt = E_scheme;
  JetDefinition *jetDef_wta = new JetDefinition(antikt_algorithm, Rparam, recombScheme_akt_wta, strategy);
  JetDefinition *jetDef = new JetDefinition(antikt_algorithm, Rparam, recombScheme_akt, strategy);

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

    if (!pythia.next()) continue;

    // Reset Fastjet input
    fjInputs.resize(0);

    //store particle as input for Fastjet
    for (int i = 0; i < pythia.event.size(); ++i) {
      // Final State only
      if (!pythia.event[i].isFinal()) continue;

      //Exclude neutrinos and leptons, which should only appear as final products of bosons, which should not be in the jets
      if (pythia.event[i].idAbs() == 12 || pythia.event[i].idAbs() == 14 || pythia.event[i].idAbs() == 16) continue;
      if (pythia.event[i].idAbs() == 11 || pythia.event[i].idAbs() == 13 || pythia.event[i].idAbs() == 15) continue;
      PseudoJet fj_particle = pythia.event[i];
      fjInputs.push_back(fj_particle);
    }

    TH2F *event_display = new TH2F("event_display", "Event Display", 60, -5, 5, 60, 0, 6.4);
    TH2F *axes_display = new TH2F("axes_display", "Axes Plot", 300, -5, 5, 300, 0, 6.4);
    TH2F *axes_njet_display = new TH2F("axes_njet_display", "Axes Plot (Njettiness)", 300, -5, 5, 300, 0, 6.4);
    TH2F *akt_jets_display = new TH2F("akt_jets_display", "Anti-KT Jets", 300, -5, 5, 300, 0, 6.4);
    TH2F *njet_jets_display = new TH2F("njet_jets_display", "N-jettiness Jets", 300, -5, 5, 300, 0, 6.4);

    for (int i_part = 0; i_part < fjInputs.size(); i_part++) {
      event_display->Fill(fjInputs[i_part].eta(), fjInputs[i_part].phi(), fjInputs[i_part].perp());
    }

    // Double_t ghost_perp = 0.001;
    // Double_t n_ghosts = 50;
    // for (int i_ghosts = 0; i_ghosts < n_ghosts; i_ghosts++) {
    //   for (int j_ghosts = 0; j_ghosts < n_ghosts; j_ghosts++) {
    //     Double_t ghost_eta = -5 + (double)(10/n_ghosts)*i_ghosts;
    //     Double_t ghost_phi = (double)(2*TMath::Pi()/n_ghosts)*j_ghosts;
    //     PseudoJet ghost(ghost_perp*TMath::Cos(ghost_phi),ghost_perp*TMath::Sin(ghost_phi),
    //       ghost_perp*TMath::SinH(ghost_eta),ghost_perp*TMath::CosH(ghost_eta));
    //     fjInputs.push_back(ghost);
    //   }
    // } 

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

    //FIRST WE COMPARE 1-JETTINESS TO THE HARDEST ANTI-KT JET AND QUANTIFY SIMILARITIES/DIFFERENCES

    //Create selectors for hardest jet in the event
    Selector jet_selector = SelectorNHardest(1);
    vector<PseudoJet> hardest_jet = jet_selector(centralJets);

    if (hardest_jet.size() == 1) {

      akt_jet_mass->Fill(hardest_jet[0].m());
      akt_jet_perp->Fill(hardest_jet[0].perp());
      akt_jet_phi->Fill(hardest_jet[0].phi());
      akt_jet_eta->Fill(hardest_jet[0].eta());
    }

    Selector selector = SelectorNHardest(2);
    vector<PseudoJet> hardest_twojets = selector(centralJets);
    PseudoJet big_jet(0,0,0,0);

    for (int i_jet = 0; i_jet < hardest_twojets.size(); i_jet++) {
      big_jet = join(big_jet, hardest_twojets[i_jet]);
      jet_ptspectrum_akt->Fill(hardest_twojets[i_jet].perp());
    }
    if (hardest_twojets.size() == 2) {
      jet_phispectrum_akt->Fill(abs(calcDphi(hardest_twojets[0], hardest_twojets[1])));
      zprime_invmass_akt->Fill(big_jet.m());
      zprime_invperp_akt->Fill(big_jet.perp());
      // PseudoJet second_jet = (hardest_twojets[0].perp() > hardest_twojets[1].perp()) ? hardest_twojets[1] : hardest_twojets[1];
      secondjet_mass_akt->Fill(hardest_twojets[1].m());
      secondjet_perp_akt->Fill(hardest_twojets[1].perp());
      secondjet_phi_akt->Fill(hardest_twojets[1].phi());
    }

    Selector threejet_selector = SelectorNHardest(3);
    vector<PseudoJet> hardest_threejets = threejet_selector(centralJets);

    PseudoJet third_jet = hardest_threejets[hardest_threejets.size() - 1];
    PseudoJet threejet_big_jet_akt(0,0,0,0);

    if ((third_jet.delta_R(hardest_twojets[0]) < 2*Rparam || third_jet.delta_R(hardest_twojets[1]) < 2*Rparam)) {
      threejet_big_jet_akt = join(big_jet, third_jet);
    }
    else threejet_big_jet_akt = big_jet;

    // for (int i_jet = 0; i_jet < hardest_threejets.size(); i_jet++) {
    //   threejet_big_jet_akt = join(threejet_big_jet_akt, hardest_threejets[i_jet]);
    // }
    threejet_invmass_akt->Fill(threejet_big_jet_akt.m());
    threejet_thirdjet_mass_akt->Fill(third_jet.m());
    threejet_thirdjet_perp_akt->Fill(third_jet.perp());
    threejet_thirdjet_phi_akt->Fill(third_jet.phi());

    //Run Fastjet algorithm
    vector <PseudoJet> inclusiveJets_wta, sortedJets_wta, centralJets_wta;
    ClusterSequence clustSeq_wta(fjInputs, *jetDef_wta);

    //Extract inclusive jets sorted by pT
    inclusiveJets_wta = clustSeq_wta.inclusive_jets(0.0);

    // Sort central jets by pT
    sortedJets_wta = sorted_by_pt(inclusiveJets_wta);

    //Cut on jets within eta range
    centralJets_wta = eta_selector(sortedJets_wta);

    //FIRST WE COMPARE 1-JETTINESS TO THE HARDEST ANTI-KT JET

    //Create selectors for hardest jet in the event
    vector<PseudoJet> hardest_jet_wta = jet_selector(centralJets_wta);

    if (hardest_jet_wta.size() == 1) {
      PseudoJet new_jet = join(hardest_jet_wta[0].constituents());
      akt_wta_jet_mass->Fill(new_jet.m());
      akt_wta_jet_perp->Fill(new_jet.perp());
      akt_wta_jet_phi->Fill(new_jet.phi());
      akt_wta_jet_eta->Fill(new_jet.eta());
    }

    vector<PseudoJet> hardest_twojets_wta = selector(centralJets_wta);
    PseudoJet big_jet_wta(0,0,0,0);

    for (int i_jet = 0; i_jet < hardest_twojets_wta.size(); i_jet++) {
      big_jet_wta = join(big_jet_wta, hardest_twojets_wta[i_jet]);
      jet_ptspectrum_akt_wta->Fill(hardest_twojets_wta[i_jet].perp());
    }
    if (hardest_twojets_wta.size() == 2) {
      jet_phispectrum_akt_wta->Fill(abs(calcDphi(hardest_twojets_wta[0], hardest_twojets_wta[1])));
      zprime_invmass_akt_wta->Fill(big_jet_wta.m());
      zprime_invperp_akt_wta->Fill(big_jet_wta.perp());
      // PseudoJet second_jet_wta = (hardest_twojets_wta[0].perp() > hardest_twojets_wta[1].perp()) ? hardest_twojets_wta[1] : hardest_twojets_wta[1];
      PseudoJet second_jet_wta = hardest_twojets_wta[1];
      second_jet_wta = join(second_jet_wta.constituents());
      secondjet_mass_akt_wta->Fill(second_jet_wta.m());
      secondjet_perp_akt_wta->Fill(second_jet_wta.perp());
      secondjet_phi_akt_wta->Fill(second_jet_wta.phi());
    }

    vector<PseudoJet> hardest_threejets_wta = threejet_selector(centralJets_wta);

    PseudoJet threejet_big_jet_akt_wta(0,0,0,0);
    PseudoJet third_jet_wta = hardest_threejets_wta[hardest_threejets_wta.size() - 1];
    third_jet_wta = join(third_jet_wta.constituents());

    if ((third_jet_wta.delta_R(hardest_twojets_wta[0]) < 2*Rparam || third_jet_wta.delta_R(hardest_twojets_wta[1]) < 2*Rparam)) {
      threejet_big_jet_akt_wta = join(big_jet_wta, third_jet_wta);
    }
    else threejet_big_jet_akt_wta = big_jet_wta;

    // for (int i_jet = 0; i_jet < hardest_threejets.size(); i_jet++) {
    //   threejet_big_jet_akt_wta = join(threejet_big_jet_akt_wta, hardest_threejets[i_jet]);
    // }
    threejet_invmass_akt_wta->Fill(threejet_big_jet_akt_wta.m());    
    threejet_thirdjet_mass_akt_wta->Fill(third_jet_wta.m());
    threejet_thirdjet_perp_akt_wta->Fill(third_jet_wta.perp());
    threejet_thirdjet_phi_akt_wta->Fill(third_jet_wta.phi());

    for (unsigned int B = 0; B < betalist.size(); B++) {

      double beta = betalist[B];
      double power = (double)1/beta;
      // double delta;

      const JetDefinition::Recombiner *recombScheme;
      if (beta > 1) recombScheme = new GeneralERecombiner((double)1/(beta - 1));
      else recombScheme = new WinnerTakeAllRecombiner();

      UnnormalizedCutoffMeasure measure_function = UnnormalizedCutoffMeasure(beta, Rparam);

      AxesStruct *axes_finder;
      if (beta < 1 || beta > 3) {
        axes_finder = new AxesStruct(Manual_Axes());
      }
      else {
        axes_finder = new AxesStruct(OnePass_Manual_Axes());
      }

      // OnePass_Manual_Axes axes_finder = OnePass_Manual_Axes();

      Strategy strategy = Best;
      // const JetDefinition::Recombiner *recombScheme = new GeneralERecombiner(delta);
      // const JetDefinition::Recombiner *recombScheme = new WinnerTakeAllRecombiner();
      JetDefinition *jetDef = new JetDefinition(genkt_algorithm, Rparam, power, recombScheme, strategy);
      ClusterSequence clustSeq(fjInputs, *jetDef);
      vector<PseudoJet> exclusive_1jet_start = clustSeq.exclusive_jets(2);
      vector<PseudoJet> exclusive_1jet = findMinAxes(fjInputs, exclusive_1jet_start, 1, beta, Rparam);

      // OnePass_WTA_KT_Axes axes_finder = OnePass_WTA_KT_Axes();

      NjettinessPlugin njet_plugin_1jet(1, axes_finder->def(), measure_function);
      JetDefinition njet_def_1jet(&njet_plugin_1jet);
      njet_plugin_1jet.setAxes(exclusive_1jet);
      ClusterSequence njet_cluster_1jet(fjInputs, njet_def_1jet);
      const NjettinessExtras *extras_1jet = njettiness_extras(njet_cluster_1jet);
      vector<PseudoJet> njet_1jet = extras_1jet->jets();

      NjettinessPlugin njet_plugin_manual_1jet(1, OnePass_Manual_Axes(), measure_function);
      JetDefinition njet_def_manual_1jet(&njet_plugin_manual_1jet);

      // if (B == 0) {
      //   for (int a = 0; a < njet_1jet.size(); a++) {
      //     axes_njet_display->Fill(njet_1jet[a].eta(), njet_1jet[a].phi());
      //     vector<PseudoJet> constituents = njet_1jet[a].constituents();
      //     for (int i_const = 0; i_const < constituents.size(); i_const++) {
      //       njet_jets_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
      //     }          
      //   }
      // }

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

          jetmass_diff->Fill(hardest_jet[0].m() - njet_1jet[0].m());
          jetperp_diff->Fill(hardest_jet[0].perp() - njet_1jet[0].perp());
          jetphi_diff->Fill(abs(calcDphi(hardest_jet[0], njet_1jet[0])));
        }

        if (hardest_jet_wta.size() == 1) {
          njet_plugin_manual_1jet.setAxes(hardest_jet_wta);
          ClusterSequence njet_cluster_manual_1jet(fjInputs, njet_def_manual_1jet);
          const NjettinessExtras *extras_manual_1jet = njettiness_extras(njet_cluster_manual_1jet);

          TH1* tau1_diff_wta = (TH1*)tau1_diff_wta_hists.At(B);
          tau1_diff_wta->Fill(extras_manual_1jet->totalTau() - extras_1jet->totalTau());

          TH1* jetmass_diff_wta = (TH1*)jetmass_diff_wta_hists.At(B);
          TH1* jetperp_diff_wta = (TH1*)jetperp_diff_wta_hists.At(B);
          TH1* jetphi_diff_wta = (TH1*)jetphi_diff_wta_hists.At(B);

          PseudoJet new_jet = join(hardest_jet_wta[0].constituents());

          jetmass_diff_wta->Fill(new_jet.m() - njet_1jet[0].m());
          jetperp_diff_wta->Fill(new_jet.perp() - njet_1jet[0].perp());
          jetphi_diff_wta->Fill(abs(calcDphi(new_jet, njet_1jet[0])));
        }
      }

      // NEXT WE LOOK AT 2-JETTINESS
      vector<PseudoJet> exclusive_2jets_start = clustSeq.exclusive_jets(3);
      vector<PseudoJet> exclusive_2jets = findMinAxes(fjInputs, exclusive_2jets_start, 2, beta, Rparam);

      // vector<PseudoJet> exclusive_2jets = clustSeq.exclusive_jets(2);
      //Use Njettiness algorithm for comparison
      NjettinessPlugin njet_plugin(2, axes_finder->def(), measure_function);
      JetDefinition njet_def(&njet_plugin);
      njet_plugin.setAxes(exclusive_2jets);
      ClusterSequence njet_cluster(fjInputs, njet_def);
      const NjettinessExtras *extras = njettiness_extras(njet_cluster);
      vector<PseudoJet> njet_jets = extras->jets();

      PseudoJet big_jet_njet(0,0,0,0);
      for (int i_jet = 0; i_jet < njet_jets.size(); i_jet++) {
        big_jet_njet = join(big_jet_njet, njet_jets[i_jet]);
        TH1* jet_ptspectrum_njet = (TH1*)jet_ptspectrum_njet_hists.At(B);
        jet_ptspectrum_njet->Fill(njet_jets[i_jet].perp());
      }
      if (njet_jets.size() == 2) {
        TH1* jet_phispectrum_njet = (TH1*)jet_phispectrum_njet_hists.At(B);
        TH1* twojet_firstjet_mass_njet = (TH1*)twojet_firstjet_mass_njet_hists.At(B);
        TH1* twojet_firstjet_perp_njet = (TH1*)twojet_firstjet_perp_njet_hists.At(B);
        TH1* twojet_firstjet_phi_njet = (TH1*)twojet_firstjet_phi_njet_hists.At(B);
        TH1* twojet_secondjet_mass_njet = (TH1*)twojet_secondjet_mass_njet_hists.At(B);
        TH1* twojet_secondjet_perp_njet = (TH1*)twojet_secondjet_perp_njet_hists.At(B);
        TH1* twojet_secondjet_phi_njet = (TH1*)twojet_secondjet_phi_njet_hists.At(B);

        jet_phispectrum_njet->Fill(abs(calcDphi(njet_jets[0], njet_jets[1])));
        vector<PseudoJet> sorted_njet_jets = sorted_by_pt(njet_jets);

        // PseudoJet hardest_jet = (njet_jets[0].delta_R(njet_1jet[0]) > njet_jets[1].delta_R(njet_1jet[0])) ? njet_jets[1] : njet_jets[0];
        // PseudoJet second_jet = (njet_jets[0].delta_R(njet_1jet[0]) > njet_jets[1].delta_R(njet_1jet[0])) ? njet_jets[0] : njet_jets[1];
        twojet_firstjet_mass_njet->Fill(sorted_njet_jets[0].m());
        twojet_firstjet_perp_njet->Fill(sorted_njet_jets[0].perp());
        twojet_firstjet_phi_njet->Fill(sorted_njet_jets[0].phi());
        twojet_secondjet_mass_njet->Fill(sorted_njet_jets[1].m());
        twojet_secondjet_perp_njet->Fill(sorted_njet_jets[1].perp());
        twojet_secondjet_phi_njet->Fill(sorted_njet_jets[1].phi());

      }

      TH1* zprime_invmass_njet = (TH1*)zprime_invmass_njet_hists.At(B);
      TH1* zprime_invperp_njet = (TH1*)zprime_invperp_njet_hists.At(B);

      zprime_invmass_njet->Fill(big_jet_njet.m());
      zprime_invperp_njet->Fill(big_jet_njet.perp());

      NjettinessPlugin njet_plugin_manual(2, OnePass_Manual_Axes(), measure_function);
      JetDefinition njet_def_manual(&njet_plugin_manual);
      if (hardest_twojets.size() == 2) {
        njet_plugin_manual.setAxes(hardest_twojets);
        ClusterSequence njet_cluster_manual(fjInputs, njet_def_manual);
        const NjettinessExtras *extras_manual = njettiness_extras(njet_cluster_manual);
        vector<PseudoJet> njet_jets_manual = extras_manual->jets();


        TH1* tau2_diff = (TH1*)tau2_diff_hists.At(B);
        tau2_diff->Fill(extras_manual->totalTau() - extras->totalTau());

        TH1* jet_phispectrum_diff = (TH1*)jet_phispectrum_diff_hists.At(B);
        TH1* twojet_bothjets_perp_diff = (TH1*)twojet_bothjets_perp_diff_hists.At(B);
        TH1* twojet_bothjets_phi_diff = (TH1*)twojet_bothjets_phi_diff_hists.At(B);

        PseudoJet closerJet0 = (njet_jets[0].delta_R(hardest_twojets[0]) < njet_jets[1].delta_R(hardest_twojets[0])) ? njet_jets[0] : njet_jets[1];
        PseudoJet closerJet1 = (njet_jets[0].delta_R(hardest_twojets[0]) < njet_jets[1].delta_R(hardest_twojets[0])) ? njet_jets[1] : njet_jets[0];

        twojet_bothjets_perp_diff->Fill(hardest_twojets[0].perp() - closerJet0.perp());
        twojet_bothjets_perp_diff->Fill(hardest_twojets[1].perp() - closerJet1.perp());
        twojet_bothjets_phi_diff->Fill(calcDphi(hardest_twojets[0], closerJet0));
        twojet_bothjets_phi_diff->Fill(calcDphi(hardest_twojets[1], closerJet1));

        double njet_phi_diff = abs(calcDphi(njet_jets[0], njet_jets[1]));
        double akt_phi_diff = abs(calcDphi(hardest_twojets[0], hardest_twojets[1]));
        double akt_wta_phi_diff = abs(calcDphi(hardest_twojets_wta[0], hardest_twojets_wta[1]));
        jet_phispectrum_diff->Fill(abs(calcDphi(njet_phi_diff, akt_phi_diff)));

      }

      TH1* zprime_invmass_diff = (TH1*)zprime_invmass_diff_hists.At(B);
      zprime_invmass_diff->Fill(big_jet.m() - big_jet_njet.m());

      // THEN WE LOOK AT 3-JETTINESS

      vector<PseudoJet> exclusive_3jets = clustSeq.exclusive_jets(3);
      NjettinessPlugin njet_plugin_threejets(3, axes_finder->def(), measure_function);
      JetDefinition njet_def_threejets(&njet_plugin_threejets);
      njet_plugin_threejets.setAxes(exclusive_3jets);
      ClusterSequence njet_cluster_threejets(fjInputs, njet_def_threejets);
      const NjettinessExtras *extras_threejets = njettiness_extras(njet_cluster_threejets);
      vector<PseudoJet> njet_threejets = extras_threejets->jets();

      PseudoJet threejet_big_jet_njet(0,0,0,0);
      vector<PseudoJet> sorted_njet_threejets = sorted_by_pt(njet_threejets);
      PseudoJet third_jet_njet = sorted_njet_threejets[2];

      if ((third_jet_njet.delta_R(njet_jets[0]) < 2*Rparam || third_jet_njet.delta_R(njet_jets[1]) < 2*Rparam)) {
        threejet_big_jet_njet = join(big_jet_njet, third_jet_njet);
      }
      else threejet_big_jet_njet = big_jet_njet;

      // for (int i_jet = 0; i_jet < njet_threejets.size(); i_jet++) {
      //   threejet_big_jet_njet = join(threejet_big_jet_njet, njet_threejets[i_jet]);
      // }
      TH1* threejet_invmass_njet = (TH1*)threejet_invmass_njet_hists.At(B);
      threejet_invmass_njet->Fill(threejet_big_jet_njet.m());

      // TH1* threejet_firstjet_mass_njet = (TH1*)threejet_firstjet_mass_njet_hists.At(B);
      // TH1* threejet_firstjet_perp_njet = (TH1*)threejet_firstjet_perp_njet_hists.At(B);
      // TH1* threejet_firstjet_phi_njet = (TH1*)threejet_firstjet_phi_njet_hists.At(B);
      // TH1* threejet_secondjet_mass_njet = (TH1*)threejet_secondjet_mass_njet_hists.At(B);
      // TH1* threejet_secondjet_perp_njet = (TH1*)threejet_secondjet_perp_njet_hists.At(B);
      // TH1* threejet_secondjet_phi_njet = (TH1*)threejet_secondjet_phi_njet_hists.At(B);
      TH1* threejet_thirdjet_mass_njet = (TH1*)threejet_thirdjet_mass_njet_hists.At(B);
      TH1* threejet_thirdjet_perp_njet = (TH1*)threejet_thirdjet_perp_njet_hists.At(B);
      TH1* threejet_thirdjet_phi_njet = (TH1*)threejet_thirdjet_phi_njet_hists.At(B);
      TH1* threejet_thirdjet_mass_diff = (TH1*)threejet_thirdjet_mass_diff_hists.At(B);
      TH1* threejet_thirdjet_perp_diff = (TH1*)threejet_thirdjet_perp_diff_hists.At(B);
      TH1* threejet_thirdjet_phi_diff = (TH1*)threejet_thirdjet_phi_diff_hists.At(B);
      TH1* threejet_thirdjet_mass_diff_wta = (TH1*)threejet_thirdjet_mass_diff_wta_hists.At(B);
      TH1* threejet_thirdjet_perp_diff_wta = (TH1*)threejet_thirdjet_perp_diff_wta_hists.At(B);
      TH1* threejet_thirdjet_phi_diff_wta = (TH1*)threejet_thirdjet_phi_diff_wta_hists.At(B);

      // threejet_firstjet_mass_njet->Fill(sorted_njet_threejets[0].m());
      // threejet_firstjet_perp_njet->Fill(sorted_njet_threejets[0].perp());
      // threejet_firstjet_phi_njet->Fill(sorted_njet_threejets[0].phi());
      // threejet_secondjet_mass_njet->Fill(sorted_njet_threejets[1].m());
      // threejet_secondjet_perp_njet->Fill(sorted_njet_threejets[1].perp());
      // threejet_secondjet_phi_njet->Fill(sorted_njet_threejets[1].phi());
      threejet_thirdjet_mass_njet->Fill(sorted_njet_threejets[2].m());
      threejet_thirdjet_perp_njet->Fill(sorted_njet_threejets[2].perp());
      threejet_thirdjet_phi_njet->Fill(sorted_njet_threejets[2].phi());

      threejet_thirdjet_mass_diff->Fill(third_jet.m() - sorted_njet_threejets[2].m());
      threejet_thirdjet_perp_diff->Fill(third_jet.perp() - sorted_njet_threejets[2].perp());
      threejet_thirdjet_phi_diff->Fill(calcDphi(third_jet, sorted_njet_threejets[2]));
      threejet_thirdjet_mass_diff_wta->Fill(third_jet_wta.m() - sorted_njet_threejets[2].m());
      threejet_thirdjet_perp_diff_wta->Fill(third_jet_wta.perp() - sorted_njet_threejets[2].perp());
      threejet_thirdjet_phi_diff_wta->Fill(calcDphi(third_jet_wta, sorted_njet_threejets[2]));
      
      NjettinessPlugin njet_plugin_manual_3jets(3, OnePass_Manual_Axes(), measure_function);
      JetDefinition njet_def_manual_3jets(&njet_plugin_manual_3jets);
      if (hardest_threejets.size() == 3) {
        njet_plugin_manual_3jets.setAxes(hardest_threejets);
        ClusterSequence njet_cluster_manual_3jets(fjInputs, njet_def_manual_3jets);
        const NjettinessExtras *extras_manual_3jets = njettiness_extras(njet_cluster_manual_3jets);
        vector<PseudoJet> njet_jets_manual_3jets = extras_manual_3jets->jets();

        TH1* tau3_diff = (TH1*)tau3_diff_hists.At(B);
        tau3_diff->Fill(extras_manual_3jets->totalTau() - extras_threejets->totalTau());
      }
      NjettinessPlugin njet_plugin_manual_wta_3jets(3, OnePass_Manual_Axes(), measure_function);
      JetDefinition njet_def_manual_wta_3jets(&njet_plugin_manual_wta_3jets);
      if (hardest_threejets.size() == 3) {
        njet_plugin_manual_wta_3jets.setAxes(hardest_threejets);
        ClusterSequence njet_cluster_manual_wta_3jets(fjInputs, njet_def_manual_wta_3jets);
        const NjettinessExtras *extras_manual_wta_3jets = njettiness_extras(njet_cluster_manual_wta_3jets);
        vector<PseudoJet> njet_jets_manual_wta_3jets = extras_manual_wta_3jets->jets();

        TH1* tau3_diff_wta = (TH1*)tau3_diff_wta_hists.At(B);
        tau3_diff_wta->Fill(extras_manual_wta_3jets->totalTau() - extras_threejets->totalTau());
      }
    }

    // event_display->SetStats(0);
    // axes_display->SetStats(0);
    // axes_njet_display->SetStats(0);
    // akt_jets_display->SetStats(0);
    // njet_jets_display->SetStats(0);

    // TCanvas *display = new TCanvas("display", "Event Display", 1093, 700);
    // display->cd();
    // display->SetFixedAspectRatio();
    // event_display->GetXaxis()->SetTitle("#eta");
    // event_display->GetYaxis()->SetTitle("#phi");
    // event_display->SetFillColor(kBlack);
    // event_display->SetLineColor(kBlack);
    // event_display->SetLineWidth(1);
    // event_display->Draw("box");

    // for (int a = 0; a < hardest_jet.size(); a++) {
    //   axes_display->Fill(hardest_jet[a].eta(), hardest_jet[a].phi());
    //   vector<PseudoJet> constituents = hardest_jet[a].constituents();
    //   for (int i_const = 0; i_const < constituents.size(); i_const++) {
    //     akt_jets_display->Fill(constituents[i_const].eta(), constituents[i_const].phi());
    //   }          
    // }
    // axes_display->SetMarkerStyle(3);
    // axes_display->SetMarkerSize(3);
    // axes_display->SetMarkerColor(kBlue);
    // akt_jets_display->SetMarkerStyle(21);
    // akt_jets_display->SetMarkerSize(0.5);
    // akt_jets_display->SetMarkerColor(kBlue);
    // axes_display->Draw("SAMES");
    // akt_jets_display->Draw("SAMES");

    // axes_njet_display->SetMarkerStyle(3);
    // axes_njet_display->SetMarkerSize(3);
    // axes_njet_display->SetMarkerColor(kRed);
    // njet_jets_display->SetMarkerStyle(21);
    // njet_jets_display->SetMarkerSize(0.5);
    // njet_jets_display->SetMarkerColor(kRed);
    // axes_njet_display->Draw("SAMES");
    // njet_jets_display->Draw("SAMES");

    // display->Write();

    //   for (int a = 0; a < hardest_jet.size(); a++) {
    //     cout << hardest_jet[a].eta() << " " << hardest_jet[a].phi() << " ";
    //   }
    //   cout << endl;
    //   for (int a = 0; a < hardest_jet.size(); a++) {
    //     cout << njet_jets[a].eta() << " " << njet_jets[a].phi() << " ";
    //   }
    //   cout << endl << endl;

    // ostringstream ss;
    // ss << i_event;
    // TString title;
    // if (i_sample == 0) title = "ttbar_display" + ss.str() + "_6jets_test.eps";
    // if (i_sample == 1) title = "dijets_display" + ss.str() + "_6jets_test.eps";

    // display->Print(title, "eps");

    delete event_display;
    delete axes_display;
    delete axes_njet_display;
    delete akt_jets_display;
    delete njet_jets_display;
    // delete display;                                                                                                                                                                                                                                                                                                                                                         

  }

 // pythia.stat();

  akt_jet_mass->Write();
  akt_jet_perp->Write();
  akt_jet_phi->Write();
  akt_jet_eta->Write();
  zprime_invmass_akt->Write();
  jet_ptspectrum_akt->Write();
  jet_phispectrum_akt->Write();

  akt_wta_jet_mass->Write();
  akt_wta_jet_perp->Write();
  akt_wta_jet_phi->Write();
  akt_wta_jet_eta->Write();
  zprime_invmass_akt_wta->Write();
  jet_ptspectrum_akt_wta->Write();
  jet_phispectrum_akt_wta->Write();

  secondjet_mass_akt->Write();
  secondjet_perp_akt->Write();
  secondjet_phi_akt->Write();

  threejet_invmass_akt->Write();
  threejet_thirdjet_mass_akt->Write();
  threejet_thirdjet_perp_akt->Write();
  threejet_thirdjet_phi_akt->Write();

  secondjet_mass_akt_wta->Write();
  secondjet_perp_akt_wta->Write();
  secondjet_phi_akt_wta->Write();

  threejet_invmass_akt_wta->Write();
  threejet_thirdjet_mass_akt_wta->Write();
  threejet_thirdjet_perp_akt_wta->Write();
  threejet_thirdjet_phi_akt_wta->Write();

  TCanvas *mass_compare_2jets = new TCanvas("mass_compare_2jets", "mass_compare_2jets", 600, 600);
  mass_compare_2jets->cd();
  double zprime_invmass_akt_scale = 1/zprime_invmass_akt->Integral(0, zprime_invmass_akt->GetNbinsX() + 1);
  zprime_invmass_akt->Scale(zprime_invmass_akt_scale);
  zprime_invmass_akt->SetLineColor(kBlack);
  zprime_invmass_akt->SetTitle("Dijet invariant mass");
  zprime_invmass_akt->GetXaxis()->SetTitle("m_{jj} (GeV)");
  zprime_invmass_akt->GetYaxis()->SetTitle("Relative Occurrence");
  zprime_invmass_akt->GetYaxis()->SetTitleOffset(1.5);
  zprime_invmass_akt->SetMinimum(0.0);
  double max_val = zprime_invmass_akt->GetMaximum();

  double zprime_invmass_akt_wta_scale = 1/zprime_invmass_akt_wta->Integral(0, zprime_invmass_akt_wta->GetNbinsX() + 1);
  zprime_invmass_akt_wta->Scale(zprime_invmass_akt_wta_scale);
  zprime_invmass_akt_wta->SetLineColor(28);
  if (zprime_invmass_akt_wta->GetMaximum() > max_val) max_val = zprime_invmass_akt_wta->GetMaximum();

  zprime_invmass_akt->Draw();
  zprime_invmass_akt_wta->Draw("SAMES");
  TLegend *leg_mass_2jets = new TLegend(0.12, 0.5, 0.5, 0.88);
  leg_mass_2jets->SetFillColor(kWhite);
  leg_mass_2jets->SetLineColor(kWhite);
  leg_mass_2jets->AddEntry(zprime_invmass_akt, "Anti-KT (E-scheme)", "L");
  leg_mass_2jets->AddEntry(zprime_invmass_akt_wta, "Anti-KT (WTA)", "L");

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* zprime_invmass_njet = (TH1*)zprime_invmass_njet_hists.At(B);
    double zprime_invmass_njet_scale = 1/zprime_invmass_njet->Integral(0, zprime_invmass_njet->GetNbinsX() + 1);
    zprime_invmass_njet->Scale(zprime_invmass_akt_scale);
    if (B == 0) zprime_invmass_njet->SetLineColor(kYellow);
    if (B == 1) zprime_invmass_njet->SetLineColor(kGreen);
    if (B == 2) zprime_invmass_njet->SetLineColor(kRed);
    if (B == 3) zprime_invmass_njet->SetLineColor(kBlue);
    zprime_invmass_njet->Draw("SAMES");
    leg_mass_2jets->AddEntry(zprime_invmass_njet, "#beta = " + (TString)ss.str());
    if (zprime_invmass_njet->GetMaximum() > max_val) max_val = zprime_invmass_njet->GetMaximum();
  }

  zprime_invmass_akt->SetMaximum(1.2*max_val);
  leg_mass_2jets->Draw("SAMES");
  mass_compare_2jets->Write();
  mass_compare_2jets->Print("zprime_mass_compare_2jets.eps", "eps");


  TCanvas *invperp_compare_2jets = new TCanvas("invperp_compare_2jets", "invperp_compare_2jets", 600, 600);
  invperp_compare_2jets->cd();
  double zprime_invperp_akt_scale = 1/zprime_invperp_akt->Integral(0, zprime_invperp_akt->GetNbinsX() + 1);
  zprime_invperp_akt->Scale(zprime_invperp_akt_scale);
  zprime_invperp_akt->SetLineColor(kBlack);
  zprime_invperp_akt->SetTitle("Dijet invariant invperp");
  zprime_invperp_akt->GetXaxis()->SetTitle("m_{jj} (GeV)");
  zprime_invperp_akt->GetYaxis()->SetTitle("Relative Occurrence");
  zprime_invperp_akt->GetYaxis()->SetTitleOffset(1.5);
  zprime_invperp_akt->SetMinimum(0.0);
  double max_val_invperp = zprime_invperp_akt->GetMaximum();

  double zprime_invperp_akt_wta_scale = 1/zprime_invperp_akt_wta->Integral(0, zprime_invperp_akt_wta->GetNbinsX() + 1);
  zprime_invperp_akt_wta->Scale(zprime_invperp_akt_wta_scale);
  zprime_invperp_akt_wta->SetLineColor(28);
  if (zprime_invperp_akt_wta->GetMaximum() > max_val_invperp) max_val_invperp = zprime_invperp_akt_wta->GetMaximum();

  zprime_invperp_akt->Draw();
  zprime_invperp_akt_wta->Draw("SAMES");
  TLegend *leg_invperp_2jets = new TLegend(0.12, 0.5, 0.5, 0.88);
  leg_invperp_2jets->SetFillColor(kWhite);
  leg_invperp_2jets->SetLineColor(kWhite);

  leg_invperp_2jets->AddEntry(zprime_invperp_akt, "Anti-KT (E-scheme)", "L");
  leg_invperp_2jets->AddEntry(zprime_invperp_akt_wta, "Anti-KT (WTA)", "L");

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* zprime_invperp_njet = (TH1*)zprime_invperp_njet_hists.At(B);
    double zprime_invperp_njet_scale = 1/zprime_invperp_njet->Integral(0, zprime_invperp_njet->GetNbinsX() + 1);
    zprime_invperp_njet->Scale(zprime_invperp_akt_scale);
    if (B == 0) zprime_invperp_njet->SetLineColor(kYellow);
    if (B == 1) zprime_invperp_njet->SetLineColor(kGreen);
    if (B == 2) zprime_invperp_njet->SetLineColor(kRed);
    if (B == 3) zprime_invperp_njet->SetLineColor(kBlue);
    zprime_invperp_njet->Draw("SAMES");
    leg_invperp_2jets->AddEntry(zprime_invperp_njet, "#beta = " + (TString)ss.str());
    if (zprime_invperp_njet->GetMaximum() > max_val_invperp) max_val_invperp = zprime_invperp_njet->GetMaximum();
  }

  zprime_invperp_akt->SetMaximum(1.2*max_val_invperp);
  leg_invperp_2jets->Draw("SAMES");
  invperp_compare_2jets->Write();
  invperp_compare_2jets->Print("zprime_invperp_compare_2jets.eps", "eps");


  TCanvas *perp_compare_2jets = new TCanvas("perp_compare_2jets", "perp_compare_2jets", 600, 600);
  perp_compare_2jets->cd();
  double jet_ptspectrum_akt_scale = 1/jet_ptspectrum_akt->Integral(0, jet_ptspectrum_akt->GetNbinsX() + 1);
  jet_ptspectrum_akt->Scale(jet_ptspectrum_akt_scale);
  jet_ptspectrum_akt->SetLineColor(kBlack);
  jet_ptspectrum_akt->SetTitle("pT of 2 hardest jets");
  jet_ptspectrum_akt->GetXaxis()->SetTitle("pT_{j} (GeV)");
  jet_ptspectrum_akt->GetYaxis()->SetTitle("Relative Occurrence");
  jet_ptspectrum_akt->GetYaxis()->SetTitleOffset(1.5);
  jet_ptspectrum_akt->SetMinimum(0.0);
  double max_val_perp_2jets = jet_ptspectrum_akt->GetMaximum();

  double jet_ptspectrum_akt_wta_scale = 1/jet_ptspectrum_akt_wta->Integral(0, jet_ptspectrum_akt_wta->GetNbinsX() + 1);
  jet_ptspectrum_akt_wta->Scale(jet_ptspectrum_akt_wta_scale);
  jet_ptspectrum_akt_wta->SetLineColor(28);
  if (jet_ptspectrum_akt_wta->GetMaximum() > max_val_perp_2jets) max_val_perp_2jets = jet_ptspectrum_akt_wta->GetMaximum();

  jet_ptspectrum_akt->Draw();
  jet_ptspectrum_akt_wta->Draw("SAMES");
  TLegend *leg_perp_2jets = new TLegend(0.5, 0.5, 0.88, 0.88);
  leg_perp_2jets->SetFillColor(kWhite);
  leg_perp_2jets->SetLineColor(kWhite);
  leg_perp_2jets->AddEntry(jet_ptspectrum_akt, "Anti-KT (E-scheme)", "L");
  leg_perp_2jets->AddEntry(jet_ptspectrum_akt_wta, "Anti-KT (WTA)", "L");

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* jet_ptspectrum_njet = (TH1*)jet_ptspectrum_njet_hists.At(B);
    double jet_ptspectrum_njet_scale = 1/jet_ptspectrum_njet->Integral(0, jet_ptspectrum_njet->GetNbinsX() + 1);
    jet_ptspectrum_njet->Scale(jet_ptspectrum_akt_scale);
    if (B == 0) jet_ptspectrum_njet->SetLineColor(kYellow);
    if (B == 1) jet_ptspectrum_njet->SetLineColor(kGreen);
    if (B == 2) jet_ptspectrum_njet->SetLineColor(kRed);
    if (B == 3) jet_ptspectrum_njet->SetLineColor(kBlue);
    jet_ptspectrum_njet->Draw("SAMES");
    leg_perp_2jets->AddEntry(jet_ptspectrum_njet, "#beta = " + (TString)ss.str());
    if (jet_ptspectrum_njet->GetMaximum() > max_val_perp_2jets) max_val_perp_2jets = jet_ptspectrum_njet->GetMaximum();
  }

  jet_ptspectrum_akt->SetMaximum(1.2*max_val_perp_2jets);
  leg_perp_2jets->Draw("SAMES");
  perp_compare_2jets->Write();
  perp_compare_2jets->Print("zprime_perp_compare_2jets.eps", "eps");


  TCanvas *phi_compare_2jets = new TCanvas("phi_compare_2jets", "phi_compare_2jets", 600, 600);
  phi_compare_2jets->cd();
  double jet_phispectrum_akt_scale = 1/jet_phispectrum_akt->Integral(0, jet_phispectrum_akt->GetNbinsX() + 1);
  jet_phispectrum_akt->Scale(jet_phispectrum_akt_scale);
  jet_phispectrum_akt->SetLineColor(kBlack);
  jet_phispectrum_akt->SetTitle("Azimuthal difference between 2 hardest jets");
  jet_phispectrum_akt->GetXaxis()->SetTitle("#phi_{j} (GeV)");
  jet_phispectrum_akt->GetYaxis()->SetTitle("Relative Occurrence");
  jet_phispectrum_akt->GetYaxis()->SetTitleOffset(1.5);
  jet_phispectrum_akt->SetMinimum(0.0);
  double max_val_phi_2jets = jet_phispectrum_akt->GetMaximum();

  double jet_phispectrum_akt_wta_scale = 1/jet_phispectrum_akt_wta->Integral(0, jet_phispectrum_akt_wta->GetNbinsX() + 1);
  jet_phispectrum_akt_wta->Scale(jet_phispectrum_akt_wta_scale);
  jet_phispectrum_akt_wta->SetLineColor(28);
  if (jet_phispectrum_akt_wta->GetMaximum() > max_val_phi_2jets) max_val_phi_2jets = jet_phispectrum_akt_wta->GetMaximum();

  jet_phispectrum_akt->Draw();
  jet_phispectrum_akt_wta->Draw("SAMES");
  TLegend *leg_phi_2jets = new TLegend(0.12, 0.5, 0.5, 0.88);
  leg_phi_2jets->SetFillColor(kWhite);
  leg_phi_2jets->SetLineColor(kWhite);
  leg_phi_2jets->AddEntry(jet_phispectrum_akt, "Anti-KT (E-scheme)", "L");
  leg_phi_2jets->AddEntry(jet_phispectrum_akt_wta, "Anti-KT (WTA)", "L");

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* jet_phispectrum_njet = (TH1*)jet_phispectrum_njet_hists.At(B);
    double jet_phispectrum_njet_scale = 1/jet_phispectrum_njet->Integral(0, jet_phispectrum_njet->GetNbinsX() + 1);
    jet_phispectrum_njet->Scale(jet_phispectrum_akt_scale);
    if (B == 0) jet_phispectrum_njet->SetLineColor(kYellow);
    if (B == 1) jet_phispectrum_njet->SetLineColor(kGreen);
    if (B == 2) jet_phispectrum_njet->SetLineColor(kRed);
    if (B == 3) jet_phispectrum_njet->SetLineColor(kBlue);
    jet_phispectrum_njet->Draw("SAMES");
    leg_phi_2jets->AddEntry(jet_phispectrum_njet, "#beta = " + (TString)ss.str());
    if (jet_phispectrum_njet->GetMaximum() > max_val_phi_2jets) max_val_phi_2jets = jet_phispectrum_njet->GetMaximum();
  }

  jet_phispectrum_akt->SetMaximum(1.2*max_val_phi_2jets);
  leg_phi_2jets->SetFillColor(kWhite);
  leg_phi_2jets->SetLineColor(kWhite);
  leg_phi_2jets->Draw("SAMES");
  phi_compare_2jets->Write();
  phi_compare_2jets->Print("zprime_phi_compare_2jets.eps", "eps");


  // TCanvas *mass_compare_secondjet = new TCanvas("mass_compare_secondjet", "mass_compare_secondjet", 600, 600);
  // mass_compare_secondjet->cd();
  // double secondjet_mass_akt_scale = 1/secondjet_mass_akt->Integral(0, secondjet_mass_akt->GetNbinsX() + 1);
  // secondjet_mass_akt->Scale(secondjet_mass_akt_scale);
  // secondjet_mass_akt->SetLineColor(kBlack);
  // secondjet_mass_akt->SetTitle("Mass of second jet");
  // secondjet_mass_akt->GetXaxis()->SetTitle("m_{j} (GeV)");
  // secondjet_mass_akt->GetYaxis()->SetTitle("Relative Occurrence");
  // secondjet_mass_akt->GetYaxis()->SetTitleOffset(1.5);
  // secondjet_mass_akt->SetMinimum(0.0);
  // double max_val_mass_secondjet = secondjet_mass_akt->GetMaximum();

  // double secondjet_mass_akt_wta_scale = 1/secondjet_mass_akt_wta->Integral(0, secondjet_mass_akt_wta->GetNbinsX() + 1);
  // secondjet_mass_akt_wta->Scale(secondjet_mass_akt_wta_scale);
  // secondjet_mass_akt_wta->SetLineColor(28);
  // if (secondjet_mass_akt_wta->GetMaximum() > max_val_mass_secondjet) max_val_mass_secondjet = secondjet_mass_akt_wta->GetMaximum();

  // secondjet_mass_akt->Draw();
  // secondjet_mass_akt_wta->Draw("SAMES");
  // TLegend *leg_mass_secondjet = new TLegend(0.5, 0.5, 0.9, 0.9);
  // leg_mass_secondjet->AddEntry(secondjet_mass_akt, "Anti-KT (E-scheme)", "L");
  // leg_mass_secondjet->AddEntry(secondjet_mass_akt_wta, "Anti-KT (WTA)", "L");

  // for (int B = 0; B < n_betas; B++) {

  //   double beta = betalist[B];
  //   ostringstream ss;
  //   ss << beta;

  //   TH1* secondjet_mass_njet = (TH1*)secondjet_mass_njet_hists.At(B);
  //   double secondjet_mass_njet_scale = 1/secondjet_mass_njet->Integral(0, secondjet_mass_njet->GetNbinsX() + 1);
  //   secondjet_mass_njet->Scale(secondjet_mass_akt_scale);
  //   if (B == 0) secondjet_mass_njet->SetLineColor(kYellow);
  //   if (B == 1) secondjet_mass_njet->SetLineColor(kGreen);
  //   if (B == 2) secondjet_mass_njet->SetLineColor(kRed);
  //   if (B == 3) secondjet_mass_njet->SetLineColor(kBlue);
  //   secondjet_mass_njet->Draw("SAMES");
  //   leg_mass_secondjet->AddEntry(secondjet_mass_njet, "#beta = " + (TString)ss.str());
  //   if (secondjet_mass_njet->GetMaximum() > max_val_mass_secondjet) max_val_mass_secondjet = secondjet_mass_njet->GetMaximum();
  // }

  // secondjet_mass_akt->SetMaximum(1.2*max_val_mass_secondjet);
  // leg_mass_secondjet->Draw("SAMES");
  // mass_compare_secondjet->Write();
  // mass_compare_secondjet->Print("mass_compare_secondjet.eps", "eps");


  // TCanvas *perp_compare_secondjet = new TCanvas("perp_compare_secondjet", "perp_compare_secondjet", 600, 600);
  // perp_compare_secondjet->cd();
  // double secondjet_perp_akt_scale = 1/secondjet_perp_akt->Integral(0, secondjet_perp_akt->GetNbinsX() + 1);
  // secondjet_perp_akt->Scale(secondjet_perp_akt_scale);
  // secondjet_perp_akt->SetLineColor(kBlack);
  // secondjet_perp_akt->SetTitle("pT of second jet");
  // secondjet_perp_akt->GetXaxis()->SetTitle("pT_{j} (GeV)");
  // secondjet_perp_akt->GetYaxis()->SetTitle("Relative Occurrence");
  // secondjet_perp_akt->GetYaxis()->SetTitleOffset(1.5);
  // secondjet_perp_akt->SetMinimum(0.0);
  // double max_val_perp_secondjet = secondjet_perp_akt->GetMaximum();

  // double secondjet_perp_akt_wta_scale = 1/secondjet_perp_akt_wta->Integral(0, secondjet_perp_akt_wta->GetNbinsX() + 1);
  // secondjet_perp_akt_wta->Scale(secondjet_perp_akt_wta_scale);
  // secondjet_perp_akt_wta->SetLineColor(28);
  // if (secondjet_perp_akt_wta->GetMaximum() > max_val_perp_secondjet) max_val_perp_secondjet = secondjet_perp_akt_wta->GetMaximum();

  // secondjet_perp_akt->Draw();
  // secondjet_perp_akt_wta->Draw("SAMES");
  // TLegend *leg_perp_secondjet = new TLegend(0.5, 0.5, 0.9, 0.9);
  // leg_perp_secondjet->AddEntry(secondjet_perp_akt, "Anti-KT (E-scheme)", "L");
  // leg_perp_secondjet->AddEntry(secondjet_perp_akt_wta, "Anti-KT (WTA)", "L");

  // for (int B = 0; B < n_betas; B++) {

  //   double beta = betalist[B];
  //   ostringstream ss;
  //   ss << beta;

  //   TH1* secondjet_perp_njet = (TH1*)secondjet_perp_njet_hists.At(B);
  //   double secondjet_perp_njet_scale = 1/secondjet_perp_njet->Integral(0, secondjet_perp_njet->GetNbinsX() + 1);
  //   secondjet_perp_njet->Scale(secondjet_perp_akt_scale);
  //   if (B == 0) secondjet_perp_njet->SetLineColor(kYellow);
  //   if (B == 1) secondjet_perp_njet->SetLineColor(kGreen);
  //   if (B == 2) secondjet_perp_njet->SetLineColor(kRed);
  //   if (B == 3) secondjet_perp_njet->SetLineColor(kBlue);
  //   secondjet_perp_njet->Draw("SAMES");
  //   leg_perp_secondjet->AddEntry(secondjet_perp_njet, "#beta = " + (TString)ss.str());
  //   if (secondjet_perp_njet->GetMaximum() > max_val_perp_secondjet) max_val_perp_secondjet = secondjet_perp_njet->GetMaximum();
  // }

  // secondjet_perp_akt->SetMaximum(1.2*max_val_perp_secondjet);
  // leg_perp_secondjet->Draw("SAMES");
  // perp_compare_secondjet->Write();
  // perp_compare_secondjet->Print("perp_compare_secondjet.eps", "eps");


  // TCanvas *phi_compare_secondjet = new TCanvas("phi_compare_secondjet", "phi_compare_secondjet", 600, 600);
  // phi_compare_secondjet->cd();
  // double secondjet_phi_akt_scale = 1/secondjet_phi_akt->Integral(0, secondjet_phi_akt->GetNbinsX() + 1);
  // secondjet_phi_akt->Scale(secondjet_phi_akt_scale);
  // secondjet_phi_akt->SetLineColor(kBlack);
  // secondjet_phi_akt->SetTitle("Azimuth of second jet");
  // secondjet_phi_akt->GetXaxis()->SetTitle("#phi_{j} (GeV)");
  // secondjet_phi_akt->GetYaxis()->SetTitle("Relative Occurrence");
  // secondjet_phi_akt->GetYaxis()->SetTitleOffset(1.5);
  // secondjet_phi_akt->SetMinimum(0.0);
  // double max_val_phi_secondjet = secondjet_phi_akt->GetMaximum();

  // double secondjet_phi_akt_wta_scale = 1/secondjet_phi_akt_wta->Integral(0, secondjet_phi_akt_wta->GetNbinsX() + 1);
  // secondjet_phi_akt_wta->Scale(secondjet_phi_akt_wta_scale);
  // secondjet_phi_akt_wta->SetLineColor(28);
  // if (secondjet_phi_akt_wta->GetMaximum() > max_val_phi_secondjet) max_val_phi_secondjet = secondjet_phi_akt_wta->GetMaximum();

  // secondjet_phi_akt->Draw();
  // secondjet_phi_akt_wta->Draw("SAMES");
  // TLegend *leg_phi_secondjet = new TLegend(0.5, 0.5, 0.9, 0.9);
  // leg_phi_secondjet->AddEntry(secondjet_phi_akt, "Anti-KT (E-scheme)", "L");
  // leg_phi_secondjet->AddEntry(secondjet_phi_akt_wta, "Anti-KT (WTA)", "L");

  // for (int B = 0; B < n_betas; B++) {

  //   double beta = betalist[B];
  //   ostringstream ss;
  //   ss << beta;

  //   TH1* secondjet_phi_njet = (TH1*)secondjet_phi_njet_hists.At(B);
  //   double secondjet_phi_njet_scale = 1/secondjet_phi_njet->Integral(0, secondjet_phi_njet->GetNbinsX() + 1);
  //   secondjet_phi_njet->Scale(secondjet_phi_akt_scale);
  //   if (B == 0) secondjet_phi_njet->SetLineColor(kYellow);
  //   if (B == 1) secondjet_phi_njet->SetLineColor(kGreen);
  //   if (B == 2) secondjet_phi_njet->SetLineColor(kRed);
  //   if (B == 3) secondjet_phi_njet->SetLineColor(kBlue);
  //   secondjet_phi_njet->Draw("SAMES");
  //   leg_phi_secondjet->AddEntry(secondjet_phi_njet, "#beta = " + (TString)ss.str());
  //   if (secondjet_phi_njet->GetMaximum() > max_val_phi_secondjet) max_val_phi_secondjet = secondjet_phi_njet->GetMaximum();
  // }

  // secondjet_phi_akt->SetMaximum(1.2*max_val_phi_secondjet);
  // leg_phi_secondjet->Draw("SAMES");
  // phi_compare_secondjet->Write();
  // phi_compare_secondjet->Print("phi_compare_secondjet.eps", "eps");

  TCanvas *massdiff_compare_2jets = new TCanvas("massdiff_compare_2jets", "massdiff_compare_2jets", 600, 600);
  massdiff_compare_2jets->cd();
  massdiff_compare_2jets->SetLogy();
  double max_val_massdiff_2jets = 0;
  TLegend *leg_massdiff_2jets = new TLegend(0.55, 0.55, 0.88, 0.88);
  leg_massdiff_2jets->SetFillColor(kWhite);
  leg_massdiff_2jets->SetLineColor(kWhite);


  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* zprime_invmass_diff = (TH1*)zprime_invmass_diff_hists.At(B);
    double zprime_invmass_diff_scale = 1/zprime_invmass_diff->Integral(0, zprime_invmass_diff->GetNbinsX() + 1);
    zprime_invmass_diff->Scale(zprime_invmass_diff_scale);
    if (B == 0) zprime_invmass_diff->SetLineColor(kYellow);
    if (B == 1) zprime_invmass_diff->SetLineColor(kGreen);
    if (B == 2) zprime_invmass_diff->SetLineColor(kRed);
    if (B == 3) zprime_invmass_diff->SetLineColor(kBlue);

    if (B == 0) zprime_invmass_diff->Draw();
    else zprime_invmass_diff->Draw("SAMES");
    zprime_invmass_diff->SetTitle("Difference in invariant mass from anti-KT (E-scheme) and 2-jettiness");
    zprime_invmass_diff->GetXaxis()->SetTitle("m_{jj}_{akt} - m_{jj}_{njet} (GeV)");
    zprime_invmass_diff->GetYaxis()->SetTitle("Relative Occurrence");
    zprime_invmass_diff->GetYaxis()->SetTitleOffset(1.5);

    leg_massdiff_2jets->AddEntry(zprime_invmass_diff, "#beta = " + (TString)ss.str());
    if (zprime_invmass_diff->GetMaximum() > max_val_massdiff_2jets) {
      max_val_massdiff_2jets = zprime_invmass_diff->GetMaximum();
      zprime_invmass_diff->SetMaximum(1.2*max_val_massdiff_2jets);
    }
  }

  leg_massdiff_2jets->Draw("SAMES");
  massdiff_compare_2jets->Write();
  massdiff_compare_2jets->Print("zprime_massdiff_compare_2jets.eps", "eps");


  TCanvas *perpdiff_compare_2jets = new TCanvas("perpdiff_compare_2jets", "perpdiff_compare_2jets", 600, 600);
  perpdiff_compare_2jets->cd();
  perpdiff_compare_2jets->SetLogy();
  double max_val_perpdiff_2jets = 0;
  TLegend *leg_perpdiff_2jets = new TLegend(0.5, 0.5, 0.88, 0.88);
  leg_perpdiff_2jets->SetFillColor(kWhite);
  leg_perpdiff_2jets->SetLineColor(kWhite);

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* twojet_bothjets_perp_diff = (TH1*)twojet_bothjets_perp_diff_hists.At(B);
    double twojet_bothjets_perp_diff_scale = 1/twojet_bothjets_perp_diff->Integral(0, twojet_bothjets_perp_diff->GetNbinsX() + 1);
    twojet_bothjets_perp_diff->Scale(twojet_bothjets_perp_diff_scale);
    if (B == 0) twojet_bothjets_perp_diff->SetLineColor(kYellow);
    if (B == 1) twojet_bothjets_perp_diff->SetLineColor(kGreen);
    if (B == 2) twojet_bothjets_perp_diff->SetLineColor(kRed);
    if (B == 3) twojet_bothjets_perp_diff->SetLineColor(kBlue);

    if (B == 0) twojet_bothjets_perp_diff->Draw();
    else twojet_bothjets_perp_diff->Draw("SAMES");
    twojet_bothjets_perp_diff->SetTitle("Difference in pT between 2 hardest jets from anti-KT (E-scheme) and 2-jettiness");
    twojet_bothjets_perp_diff->GetXaxis()->SetTitle("pT_{j}_{akt} - pT_{j}_{njet} (GeV)");
    twojet_bothjets_perp_diff->GetYaxis()->SetTitle("Relative Occurrence");
    twojet_bothjets_perp_diff->GetYaxis()->SetTitleOffset(1.5);

    leg_perpdiff_2jets->AddEntry(twojet_bothjets_perp_diff, "#beta = " + (TString)ss.str());
    if (twojet_bothjets_perp_diff->GetMaximum() > max_val_perpdiff_2jets) {
      max_val_perpdiff_2jets = twojet_bothjets_perp_diff->GetMaximum();
      twojet_bothjets_perp_diff->SetMaximum(1.2*max_val_perpdiff_2jets);
    }
  }

  leg_perpdiff_2jets->Draw("SAMES");
  perpdiff_compare_2jets->Write();
  perpdiff_compare_2jets->Print("zprime_perpdiff_compare_2jets.eps", "eps");


  TCanvas *phidiff_compare_2jets = new TCanvas("phidiff_compare_2jets", "phidiff_compare_2jets", 600, 600);
  phidiff_compare_2jets->cd();
  phidiff_compare_2jets->SetLogy();
  double max_val_phidiff_2jets = 0;
  TLegend *leg_phidiff_2jets = new TLegend(0.5, 0.5, 0.88, 0.88);
  leg_phidiff_2jets->SetFillColor(kWhite);
  leg_phidiff_2jets->SetLineColor(kWhite);

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* twojet_bothjets_phi_diff = (TH1*)twojet_bothjets_phi_diff_hists.At(B);
    double twojet_bothjets_phi_diff_scale = 1/twojet_bothjets_phi_diff->Integral(0, twojet_bothjets_phi_diff->GetNbinsX() + 1);
    twojet_bothjets_phi_diff->Scale(twojet_bothjets_phi_diff_scale);
    if (B == 0) twojet_bothjets_phi_diff->SetLineColor(kYellow);
    if (B == 1) twojet_bothjets_phi_diff->SetLineColor(kGreen);
    if (B == 2) twojet_bothjets_phi_diff->SetLineColor(kRed);
    if (B == 3) twojet_bothjets_phi_diff->SetLineColor(kBlue);

    if (B == 0) twojet_bothjets_phi_diff->Draw();
    else twojet_bothjets_phi_diff->Draw("SAMES");
    twojet_bothjets_phi_diff->SetTitle("Difference in pT between 2 hardest jets from anti-KT (E-scheme) and 2-jettiness");
    twojet_bothjets_phi_diff->GetXaxis()->SetTitle("pT_{j}_{akt} - pT_{j}_{njet} (GeV)");
    twojet_bothjets_phi_diff->GetYaxis()->SetTitle("Relative Occurrence");
    twojet_bothjets_phi_diff->GetYaxis()->SetTitleOffset(1.5);

    leg_phidiff_2jets->AddEntry(twojet_bothjets_phi_diff, "#beta = " + (TString)ss.str());
    if (twojet_bothjets_phi_diff->GetMaximum() > max_val_phidiff_2jets) {
      max_val_phidiff_2jets = twojet_bothjets_phi_diff->GetMaximum();
      twojet_bothjets_phi_diff->SetMaximum(1.2*max_val_phidiff_2jets);
    }
  }

  leg_phidiff_2jets->Draw("SAMES");
  phidiff_compare_2jets->Write();
  phidiff_compare_2jets->Print("zprime_phidiff_compare_2jets.eps", "eps");


  // TCanvas *phidiff_compare_2jets = new TCanvas("phidiff_compare_2jets", "phidiff_compare_2jets", 600, 600);
  // phidiff_compare_2jets->cd();
  // phidiff_compare_2jets->SetLogy();
  // double max_val_phidiff_2jets = 0;
  // TLegend *leg_phidiff_2jets = new TLegend(0.5, 0.5, 0.88, 0.88);
  // leg_phidiff_2jets->SetFillColor(kWhite);
  // leg_phidiff_2jets->SetLineColor(kWhite);

  // for (int B = 0; B < n_betas; B++) {

  //   double beta = betalist[B];
  //   ostringstream ss;
  //   ss << beta;

  //   TH1* jet_phispectrum_diff = (TH1*)jet_phispectrum_diff_hists.At(B);
  //   double jet_phispectrum_diff_scale = 1/jet_phispectrum_diff->Integral(0, jet_phispectrum_diff->GetNbinsX() + 1);
  //   jet_phispectrum_diff->Scale(jet_phispectrum_diff_scale);
  //   if (B == 0) jet_phispectrum_diff->SetLineColor(kYellow);
  //   if (B == 1) jet_phispectrum_diff->SetLineColor(kGreen);
  //   if (B == 2) jet_phispectrum_diff->SetLineColor(kRed);
  //   if (B == 3) jet_phispectrum_diff->SetLineColor(kBlue);
  //   jet_phispectrum_diff->Draw("SAMES");
  //   jet_phispectrum_diff->SetTitle("Difference in Dijet Azimuthal Distance for Anti-KT (E-scheme) and 2-jettiness");
  //   jet_phispectrum_diff->GetXaxis()->SetTitle("|#phi_{jj}_{akt} - #phi_{jj}_{njet}| (GeV)");
  //   jet_phispectrum_diff->GetYaxis()->SetTitle("Relative Occurrence");
  //   jet_phispectrum_diff->GetYaxis()->SetTitleOffset(1.5);

  //   leg_phidiff_2jets->AddEntry(jet_phispectrum_diff, "#beta = " + (TString)ss.str());
  //   if (jet_phispectrum_diff->GetMaximum() > max_val_phidiff_2jets) {
  //     max_val_phidiff_2jets = jet_phispectrum_diff->GetMaximum();
  //     jet_phispectrum_diff->SetMaximum(1.2*max_val_phidiff_2jets);
  //   }
  // }

  // leg_phidiff_2jets->Draw("SAMES");
  // phidiff_compare_2jets->Write();
  // phidiff_compare_2jets->Print("zprime_phidiff_compare_2jets.eps", "eps");

  TCanvas *tau2diff_compare_2jets = new TCanvas("tau2diff_compare_2jets", "tau2diff_compare_2jets", 600, 600);
  tau2diff_compare_2jets->cd();
  tau2diff_compare_2jets->SetLogy();
  double max_val_tau2diff = 0;
  TLegend *leg_tau2diff_2jets = new TLegend(0.5, 0.5, 0.9, 0.9);
  leg_tau2diff_2jets->SetFillColor(kWhite);
  leg_tau2diff_2jets->SetLineColor(kWhite);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* tau2_diff = (TH1*)tau2_diff_hists.At(B);
    double tau2_diff_scale = 1/tau2_diff->Integral(0, tau2_diff->GetNbinsX() + 1);
    tau2_diff->Scale(tau2_diff_scale);
    if (B == 0) tau2_diff->SetLineColor(kYellow);
    if (B == 1) tau2_diff->SetLineColor(kGreen);
    if (B == 2) tau2_diff->SetLineColor(kRed);
    if (B == 3) tau2_diff->SetLineColor(kBlue);
    tau2_diff->SetTitle("#tau_{2} difference of anti-kT (E-scheme) and 2-jettiness jets");
    tau2_diff->GetXaxis()->SetTitle("#tau_{2}_{akt} - #tau_{2}_{njet}");
    tau2_diff->GetYaxis()->SetTitle("Relative Occurrence");
    tau2_diff->GetYaxis()->SetTitleOffset(1.5);
    tau2_diff->Draw("SAMES");
    leg_tau2diff_2jets->AddEntry(tau2_diff, "#beta = " + (TString)ss.str());
    if (tau2_diff->GetMaximum() > max_val_tau2diff) {
      max_val_tau2diff = tau2_diff->GetMaximum();
      tau2_diff->SetMaximum(1.2*max_val_tau2diff);
    }
  }

  leg_tau2diff_2jets->Draw("SAMES");
  tau2diff_compare_2jets->Write();
  tau2diff_compare_2jets->Print("zprime_tau2diff_compare_2jets.eps", "eps");


  TCanvas *mass_compare_1jet = new TCanvas("mass_compare_1jet", "mass_compare_1jet", 600, 600);
  mass_compare_1jet->cd();

  double akt_jet_mass_scale = 1/akt_jet_mass->Integral(0, akt_jet_mass->GetNbinsX() + 1);
  akt_jet_mass->Scale(akt_jet_mass_scale);
  akt_jet_mass->SetLineColor(kBlack);
  akt_jet_mass->SetTitle("Mass of Individual Jet");
  akt_jet_mass->GetXaxis()->SetTitle("m_{j} (GeV)");
  akt_jet_mass->GetYaxis()->SetTitle("Relative Occurrence");
  akt_jet_mass->GetYaxis()->SetTitleOffset(1.5);
  akt_jet_mass->SetMinimum(0.0);
  double max_val_mass = akt_jet_mass->GetMaximum();

  double akt_wta_jet_mass_scale = 1/akt_wta_jet_mass->Integral(0, akt_jet_mass->GetNbinsX() + 1);
  akt_wta_jet_mass->Scale(akt_wta_jet_mass_scale);
  akt_wta_jet_mass->SetLineColor(28);
  if (akt_wta_jet_mass->GetMaximum() > max_val) max_val = akt_wta_jet_mass->GetMaximum();

  akt_jet_mass->Draw();
  akt_wta_jet_mass->Draw("SAMES");
  TLegend *leg_mass_1jet = new TLegend(0.5, 0.5, 0.88, 0.88);
  leg_mass_1jet->SetFillColor(kWhite);
  leg_mass_1jet->SetLineColor(kWhite);
  leg_mass_1jet->AddEntry(akt_jet_mass, "Anti-KT (E-scheme)", "L");
  leg_mass_1jet->AddEntry(akt_wta_jet_mass, "Anti-KT (WTA)", "L");

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* njet_jet_mass = (TH1*)njet_jet_mass_hists.At(B);
    double njet_jet_mass_scale = 1/njet_jet_mass->Integral(0, njet_jet_mass->GetNbinsX() + 1);
    njet_jet_mass->Scale(njet_jet_mass_scale);
    if (B == 0) njet_jet_mass->SetLineColor(kYellow);
    if (B == 1) njet_jet_mass->SetLineColor(kGreen);
    if (B == 2) njet_jet_mass->SetLineColor(kRed);
    if (B == 3) njet_jet_mass->SetLineColor(kBlue);
    njet_jet_mass->Draw("SAMES");
    leg_mass_1jet->AddEntry(njet_jet_mass, "#beta = " + (TString)ss.str());
    if (njet_jet_mass->GetMaximum() > max_val_mass) max_val_mass = njet_jet_mass->GetMaximum();
  }

  akt_jet_mass->SetMaximum(1.2*max_val_mass);
  leg_mass_1jet->Draw("SAMES");
  mass_compare_1jet->Write();
  mass_compare_1jet->Print("zprime_mass_compare_1jet.eps", "eps");

  TCanvas *perp_compare_1jet = new TCanvas("perp_compare_1jet", "perp_compare_1jet", 600, 600);
  perp_compare_1jet->cd();
  double akt_jet_perp_scale = 1/akt_jet_perp->Integral(0, akt_jet_perp->GetNbinsX() + 1);
  akt_jet_perp->Scale(akt_jet_perp_scale);
  akt_jet_perp->SetLineColor(kBlack);
  akt_jet_perp->SetTitle("pT of Individual Jet");
  akt_jet_perp->GetXaxis()->SetTitle("pT_{j} (GeV)");
  akt_jet_perp->GetYaxis()->SetTitle("Relative Occurrence");
  akt_jet_perp->GetYaxis()->SetTitleOffset(1.6);
  akt_jet_perp->SetMinimum(0.0);
  double max_val_perp = akt_jet_perp->GetMaximum();

  double akt_wta_jet_perp_scale = 1/akt_wta_jet_perp->Integral(0, akt_wta_jet_perp->GetNbinsX() + 1);
  akt_wta_jet_perp->Scale(akt_wta_jet_perp_scale);
  akt_wta_jet_perp->SetLineColor(28);
  if (akt_wta_jet_perp->GetMaximum() > max_val) max_val = akt_wta_jet_perp->GetMaximum();

  akt_jet_perp->Draw();
  akt_wta_jet_perp->Draw("SAMES");
  TLegend *leg_perp_1jet = new TLegend(0.52, 0.52, 0.88, 0.88);
  leg_perp_1jet->SetFillColor(kWhite);
  leg_perp_1jet->SetLineColor(kWhite);
  leg_perp_1jet->AddEntry(akt_jet_perp, "Anti-KT (E-scheme)", "L");
  leg_perp_1jet->AddEntry(akt_wta_jet_perp, "Anti-KT (WTA)", "L");


  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* njet_jet_perp = (TH1*)njet_jet_perp_hists.At(B);
    double njet_jet_perp_scale = 1/njet_jet_perp->Integral(0, njet_jet_perp->GetNbinsX() + 1);
    njet_jet_perp->Scale(njet_jet_perp_scale);
    if (B == 0) njet_jet_perp->SetLineColor(kYellow);
    if (B == 1) njet_jet_perp->SetLineColor(kGreen);
    if (B == 2) njet_jet_perp->SetLineColor(kRed);
    if (B == 3) njet_jet_perp->SetLineColor(kBlue);
    njet_jet_perp->Draw("SAMES");
    leg_perp_1jet->AddEntry(njet_jet_perp, "#beta = " + (TString)ss.str());
    if (njet_jet_perp->GetMaximum() > max_val_perp) max_val_perp = njet_jet_perp->GetMaximum();
  }

  akt_jet_perp->SetMaximum(1.2*max_val_perp);
  leg_perp_1jet->Draw("SAMES");
  perp_compare_1jet->Write();
  perp_compare_1jet->Print("zprime_perp_compare_1jet.eps", "eps");


  TCanvas *phi_compare_1jet = new TCanvas("phi_compare_1jet", "phi_compare_1jet", 600, 600);
  phi_compare_1jet->cd();
  double akt_jet_phi_scale = 1/akt_jet_phi->Integral(0, akt_jet_phi->GetNbinsX() + 1);
  akt_jet_phi->Scale(akt_jet_phi_scale);
  akt_jet_phi->SetLineColor(kBlack);
  akt_jet_phi->SetTitle("Azimuth of Individual Jet");
  akt_jet_phi->GetXaxis()->SetTitle("#phi_{j}");
  akt_jet_phi->GetYaxis()->SetTitle("Relative Occurrence");
  akt_jet_phi->GetYaxis()->SetTitleOffset(1.6);
  akt_jet_phi->SetMinimum(0.0);
  double max_val_phi = akt_jet_phi->GetMaximum();

  double akt_wta_jet_phi_scale = 1/akt_wta_jet_phi->Integral(0, akt_jet_phi->GetNbinsX() + 1);
  akt_wta_jet_phi->Scale(akt_wta_jet_phi_scale);
  akt_wta_jet_phi->SetLineColor(28);
  if (akt_wta_jet_phi->GetMaximum() > max_val) max_val = akt_wta_jet_phi->GetMaximum();

  akt_jet_phi->Draw();
  akt_wta_jet_phi->Draw("SAMES");
  TLegend *leg_phi_1jet = new TLegend(0.12, 0.12, 0.5, 0.5);
  leg_phi_1jet->SetFillColor(kWhite);
  leg_phi_1jet->SetLineColor(kWhite);
  leg_phi_1jet->AddEntry(akt_jet_phi, "Anti-KT (E-scheme)", "L");
  leg_phi_1jet->AddEntry(akt_wta_jet_phi, "Anti-KT (WTA)", "L");

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* njet_jet_phi = (TH1*)njet_jet_phi_hists.At(B);
    double njet_jet_phi_scale = 1/njet_jet_phi->Integral(0, njet_jet_phi->GetNbinsX() + 1);
    njet_jet_phi->Scale(njet_jet_phi_scale);
    if (B == 0) njet_jet_phi->SetLineColor(kYellow);
    if (B == 1) njet_jet_phi->SetLineColor(kGreen);
    if (B == 2) njet_jet_phi->SetLineColor(kRed);
    if (B == 3) njet_jet_phi->SetLineColor(kBlue);
    njet_jet_phi->SetTitle("phi of 1-jettiness jet");
    njet_jet_phi->Draw("SAMES");
    leg_phi_1jet->AddEntry(njet_jet_phi, "#beta = " + (TString)ss.str());
    if (njet_jet_phi->GetMaximum() > max_val_phi) max_val_phi = njet_jet_phi->GetMaximum();
  }

  akt_jet_phi->SetMaximum(1.1*max_val_phi);
  leg_phi_1jet->Draw("SAMES");
  phi_compare_1jet->Write();
  phi_compare_1jet->Print("zprime_phi_compare_1jet.eps", "eps");

  TCanvas *eta_compare_1jet = new TCanvas("eta_compare_1jet", "eta_compare_1jet", 600, 600);
  eta_compare_1jet->cd();
  double akt_jet_eta_scale = 1/akt_jet_eta->Integral(0, akt_jet_eta->GetNbinsX() + 1);
  akt_jet_eta->Scale(akt_jet_eta_scale);
  akt_jet_eta->SetLineColor(kBlack);
  akt_jet_eta->SetTitle("Rapidity of Individual Jet");
  akt_jet_eta->GetXaxis()->SetTitle("#eta_{j}");
  akt_jet_eta->GetYaxis()->SetTitle("Relative Occurrence");
  akt_jet_eta->GetYaxis()->SetTitleOffset(1.6);
  akt_jet_eta->SetMinimum(0.0);
  double max_val_eta = akt_jet_eta->GetMaximum();

  double akt_wta_jet_eta_scale = 1/akt_wta_jet_eta->Integral(0, akt_jet_eta->GetNbinsX() + 1);
  akt_wta_jet_eta->Scale(akt_wta_jet_eta_scale);
  akt_wta_jet_eta->SetLineColor(28);
  if (akt_wta_jet_eta->GetMaximum() > max_val) max_val = akt_wta_jet_eta->GetMaximum();

  akt_jet_eta->Draw();
  akt_wta_jet_eta->Draw("SAMES");
  TLegend *leg_eta_1jet = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg_eta_1jet->SetFillColor(kWhite);
  leg_eta_1jet->SetLineColor(kWhite);
  leg_eta_1jet->AddEntry(akt_jet_eta, "Anti-KT (E-scheme)", "L");
  leg_eta_1jet->AddEntry(akt_wta_jet_eta, "Anti-KT (WTA)", "L");

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* njet_jet_eta = (TH1*)njet_jet_eta_hists.At(B);
    double njet_jet_eta_scale = 1/njet_jet_eta->Integral(0, njet_jet_eta->GetNbinsX() + 1);
    njet_jet_eta->Scale(njet_jet_eta_scale);
    if (B == 0) njet_jet_eta->SetLineColor(kYellow);
    if (B == 1) njet_jet_eta->SetLineColor(kGreen);
    if (B == 2) njet_jet_eta->SetLineColor(kRed);
    if (B == 3) njet_jet_eta->SetLineColor(kBlue);
    njet_jet_eta->SetTitle("eta of 1-jettiness jet");
    njet_jet_eta->Draw("SAMES");
    leg_eta_1jet->AddEntry(njet_jet_eta, "#beta = " + (TString)ss.str());
    if (njet_jet_eta->GetMaximum() > max_val_eta) max_val_eta = njet_jet_eta->GetMaximum();
  }

  akt_jet_eta->SetMaximum(1.2*max_val_eta);
  leg_eta_1jet->Draw("SAMES");
  eta_compare_1jet->Write();
  eta_compare_1jet->Print("zprime_eta_compare_1jet.eps", "eps");

  TCanvas *massdiff_compare_1jet = new TCanvas("massdiff_compare_1jet", "massdiff_compare_1jet", 600, 600);
  massdiff_compare_1jet->cd();
  massdiff_compare_1jet->SetLogy();
  double max_val_massdiff = 0;
  TLegend *leg_massdiff_1jet = new TLegend(0.5, 0.5, 0.88, 0.88);
  leg_massdiff_1jet->SetFillColor(kWhite);
  leg_massdiff_1jet->SetLineColor(kWhite);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* jetmass_diff = (TH1*)jetmass_diff_hists.At(B);
    double jetmass_diff_scale = 1/jetmass_diff->Integral(0, jetmass_diff->GetNbinsX() + 1);
    jetmass_diff->Scale(jetmass_diff_scale);
    if (B == 0) jetmass_diff->SetLineColor(kYellow);
    if (B == 1) jetmass_diff->SetLineColor(kGreen);
    if (B == 2) jetmass_diff->SetLineColor(kRed);
    if (B == 3) jetmass_diff->SetLineColor(kBlue);
    jetmass_diff->SetTitle("Mass difference of hardest anti-kT (E-scheme) and 1-jettiness jets");
    jetmass_diff->GetXaxis()->SetTitle("m_{akt} - m_{njet}");
    jetmass_diff->GetYaxis()->SetTitle("Relative Occurrence");
    jetmass_diff->GetYaxis()->SetTitleOffset(1.5);
    jetmass_diff->Draw("SAMES");
    leg_massdiff_1jet->AddEntry(jetmass_diff, "#beta = " + (TString)ss.str());
    if (jetmass_diff->GetMaximum() > max_val_massdiff) {
      max_val_massdiff = jetmass_diff->GetMaximum();
      jetmass_diff->SetMaximum(1.2*max_val_massdiff);
    }
  }

  leg_massdiff_1jet->Draw("SAMES");
  massdiff_compare_1jet->Write();
  massdiff_compare_1jet->Print("zprime_massdiff_compare_1jet.eps", "eps");

  TCanvas *perpdiff_compare_1jet = new TCanvas("perpdiff_compare_1jet", "perpdiff_compare_1jet", 600, 600);
  perpdiff_compare_1jet->cd();
  perpdiff_compare_1jet->SetLogy();
  double max_val_perpdiff = 0;
  TLegend *leg_perpdiff_1jet = new TLegend(0.5, 0.5, 0.88, 0.88);
  leg_perpdiff_1jet->SetFillColor(kWhite);
  leg_perpdiff_1jet->SetLineColor(kWhite);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* jetperp_diff = (TH1*)jetperp_diff_hists.At(B);
    double jetperp_diff_scale = 1/jetperp_diff->Integral(0, jetperp_diff->GetNbinsX() + 1);
    jetperp_diff->Scale(jetperp_diff_scale);
    if (B == 0) jetperp_diff->SetLineColor(kYellow);
    if (B == 1) jetperp_diff->SetLineColor(kGreen);
    if (B == 2) jetperp_diff->SetLineColor(kRed);
    if (B == 3) jetperp_diff->SetLineColor(kBlue);
    jetperp_diff->SetTitle("pT difference of hardest anti-kT (E-scheme) and 1-jettiness jets");
    jetperp_diff->GetXaxis()->SetTitle("pT_{akt} - pT_{njet}");
    jetperp_diff->GetYaxis()->SetTitle("Relative Occurrence");
    jetperp_diff->GetYaxis()->SetTitleOffset(1.5);

    if (B == 0) jetperp_diff->Draw();
    else jetperp_diff->Draw("SAMES");

    leg_perpdiff_1jet->AddEntry(jetperp_diff, "#beta = " + (TString)ss.str());
    if (jetperp_diff->GetMaximum() > max_val_perpdiff) {
      max_val_perpdiff = jetperp_diff->GetMaximum();
      jetperp_diff->SetMaximum(1.2*max_val_perpdiff);
    }
  }

  leg_perpdiff_1jet->Draw("SAMES");
  perpdiff_compare_1jet->Write();
  perpdiff_compare_1jet->Print("zprime_perpdiff_compare_1jet.eps", "eps");

  TCanvas *phidiff_compare_1jet = new TCanvas("phidiff_compare_1jet", "phidiff_compare_1jet", 600, 600);
  phidiff_compare_1jet->cd();
  phidiff_compare_1jet->SetLogy();
  double max_val_phidiff = 0;
  TLegend *leg_phidiff_1jet = new TLegend(0.3, 0.5, 0.7, 0.88);
  leg_phidiff_1jet->SetFillColor(kWhite);
  leg_phidiff_1jet->SetLineColor(kWhite);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* jetphi_diff = (TH1*)jetphi_diff_hists.At(B);
    double jetphi_diff_scale = 1/jetphi_diff->Integral(0, jetphi_diff->GetNbinsX() + 1);
    jetphi_diff->Scale(jetphi_diff_scale);
    if (B == 0) jetphi_diff->SetLineColor(kYellow);
    if (B == 1) jetphi_diff->SetLineColor(kGreen);
    if (B == 2) jetphi_diff->SetLineColor(kRed);
    if (B == 3) jetphi_diff->SetLineColor(kBlue);
    jetphi_diff->SetTitle("#phi difference of hardest anti-kT (E-scheme) and 1-jettiness jets");
    jetphi_diff->GetXaxis()->SetTitle("#phi_{akt} - #phi_{njet}");
    jetphi_diff->GetYaxis()->SetTitle("Relative Occurrence");
    jetphi_diff->GetYaxis()->SetTitleOffset(1.5);

    if (B == 0) jetphi_diff->Draw();
    else jetphi_diff->Draw("SAMES");
    leg_phidiff_1jet->AddEntry(jetphi_diff, "#beta = " + (TString)ss.str());
    if (jetphi_diff->GetMaximum() > max_val_phidiff) {
      max_val_phidiff = jetphi_diff->GetMaximum();
      jetphi_diff->SetMaximum(1.2*max_val_phidiff);
    }
  }

  leg_phidiff_1jet->Draw("SAMES");
  phidiff_compare_1jet->Write();
  phidiff_compare_1jet->Print("zprime_phidiff_compare_1jet.eps", "eps");

  TCanvas *tau1diff_compare_1jet = new TCanvas("tau1diff_compare_1jet", "tau1diff_compare_1jet", 600, 600);
  tau1diff_compare_1jet->cd();
  tau1diff_compare_1jet->SetLogy();
  double max_val_tau1diff = 0;
  TLegend *leg_tau1diff_1jet = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_tau1diff_1jet->SetFillColor(kWhite);
  leg_tau1diff_1jet->SetLineColor(kWhite);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* tau1_diff = (TH1*)tau1_diff_hists.At(B);
    double tau1_diff_scale = 1/tau1_diff->Integral(0, tau1_diff->GetNbinsX() + 1);
    tau1_diff->Scale(tau1_diff_scale);
    if (B == 0) tau1_diff->SetLineColor(kYellow);
    if (B == 1) tau1_diff->SetLineColor(kGreen);
    if (B == 2) tau1_diff->SetLineColor(kRed);
    if (B == 3) tau1_diff->SetLineColor(kBlue);
    tau1_diff->SetTitle("#tau_{1} difference of hardest anti-kT (E-scheme) and 1-jettiness jets");
    tau1_diff->GetXaxis()->SetTitle("#tau_{1}_{akt} - #tau_{1}_{njet}");
    tau1_diff->GetYaxis()->SetTitle("Relative Occurrence");
    tau1_diff->GetYaxis()->SetTitleOffset(1.5);    tau1_diff->Draw("SAMES");
    leg_tau1diff_1jet->AddEntry(tau1_diff, "#beta = " + (TString)ss.str());
    if (tau1_diff->GetMaximum() > max_val_tau1diff) {
      max_val_tau1diff = tau1_diff->GetMaximum();
      tau1_diff->SetMaximum(1.2*max_val_tau1diff);
    }
  }

  leg_tau1diff_1jet->Draw("SAMES");
  tau1diff_compare_1jet->Write();
  tau1diff_compare_1jet->Print("zprime_tau1diff_compare_1jet.eps", "eps");

  TCanvas *massdiff_wta_compare_1jet = new TCanvas("massdiff_wta_compare_1jet", "massdiff_wta_compare_1jet", 600, 600);
  massdiff_wta_compare_1jet->cd();
  massdiff_wta_compare_1jet->SetLogy();
  double max_val_massdiff_wta = 0;
  TLegend *leg_massdiff_wta_1jet = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_massdiff_wta_1jet->SetFillColor(kWhite);
  leg_massdiff_wta_1jet->SetLineColor(kWhite);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* jetmass_diff_wta = (TH1*)jetmass_diff_wta_hists.At(B);
    double jetmass_diff_wta_scale = 1/jetmass_diff_wta->Integral(0, jetmass_diff_wta->GetNbinsX() + 1);
    jetmass_diff_wta->Scale(jetmass_diff_wta_scale);
    if (B == 0) jetmass_diff_wta->SetLineColor(kYellow);
    if (B == 1) jetmass_diff_wta->SetLineColor(kGreen);
    if (B == 2) jetmass_diff_wta->SetLineColor(kRed);
    if (B == 3) jetmass_diff_wta->SetLineColor(kBlue);
    jetmass_diff_wta->SetTitle("Mass difference of hardest anti-kT (WTA) and 1-jettiness jets");
    jetmass_diff_wta->GetXaxis()->SetTitle("m_{akt} - m_{njet}");
    jetmass_diff_wta->GetYaxis()->SetTitle("Relative Occurrence");
    jetmass_diff_wta->GetYaxis()->SetTitleOffset(1.5);
    jetmass_diff_wta->Draw("SAMES");
    leg_massdiff_wta_1jet->AddEntry(jetmass_diff_wta, "#beta = " + (TString)ss.str());
    if (jetmass_diff_wta->GetMaximum() > max_val_massdiff_wta) {
      max_val_massdiff_wta = jetmass_diff_wta->GetMaximum();
      jetmass_diff_wta->SetMaximum(1.2*max_val_massdiff_wta);
    }
  }

  leg_massdiff_wta_1jet->Draw("SAMES");
  massdiff_wta_compare_1jet->Write();
  massdiff_wta_compare_1jet->Print("zprime_massdiff_wta_compare_1jet.eps", "eps");

  TCanvas *perpdiff_wta_compare_1jet = new TCanvas("perpdiff_wta_compare_1jet", "perpdiff_wta_compare_1jet", 600, 600);
  perpdiff_wta_compare_1jet->cd();
  perpdiff_wta_compare_1jet->SetLogy();
  double max_val_perpdiff_wta = 0;
  TLegend *leg_perpdiff_wta_1jet = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_perpdiff_wta_1jet->SetFillColor(kWhite);
  leg_perpdiff_wta_1jet->SetLineColor(kWhite);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* jetperp_diff_wta = (TH1*)jetperp_diff_wta_hists.At(B);
    double jetperp_diff_wta_scale = 1/jetperp_diff_wta->Integral(0, jetperp_diff_wta->GetNbinsX() + 1);
    jetperp_diff_wta->Scale(jetperp_diff_wta_scale);
    if (B == 0) jetperp_diff_wta->SetLineColor(kYellow);
    if (B == 1) jetperp_diff_wta->SetLineColor(kGreen);
    if (B == 2) jetperp_diff_wta->SetLineColor(kRed);
    if (B == 3) jetperp_diff_wta->SetLineColor(kBlue);
    jetperp_diff_wta->SetTitle("pT difference of hardest anti-kT (WTA) and 1-jettiness jets");
    jetperp_diff_wta->GetXaxis()->SetTitle("pT_{akt} - pT_{njet}");
    jetperp_diff_wta->GetYaxis()->SetTitle("Relative Occurrence");
    jetperp_diff_wta->GetYaxis()->SetTitleOffset(1.5);
    jetperp_diff_wta->Draw("SAMES");
    leg_perpdiff_wta_1jet->AddEntry(jetperp_diff_wta, "#beta = " + (TString)ss.str());
    if (jetperp_diff_wta->GetMaximum() > max_val_perpdiff_wta) {
      max_val_perpdiff_wta = jetperp_diff_wta->GetMaximum();
      jetperp_diff_wta->SetMaximum(1.2*max_val_perpdiff_wta);
    }
  }

  leg_perpdiff_wta_1jet->Draw("SAMES");
  perpdiff_wta_compare_1jet->Write();
  perpdiff_wta_compare_1jet->Print("zprime_perpdiff_wta_compare_1jet.eps", "eps");

  TCanvas *phidiff_wta_compare_1jet = new TCanvas("phidiff_wta_compare_1jet", "phidiff_wta_compare_1jet", 600, 600);
  phidiff_wta_compare_1jet->cd();
  phidiff_wta_compare_1jet->SetLogy();
  double max_val_phidiff_wta = 0;
  TLegend *leg_phidiff_wta_1jet = new TLegend(0.7, 0.7, 0.9, 0.9);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* jetphi_diff_wta = (TH1*)jetphi_diff_wta_hists.At(B);
    double jetphi_diff_wta_scale = 1/jetphi_diff_wta->Integral(0, jetphi_diff_wta->GetNbinsX() + 1);
    jetphi_diff_wta->Scale(jetphi_diff_wta_scale);
    if (B == 0) jetphi_diff_wta->SetLineColor(kYellow);
    if (B == 1) jetphi_diff_wta->SetLineColor(kGreen);
    if (B == 2) jetphi_diff_wta->SetLineColor(kRed);
    if (B == 3) jetphi_diff_wta->SetLineColor(kBlue);
    jetphi_diff_wta->SetTitle("#phi difference of hardest anti-kT (WTA) and 1-jettiness jets");
    jetphi_diff_wta->GetXaxis()->SetTitle("#phi_{akt} - #phi_{njet}");
    jetphi_diff_wta->GetYaxis()->SetTitle("Relative Occurrence");
    jetphi_diff_wta->GetYaxis()->SetTitleOffset(1.5);
    jetphi_diff_wta->Draw("SAMES");
    leg_phidiff_wta_1jet->AddEntry(jetphi_diff_wta, "#beta = " + (TString)ss.str());
    if (jetphi_diff_wta->GetMaximum() > max_val_phidiff_wta) {
      max_val_phidiff_wta = jetphi_diff_wta->GetMaximum();
      jetphi_diff_wta->SetMaximum(1.2*max_val_phidiff_wta);
    }
  }

  leg_phidiff_wta_1jet->Draw("SAMES");
  phidiff_wta_compare_1jet->Write();
  phidiff_wta_compare_1jet->Print("zprime_phidiff_wta_compare_1jet.eps", "eps");

  TCanvas *tau1diff_wta_compare_1jet = new TCanvas("tau1diff_wta_compare_1jet", "tau1diff_wta_compare_1jet", 600, 600);
  tau1diff_wta_compare_1jet->cd();
  tau1diff_wta_compare_1jet->SetLogy();
  double max_val_tau1diff_wta = 0;
  TLegend *leg_tau1diff_wta_1jet = new TLegend(0.7, 0.7, 0.9, 0.9);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* tau1_diff_wta = (TH1*)tau1_diff_wta_hists.At(B);
    double tau1_diff_wta_scale = 1/tau1_diff_wta->Integral(0, tau1_diff_wta->GetNbinsX() + 1);
    tau1_diff_wta->Scale(tau1_diff_wta_scale);
    if (B == 0) tau1_diff_wta->SetLineColor(kYellow);
    if (B == 1) tau1_diff_wta->SetLineColor(kGreen);
    if (B == 2) tau1_diff_wta->SetLineColor(kRed);
    if (B == 3) tau1_diff_wta->SetLineColor(kBlue);
    tau1_diff_wta->SetTitle("#tau_{1} difference of hardest anti-kT (WTA) and 1-jettiness jets");
    tau1_diff_wta->GetXaxis()->SetTitle("#tau_{1}_{akt} - #tau_{1}_{njet}");
    tau1_diff_wta->GetYaxis()->SetTitle("Relative Occurrence");
    tau1_diff_wta->GetYaxis()->SetTitleOffset(1.5);    tau1_diff_wta->Draw("SAMES");
    leg_tau1diff_wta_1jet->AddEntry(tau1_diff_wta, "#beta = " + (TString)ss.str());
    if (tau1_diff_wta->GetMaximum() > max_val_tau1diff_wta) {
      max_val_tau1diff_wta = tau1_diff_wta->GetMaximum();
      tau1_diff_wta->SetMaximum(1.2*max_val_tau1diff_wta);
    }
  }

  leg_tau1diff_wta_1jet->Draw("SAMES");
  tau1diff_wta_compare_1jet->Write();
  tau1diff_wta_compare_1jet->Print("zprime_tau1diff_wta_compare_1jet.eps", "eps");


  TCanvas *invmass_compare_3jets = new TCanvas("invmass_compare_3jets", "invmass_compare_3jets", 600, 600);
  invmass_compare_3jets->cd();
  double threejet_invmass_akt_scale = 1/threejet_invmass_akt->Integral(0, threejet_invmass_akt->GetNbinsX() + 1);
  threejet_invmass_akt->Scale(threejet_invmass_akt_scale);
  threejet_invmass_akt->SetLineColor(kBlack);
  threejet_invmass_akt->SetTitle("3-jet invariant mass");
  threejet_invmass_akt->GetXaxis()->SetTitle("m_{jj} (GeV)");
  threejet_invmass_akt->GetYaxis()->SetTitle("Relative Occurrence");
  threejet_invmass_akt->GetYaxis()->SetTitleOffset(1.5);
  threejet_invmass_akt->SetMinimum(0.0);
  double max_val_invmass_3jets = threejet_invmass_akt->GetMaximum();

  double threejet_invmass_akt_wta_scale = 1/threejet_invmass_akt_wta->Integral(0, threejet_invmass_akt_wta->GetNbinsX() + 1);
  threejet_invmass_akt_wta->Scale(threejet_invmass_akt_wta_scale);
  threejet_invmass_akt_wta->SetLineColor(28);
  if (threejet_invmass_akt_wta->GetMaximum() > max_val_invmass_3jets) max_val_invmass_3jets = threejet_invmass_akt_wta->GetMaximum();

  threejet_invmass_akt->Draw();
  threejet_invmass_akt_wta->Draw("SAMES");
  TLegend *leg_invmass_3jets = new TLegend(0.12, 0.5, 0.5, 0.88);
  leg_invmass_3jets->SetFillColor(kWhite);
  leg_invmass_3jets->SetLineColor(kWhite);

  leg_invmass_3jets->AddEntry(threejet_invmass_akt, "Anti-KT (E-scheme)", "L");
  leg_invmass_3jets->AddEntry(threejet_invmass_akt_wta, "Anti-KT (WTA)", "L");

  for (int B = 0; B < n_betas; B++) {

    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* threejet_invmass_njet = (TH1*)threejet_invmass_njet_hists.At(B);
    double threejet_invmass_njet_scale = 1/threejet_invmass_njet->Integral(0, threejet_invmass_njet->GetNbinsX() + 1);
    threejet_invmass_njet->Scale(threejet_invmass_akt_scale);
    if (B == 0) threejet_invmass_njet->SetLineColor(kYellow);
    if (B == 1) threejet_invmass_njet->SetLineColor(kGreen);
    if (B == 2) threejet_invmass_njet->SetLineColor(kRed);
    if (B == 3) threejet_invmass_njet->SetLineColor(kBlue);
    threejet_invmass_njet->Draw("SAMES");
    leg_invmass_3jets->AddEntry(threejet_invmass_njet, "#beta = " + (TString)ss.str());
    if (threejet_invmass_njet->GetMaximum() > max_val_invmass_3jets) max_val_invmass_3jets = threejet_invmass_njet->GetMaximum();
  }

  threejet_invmass_akt->SetMaximum(1.2*max_val_invmass_3jets);
  leg_invmass_3jets->Draw("SAMES");
  invmass_compare_3jets->Write();
  invmass_compare_3jets->Print("zprime_invmass_compare_3jets.eps", "eps");

  TCanvas *mass_compare_3jets = new TCanvas("mass_compare_3jets", "mass_compare_3jets", 600, 600);
  mass_compare_3jets->cd();

  double threejet_thirdjet_mass_akt_scale = 1/threejet_thirdjet_mass_akt->Integral(0, threejet_thirdjet_mass_akt->GetNbinsX() + 1);
  threejet_thirdjet_mass_akt->Scale(threejet_thirdjet_mass_akt_scale);
  threejet_thirdjet_mass_akt->SetLineColor(kBlack);
  threejet_thirdjet_mass_akt->SetTitle("Mass of Third Jet");
  threejet_thirdjet_mass_akt->GetXaxis()->SetTitle("m_{j} (GeV)");
  threejet_thirdjet_mass_akt->GetYaxis()->SetTitle("Relative Occurrence");
  threejet_thirdjet_mass_akt->GetYaxis()->SetTitleOffset(1.5);
  threejet_thirdjet_mass_akt->SetMinimum(0.0);
  double max_val_mass_3jets = threejet_thirdjet_mass_akt->GetMaximum();

  double threejet_thirdjet_mass_akt_wta_scale = 1/threejet_thirdjet_mass_akt_wta->Integral(0, threejet_thirdjet_mass_akt->GetNbinsX() + 1);
  threejet_thirdjet_mass_akt_wta->Scale(threejet_thirdjet_mass_akt_wta_scale);
  threejet_thirdjet_mass_akt_wta->SetLineColor(28);
  if (threejet_thirdjet_mass_akt_wta->GetMaximum() > max_val) max_val = threejet_thirdjet_mass_akt_wta->GetMaximum();

  threejet_thirdjet_mass_akt->Draw();
  threejet_thirdjet_mass_akt_wta->Draw("SAMES");
  TLegend *leg_mass_3jets = new TLegend(0.5, 0.5, 0.88, 0.88);
  leg_mass_3jets->SetFillColor(kWhite);
  leg_mass_3jets->SetLineColor(kWhite);
  leg_mass_3jets->AddEntry(threejet_thirdjet_mass_akt, "Anti-KT (E-scheme)", "L");
  leg_mass_3jets->AddEntry(threejet_thirdjet_mass_akt_wta, "Anti-KT (WTA)", "L");

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* threejet_thirdjet_mass_njet = (TH1*)threejet_thirdjet_mass_njet_hists.At(B);
    double threejet_thirdjet_mass_njet_scale = 1/threejet_thirdjet_mass_njet->Integral(0, threejet_thirdjet_mass_njet->GetNbinsX() + 1);
    threejet_thirdjet_mass_njet->Scale(threejet_thirdjet_mass_njet_scale);
    if (B == 0) threejet_thirdjet_mass_njet->SetLineColor(kYellow);
    if (B == 1) threejet_thirdjet_mass_njet->SetLineColor(kGreen);
    if (B == 2) threejet_thirdjet_mass_njet->SetLineColor(kRed);
    if (B == 3) threejet_thirdjet_mass_njet->SetLineColor(kBlue);
    threejet_thirdjet_mass_njet->Draw("SAMES");
    leg_mass_3jets->AddEntry(threejet_thirdjet_mass_njet, "#beta = " + (TString)ss.str());
    if (threejet_thirdjet_mass_njet->GetMaximum() > max_val_mass_3jets) max_val_mass_3jets = threejet_thirdjet_mass_njet->GetMaximum();
  }

  threejet_thirdjet_mass_akt->SetMaximum(1.2*max_val_mass_3jets);
  leg_mass_3jets->Draw("SAMES");
  mass_compare_3jets->Write();
  mass_compare_3jets->Print("zprime_mass_compare_3jets.eps", "eps");

  TCanvas *perp_compare_3jets = new TCanvas("perp_compare_3jets", "perp_compare_3jets", 600, 600);
  perp_compare_3jets->cd();

  double threejet_thirdjet_perp_akt_scale = 1/threejet_thirdjet_perp_akt->Integral(0, threejet_thirdjet_perp_akt->GetNbinsX() + 1);
  threejet_thirdjet_perp_akt->Scale(threejet_thirdjet_perp_akt_scale);
  threejet_thirdjet_perp_akt->SetLineColor(kBlack);
  threejet_thirdjet_perp_akt->SetTitle("pT of Third Jet");
  threejet_thirdjet_perp_akt->GetXaxis()->SetTitle("pT_{j} (GeV)");
  threejet_thirdjet_perp_akt->GetYaxis()->SetTitle("Relative Occurrence");
  threejet_thirdjet_perp_akt->GetYaxis()->SetTitleOffset(1.5);
  threejet_thirdjet_perp_akt->SetMinimum(0.0);
  double max_val_perp_3jets = threejet_thirdjet_perp_akt->GetMaximum();

  double threejet_thirdjet_perp_akt_wta_scale = 1/threejet_thirdjet_perp_akt_wta->Integral(0, threejet_thirdjet_perp_akt->GetNbinsX() + 1);
  threejet_thirdjet_perp_akt_wta->Scale(threejet_thirdjet_perp_akt_wta_scale);
  threejet_thirdjet_perp_akt_wta->SetLineColor(28);
  if (threejet_thirdjet_perp_akt_wta->GetMaximum() > max_val) max_val = threejet_thirdjet_perp_akt_wta->GetMaximum();

  threejet_thirdjet_perp_akt->Draw();
  threejet_thirdjet_perp_akt_wta->Draw("SAMES");
  TLegend *leg_perp_3jets = new TLegend(0.5, 0.5, 0.88, 0.88);
  leg_perp_3jets->SetFillColor(kWhite);
  leg_perp_3jets->SetLineColor(kWhite);
  leg_perp_3jets->AddEntry(threejet_thirdjet_perp_akt, "Anti-KT (E-scheme)", "L");
  leg_perp_3jets->AddEntry(threejet_thirdjet_perp_akt_wta, "Anti-KT (WTA)", "L");

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* threejet_thirdjet_perp_njet = (TH1*)threejet_thirdjet_perp_njet_hists.At(B);
    double threejet_thirdjet_perp_njet_scale = 1/threejet_thirdjet_perp_njet->Integral(0, threejet_thirdjet_perp_njet->GetNbinsX() + 1);
    threejet_thirdjet_perp_njet->Scale(threejet_thirdjet_perp_njet_scale);
    if (B == 0) threejet_thirdjet_perp_njet->SetLineColor(kYellow);
    if (B == 1) threejet_thirdjet_perp_njet->SetLineColor(kGreen);
    if (B == 2) threejet_thirdjet_perp_njet->SetLineColor(kRed);
    if (B == 3) threejet_thirdjet_perp_njet->SetLineColor(kBlue);
    threejet_thirdjet_perp_njet->Draw("SAMES");
    leg_perp_3jets->AddEntry(threejet_thirdjet_perp_njet, "#beta = " + (TString)ss.str());
    if (threejet_thirdjet_perp_njet->GetMaximum() > max_val_perp_3jets) max_val_perp_3jets = threejet_thirdjet_perp_njet->GetMaximum();
  }

  threejet_thirdjet_perp_akt->SetMaximum(1.2*max_val_perp_3jets);
  leg_perp_3jets->Draw("SAMES");
  perp_compare_3jets->Write();
  perp_compare_3jets->Print("zprime_perp_compare_3jets.eps", "eps");

  TCanvas *phi_compare_3jets = new TCanvas("phi_compare_3jets", "phi_compare_3jets", 600, 600);
  phi_compare_3jets->cd();

  double threejet_thirdjet_phi_akt_scale = 1/threejet_thirdjet_phi_akt->Integral(0, threejet_thirdjet_phi_akt->GetNbinsX() + 1);
  threejet_thirdjet_phi_akt->Scale(threejet_thirdjet_phi_akt_scale);
  threejet_thirdjet_phi_akt->SetLineColor(kBlack);
  threejet_thirdjet_phi_akt->SetTitle("Azimuth of Third Jet");
  threejet_thirdjet_phi_akt->GetXaxis()->SetTitle("#phi_{j} (GeV)");
  threejet_thirdjet_phi_akt->GetYaxis()->SetTitle("Relative Occurrence");
  threejet_thirdjet_phi_akt->GetYaxis()->SetTitleOffset(1.5);
  threejet_thirdjet_phi_akt->SetMinimum(0.0);
  double max_val_phi_3jets = threejet_thirdjet_phi_akt->GetMaximum();

  double threejet_thirdjet_phi_akt_wta_scale = 1/threejet_thirdjet_phi_akt_wta->Integral(0, threejet_thirdjet_phi_akt->GetNbinsX() + 1);
  threejet_thirdjet_phi_akt_wta->Scale(threejet_thirdjet_phi_akt_wta_scale);
  threejet_thirdjet_phi_akt_wta->SetLineColor(28);
  if (threejet_thirdjet_phi_akt_wta->GetMaximum() > max_val) max_val = threejet_thirdjet_phi_akt_wta->GetMaximum();

  threejet_thirdjet_phi_akt->Draw();
  threejet_thirdjet_phi_akt_wta->Draw("SAMES");
  TLegend *leg_phi_3jets = new TLegend(0.5, 0.5, 0.8, 0.8);
  leg_phi_3jets->SetFillColor(kWhite);
  leg_phi_3jets->SetLineColor(kWhite);
  leg_phi_3jets->AddEntry(threejet_thirdjet_phi_akt, "Anti-KT (E-scheme)", "L");
  leg_phi_3jets->AddEntry(threejet_thirdjet_phi_akt_wta, "Anti-KT (WTA)", "L");

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* threejet_thirdjet_phi_njet = (TH1*)threejet_thirdjet_phi_njet_hists.At(B);
    double threejet_thirdjet_phi_njet_scale = 1/threejet_thirdjet_phi_njet->Integral(0, threejet_thirdjet_phi_njet->GetNbinsX() + 1);
    threejet_thirdjet_phi_njet->Scale(threejet_thirdjet_phi_njet_scale);
    if (B == 0) threejet_thirdjet_phi_njet->SetLineColor(kYellow);
    if (B == 1) threejet_thirdjet_phi_njet->SetLineColor(kGreen);
    if (B == 2) threejet_thirdjet_phi_njet->SetLineColor(kRed);
    if (B == 3) threejet_thirdjet_phi_njet->SetLineColor(kBlue);
    threejet_thirdjet_phi_njet->Draw("SAMES");
    leg_phi_3jets->AddEntry(threejet_thirdjet_phi_njet, "#beta = " + (TString)ss.str());
    if (threejet_thirdjet_phi_njet->GetMaximum() > max_val_phi_3jets) max_val_phi_3jets = threejet_thirdjet_phi_njet->GetMaximum();
  }

  threejet_thirdjet_phi_akt->SetMaximum(1.2*max_val_phi_3jets);
  leg_phi_3jets->Draw("SAMES");
  phi_compare_3jets->Write();
  phi_compare_3jets->Print("zprime_phi_compare_3jets.eps", "eps");

  TCanvas *massdiff_compare_3jets = new TCanvas("massdiff_compare_3jets", "massdiff_compare_3jets", 600, 600);
  massdiff_compare_3jets->cd();
  massdiff_compare_3jets->SetLogy();
  double max_val_massdiff_3jets = 0;
  TLegend *leg_massdiff_3jets = new TLegend(0.5, 0.5, 0.88, 0.88);
  leg_massdiff_3jets->SetFillColor(kWhite);
  leg_massdiff_3jets->SetLineColor(kWhite);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* threejet_thirdjet_mass_diff = (TH1*)threejet_thirdjet_mass_diff_hists.At(B);
    double threejet_thirdjet_mass_diff_scale = 1/threejet_thirdjet_mass_diff->Integral(0, threejet_thirdjet_mass_diff->GetNbinsX() + 1);
    threejet_thirdjet_mass_diff->Scale(threejet_thirdjet_mass_diff_scale);
    if (B == 0) threejet_thirdjet_mass_diff->SetLineColor(kYellow);
    if (B == 1) threejet_thirdjet_mass_diff->SetLineColor(kGreen);
    if (B == 2) threejet_thirdjet_mass_diff->SetLineColor(kRed);
    if (B == 3) threejet_thirdjet_mass_diff->SetLineColor(kBlue);
    threejet_thirdjet_mass_diff->SetTitle("Mass difference of third hardest anti-kT (E-scheme) and 3-jettiness jet");
    threejet_thirdjet_mass_diff->GetXaxis()->SetTitle("m_{akt} - m_{njet}");
    threejet_thirdjet_mass_diff->GetYaxis()->SetTitle("Relative Occurrence");
    threejet_thirdjet_mass_diff->GetYaxis()->SetTitleOffset(1.5);
    
    if (B == 0) threejet_thirdjet_mass_diff->Draw();
    else threejet_thirdjet_mass_diff->Draw("SAMES");
    leg_massdiff_3jets->AddEntry(threejet_thirdjet_mass_diff, "#beta = " + (TString)ss.str());
    if (threejet_thirdjet_mass_diff->GetMaximum() > max_val_massdiff_3jets) {
      max_val_massdiff_3jets = threejet_thirdjet_mass_diff->GetMaximum();
      threejet_thirdjet_mass_diff->SetMaximum(1.2*max_val_massdiff_3jets);
    }
  }

  leg_massdiff_3jets->Draw("SAMES");
  massdiff_compare_3jets->Write();
  massdiff_compare_3jets->Print("zprime_massdiff_compare_3jets.eps", "eps");

  TCanvas *perpdiff_compare_3jets = new TCanvas("perpdiff_compare_3jets", "perpdiff_compare_3jets", 600, 600);
  perpdiff_compare_3jets->cd();
  perpdiff_compare_3jets->SetLogy();
  double max_val_perpdiff_3jets = 0;
  TLegend *leg_perpdiff_3jets = new TLegend(0.6, 0.6, 0.88, 0.88);
  leg_perpdiff_3jets->SetFillColor(kWhite);
  leg_perpdiff_3jets->SetLineColor(kWhite);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* threejet_thirdjet_perp_diff = (TH1*)threejet_thirdjet_perp_diff_hists.At(B);
    double threejet_thirdjet_perp_diff_scale = 1/threejet_thirdjet_perp_diff->Integral(0, threejet_thirdjet_perp_diff->GetNbinsX() + 1);
    threejet_thirdjet_perp_diff->Scale(threejet_thirdjet_perp_diff_scale);
    if (B == 0) threejet_thirdjet_perp_diff->SetLineColor(kYellow);
    if (B == 1) threejet_thirdjet_perp_diff->SetLineColor(kGreen);
    if (B == 2) threejet_thirdjet_perp_diff->SetLineColor(kRed);
    if (B == 3) threejet_thirdjet_perp_diff->SetLineColor(kBlue);
    threejet_thirdjet_perp_diff->SetTitle("pT difference of third hardest anti-kT (E-scheme) and 3-jettiness jet");
    threejet_thirdjet_perp_diff->GetXaxis()->SetTitle("pT_{akt} - pT_{njet}");
    threejet_thirdjet_perp_diff->GetYaxis()->SetTitle("Relative Occurrence");
    threejet_thirdjet_perp_diff->GetYaxis()->SetTitleOffset(1.5);
    if (B == 0) threejet_thirdjet_perp_diff->Draw();
    else threejet_thirdjet_perp_diff->Draw("SAMES");
    leg_perpdiff_3jets->AddEntry(threejet_thirdjet_perp_diff, "#beta = " + (TString)ss.str());
    if (threejet_thirdjet_perp_diff->GetMaximum() > max_val_perpdiff_3jets) {
      max_val_perpdiff_3jets = threejet_thirdjet_perp_diff->GetMaximum();
      threejet_thirdjet_perp_diff->SetMaximum(1.2*max_val_perpdiff_3jets);
    }
  }

  leg_perpdiff_3jets->Draw("SAMES");
  perpdiff_compare_3jets->Write();
  perpdiff_compare_3jets->Print("zprime_perpdiff_compare_3jets.eps", "eps");

  TCanvas *phidiff_compare_3jets = new TCanvas("phidiff_compare_3jets", "phidiff_compare_3jets", 600, 600);
  phidiff_compare_3jets->cd();
  phidiff_compare_3jets->SetLogy();
  double max_val_phidiff_3jets = 0;
  TLegend *leg_phidiff_3jets = new TLegend(0.5, 0.5, 0.88, 0.88);
  leg_phidiff_3jets->SetFillColor(kWhite);
  leg_phidiff_3jets->SetLineColor(kWhite);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* threejet_thirdjet_phi_diff = (TH1*)threejet_thirdjet_phi_diff_hists.At(B);
    double threejet_thirdjet_phi_diff_scale = 1/threejet_thirdjet_phi_diff->Integral(0, threejet_thirdjet_phi_diff->GetNbinsX() + 1);
    threejet_thirdjet_phi_diff->Scale(threejet_thirdjet_phi_diff_scale);
    if (B == 0) threejet_thirdjet_phi_diff->SetLineColor(kYellow);
    if (B == 1) threejet_thirdjet_phi_diff->SetLineColor(kGreen);
    if (B == 2) threejet_thirdjet_phi_diff->SetLineColor(kRed);
    if (B == 3) threejet_thirdjet_phi_diff->SetLineColor(kBlue);
    threejet_thirdjet_phi_diff->SetTitle("Azimuthal difference of third hardest anti-kT (E-scheme) and 3-jettiness jet");
    threejet_thirdjet_phi_diff->GetXaxis()->SetTitle("|#phi_{akt} - #phi_{njet}|");
    threejet_thirdjet_phi_diff->GetYaxis()->SetTitle("Relative Occurrence");
    threejet_thirdjet_phi_diff->GetYaxis()->SetTitleOffset(1.5);

    if (B == 0) threejet_thirdjet_phi_diff->Draw();
    else threejet_thirdjet_phi_diff->Draw("SAMES");
    leg_phidiff_3jets->AddEntry(threejet_thirdjet_phi_diff, "#beta = " + (TString)ss.str());
    if (threejet_thirdjet_phi_diff->GetMaximum() > max_val_phidiff_3jets) {
      max_val_phidiff_3jets = threejet_thirdjet_phi_diff->GetMaximum();
      threejet_thirdjet_phi_diff->SetMaximum(1.2*max_val_phidiff_3jets);
    }
  }

  leg_phidiff_3jets->Draw("SAMES");
  phidiff_compare_3jets->Write();
  phidiff_compare_3jets->Print("zprime_phidiff_compare_3jets.eps", "eps");

  TCanvas *tau3diff_compare_3jets = new TCanvas("tau3diff_compare_3jets", "tau3diff_compare_3jets", 600, 600);
  tau3diff_compare_3jets->cd();
  tau3diff_compare_3jets->SetLogy();
  double max_val_tau3diff = 0;
  TLegend *leg_tau3diff_3jets = new TLegend(0.5, 0.5, 0.88, 0.88);
  leg_tau3diff_3jets->SetFillColor(kWhite);
  leg_tau3diff_3jets->SetLineColor(kWhite);


  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* tau3_diff = (TH1*)tau3_diff_hists.At(B);
    double tau3_diff_scale = 1/tau3_diff->Integral(0, tau3_diff->GetNbinsX() + 1);
    tau3_diff->Scale(tau3_diff_scale);
    if (B == 0) tau3_diff->SetLineColor(kYellow);
    if (B == 1) tau3_diff->SetLineColor(kGreen);
    if (B == 2) tau3_diff->SetLineColor(kRed);
    if (B == 3) tau3_diff->SetLineColor(kBlue);
    tau3_diff->SetTitle("#tau_{3} difference of 3 hardest anti-kT (E-scheme) and 3-jettiness jets");
    tau3_diff->GetXaxis()->SetTitle("#tau_{3}_{akt} - #tau_{3}_{njet}");
    tau3_diff->GetYaxis()->SetTitle("Relative Occurrence");
    tau3_diff->GetYaxis()->SetTitleOffset(1.5);    tau3_diff->Draw("SAMES");
    leg_tau3diff_3jets->AddEntry(tau3_diff, "#beta = " + (TString)ss.str());
    if (tau3_diff->GetMaximum() > max_val_tau3diff) {
      max_val_tau3diff = tau3_diff->GetMaximum();
      tau3_diff->SetMaximum(1.2*max_val_tau3diff);
    }
  }

  leg_tau3diff_3jets->Draw("SAMES");
  tau3diff_compare_3jets->Write();
  tau3diff_compare_3jets->Print("zprime_tau3diff_compare_3jets.eps", "eps");


  TCanvas *massdiff_wta_compare_3jets = new TCanvas("massdiff_wta_compare_3jets", "massdiff_wta_compare_3jets", 600, 600);
  massdiff_wta_compare_3jets->cd();
  massdiff_wta_compare_3jets->SetLogy();
  double max_val_massdiff_wta_3jets = 0;
  TLegend *leg_massdiff_wta_3jets = new TLegend(0.7, 0.7, 0.9, 0.9);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* threejet_thirdjet_mass_diff_wta = (TH1*)threejet_thirdjet_mass_diff_wta_hists.At(B);
    double threejet_thirdjet_mass_diff_wta_scale = 1/threejet_thirdjet_mass_diff_wta->Integral(0, threejet_thirdjet_mass_diff_wta->GetNbinsX() + 1);
    threejet_thirdjet_mass_diff_wta->Scale(threejet_thirdjet_mass_diff_wta_scale);
    if (B == 0) threejet_thirdjet_mass_diff_wta->SetLineColor(kYellow);
    if (B == 1) threejet_thirdjet_mass_diff_wta->SetLineColor(kGreen);
    if (B == 2) threejet_thirdjet_mass_diff_wta->SetLineColor(kRed);
    if (B == 3) threejet_thirdjet_mass_diff_wta->SetLineColor(kBlue);
    threejet_thirdjet_mass_diff_wta->SetTitle("Mass difference of third hardest anti-kT (WTA) and 3-jettiness jet");
    threejet_thirdjet_mass_diff_wta->GetXaxis()->SetTitle("m_{akt} - m_{njet}");
    threejet_thirdjet_mass_diff_wta->GetYaxis()->SetTitle("Relative Occurrence");
    threejet_thirdjet_mass_diff_wta->GetYaxis()->SetTitleOffset(1.5);
    threejet_thirdjet_mass_diff_wta->Draw("SAMES");
    leg_massdiff_wta_3jets->AddEntry(threejet_thirdjet_mass_diff_wta, "#beta = " + (TString)ss.str());
    if (threejet_thirdjet_mass_diff_wta->GetMaximum() > max_val_massdiff_wta_3jets) {
      max_val_massdiff_wta_3jets = threejet_thirdjet_mass_diff_wta->GetMaximum();
      threejet_thirdjet_mass_diff_wta->SetMaximum(1.2*max_val_massdiff_wta_3jets);
    }
  }

  leg_massdiff_wta_3jets->Draw("SAMES");
  massdiff_wta_compare_3jets->Write();
  massdiff_wta_compare_3jets->Print("zprime_massdiff_wta_compare_3jets.eps", "eps");

  TCanvas *perpdiff_wta_compare_3jets = new TCanvas("perpdiff_wta_compare_3jets", "perpdiff_wta_compare_3jets", 600, 600);
  perpdiff_wta_compare_3jets->cd();
  perpdiff_wta_compare_3jets->SetLogy();
  double max_val_perpdiff_wta_3jets = 0;
  TLegend *leg_perpdiff_wta_3jets = new TLegend(0.7, 0.7, 0.9, 0.9);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* threejet_thirdjet_perp_diff_wta = (TH1*)threejet_thirdjet_perp_diff_wta_hists.At(B);
    double threejet_thirdjet_perp_diff_wta_scale = 1/threejet_thirdjet_perp_diff_wta->Integral(0, threejet_thirdjet_perp_diff_wta->GetNbinsX() + 1);
    threejet_thirdjet_perp_diff_wta->Scale(threejet_thirdjet_perp_diff_wta_scale);
    if (B == 0) threejet_thirdjet_perp_diff_wta->SetLineColor(kYellow);
    if (B == 1) threejet_thirdjet_perp_diff_wta->SetLineColor(kGreen);
    if (B == 2) threejet_thirdjet_perp_diff_wta->SetLineColor(kRed);
    if (B == 3) threejet_thirdjet_perp_diff_wta->SetLineColor(kBlue);
    threejet_thirdjet_perp_diff_wta->SetTitle("pT difference of third hardest anti-kT (WTA) and 3-jettiness jet");
    threejet_thirdjet_perp_diff_wta->GetXaxis()->SetTitle("pT_{akt} - pT_{njet}");
    threejet_thirdjet_perp_diff_wta->GetYaxis()->SetTitle("Relative Occurrence");
    threejet_thirdjet_perp_diff_wta->GetYaxis()->SetTitleOffset(1.5);
    threejet_thirdjet_perp_diff_wta->Draw("SAMES");
    leg_perpdiff_wta_3jets->AddEntry(threejet_thirdjet_perp_diff_wta, "#beta = " + (TString)ss.str());
    if (threejet_thirdjet_perp_diff_wta->GetMaximum() > max_val_perpdiff_wta_3jets) {
      max_val_perpdiff_wta_3jets = threejet_thirdjet_perp_diff_wta->GetMaximum();
      threejet_thirdjet_perp_diff_wta->SetMaximum(1.2*max_val_perpdiff_wta_3jets);
    }
  }

  leg_perpdiff_wta_3jets->Draw("SAMES");
  perpdiff_wta_compare_3jets->Write();
  perpdiff_wta_compare_3jets->Print("zprime_perpdiff_wta_compare_3jets.eps", "eps");

  TCanvas *phidiff_wta_compare_3jets = new TCanvas("phidiff_wta_compare_3jets", "phidiff_wta_compare_3jets", 600, 600);
  phidiff_wta_compare_3jets->cd();
  phidiff_wta_compare_3jets->SetLogy();
  double max_val_phidiff_wta_3jets = 0;
  TLegend *leg_phidiff_wta_3jets = new TLegend(0.7, 0.7, 0.9, 0.9);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* threejet_thirdjet_phi_diff_wta = (TH1*)threejet_thirdjet_phi_diff_wta_hists.At(B);
    double threejet_thirdjet_phi_diff_wta_scale = 1/threejet_thirdjet_phi_diff_wta->Integral(0, threejet_thirdjet_phi_diff_wta->GetNbinsX() + 1);
    threejet_thirdjet_phi_diff_wta->Scale(threejet_thirdjet_phi_diff_wta_scale);
    if (B == 0) threejet_thirdjet_phi_diff_wta->SetLineColor(kYellow);
    if (B == 1) threejet_thirdjet_phi_diff_wta->SetLineColor(kGreen);
    if (B == 2) threejet_thirdjet_phi_diff_wta->SetLineColor(kRed);
    if (B == 3) threejet_thirdjet_phi_diff_wta->SetLineColor(kBlue);
    threejet_thirdjet_phi_diff_wta->SetTitle("Azimuthal difference of third hardest anti-kT (E-scheme) and 3-jettiness jet");
    threejet_thirdjet_phi_diff_wta->GetXaxis()->SetTitle("|#phi_{akt} - #phi_{njet}|");
    threejet_thirdjet_phi_diff_wta->GetYaxis()->SetTitle("Relative Occurrence");
    threejet_thirdjet_phi_diff_wta->GetYaxis()->SetTitleOffset(1.5);
    threejet_thirdjet_phi_diff_wta->Draw("SAMES");
    leg_phidiff_wta_3jets->AddEntry(threejet_thirdjet_phi_diff_wta, "#beta = " + (TString)ss.str());
    if (threejet_thirdjet_phi_diff_wta->GetMaximum() > max_val_phidiff_wta_3jets) {
      max_val_phidiff_wta_3jets = threejet_thirdjet_phi_diff_wta->GetMaximum();
      threejet_thirdjet_phi_diff_wta->SetMaximum(1.2*max_val_phidiff_wta_3jets);
    }
  }

  leg_phidiff_wta_3jets->Draw("SAMES");
  phidiff_wta_compare_3jets->Write();
  phidiff_wta_compare_3jets->Print("zprime_phidiff_wta_compare_3jets.eps", "eps");

  TCanvas *tau3diff_wta_compare_3jets = new TCanvas("tau3diff_wta_compare_3jets", "tau3diff_wta_compare_3jets", 600, 600);
  tau3diff_wta_compare_3jets->cd();
  tau3diff_wta_compare_3jets->SetLogy();
  double max_val_tau3diff_wta = 0;
  TLegend *leg_tau3diff_wta_3jets = new TLegend(0.7, 0.7, 0.9, 0.9);

  for (int B = 0; B < betalist.size(); B++) {
    double beta = betalist[B];
    ostringstream ss;
    ss << beta;

    TH1* tau3_diff_wta = (TH1*)tau3_diff_wta_hists.At(B);
    double tau3_diff_wta_scale = 1/tau3_diff_wta->Integral(0, tau3_diff_wta->GetNbinsX() + 1);
    tau3_diff_wta->Scale(tau3_diff_wta_scale);
    if (B == 0) tau3_diff_wta->SetLineColor(kYellow);
    if (B == 1) tau3_diff_wta->SetLineColor(kGreen);
    if (B == 2) tau3_diff_wta->SetLineColor(kRed);
    if (B == 3) tau3_diff_wta->SetLineColor(kBlue);
    tau3_diff_wta->SetTitle("#tau_{3} difference of hardest anti-kT (WTA) and 3-jettiness jet");
    tau3_diff_wta->GetXaxis()->SetTitle("#tau_{3}_{akt} - #tau_{3}_{njet}");
    tau3_diff_wta->GetYaxis()->SetTitle("Relative Occurrence");
    tau3_diff_wta->GetYaxis()->SetTitleOffset(1.5);
    tau3_diff_wta->Draw("SAMES");
    leg_tau3diff_wta_3jets->AddEntry(tau3_diff_wta, "#beta = " + (TString)ss.str());
    if (tau3_diff_wta->GetMaximum() > max_val_tau3diff_wta) {
      max_val_tau3diff_wta = tau3_diff_wta->GetMaximum();
      tau3_diff_wta->SetMaximum(1.2*max_val_tau3diff_wta);
    }
  }

  leg_tau3diff_wta_3jets->Draw("SAMES");
  tau3diff_wta_compare_3jets->Write();
  tau3diff_wta_compare_3jets->Print("zprime_tau3diff_wta_compare_3jets.eps", "eps");



  // double zprime_invmass_akt_scale = 1/zprime_invmass_akt->Integral();
  // double zprime_invmass_njet_beta1_scale = 1/zprime_invmass_njet_beta1->Integral();
  // double zprime_invmass_njet_beta2_scale = 1/zprime_invmass_njet_beta2->Integral();
  // zprime_invmass_akt->Scale(zprime_invmass_akt_scale);
  // zprime_invmass_njet_beta1->Scale(zprime_invmass_njet_beta1_scale);
  // zprime_invmass_njet_beta2->Scale(zprime_invmass_njet_beta2_scale);

  // TCanvas *mass_compare = new TCanvas("mass_compare", "Mass Comparison", 600, 600);
  // mass_compare->cd();
  // zprime_invmass_akt->SetLineColor(kRed);
  // zprime_invmass_njet_beta1->SetLineColor(kBlue);
  // zprime_invmass_njet_beta2->SetLineColor(kGreen);
  // zprime_invmass_akt->Draw();
  // zprime_invmass_njet_beta1->Draw("same");
  // zprime_invmass_njet_beta2->Draw("same");
  // TLegend *leg_mass = new TLegend(0.1, 0.7, 0.3, 0.9);
  // leg_mass->AddEntry(zprime_invmass_akt, "Anti-KT", "L");
  // leg_mass->AddEntry(zprime_invmass_njet_beta1, "Njettiness (#beta = 1)", "L");
  // leg_mass->AddEntry(zprime_invmass_njet_beta2, "Njettiness (#beta = 2)", "L");
  // leg_mass->Draw();
  // mass_compare->Write();

  // TCanvas *perp_compare = new TCanvas("perp_compare", "perp Comparison", 600, 600);
  // perp_compare->cd();
  // jet_ptspectrum_akt->SetLineColor(kRed);
  // jet_ptspectrum_njet->SetLineColor(kBlue);
  // jet_ptspectrum_akt->Draw();
  // jet_ptspectrum_njet->Draw("same");
  // TLegend *leg_perp = new TLegend(0.1, 0.7, 0.3, 0.9);
  // leg_perp->AddEntry(jet_ptspectrum_akt, "Anti-KT", "L");
  // leg_perp->AddEntry(jet_ptspectrum_njet, "Njettiness", "L");
  // leg_perp->Draw();
  // perp_compare->Write();

  // TCanvas *phi_compare = new TCanvas("phi_compare", "phi Comparison", 600, 600);
  // phi_compare->cd();
  // jet_phispectrum_akt->SetLineColor(kRed);
  // jet_phispectrum_njet->SetLineColor(kBlue);
  // jet_phispectrum_akt->Draw();
  // jet_phispectrum_njet->Draw("same");
  // TLegend *leg_phi = new TLegend(0.1, 0.7, 0.3, 0.9);
  // leg_phi->AddEntry(jet_phispectrum_akt, "Anti-KT", "L");
  // leg_phi->AddEntry(jet_phispectrum_njet, "Njettiness", "L");
  // leg_phi->Draw();
  // phi_compare->Write();


  //Prevent memory leaks  
  delete outFile;
  delete jetDef;
  // delete recombScheme;
  return 0;
}