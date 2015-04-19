/* 
This is the primary header file for the Global Soft Drop Algorithm, originally created by TJ Wilkason.

The Global Soft Drop Algorithm ("DropCone") utilizes the principle of Soft Drop (arXiv:1402.2657) over the whole event. 
It runs by first clustering all particles in the event into a single particle through a kT-like algorithm and a WTA recombiner.
The algorithm then uses a priority queue to recursively unwind the clustering sequence, dropping all particles 

*/

// COMMENT YOUR DAMN CODE BEFORE YOU FORGET

#ifndef _GlobalSoftDrop_hh_
#define _GlobalSoftDrop_hh_

//FASTJET HEADERS
// #include "FastJet3.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/WinnerTakeAllRecombiner.hh"

#include <cmath>
#include <queue>
#include <unordered_map>

using namespace std;
using namespace fastjet;

class GlobalSoftDrop;

struct BranchFitness {
	PseudoJet particle;
	double fitness;
	double zfrac;
};

class CompareJetBranch {
public:
	// Returns true if b1 is lower in the queue than b2
	bool operator()(BranchFitness b1, BranchFitness b2) {
		return (b1.fitness < b2.fitness);
	}
};

class BranchQueue {
public:
	BranchQueue() {}
	~BranchQueue() {delete seedJet_clustSeq;}
	// vector<PseudoJet> runBranchQueue(PseudoJet leading_particle, int njets, double alpha, double zcut);
	vector<PseudoJet> runBranchQueue(std::vector<fastjet::PseudoJet> input_particles, fastjet::JetAlgorithm jet_algorithm, int njets, double alpha, double zcut);

private:
	PseudoJet seedJet;
	ClusterSequence* seedJet_clustSeq;

	// step 1: cluster all particles in the event into a single "jet" through WTA recombination
	void initialCluster(std::vector<fastjet::PseudoJet> input_particles, fastjet::JetAlgorithm jet_algorithm);

	vector<BranchFitness> splitBranch(BranchFitness branch, double alpha);
	vector<PseudoJet> queue2Vector(priority_queue<BranchFitness, vector<BranchFitness>, CompareJetBranch> myQueue);
};

class GlobalSoftDrop {
public:
	GlobalSoftDrop(vector<PseudoJet> input_particles, 
					JetAlgorithm jet_algorithm,
					int njets, double alpha, double zcut) 
	: _input_particles(input_particles), _jet_algorithm(jet_algorithm),
	_njets(njets), _alpha(alpha), _zcut(zcut) {}

	void runAlgorithm() {
		findAxes();
		calculateRadii();
		clusterJets();
	}

	vector<PseudoJet> getAxes() { return myAxes; }
	vector<PseudoJet> getJets() { return myJets; }
	double getRadius() { return jet_radius; }

private:
	vector<PseudoJet> _input_particles;
	JetAlgorithm _jet_algorithm;
	int _njets;
	double _alpha;
	double _zcut;

	// PseudoJet seedJet;
	// ClusterSequence _seedJet_clustSeq;
	vector<PseudoJet> myAxes;
	vector<PseudoJet> myJets;
	double jet_radius;

	// step 1: use the priority queue to find the desired number of axes in the event
	void findAxes();
	// step 2: use the distance between the axes to calculate the desired radius for each of the jets
	void calculateRadii();
	// step 3: cluster particles into jets through nearest-neighbor clustering (auto-tessellation)
	void clusterJets();
};

#endif  // _GlobalSofDrop_hh__