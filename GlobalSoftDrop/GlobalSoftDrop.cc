#include "GlobalSoftDrop.hh" 

using namespace std;
using namespace fastjet;


vector<BranchFitness> BranchQueue::splitBranch(BranchFitness branch, double alpha) {
	PseudoJet particle = branch.particle;
	PseudoJet parent0, parent1;
	particle.has_parents(parent0, parent1);

	BranchFitness parent0Branch, parent1Branch;

	if (parent0.has_structure() && parent1.has_structure()) {

		// double angle_distance = pow((pow(parent0.eta() - parent1.eta(), 2) + pow(parent0.phi() - parent1.phi(), 2)), 0.5);
		double angle_distance = parent0.delta_R(parent1);
		double parent0frac = (parent0.perp() / (parent0.perp() + parent1.perp()))*pow(angle_distance, alpha);
		double parent1frac = (parent1.perp() / (parent0.perp() + parent1.perp()))*pow(angle_distance, alpha);

		double parent0Fitness = parent0.exclusive_subdmerge(1);
		double parent1Fitness = parent1.exclusive_subdmerge(1);

		parent0Branch = {parent0, parent0Fitness, parent0frac};
		parent1Branch = {parent1, parent1Fitness, parent1frac};
	}

	vector<BranchFitness> parentBranches;
	parentBranches.push_back(parent0Branch);
	parentBranches.push_back(parent1Branch);

	return parentBranches;
}

vector<PseudoJet> BranchQueue::queue2Vector(priority_queue<BranchFitness, vector<BranchFitness>, CompareJetBranch> myQueue) {
  vector<PseudoJet> myVector;
  int size = myQueue.size();
  for (int i = 0; i < size; i++) {
    BranchFitness tempBlob = myQueue.top();
    myVector.push_back(tempBlob.particle);
    myQueue.pop();
  }
  for (int i = 0; i < size; i++) {
    BranchFitness tempBlob = {myVector[i], myVector[i].exclusive_subdmerge(1)};
    myQueue.push(tempBlob);
  }
  return myVector;
}

// void GlobalSoftDrop::initialCluster() {
void BranchQueue::initialCluster(std::vector<fastjet::PseudoJet> input_particles, fastjet::JetAlgorithm jet_algorithm) {

	double Rparam = JetDefinition::max_allowable_R;
	Strategy strategy = Best;
	JetDefinition::Recombiner *recombScheme = new fastjet::contrib::WinnerTakeAllRecombiner();
	JetDefinition *jetDef = new JetDefinition(jet_algorithm, Rparam, recombScheme, strategy);

	seedJet_clustSeq = new ClusterSequence(input_particles, *jetDef);
	vector <PseudoJet> inclusiveJets = seedJet_clustSeq->inclusive_jets(0.0);
	seedJet = inclusiveJets[0];
}

vector<PseudoJet> BranchQueue::runBranchQueue(std::vector<fastjet::PseudoJet> input_particles, fastjet::JetAlgorithm jet_algorithm, int njets, double alpha, double zcut) {
	
	initialCluster(input_particles, jet_algorithm);

	priority_queue<BranchFitness, vector<BranchFitness>, CompareJetBranch> branchQueue;

	double firstBranchfitness = seedJet.exclusive_subdmerge(1);

	BranchFitness firstBranch = {seedJet, firstBranchfitness, 0.0};
	branchQueue.push(firstBranch);

	vector<PseudoJet> branchAxes;
	BranchFitness currentBranch;

	vector<PseudoJet> tempAxes;
	int max_axes_size = 0;

	double Rmax = std::numeric_limits<double>::max();

	while (branchQueue.size() < (unsigned)njets) {

		if (branchQueue.size() == 0) {
			cerr << "Error: branch queue empty." << endl;
			break;
		}

		if (branchQueue.size() >= (unsigned)max_axes_size) {
			max_axes_size = branchQueue.size();
			tempAxes = queue2Vector(branchQueue);
      		// branchQueue = vector2Queue(tempAxes);
		}

		currentBranch = branchQueue.top();
		branchQueue.pop();

		vector<BranchFitness> parentBranches = splitBranch(currentBranch, alpha);
		if (parentBranches[0].particle.has_structure() && parentBranches[1].particle.has_structure()) {
			double parent0splitting = parentBranches[0].particle.delta_R(firstBranch.particle);
			double parent1splitting = parentBranches[1].particle.delta_R(firstBranch.particle);
			if (parentBranches[0].zfrac < zcut || parent0splitting > Rmax) {
				branchQueue.push(parentBranches[1]);
			}
			else if (parentBranches[1].zfrac < zcut || parent1splitting > Rmax) {
				branchQueue.push(parentBranches[0]);
			}
			else {
				branchQueue.push(parentBranches[0]);
				branchQueue.push(parentBranches[1]);
			}
		}
	} 

	int n_found = branchQueue.size();
	if (n_found == 0) branchAxes = tempAxes;

	for (int i = 0; i < n_found; i++) {
		BranchFitness final_branch = branchQueue.top();    
		branchAxes.push_back(final_branch.particle);
		branchQueue.pop();
	}
	return branchAxes;
}

void GlobalSoftDrop::findAxes() {
	vector<PseudoJet> final_axes;
	BranchQueue *branch_queue = new BranchQueue();
	// final_axes = branch_queue->runBranchQueue(seedJet, _seedJet_clustSeq, _njets, _alpha, _zcut);
	final_axes = branch_queue->runBranchQueue(_input_particles, _jet_algorithm, _njets, _alpha, _zcut);

	myAxes = final_axes;
	delete branch_queue;
}

void GlobalSoftDrop::calculateRadii() {

  	std::string bitmask(2, 1); // K leading 1's
  	bitmask.resize(myAxes.size(), 0); // N-K trailing 0's

  	double avg_distance = 0;
  	double max_distance = 0;
  	int axis_counter = 0;

  	// DECIDE ON THE DEFAULT VALUE -- TJW
  	if (_njets == 1) jet_radius = 0.5;
  	else {

	  	do {
	  		vector<int> axis_indices;
			for (unsigned int i = 0; i < myAxes.size(); ++i) { // [0..N-1] integers
				if (bitmask[i]) axis_indices.push_back(i);
			}
			if (axis_indices.size() > 1) {
				double distance = myAxes[axis_indices[0]].delta_R(myAxes[axis_indices[1]]);
				avg_distance += distance;
				if (distance > max_distance) max_distance = distance;

				axis_counter++;
			}
		} while (std::prev_permutation(bitmask.begin(), bitmask.end()));

		avg_distance = avg_distance/(2*axis_counter); 
		jet_radius = avg_distance;
	}
	// jet_radius = 0.5;
}

void GlobalSoftDrop::clusterJets() {

	// JetDefinition::Recombiner *recombScheme = new fastjet::contrib::WinnerTakeAllRecombiner();

	vector<PseudoJet> final_jets;

  	// INITIALIZE JETS
	for (unsigned int i = 0; i < myAxes.size(); i++) {
		PseudoJet empty_jet(0,0,0,0);
		final_jets.push_back(empty_jet);
	}

  	// TAKE INITIAL PARTICLES AND CLUSTER THEM ACCORDING TO NEAREST AXIS
	for (unsigned int i_particles = 0; i_particles < _input_particles.size(); i_particles++) {
		double min_axis_distance = std::numeric_limits<int>::max();
		int min_axis_index = -1;
		for (unsigned int i_axis = 0; i_axis < myAxes.size(); i_axis++) {
			double axis_distance = _input_particles[i_particles].delta_R(myAxes[i_axis]);
			if ((axis_distance < min_axis_distance) && (axis_distance < jet_radius)) {
				min_axis_distance = axis_distance;
				min_axis_index = i_axis;
			}
		}
		if (min_axis_index != -1) {
			// final_jets[min_axis_index] = join(final_jets[min_axis_index], _input_particles[i_particles], *recombScheme);
			final_jets[min_axis_index] = join(final_jets[min_axis_index], _input_particles[i_particles]);
		}
	}

	myJets = final_jets;
}