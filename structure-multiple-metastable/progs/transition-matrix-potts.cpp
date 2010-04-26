// Calculate the transition matrix occurring in the master equation for
// small Potts systems

// Started 15/08/2006

#include <iostream>
#include <vector>
#include <cmath>

#include <fstream>
#include <sstream>

#include "IntegerValuedLattice.H"
#include "ParameterFileParser.h"

using namespace std;

IntegerValuedLattice *Lattice;

double **interactions;
double *h;

int *magnetisation;
double interaction_energy;

double delta_H;
double sum_off_diag;

ofstream outfile;
ofstream labelfile;

bool make_diagonal = false;
// This determines the order in which configurations are generated
// If true, it generates configs in an order to make the matrix as
//   diagonal as possible
// If false, it generates them in numerical order 

bool symmetric = true;
// if true, generates the symmetric version of the matrix
// otherwise generates the normal Glauber dynamics version
// (c.f. R. Melin, J Physique (France) I, 1996)

// note that the symmetric version is *not* a stochastic matrix,
// i.e. does not have row sums equal to 1
// since the diagonal elements are the SAME as in the original
// Glauber matrix

// must account for this in the calculation


int glut_simulation_window = 0;
int q;   // no. of spin states allowed for each spin (Ising=2)

int Lx;
int Ly;

double beta;

int z;  // lattice coord number

int N;
int no_configs;

vector<int> labels;   // ordered array of labels
int* location;        // location of each integer in the label array
int max_power; 
int mask;

double* Hamiltonian;


//double** transition_matrix;



double calc_energy() {
	// calculate total energy of a configuration

	interaction_energy = 0;
	double field_energy = 0.0;
	
	for (int i=0; i<N; i++) {
		
		int value = Lattice->sites[i];
		
		for (int j=0; j<z; j++) {
			
			int neighbour = Lattice->neighbours[i][j];
			int neighbour_value = Lattice->sites[neighbour];
			
			interaction_energy -= interactions[value][neighbour_value];
			
		}
		
		field_energy -= h[value];
		
	}
	
	return (0.5*interaction_energy + field_energy);

		// double counted all bonds
}


void calc_magnetisations() {

	int value;

	for (int j=0; j<q; j++) {
		magnetisation[j] = 0;
	}

	for (int i=0; i<N; i++) {
		value = Lattice->sites[i];
		magnetisation[value] ++;
	}

}
  


double calc_transition_prob(int spin, int new_value) {
	
	delta_H = 0.0;

	int neighbour;
	int neighbour_value;
	
	int old_value = Lattice->sites[spin];
	
	double transition_prob = 0.0;
	
	
	for (int i=0; i<z; i++) {  // loop through neighbours
		neighbour = Lattice->neighbours[spin][i];
		neighbour_value = Lattice->sites[neighbour];
		
		delta_H += 
			interactions[old_value][neighbour_value] 
			- interactions[new_value][neighbour_value];
		
	}
	
	delta_H += h[old_value] - h[new_value];
	
	if (delta_H < 0.0) transition_prob = 1.0;
	else transition_prob = exp(-beta * delta_H);

	transition_prob /= N * (q-1);  // proposal probability
	
	sum_off_diag += transition_prob;
	

	if (symmetric) {
		transition_prob *= exp(0.5*beta*delta_H);
	}

	return transition_prob;

}





void initialise_lattice(int config_number) {
	// initialise the lattice with the configuration config_number treated as ternary

	using namespace std;

	mask = max_power;

	//cout << "\n" << "Initialising lattice with configuration " << config_number << "\n";

	for (int position = 0; position < N; position++) {  // position in lattice
		//cout << "mask = " << mask 
		//		 << "; config_number = " << config_number
		//		 << "; position=" << position << "\n";
		
		int current_value = config_number / mask;
	
		Lattice->sites[position] = current_value;   // ternary digit  
		//cout << "Setting to " << config_number / mask << "\n";

		config_number -= current_value*mask;
		mask /= q;

	}

	//Lattice->PlotGridRectangular();

}



void find_transitions(int config_number) {
	// find transitions and their probabilities from configuration "config_number"

	// note that transition from a state to itself must be 1 minus the sum of the others

	using namespace std;

	sum_off_diag = 0.0;

	
	//cout << "\n" << "Find transitions from " << config_number << "; location is " << location[config_number] << "\n";

	initialise_lattice(config_number);

	Hamiltonian[config_number] = -beta * calc_energy();

	calc_magnetisations();

	if (!make_diagonal) {

		labelfile << config_number;

		for (int i=0; i<q-1; i++) {
			labelfile << "\t" << magnetisation[i];
		}
		
		labelfile << "\t" << interaction_energy << "\n";
	}


	int current_config_number = config_number;

	if (make_diagonal) {
		if (location[config_number] < 0) {
			labels.push_back(config_number);
			location[config_number] = labels.size() - 1;
		}
	}


	mask = max_power;

	for (int spin_to_flip = 0; spin_to_flip < N; spin_to_flip++) {

		int current_value = current_config_number / mask;
		// just used for extracting digits to base q

		current_config_number -= (current_value * mask);

		for (int i=1; i<=q-1; i++) {
			int new_value = (current_value + i) % q;

			double transition_prob = calc_transition_prob(spin_to_flip, new_value);

			int new_config = config_number + (new_value - current_value)*mask;

			if (make_diagonal) {
				if (location[new_config] < 0 ) {
					labels.push_back(new_config);
					location[new_config] = labels.size() - 1;
				}
			}


			//transition_prob /= N*(q-1);   // proposal probability
			
			
			//cout << "Transition " << location[config_number] << " -> " << location[new_config] 
			//		 << "; prob. " << transition_prob << "\n";
			
			
			//transition_matrix[location[config_number]][location[new_config]] = transition_prob;
			
			//sum_off_diag += transition_prob;
			
			//cout << "Transtion " << config_number << " -> " << new_config << ", prob. "
			//			 << transition_prob << "\n";
			
			if (make_diagonal) {
				outfile << location[config_number]+1 << "\t" 
								<< location[new_config]+1 << "\t" 
								<< transition_prob << "\n";
			}

			else {
				outfile << config_number + 1 << "\t"
								<< new_config+1 << "\t"
								<< transition_prob << "\n";
			}

			
			//transition_matrix[location[config_number]][location[config_number]] = 1.0 - sum_off_diag;
			
			



		}
		
		mask = mask / q;


	}

	//cout << "Transition " << location[config_number] << " -> " << location[config_number] 
	//		 << "; prob. " << 1.0 - sum_off_diag;

	if (make_diagonal) {
		outfile << location[config_number]+1 << "\t" 
						<< location[config_number]+1 << "\t" 
						<< 1.0-sum_off_diag << "\n";
	}

	else {
		outfile << config_number + 1 << "\t"
						<< config_number + 1 << "\t"
						<< 1.0-sum_off_diag << "\n";
	}

	//cout << location[config_number]+1 << "\t" << sum_off_diag << "\n";

}


	void readParameters (string filename) {

		ParameterFileParser file (filename);

		int file_q;
		
		file.parse (file_q);
		if (file_q != q) {
			cout << "\n q from parameter file does not match; not reading param file!\n\n";
			return;
		}
		
		file.parse(beta);

		for (int i=0; i<q; i++) {
			file.parse(h[i]);
		}
		
		string string1, string2;
		
		while (!file.eof()) {
		
			file.parse(string1, string2);
		
			int i = int(string1[0]) - int('0');
			int j = int(string1[1]) - int('0');
			
			double J = atof(string2.c_str());
			
			cout << "Updating interaction " << i << ", " << j << " to value " << J << "\n";
			interactions[i][j] = J;
		}
	}



void readCommandLineArguments(int argc, char *argv[]) {
	
	if (argc < 6) {
		cout << "Not enough parameters!\n" 
				 << "Syntax: transition-matrix-potts q Lx Ly parameter-file symmetric [b=] [h1=] ... [hq=] [J00=] ..." << "\n";
		cout << "Later parameters override earlier ones" << "\n";
		exit(1);
  }

	cout << "argc = " << argc << "\n";

	q = atoi(argv[1]);
  Lx = atoi(argv[2]);
	Ly = atoi(argv[3]);

		
	/// need the above parameters to define arrays of correct size:
	
	interactions = new double*[q];
	for (int i=0; i<q; i++) {
		interactions[i] = new double[q];
		
		for (int j=0; j<q; j++) {
			interactions[i][j] = 0.0;
		}

		interactions[i][i] = 1.0;
		
	}

	magnetisation = new int[q];

	
	h = new double[q]; 
	for (int i=0; i<q; i++) {
		h[i] = 0.0;
	}
	

	
	/// carry on processing command line params
	
	readParameters(argv[4]);
	
	symmetric = bool( atoi(argv[5]) );
	
	cout <<  "Using parameters:\n"
			 << "q = " << q << "\n"
			 << "size = " << Lx << " x " << Ly << "\n"
			 << "parameter file: " << argv[4] << "\n"
			 << "symmetric: " << symmetric << "\n"; 
			

	


	int j, k, l;
	
	for (int i=6; i<argc; i++) {
		cout << "Reading parameter " << i << "\n";

		switch(argv[i][0]) {

		case 'b':
			beta = atof(argv[i]+2);
			break;

		case 'h':
			j = int(argv[i][1]) - int('0');
			h[j] = atof(argv[i]+3);
			break;

		case 'J':
			k = int(argv[i][1]) - int('0');
			l = int(argv[i][2]) - int('0');
			interactions[k][l] = atof(argv[i]+4);
			break;

		//default:
		//	readParameterFile(argv[i]);
		//	break;
		}
	}


}



void outputParametersAsComments(ostream& paramfile) {
  // output params to file with marker at beginning of each line
	// use '%' for MatrixMarket format
	// Can set gnuplot so uses '%' for comment as well:
	// set datafile commentschars "#%"

	paramfile <<  "%" << "\n";
	paramfile << "% Generalised Potts model with parameters:\n";

  paramfile << "% " << q << "\tq\n";
  paramfile << "% " << Lx << "\tLx\n";
	paramfile << "% " << Ly << "\tLy\n";
  paramfile << "% " << beta << "\tbeta\n";
  paramfile << "% " << "\n";
	
  for (int i=0; i<q; i++) {
		paramfile << "% " << h[i] << "\t" <<"h[" << i << "]\n";
  }
	
  paramfile << "%\n";
  
  for (int i=0; i<q; i++) {
		for (int j=i; j<q; j++) {
			paramfile << "% " << interactions[i][j] << "\t"
								<< "J(" << i << "," << j << ")\n";
	  
	  
	  
		}
  }

  paramfile << "%\n";
  
}



int main(int argc, char *argv[]) {


	using namespace std;

	// command line arguments:
	// q Lx Ly parameters.txt
	
	

	readCommandLineArguments(argc, argv);

	Lattice = new IntegerValuedLattice(glut_simulation_window,
																		 0,
																		 Lx, Ly, q+1);

	N = Lattice->N;
	z = Lattice->z;

	outfile.open("output-parameters.txt");
	outputParametersAsComments(outfile);
	outfile.close();


	/*

	// Ising for comparison:
	q=2;
	interactions[0][0] = interactions[1][1] = 1.0;
	interactions[0][1] = interactions[1][0] = -1.0;

	h[0] = -0.5;
	h[1] = 0.5;

	//beta = 1.0;
	
	
	


	//q = 3;
	

	
	
	// Sequence of 3 metastable states 1->2->3 :
	
	interactions[0][0] = interactions[1][1] = interactions[2][2] = 1.0;
	interactions[0][1] = interactions[1][2] = 0.1;
	interactions[0][2] = -1.0;

	h[0] = 0.0;
	h[1] = 0.5;
	h[2] = 0.8;
	
	

	
	
	
	// Competing states:
	interactions[0][0] = interactions[1][1] = interactions[2][2] = 1.0;
	interactions[0][1] = interactions[0][2] = 0.1;
	interactions[1][2] = -1.0;

	h[0] = 0.0;
	h[1] = 0.5;
	h[2] = 0.5;
	


	
	// Competing states with positive 1--2 interaction:
	interactions[0][0] = interactions[1][1] = interactions[2][2] = 1.0;
	interactions[0][1] = interactions[0][2] = 0.1;
	interactions[1][2] = 0.1;

	h[0] = 0.0;
	h[1] = 0.5;
	h[2] = 1.0;
	
	*/
	

	cout << "Interaction matrix: " << "\n";
	// symmetrise interaction matrix:
	for (int i=0; i<q; i++) {
		for (int j=0; j<i; j++) {
			interactions[i][j] = interactions[j][i];
		}

		for (int j=0; j<q; j++) {
			cout << interactions[i][j] << "\t";
		}

		cout << "\n";

	}
	



	cout << "beta = " << beta << "\n";

	no_configs = 1;
	for (int i=0; i<N; i++) {
		no_configs = no_configs * q;
	}


	location = new int[no_configs];
	for (int i=0; i<no_configs; i++) {
		location[i] = -1;  // not in label list yet
	}

	Hamiltonian = new double[no_configs];

	
	/*
	transition_matrix = new double*[no_configs];
	for (int i=0; i<no_configs; i++) {
		transition_matrix[i] = new double[no_configs];
		for (int j=0; j<no_configs; j++) {
			transition_matrix[i][j] = 0.0;
		}
	}
	*/


	max_power = no_configs / q;    // for masking

	cout << "\n" << N << " spins in lattice; " << no_configs << " configurations" << "\n";


	labels.push_back(no_configs-1);
	location[no_configs-1] = 0;

	//for (int i=no_configs-1; i>=0; i--) {
	//	find_transitions(i);         // find_transitions generates new entries in labels!
	//}
	

	string filename;
	
	if (symmetric) {
		filename = "potts-matrix-symmetric.dat";
	}

	else {
		filename = "potts-matrix.dat";
	}

	outfile.open(filename.c_str()); 
	outfile.precision(17);

	labelfile.open("labels-potts.dat");

	cout << "Writing to file: " << filename << " in sparse (MatrixMarket) form" << "\n";
 
	outfile << "%%MatrixMarket matrix coordinate real general" << "\n";
	outfile << "% Sparse matrix stored in MatrixMarket format" << "\n";

	outputParametersAsComments(outfile);

	outfile << no_configs << "\t" << no_configs << "\t" 
					<< no_configs*(N*(q-1)+1) << "\n";
	
	// the last number is the number of transitions, i.e. non-zero entries:
	// each row (of which there are no_configs) has (N*(q-1))+1 transitions:
	// flipping each of N spins to each of q-1 other states, 
	// or flipping no spins and remaining the same


	if (make_diagonal) {
		for (int i=0; i<no_configs; i++) {
			if (i%5000==0)
				cout << i << "\n";
			find_transitions(labels[i]);   // do in order of label generation!
		}
	}

	else {
		for (int i=0; i<no_configs; i++) {
			if (i%5000==0) 
				cout << i << "\n";
			find_transitions(i);
		}
	}

	cout << "Label array: " << "\n";

	for (int i=0; i<labels.size(); i++) {
		cout << i << ": " << labels[i] << "\n";
	}


	// calculate log Z and output P_0^{-1/2}

	double f;
	double f_max = -100000.0;
	double log_Z = 0.0;

	for (int i=0; i<no_configs; i++) {
		f = Hamiltonian[i];
		if (f > f_max) f_max = f;
	}
	
	for (int i=0; i<no_configs; i++) {
		f = Hamiltonian[i];
		log_Z += exp(f - f_max);
	}

	log_Z = f_max + log(log_Z);

	filename = "eqm-distn.dat";
	ofstream eqm_file(filename.c_str());
	cout << "Writing eqm distn to " << filename << "\n";
 
	eqm_file.precision(17);

	eqm_file << "# log[(P_0)_i ^ {1/2}]" << "\n";
	for (int i=0; i<no_configs; i++) {
		f = Hamiltonian[i];

		//eqm_file << 0.5*(f - log_Z) << "\n";

		if (symmetric) {
			eqm_file << exp(0.5 * (f - log_Z)) << "\n";
		}

		else {
			eqm_file << exp(f - log_Z) << "\n";
		}
	}
	
	


	
	/*
	cout << "\n" << "Transition matrix: " << "\n";
	for (int i=0; i<no_configs; i++) {
			cout << i << ": \t";

		for (int j=0; j<no_configs; j++) {
			if (transition_matrix[i][j] > 1.0E-8) cout << "1";
			else cout << ".";

			cout << " ";

			// outfile << transition_matrix[i][j] << " ";

			if (transition_matrix[i][j] > 1.0E-14) {
				outfile << i+1 << "\t" << j+1 << "\t" << transition_matrix[i][j] << "\n";
				// Mathematica needs arrays starting at 1, hence the +1
			}

		}

		cout << "\n";
		//outfile << "\n";
	}
	*/

}
