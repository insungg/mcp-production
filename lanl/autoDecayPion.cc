#include <iostream>
#include <fstream>

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TCanvas.h"

using namespace std;
using namespace ROOT::Math;

// constants
const double mpi   		= 134.9768; // MeV, charged pion mass
const double mp    		= 938.272;
const double alpha 		= 0.0072973526; // 1/137
const double PI    		= 3.14159265358979;
const double conv2rad 	= PI/180.0;

double envelope = 150; // enveloping function of ddBrPi2gxx for rejection sampling

// model parameter, like coupling constants
const double epsilon   = 1.0;
const double BrPi2gg   = 1e10; // Br(pion -> gamam gamma) 
double mchi;


// store particle information
struct Particle {
	double mass;
	double flightTime;
	TLorentzVector momentum;
	TLorentzVector startPosition;
	TLorentzVector endPosition;
};

//  sqrt(Kallen function(x,y,z))/2x = lambda, it defines momentum in two-body decay
double lambda(double x, double y, double z) {
	return sqrt(pow(x,4) + pow(y,4) + pow(z,4) -2*pow(x*y,2) - 2*pow(y*z,2) - 2*pow(z*x, 2))/(2*x);
}

// consider only massless dark vector mediator, so no on-shell contribution
double ddBrPi2gxx(double s, double theta) {
	return sin(theta) * pow(epsilon, 2) * alpha / (4*PI*s) * pow(1-s/pow(mpi, 2), 3) * sqrt(1-4*pow(mchi, 2)/s) * (1 - 4*pow(mchi, 2)/s) * pow(sin(theta), 2) * BrPi2gg;
}


void autoDecayPion() { // mass of chi in MeV
	// initialize
	ifstream pionTxt("Pion_Production.txt");
	if (!pionTxt.is_open()) {
		cout << "Failed to open pion data" << endl;
		return -1;
	}
	ifstream massTxt("mass3.txt");
	if (!massTxt.is_open()) {
		cout << "Failed to open mass data" << endl;
		return -1;
	}
	ofstream ageoTxt("ageo_e73.txt");
	if (!ageoTxt.is_open()) {
		cout << "Failed to open ageo data" << endl;
		return -1;
	}

	vector<vector<double>> records;
	int nMass = 22;
	string line;
	TRandom3 rand;
	for (int i = 0; i < nMass; i++) {
		if (!getline(massTxt, line)) break;
		mchi = atof(line.c_str());
		double PX, PY, PZ, PP, M, E; 
		int hitCounter10m = 0;
		int hitCounter35m = 0;
		int hitCounter60m = 0;
		int hitCounter100m = 0;


		// loop through pion data 
		double momP, momTheta, momPhi;
		while (pionTxt >> momP >> momTheta >> momPhi) {
			// sample s(=off-shell photon mass) and theta(angle between mcp momentum 
			// in the rest frame of off-shell photon and its z axis
			double s, theta;
			while (true) {
				s 	  = rand.Uniform(4.0*pow(mchi, 2), pow(mpi, 2));
				theta = rand.Uniform(0.0, PI); 
				double y = ddBrPi2gxx(s, theta);
				// cout << y << endl;
				if (y >= rand.Uniform(0, 1) * envelope) {
					// renew envelope when y > envelope
					if (y > envelope)
						envelope = y;
					break;
				}
			}

			// two-body decay off-shell photon V -> xxbar in the rest frame of V
			double P   = lambda(sqrt(s), mchi, mchi);
			double phi = rand.Uniform(0, 2.0*PI);
			Particle mcp1; 
			mcp1.momentum.SetPxPyPzE(P*sin(theta)*cos(phi), P*sin(theta)*sin(phi), P*cos(theta), sqrt(P*P + mchi*mchi));
			// cout << mcp1.momentum.Px() << " " << mcp1.momentum.Py() << " " << mcp1.momentum.Pz() << endl;
			Particle mcp2;
			mcp2.momentum.SetPxPyPzE(-P*sin(theta)*cos(phi), -P*sin(theta)*sin(phi), -P*cos(theta), sqrt(P*P + mchi*mchi));

			// boost V from its rest frame along z direction
			double vP    = lambda(mpi, sqrt(s), 0); // pion -> V + gamma, mass of gamma = 0
			double vBeta = vP / sqrt(vP*vP + s); 
			mcp1.momentum.Boost(0, 0, vBeta);
			mcp2.momentum.Boost(0, 0, vBeta);
			// pion -> V + gamma is an isotropic two-body decay processs
			double vTheta = rand.Uniform(0.0, PI);
			double vPhi   = rand.Uniform(0.0, 2*PI);
			// rotate z axis of V back to pion rest frame
			mcp1.momentum.RotateZ(vTheta);
			mcp1.momentum.RotateY(vPhi);
			mcp2.momentum.RotateZ(vTheta);
			mcp2.momentum.RotateY(vPhi);
			// boost pion rest frame to lab frame
			double momE  = sqrt(pow(momP, 2) + pow(mpi, 2));
			double momBx = momP*sin(momTheta)*cos(momPhi) / momE;
			double momBy = momP*sin(momTheta)*sin(momPhi) / momE;
			double momBz = momP*cos(momTheta) / momE;
			mcp1.momentum.Boost(momBx, momBy, momBz);
			mcp2.momentum.Boost(momBx, momBy, momBz);

			PX = mcp1.momentum.Px();
			PY = mcp1.momentum.Py();
			PZ = mcp1.momentum.Pz();
			PP = mcp1.momentum.P();
			M  = mcp1.momentum.M();
			E  = mcp1.momentum.E();
			
			double bx, by, bz, scale, y, z;
			bx = PX/E; by = PY/E; bz = PZ/E;
			scale = 10.0/bx; // time to arriave at the target multiplied by c
			
			y = by * scale; // by * c * scale / c = by * scale
			z = bz * scale; 

			if (abs(y) < 0.25 && abs(z) < 0.2)
				hitCounter10m++;

		
			bx = PX/E; by = PY/E; bz = PZ/E;
			scale = 35.0/bx; // time to arriave at the target multiplied by c
			
			y = by * scale; // by * c * scale / c = by * scale
			z = bz * scale; 

			if (abs(y) < 0.25 && abs(z) < 0.2)
				hitCounter35m++;

			bx = PX/E; by = PY/E; bz = PZ/E;
			scale = 60.0/bx; // time to arriave at the target multiplied by c
			
			y = by * scale; // by * c * scale / c = by * scale
			z = bz * scale; 

			if (abs(y) < 0.25 && abs(z) < 0.2)
				hitCounter60m++;

			bx = PX/E; by = PY/E; bz = PZ/E;
			scale = 100.0/bx; // time to arriave at the target multiplied by c
			
			y = by * scale; // by * c * scale / c = by * scale
			z = bz * scale; 

			if (abs(y) < 0.25 && abs(z) < 0.2)
				hitCounter100m++;
		}
		vector<double> ageo; ageo.push_back(mchi);
		ageo.push_back(hitCounter10m);
		ageo.push_back(hitCounter35m);
		ageo.push_back(hitCounter60m);
		ageo.push_back(hitCounter100m);
		records.push_back(ageo);
		cout << i << " th mass point done" << endl;
		pionTxt.clear();
		pionTxt.seekg(0, std::ios::beg); // Set file position to the beginning
	}

	for (const auto &row : records) {
		for (double value : row)
			ageoTxt << value << " ";
		ageoTxt << endl;
	}
	pionTxt.close();
	massTxt.close();
	ageoTxt.close();

	return 0;
}
