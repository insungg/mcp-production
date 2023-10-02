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
const double mchi      = 10.0; // MeV, mcp mass
const double BrPi2gg   = 1e9; // Br(pion -> gamam gamma) 


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


void decayPion() {
	// initialize
	ifstream input("Pion_Production.txt");
	if (!input.is_open()) {
		cout << "Failed to open pion data" << endl;
		return -1;
	}
	TFile *output = new TFile("mcp-production.root", "RECREATE");
	TTree *tree = new TTree("mcp", "mcp");
	double PX, PY, PZ, PP, M; 
	tree->Branch("Px", &PX);
	tree->Branch("Py", &PY);
	tree->Branch("Pz", &PZ);
	tree->Branch("PP", &PP);
	tree->Branch("M",  &M );
	TRandom3 rand;



	// loop through pion data 
	double momP, momTheta, momPhi;
	while (input >> momP >> momTheta >> momPhi) {
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
		double P 	   = lambda(sqrt(s), mchi, mchi);
		double phi = rand.Uniform(0, 2.0*PI);
		Particle mcp1; 
		mcp1.momentum.SetPxPyPzE(P*sin(theta)*cos(phi), P*sin(theta)*sin(phi), P*cos(theta), sqrt(P*P + mchi*mchi));
		// cout << mcp1.momentum.Px() << " " << mcp1.momentum.Py() << " " << mcp1.momentum.Pz() << endl;
		Particle mcp2;
		mcp2.momentum.SetPxPyPzE(-P*sin(theta)*cos(phi), -P*sin(theta)*sin(phi), -P*cos(theta), sqrt(P*P + mchi*mchi));

		// boost V from its rest frame along z direction
		double vP    = lambda(mpi, sqrt(s), 0); // pion -> V + gamma, mass of gamma = 0
		double vBeta = vP / sqrt(vP*vP + s); 
		// cout << "vBeta : " << vBeta << endl;
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
		tree->Fill();

		PX = mcp2.momentum.Px();
		PY = mcp2.momentum.Py();
		PZ = mcp2.momentum.Pz();
		PP = sqrt(pow(PX, 2) + pow(PY, 2) + pow(PZ, 2));
		M  = mcp2.momentum.M();
		tree->Fill();
		
	}
	
	
	output->Write();
	output->Close();

	return 0;
}
