#include <iostream>
#include <cmath>

#include "TRandom3.h"
#include "TH2D.h"
#include "TCanvas.h"

using namespace std;

// constants
const double mpi = 134.9768; // MeV, neutral pion mass
const double mp  = 938.272;
const double PI = 3.14159265358979;
const double conv2rad = PI/180.0;
// NOTE : The reference paper uses degree, not radian

// beam property
const int Z = 6; // target proton number, Z = 6 for carbon target
const double Tp = 800.0; // energy of proton beam, 800 MeV for LANSCE
const double TpiMax = Tp - mpi; // maximum kinetic energy of pion produced


// distrbution parameters
const int TRIALS = 10000000;
const double knots[]  = {0, 0, 0, 30, 70, 180, 180, 180};
const double B = 25.0;

double envelope = 10; // enveloping function of dSigma used for rejection sampling
double TA, sigmaA, NormZ;
double bsplineCoeff[5];


// set parpameters of distribution 
void initialize();
// define Bspline function
double Bspline(int idx, int order, double theta);
// Burman-Smith distribution, theta in radian and Tpi in MeV
double ddSigma(double Tpi, double theta);


// main program
void generatePion() {
	ofstream output("Pion_Production.txt");
	if (!output.is_open()) {
		cout << "Failed to create output file" << endl;
		return -1;
	}

	// initialize
	initialize();
	TRandom3 rand;
	vector<vector<double>> records;
	TH2D *hist = new TH2D("hist", "hist", 1000, 0, 1000, 100, -1, 1);
	
	// sample p (magnitude of 3-momentum of pion), theta, phi
	for (int i = 0; i < TRIALS; i++) {
		vector<double> momentum;
		// rejection sampling
		while (true) {
			double Tpi    = rand.Uniform(0.0, TpiMax);
			double theta  = rand.Uniform(0.0, 180.0);

			double y = ddSigma(Tpi, theta);
			// accept when y/envelope > u where u ~ Unif(0, 1)
			if ( y >= rand.Uniform(0, 1) * envelope) {
				// renew envelope when y > envelope
				if ( y > envelope )
					envelope = y;

				double P   = sqrt( pow(Tpi + mpi, 2) - pow(mpi, 2) ); // pion momentum in MeV/c
				double phi = rand.Uniform(0.0, 2.0*PI);
				theta = theta * conv2rad;

				momentum.push_back(P);
				momentum.push_back(cos(theta));
				momentum.push_back(phi);

				hist->Fill(P, cos(theta));
				// cout << "Y : " << y << " in " << i << endl;
				break;
			}
		}
		records.push_back(momentum);
		if (i % 10000 == 0) cout << i << " th job done" << endl;
	}

	// save momentum distribution
	for (const auto &row : records) {
		for (double value : row) 
			output << value << " ";
		output << endl;
	}
	output.close();


	// draw histograms
	TCanvas *c1 = new TCanvas(); c1->cd();
	hist->GetXaxis()->SetTitle("Momentum [MeV/c]");
	hist->GetYaxis()->SetTitle("Cos(theta)");
	hist->Draw("colz");

	TCanvas *c2 = new TCanvas(); c2->cd();
	// TH1D *projX = (TH1D *) hist->ProjectX("projX");



	return 0;
}


void initialize() {
	double TA585, TA730, sigmaA585, sigmaA730;
	if(Z == 1){
        TA730 = 53.0;
        TA585 = 64.5;
        sigmaA730 = 127.0;
        sigmaA585 = 155.0;
    }
    else if(Z <= 8){
        TA730 = 34.2;
        TA585 = 28.9;
        sigmaA730 = 150.0;
        sigmaA585 = 130.0;
    }
    else if(Z < 92){
        TA730 = 29.9;
        TA585 = 26.0;
        sigmaA730 = 166.0;
        sigmaA585 = 135.0;
    }
	else {
		cout << "Z should be less than 93" << endl;
		return -1;
	}
	TA 	  = ( TA730 * (Tp - 585) - TA585 * (Tp - 730) ) / (730 - 585);
	sigmaA = ( sigmaA730 * (Tp - 585) - sigmaA585 * (Tp-730) ) / (730 - 585);
	

	NormZ = 0;
	double normzc[] = {0.8851, -0.1015, 0.1459, -0.0265};
	for (int i = 0; i < 4; i++) 
		NormZ += normzc[i] * pow(log(Z), i) * pow(Z, 0.3333);

	bsplineCoeff[0] = std::min( 27.0 - 4.0 * pow( (730.0-Tp)/(730.0-585.0), 2 ), 27.0 );
	bsplineCoeff[1] = 18.2; 
	bsplineCoeff[2] = 8.0;
	bsplineCoeff[3] = 13.0 + (Z - 12.0) / 10.0;
	bsplineCoeff[4] = 9.0 + (Z-12.0)/10.0 - (Tp - 685.0)/20.0;
	
}


// see Cox-deBoor's algorithm which recursively defines Bspline
double Bspline(int idx, int order, double theta) {
	if (order == 0) {
		if ( theta >= knots[idx] && theta <= knots[idx+1]) return 1.0;
		else return 0.0;
	}

	double total = 0;
	if (knots[idx + order] != knots[idx]) 
		total += (theta - knots[idx]) / (knots[idx + order] - knots[idx]) * Bspline(idx, order - 1, theta);
	if (knots[idx + order + 1] != knots[idx + 1])
		total += (knots[idx + order + 1] - theta) / (knots[idx + order + 1] - knots[idx + 1]) * Bspline(idx + 1, order - 1, theta);

	return total;
}


double ddSigma(double Tpi, double theta) { // note : Tpi in MeV, theta in degree
	double Tbar  = 48 + 330 * exp(-theta/TA);
	double sigma = sigmaA * exp(-theta / 85.0); 
	double TF;
	if (Z == 1) {
		// hydrogen is not implemented yet
	}
	else
		TF = Tp - 140 - 2*B;

	double amp = 0;
	for (int i = 0; i < 5; i++) {
		amp += bsplineCoeff[i] * Bspline(i, 2, theta);
	}
	amp = amp * NormZ;


	// d/dtheta = d/dcos X dcos/dtheta = sin d/dcos 
	return sin(theta*conv2rad) * amp * exp(-1.0 * pow( (Tbar - Tpi) / sqrt(2) / sigma, 2) ) / (1.0 + exp((Tpi-TF) / B)); 
}

