// This Pythia script generates p-p fixed target collision events 
// and collect the momentum of milicharged particles
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <mutex>
#include <thread>
#include <string.h>

#include "TTree.h"
#include "Pythia8/Pythia.h"
#include "Pythia8/PythiaParallel.h"

using namespace std;
using namespace Pythia8;

#define NPARTICLE 4


// usage : submet [seed] [nCores] [nEvents]
int main(int argc, char *argv[])
{
	// check argv is valid
    if (argc < 4) { 
		cout << "REQUIRED FORMAT : submet [seed] [nMass] [nJobs]" << endl; 
        return -1;
    }

	// open mass point
	ifstream masstxt("mass2.txt");
	if (!masstxt.is_open()) {
		cout << "Failed to open mass points" << endl;
		return -1;
	}
	// create output file
	ofstream ageotxt("ageo_numi_high4.txt");
	if (!ageotxt.is_open()) {
		cout << "Failed to create the ageo.txt" << endl;
		return -1;
	}
	ofstream gentxt("gen_numi_high4.txt");
	if (!gentxt.is_open()) {
		cout << "Failed to create the gen.txt" << endl;
		return -1;
	}
	


	// define parameters
    long long  seed = strtol(argv[1], NULL, 10);
    long long  nMass = strtol(argv[2], NULL, 10);
    long long  nJobs = strtol(argv[3], NULL, 10);
	// define records
	vector< vector<double> > ageoRecords;
	vector< vector<double> > genRecords;

	string line;
    for (long long i = 0; i < nMass; i++) {
		// stop when hit the EOF
		if (!getline(masstxt, line)) break;

		long long genCounter[NPARTICLE];
		long long hitCounter[NPARTICLE];
		mutex counterMutex;
		
		for (int j = 0; j < NPARTICLE; j++) {
			genCounter[j] = 1;
			hitCounter[j] = 0;
		}


		// define pythia for multithreading
		PythiaParallel pythia;
		// read beam config and mcp config
		pythia.readFile("beam.config");
		pythia.readFile("momentum.config");
		// set mcp mass
        double mass = atof(line.c_str()); 
		// string defineMass = "31:m0="; defineMass += to_string(mass);	
		// pythia.readString(defineMass);
		pythia.readString(Form("31:m0=%lf", mass)); 
		// set random seed	
		// string defineSeed = "Random:seed="; defineSeed += to_string(seed);
    	pythia.readString("Random:setSeed=on");  
		// pythia.readString(defineSeed);
    	pythia.readString( Form("Random:seed=%d", (int ) (seed)) );  
		// let asyncrhonous process, otherwise the program will slow down
		pythia.readString("Parallelism:processAsync = on");
		pythia.readString("Parallelism:numThreads=100");
		// string defineJobs = "Main:numberOfEvents=";
		// defineJobs += to_string(nJobs); pythia.readString(defineJobs);
		pythia.readString( Form("Main:numberOfEvents = %lld", nJobs) );

		pythia.init();
		pythia.particleData.list(31);

		pythia.run([&] (Pythia *pythiaPtr) {
			int particleId[NPARTICLE] = {111, 221, 443, 553};
			for (long long j = 0; j < pythiaPtr->event.size(); j++) {
				if ( abs(pythiaPtr->event.at(j).id()) == 31) {
					double p  = pythiaPtr->event.at(j).pAbs();
					double pt = pythiaPtr->event.at(j).pT();
					double px = pythiaPtr->event.at(j).px();
					double py = pythiaPtr->event.at(j).py();
					double pz = pythiaPtr->event.at(j).pz();
					double e  = pythiaPtr->event.at(j).e();
					double m  = pythiaPtr->event.at(j).m();

					double x, y;
					double dump = 1040;
					x = px/e * dump; y = py/e * dump;

					int mom1_id = pythiaPtr->event.at(j).mother1();
					int mom1;
					if (mom1_id) mom1 = pythiaPtr->event.at(mom1_id).id(); else mom1 = 0;


					for (int k = 2; k < NPARTICLE - 1; k++) {
						if (mom1 == particleId[k]) {
							std::lock_guard<mutex> lock(counterMutex);
							genCounter[k]++;
							if ( abs(x) < 0.5 && abs(y) < 0.5)
								hitCounter[k]++;
						}
					}
				}
			}
		// The mutex will be released automatically when the lock_gaurd goes out of scope.
		});


		vector<double> gen; gen.push_back(mass);
		for (int j = 0; j < NPARTICLE; j++)
			gen.push_back(genCounter[j]);

		vector<double> ageo; ageo.push_back(mass);
		for (int j = 0; j < NPARTICLE; j++)
			ageo.push_back(hitCounter[j]);

		genRecords.push_back(gen);
		ageoRecords.push_back(ageo);

		// write root file
		// output->Write();
		// output->Close();
		
		cout << i << " th mass is done" << endl;
    }





	// sort by mass before save
	sort(ageoRecords.begin(), ageoRecords.end(), 
			[] (const vector<double> &a, const vector<double> &b) { 
				return a[0] < b[0];
			}
	);
	sort(genRecords.begin(), genRecords.end(), 
			[] (const vector<double> &a, const vector<double> &b) { 
				return a[0] < b[0];
			}
	);
	// write final result to ageo.txt
	for (const auto &row : ageoRecords) {
		for (double value : row) 
			ageotxt << value << " ";
		ageotxt << endl;
	}
	for (const auto &row : genRecords) {
		for (double value : row) 
			gentxt << value << " ";
		gentxt << endl;
	}

	// close file streams
	masstxt.close();
	ageotxt.close();
	gentxt.close();

    return 0;
}
