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
#define PARCUT 2 


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
	ofstream ageo_numi_txt("ageo_numi_h3.txt");
	if (!ageo_numi_txt.is_open()) {
		cout << "Failed to create the ageo_numi.txt" << endl;
		return -1;
	}
	ofstream ageo_dune_txt("ageo_dune_h3.txt");
	if (!ageo_dune_txt.is_open()) {
		cout << "Failed to create the ageo_dune.txt" << endl;
		return -1;
	}
	ofstream gentxt("gen_h3.txt");
	if (!gentxt.is_open()) {
		cout << "Failed to create the gen.txt" << endl;
		return -1;
	}
	


	// define parameters
    long long  seed = strtol(argv[1], NULL, 10);
    long long  nMass = strtol(argv[2], NULL, 10);
    long long  nJobs = strtol(argv[3], NULL, 10);
	// define records
	vector< vector<double> > ageoRecords_numi;
	vector< vector<double> > ageoRecords_dune;
	vector< vector<double> > genRecords;

	string line;
    for (long long i = 0; i < nMass; i++) {
		// stop when hit the EOF
		if (!getline(masstxt, line)) break;

		long long genCounter[NPARTICLE];
		long long hitCounter_numi[NPARTICLE];
		long long hitCounter_dune[NPARTICLE];
		mutex counterMutex;
		
		for (int j = 0; j < NPARTICLE; j++) {
			genCounter[j] = 1;
			hitCounter_numi[j] = 0;
			hitCounter_dune[j] = 0;
		}

		// define pythia for multithreading
		PythiaParallel pythia;
		// read beam config and mcp config
		pythia.readFile("beam.config");
		pythia.readFile("momentum.config");
		// set mcp mass
        double mass = atof(line.c_str()); 
		pythia.readString(Form("31:m0=%lf", mass)); 
    	pythia.readString("Random:setSeed=on");  
    	pythia.readString( Form("Random:seed=%d", (int ) (seed)) );  
		// let asyncrhonous process, otherwise the program will slow down
		pythia.readString("Parallelism:processAsync = on");
		pythia.readString("Parallelism:numThreads=100");
		pythia.readString( Form("Main:numberOfEvents = %lld", nJobs) );

		pythia.init();
		pythia.particleData.list(31);

		pythia.run([&] (Pythia *pythiaPtr) {
			int particleId[NPARTICLE] = {111, 221, 443, 553};
			for (long long j = 0; j < pythiaPtr->event.size(); j++) {
				if ( abs(pythiaPtr->event.at(j).id()) == 31) {
					double pt = pythiaPtr->event.at(j).pT();
					double px = pythiaPtr->event.at(j).px();
					double py = pythiaPtr->event.at(j).py();
					double pz = pythiaPtr->event.at(j).pz();
					double e  = pythiaPtr->event.at(j).e();
					double m  = pythiaPtr->event.at(j).m();

					double scale_x = px/pz;
					double scale_y = py/pz;

					double x_numi, y_numi;
					double x_dune, y_dune;
					double dump_numi = 1040;
					double dump_dune = 574;
					x_numi = scale_x * dump_numi; y_numi = scale_y * dump_numi;
					x_dune = scale_x * dump_dune; y_dune = scale_y * dump_dune;

					int mom1_id = pythiaPtr->event.at(j).mother1();
					int mom1;
					if (mom1_id) mom1 = pythiaPtr->event.at(mom1_id).id(); else mom1 = 0;


					for (int k = PARCUT; k < NPARTICLE; k++) {
						if (mom1 == particleId[k]) {
							std::lock_guard<mutex> lock(counterMutex);
							genCounter[k]++;
							if ( abs(x_numi) < 0.5 && abs(y_numi) < 0.5)
								hitCounter_numi[k]++;
							if ( abs(x_dune) < 0.5 && abs(y_dune) < 0.5)
								hitCounter_dune[k]++;
						}
					}
				}
			}
		// The mutex will be released automatically when the lock_gaurd goes out of scope.
		});


		vector<double> gen; gen.push_back(mass);
		for (int j = PARCUT; j < NPARTICLE; j++)
			gen.push_back(genCounter[j]);

		vector<double> ageo_numi; ageo_numi.push_back(mass);
		for (int j = PARCUT; j < NPARTICLE; j++)
			ageo_numi.push_back(hitCounter_numi[j] * 1.0 / genCounter[j]);

		vector<double> ageo_dune; ageo_dune.push_back(mass);
		for (int j = PARCUT; j < NPARTICLE; j++)
			ageo_dune.push_back(hitCounter_dune[j] * 1.0 / genCounter[j]);

		genRecords.push_back(gen);
		ageoRecords_numi.push_back(ageo_numi);
		ageoRecords_dune.push_back(ageo_dune);
		
		cout << i << " th mass is done" << endl;
    }





	// sort by mass before save
	sort(ageoRecords_numi.begin(), ageoRecords_numi.end(), 
			[] (const vector<double> &a, const vector<double> &b) { 
				return a[0] < b[0];
			}
	);
	sort(ageoRecords_dune.begin(), ageoRecords_dune.end(), 
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
	for (const auto &row : ageoRecords_numi) {
		for (double value : row) 
			ageo_numi_txt << value << " ";
		ageo_numi_txt << endl;
	}
	for (const auto &row : ageoRecords_dune) {
		for (double value : row) 
			ageo_dune_txt << value << " ";
		ageo_dune_txt << endl;
	}
	for (const auto &row : genRecords) {
		for (double value : row) 
			gentxt << value << " ";
		gentxt << endl;
	}

	// close file streams
	masstxt.close();
	ageo_numi_txt.close();
	ageo_dune_txt.close();
	gentxt.close();

    return 0;
}
