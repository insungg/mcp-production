// This Pythia script generates p-p fixed target collision events 
// and collect the momentum of milicharged particles
#include <iostream>
#include <stdlib.h>
#include <fstream>

using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TThread.h"
#include "Pythia8/Pythia.h"


// define ageo vector for record
vector< vector<double> > ageoRecords;
vector< vector<double> > genRecords;

// multi-threading callback function
void* handler(void *t_argv); 

// usage : submet [seed] [nCores] [nEvents]
int main(int argc, char *argv[])
{
	// check argv is valid
    if (argc < 4) { 
		cout << "REQUIRED FORMAT : submet [seed] [nCores] [nJobs]" << endl; 
        return -1;
    }

	// open mass point
	ifstream masstxt("mass.txt");
	if (!masstxt.is_open()) {
		cout << "Failed to open mass points" << endl;
		return -1;
	}
	// create output file
	ofstream ageotxt("ageo.txt");
	if (!ageotxt.is_open()) {
		cout << "Failed to create the ageo.txt" << endl;
		return -1;
	}
	ofstream gentxt("gen.txt");
	if (!gentxt.is_open()) {
		cout << "Failed to create the ageo.txt" << endl;
		return -1;
	}

	// define thread parameters
    long long  seed = strtol(argv[1], NULL, 10);
    long long  nCores = strtol(argv[2], NULL, 10);
    long long  nJobs = strtol(argv[3], NULL, 10);
	cout << "seed : " << seed << "    nCores : " << nCores << "    nJobs : " << nJobs << endl; 

    // define threads for multi-threading 
    TThread* th[nCores];
    Double_t thread_argv[3];
    thread_argv[0] = seed;
    thread_argv[1] = nJobs;

    // throw jobs to threads
	string line;
    for (int i = 0; i < nCores; i++) {
		// stop when hit the EOF
		if (!getline(masstxt, line)) break;

		// define thread index
		thread_argv[2] = i;
		// define mcp mass
        thread_argv[3] = atof(line.c_str());

		cout << thread_argv[3] << endl;

		// call handler
        th[i] = new TThread(Form("th%d", i), handler, (void *) thread_argv);
        th[i]->Run();

		// wait until callbalck initialization finish
		sleep(1.5); 
    }

    // join works when each josb is finished
    for(int i = 0; i < nCores; i++) { 
        th[i]->Join();
		delete th[i];
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



void* handler(void *t_argv) 
{
    Double_t *parameters = (Double_t *) t_argv; 

    Int_t seed  = parameters[0];
    Int_t nJobs = parameters[1];
	Int_t ith   = parameters[2];

    Double_t mass  = parameters[3];
	
	cout << seed << " " << nJobs << " " << ith << " " << mass << endl;

    // define root file
    TFile* output = new TFile( 
			Form( "/data3/ishwang/pythia_rootfiles/ageo_seed%d_mass%.3f_%dth.root", seed, mass, ith) , "RECREATE" );
    
    // define kinematic variables
    Int_t id, status;
    Double_t px, py, pz, pt, p, m, e;
    Double_t bx, by;
	Double_t  x,  y;
	Int_t mom1, mom2, dau1, dau2;

    Double_t dump_length = 1040; // in metre
	Double_t c = 299792458; // speed of light in m/s
    Double_t flight_time = dump_length / c;


	// define tags
	Int_t   nParticle = 4;
	TString particleName[nParticle] = {"pi", "eta", "jpsi", "upsilon"};
	Int_t   particleId[nParticle]   = {111, 221, 443, 553}; 


    // create branches for mother mesons
    TTree *tree[nParticle];
	for (int i = 0; i < nParticle; i++) {
		tree[i] = new TTree(particleName[i], particleName[i]);
    	tree[i]->Branch("id", &id);
    	tree[i]->Branch("px", &px);
    	tree[i]->Branch("py", &py);
    	tree[i]->Branch("pz", &pz);
    	tree[i]->Branch("pt", &pt);
    	tree[i]->Branch("p",   &p);
    	tree[i]->Branch("x",   &x);
    	tree[i]->Branch("y",   &y);
    	tree[i]->Branch("m",   &m);
		tree[i]->Branch("mom1", &mom1);
		tree[i]->Branch("mom2", &mom2);
		tree[i]->Branch("dau1", &dau1);
		tree[i]->Branch("dau2", &dau2);
		tree[i]->Branch("status", &status);
	}
	
	/*
	// create branches for total mcp contribution
	TTree *mcp = new TTree("mcp", "mcp");
    mcp->Branch("id", &id);
   	mcp->Branch("px", &px);
   	mcp->Branch("py", &py);
   	mcp->Branch("pz", &pz);
   	mcp->Branch("pt", &pt);
   	mcp->Branch("p",   &p);
   	mcp->Branch("x",   &x);
   	mcp->Branch("y",   &y);
   	mcp->Branch("m",   &m);
	mcp->Branch("mom1", &mom1);
	mcp->Branch("mom2", &mom2);
	mcp->Branch("dau1", &dau1);
	mcp->Branch("dau2", &dau2);
	mcp->Branch("status", &status);
	*/ 

	cout << "root was successfully initialized " << ith << " threads" << endl;



    // create pythia generator
    Pythia8::Pythia pythia;
    
    // read beam properties
    pythia.readFile("beam.config");
    pythia.readFile("momentum.config");
	pythia.readString(Form("%s%lf", "31:m0=", mass));
    
    // define random seed
    pythia.readString("Random:setSeed=on");  
    pythia.readString( Form("Random:seed=%d", (int ) (seed) ) );  

	// init pythia
    pythia.init();
	cout << "pythia was successfully initialized " << ith << " threads" << endl;

	// check if mcp is defined
    pythia.particleData.list(31); 

	// define hit counter
	Int_t genCount[nParticle];
	Int_t hitCount[nParticle];
	for (Int_t i = 0; i < nParticle; i++) {
		genCount[i] = 1;
		hitCount[i] = 0;
	}

    // event generation  
    for (long long i = 0; i < nJobs; i ++) {
        if (!pythia.next()) continue; // skip when generation failed 
        
		// loop through event
		for (long long j = 0; j < pythia.event.size(); j++) {
			// detect mcp
			if ( abs(pythia.event.at(j).id()) == 31 ) {
				id = pythia.event.at(j).id();
    		    px = pythia.event.at(j).px();
    		    py = pythia.event.at(j).py();
    		    pz = pythia.event.at(j).pz();
				pt = pythia.event.at(j).pT();
				p  = pythia.event.at(j).pAbs();
    	    	m  = pythia.event.at(j).m();
                e  = pythia.event.at(j).e();

				status = pythia.event.at(j).status();

                Int_t mom1_id = pythia.event.at(j).mother1();
                Int_t mom2_id = pythia.event.at(j).mother2();
				Int_t dau1_id = pythia.event.at(j).daughter1();
				Int_t dau2_id = pythia.event.at(j).daughter2();

				if (mom1_id) mom1 = pythia.event.at(mom1_id).id(); else mom1 = 0;
				if (mom2_id) mom2 = pythia.event.at(mom2_id).id(); else mom2 = 0;
				if (dau1_id) dau1 = pythia.event.at(dau1_id).id(); else dau1 = 0;
				if (dau2_id) dau2 = pythia.event.at(dau2_id).id(); else dau2 = 0;


                bx = px/e; by = py/e;

                x  = bx * dump_length;
                y  = by * dump_length;
				
				// mcp->Fill();

				// check if mcp hit target
				for (int k = 0; k < nParticle; k++) {
					if (particleId[k] == mom1) {
						genCount[k]++;
						tree[k]->Fill();
						if (abs(x) < 0.5 && abs(y) < 0.5) // detector face = 1m x 1m
							hitCount[k]++;
					}
				}
			}  
		}
    }  

	// define hit ratio (=ageo)
	vector<double>  ageo; ageo.push_back(mass);
	for (int i = 0; i < nParticle; i++)
		ageo.push_back( hitCount[i] * 1.0 / genCount[i] ); // convert int to double
	vector<double>  gen; gen.push_back(mass);
	for (int i = 0; i < nParticle; i++)
		gen.push_back( genCount[i] );
	
	// store final data at ith row
	// ageoRecords.insert(ageoRecords.begin() + ith, ageo);
	ageoRecords.push_back(ageo);
	genRecords.push_back(gen);


    output->Write();
    output->Close();

	cout << ith << " core  was successfully finished job" << endl;
}

