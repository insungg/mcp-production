// This Pythia script generates p-p fixed target collision events 
// and counts number of specified mesons to meson production per P.O.T
#include <iostream>
#include <stdlib.h>

#include "TFile.h"
#include "TTree.h"
#include "TThread.h"
#include "Pythia8/Pythia.h"

using namespace std;


// multi-threading callback function
void* handler(void *t_argv); 


// usage : submet [seed] [nCores] [nEvents]
int main(int argc, char *argv[])
{
    if (argc < 4) { 
        cout << "REQUIRED FORMAT : submet [seed] [nCores] [nEvents]" << endl; 
        return -1;
    }

    // define parameters
    long long  seed = strtol(argv[1], NULL, 10);
    long long  nCores = strtol(argv[2], NULL, 10);
    long long  nEvents = strtol(argv[3], NULL, 10);
    long long  nJobs = nEvents / nCores;
    cout << "seed : " << seed << "    nCores : " << nCores << "    nJobs : " << nJobs << endl; 

    // define threads for multi-threading 
    TThread* th[nCores];
    long long thread_argv[3];
    thread_argv[0] = seed;
    thread_argv[1] = nJobs;

    // throw jobs to threads
    for (int i = 0; i < nCores; i++) {
        thread_argv[2] = i;
        th[i] = new TThread(Form("th%d", i), handler, (void *) thread_argv);
        th[i]->Run();
		sleep(1); // wait until callbalck init finish
    }

    // join works when each josbs are finished
    for(int i = 0; i < nCores; i++) 
        th[i]->Join();

    return 0;
}



void* handler(void *t_argv) 
{
     long long *parameters = (long long *) t_argv; 
    
     long long seed = parameters[0];
     long long nJobs = parameters[1];
     long long ith = parameters[2];

    // define root file
    TFile* output = new TFile( Form( "/data3/ishwang/pythia_rootfiles/cm_fast_seed%d_t%d.root", seed, ith) , "RECREATE" );
    
    // define tags
    int nParticle = 13;
    TString particleName[nParticle] = {"Jpsi", "Upsilon(1s)", "Upsilon(11D)","Upsilon(2s)","Upsilon(12D)","Upsilon(3s)","Upsilon(4s)", "Upsilon(10860)","Upsilon(11020)", "Upsilon(22D)", "Upsilon(22D)", "Upsilon(31D)", "Upsilon(32D)"};
    int particleId[nParticle] = {443, 553, 30553, 100553, 130553, 200553, 300553, 9000553, 9010553, 20555, 120555, 557, 100557};
    
    // define kinematic variables
    int id;
    double px, py, pz, e, m;
	double dump_length = 1040;
	double bx, by, x, y;
    
    // create branches
    TTree *tree[nParticle];
    for (int i = 0; i < nParticle; i++) { 
        tree[i] = new TTree(particleName[i], particleName[i]);
        tree[i]->Branch("id",      &id);
        tree[i]->Branch("px",      &px);
        tree[i]->Branch("py",      &py);
        tree[i]->Branch("pz",      &pz);
        tree[i]->Branch("m" ,       &m);
		tree[i]->Branch("x" ,       &x);
		tree[i]->Branch("y" ,       &y);
    }


    // create pythia generator
    Pythia8::Pythia pythia;
    
    // read beam properties
    pythia.readFile("beam.config");
    
    // define random seed
    pythia.readString("Random:setSeed = on"); // use random seed 
    pythia.readString( Form("Random:seed = %d", seed + ith) ); // + ith to prevent redundant event generation 

    pythia.init();

    // event generation  
    for (int i = 0; i < nJobs; i ++) {
        if (!pythia.next()) continue; // skip when generation failed 
        
        // event loop 
        for (int j = 0; j < pythia.event.size(); j++) {		
	    	for (int k = 0; k < nParticle; k++) {
				if ( abs(pythia.event.at(j).id()) == particleId[k] ) {
					id = pythia.event.at(j).id();
            		px = pythia.event.at(j).px();
            		py = pythia.event.at(j).py();
            		pz = pythia.event.at(j).pz();
            		m  = pythia.event.at(j).m();
					
					e  =  pythia.event.at(j).e();
					bx = px/e; by = py/e;
					x  = bx * dump_length;
					y  = by * dump_length;
            		
					tree[k]->Fill();
				}
	    	}
        }
    }  
 
    output->Write();
    output->Close();
}

