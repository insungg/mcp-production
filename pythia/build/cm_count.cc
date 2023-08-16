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
    int seed = strtol(argv[1], NULL, 10);
    int nCores = strtol(argv[2], NULL, 10);
    int nEvents = strtol(argv[3], NULL, 10);
    int nJobs = nEvents / nCores;
    cout << "seed : " << seed << "    nCores : " << nCores << "    nJobs : " << nJobs << endl; 

    // define threads for multi-threading 
    TThread* th[nCores];
    int thread_argv[3];
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
    int *parameters = (int *) t_argv; 
    
    int seed = parameters[0];
    int nJobs = parameters[1];
    int ith = parameters[2];

    // define root file
    TFile* output = new TFile( Form( "/data/pythia_rootfiles/cm_count_seed%d_t%d.root", seed, ith) , "RECREATE" );
    
    // define tags
    int nParticle = 7;
    TString particleName[nParticle] = {"Pion", "Eta", "Rho", "Omega", "Phi", "Jpsi", "Upsilon"};
    int particleId[nParticle] = {111, 221, 113, 223, 333, 443, 553};  
    
    // define kinematic variables
    int id;
    double p, pt, pz, m;
    
    // create branches
    TTree *tree[nParticle];
    for (int i = 0; i < nParticle; i++) { 
        tree[i] = new TTree(particleName[i], particleName[i]);
        tree[i]->Branch("id",      &id);
        tree[i]->Branch("p" ,       &p);
        tree[i]->Branch("pt",      &pt);
        tree[i]->Branch("pz",      &pz);
        tree[i]->Branch("m" ,       &m);
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
            		p  = pythia.event.at(j).pAbs();
            		pt = pythia.event.at(j).pT();
            		pz = pythia.event.at(j).pz();
            		m  = pythia.event.at(j).m();
            		
					tree[k]->Fill();
				}
	    	}
        }
    }  
 
    output->Write();
    output->Close();
}

