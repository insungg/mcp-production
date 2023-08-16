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
    long long seed = strtol(argv[1], NULL, 10);
    long long nCores = strtol(argv[2], NULL, 10);
    long long nEvents = strtol(argv[3], NULL, 10);
    long long nJobs = nEvents / nCores;
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
    TFile* output = new TFile( Form( "/mnt/ntfs1/pythia_rootifles/upsilon_seed%d_t%d.root", (int) seed, (int) ith) , "RECREATE" );
    
    
    // define kinematic variables
    int id;
    double px, py, pz, m;
    
    // create branches
    TTree *tree = new TTree("Upsilon", "Upsilon");
    tree->Branch("id",      &id);
    tree->Branch("px",      &px);
    tree->Branch("py",      &py);
    tree->Branch("pz",      &pz);
    tree->Branch("m" ,       &m);


    // create pythia generator
    Pythia8::Pythia pythia;
    
    // read beam properties
    pythia.readFile("beam.config");
    
    // define random seed
    pythia.readString("Random:setSeed = on"); // use random seed 
    pythia.readString( Form("Random:seed = %d", (int) seed + (int) ith) ); // + ith to prevent redundant event generation 

    pythia.init();

    // event generation  
    for (int i = 0; i < nJobs; i ++) {
        if (!pythia.next()) continue; // skip when generation failed 
        
        // event loop 
        for (int j = 0; j < pythia.event.size(); j++) {		
		if ( abs(pythia.event.at(j).id()) == 553 ) {
			id = pythia.event.at(j).id();
            		px = pythia.event.at(j).px();
            		py = pythia.event.at(j).py();
            		pz = pythia.event.at(j).pz();
            		m  = pythia.event.at(j).m();
            		
			tree->Fill();
		}
	}
    }
      
 
    output->Write();
    output->Close();
}

