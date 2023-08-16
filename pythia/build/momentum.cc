// This Pythia script generates p-p fixed target collision events 
// and collect the momentum of milicharged particles
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
		sleep(1); // wait until initialization finish
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
    TFile* output = new TFile( Form( "/mnt/ntfs1/pythia_rootifles/momentum_seed%d_t%d.root", seed, ith) , "RECREATE" );
    
    // define kinematic variables
    int id;
    double px, py, pz, pt, p, m;
    double angle;

    double dump_length = 110;
    
    // create branches
    TTree *tree = new TTree( "mCP" , "mCP");
    tree->Branch("id", &id);
    tree->Branch("pz", &pz);
    tree->Branch("pt", &pt);
    tree->Branch("p",   &p);
    tree->Branch("m",   &m);
    tree->Branch("angle", &angle);
    

    TH2D *board = new TH2D();

    // create pythia generator
    Pythia8::Pythia pythia;
    
    // read beam properties
    pythia.readFile("beam.config");
    pythia.readFile("momentum.config");
    
    // define random seed
    pythia.readString("Random:setSeed = on"); // use random seed 
    pythia.readString(Form("Random:seed = %d", seed + ith)); // + ith to prevent redundant event generation 

    pythia.init();


    // pythia.particleData.list(31); // check if mcp is defined

    // event generation  
    for (int i = 0; i < nJobs; i ++) {
        if (!pythia.next()) continue; // skip when generation failed 
        
	for (int j = 0; j < pythia.event.size(); j++) {
		if ( abs(pythia.event.at(j).id()) == 31 ) {
			id = pythia.event.at(j).id();
        		px = pythia.event.at(j).px();
        		py = pythia.event.at(j).py();
        		pz = pythia.event.at(j).pz();
			
			p = sqrt( pow(px,2) + pow(py,2) + pow(pz,2)  );
			pt = sqrt( pow(px,2) + pow(py,2)  );

        		m  = pythia.event.at(j).m();
        	
			tree->Fill();
		}  

	}
    }  
 
    output->Write();
    output->Close();
}

