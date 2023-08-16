// This Pythia script generates p-p fixed target collision events 
// and collect the momentum of milicharged particles
#include <iostream>
#include <stdlib.h>

#include "TFile.h"
#include "TTree.h"
#include "TThread.h"
#include "Pythia8/Pythia.h"

// multi-threading callback function
void* handler(void *t_argv); 


// usage : submet [seed] [nCores] [nEvents]
int main(int argc, char *argv[])
{
    if (argc < 4) { 
		std::cout << "REQUIRED FORMAT : submet [seed] [nCores] [nEvents]" << std::endl; 
        return -1;
    }

    long long  seed = strtol(argv[1], NULL, 10);
    long long  nCores = strtol(argv[2], NULL, 10);
    long long  nEvents = strtol(argv[3], NULL, 10);
    long long  nJobs = nEvents / nCores;
	std::cout << "seed : " << seed << "    nCores : " << nCores << "    nJobs : " << nJobs << std::endl; 

    
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
		sleep(1); // wait until callbalck initialization finish
    }

    // join works when each josbs are finished
    for(int i = 0; i < nCores; i++) 
        th[i]->Join();

    return 0;
}



void* handler(void *t_argv) 
{
    long long  *parameters = (long long  *) t_argv; 
    
    long long seed = parameters[0];
    long long nJobs = parameters[1];
    long long ith = parameters[2];

    // define root file
    TFile* output = new TFile( Form( "/data/pythia_rootfiles/jpsi_seed%d_t%d.root", seed, ith) , "RECREATE" );
    
    // define kinematic variables
    int id, status;
    double px, py, pz, pt, p, m, e;
    double bx, by;
	double x, y;
	int mom1, mom2, daught1, daught2;

    double dump_length = 1040; // in metre
	double c = 299792458; // speed of light in m/s
    double flight_time = dump_length / c;

	
	int nFinals;
    
    // create branches
    TTree *tree = new TTree( "mCP" , "mCP");
    tree->Branch("id", &id);
    tree->Branch("status", &status);
    tree->Branch("px", &px);
    tree->Branch("py", &py);
    tree->Branch("pz", &pz);
    tree->Branch("pt", &pt);
    tree->Branch("p",   &p);
    tree->Branch("x",   &x);
    tree->Branch("y",   &y);
    tree->Branch("m",   &m);
	tree->Branch("mom1", &mom1);
	tree->Branch("mom2", &mom2);
	tree->Branch("daughter1", &daught1);
	tree->Branch("daughter2", &daught2);

	TTree *jpsi = new TTree( "jpsi", "jpsi");
	jpsi->Branch("nFinals", &nFinals);
	jpsi->Branch("id", &id);
	jpsi->Branch("status", &status);
	jpsi->Branch("mom1", &mom1);
	jpsi->Branch("mom2", &mom2);
	jpsi->Branch("daughter1", &daught1);
	jpsi->Branch("daughter2", &daught2);
	
    

	std::cout << "root was successfully initialized " << ith << " threads" << std::endl;



    // create pythia generator
    Pythia8::Pythia pythia;
    
    // read beam properties
    pythia.readFile("beam.config");
    pythia.readFile("momentum.config");
    
    // define random seed
    pythia.readString("Random:setSeed=on"); // use random seed 
    pythia.readString( Form("Random:seed=%d", (int ) (seed + ith)) ); // + ith to prevent redundant event generation 

    pythia.init();

	std::cout << "pythia was successfully initialized " << ith << " threads" << std::endl;


    pythia.particleData.list(31); // check if mcp is defined

    // event generation  
    for (long long i = 0; i < nJobs; i ++) {
        if (!pythia.next()) continue; // skip when generation failed 
        
		bool flag = false;
		for (long long j = 0; j < pythia.event.size(); j++) {
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
                mom1 = pythia.event.at(j).mother1();
                mom2 = pythia.event.at(j).mother2();
				daught1 = pythia.event.at(j).daughter1();
				daught2 = pythia.event.at(j).daughter2();

				if ( abs(mom1) == 443) flag = true;

                bx = px/e; by = py/e;

                x  = bx * dump_length;
                y  = by * dump_length;

				tree->Fill();
			}  
			if ( abs(pythia.event.at(j).id()) == 443) flag = true;

		}

		if (flag) {
			for (long long j = 0; j < pythia.event.size(); j++) {
				id = pythia.event.at(j).id();
				status = pythia.event.at(j).status();
                mom1 = pythia.event.at(j).mother1();
                mom2 = pythia.event.at(j).mother2();
				daught1 = pythia.event.at(j).daughter1();
				daught2 = pythia.event.at(j).daughter2();

				jpsi->Fill();
			}
			nFinals = pythia.event.size();
			jpsi->Fill();
		}

		if (i % 10000000 == 0) std::cout << ith << " finished " << i << " jobs" << std::endl;
    }  

    output->Write();
    output->Close();

	std::cout << ith << " core  was successfully finished job" << std::endl;
}

