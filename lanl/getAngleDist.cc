using namespace std;
void getAngleDist() {
	TFile *input = new TFile("mcp-production.root");
	TTree *tree  = (TTree *) input->Get("mcp");

	double px, py, pz;
	double m = 10.0;
	double r = 30.0;
	tree->SetBranchAddress("Px", &px);
	tree->SetBranchAddress("Py", &py);
	tree->SetBranchAddress("Pz", &pz);

    TH2D* hist = new TH2D("hist", "hist", 100, 0, TMath::Pi(), 100, -1, 1);

	for (int i = 0; i < tree->GetEntries(); i++) {
		tree->GetEntry(i);

        TVector3 v1(px, py, pz);
		
        // double p = sqrt(px*px + py+py + pz + pz);
        double phi = v1.Phi();
        double cos = v1.CosTheta();

        hist->Fill(phi, cos);
	}


    hist->GetXaxis()->SetTitle("phi");
    hist->GetYaxis()->SetTitle("cos(theta)");
    hist->SetStats(0);
    hist->Draw("colz");
    

	return 0;
}
