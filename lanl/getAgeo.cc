using namespace std;
void getAgeo() {
	TFile *input = new TFile("mcp-production.root");
	TTree *tree  = (TTree *) input->Get("mcp");

	double px, py, pz;
	double m = 10.0;
	double r = 10.0;
	tree->SetBranchAddress("Px", &px);
	tree->SetBranchAddress("Py", &py);
	tree->SetBranchAddress("Pz", &pz);

	int counter_10m = 0;
	int counter_35m = 0;
	int counter_60m = 0;
	int counter_100m = 0;

	for (int i = 0; i < tree->GetEntries(); i++) {
		tree->GetEntry(i);
		
		double bx, by, bz;
		double e = sqrt(pow(px,2) + pow(py, 2) + pow(pz, 2) + pow(m, 2));
		bx = px/e; by = py/e; bz = pz/e;

		// cout << "bx : " << bx << " e : " << e << endl;

		double scale = r/bx;
		// cout << scale << endl;

		double y = by * scale;
		double z = bz * scale;

		// cout << "y : " << y << " z : " << z << endl;

		if (abs(y) <= 0.25 && abs(z) <= 0.2)
			counter_10m++;

        scale = 35.0/bx;
        y = by * scale;
        z = bz * scale;
		if (abs(y) <= 0.25 && abs(z) <= 0.2)
			counter_35m++;

        scale = 60.0/bx;
        y = by * scale;
        z = bz * scale;
		if (abs(y) <= 0.25 && abs(z) <= 0.2)
			counter_60m++;

        scale = 100.0/bx;
        y = by * scale;
        z = bz * scale;
		if (abs(y) <= 0.25 && abs(z) <= 0.2)
			counter_100m++;
	}

	cout << counter_10m << endl;
	cout << counter_35m << endl;
	cout << counter_60m << endl;
	cout << counter_100m << endl;

	return 0;
}
