using namespace std;
void getAgeo() {
	TFile *input = new TFile("mcp-production.root");
	TTree *tree  = (TTree *) input->Get("mcp");

	double px, py, pz;
	double m = 10.0;
	double r = 30.0;
	tree->SetBranchAddress("Px", &px);
	tree->SetBranchAddress("Py", &py);
	tree->SetBranchAddress("Pz", &pz);

	int counter_vert = 0;
	int counter_hori = 0;

	for (int i = 0; i < tree->GetEntries(); i++) {
		tree->GetEntry(i);
		
		double bx, by, bz;
		double e = sqrt(pow(px,2) + pow(py, 2) + pow(pz, 2) + pow(m, 2));
		bx = px/e; by = py/e; bz = pz/e;

		cout << "bx : " << bx << " e : " << e << endl;

		double scale = r/bx;
		cout << scale << endl;

		double y = by * scale;
		double z = bz * scale;

		cout << "y : " << y << " z : " << z << endl;

		if (abs(y) <= 0.5 && abs(z) <= 0.5)
			counter_vert++;

		if (abs(y) <= 1 && abs(z) <= 1)
			counter_hori++;
	}

	cout << counter_vert << endl;
	cout << counter_hori << endl;

	return 0;
}
