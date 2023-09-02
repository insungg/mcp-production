#include <fstream>
#include <TTree>

int main() {
	ofstreami output("autoMad.txt");
	if (!output.is_open()) {
		cout << "Failed to make output file" << endl; 
		return -1;
	}

	



}
