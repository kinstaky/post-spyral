#include <iostream>
#include <filesystem>
#include <fstream>

#include <TROOT.h>
#include <TCutG.h>


int main(int argc, char **argv) {
	if (argc != 2) {
		std::cout << "Usage: " << argv[0] << " <cut_file>" << std::endl;
		return -1;
	}

	gROOT->Macro(argv[1]);

	std::filesystem::path cut_file(argv[1]);
	std::filesystem::path cut_data_file = cut_file.replace_extension(".txt");

	TCutG *cut = (TCutG *)gROOT->FindObject("CUTG");
	if (!cut) {
		std::cout << "Error: cut not found" << std::endl;
		return -1;
	}

	std::ofstream fout(cut_data_file);
	double *x = cut->GetX();
	double *y = cut->GetY();
	fout << cut->GetN() << "\n";
	for (int i = 0; i < cut->GetN(); ++i) {
		fout << x[i] << " " << y[i] << "\n";
	}
	fout.close();
	return 0;
}