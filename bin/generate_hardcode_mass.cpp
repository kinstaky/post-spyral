#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "include/nuclear_data.h"
#include "include/path_manager.h"

constexpr double electron_mass = 0.510998928;	// MeV/c^2
constexpr double atomic_mass_unit = 931.4940954; // MeV/c^2

void GenerateHardcodeData() {
	// get path
	std::filesystem::path exec_path = gong::PathManager::Instance().ExecPath();
	std::filesystem::path amdc_data_path = exec_path / "../assets/amdc_2020.txt";
	std::filesystem::path amdc_ion_data_path = exec_path / "../assets/amdc_ion_2020.txt";
	std::filesystem::path hardcode_data_path = exec_path / "../assets/hardcode.txt";


	// open file
	std::ifstream fin(amdc_data_path);
	if (!fin.good()) {
		std::cerr << "Error: Cannot open file "
			<< amdc_data_path.lexically_normal() << "\n";
		return;
	}
	std::ofstream ion_fout(amdc_ion_data_path);
	if (!ion_fout.good()) {
		std::cerr << "Error: Could not create file "
			<< amdc_ion_data_path.lexically_normal() << "\n";
		return;
	}
	std::ofstream hardcode_fout(hardcode_data_path);
	if (!hardcode_fout.good()) {
		std::cerr << "Error: Could not create file "
			<< hardcode_data_path.lexically_normal() << "\n";
		return;
	}
	// skip first line
	std::string line;
	std::getline(fin, line);
	ion_fout << line << "\n";
	gong::NuclearData data;
	while (fin >> data) {
		double ion_mass = data.mass - data.z * electron_mass / atomic_mass_unit;
		hardcode_fout << "{{" << data.z << ", " << data.a
			<< "}, {" << data.z << ", " << data.a
			<< ", " << std::setprecision(17) << ion_mass << ", \"" << data.name
			<< "\"}},\n";
		ion_fout << data.z << " " << data.a << " " << data.name << " "
			<< std::setprecision(17)<< ion_mass << "\n";
	}
	// close file
	fin.close();
	ion_fout.close();
	hardcode_fout.close();
}

int main(int, char **argv) {
	gong::PathManager::Instance().Initialize(argv[0]);
	GenerateHardcodeData();
	return 0;
}