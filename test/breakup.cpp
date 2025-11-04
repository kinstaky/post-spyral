#include <fstream>
#include <iostream>
#include <string>

#ifdef DEBUG
#include <TGraph.h>
#include <TFile.h>
#include <TMath.h>
#endif

#include "include/constants.h"
#include "include/particle.h"
#include "include/path_manager.h"


using gong::PathManager;
using gong::Particle;
using gong::ParticleList;
using std::filesystem::path;

bool verbose = true;

struct InputInformation {
	Particle parent;
	Particle fragment0;
	Particle fragment1;
	path cangle_langle_file_name;
	path cangle_energy_file_name;
	double angle_eps;
	double energy_eps;
};


/// @brief compare value
/// @param[in] expected expected value
/// @param[in] actual actual value
/// @param[in] angle CM angle
/// @param[in] case_num test case number
/// @param[in] message message to show
/// @return true if pass, false if fail
bool Compare(
	double expected,
	double actual,
	double eps,
	double angle,
	int case_num,
	std::string message
) {
	if (std::fabs(expected - actual) > eps) {
		if (!verbose) return false;
		std::cerr << "Compare case " << case_num << " " << message
			<< " at CM angle " << angle << " failed:\n"
			<< "  " << expected << ", " << actual
			<< ", " << std::fabs(expected - actual) << std::endl;
		return false;
	}
	return true;
}


int main(int, char **argv) {
	// initialize
	PathManager::Instance().Initialize(argv[0]);
	path data_path = PathManager::Instance().ExecPath() / "../test/data/";


	// input information
	std::vector<InputInformation> input_informations = {
		{
			Particle(
				6, 14, 280.0, ROOT::Math::XYZVector(0.0, 0.0, 1.0), 0, 20.0
			),
			Particle(
				4, 10, 0.0, ROOT::Math::XYZVector(0.0, 0.0, 1.0), 0, 2.5
			),
			Particle(2, 4),
			data_path / "14C-10Be-4He-cangle-langle.txt",
			data_path / "14C-10Be-4He-cangle-energy.txt",
			1e-3,
			1e-2
		},
		{
			Particle(
				7, 12, 200.0, ROOT::Math::XYZVector(0.0, 0.0, 1.0), 0, 10.0
			),
			Particle(6, 11),
			Particle(1, 1),
			data_path / "12N-11C-1H-cangle-langle.txt",
			data_path / "12N-11C-1H-cangle-energy.txt",
			5e-3,
			5e-3
		},
		{
			Particle(
				6, 16, 350.0, ROOT::Math::XYZVector(0.0, 0.0, 1.0), 0, 25.0
			),
			Particle(
				4, 10, 0.0, ROOT::Math::XYZVector(0.0, 0.0, 1.0), 0, 2.5
			),
			Particle(
				2, 6, 0.0, ROOT::Math::XYZVector(0.0, 0.0, 1.0), 0, 1.0
			),
			data_path / "16C-10Be-6He-cangle-langle.txt",
			data_path / "16C-10Be-6He-cangle-energy.txt",
			5e-4,
			1e-2
		}
	};

#ifdef DEBUG
	path root_file_name = data_path / "breakup.root";
	TFile opf(root_file_name.lexically_normal().c_str(), "recreate");
	std::vector<TGraph> g_cangle_frag0_langle;
	std::vector<TGraph> g_cangle_frag0_energy;
	std::vector<TGraph> g_cangle_frag1_langle;
	std::vector<TGraph> g_cangle_frag1_energy;
	for (size_t index = 0; index < input_informations.size(); ++index) {
		g_cangle_frag0_langle.push_back(TGraph());
		g_cangle_frag0_energy.push_back(TGraph());
		g_cangle_frag1_langle.push_back(TGraph());
		g_cangle_frag1_energy.push_back(TGraph());
	}
#endif

	bool pass = true;
	// loop test cases
	for (size_t index = 0; index < input_informations.size(); ++index) {
		bool angle_pass = true;
		bool energy_pass = true;
		const InputInformation &info = input_informations[index];
		// read data
		// open file
		std::ifstream angle_fin(info.cangle_langle_file_name);
		std::ifstream energy_fin(info.cangle_energy_file_name);
		// read lines without formated data
		std::string line;
		for (int i = 0; i < 6; ++i) getline(angle_fin, line);
		for (int i = 0; i < 6; ++i) getline(energy_fin, line);

		double fragment0_mass_in_u =
			info.fragment0.Mass() / gong::atomic_mass_unit;
		double fragment1_mass_in_u =
			info.fragment1.Mass() / gong::atomic_mass_unit;

		while (angle_fin.good()) {
			double cm_angle, frag0_lab_angle, frag1_lab_angle;
			double frag0_energy, frag1_energy, tmp;
			angle_fin >> cm_angle >> frag0_lab_angle >> frag1_lab_angle;
			energy_fin >> tmp >> frag0_energy >> frag1_energy;

			gong::Particle fragment0(info.fragment0);
			gong::Particle fragment1(info.fragment1);
			gong::Breakup(
				info.parent,
				cm_angle / 180.0 * gong::pi,
				fragment0, fragment1
			);
			// check result
			angle_pass &= Compare(
				frag0_lab_angle,
				fragment0.Polar() / gong::pi * 180.0,
				info.angle_eps,
				cm_angle,
				index,
				"fragment0 lab angle"
			);
			angle_pass &= Compare(
				frag1_lab_angle,
				fragment1.Polar() / gong::pi * 180.0,
				info.angle_eps,
				cm_angle,
				index,
				"fragment1 lab angle"
			);
			energy_pass &= Compare(
				frag0_energy,
				fragment0.KineticEnergy() / fragment0_mass_in_u,
				info.energy_eps,
				cm_angle,
				index,
				"fragment0 energy"
			);
			energy_pass &= Compare(
				frag1_energy,
				fragment1.KineticEnergy() / fragment1_mass_in_u,
				info.energy_eps,
				cm_angle,
				index,
				"fragment1 energy"
			);

#ifdef DEBUG
			g_cangle_frag0_langle[index].AddPoint(
				cm_angle, fragment0.Polar() * TMath::RadToDeg()
			);
			g_cangle_frag0_energy[index].AddPoint(
				cm_angle, fragment0.KineticEnergy()
			);
			g_cangle_frag1_langle[index].AddPoint(
				cm_angle, fragment1.Polar() * TMath::RadToDeg()
			);
			g_cangle_frag1_energy[index].AddPoint(
				cm_angle, fragment1.KineticEnergy()
			);
#endif
		}

		if (!angle_pass) {
			std::cout << "Case " << index << " angle pass failed.\n";
		}
		if (!energy_pass) {
			std::cout << "Case " << index << " energy pass failed.\n";
		}

		pass &= angle_pass && energy_pass;
	}

#ifdef DEBUG
	for (size_t index = 0; index < input_informations.size(); ++index) {
		g_cangle_frag0_langle[index].Write(TString::Format("gf0a%ld", index));
		g_cangle_frag0_energy[index].Write(TString::Format("gf0e%ld", index));
		g_cangle_frag1_langle[index].Write(TString::Format("gf1a%ld", index));
		g_cangle_frag1_energy[index].Write(TString::Format("gf1e%ld", index));
	}
	opf.Close();
#endif



	return pass ? 0 : 1;
}