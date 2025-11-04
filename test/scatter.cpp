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

bool verbose = false;

struct InputInformation {
	Particle beam;
	Particle target;
	double beam_excitation;
	double target_excitation;
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
			Particle(6, 16, 192.137),	// 16C in 12MeV/u
			Particle(1, 1),	// proton
			0.0,
			0.0,
			data_path / "16C-1H-cangle-langle.txt",
			data_path / "16C-1H-cangle-energy.txt",
			1e-4,
			5e-4
		},
		{
			Particle(6, 14, 350.0),	// 14C in 350MeV
			Particle(1, 2),	// deuteron
			18.0,
			0.0,
			data_path / "14C-2H-cangle-langle.txt",
			data_path / "14C-2H-cangle-energy.txt",
			1e-3,
			5e-4
		},
		{
			Particle(10, 20, 600.0),	// 20Ne in 600MeV
			Particle(4, 9),	// 9Be
			6.726,
			1.684,
			data_path / "20Ne-9Be-cangle-langle.txt",
			data_path / "20Ne-9Be-cangle-energy.txt",
			1e-1,
			1e-4
		}
	};

#ifdef DEBUG
	path root_file_name = data_path / "scatter.root";
	TFile opf(root_file_name.lexically_normal().c_str(), "recreate");
	std::vector<TGraph> g_cangle_beam_langle;
	std::vector<TGraph> g_cangle_beam_energy;
	std::vector<TGraph> g_cangle_target_langle;
	std::vector<TGraph> g_cangle_target_energy;
	for (size_t index = 0; index < input_informations.size(); ++index) {
		g_cangle_beam_langle.push_back(TGraph());
		g_cangle_beam_energy.push_back(TGraph());
		g_cangle_target_langle.push_back(TGraph());
		g_cangle_target_energy.push_back(TGraph());
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

		double scattered_beam_mass_in_u =
			info.beam.Mass() / gong::atomic_mass_unit;
		double scattered_target_mass_in_u =
			info.target.Mass() / gong::atomic_mass_unit;

		while (angle_fin.good()) {
			double cm_angle, beam_lab_angle, target_lab_angle;
			double beam_energy, target_energy, tmp;
			angle_fin >> cm_angle >> beam_lab_angle >> target_lab_angle;
			energy_fin >> tmp >> beam_energy >> target_energy;

			gong::Particle fragment0(info.beam);
			fragment0.SetExcitationEnergy(info.beam_excitation);
			gong::Particle fragment1(info.target);
			fragment1.SetExcitationEnergy(info.target_excitation);
			gong::Scatter(
				info.beam, info.target,
				cm_angle / 180.0 * gong::pi,
				fragment0, fragment1
			);
			// check result
			angle_pass &= Compare(
				beam_lab_angle,
				fragment0.Polar() / gong::pi * 180.0,
				info.angle_eps,
				cm_angle,
				index,
				"beam lab angle"
			);
			angle_pass &= Compare(
				target_lab_angle,
				fragment1.Polar() / gong::pi * 180.0,
				info.angle_eps,
				cm_angle,
				index,
				"target lab angle"
			);
			energy_pass &= Compare(
				beam_energy,
				fragment0.KineticEnergy() / scattered_beam_mass_in_u,
				info.energy_eps,
				cm_angle,
				index,
				"scattered beam energy"
			);
			energy_pass &= Compare(
				target_energy,
				fragment1.KineticEnergy() / scattered_target_mass_in_u,
				info.energy_eps,
				cm_angle,
				index,
				"scattered target energy"
			);

#ifdef DEBUG
			g_cangle_beam_langle[index].AddPoint(
				cm_angle, fragment0.Polar() * TMath::RadToDeg()
			);
			g_cangle_beam_energy[index].AddPoint(
				cm_angle, fragment0.KineticEnergy()
			);
			g_cangle_target_langle[index].AddPoint(
				cm_angle, fragment1.Polar() * TMath::RadToDeg()
			);
			g_cangle_target_energy[index].AddPoint(
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
		g_cangle_beam_langle[index].Write(TString::Format("gba%ld", index));
		g_cangle_beam_energy[index].Write(TString::Format("gbe%ld", index));
		g_cangle_target_langle[index].Write(TString::Format("gta%ld", index));
		g_cangle_target_energy[index].Write(TString::Format("gte%ld", index));
	}
	opf.Close();
#endif



	return pass ? 0 : 1;
}