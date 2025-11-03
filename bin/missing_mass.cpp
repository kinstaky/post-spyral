#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <Math/Vector3D.h>

#include "include/particle.h"
#include "include/path_manager.h"
#include "include/nuclear_data.h"
#include "include/constants.h"
#include "include/event/solved_event.h"


using ROOT::Math::XYZVector;


int main(int, char **argv) {
	gong::PathManager::Instance().Initialize(argv[0]);

	// file name
	std::filesystem::path input_file_name =
		gong::PathManager::Instance().DataPath()
		/ "InterpSolver/solved_1H.root";
	// open file
	TFile ipf(input_file_name.lexically_normal().c_str(), "read");
	// tree
	TTree *ipt = (TTree *)ipf.Get("tree");
	if (ipt == nullptr) {
		std::cerr << "Error: Get tree from file failed: "
			<< input_file_name << std::endl;
		return -1;
	}
	// event
	SolvedEvent event;
	SetupInput(ipt, event);

	// output file name
	std::filesystem::path output_file_name =
		gong::PathManager::Instance().DataPath()
		/ "spectrum/missing_mass_1H.root";
	// output ROOT file
	TFile opf(output_file_name.lexically_normal().c_str(), "recreate");
	// output tree
	TTree opt("tree", "missing mass spectrum");
	// output data
	int valid;
	double p_polar, p_azimuthal, p_kinetic;
	double c_polar, c_azimuthal, c_kinetic;
	double c_excitation_energy;
	// setup output branches
	opt.Branch("valid", &valid, "valid/I");
	opt.Branch("p_polar", &p_polar, "ptheta/D");
	opt.Branch("p_azimuthal", &p_azimuthal, "pphi/D");
	opt.Branch("p_kinetic", &p_kinetic, "pk/D");
	opt.Branch("c_polar", &c_polar, "ctheta/D");
	opt.Branch("c_azimuthal", &c_azimuthal, "cphi/D");
	opt.Branch("c_kinetic", &c_kinetic, "ck/D");
	opt.Branch("c_excitation_energy", &c_excitation_energy, "cex/D");


	catima::Material target_gas = gong::GasMaterial(
		{{1, 1, 2}},
		300.0,
		1.0
	);

	const gong::Particle target(1, 1);

	// for showing process
	long long entries = ipt->GetEntries();
	long long progress = 0;
	long long last_percentage = 0;
	// showing process
	printf("Calculating spectrum with missing mass method   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		progress += 100;
		long long percentage = progress / entries;
		if (percentage != last_percentage) {
			printf(
				"\rCalculating spectrum with missing mass method %3lld%%",
				percentage
			);
			fflush(stdout);
			last_percentage = percentage;
		}

		ipt->GetEntry(entry);

		// init
		valid = 0;

		// check vertex point
		double r2 = event.vx * event.vx + event.vy * event.vy;
		double polar_in_degree = event.polar / gong::pi * 180.0;
		if (
			r2 > 0.00062 || fabs(event.vx) < 1e-5 || fabs(event.vy) < 1e-5
			|| event.vz < 0.005 || event.vz > 0.995 || polar_in_degree > 88.0
		) {
			valid = 1;
			opt.Fill();
			continue;
		}

		// calculate beam energy at vertex point
		gong::Particle beam(6, 16, 212.89);
		double distance_before_reaction =
			sqrt(event.vx*event.vx + event.vy*event.vy + event.vz*event.vz);
		target_gas.thickness_cm(distance_before_reaction * 100.0);
		beam.LostKineticEnergy(target_gas);

		gong::Particle scattered_proton(
			1, 1,
			event.kinetic,
			ROOT::Math::XYZVector(
				cos(event.azimuthal) * sin(event.polar),
				sin(event.azimuthal) * sin(event.polar),
				cos(event.polar)
			)
		);

		gong::Particle scattered_carbon = beam + target - scattered_proton;

		// // calculate scattering proton
		// p_kinetic = event.kinetic;
		// double p_energy = p_kinetic + mass_1h;
		// double p_momentum_value = gong::MomentumFromKinetic(mass_1h, p_kinetic);
		// XYZVector p_momentum(
		// 	cos(event.azimuthal) * sin(event.polar),
		// 	sin(event.azimuthal) * sin(event.polar),
		// 	cos(event.polar)
		// );
		// p_momentum *= p_momentum_value;

		// // calculate scattering carbon
		// double c_energy = beam.Energy() + target.Energy() - p_energy;
		// XYZVector c_momentum = beam.MomentumVector() - p_momentum;
		// double c_mass = sqrt(c_energy*c_energy - c_momentum.Mag2());


		p_polar = event.polar;
		p_azimuthal = event.azimuthal;
		c_polar = scattered_carbon.Polar();
		c_azimuthal = scattered_carbon.Azimuthal();
		c_kinetic = scattered_carbon.KineticEnergy();
		c_excitation_energy = scattered_carbon.ExcitationEnergy();
		// fill
		opt.Fill();
	}
	// finish
	printf("\rCalculating spectrum with missing mass method 100%%\n");


	// save
	opf.Write();
	// close files
	ipf.Close();
	opf.Close();

	return 0;
}