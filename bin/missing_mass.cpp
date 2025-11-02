#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <Math/Vector3D.h>

#include "include/particle.h"
#include "include/path_manager.h"
#include "include/nuclear_data.h"
#include "include/constants.h"


using ROOT::Math::XYZVector;

struct SolvedEvent {
	int run, event;
	double vx, vy, vz;
	double brho, kinetic;
	double polar, azimuthal;
};

int main(int, char **argv) {
	gong::PathManager::Instance().Initialize(argv[0]);

	// file name
	TString input_file_name = TString::Format(
		"/data/a1975/InterpSolver/solved_1H.root"
	);
	// open file
	TFile ipf(input_file_name, "read");
	// tree
	TTree *ipt = (TTree *)ipf.Get("tree");
	if (ipt == nullptr) {
		std::cerr << "Error: cannot get tree " << "tree" << std::endl;
		return -1;
	}
	// event
	SolvedEvent event;
	ipt->SetBranchAddress("origin_run", &event.run);
	ipt->SetBranchAddress("origin_event", &event.event);
	ipt->SetBranchAddress("vertex_x", &event.vx);
	ipt->SetBranchAddress("vertex_y", &event.vy);
	ipt->SetBranchAddress("vertex_z", &event.vz);
	ipt->SetBranchAddress("brho", &event.brho);
	ipt->SetBranchAddress("kinetic_energy", &event.kinetic);
	ipt->SetBranchAddress("polar", &event.polar);
	ipt->SetBranchAddress("azimuthal", &event.azimuthal);


	// output file name
	TString output_file_name = TString::Format(
		"/data/a1975/spectrum/missing_mass_1H.root"
	);
	// output ROOT file
	TFile opf(output_file_name, "recreate");
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

	const catima::Material window = gong::SolidMaterial(
		{{1, 1, 14}, {12, 6, 14}, {14, 7, 4}, {16, 8, 4}},
		1422.312 * 1e-3
	);
	const catima::Material ic_gas = gong::GasMaterial(
		{{12, 6, 1}, {18, 9, 4}},
		200.0,
		0.035
	);
	catima::Material target_gas = gong::GasMaterial(
		{{1, 1, 2}},
		300.0,
		1.0
	);

	const gong::Particle target(1, 1);
	const double mass_1h = gong::GetMass(1, 1);
	const double mass_16c = gong::GetMass(6, 16);

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
		if (
			r2 > 0.00062 || fabs(event.vx) < 1e-5 || fabs(event.vy) < 1e-5
			|| event.vz < 0.005 || event.vz > 0.995 || event.polar/gong::pi*180.0 > 88.0
		) {
			valid = 1;
			opt.Fill();
			continue;
		}

		// calculate beam energy at vertex point
		gong::Particle beam(6, 16, 184.0);
		double distance_before_reaction =
			sqrt(event.vx*event.vx + event.vy*event.vy + event.vz*event.vz);
		target_gas.thickness_cm(distance_before_reaction * 100.0);
		beam.LostKineticEnergy(window)
			.LostKineticEnergy(ic_gas)
			.LostKineticEnergy(window)
			.LostKineticEnergy(window)
			.LostKineticEnergy(target_gas);

		// calculate scattering proton
		p_kinetic = event.kinetic;
		double p_energy = p_kinetic + mass_1h;
		double p_momentum_value = gong::MomentumFromKinetic(mass_1h, p_kinetic);
		XYZVector p_momentum(
			cos(event.azimuthal) * sin(event.polar),
			sin(event.azimuthal) * sin(event.polar),
			cos(event.polar)
		);
		p_momentum *= p_momentum_value;

		// calculate scattering carbon
		double c_energy = beam.Energy() + target.Energy() - p_energy;
		XYZVector c_momentum = beam.MomentumVector() - p_momentum;
		double c_mass = sqrt(c_energy*c_energy - c_momentum.Mag2());


		p_polar = event.polar;
		p_azimuthal = event.azimuthal;
		c_polar = c_momentum.Theta();
		c_azimuthal = c_momentum.Phi();
		c_kinetic = c_energy - c_mass;
		c_excitation_energy = c_mass - mass_16c;
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