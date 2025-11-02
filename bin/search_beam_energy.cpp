#include <iostream>
#include <fstream>
#include <algorithm>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TTree.h>
#include <TString.h>
#include <TCutG.h>
#include <TMultiGraph.h>
#include <TMath.h>

#include "include/constants.h"
#include "include/particle.h"
#include "include/path_manager.h"
#include "include/event/solved_event.h"

int main(int, char **argv) {
	gong::PathManager::Instance().Initialize(argv[0]);

	// input file name
	std::filesystem::path input_file_name =
		gong::PathManager::Instance().DataPath()
		/ "InterpSolver/solved_1H.root";
	// open input file
	TFile ipf(input_file_name.lexically_normal().c_str(), "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from file failed: "
			<< input_file_name << std::endl;
		return -1;
	}
	// solved event
	SolvedEvent event;
	// setup input branches
	SetupInput(ipt, event);


	// read cut
	std::filesystem::path cut_file_name =
		gong::PathManager::Instance().DataPath()
		/ "cuts/C16-p-elastic-Etheta.txt";
	// open file
	std::ifstream cut_fin(cut_file_name);
	if (!cut_fin) {
		std::cerr << "Error: Open cut file failed: "
			<< cut_file_name << std::endl;
		return -1;
	}
	int cut_points = 0;
	std::vector<double> cut_x, cut_y;
	cut_fin >> cut_points;
	for (int i = 0; i < cut_points; ++i) {
		double x, y;
		cut_fin >> x >> y;
		cut_x.push_back(x);
		cut_y.push_back(y);
	}
	TCutG elastic_cut("elastic_cut", cut_points, cut_x.data(), cut_y.data());


	// output file name
	std::filesystem::path output_file_name =
		gong::PathManager::Instance().DataPath()
		/ "spectrum/search_beam_energy.root";
	// create output file
	TFile opf(output_file_name.lexically_normal().c_str(), "recreate");
	// E-theta for checking
	TH2F e_theta("etheta", "E-theta for proton", 1000, 0, 90, 1000, 0, 40);
	// theta distribution
	TH1F theta_distribution("theta", "theta distribution", 90, 0, 90);
	TGraph g_e_theta;
	TGraph g_difference;

	// search start
	constexpr double search_start = 170.0;
	// search stop
	constexpr double search_stop = 230.0;
	// search step
	constexpr double search_step = 1.0;
	// search number
	constexpr int search_num = (search_stop - search_start) / search_step + 1;
	double difference[search_num];
	for (int i = 0; i < search_num; ++i) difference[i] = 0.0;

	// for diluting
	int dilute_total[15] = {
		2, 5, 4, 6, 8, 10, 12, 14, 16, 18, 18, 16, 12, 7, 5
	};
	int dilute_threshold[15] = {
		1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2
	};
	int dilute_progress[15];
	for (int i = 0; i < 15; ++i) dilute_progress[i] = 0;

	// H2 gas target
	catima::Material target_gas = gong::GasMaterial(
		{{1, 1, 2}},
		300.0,
		1.0
	);
	gong::Particle target(1, 1);


	constexpr double generate_start = 150.0;
	constexpr double generate_stop = 240.0;
	constexpr double generate_step = 1.0;
	constexpr int generate_num = (generate_stop - generate_start) / generate_step + 1;
	TGraph g_e_theta_sim[generate_num];
	// generate kinematic data for different beam energy
	printf("Generating kinematic data  0/%d", generate_num);
	fflush(stdout);
	double kinetic_energy_data[generate_num][1800];
	double polar_data[generate_num][1800];
	for (int i = 0; i < generate_num; ++i) {
		double beam_energy = generate_start + i * generate_step;
		for (int j = 0; j < 1800; ++j) {
			double theta = 0.1 + j * 0.1;
			gong::Particle beam(6, 16, beam_energy);
			gong::ParticleList scattered = gong::Scatter(beam, target, theta);
			kinetic_energy_data[i][j] = scattered[1].KineticEnergy();
			polar_data[i][j] = scattered[1].Polar() * TMath::RadToDeg();
		}
		g_e_theta_sim[i] = TGraph(1800, polar_data[i], kinetic_energy_data[i]);
		printf("\rGenerating kinematic data %2d/%d", i+1, generate_num);
		fflush(stdout);
	}
	printf("\n");


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

		// check vertex point
		double r2 = event.vx * event.vx + event.vy * event.vy;
		if (
			r2 > 0.00062 || fabs(event.vx) < 1e-5 || fabs(event.vy) < 1e-5
			|| event.vz < 0.005 || event.vz > 0.995
		) {
			continue;
		}
		double polar_in_degree = event.polar / gong::pi * 180.0;
		if (polar_in_degree > 88.0) continue;
		if (!elastic_cut.IsInside(polar_in_degree, event.kinetic)) continue;

		// fill full E-theta
		e_theta.Fill(polar_in_degree, event.kinetic);

		if (int(polar_in_degree) >= 70 && int(polar_in_degree) <= 84) {
			int index = int(polar_in_degree) - 70;
			dilute_progress[index] = (dilute_progress[index] + 1) % dilute_total[index];
			if (dilute_progress[index] >= dilute_threshold[index]) continue;
		}

		// fill diluted theta
		theta_distribution.Fill(polar_in_degree);
		g_e_theta.AddPoint(polar_in_degree, event.kinetic);

		// prepare for scattering simulation
		double distance_before_reaction =
			sqrt(event.vx*event.vx + event.vy*event.vy + event.vz*event.vz);
		target_gas.thickness_cm(distance_before_reaction * 100.0);
		int search_counts = 0;
		for (
			double energy = search_start;
			energy <= search_stop;
			energy += search_step
		) {
			gong::Particle beam(6, 16, energy);
			beam.LostKineticEnergy(target_gas);
			// energy index
			int ei = int(
				(beam.KineticEnergy() - generate_start) / generate_step
			);
			if (ei < 0 || ei >= generate_num-1) {
				std::cerr << "Error: Energy index out of range: "
					<< ei << std::endl;
				return -1;
			}
			// polar index of lower energy data
			int pi_le = std::lower_bound(
				polar_data[ei],
				polar_data[ei]+1800,
				polar_in_degree,
				std::greater<double>()
			) - polar_data[ei];
			// lower bound of kinetic energy
			double kinetic_le = 0.0;
			if (pi_le == 0) {
				kinetic_le = kinetic_energy_data[ei][0];
			} else if (pi_le == 1800) {
				kinetic_le = kinetic_energy_data[ei][1799];
			} else {
				// theta in decrease order
				double theta_upper = polar_data[ei][pi_le-1];
				double theta_lower = polar_data[ei][pi_le];
				double ratio =
					(polar_in_degree - theta_lower)
					/ (theta_upper - theta_lower);
				double kinetic_lower = kinetic_energy_data[ei][pi_le-1];
				double kinetic_upper = kinetic_energy_data[ei][pi_le];
				kinetic_le = kinetic_lower * ratio + kinetic_upper * (1.0 - ratio);
			}
			// polar index for greater energy data
			int pi_ge = std::lower_bound(
				polar_data[ei+1],
				polar_data[ei+1]+1800,
				polar_in_degree,
				std::greater<double>()
			) - polar_data[ei+1];
			// upper bound of kinetic energy
			double kinetic_ge = 0.0;
			if (pi_ge == 0) {
				kinetic_ge = kinetic_energy_data[ei+1][0];
			} else if (pi_ge == 1800) {
				kinetic_ge = kinetic_energy_data[ei+1][1799];
			} else {
				// theta in decrease order
				double theta_upper = polar_data[ei+1][pi_ge-1];
				double theta_lower = polar_data[ei+1][pi_ge];
				double ratio =
					(polar_in_degree - theta_lower)
					/ (theta_upper - theta_lower);
				double kinetic_lower = kinetic_energy_data[ei+1][pi_ge-1];
				double kinetic_upper = kinetic_energy_data[ei+1][pi_ge];
				kinetic_ge = kinetic_lower * ratio + kinetic_upper * (1.0 - ratio);
			}

			double beam_energy_upper = search_start + (ei+1) * search_step;
			double e_ratio =
				(beam.KineticEnergy() - beam_energy_upper)
				/ search_step;

			double kinetic_sim = kinetic_le * e_ratio + kinetic_ge * (1.0 - e_ratio);
			difference[search_counts] += pow(kinetic_sim - event.kinetic, 2.0);

// std::cout << "Entry " << entry
// 	<< ", search energy " << energy
// 	<< ", kinetic energy " << beam.KineticEnergy()
// 	<< " , polar " << polar_in_degree
// 	<< "\n"
// 	<< "  Found beam energy index " << ei
// 	<< "\n"
// 	<< "  Polar index le " << pi_le
// 	<< ", between " << polar_data[ei][pi_le-1]
// 	<< ", " << polar_data[ei][pi_le]
// 	<< ", lower kinetic energy " << kinetic_le
// 	<< "\n"
// 	<< "  Polar index ge " << pi_ge
// 	<< ", between " << polar_data[ei+1][pi_ge-1]
// 	<< ", " << polar_data[ei+1][pi_ge]
// 	<< ", upper kinetic energy " << kinetic_ge
// 	<< "\n"
// 	<< "  Get kinetic simulated " << kinetic_sim
// 	<< ", expected " << event.kinetic
// 	<< ", difference " << difference[search_counts]
// 	<< std::endl;
			++search_counts;
		}
	}
	// finish
	printf("\rCalculating spectrum with missing mass method 100%%\n");

	// fill energy difference graph
	int dindex = 0;
	for (
		double energy = search_start;
		energy <= search_stop;
		energy += search_step
	) {
		g_difference.AddPoint(energy, difference[dindex]);
		++dindex;
	}

	// save
	// E-theta
	e_theta.Write();
	// theta distribution
	theta_distribution.Write();
	// simulated E-theta with different energy
	for (int i = 0; i < search_num; ++i) {
		g_e_theta_sim[i].Write(TString::Format("gs%d", int(i+search_start)));
	}
	// energy difference in different beam energy hypothesis
	g_difference.Write("gd");
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}