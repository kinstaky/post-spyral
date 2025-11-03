#include <iostream>

#include <cxxopts.hpp>

#include "include/path_manager.h"
#include "include/particle.h"


double BinarySearch(
	double low,
	double high,
	double target,
	const catima::Material &window,
	const catima::Material &ic_gas,
	double eps,
	int iteration
);

int main(int argc, char **argv) {
	cxxopts::Options options("backtrack_beam_energy", "Backtrack beam energy");
	options.add_options()
		("h,help", "Print help")
		("energy", "Beam energy in ATTPC upstream", cxxopts::value<double>());
	options.parse_positional({"energy"});
	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		std::cout << options.help() << std::endl;
		return 0;
	}
	if (!result.count("energy")) {
		std::cout << "Error: energy not specified\n\n"
			<< options.help() << std::endl;
		return -1;
	}

	gong::PathManager::Instance().Initialize(argv[0]);

	double beam_energy = result["energy"].as<double>();
	const catima::Material window = gong::SolidMaterial(
		{{1, 1, 14}, {12, 6, 14}, {14, 7, 4}, {16, 8, 4}},
		1422.312 * 1e-3
	);
	const catima::Material ic_gas = gong::GasMaterial(
		{{12, 6, 1}, {18, 9, 4}},
		200.0,
		0.035
	);

	double backward_energy = BinarySearch(
		beam_energy-0.5,
		beam_energy+30.0,
		beam_energy,
		window,
		ic_gas,
		1e-3,
		0
	);

	std::cout << backward_energy << std::endl;

	gong::Particle beam(6, 16, beam_energy);
	beam.LostKineticEnergy(window)
		.LostKineticEnergy(ic_gas)
		.LostKineticEnergy(window)
		.LostKineticEnergy(window);

	return 0;
}


double BinarySearch(
	double low,
	double high,
	double target,
	const catima::Material &window,
	const catima::Material &ic_gas,
	double eps,
	int iteration
) {

	double mid = (low + high) / 2.0;
	if (iteration == 1000) return mid;
	if (high - low < eps) return mid;
	// low
	double low_value = gong::Particle(6, 16, low)
		.LostKineticEnergy(window)
		.LostKineticEnergy(ic_gas)
		.LostKineticEnergy(window)
		.LostKineticEnergy(window)
		.KineticEnergy();
	double mid_value = gong::Particle(6, 16, mid)
		.LostKineticEnergy(window)
		.LostKineticEnergy(ic_gas)
		.LostKineticEnergy(window)
		.LostKineticEnergy(window)
		.KineticEnergy();
	double high_value = gong::Particle(6, 16, high)
		.LostKineticEnergy(window)
		.LostKineticEnergy(ic_gas)
		.LostKineticEnergy(window)
		.LostKineticEnergy(window)
		.KineticEnergy();

	if (fabs(low_value - target) < eps) return low;
	if (fabs(mid_value - target) < eps) return mid;
	if (fabs(high_value - target) < eps) return high;
	if (mid_value < target) {
		return BinarySearch(mid, high, target, window, ic_gas, eps, iteration+1);
	} else {
		return BinarySearch(low, mid, target, window, ic_gas, eps, iteration+1);
	}
}
