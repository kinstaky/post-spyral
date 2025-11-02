#include <catima/catima.h>
#include "include/nuclear_data.h"
#include "include/path_manager.h"
#include "include/particle.h"

constexpr double energy_output[5] = {
	194.11610256329004,
	190.26862653596302,
	188.33812261604783,
	186.39136051917055,
	175.44198856910415
};
constexpr double eps = 1e-3;

bool verbose = false;

int main(int, char **argv) {
	// initialize path manager
	gong::PathManager::Instance().Initialize(argv[0]);

	catima::Material window_material = gong::SolidMaterial(
		{{1, 1, 14}, {12, 6, 14}, {14, 7, 4}, {16, 8, 4}},
		1422.312 * 1e-3
	);
	catima::Material ic_gas_material = gong::GasMaterial(
		{{12, 6, 1}, {18, 9, 4}},
		200.0,
		0.035
	);
	catima::Material target_material = gong::GasMaterial(
		{{1, 1, 2}},
		300.0,
		1.0
	);


	gong::Particle carbon16(6, 16, 196.0);
	carbon16.LostKineticEnergy(window_material);
	double energy_ic = carbon16.KineticEnergy();
	if (verbose) {
		std::cout << "Energy ic difference = "
			<< fabs(energy_ic - energy_output[0]) << std::endl;
	}
	if (fabs(energy_ic - energy_output[0]) > eps) {
		std::cerr << "Energy ic difference = "
			<< fabs(energy_ic - energy_output[0]) << std::endl;
		return 1;
	}
	// kinetic energy after the IC gas
	carbon16.LostKineticEnergy(ic_gas_material);
	double energy_ic_exit = carbon16.KineticEnergy();
	if (verbose) {
		std::cout << "Energy ic exit difference = "
			<< fabs(energy_ic_exit - energy_output[1]) << std::endl;
	}
	if (fabs(energy_ic_exit - energy_output[1]) > eps) {
		std::cerr << "Energy ic exit difference = "
			<< fabs(energy_ic_exit - energy_output[1]) << std::endl;
		return 1;
	}
	// kinetic energy after the second IC window
	carbon16.LostKineticEnergy(window_material);
	double energy_post_ic = carbon16.KineticEnergy();
	if (verbose) {
		std::cout << "Energy post ic difference = "
			<< fabs(energy_post_ic - energy_output[2]) << std::endl;
	}
	if (fabs(energy_post_ic - energy_output[2]) > eps) {
		std::cerr << "Energy post ic difference = "
			<< fabs(energy_post_ic - energy_output[2]) << std::endl;
		return 1;
	}
	// kinetic energy after the AT-TPC window
	carbon16.LostKineticEnergy(window_material);
	double energy_start = carbon16.KineticEnergy();
	if (verbose) {
		std::cout << "Energy attpc difference = "
		<< fabs(energy_start - energy_output[3]) << std::endl;
	}
	if (fabs(energy_start - energy_output[3]) > eps) {
		std::cerr << "Energy attpc difference = "
			<< fabs(energy_start - energy_output[3]) << std::endl;
		return 1;
	}
	// kinetic energy after the AT-TPC gas
	carbon16.LostKineticEnergy(target_material);
	double energy_exit = carbon16.KineticEnergy();
	if (verbose) {
		std::cout << "Energy attpc exit difference = "
		<< fabs(energy_exit - energy_output[4]) << std::endl;
	}
	if (fabs(energy_exit - energy_output[4]) > eps) {
		std::cerr << "Energy attpc exit difference = "
			<< fabs(energy_exit - energy_output[4]) << std::endl;
		return 1;
	}

	return 0;
}