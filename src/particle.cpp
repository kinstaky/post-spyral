#include "include/particle.h"
#include "include/nuclear_data.h"
#include "include/constants.h"

namespace gong {

Particle::Particle(
	int z,
	int a,
	double kinetic,
	const ROOT::Math::XYZVector &direction,
	int q,
	double excitation
) :
	z_(z),
	a_(a),
	q_(q),
	excitation_(excitation),
	kinetic_(kinetic),
	direction_(direction.Unit())
{
	// get mass
	mass_ = GetMass(z_, a_);
	// calculate energy
	energy_ = mass_ + excitation_ + kinetic_;
	// calculate momentum
	momentum_ = MomentumFromKinetic(mass_ + excitation_, kinetic_);
}


Particle& Particle::SetMass0(double mass) {
	mass_ = mass;
	// update energy
	energy_ = mass_ + excitation_ + kinetic_;
	// update momentum
	momentum_ = MomentumFromKinetic(mass_ + excitation_, kinetic_);
	return *this;
}


Particle& Particle::LostKineticEnergy(catima::Material material) {
	material.thickness(material.thickness() / cos(direction_.Theta()));
	double mass_in_u = (mass_ + excitation_) / atomic_mass_unit;
	catima::Projectile projectile(mass_in_u, z_, 0, kinetic_ / mass_in_u);
	// calculate kinetic energy lost in the material
	kinetic_ = catima::energy_out(projectile, material) * mass_in_u;
	// update energy
	energy_ = mass_ + excitation_ + kinetic_;
	// update momentum
	momentum_ = MomentumFromKinetic(mass_ + excitation_, kinetic_);
	return *this;
}



Particle& Particle::AddKineticEnergy(double energy) {
	kinetic_ += energy;
	// update energy
	energy_ = mass_ + excitation_ + kinetic_;
	// update momentum
	momentum_ = MomentumFromKinetic(mass_ + excitation_, kinetic_);
	return *this;
}


Particle& Particle::SetKineticEnergy(double energy) {
	kinetic_ = energy;
	// update energy
	energy_ = mass_ + excitation_ + kinetic_;
	// update momentum
	momentum_ = MomentumFromKinetic(mass_ + excitation_, kinetic_);
	return *this;
}


Particle& Particle::SetExcitationEnergy(double excitation) {
	excitation_ = excitation;
	// update energy
	energy_ = mass_ + excitation_ + kinetic_;
	// update momentum
	momentum_ = MomentumFromKinetic(mass_ + excitation_, kinetic_);
	return *this;
}


Particle& Particle::SetDirection(const ROOT::Math::XYZVector &direction) {
	direction_ = direction.Unit();
	return *this;
}


Particle& Particle::SetMomentum(const ROOT::Math::XYZVector &momentum) {
	direction_ = momentum.Unit();
	momentum_ = momentum.R();
	// update kinetic energy
	kinetic_ = KineticFromMomentum(mass_ + excitation_, momentum_);
	// update energy
	energy_ = mass_ + excitation_ + kinetic_;
	return *this;
}


Particle& Particle::SetMomentum(double momentum) {
	momentum_ = momentum;
	// update kinetic energy
	kinetic_ = KineticFromMomentum(mass_ + excitation_, momentum_);
	// update energy
	energy_ = mass_ + excitation_ + kinetic_;
	return *this;
}


catima::Material SolidMaterial(
	std::initializer_list<std::array<double,3>> compound,
	double thickness
) {
	catima::Material material;
	for (auto &element : compound) {
		material.add_element(
			gong::GetMassInUnit(element[1], element[0]),
			element[1],
			element[2]
		);
	}
	material.thickness(thickness * 1e-3);
	return material;
}


catima::Material GasMaterial(
	std::initializer_list<std::array<double,3>> compound,
	double pressure,
	double length
) {
	catima::Material material;
	int molar_mass = 0;	// g/mol
	for (auto &element : compound) {
		material.add_element(
			gong::GetMassInUnit(element[1], element[0]),
			element[1],
			element[2]
		);
		molar_mass += element[0] * element[2];
	}
	// convert to Pa(kg*m^-1*s^-2), 0.0075006 = 760 / 1.01325e5
	pressure /= 0.0075006;
	// rho = g/mol * Pa / (J*K-1*mol-1*K) * (cm^3/m^3)
	//     = g/mol * kg*m^-1*s^-2 / (kg*m^2*s^-2*mol^-1*K^-1 * K) * (cm^3/m^3)
	material.density(molar_mass * pressure / (gas_constant * 293.15) * 1e-6);
	material.thickness_cm(length * 100.0);
	return material;
}


ParticleList Scatter(
	const Particle &beam,
	const Particle &target,
	double angle,
	double beam_excitation,
	double target_excitation
) {
	ParticleList result;
	result.push_back(beam);
	result[0].SetExcitationEnergy(beam_excitation);
	result.push_back(target);
	result[1].SetExcitationEnergy(target_excitation);
	double beta_mass_center = beam.Momentum() / (beam.Energy() + target.Mass());
	double gamma_mass_center =
		1.0 / sqrt(1.0 - beta_mass_center * beta_mass_center);
	double reaction_energy = sqrt(
		(beam.Energy() + target.Mass() + beam.Momentum())
		* (beam.Energy() + target.Mass() - beam.Momentum())
	);
	// momentum of exit particle or recoil particle in center of mass frame
	double exit_momentum_center =
		sqrt(
			(reaction_energy - result[0].Mass() - result[1].Mass())
			* (reaction_energy - result[0].Mass() + result[1].Mass())
			* (reaction_energy + result[0].Mass() - result[1].Mass())
			* (reaction_energy + result[0].Mass() + result[1].Mass())
		) / (2.0 * reaction_energy);
	// exit momentum parallel and vertical part
	double exit_momentum_center_parallel = exit_momentum_center * cos(angle);
	double exit_momentum_center_vertical = exit_momentum_center * sin(angle);
	// exit energy in c.m.
	double exit_energy_center = sqrt(
		exit_momentum_center * exit_momentum_center
		+ result[0].Mass() * result[0].Mass()
	);
	// exit energy in lab frame
	double exit_energy = gamma_mass_center * exit_energy_center
		+ gamma_mass_center * beta_mass_center * exit_momentum_center_parallel;
	// exit momentum parallel part in lab frame
	double exit_momentum_parallel = gamma_mass_center * exit_momentum_center_parallel
		+ gamma_mass_center * beta_mass_center * exit_energy_center;
	// exit momentum vertical part in lab frame
	double exit_momentum_vertical = exit_momentum_center_vertical;
	// calculate exit angle in lab frame
	double exit_angle = fabs(atan(exit_momentum_vertical / exit_momentum_parallel));
	exit_angle = exit_momentum_parallel > 0 ?
		exit_angle : pi - exit_angle;
	result[0].SetKineticEnergy(exit_energy - result[0].Mass());
	result[0].SetDirection(ROOT::Math::XYZVector(
		sin(exit_angle), 0.0, cos(exit_angle)
	));

	// p of recoil particle in center of mass frame
	double recoil_momentum_center = exit_momentum_center;
	// recoil momentum parallel and vertical part
	double recoil_momentum_center_parallel = -recoil_momentum_center * cos(angle);
	double recoil_momentum_center_vertical = -recoil_momentum_center * sin(angle);
	// recoil energy in c.m.
	double recoil_energy_center = sqrt(
		recoil_momentum_center * recoil_momentum_center
		+ result[1].Mass() * result[1].Mass()
	);
	// recoil energy in lab frame
	double recoil_energy = gamma_mass_center * recoil_energy_center
		+ gamma_mass_center * beta_mass_center * recoil_momentum_center_parallel;
	// recoil momentum parallel part in lab frame
	double recoil_momentum_parallel =
		gamma_mass_center * recoil_momentum_center_parallel
		+ gamma_mass_center * beta_mass_center * recoil_energy_center;
	// recoil momentum vertical part in lab frame
	double recoil_momentum_vertical = recoil_momentum_center_vertical;
	// recoil angle in lab frame
	double recoil_angle =
		fabs(atan(recoil_momentum_vertical / recoil_momentum_parallel));
	recoil_angle = recoil_momentum_parallel > 0 ?
		recoil_angle : pi - recoil_angle;
	result[1].SetKineticEnergy(recoil_energy - result[1].Mass());
	result[1].SetDirection(ROOT::Math::XYZVector(
		sin(recoil_angle), 0.0, cos(recoil_angle)
	));

	return result;
}


};	//	namespace gong