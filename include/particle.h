#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include <string>
#include <catima/catima.h>
#include <Math/Vector3D.h>

namespace gong {

class Particle {
public:
	Particle(
		int z,
		int a,
		double kinetic = 0.0,
		const ROOT::Math::XYZVector &direction = ROOT::Math::XYZVector(0, 0, 1),
		int q = 0,
		double excitation = 0.0
	);


	/// @brief default copy constructor
	/// @param other other particle
	Particle(const Particle &other) = default;


	/// @brief default destructor
	~Particle() = default;


	/// @brief set mass at ground state in MeV/c^2
	/// @param mass mass at ground state in MeV/c^2
	/// @returns reference to self
	Particle& SetMass0(double mass);


	/// @brief get mass at ground state in MeV/c^2
	/// @returns mass at ground state in MeV/c^2
	inline double Mass0() const { return mass_; }


	/// @brief get mass (including excitation energy) in MeV/c^2
	/// @returns mass (including excitation energy) in MeV/c^2
	inline double Mass() const { return mass_ + excitation_; }


	/// @brief calculate kinetic energy lost in a material
	/// @param material material in catima
	/// @returns reference to self
	Particle& LostKineticEnergy(catima::Material material);


	/// @brief add kinetic energy in MeV
	/// @param[in] energy energy in MeV
	/// @returns reference to self
	Particle& AddKineticEnergy(double energy);


	/// @brief set kinetic energy in MeV
	/// @param[in] energy energy in MeV
	/// @returns reference to self
	Particle& SetKineticEnergy(double energy);


	/// @brief get kinetic energy in MeV
	/// @returns kinetic energy in MeV
	inline double KineticEnergy() const { return kinetic_; }


	/// @brief get kinetic energy in MeV
	/// @returns kinetic energy in MeV
	inline double Energy() const { return energy_; }


	/// @brief set excitation energy in MeV
	/// @param[in] energy excitation energy in MeV
	/// @returns reference to self
	Particle& SetExcitationEnergy(double energy);


	/// @brief get excitation energy in MeV
	/// @returns excitation energy in MeV
	inline double ExcitationEnergy() const { return excitation_; }


	/// @brief set direction
	/// @param[in] direction direction vector
	/// @returns reference to self
	Particle& SetDirection(const ROOT::Math::XYZVector &direction);


	/// @brief get direction vector
	/// @returns direction vector
	inline ROOT::Math::XYZVector Direction() const { return direction_; }

	/// @brief get direction angle in radian
	/// @returns direction angle in radian
	inline double Polar() const { return direction_.Theta(); }


	/// @brief get direction angle in radian
	/// @returns direction angle in radian
	inline double Azimuthal() const { return direction_.Phi(); }


	/// @brief set momentum in MeV/c
	/// @param[in] momentum momentum vector in MeV/c
	/// @returns reference to self
	Particle& SetMomentum(const ROOT::Math::XYZVector &momentum);


	/// @brief set momentum in MeV/c
	/// @param[in] momentum momentum in MeV/c
	/// @returns reference to self
	Particle& SetMomentum(double momentum);


	/// @brief get momentum value in MeV/c
	/// @return momentum value in MeV/c
	inline double Momentum() const { return momentum_; }


	/// @brief get momentum in MeV/c
	/// @returns momentum in MeV/c
	inline ROOT::Math::XYZVector MomentumVector() const { return momentum_ * direction_; }


private:
	// atomic number
	int z_;
	// atomic mass
	int a_;
	// charge, 0 for ion
	int q_;
	// mass in MeV/c^2
	double mass_;
	// excitation energy in MeV
	double excitation_;
	// energy in MeV
	double energy_;
	// kinetic energy in MeV
	double kinetic_;
	// momentum value in MeV/c
	double momentum_;
	// momentum in MeV/c
	ROOT::Math::XYZVector direction_;
};


/// @brief get momentum from kinematic energy
/// @param[in] mass mass of particle
/// @param[in] kinetic kinetic energy of particle
/// @returns momentum of particle
inline double MomentumFromKinetic(double mass, double kinetic) {
	return sqrt((2.0 * mass + kinetic) * kinetic);
}


/// @brief get kinetic energy from momentum
/// @param[in] mass mass of particle
/// @param[in] momentum momentum of particle
/// @returns kinetic energy of particle
inline double KineticFromMomentum(double mass, double momentum) {
	return sqrt(pow(momentum, 2.0) + pow(mass, 2.0)) - mass;
}


/// @brief create gas material
/// @param compound compound in catima format, {a, z, stoichiometric}
/// @param pressure pressure in Torr
/// @param length length in meters(m)
/// @returns gas material in catima format
/// @code{.cpp}
/// catima::Material h2 = GasMaterial({{1,1,2}}, 200.0, 1.0);
/// @endcode
catima::Material GasMaterial(
	std::initializer_list<std::array<double,3>> compound,
	double pressure,
	double length
);


/// @brief create solid material
/// @param compound compound in catima format, {a, z, stoichiomatric}
/// @param thickness thickness in mg/cm^2
/// @returns solid material in catima format
/// Example usage:
/// @code{.cpp}
/// catima::Material CD2 = SolidMaterial({{12,6,1}, {2,1,2}}, 9.53e-3);
/// @endcode
catima::Material SolidMaterial(
	std::initializer_list<std::array<double,3>> compound,
	double thickness
);


using ParticleList = std::vector<Particle>;



/// @brief simulate beam scatter on target
/// @param[in] beam beam particle
/// @param[in] target target particle
/// @param[in] angle scattering angle of beam in radian in center of mass coordinate
/// @param[in] beam_excitation beam excitation energy in MeV
/// @param[in] target_excitation target excitation energy in MeV
/// @returns list of scattered particles
///		index 0: scattered beam, index 1: scattered target
ParticleList Scatter(
	const Particle &beam,
	const Particle &target,
	double angle,
	double beam_excitation = 0.0,
	double target_excitation = 0.0
);


};	//	namspace gong

#endif	// __PARTICLE_H__