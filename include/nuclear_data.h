#ifndef __NUCLEAR_DATA_H__
#define __NUCLEAR_DATA_H__

#include <string>
#include <iostream>
#ifdef USE_HARDCODE_DATA
#include <map>
#endif

namespace gong {

struct NuclearData {
	int z;
	int a;
	double mass;
	std::string name;
};


inline std::istream &operator>>(std::istream &is, NuclearData &nuclear) {
	is >> nuclear.z >> nuclear.a >> nuclear.name >> nuclear.mass;
	return is;
}

inline std::ostream &operator<<(std::ostream &os, const NuclearData &nuclear) {
	os << nuclear.z << " " << nuclear.a << " " << nuclear.name << " " << nuclear.mass;
	return os;
}


/// @brief get nuclear data from atomic mass number and atomic number
/// @param[in] z atomic number
/// @param[in] a atomic mass number
/// @param[in] q charge
/// @returns nuclear data
NuclearData GetNuclear(int z, int a, int q = 0);


/// @brief get mass from atomic mass number and atomic number
/// @param[in] z atomic number
/// @param[in] a atomic mass number
/// @param[in] q charge
/// @returns mass in MeV/c^2
double GetMass(int z, int a, int q = 0);


/// @brief get mass in unit of atomic mass unit
/// @param[in] z atomic number
/// @param[in] a atomic mass number
/// @param[in] q charge
/// @returns mass in atomic mass unit
double GetMassInUnit(int z, int a, int q = 0);

#ifdef USE_HARDCODE_DATA
extern const std::map<std::pair<int, int>, NuclearData> nuclear_data_map;
#endif

}	// namespace gong

#endif // __NUCLEAR_DATA_H__
