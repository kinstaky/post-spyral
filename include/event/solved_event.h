#ifndef __SOLVED_EVENT_H__
#define __SOLVED_EVENT_H__

#include <TTree.h>

struct SolvedEvent {
	int run, event;
	double vx, vy, vz;
	double brho, kinetic;
	double polar, azimuthal;
};

/// @brief Setup input branches for tree
/// @param[in] tree pointer to input tree
/// @param[in] event solved event
/// @returns 0 if success
int SetupInput(TTree* tree, SolvedEvent &event);


/// @brief Setup output branches for tree
/// @param[in] tree pointer to output tree
/// @param[in] event solved event
/// @returns 0 if success
int SetupOutput(TTree* tree, SolvedEvent &event);

#endif	// 	__SOLVED_EVENT_H__