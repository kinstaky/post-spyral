#include "include/event/solved_event.h"

int SetupInput(TTree* tree, SolvedEvent &event) {
	tree->SetBranchAddress("origin_run", &event.run);
	tree->SetBranchAddress("origin_event", &event.event);
	tree->SetBranchAddress("vertex_x", &event.vx);
	tree->SetBranchAddress("vertex_y", &event.vy);
	tree->SetBranchAddress("vertex_z", &event.vz);
	tree->SetBranchAddress("brho", &event.brho);
	tree->SetBranchAddress("kinetic_energy", &event.kinetic);
	tree->SetBranchAddress("polar", &event.polar);
	tree->SetBranchAddress("azimuthal", &event.azimuthal);
	return 0;
}


int SetupOutput(TTree* tree, SolvedEvent &event) {
	tree->Branch("origin_run", &event.run, "run/I");
	tree->Branch("origin_event", &event.event, "evt/I");
	tree->Branch("vertex_x", &event.vx, "vx/D");
	tree->Branch("vertex_y", &event.vy, "vy/D");
	tree->Branch("vertex_z", &event.vz, "vz/D");
	tree->Branch("brho", &event.brho, "brho/D");
	tree->Branch("kinetic_energy", &event.kinetic, "ek/D");
	tree->Branch("polar", &event.polar, "theta/D");
	tree->Branch("azimuthal", &event.azimuthal, "phi/D");
	return 0;
}
