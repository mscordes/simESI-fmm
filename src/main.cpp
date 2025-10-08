#include "Core.h"		

#include <iostream>
#include <stdexcept>

using namespace std;

/******************************************************************************
 *                                                                            *
 *                                 simESI-fmm                                 *
 *                                                                            *
 *   Calculates the effect of atmospheric interaction on ionic temperature    *
 *   per original IICT while incorporating the effect of simultaneous         *
 *   collisions, and inelastic collisions with water.                         *
 *                                                                            *
 *   Copyright 2025 Elyssia S. Gallagher, Michael Cordes,                     *
 *   and Baylor University                                                    *
 *                                                                            *
 *   Created by Michael Cordes at the Gallagher Lab, Baylor University.       *
 *                                                                            *
 ******************************************************************************/

int main(int argc, char** argv) {
	try {
		run(argc, argv);
	}
	catch (const runtime_error& e) {
		if (string(e.what()) == "\n\nHelp requested")
			return 0;
		cerr << "\n\nRuntime error: " << e.what() << endl;
		return 1;
	}
	catch (...) {
		cerr << "\n\nUnknown error occurred." << endl;
		return 1;
	}
	return 0; 
}