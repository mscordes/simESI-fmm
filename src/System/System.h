#pragma once
#include "State.h"
#include "Parameters.h"

struct Writers;  // Forward declaration

struct RunFlags {
	bool createRunFile : 1;
	bool exchanges : 1;
	bool protExchanges : 1;

	void reset() {
		createRunFile = false;
		exchanges = false;
		protExchanges = false;
	}
};

struct System {
	Parameters params;
	State state;

	System() = default;

	void run();	
	void initialize();
	void doStep(
		int step, 
		Writers& writers, 
		RunFlags& flags, 
		const chrono::steady_clock::time_point& tInit, 
		int& numRestarts, 
		int& numContSteps);
};