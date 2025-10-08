#pragma once
#include <chrono>

using namespace std;

inline auto now() { 
	return chrono::steady_clock::now(); 
}

inline double seconds(auto&& duration) {
	return chrono::duration<double>(duration).count();
}