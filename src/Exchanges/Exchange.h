#pragma once
#include "Atom.h"
#include "Clustering.h"
#include "CoordinateManip.h"
#include "State.h"

#include <ostream>
#include <optional>

struct Exchange {
	shared_ptr<Atom> hydrogen;
	shared_ptr<Atom> acceptor;
	int step;
	int hop;
	int clusterID;
	int nearWaters; 
	double ph;
	double temp;
	optional<double> h_pka;
	optional<double> h_gpb;
	optional<double> a_pka;
	optional<double> a_gpb;
	double energy;      
	bool isNterm{ false };
	bool isCterm{ false };
	bool isGrotthuss{ false };

	Exchange() = default;

	Exchange(
		const State& state,
		const shared_ptr<Atom>& hydrogen,
		const shared_ptr<Atom>& acceptor,
		int step, int hop,
  		const Cluster& cluster,
		const vector<Vec3D>& nonWaterCoords,
		const vector<double>& nonWaterCharges,
		const unordered_map<string, vector<Vec3D>>& sCoords,
		const unordered_map<string, vector<double>>& sCharges,
		const vector<WaterQuad>& waterCoords,
	  	double temp);

	Exchange( // For intramolecular (mobile proton) protein transfers
		const State& state,
		const shared_ptr<Atom>& hydrogen,
		const shared_ptr<Atom>& acceptor,
		double energy, 
		int step, 
		double temp);

	void print(ostream& os) const;
	double calcExchangeEnergy() const;
	double calcGrotthussEnergy(
		const State& state,
		const vector<Vec3D>& nonWaterCoords,
		const vector<double>& nonWaterCharges,
		const unordered_map<string, vector<Vec3D>>& sCoords,
		const unordered_map<string, vector<double>>& sCharges,
		const vector<WaterQuad>& waterCoords) const;
};