#pragma once
#include "State.h"
#include "Parameters.h"
#include "System.h"

struct WaterQuad { Vec3D O, H1, H2, M; };

struct WaterInfo {
	vector<shared_ptr<Residue>>	waters;
	vector<Vec3D>				Ocoords;
	vector<Vec3D>				Hcoords;
	vector<WaterQuad>			coords;
};

vector<Vec3D> getAllCoords(const State& state);
double getProteinMass(const State& state);
vector<Vec3D> getProteinCoords(const State& state);
vector<Vec3D> getProteinCarbonCoords(const State& state);
vector<Vec3D> getProteinSkeletonCoords(const State& state);
shared_ptr<Atom> findClosestWaterO(double& dist, const Vec3D& ref, const WaterInfo& waterInfo);
shared_ptr<Atom> findClosestWaterH(double& dist, const Vec3D& ref, const WaterInfo& waterInfo);
void centerState(State& state, const Vec3D& boxD);
void centerState(State& state, const double boxScalar);
Vec3D centerSmallBox(State& state);
void carveDroplet(State& state, const double envelope = Parameters::Get().getDropletSize());
shared_ptr<Residue> parseResidueGRO(const string& groFile);
shared_ptr<Residue> cloneResidue(const Residue& src);
void insertMolecIntoDroplet(
	const string& resType, 
	const unordered_map<string, int>& tooSeed,
	State& state, 
	const vector<Vec3D>& protCoords, 
	vector<Vec3D>& coords, 
	const Vec3D& boxD, 
	double minSep = 0.30);
void insertMolec(const Vec3D& newLoc, const shared_ptr<Residue>& templateResidue, State& state);
void fixDisulfides(State& state, RunFlags& flags);
WaterInfo getWaterInfo(const State& state);
unordered_map<string, vector<Vec3D>> getSoluteCoords(const State& state);
unordered_set<shared_ptr<Residue>> getGaseousWaters(const State& state);