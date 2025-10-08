#pragma once
#include "Vec3D.h"
#include <array>
#include <vector>

double stdev(const vector<double>& v);
inline double distance(const Vec3D& a, const Vec3D& b) noexcept;
vector<vector<double>> cdist(const vector<Vec3D>& v);
double getAngle(const Vec3D& c1, const Vec3D& c2, const Vec3D& c3);
Vec3D setBondLength(const Vec3D& ref, const Vec3D& coord, double l);
Vec3D getCenter(const vector<Vec3D>& coords);
vector<Vec3D> sampleMaxwell(int n, double temperature, char element);
double convertToRadians(double degrees);
Vec3D planeNormalVector(const Vec3D& p1, const Vec3D& p2, const Vec3D& p3);
Vec3D rotateBond(const Vec3D& p1, const Vec3D& p2, Vec3D& p3, double theta);
array<Vec3D, 3> rotationMatrix(const Vec3D& axis, double theta);
Vec3D applyRotationMatrix(const array<Vec3D, 3>& R, const Vec3D& v);
Vec3D bisector(const Vec3D& a, const Vec3D& b);
bool MCMC(double E, double temp);
double computeCoulomb(const vector<vector<double>>& distances, const vector<double>& charges);