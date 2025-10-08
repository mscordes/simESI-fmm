#include "MathUtils.h"
#include "Constants.h"
#include "RNG.h"

#include <cmath>
#include <numeric>

double stdev(const vector<double>& v) {
    if (v.empty()) return 0.0;
    const double mean = accumulate(v.begin(), v.end(), 0.0) / v.size();
    double sum = 0.0;

    for (double val : v) {
        const double diff = val - mean;
        sum += diff * diff;
    }
    return sqrt(sum / v.size());
}

inline double distance(const Vec3D& a, const Vec3D& b) noexcept {
	const double dx = a.x - b.x;
	const double dy = a.y - b.y;
	const double dz = a.z - b.z;
	return sqrt(dx * dx + dy * dy + dz * dz);
}

vector<vector<double>> cdist(const vector<Vec3D>& v) {
    size_t size = v.size();
    vector<vector<double>> distances(size, vector<double>(size, 0.0));

    for (size_t i = 0; i < size; ++i)
        for (size_t j = i + 1; j < size; ++j)
            distances[i][j] = distances[j][i] = v[i].distance(v[j]);

    return distances;
}

double getAngle(const Vec3D& c1, const Vec3D& c2, const Vec3D& c3) {
	Vec3D u = c1 - c2;
	double unorm = u.magnitude();
	Vec3D v = c2 - c3;
	double vnorm = v.magnitude();
	return acos(u.dot(v) / (unorm * vnorm)) * 57.2958;
}

Vec3D setBondLength(const Vec3D& ref, const Vec3D& coord, double l) {
	return ref + (coord - ref).normalize() * l;
}

Vec3D getCenter(const vector<Vec3D>& coords) {
	Vec3D sum;
	for (const auto& c : coords)
		sum = sum + c;
	return sum / coords.size();
}

// Sample n velocities from a Maxwell-Boltzmann distribution in nm/ps
vector<Vec3D> sampleMaxwell(int n, double temperature, char element) {
	double mass = Constants::atomMasses.at(element);
	double sigma = sqrt(((Constants::R / 1000) * temperature) / mass);
	normal_distribution<double> dist(0.0, sigma);

	vector<Vec3D> velocities;
	velocities.reserve(n);

	for (int i = 0; i < n; ++i)
		velocities.emplace_back(Vec3D{ dist(RNG::gen), dist(RNG::gen), dist(RNG::gen) });

	return velocities;
}

double convertToRadians(double degrees) {
    return degrees * Constants::pi / 180.0;
}

Vec3D planeNormalVector(const Vec3D& p1, const Vec3D& p2, const Vec3D& p3) {
    Vec3D u = p2 - p1;
    Vec3D v = p3 - p1;
    return u.cross(v).normalize();
}

Vec3D rotateBond(const Vec3D& p1, const Vec3D& p2, Vec3D& p3, double theta) {
    Vec3D axis = planeNormalVector(p1, p2, p3);
    Vec3D u = p2 - p1;
    Vec3D v = p2 - p3;

    double thetaOld = acos(clamp(u.dot(v) / (u.magnitude() * v.magnitude()), -1.0, 1.0));
    double thetaNew = convertToRadians(theta);
    double dTheta = thetaOld - thetaNew;

    // Rodrigues rotation
    Vec3D v_rot = v * cos(dTheta)
        + axis.cross(v) * sin(dTheta)
        + axis * axis.dot(v) * (1 - cos(dTheta));

    return p2 + v_rot; // new p3
}

array<Vec3D, 3> rotationMatrix(const Vec3D& axis, double theta) {
    Vec3D n = axis.normalize();
    double half = theta * 0.5;
    double s = sin(half);
    double a = cos(half);
    double b = n.x * s;
    double c = n.y * s;
    double d = n.z * s;

    double aa = a * a, bb = b * b, cc = c * c, dd = d * d;
    double bc = b * c, ad = a * d, ac = a * c, ab = a * b, bd = b * d, cd = c * d;

    return {
        Vec3D{aa + bb - cc - dd, 2 * (bc + ad),     2 * (bd - ac)},
        Vec3D{2 * (bc - ad),      aa + cc - bb - dd, 2 * (cd + ab)},
        Vec3D{2 * (bd + ac),      2 * (cd - ab),     aa + dd - bb - cc}
    };
}

Vec3D applyRotationMatrix(const array<Vec3D, 3>& R, const Vec3D& v) {
    return Vec3D{
        R[0].x * v.x + R[0].y * v.y + R[0].z * v.z,
        R[1].x * v.x + R[1].y * v.y + R[1].z * v.z,
        R[2].x * v.x + R[2].y * v.y + R[2].z * v.z
    };
}

Vec3D bisector(const Vec3D& a, const Vec3D& b) {
    Vec3D a_norm = a.normalize();
    Vec3D b_norm = b.normalize();
    Vec3D bis = a_norm + b_norm;
    return bis.normalize();
}

// Metropolis Criterion 
bool MCMC(double E, double temp) {
    if (E > 100) return false;
    else if (E < 0) return true;

    double RT = (Constants::R / 1000.0) * temp; // kJ/mol
    double boltzmann = exp(-E / RT);
    if (boltzmann > 1.0) boltzmann = 1.0;

    static uniform_real_distribution<double> d(0.0, 1.0);
    return d(RNG::gen) < boltzmann;
}

double computeCoulomb(
    const vector<vector<double>>& distances, 
    const vector<double>& charges) 
{
    static const double ke = Constants::ke;
    double e = 0.0;
    size_t size = charges.size();
    for (size_t i = 0; i < size; ++i)
        for (size_t j = i + 1; j < size; ++j)
            e += (charges[i] * charges[j]) / distances[i][j];
    return e * ke; // kJ / mol
}
