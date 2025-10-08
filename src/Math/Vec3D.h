#pragma once
#include <cmath>

using namespace std;

struct Vec3D 
{
    double x, y, z;

    Vec3D() : x(0), y(0), z(0) {}
    Vec3D(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec3D operator+(const Vec3D& other) const {
        return Vec3D(x + other.x, y + other.y, z + other.z);
    }

    Vec3D operator+(double scalar) const {
        return Vec3D(x + scalar, y + scalar, z + scalar);
    }

    Vec3D operator-(const Vec3D& other) const {
        return Vec3D(x - other.x, y - other.y, z - other.z);
    }

    Vec3D operator-(double scalar) const {
        return Vec3D(x - scalar, y - scalar, z - scalar);
    }

    Vec3D operator*(double scalar) const {
        return Vec3D(x * scalar, y * scalar, z * scalar);
    }

    Vec3D operator/(double scalar) const {
        return Vec3D(x / scalar, y / scalar, z / scalar);
    }

    double dot(const Vec3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3D cross(const Vec3D& other) const {
        return Vec3D(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }

    double magnitude() const {
        return sqrt(x * x + y * y + z * z);
    }

    Vec3D normalize() const {
        double mag = magnitude();
        if (mag == 0) return Vec3D(0, 0, 0);
        return Vec3D(x / mag, y / mag, z / mag);
    }

    inline double distance_sq(const Vec3D& other) const noexcept {
        const double dx = x - other.x;
        const double dy = y - other.y;
        const double dz = z - other.z;
        return dx * dx + dy * dy + dz * dz;
    }

    inline double distance(const Vec3D& other) const noexcept {
        return sqrt(distance_sq(other));
    }
};