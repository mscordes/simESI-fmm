#include "Core.h"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <random>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>
#include <memory>

// Ensure proper geometry of newly formed H3O+ and water

namespace Core {

    float convertToRadians(const float& degrees) {
        return degrees * M_PI / 180.0f;
    }

    // Angle betwen three points in degrees
    float getAngle(const std::array<float, 3>& p1, const std::array<float, 3>& p2, const std::array<float, 3>& p3) {
        // Calculate vectors u and v
        std::array<float, 3> u = subtractVectors(p1, p2);
        float unorm = std::sqrt(dotProduct(u, u));

        std::array<float, 3> v = subtractVectors(p2, p3);
        float vnorm = std::sqrt(dotProduct(v, v));

        return std::acos(dotProduct(u, v) / (unorm * vnorm)) * 57.2958f;
    }

    std::array<float, 3> addVectors(const std::array<float, 3>& u, const std::array<float, 3>& v) {
        return { u[0] + v[0], u[1] + v[1], u[2] + v[2] };
    }

	std::array<float, 3> subtractVectors(const std::array<float, 3>& u, const std::array<float, 3>& v) {
		return { u[0] - v[0], u[1] - v[1], u[2] - v[2] };
	}

    std::array<float, 3> scaleVector(const std::array<float, 3>& vec, const float& scalar) {
        return { vec[0] * scalar,  vec[1] * scalar,  vec[2] * scalar };
    }

    std::array<float, 3> crossProduct(const std::array<float, 3>& u, const std::array<float, 3>& v) {

        return {
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0]
        };
    }

    // Normalize a vector
	std::array<float, 3> normalizeVector(const std::array<float, 3>& vec) {

		float norm = std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
		if (norm == 0) {
			return { 0.0f, 0.0f, 0.0f };
		}

		return { vec[0] / norm, vec[1] / norm, vec[2] / norm };
	}

	// Dot product for 3x3 (rotation) matrix with 3D vector
    float dotProduct(const std::array<float, 3>& u, const std::array<float, 3>& v) {

        float dot = 0.0f;
        for (size_t i = 0; i < u.size(); ++i) {
            dot += u[i] * v[i];
        }
        return dot;
    }

    // Function to calculate the unit normal vector to a plane defined by four points
    std::array<float, 3> planeNormalVector(const std::array<float, 3>& p1, const std::array<float, 3>& p2,
        const std::array<float, 3>& p3, const std::array<float, 3>& p4) {

        // Calculate vectors u and v
		std::array<float, 3> u = subtractVectors(p1, p2);
		std::array<float, 3> v = subtractVectors(p3, p4);

        // Compute the cross product of u and v
        std::array<float, 3> normal = crossProduct(u, v);

        // Compute the norm of the normal vector
        float normal_norm = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

        // Normalize the normal vector to get the unit normal
        if (normal_norm == 0) {
            std::cerr << "The input vectors are parallel or invalid; normal vector cannot be computed.";
            exit(1);
        }

        return { normal[0] / normal_norm, normal[1] / normal_norm, normal[2] / normal_norm };
    }

    // Defines plane of rotation between three points and rotates p3 to desired angle.
    std::array<float, 3> rotateBond(const std::array<float, 3>& p1, const std::array<float, 3>& p2,
        std::array<float, 3>& p3, const float& theta) {

        // Define the plane normal vector
        std::array<float, 3> k = planeNormalVector(p1, p2, p2, p3);

        // Calculate vectors u and v
        std::array<float, 3> u = subtractVectors(p2, p1);
        float unorm = std::sqrt(dotProduct(u, u));

        std::array<float, 3> v = subtractVectors(p2, p3);
        float vnorm = std::sqrt(dotProduct(v, v));

        // Calculate old and new angles
        float theta_old = std::acos(dotProduct(u, v) / (unorm * vnorm));
        float theta_new = convertToRadians(theta);
        float theta_final = theta_old - theta_new;

        // Rodriguez rotation formula
        std::array<float, 3> v_rot = addVectors(
            scaleVector(v, std::cos(theta_final)),
            addVectors(
                scaleVector(crossProduct(k, v), std::sin(theta_final)),
                scaleVector(k, dotProduct(k, v) * (1 - std::cos(theta_final)))
            )
        );

        // Calculate the new position of p3
        std::array<float, 3> new_p3 = subtractVectors(p2, v_rot);
        return new_p3;
    }


    // Function to calculate a 3x3 rotation matrix
    std::vector<std::array<float, 3>> rotationMatrix(const std::array<float, 3>&axis, const float& theta) {

        // Normalize the axis vector
        float axis_norm = std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
        if (axis_norm == 0) {
            std::cerr << "Axis vector cannot be zero." << std::endl;
            exit(1);
        }

        std::array<float, 3> normalized_axis = { axis[0] / axis_norm, axis[1] / axis_norm, axis[2] / axis_norm };

        // Precompute terms for the quaternion representation
        float a = std::cos(theta / 2.0);
        float b = -normalized_axis[0] * std::sin(theta / 2.0f);
        float c = -normalized_axis[1] * std::sin(theta / 2.0f);
        float d = -normalized_axis[2] * std::sin(theta / 2.0f);

        // Calculate intermediate terms
        float aa = a * a, bb = b * b, cc = c * c, dd = d * d;
        float bc = b * c, ad = a * d, ac = a * c, ab = a * b, bd = b * d, cd = c * d;

        // Build and return the rotation matrix
        return {
            {aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)},
            {2.0f * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)},
            {2.0f * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc}
        };
    }

    std::array<float, 3> dotProductforRot(const std::vector<std::array<float, 3>>& rotationMatrix, const std::array<float, 3>& v) {
		std::array<float, 3> dot = { 0.0f, 0.0f, 0.0f };

        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                dot[i] += rotationMatrix[i][j] * v[j]; 
            }
        }
        return dot;
    }

    // Proper geometry of newly formed H3O+ 
	void H3O_transform(std::vector<std::shared_ptr<Atom>>& atoms) {
        std::vector<std::array<float, 3>> coords;
		for (const auto& atom : atoms) {
			coords.push_back(atom->coord);
		}

        // Make sure new H-bond lenght 0.09686 nm
        coords[4] = setBondLength(coords[0], coords[4], 0.09686f);

        // Make MOH bond 74.4* where M is the virtual site. Keep atoms in same initial plane, only move M atom.
		coords[3] = rotateBond(coords[4], coords[0], coords[3], 74.4f);
        coords[3] = setBondLength(coords[0], coords[3], 0.015f);

        // Delete two non-exchanged H's and replace by rotating exchanged H 120 degrees 2x
		std::array<float, 3> v = subtractVectors(coords[0], coords[4]);
		std::array<float, 3> axis = subtractVectors(coords[0], coords[3]);
        std::array<float, 3> k = normalizeVector(axis);
        coords[1] = subtractVectors(coords[0], dotProductforRot(rotationMatrix(axis, convertToRadians(120.0f)), v));
        coords[2] = subtractVectors(coords[0], dotProductforRot(rotationMatrix(axis, convertToRadians(240.0f)), v));
		coords[1] = setBondLength(coords[0], coords[1], 0.09686f);
		coords[2] = setBondLength(coords[0], coords[2], 0.09686f);

		// Set final coordinates
        for (size_t i = 0; i < atoms.size(); ++i) {
			atoms[i]->coord = coords[i];
        }
	}

    // Proper geometry of newly formed water
    void water_transform(std::vector<std::shared_ptr<Atom>>& atoms) {
        std::vector<std::array<float, 3>> coords;
        for (const auto& atom : atoms) {
            coords.push_back(atom->coord);
        }

        // Set length of O-H bonds to 0.09686 nm
        coords[1] = setBondLength(coords[0], coords[1], 0.09686f); // HW1
        coords[2] = setBondLength(coords[0], coords[2], 0.09686f); // HW2

        // Defines plane of HOH and sets virtual site bisecting two H's in plane at bond length of 0.015nm
        std::array<float, 3> u = subtractVectors(coords[1], coords[0]);
		std::array<float, 3> v = subtractVectors(coords[2], coords[0]);
		std::array<float, 3> un = normalizeVector(u);
		std::array<float, 3> vn = normalizeVector(v);
		std::array<float, 3> bisector = addVectors(un, addVectors(vn, coords[0]));
		coords[3] = setBondLength(coords[0], bisector, 0.015f); 
        
        // Sets water bond angle to 104.5
		coords[1] = rotateBond(coords[2], coords[0], coords[1], 104.5f);
		coords[1] = setBondLength(coords[0], coords[1], 0.09686f);

        // Set final coordinates
        for (size_t i = 0; i < atoms.size(); ++i) {
            atoms[i]->coord = coords[i];
        }
    }

}