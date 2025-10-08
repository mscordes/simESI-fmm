#include "Atom.h"
#include "Residue.h"
#include <iostream>

#include <utility>

Atom::Atom(
    string name_,
    char element_,
    Vec3D coord_,
    Vec3D velocity_,
    shared_ptr<Residue> parent_)
    : name(move(name_)),
    element(element_),
    coord(move(coord_)),
    velocity(move(velocity_)),
    parent(parent_)
{
}

void Atom::print() const {
	auto p = parent.lock();
    cout << "  Atom: " << name << ", Residue: " << p->ID << ", Chain: " << p->chain
        << ", Coord: [" << coord.x << "," << coord.y << "," << coord.z << "]"
        << ", Velocity: [" << velocity.x << "," << velocity.y << "," << velocity.z << "]\n";
}