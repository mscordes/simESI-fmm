#pragma once
#include "Vec3D.h"
#include <memory>
#include <string>

using namespace std;

struct Residue; // Forward declaration of Residue

struct Atom {
    string name;  
    char   element;
    Vec3D  coord;             
    Vec3D  velocity;   
    weak_ptr<Residue> parent;   

    Atom() = default;

    Atom(string name,
        char element,
        Vec3D coord,
        Vec3D velocity,
        shared_ptr<Residue> parent);

    void print() const;
};