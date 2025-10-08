#include "Transforms.h"
#include "Atom.h"
#include "Constants.h"
#include "MathUtils.h"

#include <array>

// Proper geometry of newly formed H3O+ 
void hydroniumTransform(shared_ptr<Residue>& hydronium) {
    Vec3D&  O = hydronium->atoms[0]->coord;
    Vec3D& H1 = hydronium->atoms[1]->coord;
    Vec3D& H2 = hydronium->atoms[2]->coord;
    Vec3D&  M = hydronium->atoms[3]->coord;
    Vec3D& H3 = hydronium->atoms[4]->coord;

    H3 = setBondLength(O, H3, 0.09686);

    Vec3D center = getCenter({ O, H1, H2, H3 });
    Vec3D plane = planeNormalVector(H1, H2, H3);
    if ((center - O).dot(plane) > 0) plane = plane * -1;

    M = rotateBond(O, H3, plane, 74.6);
    M = setBondLength(O, M, 0.015);

    Vec3D v = H3 - O;               
    Vec3D axis = (M - O).normalize();  

    H1 = O + applyRotationMatrix(rotationMatrix(axis, convertToRadians(120.0)), v);
    H2 = O + applyRotationMatrix(rotationMatrix(axis, convertToRadians(240.0)), v);

    H1 = setBondLength(O, H1, 0.09686);
    H2 = setBondLength(O, H2, 0.09686);
}

// Proper geometry of newly formed water
void waterTransform(shared_ptr<Residue>& water) {
    Vec3D&  O = water->atoms[0]->coord;
    Vec3D& H1 = water->atoms[1]->coord;
    Vec3D& H2 = water->atoms[2]->coord;
    Vec3D&  M = water->atoms[3]->coord;

    H1 = setBondLength(O, H1, 0.09572); 

    Vec3D center = getCenter({ O, H1, H2});
    Vec3D bs = bisector(O - H1, O - H2);
    if ((center - O).dot(bs) > 0) bs = bs * -1;

    M = rotateBond(O, H1, bs, 52.5);
    M = setBondLength(O, M, 0.01546);

    Vec3D v = H1 - O;
    Vec3D axis = (M - O).normalize();

    H2 = O + applyRotationMatrix(rotationMatrix(axis, convertToRadians(180.0)), v);
    H2 = setBondLength(O, H2, 0.09572);
}