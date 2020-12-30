#ifndef RAY_H
#define RAY_H

#include "vec.h"

class Ray {
public:
    Vec3 origin;
    Vec3 direction;
    Ray(Vec3 o = 0, Vec3 d = 0) { origin = o; direction = d; }

    Vec3 point_at(float t) const { return origin + direction * t; }
};

#endif