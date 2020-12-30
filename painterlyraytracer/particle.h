#ifndef PARTICLE_H
#define PARTICLE_H

#include "vec.h"

class PaintParticle {
public:
    float x, y;
    Vec3 color, direction, normal;
    float depth;
    int size, edge;
    PaintParticle(float x0, float y0, Vec3 &c, Vec3 &dir, Vec3 &n, float dep, int s, int e) {
        x = x0;
        y = y0;
        c = color;
        direction = dir;
        normal = n;
        depth = dep;
        size = s;
        edge = e;
    }
};

class StrokeParticle {
public:
    float x, y;
    Vec3 color;
    int size;
    StrokeParticle(float x0, float y0, Vec3 &c, int s) {
        x = x0;
        y = y0;
        c = color;
        size = s;
    }
};

#endif