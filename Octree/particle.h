#ifndef PARTICLE_H
#define PARTICLE_H

#include "vec3.h"

class PaintParticle {
public:
    float x, y;
    vec3 color, direction, normal;
    float depth;
    int size, edge;
    PaintParticle(float x0, float y0, vec3 &c, vec3 &dir, vec3 &n, float dep, int s, int e) {
        x = x0;
        y = y0;
        color = c;
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
    vec3 color;
    int size;
    StrokeParticle(float x0, float y0, vec3 &c, int s) {
        x = x0;
        y = y0;
        color = c;
        size = s;
    }
};

#endif
