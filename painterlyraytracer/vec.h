#ifndef VEC_H
#define VEC_H

#include <cstdio>
#include <math.h>

class Vec3 {
public:
    float x, y, z;
    Vec3(float x0 = 0, float y0 = 0, float z0 = 0) { x = x0; y = y0; z = z0; }
    Vec3 operator+(const Vec3 &v2) const { return Vec3(x + v2.x, y + v2.y, z + v2.z); }
    Vec3 operator-(const Vec3 &v2) const { return Vec3(x - v2.x, y - v2.y, z - v2.z); }
    Vec3 operator*(float s) const { return Vec3(x * s, y * s, z * s); }
    Vec3 operator/(float s) const { return Vec3(x / s, y / s, z / s); }

    float mag() const {
        float squared_sum = x * x + y * y + z * z;
        return sqrt(squared_sum);
    }

    Vec3 norm() const { 
        return *this/mag();
    }
    
    float dot(const Vec3 &v2) {
        return x * v2.x + y * v2.y + z * v2.z;
    }

    Vec3 cross(const Vec3 &v2) {
        float rx = (y * v2.z) - (z * v2.y);
        float ry = (z * v2.x) - (x * v2.z);
        float rz = (x * v2.y) - (y * v2.x);
        return Vec3(rx, ry, rz);
    }

    Vec3 multVbyV(const Vec3 &v2) {
        return Vec3(x * v2.x, y * v2.y, z * v2.z);
    }

    void print() {
        printf("(%f, %f, %f) \n", x, y, z);
    }
};

#endif