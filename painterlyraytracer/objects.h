#ifndef OBJECTS_H
#define OBJECTS_H

#include "utilities.h"

#define EPSILON 1e-7

class Intersection {
public:
    float distance;
    Vec3 point, normal;
    int edge;
    Intersection(float d = 0, Vec3 p = 0, Vec3 n = 0, int e = 0) { 
        distance = d; 
        point = p; 
        normal = n; 
        edge = e; 
    }

    void set_attributes(float t, const Ray &ray, Vec3 center) {
        distance = t;
        point = ray.point_at(t);
        normal = (point - center).norm();

        set_normal(ray);
        check_edge(ray);
    }

    void set_normal(const Ray &ray) {
        if (normal.dot(ray.direction) > 0) { normal = normal * -1; }
    }

    void check_edge(const Ray &ray) {
        if (normal.dot(ray.direction) < EPSILON) { 
            edge = 1; 
        } else {
            edge = 0;
        }
    }
};

class SceneObj {
public:
    int material, type;
    Vec3 color;
    Vec3 emission;
    virtual bool intersect(const Ray &ray, Intersection &intersection) const = 0;
};

class Sphere : public SceneObj {
public:
    Vec3 center;
    float radius;
    Sphere(Vec3 c = 0, float r = 0) { center = c; radius = r; }
    bool intersect(const Ray &ray, Intersection &intersection) const override {
        Vec3 co = center - ray.origin;
        float e = co.dot(ray.direction);
        if (e < 0) { return false; }

        float h2 = co.dot(co) - e*e;
        if (h2 > radius*radius) { return false; }

        float f = sqrt(radius*radius - h2);

        // calculate distance
        float t = 0;
        float t1 = e - f;
        float t2 = e + f;

        if (t1 < 0 && t2 < 0) { return false; }

        if (t1 < t2) {
            if (t1 < 0) {
                t = t2;
            } else {
                t = t1;
            }
        } else {
            t = t2;
        }

        // set intersection attributes (point, surface normal, etc.)
        intersection = Intersection();
        intersection.set_attributes(t, ray, center); // this won't work for triangles, planes

        return true;
    }
};

#endif