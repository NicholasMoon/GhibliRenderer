#ifndef LIGHTS_H
#define LIGHTS_H

#include "utilities.h"
#include "objects.h"
// #include <iostream>

class Light {
public:
    Vec3 position;
    Vec3 color;
    virtual Vec3 direction(const Intersection &intersection) const = 0;
    virtual double distance(const Intersection &intersection) const = 0;
    virtual double intensity(const Intersection &intersection) const = 0;
};

class Sun : public Light {
public:
    Sun(Vec3 p = 0, Vec3 c = 0) { position = p; color = c; }
    Vec3 direction(const Intersection &intersection) const override {
        return position.norm();
    }

    double distance(const Intersection &intersection) const override {
        return std::numeric_limits<double>::max();
    }

    double intensity(const Intersection &intersection) const override {
        return 1;
    }
};

class Bulb : public Light {
public:
    Bulb(Vec3 p = 0, Vec3 c = 0) { position = p; color = c; }
    Vec3 direction(const Intersection &intersection) const override {
        return (position - intersection.point).norm();
    }

    double distance(const Intersection &intersection) const override {
        double d = (position - intersection.point).mag();
        return d;
    }

    double intensity(const Intersection &intersection) const override {
        double d = distance(intersection);
        return 1 / (d*d);
    }
};

#endif