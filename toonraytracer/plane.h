// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// plane.h

#ifndef PLANE_H
#define PLANE_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <limits>

#include "object.h"
#include "vec3.h"
#include "ray.h"
#include "light.h"

class plane : public object {
 public:
  plane(double A, double B, double C, double D, double color[3]);
  plane(double A, double B, double C, double D, double color[3], material *mat, int objectID);
  bool hit(ray *incoming_ray, std::vector<object*> &objects, std::vector<light*> &lights, double color[4], double &distance);
  bool shadowHit(ray *incoming_ray, light* target_light, double &distance);
  vec3 getColor();
  vec3 getNormal(double x, double y, double z);
  ~plane();
  double A, B, C, D;
  vec3 p;
  vec3 n;
  double c[3];
 private:
  
};

#endif

