// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// sphere.h

#ifndef SPHERE_H
#define SPHERE_H

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

class sphere : public object {
 public:
  sphere(double x, double y, double z, double r, double c[3]);
  sphere(double x, double y, double z, double r, double c[3], material *mat, int objectID);
  sphere(double x, double y, double z, double r, double c[3], material *mat, int objectID, int object_type);
  sphere(double x, double y, double z, double r, double c[3], double e[3], material *mat, int objectID, int object_type);
  bool hit(ray *incoming_ray, std::vector<object*> &objects, std::vector<light*> &lights, double color[4], double &distance);
  bool shadowHit(ray *incoming_ray, light* target_light, double &distance);
  vec3 getColor();
  vec3 getEmission();
  vec3 getNormal(double x, double y, double z, int flat);
  bool in_bounding_box(AABB *bounding_box);
  void updateWorldBoundaries(vec3 &min_coordinates, vec3 &max_coordinates);
  ~sphere();
  double x, y, z, r;
  double c[3];
  double e[3];

 private:
  
};

#endif

