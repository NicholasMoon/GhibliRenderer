// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// object.h

#ifndef OBJECT_H
#define OBJECT_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

class ray;
class AABB;
#include "vec3.h"
#include "material.h"
#include "ray.h"
#include "light.h"
#include "AABB.h"

class object {
 public:
  virtual bool hit(ray *incoming_ray, std::vector<object*> &objects, std::vector<light*> &lights, double color[4], double &distance) { std::cout << "wrong hit" << std::endl; return 0; };
  virtual bool shadowHit(ray *incoming_ray, light* target_light, double &distance) { std::cout << "wrong hit" << std::endl; return 0; };
  virtual vec3 getColor() { return vec3(0,0,0); };
  virtual vec3 getNormal(double x, double y, double z, int flat) { return vec3(0,0,0); };
  virtual bool in_bounding_box(AABB *bounding_box) { return true; };
  virtual void updateWorldBoundaries(vec3 &min_coordinates, vec3 &max_coordinates) { return; };
  material *mat;
  int objectID;
  int object_type;
 private:
  
};

#endif

