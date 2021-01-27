// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// tri.h

#ifndef TRI_H
#define TRI_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <limits>

#include "object.h"
#include "vec3.h"
#include "vertex.h"
#include "ray.h"
#include "light.h"

class tri : public object {
 public:
  tri(vertex *v1, vertex *v2, vertex *v3, double c[3]);
  tri(vertex *v1, vertex *v2, vertex *v3, double c[3], material *mat, int objectID);
  tri(vertex *v1, vertex *v2, vertex *v3, double c[3], material *mat, int objectID, int object_type);
  bool hit(ray *incoming_ray, std::vector<object*> &objects, std::vector<light*> &lights, double color[4], double &distance);
  bool shadowHit(ray *incoming_ray, light* target_light, double &distance);
  vec3 getColor();
  vec3 getNormal(double x, double y, double z, int flat);
  ~tri();
  vertex v1, v2, v3;
  vec3 n;
  vec3 edge1, edge2;
  vec3 e2, e3;
  vec3 barycentric;
  double c[3];
 private:
  
};

#endif

