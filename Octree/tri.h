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

double Max_Double(double a, double b, double c);
double Min_Double(double a, double b, double c);

class tri : public object {
 public:
  tri(vertex *v1, vertex *v2, vertex *v3, double c[3]);
  tri(vertex *v1, vertex *v2, vertex *v3, double c[3], material *mat, int objectID);
  tri(vertex *v1, vertex *v2, vertex *v3, double c[3], material *mat, int objectID, int object_type);
  tri(vertex *v1, vertex *v2, vertex *v3, double c[3], double e[3], material *mat, int objectID, int object_type); 
  bool hit(ray *incoming_ray, std::vector<object*> &objects, std::vector<light*> &lights, double color[4], double &distance);
  bool shadowHit(ray *incoming_ray, light* target_light, double &distance);
  vec3 getColor();
  vec3 getEmission();
  vec3 getNormal(double x, double y, double z, int flat);
  bool SAT_projection(vec3 &axis, vec3 &center, vec3 &extents, vec3 &box_normal_x, vec3 &box_normal_y, vec3 &box_normal_z, vec3 &vc1, vec3 &vc2, vec3 &vc3);
  bool in_bounding_box(AABB *bounding_box);
  void updateWorldBoundaries(vec3 &min_coordinates, vec3 &max_coordinates);
  ~tri();
  vertex v1, v2, v3;
  vec3 n;
  vec3 edge1, edge2;
  vec3 e2, e3;
  vec3 barycentric;
  double c[3];
  double e[3];
 private:
  
};

#endif

