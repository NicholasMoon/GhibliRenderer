// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// ray.h

#ifndef RAY_H
#define RAY_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <limits>
#include <random>


class object;
#include "vec3.h"
#include "object.h"
#include "light.h"

class ray {
 public:
  ray(double ox, double oy, double oz, double dx, double dy, double dz) : origin(ox,oy,oz), direction(dx,dy,dz) { };
  bool cast(std::vector<object*> &objects, std::vector<light*> &lights, double color[4], int bounces, int lastObject, std::default_random_engine &generator, int primary_ray, int x, int y, int width, double *depthMap, int &primary_objID, vec3 &hit_normal, int &object_type, std::vector<int> &hit_list, int &shadowed, int flat, int light_samples, int indirect_samples, int indirect_bounces);
  double detect_edge(std::vector<object*> &objects, std::vector<light*> &lights, double color[4], int bounces, int lastObject, std::default_random_engine &generator, int &stencil_objID, int &object_type, std::vector<int> &hit_list);
  double castLight(std::vector<object*> &objects, light *targetLight, double distance);
  ~ray();
  vec3 origin;
  vec3 direction;
 private:
  
};

#endif