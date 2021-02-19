// octree.h

#ifndef AABB_H
#define AABB_H

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

class AABB {
 public:
  AABB(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z): min_coordinates(min_x, min_y, min_z), max_coordinates(max_x, max_y, max_z) { };
  bool intersect(ray *r);
  vec3 min_coordinates;
  vec3 max_coordinates;
  ~AABB();
 private:
  
};

#endif
