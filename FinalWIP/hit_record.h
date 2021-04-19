// scene.h

#ifndef HITRECORD_H
#define HITRECORD_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <limits>
#include <random>

class object;
class vec3;
class ray;
#include "vec3.h"
#include "ray.h"
#include "object.h"
#include "scene.h"

class HitRecord {
 public:
  HitRecord(int x, int y);
  void setAttributes(double distance, vec3 &hit_point, vec3 &hit_normal);
  ~HitRecord();

  vec3 hit_point, primary_hit_normal, hit_normal;
  double distance;
  int primary_ray;
  int x, y;
  int lastObject, primary_objID, object_type;
  int shadowed;
  std::vector<int> hit_list;
 private:
  
};

#endif
