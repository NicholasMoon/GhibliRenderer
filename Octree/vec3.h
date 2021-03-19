// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// vec3.h

#ifndef VEC3_H
#define VEC3_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

class vec3 {
 public:
  vec3(double x, double y, double z);
  vec3();
  vec3(const vec3 &v);
  vec3 subtract(const vec3 &v);
  void subtract(const vec3 &v, vec3 &result);
  void subtractAndSave(const vec3 &v);
  double dot(vec3 &v);
  vec3 cross(vec3 &v);
  void normalize();
  double magnitude();
  void reverse();
  double distance(vec3 &v);
  void printVec();
  ~vec3();
  double x, y, z;
 private:
  
};

#endif

