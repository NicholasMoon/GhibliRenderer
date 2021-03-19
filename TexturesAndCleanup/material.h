// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// material.h

#ifndef MATERIAL_H
#define MATERIAL_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <limits>

#include "vec3.h"

class material {
 public:
  material();
  material(vec3 &shininess, vec3 &transparency, double ior, double roughness, double eccentricity);
  material(vec3 &diffuse_color, vec3 &ambient, vec3 &diffuse, vec3 &shininess, vec3 &transparency, double ior, double roughness, double eccentricity, std::string diffuse_map);
  ~material();
  vec3 diffuse_color;
  vec3 ambient;
  vec3 diffuse;
  vec3 shininess;
  vec3 transparency;
  double ior;
  double roughness;
  double eccentricity;
  std::string diffuse_map;
 private:
  
};

#endif
