// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// vertex.h

#ifndef VERTEX_H
#define VERTEX_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include "vec3.h"

class vertex {
 public:
  vertex(vec3 &xyz, vec3 &normal, vec3 &color);
  vertex();
  ~vertex();
  vec3 xyz;
  vec3 normal;
  vec3 color;
 private:
  
};

#endif

