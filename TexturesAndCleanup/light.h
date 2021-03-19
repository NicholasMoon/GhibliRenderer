// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// light.h

#ifndef LIGHT_H
#define LIGHT_H

#include <math.h>
#include "vec3.h"

class light {
 public:
  light(int type, double x, double y, double z, double c[3]);
  light(int type, vec3 &p0, vec3 &p1, vec3 &p2, vec3 &p3, double c[3]);
  ~light();
  double x, y, z;
  vec3 p0;
  vec3 p1;
  vec3 p2;
  vec3 p3;
  vec3 light_center;
  double c[3];
  int type;
 private:
  
};

#endif

