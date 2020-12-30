// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// light.h

#ifndef LIGHT_H
#define LIGHT_H

#include <math.h>

class light {
 public:
  light(int type, double x, double y, double z, double c[3]);
  ~light();
  double x, y, z;
  double c[3];
  int type;
 private:
  
};

#endif

