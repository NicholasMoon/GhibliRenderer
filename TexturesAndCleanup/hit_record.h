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

#include "vec3.h"
#include "object.h"
#include "light.h"

class HitRecord {
 public:
  HitRecord(char** inputFile);
  ~HitRecord();

  double color[4];
 private:
  
};

#endif
