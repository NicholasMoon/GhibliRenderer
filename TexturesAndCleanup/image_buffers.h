// scene.h

#ifndef IMAGEBUFFERS_H
#define IMAGEBUFFERS_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <limits>
#include <random>

class ImageBuffers {
 public:
  ImageBuffers(int dimension);
  ~ImageBuffers();


  double *depthMap;
  double *edgeMap;
  int *objectTypeMap;
  int *objectBoundaryMap;
  int *objectIDMap;
  int *shadowMap;
  float *paintMap; // paint density values
 private:
  
};

#endif
