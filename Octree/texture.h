#ifndef TEXTURE_H
#define TEXTURE_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include "CImg.h"
using namespace cimg_library;

#include "vec3.h"
#include "material.h"
#include "utilities.h"

class texture {
 public:
  virtual vec3 getColor(double u, double v, vec3 &hit_point) { return vec3(0,0,0); };
 private:
  
};

class color_texture : public texture {
 public:
  color_texture() {};
  color_texture(vec3 &color) : color(color) {};
  vec3 getColor(double u, double v, vec3 &hit_point);
  vec3 color;
 private:
 
};

class image_texture : public texture {
 public:
  image_texture() {};
  image_texture(std::string file);
  vec3 getColor(double u, double v, vec3 &hit_point);
  std::string file;
  CImg<float> image;
  int width, height;
 private:
  
};

#endif