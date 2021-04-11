#ifndef IMAGETEXTURE_H
#define IMAGETEXTURE_H

#include "texture.h"

class image_texture : public texture {
 public:
  image_texture() {};
  image_texture(std::string file);
  ~image_texture();
  vec3 getColor(double u, double v, vec3 &hit_point);
  std::string file;
  CImg<float> image;
  int width, height;
 private:
  
};

#endif
