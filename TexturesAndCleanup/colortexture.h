#ifndef COLORTEXTURE_H
#define COLORTEXTURE_H

#include "texture.h"

class color_texture : public texture {
 public:
  color_texture() {};
  color_texture(vec3 &color) : color(color) {};
  ~color_texture();
  vec3 getColor(double u, double v, vec3 &hit_point);
  vec3 color;
 private:
 
};

#endif
