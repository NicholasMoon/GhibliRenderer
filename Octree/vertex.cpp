// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// vertex.cpp

#include "vertex.h"


vertex::vertex(vec3 &xyz, vec3 &uv, vec3 &normal, vec3 &color) {
  this->xyz = xyz;
  this->uv = uv;
  this->normal = normal;
  this->color = color;
}

vertex::vertex() {
  vec3 origin(0,0,0);
  this->xyz = origin;
  this->uv = origin;
  this->normal = origin;
  this->color = origin;
}

vertex::~vertex() {

}
