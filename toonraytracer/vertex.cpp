// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// vertex.cpp

#include "vertex.h"


vertex::vertex(vec3 &xyz, vec3 &normal, vec3 &color) {
  this->xyz = xyz;
  this->normal = normal;
  this->color = color;
}

vertex::vertex() {
  vec3 origin(0,0,0);
  this->xyz = origin;
  this->normal = origin;
  this->color = origin;
}

vertex::~vertex() {

}
