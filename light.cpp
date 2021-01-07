// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// light.cpp

#include "light.h"


light::light(int type, double x, double y, double z, double c[3]) {
  this->type = type;
  if (type == 1) {
	  double magnitude = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
	  this->x = x / magnitude;
	  this->y = y / magnitude;
	  this->z = z / magnitude;
  }
  else {
	this->x = x;
	this->y = y;
	this->z = z;
  }
  this->c[0] = c[0];
  this->c[1] = c[1];
  this->c[2] = c[2];
}

light::~light() {

}
