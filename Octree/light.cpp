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

light::light(int type, vec3 &p0, vec3 &p1, vec3 &p2, vec3 &p3, double c[3]) {
	this->type = type;
	this->p0 = p0;
	this->p1 = p1;
	this->p2 = p2;
	this->p3 = p3;
	this->c[0] = c[0];
  	this->c[1] = c[1];
  	this->c[2] = c[2];
	double lcx,lcy,lcz = 0;
	vec3 lc((p0.x + p1.x + p2.x + p3.x) / 4, (p0.y + p1.y + p2.y + p3.y) / 4, (p0.z + p1.z + p2.z + p3.z) / 4);
	this->light_center = lc;
}

light::~light() {

}
