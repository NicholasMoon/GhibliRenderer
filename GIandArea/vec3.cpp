// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// vec3.cpp

#include "vec3.h"


vec3::vec3(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

vec3::vec3() {
  this->x = 0;
  this->y = 0;
  this->z = 0;
}

vec3::vec3(const vec3 &v) {
  this->x = v.x;
  this->y = v.y;
  this->z = v.z;
}

vec3 vec3::subtract(const vec3 &v) {
  vec3 result(0,0,0);
  result.x = this->x - v.x;
  result.y = this->y - v.y;
  result.z = this->z - v.z;
  return result;
}

double vec3::dot(vec3 &v) {
  return this->x * v.x + this->y * v.y + this->z * v.z;
}

vec3 vec3::cross(vec3 &v) {
  return vec3(this->y * v.z - this->z * v.y, this->z * v.x - this->x * v.z, this->x * v.y - this->y * v.x);
}

void vec3::normalize() {
  double magnitude = this->magnitude();
  this->x /= magnitude;
  this->y /= magnitude;
  this->z /= magnitude;
}

double vec3::magnitude() {
  return sqrt(pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2));
}

void vec3::reverse() {
  this->x *= -1;
  this->y *= -1;
  this->z *= -1;
}

double vec3::distance(vec3 &v) {
  return sqrt(pow(v.x - this->x, 2) + pow(v.y - this->y, 2) + pow(v.z - this->z, 2));
}

vec3::~vec3() {

}
