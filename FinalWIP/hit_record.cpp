// octreenode.cpp

#include "hit_record.h"

HitRecord::HitRecord(int x, int y) {
    this->hit_point = vec3(0,0,0);
    this->primary_hit_normal = vec3(0,0,0);
    this->hit_normal = vec3(0,0,0);
    this->distance = 0.0;
    this->primary_ray = -1;
    this->x = x;
    this->y = y;
    this->lastObject = -1;
    this->primary_objID = -1;
    this->object_type = 0;
    this->shadowed = 0;
}

void HitRecord::setAttributes(double distance, vec3 &hit_point, vec3 &hit_normal) {
    this->distance = distance;
    this->hit_point = hit_point;
    this->primary_hit_normal = hit_normal;
    // this->hit_normal = hit_normal;
}

HitRecord::~HitRecord() {

}
