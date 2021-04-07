// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// ray.h

#ifndef RAY_H
#define RAY_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <limits>
#include <random>


class object;
class OctreeNode;
class Scene;
class HitRecord;
#include "vec3.h"
#include "object.h"
#include "light.h"
#include "octreenode.h"
#include "scene.h"
#include "image_buffers.h"
#include "hit_record.h"


class ray {
 public:
  ray(double ox, double oy, double oz, double dx, double dy, double dz) : origin(ox,oy,oz), direction(dx,dy,dz) { 
	inverse_direction.x = 1.0 / dx;
	inverse_direction.y = 1.0 / dy;
	inverse_direction.z = 1.0 / dz;
	dir_sign[0] = (inverse_direction.x < 0);
	dir_sign[1] = (inverse_direction.y < 0);
	dir_sign[2] = (inverse_direction.z < 0);
  };

  // bool cast(Scene *theScene, ImageBuffers *imageBuffers, double color[4], int bounces, int lastObject, int primary_ray, int x, int y, int &primary_objID, vec3 &hit_normal, int &object_type, std::vector<int> &hit_list, int &shadowed, int indirect_bounces);
  bool cast(Scene *theScene, ImageBuffers *imageBuffers, HitRecord *hitRecord, double color[4], int bounces, int indirect_bounces);
  double detect_edge(Scene *theScene, double color[4], int bounces, int lastObject, int &stencil_objID, int &object_type, std::vector<int> &hit_list);
  double castLight(std::vector<object*> &objects, light *targetLight, double distance);
  ~ray();
  vec3 origin;
  vec3 direction;
  vec3 inverse_direction;
  int dir_sign[3];
 private:
  
};

#endif
