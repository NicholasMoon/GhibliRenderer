// octree.h

#ifndef OCTREENODE_H
#define OCTREENODE_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <limits>
#include <random>


class object;
#include "vec3.h"
#include "object.h"
#include "light.h"

class OctreeNode {
 public:
  OctreeNode();
  OctreeNode(int num_obj_cutoff, int level, int position);
  OctreeNode(AABB *bounding_box, int num_obj_cutoff, int level, int position);
  void addChildNode(int position, std::vector<object*> &objects);
  void buildOctree(std::vector<object*> &objects);
  bool traverseOctree(ray *r, std::vector<object*> *objects);
  void getMinCoordinates(int position, vec3 &min_coordinates);
  void getMaxCoordinates(int position, vec3 &max_coordinates);
  void printOctree();
  ~OctreeNode();
  OctreeNode *children[8];
  char child_valid[8] = {0};
  std::vector<object*> node_objects;
  AABB *bounding_box;
  int num_obj_cutoff;
  int level;
  int position;
 private:
  
};

#endif
