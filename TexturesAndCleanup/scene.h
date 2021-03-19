// scene.h

#ifndef SCENE_H
#define SCENE_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <limits>
#include <random>
#include <map>

#include "vec3.h"
#include "vertex.h"
#include "object.h"
#include "plane.h"
#include "sphere.h"
#include "tri.h"
#include "imagetexture.h"
#include "colortexture.h"
#include "texture.h"


class Scene {
 public:
  Scene(char** argv);
  void readMaterialFile(std::string material_file, std::map<std::string, material *> &materials, std::vector<std::string> &texture_files);
  void initTextureMaps(std::vector<std::string> &texture_files, std::map<std::string, image_texture *> &texture_maps);
  void PrintScene();
  ~Scene();


  int width, height;
  int dimension;
  int bounces;
  double environmentColor[4];
  vec3 eye;
  vec3 forward;
  vec3 right;
  vec3 up;
  int stencil_ring_samples;
  double stencil_rings;
  double outline_cutoff;
  int shoot;
  int spp;
  int output_depth;
  int octree_box_limit_x8;
  int max_octree_depth;
  double e_slope;
  double stencil_radius;
  int light_samples;
  int indirect_samples;
  int indirect_bounces;
  int flat;
  int environment;
  std::ifstream inputFile;
  std::string outputFileName;

  std::default_random_engine generator;

  OctreeNode *octree;
  std::vector<light*> lights;
 private:
  
};

#endif
