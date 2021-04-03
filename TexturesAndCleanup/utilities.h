#ifndef UTILITIES_H
#define UTILITIES_H

#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>
#include <memory>
#include <ctime>

#include "vec3.h"
#include "object.h"
#include "CImg.h"
#include "scene.h"

using namespace cimg_library;

float random_number();
double distance(double x, double x0, double y, double y0);
double sample_parabola(double a, double x);
void alpha_composite(vec3 &c1, vec3 &c2, double a1, double a2, vec3 &final_color, double &final_alpha);
float clamp(float color, float min, float max);
void write_color(CImg<float> &img, int x, int y, vec3 color, int alpha);
double printProgress(int pixnum, int totalpixels, double milestone, time_t start);
void getWorldBoundaries(vec3 &min_coordinates, vec3 &max_coordinates, std::vector<object*> &objects);
double Max_Double(double a, double b, double c);
double Min_Double(double a, double b, double c);
double clip(double color);
double maxItoD(int x, int y);
vec3 cameraDirection(double xi, double yi, Scene *theScene);
void getDrawingGradient(vec3 hit_normal, double *old_color, double *light_x_color, double *light_y_color, vec3 &final_gradient, int pixnum);


#endif

