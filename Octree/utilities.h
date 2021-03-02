#ifndef UTILITIES_H
#define UTILITIES_H

#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>
#include <memory>

#include "vec3.h"
#include "ray.h"
#include "CImg.h"

using namespace cimg_library;

float random_number();
double distance(double x, double x0, double y, double y0);
double sample_parabola(double a, double x);
// void alpha_composite(vec3 &c1, vec3 &c2, double a1, double a2, vec3 &final_color, double &final_alpha);
float clamp(float color, float min, float max);
void write_color(CImg<float> &img, int x, int y, vec3 color, int alpha);

#endif

