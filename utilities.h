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

std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0.0, 1.0);

float random_number() {
    return distribution(generator);
}

float clamp(float color, float min, float max) {
    float c = color;
    if (c >= max) {
        c = max;
    } else if (c <= min) {
        c = min;
    }

    return c;
}

void write_color(CImg<float> &img, int x, int y, vec3 color, int alpha) {
    img(x, y, 0, 0) = clamp(color.x, 0, 1) * 255;
    img(x, y, 0, 1) = clamp(color.y, 0, 1) * 255;
    img(x, y, 0, 2) = clamp(color.z, 0, 1) * 255;
    img(x, y, 0, 3) = alpha;
}

#endif
