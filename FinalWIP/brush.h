#ifndef BRUSH_H
#define BRUSH_H

#include <algorithm>
#include "utilities.h"
#include "style.h"
#include "stroke.h"

class brush {
public:
    brush(int radius, int falloff);
    void create_mask();
    double get_mask_value(std::tuple<int, int> offset);
    double constant(std::tuple<int, int> offset);
    double linear(std::tuple<int, int> offset);
    // double inverse_square(std::tuple<int, int> offset);
    // double smoothstep(std::tuple<int, int> offset);
    void paint(stroke *stroke, ImageBuffers *imageBuffers, CImg<float> *img);
    ~brush();
    int radius, mask_size, falloff;
    double **mask;
private:

};

#endif
