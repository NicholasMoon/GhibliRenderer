#ifndef BRUSH_H
#define BRUSH_H

#include <algorithm>
#include "utilities.h"
#include "style.h"
#include "stroke.h"

class brush {
public:
    brush(int radius);
    // brush(int radius, int falloff);
    // double distance(double x, double x0, double y, double y0);
    virtual double falloff(std::tuple<int, int> offset) = 0;
    void create_mask();
    // double get_mask_value(std::tuple<int, int> offset);
    // double constant(std::tuple<int, int> offset);
    // double linear(std::tuple<int, int> offset);
    // double inverse_square(std::tuple<int, int> offset);
    // double smoothstep(std::tuple<int, int> offset);
    void paint(stroke *stroke, float *paintBuffer, int *objectTypeMap, int *objectBoundaryMap, int *objectIDMap, int currObjectBoundary, int currObjectID, CImg<float> *img);
    ~brush();
    int radius, mask_size;
    double **mask;
private:

};

class constbrush : public brush {
public:
    using brush::brush;
    double falloff(std::tuple<int, int> offset);

private:

};

#endif
