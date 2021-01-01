#ifndef STYLE_H
#define STYLE_H

#include <vector>

class Style {
public:
    int brush_size; // change to array when you want multiple sizes
    int curvature_filter;
    int min_stroke_length, max_stroke_length;
};

class Pointillist : public Style {
public:
    Pointillist() {
        brush_size = 5;
        curvature_filter = 1;
        min_stroke_length = 1;
        max_stroke_length = 1;
    }
};

class Expressionist : public Style {
public:
    Expressionist() {
        brush_size = 2;
        curvature_filter = 1;
        min_stroke_length = 1;
        max_stroke_length = 4;
    }
};

#endif
