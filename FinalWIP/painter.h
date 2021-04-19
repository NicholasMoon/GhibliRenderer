#ifndef PAINTER_H
#define PAINTER_H

#include <algorithm>
#include "utilities.h"
#include "style.h"
#include "particle.h"
#include "brush.h"
#include "stroke.h"

class painter {
public:
    painter(int width, int height, std::vector<stroke *> bottom, std::vector<stroke *> middle, std::vector<stroke *> top);
    void paint_layer(std::vector<stroke *> layer, brush *brush, ImageBuffers *imageBuffers, CImg<float> *img);
    void paint(ImageBuffers *imageBuffers, CImg<float> *img);
    brush *choose_brush();
    ~painter();
    int canvasWidth, canvasHeight;
    std::vector<stroke *> bottomLayer, middleLayer, topLayer;
    brush *brushes[3];
private:

};

#endif
