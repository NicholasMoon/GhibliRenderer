#include "painter.h"


painter::painter(int width, int height, std::vector<stroke *> bottom, std::vector<stroke *> middle, std::vector<stroke *> top) {
    this->canvasWidth = width;
    this->canvasHeight = height;
    this->bottomLayer = bottom;
    this->middleLayer = middle;
    this->topLayer = top;

    // initialize based on resolution
    brush *brush_small = new brush(1, 2);
    brush_small->create_mask();

    brush *brush_medium = new brush(3, 2);
    brush_medium->create_mask();

    brush *brush_large = new brush(6, 2);
    brush_large->create_mask();

    this->brushes[0] = brush_large;
    this->brushes[1] = brush_medium;
    this->brushes[2] = brush_small;
}



void painter::paint_layer(std::vector<stroke *> layer, brush *brush, ImageBuffers *imageBuffers, CImg<float> *img) {
    std::shuffle(layer.begin(), layer.end(), std::random_device()); 
    
    for (int s=0; s < layer.size(); s++) {
        brush->paint(layer.at(s), imageBuffers, img);
    }
    
}

void painter::paint(ImageBuffers *imageBuffers, CImg<float> *img) {
    std::vector<stroke *> current_layer;
    brush *current_brush;
    for (int i=0; i < 3; i++) {
        if (i == 0) {
            current_layer = this->bottomLayer; 
            current_brush = this->brushes[0];
        }
        else if (i == 1) {
            current_layer = this->middleLayer; 
            current_brush = this->brushes[1];
        }
        else if (i == 2) {
            current_layer = this->topLayer; 
            current_brush = this->brushes[2];
        }
        paint_layer(current_layer, current_brush, imageBuffers, img);
    }  
}

painter::~painter() {

}