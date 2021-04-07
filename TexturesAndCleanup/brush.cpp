#include "brush.h"

// at init time, initialize set of 3 brushes - small, medium, large (do based on image size)
// change depending on area and distance from the camera
// make it part of style.h

brush::brush(int radius, int falloff) {
    this->radius = radius;
    this->mask_size = (radius * 2) + 1;
    this->falloff = falloff;

    this->mask = new double*[this->mask_size];
    for (int i=0; i < this->mask_size; i++) {
        this->mask[i] = new double[this->mask_size];
    }
}

void brush::create_mask() {
    int radius = this->radius;
    for (int i=-radius; i < radius+1; i++) {
        for (int j=-radius; j < radius+1; j++) {
            auto offset = std::make_tuple (i, j);
            double mask_val = get_mask_value(offset); 
            this->mask[i+radius][j+radius] = mask_val;
        }
    }
}

double brush::get_mask_value(std::tuple<int, int> offset) {
    double mask_val = 0;
    switch(this->falloff) {
        case 1:
            mask_val = constant(offset);
            break;
        case 2:
            mask_val = linear(offset);
            break;
        default:
            mask_val = constant(offset);
    }
    return mask_val;
}

// Falloff functions - return alpha of color
double brush::constant(std::tuple<int, int> offset) {
    int x = (0 - std::get<0>(offset));
    int x2 = x * x;
    int y = (0 - std::get<1>(offset));
    int y2 = y * y;
    if (x2 + y2 <= this->radius * this->radius) return 1; 
    return 0;
}

double brush::linear(std::tuple<int, int> offset) {
    double dist = distance(0, std::get<0>(offset), 0, std::get<1>(offset));
    if (dist <= this->radius) {
        double y = 1 - (1 / (double)this->radius) * dist;
        return y;
    }
    return 0;
}

// double brush::inverse_square(std::tuple<int, int> offset) {
//     double dist = distance(0, std::get<0>(offset), 0, std::get<1>(offset));
//     if (dist <= this->radius) {
//         if (dist == 0) return 1;
//         double y = 1 / (dist * dist);
//         return y; 
//     }
//     return 0;
// }

// void brush::paint(stroke *stroke, float *paintMap, int *objectTypeMap, int *objectBoundaryMap, int *objectIDMap, int primaryObjectBoundary, int currObjectID, CImg<float> *img) {
void brush::paint(stroke *stroke, ImageBuffers *imageBuffers, CImg<float> *img) {  
    for (int i=0; i < stroke->points.size(); i++) { // iterate through all points in stroke
        auto curr_stroke_point = stroke->points.at(i);
        vec3 color = stroke->color; 

        int x = std::get<0>(curr_stroke_point);
        int y = std::get<1>(curr_stroke_point);
        // Paint surrounding pixels
        for (int j=-this->radius; j < this->radius+1; j++) {
            for (int k=-this->radius; k < this->radius+1; k++) {
                int curr_x = x + j;
                int curr_y = y + k;
                if (curr_x < 0 || curr_x >= img->width()) continue;
                if (curr_y < 0 || curr_y >= img->height()) continue;
                if (imageBuffers->objectTypeMap[curr_y * img->width() + curr_x]) continue; 
                if (imageBuffers->objectBoundaryMap[curr_y * img->width() + curr_x] != stroke->primaryObjectBoundary) continue; 
                if (imageBuffers->objectIDMap[curr_y * img->width() + curr_x] != stroke->currObjectID) continue;
                if (imageBuffers->environmentMap[curr_y * img->width() + curr_x]) continue;
                double paint_num = this->mask[j+radius][k+radius];
                if (paint_num) {
                    vec3 mix_color, final_color;
                    double mix_alpha, final_alpha;
                    final_color = color;
                    final_alpha = 255;
                    imageBuffers->paintMap[curr_y * img->width() + curr_x] = 1;
                    if (imageBuffers->paintMap[curr_y * img->width() + curr_x] == 0) { // mix canvas color with paint color
                        double r = ((*img)(curr_x, curr_y, 0, 0) / 255);
                        double g = ((*img)(curr_x, curr_y, 0, 1) / 255);
                        double b = ((*img)(curr_x, curr_y, 0, 2) / 255);
                        mix_color = vec3(r, g, b);
                        mix_alpha = 1 - paint_num;
                        alpha_composite(stroke->color, mix_color, paint_num, mix_alpha, final_color, final_alpha);
                        imageBuffers->paintMap[curr_y * img->width() + curr_x] = paint_num;
                    } 
                    else {
                        double r = ((*img)(curr_x, curr_y, 0, 0) / 255);
                        double g = ((*img)(curr_x, curr_y, 0, 1) / 255);
                        double b = ((*img)(curr_x, curr_y, 0, 2) / 255);
                        mix_color = vec3(r, g, b);
                        mix_alpha = imageBuffers->paintMap[curr_y * img->width() + curr_x];
                        if (paint_num > mix_alpha) {
                            mix_alpha = 1 - paint_num;
                            imageBuffers->paintMap[curr_y * img->width() + curr_x] = paint_num;
                        } 
                        alpha_composite(stroke->color, mix_color, paint_num, mix_alpha, final_color, final_alpha);
                    }
                    write_color(*img, curr_x, curr_y, final_color, final_alpha); 
                }
            }
        }
    }
}

brush::~brush() {

}