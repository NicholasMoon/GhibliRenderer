#include "brush.h"

// at init time, initialize set of 3 brushes - small, medium, large (do based on image size)
// change depending on area and distance from the camera
// make it part of style.h

brush::brush(int radius) {
    this->radius = radius;
    this->mask_size = (radius * 2) + 1;
    // this->falloff = falloff;
    // this->img = i;

    this->mask = new double*[this->mask_size];
    for (int i=0; i < this->mask_size; i++) {
        this->mask[i] = new double[this->mask_size];
    }
}

// ~Brush() {
//     for (int i=0; i < mask_size; i++) {
//         delete [] mask[i];
//     }
//     delete [] mask;
// }

void brush::create_mask() {
    int radius = this->radius;
    for (int i=-radius; i < radius+1; i++) {
        for (int j=-radius; j < radius+1; j++) {
            auto offset = std::make_tuple (i, j);
            double mask_val = this->falloff(offset); 
            this->mask[i+radius][j+radius] = mask_val;
        }
    }
}

// double brush::get_mask_value(std::tuple<int, int> offset) {
//     // double mask_val = 0;
//     // switch(falloff) {
//     //     case 1:
//     //         mask_val = constant(offset);
//     //         break;
//     //     case 2:
//     //         mask_val = linear(offset);
//     //         break;
//     //     case 3:
//     //         mask_val = inverse_square(offset);
//     //         break;
//     //     case 4:
//     //         mask_val = smoothstep(offset);
//     //         break;
//     //     default:
//     //         mask_val = constant(offset);
//     // }
//     double mask_val = brush::falloff(offset);
//     return mask_val;
// }

// Falloff functions - return alpha of color
// double brush::constant(std::tuple<int, int> offset) {
//     int x = (0 - std::get<0>(offset));
//     int x2 = x * x;
//     int y = (0 - std::get<1>(offset));
//     int y2 = y * y;
//     if (x2 + y2 <= this->radius * this->radius) return 1; 
//     return 0;
// }

// double brush::linear(std::tuple<int, int> offset) {
//     double dist = distance(0, std::get<0>(offset), 0, std::get<1>(offset));
//     if (dist <= this->radius) {
//         double y = 1 - (1 / (double)this->radius) * dist;
//         // std::cout << "dist: " << dist << " y: " << y << std::endl;
//         return y;
//     }
//     return 0;
// }

// double brush::inverse_square(std::tuple<int, int> offset) {
//     double dist = distance(0, std::get<0>(offset), 0, std::get<1>(offset));
//     if (dist <= this->radius) {
//         if (dist == 0) return 1;
//         double y = 1 / (dist * dist);
//         return y; 
//     }
//     return 0;
// }

// double brush::smoothstep(std::tuple<int, int> offset) {
//     int x = (0 - std::get<0>(offset));
//     int x2 = x * x;
//     int y = (0 - std::get<1>(offset));
//     int y2 = y * y;
//     if (x2 + y2 <= this->radius * this->radius) return 1; 
//     return 0;
// }

void brush::paint(stroke *stroke, float *paintMap, int *objectTypeMap, int *objectBoundaryMap, int *objectIDMap, int primaryObjectBoundary, int currObjectID, CImg<float> *img) { 
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
                // std::cout << "i: " << i << " j: " << j << std::endl;
                if (curr_x < 0 || curr_x >= img->width()) continue;
                if (curr_y < 0 || curr_y >= img->height()) continue;
                if (objectTypeMap[curr_y * img->width() + curr_x]) continue; 
                if (objectBoundaryMap[curr_y * img->width() + curr_x] != primaryObjectBoundary) continue; 
                if (objectIDMap[curr_y * img->width() + curr_x] != currObjectID) continue;
                double paint_num = this->mask[j+radius][k+radius];
                if (paint_num) {
                    vec3 mix_color, final_color;
                    double mix_alpha, final_alpha;
                    final_color = color;
                    final_alpha = 255;
                    // if (paintMap[curr_y * img->width() + curr_x] == 0) {
                    //     mix_color = vec3(1, 1, 1);
                    //     mix_alpha = 1 - paint_num;
                    //     alpha_composite(stroke.color, mix_color, paint_num, mix_alpha, final_color, final_alpha);
                    //     paintMap[curr_y * img->width() + curr_x] = 1;
                    // } 
                    // else {
                    //     double r = ((*img)(curr_x, curr_y, 0, 0) / 255);
                    //     double g = ((*img)(curr_x, curr_y, 0, 1) / 255);
                    //     double b = ((*img)(curr_x, curr_y, 0, 2) / 255);
                    //     mix_color = vec3(r, g, b);
                    //     mix_alpha = ((*img)(curr_x, curr_y, 0, 3) / 255);
                    //     if (paint_num > mix_alpha) {
                    //         mix_alpha = 1 - paint_num;
                    //     } else {
                    //         paint_num = 1 - mix_alpha;
                    //     }
                    //     alpha_composite(stroke.color, mix_color, paint_num, mix_alpha, final_color, final_alpha);
                    // }
                    write_color(*img, curr_x, curr_y, final_color, final_alpha); 
                }
            }
        }
    }
}

brush::~brush() {

}

double constbrush::falloff(std::tuple<int, int> offset) {
    int x = (0 - std::get<0>(offset));
    int x2 = x * x;
    int y = (0 - std::get<1>(offset));
    int y2 = y * y;
    if (x2 + y2 <= this->radius * this->radius) return 1; 
    return 0;
}