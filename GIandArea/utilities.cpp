#include "utilities.h"

std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0.0, 1.0);

float random_number() {
    return distribution(generator);
}

double distance(double x, double x0, double y, double y0) {
    double x2 = (x - x0) * (x - x0);
    double y2 = (y - y0) * (y - y0);
    return sqrt(x2 + y2);
}

double sample_parabola(double a, double x) { 
    return a * (x * x);
}

// void alpha_composite(vec3 &c1, vec3 &c2, double a1, double a2, vec3 &final_color, double &final_alpha) {
//     final_alpha = a1 + a2 * (1-a1);
//     final_color = (c1 * a1) + c2 * (1.0 - a1);
//     // final_color.print();
//     // std::cout << final_alpha << std::endl;
// }

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
