// #ifndef STROKE_H
// #define STROKE_H

// #include "utilities.h"
// #include <vector>

// class Stroke {
// public:
//     double x0, y0;
//     vec3 color, direction, normal;
//     double depth;
//     std::vector<std::tuple<double, double>> points;
//     double curvature_filter; // higher with larger gradient
//     int min_length, max_length; // change depending on distance from camera and position (use smaller strokes for edges) and image size

//     Stroke(double x_0, double y_0, vec3 c, vec3 dir, vec3 norm, double dep) {
//         x0 = x_0;
//         y0 = y_0;
//         color = c;
//         direction = dir;
//         normal = norm;
//         depth = dep;
//     }

//     void create(vec3 gradient, double min_x, double min_y, double max_x, double max_y) {
//         std::tuple<double, double> s0 = std::make_tuple (x0, y0);
//         points.push_back(s0);

//         double x = x0;
//         double y = y0;
//         double lastDx = 0;
//         double lastDy = 0;

//         for (int i=1; i < max_length; i++) {
//             // Get unit vector of gradient
//             double gx = gradient.x;
//             double gy = -gradient.y;
            
//             if (gx == 0 && gy == 0) continue;
//             double dx = -gy;
//             double dy = gx;
//             // normal directions: (-gy, gx) or (gy, -gx)

//             if ((lastDx * dx + lastDy * dy) < 0) {
//                 dx = -dx;
//                 dy = -dy;
//             }
//             // dx = curvature_filter * dx + (1 - curvature_filter) * lastDx;
//             // dy = curvature_filter * dy + (1 - curvature_filter) * lastDy;
//             // dx = dx / sqrt(dx * dx + dy * dy);
//             // dy = dy / sqrt(dx * dx + dy * dy);

//             // float yDistort = Vec3(dx, dy, 0).mag() * style.curvature_filter * (x - 0.5) * (x - 0.5);

//             // Filter stroke direction
//             x = x + dx;
//             y = y + dy;
//             lastDx = dx;
//             lastDy = dy;

//             if (x < min_x || x >= max_x) continue;
//             if (y < min_y || y >= max_y) continue;

//             std::tuple<double, double> s = std::make_tuple (x, y);
//             points.push_back(s);
//         }
//     }

//     void set_length(int length) { // based on pre-initalized small, medium, large strokes and reflected stroke length (length 1)
//         max_length = length;
//     }

//     void set_curvature(double curvature) { // value from 0 to 1
//         curvature_filter = curvature;
//     } 

// };

// #endif

#ifndef STROKE_H
#define STROKE_H

#include "utilities.h"
#include <vector>

class stroke {
public:
    stroke(double x_0, double y_0, vec3 col, vec3 dir, vec3 norm, double dep);
    void create(vec3 gradient, double min_x, double min_y, double max_x, double max_y);
    void set_length(int length);
    void set_curvature(double curvature);
    ~stroke();
    double x0, y0;
    vec3 color, direction, normal;
    double depth;
    std::vector<std::tuple<double, double>> points;
    double curvature_filter; // higher with larger gradient
    int min_length, max_length; // change depending on distance from camera and position (use smaller strokes for edges) and image size
private:

};

#endif
