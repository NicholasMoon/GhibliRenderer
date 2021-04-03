#include "stroke.h"

stroke::stroke(double x_0, double y_0, vec3 col, vec3 dir, vec3 norm, double dep, int primary_obj_bound, int curr_obj_id) {
    this->x0 = x_0;
    this->y0 = y_0;
    this->color = col;
    this->direction = dir;
    this->normal = norm;
    this->depth = dep;
    this->primaryObjectBoundary = primary_obj_bound;
    this->currObjectID = curr_obj_id; 
}

void stroke::create(vec3 gradient, bool inside, double min_x, double min_y, double max_x, double max_y) {
    std::tuple<double, double> s0 = std::make_tuple (x0, y0);
    this->points.push_back(s0);

    double x = x0;
    double y = y0;
    double lastDx = 0;
    double lastDy = 0;

    for (int i=1; i < max_length; i++) {
        // // std::cout << "x: " << x << "y: " << y << std::endl;
        // // Get unit vector of gradient
        // double gx = normal.x;
        // double gy = normal.y;
        // // std::cout << "gx: " << gx << " gy: " << gy << std::endl;
        // if (gx == 0 && gy == 0) continue;
        // double dx = -gy;
        // double dy = -gx;

        // Get unit vector of gradient
        double gx = gradient.x;
        double gy = -gradient.y; // rotate normal to be oriented correctly in x-y plane (positive y should go up not down)
        if (inside) gy = gradient.y;

        if (gx == 0 && gy == 0) continue;
        double dx = -gy;
        double dy = gx;
        // double dx = gx;
        // double dy = gy; 
        // normal directions: (-gy, gx) or (gy, -gx)

        if ((lastDx * dx + lastDy * dy) < 0) {
            dx = -dx;
            dy = -dy;
        }
        // dx = curvature_filter * dx + (1 - curvature_filter) * lastDx;
        // dy = curvature_filter * dy + (1 - curvature_filter) * lastDy;
        // dx = dx / sqrt(dx * dx + dy * dy);
        // dy = dy / sqrt(dx * dx + dy * dy);

        // float yDistort = Vec3(dx, dy, 0).mag() * style.curvature_filter * (x - 0.5) * (x - 0.5);

        // Filter stroke direction
        x = x + dx;
        y = y + dy;
        lastDx = dx;
        lastDy = dy;

        if (x < min_x || x >= max_x) continue;
        if (y < min_y || y >= max_y) continue;

        std::tuple<double, double> s = std::make_tuple (x, y);
        this->points.push_back(s);
    }
}

void stroke::set_length(int length) { // based on pre-initalized small, medium, large strokes and reflected stroke length (length 1)
    this->max_length = length;
}

void stroke::set_curvature(double curvature) { // value from 0 to 1
    this->curvature_filter = curvature;
} 

stroke::~stroke() {
    
}