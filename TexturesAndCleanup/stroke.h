#ifndef STROKE_H
#define STROKE_H

#include "utilities.h"
#include <vector>

class stroke {
public:
    // stroke(double x_0, double y_0, vec3 col, vec3 dir, vec3 norm, double dep, int primary_obj_bound, int curr_obj_id);
    stroke(double x_0, double y_0, vec3 col, vec3 dir, HitRecord *hit_record);
    void create(vec3 gradient, bool inside, double min_x, double min_y, double max_x, double max_y);
    void set_length(int length);
    void set_curvature(double curvature);
    ~stroke();
    double x0, y0;
    vec3 color, direction, normal;
    double depth;
    int primaryObjectBoundary, currObjectID;
    std::vector<std::tuple<double, double>> points;
    double curvature_filter; // higher with larger gradient
    int min_length, max_length; // change depending on distance from camera and position (use smaller strokes for edges) and image size
private:

};

#endif
