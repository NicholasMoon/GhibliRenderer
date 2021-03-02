// octreenode.cpp

#include "AABB.h"

bool AABB::intersect(ray *r) {
    double orig_values[3] = {r->origin.x, r->origin.y, r->origin.z};
    double dir_values[3] = {r->direction.x, r->direction.y, r->direction.z};
    double min_values[3] = {this->min_coordinates.x, this->min_coordinates.y, this->min_coordinates.z};
    double max_values[3] = {this->max_coordinates.x, this->max_coordinates.y, this->max_coordinates.z};
   /* std::cout << "Data: " << std::endl;
    std::cout << orig_values[0] << " " << orig_values[1] << " " << orig_values[2] << std::endl;
    std::cout << dir_values[0] << " " << dir_values[1] << " " << dir_values[2] << std::endl;
    std::cout << min_values[0] << " " << min_values[1] << " " << min_values[2] << std::endl;
    std::cout << max_values[0] << " " << max_values[1] << " " << max_values[2] << std::endl;*/
    double planes[3];
    double distance[3];

    int i0, i1;
    double h0, h1;

    // Check if ray starts inside box
    bool inside = true;
    if (orig_values[0] < this->min_coordinates.x || orig_values[0] > this->max_coordinates.x) inside = false;
    if (orig_values[1] < this->min_coordinates.y || orig_values[1] > this->max_coordinates.y) inside = false;
    if (orig_values[2] < this->min_coordinates.z || orig_values[2] > this->max_coordinates.z) inside = false;
    if (inside) return true;

    // Find 3 planes to test (based on sign of ray direction)
    for (int i = 0; i < 3; i++) {
        if (dir_values[i] < 0) {
            planes[i] = max_values[i];
        }
        else if (dir_values[i] > 0) {
            planes[i] = min_values[i];
        }
    }

    // Find t values
    for (int j = 0; j < 3; j++) {
        if (dir_values[j] == 0) {
            distance[j] = -1;
        }
        else {
            distance[j] = (planes[j] - orig_values[j]) / dir_values[j];
        }
    }
    // Test for intersection
    for (int k = 0; k < 3; k++) {
        i0 = (k + 1) % 3;
        i1 = (k + 2) % 3;
        //if (distance[k] < 0) return false;
	//std::cout << "dis " << distance[i0] << " " << distance[i1] << std::endl;
        //std::cout << "orig " << orig_values[i0] << " " << orig_values[i1] << std::endl;
        //std::cout << "dir " << dir_values[i0] << " " << dir_values[i1] << std::endl;
        h0 = orig_values[i0] + distance[i0] * dir_values[i0];
        h1 = orig_values[i1] + distance[i1] * dir_values[i1];
	//std::cout << h0 << " " << h1 << std::endl;
        if (h0 < min_values[i0] || h0 > max_values[i0]) return false;
        if (h1 < min_values[i1] || h1 > max_values[i1]) return false;
    }
    return true;
}

vec3 AABB::getClosestPoint(vec3 point) {
    vec3 closest_point(point.x, point.y, point.z);

    if (closest_point.x < this->min_coordinates.x) {
        closest_point.x = this->min_coordinates.x;
    }
    else if (closest_point.x > this->max_coordinates.x) {
        closest_point.x = this->max_coordinates.x;
    }

    if (closest_point.y < this->min_coordinates.y) {
        closest_point.y = this->min_coordinates.y;
    }
    else if (closest_point.y > this->max_coordinates.y) {
        closest_point.y = this->max_coordinates.y;
    }

    if (closest_point.z < this->min_coordinates.z) {
        closest_point.z = this->min_coordinates.z;
    }
    else if (closest_point.z > this->max_coordinates.z) {
        closest_point.z = this->max_coordinates.z;
    }

    return closest_point;
}

void AABB::PrintAABB() {
    std::cout << "Min: (" << this->min_coordinates.x << ", " << this->min_coordinates.y << ", " << this->min_coordinates.z << "); Max: (" << this->max_coordinates.x << ", " << this->max_coordinates.y << ", " << this->max_coordinates.z << ")" << std::endl;
}

AABB::~AABB() {

}
