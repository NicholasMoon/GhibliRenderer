// octreenode.cpp

#include "AABB.h"

bool AABB::intersect(ray *r) {
    // ray-AABB intersection https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
    double orig_values[3] = {r->origin.x, r->origin.y, r->origin.z};
    double dir_values[3] = {r->direction.x, r->direction.y, r->direction.z};
    vec3 bounds[2] = { this->min_coordinates, this->max_coordinates};
    

    double tmin, tmax, tymin, tymax, tzmin, tzmax;

    tmin = (bounds[r->dir_sign[0]].x - orig_values[0]) * r->inverse_direction.x; 
    tmax = (bounds[1 - r->dir_sign[0]].x - orig_values[0]) * r->inverse_direction.x; 
    tymin = (bounds[r->dir_sign[1]].y - orig_values[1]) * r->inverse_direction.y; 
    tymax = (bounds[1 - r->dir_sign[1]].y - orig_values[1]) * r->inverse_direction.y; 
 
    if ((tmin > tymax + 0.0001) || (tymin > tmax + 0.0001)) 
        return false; 
 
    if (tymin > tmin) 
        tmin = tymin; 
 
    if (tymax < tmax) 
        tmax = tymax; 
 
    tzmin = (bounds[r->dir_sign[2]].z - orig_values[2]) * r->inverse_direction.z; 
    tzmax = (bounds[1 - r->dir_sign[2]].z - orig_values[2]) * r->inverse_direction.z; 
 
    if ((tmin > tzmax + 0.0001) || (tzmin > tmax + 0.0001)) 
        return false; 
 
    //if (tzmin > tmin) tmin = tzmin;
    //if (tzmax < tmax) tmax = tzmax; 
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
