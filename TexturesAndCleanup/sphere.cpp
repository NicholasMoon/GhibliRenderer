// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// sphere.cpp

#include "sphere.h"


sphere::sphere(double x, double y, double z, double r, double c[3]) {
  this->x = x;
  this->y = y;
  this->z = z;
  this->r = r;
  this->c[0] = c[0];
  this->c[1] = c[1];
  this->c[2] = c[2];
}

sphere::sphere(double x, double y, double z, double r, double c[3], material *mat, int objectID) {
  this->objectID = objectID;
  this->x = x;
  this->y = y;
  this->z = z;
  this->r = r;
  this->c[0] = c[0];
  this->c[1] = c[1];
  this->c[2] = c[2];
  this->mat = mat;
}

sphere::sphere(double x, double y, double z, double r, double c[3], material *mat, int objectID, int object_type) {
  this->objectID = objectID;
  this->object_type = object_type;
  this->x = x;
  this->y = y;
  this->z = z;
  this->r = r;
  this->c[0] = c[0];
  this->c[1] = c[1];
  this->c[2] = c[2];
  this->mat = mat;
}

sphere::sphere(double x, double y, double z, double r, double c[3], double e[3], material *mat, int objectID, int object_type) { // new addition
  this->objectID = objectID;
  this->object_type = object_type;
  this->x = x;
  this->y = y;
  this->z = z;
  this->r = r;
  this->c[0] = c[0];
  this->c[1] = c[1];
  this->c[2] = c[2];
  this->e[0] = e[0];
  this->e[1] = e[1];
  this->e[2] = e[2];
  this->mat = mat;
}

bool sphere::hit(ray *incoming_ray, std::vector<object*> &objects, std::vector<light*> &lights, double color[4], double &distance) {
	vec3 rayToCenter(this->x - incoming_ray->origin.x, this->y - incoming_ray->origin.y, this->z - incoming_ray->origin.z);
	double centerDotDir = incoming_ray->direction.dot(rayToCenter);
	double midX = incoming_ray->origin.x + centerDotDir * incoming_ray->direction.x;
	double midY = incoming_ray->origin.y + centerDotDir * incoming_ray->direction.y;
	double midZ = incoming_ray->origin.z + centerDotDir * incoming_ray->direction.z;
	vec3 midpointToCenter(midX - this->x, midY - this->y, midZ - this->z);
	double PtoCDistance = midpointToCenter.dot(midpointToCenter);
	double discriminant = (this->r * this->r) - PtoCDistance;
	if (discriminant < 0) {
		return false;
	}
	double distanceFromRay = 0;
	if (rayToCenter.magnitude() >= this->r) {
		distanceFromRay = centerDotDir - sqrt(discriminant);
	}
	else {
		distanceFromRay = centerDotDir + sqrt(discriminant);
	}
	if (distanceFromRay >= distance) {
		return false;
	}
	if (distanceFromRay < 0.0000001) {
		return false;
	}
	distance = distanceFromRay;
	return true;
}

bool sphere::shadowHit(ray *incoming_ray, light* target_light, double &distance) {
	vec3 rayToCenter(this->x - incoming_ray->origin.x, this->y - incoming_ray->origin.y, this->z - incoming_ray->origin.z);
	double centerDotDir = incoming_ray->direction.dot(rayToCenter);
	double midX = incoming_ray->origin.x + centerDotDir * incoming_ray->direction.x;
	double midY = incoming_ray->origin.y + centerDotDir * incoming_ray->direction.y;
	double midZ = incoming_ray->origin.z + centerDotDir * incoming_ray->direction.z;
	vec3 midpointToCenter(midX - this->x, midY - this->y, midZ - this->z);
	double PtoCDistance = midpointToCenter.dot(midpointToCenter);
	double discriminant = (this->r * this->r) - PtoCDistance;
	if (discriminant < 0) {
		return false;
	}
	double distanceFromRay = centerDotDir - sqrt(discriminant);
	if (distanceFromRay < 0.0000001) {
		return false;
	}
	else if (distanceFromRay >= distance) {
		return false;
	}
	distance = distanceFromRay;
	return true;
}

vec3 sphere::getColor() {
	return vec3(c[0], c[1], c[2]);
}

vec3 sphere::getEmission() { 
	return vec3(e[0], e[1], e[2]);
}

vec3 sphere::getNormal(double x, double y, double z, int flat) {
	return vec3(x - this->x, y - this->y, z - this->z);
}

vec3 sphere::getTextureCoordinates(vec3 &hit_point) {
	// TODO: implement
	return vec3(0,0,0);
}

bool sphere::in_bounding_box(AABB *bounding_box) {
	// Sphere-AABB collision -> gdbooks.gitbooks.io/3dcollisions/content/Chapter2/static_sphere_aabb.html

	vec3 sphereCenter(x, y, z);
	vec3 closestPointOnBB = bounding_box->getClosestPoint(sphereCenter);
	
	vec3 differenceCenterPoint = sphereCenter.subtract(closestPointOnBB);
	double distance = differenceCenterPoint.magnitude();

	if (distance <= this->r) return true;

	return false;
}

void sphere::updateWorldBoundaries(vec3 &min_coordinates, vec3 &max_coordinates) {
	
}

sphere::~sphere() {

}
