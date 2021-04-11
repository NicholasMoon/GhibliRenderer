// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// plane.cpp

#include "plane.h"


plane::plane(double A, double B, double C, double D, double color[3]) {
  this->A = A;
  this->B = B;
  this->C = C;
  this->D = D;
  this->c[0] = color[0];
  this->c[1] = color[1];
  this->c[2] = color[2];
  if (A != 0) {
  	vec3 point(-D/A,0,0);
  	p = point;
  }
  else if (B != 0) {
	vec3 point(0,-D/B,0);
  	p = point;
  }
  else {
	vec3 point(0,0,-D/C);
  	p = point;
  }
  
  vec3 normal(A,B,C);
  normal.normalize();
  n = normal;
}

plane::plane(double A, double B, double C, double D, double color[3], material *mat, int objectID) {
  this->objectID = objectID;
  this->A = A;
  this->B = B;
  this->C = C;
  this->D = D;
  this->c[0] = color[0];
  this->c[1] = color[1];
  this->c[2] = color[2];
  this->mat = mat;
  if (A != 0) {
  	vec3 point(-D/A,0,0);
  	p = point;
  }
  else if (B != 0) {
	vec3 point(0,-D/B,0);
  	p = point;
  }
  else {
	vec3 point(0,0,-D/C);
  	p = point;
  }
  
  vec3 normal(A,B,C);
  normal.normalize();
  n = normal;
}

bool plane::hit(ray *incoming_ray, std::vector<object*> &objects, std::vector<light*> &lights, double color[4], double &distance) {
	double cosineNE = this->n.dot(incoming_ray->direction);
	if (cosineNE == 0) {
		return false;
	}
	vec3 rayToPoint(this->p.x - incoming_ray->origin.x, this->p.y - incoming_ray->origin.y, this->p.z - incoming_ray->origin.z);
	vec3 normal = this->n;
        double unitDir = incoming_ray->direction.dot(normal);
	double pointDotDir = normal.dot(rayToPoint);
	if (unitDir < 0) {
		normal.reverse();
		unitDir = incoming_ray->direction.dot(normal);
		pointDotDir = normal.dot(rayToPoint);
	}
	if (pointDotDir < 0) {
		return false;
	}
	double stepsToTake = pointDotDir / unitDir;
	double hitX = incoming_ray->origin.x + stepsToTake * incoming_ray->direction.x;
	double hitY = incoming_ray->origin.y + stepsToTake * incoming_ray->direction.y;
	double hitZ = incoming_ray->origin.z + stepsToTake * incoming_ray->direction.z;
	vec3 hitPoint(hitX, hitY, hitZ);
	if (incoming_ray->direction.x * hitX < 0 || incoming_ray->direction.y * hitY < 0 || incoming_ray->direction.z * hitZ < 0) {
		return false;
	}
	double distanceFromRay = hitPoint.magnitude();
	if (distanceFromRay >= distance) {
		return false;
	}
	distance = distanceFromRay;
	return true;
}

bool plane::shadowHit(ray *incoming_ray, light* target_light, double &distance) {
	double cosineNE = this->n.dot(incoming_ray->direction);
	if (cosineNE == 0) {
		return false;
	}
	vec3 rayToPoint(this->p.x - incoming_ray->origin.x, this->p.y - incoming_ray->origin.y, this->p.z - incoming_ray->origin.z);
	vec3 normal = this->n;
        double unitDir = incoming_ray->direction.dot(normal);
	double pointDotDir = normal.dot(rayToPoint);
	if (unitDir < 0) {
		normal.reverse();
		unitDir = incoming_ray->direction.dot(normal);
		pointDotDir = normal.dot(rayToPoint);
	}
	if (pointDotDir < 0) {
		return false;
	}
	double stepsToTake = pointDotDir / unitDir;
	double hitX = incoming_ray->origin.x + stepsToTake * incoming_ray->direction.x;
	double hitY = incoming_ray->origin.y + stepsToTake * incoming_ray->direction.y;
	double hitZ = incoming_ray->origin.z + stepsToTake * incoming_ray->direction.z;
	vec3 hitPoint(hitX, hitY, hitZ);
	if (incoming_ray->direction.x * hitX < 0 || incoming_ray->direction.y * hitY < 0 || incoming_ray->direction.z * hitZ < 0) {
		return false;
	}
	double distanceFromRay = hitPoint.magnitude();
	if (distanceFromRay >= distance) {
		return false;
	}
	if (distanceFromRay < 0.0000001) {
		return false;
	}
	distance = distanceFromRay;
	return true;
}

bool plane::connectingHit(ray *incoming_ray, double &distance) {
	return true;
}

vec3 plane::getColor() {
	
	return vec3(this->c[0], this->c[1], this->c[2]);
}

vec3 plane::getNormal(double x, double y, double z, int flat) {
	return this->n;
}

vec3 plane::getTextureCoordinates(vec3 &hit_point) {
	// TODO: implement
	return vec3(0,0,0);
}

bool plane::in_bounding_box(AABB *bounding_box) {
	// Plane-AABB collision - gdbooks.gitbooks.io/3dcollisions/content/Chapter2/static_aabb_plane.html

	vec3 boxCenter(0, 0, 0);
	boxCenter.x = (bounding_box->min_coordinates.x + bounding_box->max_coordinates.x) / 2;
	boxCenter.y = (bounding_box->min_coordinates.y + bounding_box->max_coordinates.y) / 2;
	boxCenter.z = (bounding_box->min_coordinates.z + bounding_box->max_coordinates.z) / 2;

	vec3 positiveExtent = bounding_box->max_coordinates.subtract(boxCenter);

	double halfLengthLine = std::abs(positiveExtent.x * this->n.x) + std::abs(positiveExtent.y * this->n.y) + std::abs(positiveExtent.z * this->n.z);
	double distanceToPlane = std::abs(boxCenter.dot(this->n) + this->D);

	if (distanceToPlane > halfLengthLine) return false;

	return true;
}

void plane::updateWorldBoundaries(vec3 &min_coordinates, vec3 &max_coordinates) {

}

plane::~plane() {

}
