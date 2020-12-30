// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// tri.cpp

#include "tri.h"


tri::tri(vertex *v1, vertex *v2, vertex *v3, double c[3]) {
	vertex vert1(v1->xyz, v1->normal, v1->color);
	vertex vert2(v2->xyz, v2->normal, v2->color);
	vertex vert3(v3->xyz, v3->normal, v3->color);
	this->v1 = vert1;
	this->v2 = vert2;
	this->v3 = vert3;
	this->c[0] = c[0];
	this->c[1] = c[1];
	this->c[2] = c[2];
	vec3 edge1(this->v2.xyz.x - this->v1.xyz.x, this->v2.xyz.y - this->v1.xyz.y, this->v2.xyz.z - this->v1.xyz.z);
	this->edge1 = edge1;
	vec3 edge2(this->v3.xyz.x - this->v1.xyz.x, this->v3.xyz.y - this->v1.xyz.y, this->v3.xyz.z - this->v1.xyz.z);
	this->edge2 = edge2;
	this->n = edge1.cross(edge2);
	this->n.normalize();

	this->e2 = this->edge2.cross(this->n);
	this->e3 = this->edge1.cross(this->n);

	double dot2 = this->e2.dot(this->edge1);
	double dot3 = this->e3.dot(this->edge2);

	this->e2.x /= dot2;
	this->e2.y /= dot2;
	this->e2.z /= dot2;

	this->e3.x /= dot3;
	this->e3.y /= dot3;
	this->e3.z /= dot3;
}

tri::tri(vertex *v1, vertex *v2, vertex *v3, double c[3], material *mat, int objectID) {
	this->objectID = objectID;
	vertex vert1(v1->xyz, v1->normal, v1->color);
	vertex vert2(v2->xyz, v2->normal, v2->color);
	vertex vert3(v3->xyz, v3->normal, v3->color);
	this->v1 = vert1;
	this->v2 = vert2;
	this->v3 = vert3;
	this->c[0] = c[0];
	this->c[1] = c[1];
	this->c[2] = c[2];
	this->mat = mat;
	vec3 edge1(this->v2.xyz.x - this->v1.xyz.x, this->v2.xyz.y - this->v1.xyz.y, this->v2.xyz.z - this->v1.xyz.z);
	this->edge1 = edge1;
	vec3 edge2(this->v3.xyz.x - this->v1.xyz.x, this->v3.xyz.y - this->v1.xyz.y, this->v3.xyz.z - this->v1.xyz.z);
	this->edge2 = edge2;
	this->n = edge1.cross(edge2);
	this->n.normalize();

	this->e2 = this->edge2.cross(this->n);
	this->e3 = this->edge1.cross(this->n);

	double dot2 = this->e2.dot(this->edge1);
	double dot3 = this->e3.dot(this->edge2);

	this->e2.x /= dot2;
	this->e2.y /= dot2;
	this->e2.z /= dot2;

	this->e3.x /= dot3;
	this->e3.y /= dot3;
	this->e3.z /= dot3;
}

bool tri::hit(ray *incoming_ray, std::vector<object*> &objects, std::vector<light*> &lights, double color[4], double &distance) {
	double cosineNE = this->n.dot(incoming_ray->direction);
	if (cosineNE == 0) {
		return false;
	}
	vec3 rayToPoint(this->v1.xyz.x - incoming_ray->origin.x, this->v1.xyz.y - incoming_ray->origin.y, this->v1.xyz.z - incoming_ray->origin.z);
	vec3 normal = this->n;
        double unitDir = incoming_ray->direction.dot(normal);
	double pointDotDir = normal.dot(rayToPoint);


	double stepsToTake = pointDotDir / unitDir;
	if (stepsToTake < 0) {
		return false;
	}
	double hitX = incoming_ray->origin.x + stepsToTake * incoming_ray->direction.x;
	double hitY = incoming_ray->origin.y + stepsToTake * incoming_ray->direction.y;
	double hitZ = incoming_ray->origin.z + stepsToTake * incoming_ray->direction.z;
	vec3 hitPoint(hitX, hitY, hitZ);
	
	/*if (incoming_ray->direction.x * hitX < 0 || incoming_ray->direction.y * hitY < 0 || incoming_ray->direction.z * hitZ < 0) {
		return false;
	}*/

	vec3 x1(hitX - this->v1.xyz.x, hitY - this->v1.xyz.y, hitZ - this->v1.xyz.z);

	double b2 = x1.dot(this->e2);
	double b3 = x1.dot(this->e3);
	double b1 = 1 - b2 - b3;

	if (b1 < 0 || b2 < 0 || b3 < 0 || b1 > 1 || b2 > 1 || b3 > 1) {
		return false;
	}

	double distanceFromRay = hitPoint.magnitude();
	if (distanceFromRay >= distance) {
		return false;
	}
	distance = distanceFromRay;
	this->barycentric.x = b1;
	this->barycentric.y = b2;
	this->barycentric.z = b3;
	return true;
}

bool tri::shadowHit(ray *incoming_ray, light* target_light, double &distance) {
	double cosineNE = this->n.dot(incoming_ray->direction);
	if (cosineNE == 0) {
		return false;
	}
	vec3 rayToPoint(this->v1.xyz.x - incoming_ray->origin.x, this->v1.xyz.y - incoming_ray->origin.y, this->v1.xyz.z - incoming_ray->origin.z);
	vec3 normal = this->n;
        double unitDir = incoming_ray->direction.dot(normal);
	double pointDotDir = normal.dot(rayToPoint);


	double stepsToTake = pointDotDir / unitDir;
	if (stepsToTake < 0) {
		return false;
	}
	double hitX = incoming_ray->origin.x + stepsToTake * incoming_ray->direction.x;
	double hitY = incoming_ray->origin.y + stepsToTake * incoming_ray->direction.y;
	double hitZ = incoming_ray->origin.z + stepsToTake * incoming_ray->direction.z;
	vec3 hitPoint(hitX, hitY, hitZ);
	/*if (incoming_ray->direction.x * hitX < 0 || incoming_ray->direction.y * hitY < 0 || incoming_ray->direction.z * hitZ < 0) {
		return false;
	}*/

	vec3 x1(hitX - this->v1.xyz.x, hitY - this->v1.xyz.y, hitZ - this->v1.xyz.z);

	double b2 = x1.dot(this->e2);
	double b3 = x1.dot(this->e3);
	double b1 = 1 - b2 - b3;
	
	if (b1 < 0 || b2 < 0 || b3 < 0 || b1 > 1 || b2 > 1 || b3 > 1) {
		return false;
	}

	double distanceFromRay = hitPoint.magnitude();
	if (distanceFromRay < 0.0000001) {
		return false;
	}
	if (distanceFromRay >= distance) {
		return false;
	}
	distance = distanceFromRay;
	return true;
}

vec3 tri::getColor() {
	
	return vec3(this->c[0], this->c[1], this->c[2]);
}

vec3 tri::getNormal(double x, double y, double z) {
	if (this->barycentric.x == 0 && this->barycentric.y == 0 && this->barycentric.z == 0) {
		return this->n;
	}
	else {
		x = this->barycentric.x * this->v1.normal.x + this->barycentric.y * this->v2.normal.x + this->barycentric.z * this->v3.normal.x;
		y = this->barycentric.x * this->v1.normal.y + this->barycentric.y * this->v2.normal.y + this->barycentric.z * this->v3.normal.y;
		z = this->barycentric.x * this->v1.normal.z + this->barycentric.y * this->v2.normal.z + this->barycentric.z * this->v3.normal.z;
		vec3 interpolated_normal(x,y,z);
		interpolated_normal.normalize();
		return interpolated_normal;
	}
}

tri::~tri() {

}
