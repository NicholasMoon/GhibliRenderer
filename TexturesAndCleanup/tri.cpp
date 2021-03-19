// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// tri.cpp

#include "tri.h"


tri::tri(vertex *v1, vertex *v2, vertex *v3, double c[3]) {
	vertex vert1(v1->xyz, v1->uv, v1->normal, v1->color);
	vertex vert2(v2->xyz, v2->uv, v2->normal, v2->color);
	vertex vert3(v3->xyz, v3->uv, v3->normal, v3->color);
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
	
	this->box_pointers = 0;
}

tri::tri(vertex *v1, vertex *v2, vertex *v3, double c[3], material *mat, int objectID) {
	this->objectID = objectID;
	vertex vert1(v1->xyz, v1->uv, v1->normal, v1->color);
	vertex vert2(v2->xyz, v2->uv, v2->normal, v2->color);
	vertex vert3(v3->xyz, v3->uv, v3->normal, v3->color);
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

	this->box_pointers = 0;
}

tri::tri(vertex *v1, vertex *v2, vertex *v3, double c[3], material *mat, int objectID, int object_type) {
	this->objectID = objectID;
	this->object_type = object_type;
	vertex vert1(v1->xyz, v1->uv, v1->normal, v1->color);
	vertex vert2(v2->xyz, v2->uv, v2->normal, v2->color);
	vertex vert3(v3->xyz, v3->uv, v3->normal, v3->color);
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
	
	this->box_pointers = 0;
}

tri::tri(vertex *v1, vertex *v2, vertex *v3, double c[3], double e[3], material *mat, int objectID, int object_type) {
	this->objectID = objectID;
	this->object_type = object_type;
	vertex vert1(v1->xyz, v1->uv, v1->normal, v1->color);
	vertex vert2(v2->xyz, v2->uv, v2->normal, v2->color);
	vertex vert3(v3->xyz, v3->uv, v3->normal, v3->color);
	this->v1 = vert1;
	this->v2 = vert2;
	this->v3 = vert3;
	this->c[0] = c[0];
	this->c[1] = c[1];
	this->c[2] = c[2];
	this->e[0] = e[0];
	this->e[1] = e[1];
	this->e[2] = e[2];
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

	this->box_pointers = 0;
}

tri::tri(vertex *v1, vertex *v2, vertex *v3, texture *t, double e[3], material *mat, int objectID, int object_type) { // supports texture color
	this->objectID = objectID;
	this->object_type = object_type;
	vertex vert1(v1->xyz, v1->uv, v1->normal, v1->color);
	vertex vert2(v2->xyz, v2->uv, v2->normal, v2->color);
	vertex vert3(v3->xyz, v3->uv, v3->normal, v3->color);
	this->v1 = vert1;
	this->v2 = vert2;
	this->v3 = vert3;
	this->tex = t;
	this->e[0] = e[0];
	this->e[1] = e[1];
	this->e[2] = e[2];
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

	this->box_pointers = 0;
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

	if (b1 < -0.0000001 || b2 < -0.0000001 || b3 < -0.0000001 || b1 > 1.0000001 || b2 > 1.0000001 || b3 > 1.0000001) {
		return false;
	}
	double travelX = stepsToTake * incoming_ray->direction.x;
	double travelY = stepsToTake * incoming_ray->direction.y;
	double travelZ = stepsToTake * incoming_ray->direction.z;
	vec3 travelPoint(travelX, travelY, travelZ);
	double distanceFromRay = travelPoint.magnitude();
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
	
	if (b1 < -0.0000001 || b2 < -0.0000001 || b3 < -0.0000001 || b1 > 1.0000001 || b2 > 1.0000001 || b3 > 1.0000001) {
		return false;
	}

	double travelX = stepsToTake * incoming_ray->direction.x;
	double travelY = stepsToTake * incoming_ray->direction.y;
	double travelZ = stepsToTake * incoming_ray->direction.z;
	vec3 travelPoint(travelX, travelY, travelZ);
	double distanceFromRay = travelPoint.magnitude();
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

vec3 tri::getColor(double u, double v, vec3 &hit_point) {
	
	return this->tex->getColor(u, v, hit_point);
}

vec3 tri::getEmission() { 
	return vec3(e[0], e[1], e[2]);
}

vec3 tri::getNormal(double x, double y, double z, int flat) {
	if (flat) {
		return this->n;
	}
	else if (this->n.x == 1 || this->n.y == 1 || this->n.z == 1 || this->n.x == -1 || this->n.y == -1 || this->n.z == -1) {
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

vec3 tri::getTextureCoordinates(vec3 &hit_point) {
	vec3 x1(hit_point.x - this->v1.xyz.x, hit_point.y - this->v1.xyz.y, hit_point.z - this->v1.xyz.z);
	double b2 = x1.dot(this->e2);
	double b3 = x1.dot(this->e3);
	double b1 = 1 - b2 - b3;

	double u = this->barycentric.x * this->v1.uv.x + this->barycentric.y * this->v2.uv.x + this->barycentric.z * this->v3.uv.x;
	double v = this->barycentric.x * this->v1.uv.y + this->barycentric.y * this->v2.uv.y + this->barycentric.z * this->v3.uv.y;
	return vec3(u, v, 0);
}

bool tri::SAT_projection(vec3 &axis, vec3 &center, vec3 &extents, vec3 &box_normal_x, vec3 &box_normal_y, vec3 &box_normal_z, vec3 &vc1, vec3 &vc2, vec3 &vc3) {
	double v1_projection = vc1.dot(axis);
	double v2_projection = vc2.dot(axis);
	double v3_projection = vc3.dot(axis);

	double r = extents.x * std::abs(box_normal_x.dot(axis)) + extents.y * std::abs(box_normal_y.dot(axis)) + extents.z * std::abs(box_normal_z.dot(axis));

	if (std::max(-Max_Double(v1_projection, v2_projection, v3_projection), Min_Double(v1_projection, v2_projection, v3_projection)) > r) {
		return false;
	}
	return true;
}

bool tri::in_bounding_box(AABB *bounding_box) {
	// seperating axis theorem -> gdbooks.gitbooks.io/3dcollisions/content/Chapter4/aabb-triangle.html
	vec3 vc1(this->v1.xyz);
	vec3 vc2(this->v2.xyz);
	vec3 vc3(this->v3.xyz);
	vec3 center((bounding_box->min_coordinates.x + bounding_box->max_coordinates.x) / 2, (bounding_box->min_coordinates.y + bounding_box->max_coordinates.y) / 2, (bounding_box->min_coordinates.z + bounding_box->max_coordinates.z) / 2);
	vec3 extents(bounding_box->max_coordinates.x - center.x, bounding_box->max_coordinates.y - center.y, bounding_box->max_coordinates.z - center.z);
	
	vc1.subtractAndSave(center);
	vc2.subtractAndSave(center);
	vc3.subtractAndSave(center);


	vec3 tri_edge1;
	vc2.subtract(vc1,tri_edge1);
	vec3 tri_edge2;
	vc3.subtract(vc2,tri_edge2);
	vec3 tri_edge3;
	vc1.subtract(vc3,tri_edge3);

	

	vec3 box_normal_x(1,0,0);
	vec3 box_normal_y(0,1,0);
	vec3 box_normal_z(0,0,1);

	vec3 edge1_cross_bx, edge1_cross_by, edge1_cross_bz, edge2_cross_bx, edge2_cross_by, edge2_cross_bz, edge3_cross_bx, edge3_cross_by, edge3_cross_bz;

	edge1_cross_bx = box_normal_x.cross(tri_edge1);
	edge1_cross_by = box_normal_y.cross(tri_edge1);
	edge1_cross_bz = box_normal_z.cross(tri_edge1);

	edge2_cross_bx = box_normal_x.cross(tri_edge2);
	edge2_cross_by = box_normal_y.cross(tri_edge2);
	edge2_cross_bz = box_normal_z.cross(tri_edge2);

	edge3_cross_bx = box_normal_x.cross(tri_edge3);
	edge3_cross_by = box_normal_y.cross(tri_edge3);
	edge3_cross_bz = box_normal_z.cross(tri_edge3);

	if (!SAT_projection(edge1_cross_bx, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else if (!SAT_projection(edge1_cross_by, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else if (!SAT_projection(edge1_cross_bz, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else if (!SAT_projection(edge2_cross_bx, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else if (!SAT_projection(edge2_cross_by, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else if (!SAT_projection(edge2_cross_bz, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else if (!SAT_projection(edge3_cross_bx, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else if (!SAT_projection(edge3_cross_by, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else if (!SAT_projection(edge3_cross_bz, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else if (!SAT_projection(box_normal_x, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else if (!SAT_projection(box_normal_y, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else if (!SAT_projection(box_normal_z, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else if (!SAT_projection(this->n, center, extents, box_normal_x, box_normal_y, box_normal_z, vc1, vc2, vc3)) {
		return false;
	}
	else {
		return true;
	}
}

void tri::updateWorldBoundaries(vec3 &min_coordinates, vec3 &max_coordinates) {
	if (this->v1.xyz.x < min_coordinates.x) {
		min_coordinates.x = this->v1.xyz.x;
	}
	if (this->v2.xyz.x < min_coordinates.x) {
		min_coordinates.x = this->v2.xyz.x;
	}
	if (this->v3.xyz.x < min_coordinates.x) {
		min_coordinates.x = this->v3.xyz.x;
	}
	if (this->v1.xyz.x > max_coordinates.x) {
		max_coordinates.x = this->v1.xyz.x;
	}
	if (this->v2.xyz.x > max_coordinates.x) {
		max_coordinates.x = this->v2.xyz.x;
	}
	if (this->v3.xyz.x > max_coordinates.x) {
		max_coordinates.x = this->v3.xyz.x;
	}
	
	if (this->v1.xyz.y < min_coordinates.y) {
		min_coordinates.y = this->v1.xyz.y;
	}
	if (this->v2.xyz.y < min_coordinates.y) {
		min_coordinates.y = this->v2.xyz.y;
	}
	if (this->v3.xyz.y < min_coordinates.y) {
		min_coordinates.y = this->v3.xyz.y;
	}
	if (this->v1.xyz.y > max_coordinates.y) {
		max_coordinates.y = this->v1.xyz.y;
	}
	if (this->v2.xyz.y > max_coordinates.y) {
		max_coordinates.y = this->v2.xyz.y;
	}
	if (this->v3.xyz.y > max_coordinates.y) {
		max_coordinates.y = this->v3.xyz.y;
	}

	if (this->v1.xyz.z < min_coordinates.z) {
		min_coordinates.z = this->v1.xyz.z;
	}
	if (this->v2.xyz.z < min_coordinates.z) {
		min_coordinates.z = this->v2.xyz.z;
	}
	if (this->v3.xyz.z < min_coordinates.z) {
		min_coordinates.z = this->v3.xyz.z;
	}
	if (this->v1.xyz.z > max_coordinates.z) {
		max_coordinates.z = this->v1.xyz.z;
	}
	if (this->v2.xyz.z > max_coordinates.z) {
		max_coordinates.z = this->v2.xyz.z;
	}
	if (this->v3.xyz.z > max_coordinates.z) {
		max_coordinates.z = this->v3.xyz.z;
	}
}

tri::~tri() {

}
