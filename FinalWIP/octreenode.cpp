// octreenode.cpp

#include "octreenode.h"

OctreeNode::OctreeNode(int num_obj_cutoff, int level, int position, int num_objs, int max_depth) {
	this->num_obj_cutoff = num_obj_cutoff;
	this->level = level;
	this->position = position;
	this->num_objs = num_objs;
	this->max_depth = max_depth;
}

OctreeNode::OctreeNode(AABB *bounding_box, int num_obj_cutoff, int level, int position, int num_objs, int max_depth) {
	this->bounding_box = bounding_box;
	this->num_obj_cutoff = num_obj_cutoff;
	this->level = level;
	this->position = position;
	this->num_objs = num_objs;
	this->max_depth = max_depth;
}

bool OctreeNode::addChildren() {
	if (this->level >= this->max_depth || this->node_objects.size() < this->num_obj_cutoff) {
		//std::cout << this->level << " " << this->max_depth << std::endl;
		return false;
	}

	for (int child_number = 0; child_number < 8; child_number++) {
		vec3 child_min_coordinates(0,0,0);
		getMinCoordinates(child_number, child_min_coordinates);

		vec3 child_max_coordinates(0,0,0);
		getMaxCoordinates(child_number, child_max_coordinates);

		AABB *child_bounding_box = new AABB(child_min_coordinates.x, child_min_coordinates.y, child_min_coordinates.z, child_max_coordinates.x, child_max_coordinates.y, child_max_coordinates.z);
		OctreeNode *child = new OctreeNode(child_bounding_box, this->num_obj_cutoff, this->level + 1, child_number, this->num_objs, this->max_depth);
		for (int i = 0; i < this->node_objects.size(); i++) {
			if (this->node_objects[i]->in_bounding_box(child_bounding_box)) {
				child->node_objects.push_back(this->node_objects[i]);
			}
		}
		this->children[child_number] = child;
		this->child_valid[child_number] = 1;
	}
	for (int j = 0; j < 8; j++) {
		this->children[j]->addChildren();
	}
	this->node_objects.clear();
	return true;
}

void OctreeNode::buildOctree(std::vector<object*> &objects) {
	double max_distance = std::numeric_limits<double>::max();
	double min_distance = -max_distance;
	vec3 min_coordinates(max_distance,max_distance,max_distance);
	vec3 max_coordinates(min_distance,min_distance,min_distance);

	getWorldBoundaries(min_coordinates, max_coordinates, objects);
	min_coordinates.x -= 0.001;
	min_coordinates.y -= 0.001;
	min_coordinates.z -= 0.001;
	max_coordinates.x += 0.001;
	max_coordinates.y += 0.001;
	max_coordinates.z += 0.001;
	AABB *scene_bounding_box = new AABB(min_coordinates.x, min_coordinates.y, min_coordinates.z, max_coordinates.x, max_coordinates.y, max_coordinates.z);
	fflush(stdout);
	this->bounding_box = scene_bounding_box;
	for (int j = 0; j < objects.size(); j++) {
		this->node_objects.push_back(objects[j]);
	}
	addChildren();
}

bool OctreeNode::traverseOctree(ray *r, std::vector<std::vector<object*>> &objects) {
	int valid_children = 0;
	if (this->bounding_box->intersect(r)) {
		if (this->node_objects.size() > 0) {
			objects.push_back(this->node_objects);
			return true;
		}
		for (int k = 0; k < 8; k++) {
			valid_children += child_valid[k];
		}
		if (valid_children == 0 && this->node_objects.size() == 0) {
			return true;
		}
		for (int i = 0; i < 8; i++) {
			children[i]->traverseOctree(r, objects);
		}
		return true;
	}
	else {
		return false;
	}
	
}

void OctreeNode::getMinCoordinates(int position, vec3 &min_coordinates) {
	switch (position) {
		case (0):
			min_coordinates.x = this->bounding_box->min_coordinates.x;
			min_coordinates.y = this->bounding_box->min_coordinates.y;
			min_coordinates.z = this->bounding_box->min_coordinates.z;
			break;
		case (1):
			min_coordinates.x = (this->bounding_box->max_coordinates.x + this->bounding_box->min_coordinates.x) / 2;
			min_coordinates.y = this->bounding_box->min_coordinates.y;
			min_coordinates.z = this->bounding_box->min_coordinates.z;
			break;
		case (2):
			min_coordinates.x = (this->bounding_box->max_coordinates.x + this->bounding_box->min_coordinates.x) / 2;
			min_coordinates.y = this->bounding_box->min_coordinates.y;
			min_coordinates.z = (this->bounding_box->max_coordinates.z + this->bounding_box->min_coordinates.z) / 2;
			break;
		case (3):
			min_coordinates.x = this->bounding_box->min_coordinates.x;
			min_coordinates.y = this->bounding_box->min_coordinates.y;
			min_coordinates.z = (this->bounding_box->max_coordinates.z + this->bounding_box->min_coordinates.z) / 2;
			break;
		case (4):
			min_coordinates.x = this->bounding_box->min_coordinates.x;
			min_coordinates.y = (this->bounding_box->max_coordinates.y + this->bounding_box->min_coordinates.y) / 2;
			min_coordinates.z = this->bounding_box->min_coordinates.z;
			break;
		case (5):
			min_coordinates.x = (this->bounding_box->max_coordinates.x + this->bounding_box->min_coordinates.x) / 2;
			min_coordinates.y = (this->bounding_box->max_coordinates.y + this->bounding_box->min_coordinates.y) / 2;
			min_coordinates.z = this->bounding_box->min_coordinates.z;
			break;
		case (6):
			min_coordinates.x = (this->bounding_box->max_coordinates.x + this->bounding_box->min_coordinates.x) / 2;
			min_coordinates.y = (this->bounding_box->max_coordinates.y + this->bounding_box->min_coordinates.y) / 2;
			min_coordinates.z = (this->bounding_box->max_coordinates.z + this->bounding_box->min_coordinates.z) / 2;
			break;
		case (7):
			min_coordinates.x = this->bounding_box->min_coordinates.x;
			min_coordinates.y = (this->bounding_box->max_coordinates.y + this->bounding_box->min_coordinates.y) / 2;
			min_coordinates.z = (this->bounding_box->max_coordinates.z + this->bounding_box->min_coordinates.z) / 2;
			break;
		default:
			break;
	}
	/*min_coordinates.x -= 0.001;
	min_coordinates.y -= 0.001;
	min_coordinates.z -= 0.001;*/
}

void OctreeNode::getMaxCoordinates(int position, vec3 &max_coordinates) {
	switch (position) {
		case (0):
			max_coordinates.x = (this->bounding_box->max_coordinates.x + this->bounding_box->min_coordinates.x) / 2;
			max_coordinates.y = (this->bounding_box->max_coordinates.y + this->bounding_box->min_coordinates.y) / 2;
			max_coordinates.z = (this->bounding_box->max_coordinates.z + this->bounding_box->min_coordinates.z) / 2;
			break;
		case (1):
			max_coordinates.x = this->bounding_box->max_coordinates.x;
			max_coordinates.y = (this->bounding_box->max_coordinates.y + this->bounding_box->min_coordinates.y) / 2;
			max_coordinates.z = (this->bounding_box->max_coordinates.z + this->bounding_box->min_coordinates.z) / 2;
			break;
		case (2):
			max_coordinates.x = this->bounding_box->max_coordinates.x;
			max_coordinates.y = (this->bounding_box->max_coordinates.y + this->bounding_box->min_coordinates.y) / 2;
			max_coordinates.z = this->bounding_box->max_coordinates.z;
			break;
		case (3):
			max_coordinates.x = (this->bounding_box->max_coordinates.x + this->bounding_box->min_coordinates.x) / 2;
			max_coordinates.y = (this->bounding_box->max_coordinates.y + this->bounding_box->min_coordinates.y) / 2;
			max_coordinates.z = this->bounding_box->max_coordinates.z;
			break;
		case (4):
			max_coordinates.x = (this->bounding_box->max_coordinates.x + this->bounding_box->min_coordinates.x) / 2;
			max_coordinates.y = this->bounding_box->max_coordinates.y;
			max_coordinates.z = (this->bounding_box->max_coordinates.z + this->bounding_box->min_coordinates.z) / 2;
			break;
		case (5):
			max_coordinates.x = this->bounding_box->max_coordinates.x;
			max_coordinates.y = this->bounding_box->max_coordinates.y;
			max_coordinates.z = (this->bounding_box->max_coordinates.z + this->bounding_box->min_coordinates.z) / 2;
			break;
		case (6):
			max_coordinates.x = this->bounding_box->max_coordinates.x;
			max_coordinates.y = this->bounding_box->max_coordinates.y;
			max_coordinates.z = this->bounding_box->max_coordinates.z;
			break;
		case (7):
			max_coordinates.x = (this->bounding_box->max_coordinates.x + this->bounding_box->min_coordinates.x) / 2;
			max_coordinates.y = this->bounding_box->max_coordinates.y;
			max_coordinates.z = this->bounding_box->max_coordinates.z;
			break;
		default:
			break;
	}
	/*max_coordinates.x += 0.001;
	max_coordinates.y += 0.001;
	max_coordinates.z += 0.001;*/
}

int OctreeNode::printOctree() {
	std::cout << "Level: " << this->level << " Position: " << this->position << " Objects: " << this->node_objects.size() << " AABB: ";
	this->bounding_box->PrintAABB();
	int total_objects = this->node_objects.size();
	for (int i = 0; i < 8; i++) {
		if (this->child_valid[i]) {
			total_objects += this->children[i]->printOctree();
		}
	}
	return total_objects;
}

OctreeNode::~OctreeNode() {

}
