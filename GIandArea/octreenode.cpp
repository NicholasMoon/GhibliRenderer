// octreenode.cpp

#include "octreenode.h"

OctreeNode::OctreeNode(int num_obj_cutoff, int level, int position) {
	this->num_obj_cutoff = num_obj_cutoff;
	this->level = level;
	this->position = position;
}

OctreeNode::OctreeNode(AABB *bounding_box, int num_obj_cutoff, int level, int position) {
	this->bounding_box = bounding_box;
	this->num_obj_cutoff = num_obj_cutoff;
	this->level = level;
	this->position = position;
}

void OctreeNode::addChildNode (int position, std::vector<object*> &objects) {
	vec3 child_min_coordinates(0,0,0);
	getMinCoordinates(position, child_min_coordinates);

	vec3 child_max_coordinates(0,0,0);
	getMaxCoordinates(position, child_max_coordinates);
	//std::cout << child_min_coordinates.x << " " << child_min_coordinates.y << " " << child_min_coordinates.z << " " << child_max_coordinates.x << " " << child_max_coordinates.y << " " << child_max_coordinates.z << std::endl;

	AABB *child_bounding_box = new AABB(child_min_coordinates.x, child_min_coordinates.y, child_min_coordinates.z, child_max_coordinates.x, child_max_coordinates.y, child_max_coordinates.z);
	std::vector<object*> children_objects;
	int num_objects_in_box = 0;
	for (int i = 0; i < objects.size(); i++) {
		if (objects[i]->in_bounding_box(child_bounding_box)) {
			//std::cout << "object " << i << " is in box " << position << std::endl;
			num_objects_in_box++;
			children_objects.push_back(objects[i]);
		}
	}
	if (num_objects_in_box > this->num_obj_cutoff) {
		//std::cout << "Making New Box!" << std::endl;
		OctreeNode *child = new OctreeNode(child_bounding_box, this->num_obj_cutoff, this->level + 1, position);
		for (int i = 0; i < 8; i++) {
			child->addChildNode(i, children_objects);
		}
		this->children[position] = child;
		this->child_valid[position] = 1;
	}
	else {
		//std::cout << "Not enough for new box!" << std::endl;
		delete child_bounding_box;
		for (int i = 0; i < children_objects.size(); i++) {
			this->node_objects.push_back(children_objects[i]);
		}
	}
}

void OctreeNode::buildOctree(std::vector<object*> &objects) {
	vec3 min_coordinates(-5, -5, -5);
	vec3 max_coordinates(5, 5, 5);
	/*for (int i = 0; i < objects.size(); i++) {
		if (objects[i]->v1.x < min_coordinates.x) {
			min_coordinates.x = objects[i]->v1.x;
		}
		if (objects[i]->v2.x < min_coordinates.x) {
			min_coordinates.x = objects[i]->v2.x;
		}
		if (objects[i]->v3.x < min_coordinates.x) {
			min_coordinates.x = objects[i]->v3.x;
		}

		if (objects[i]->v1.x > max_coordinates.x) {
			max_coordinates.x = objects[i]->v1.x;
		}
		if (objects[i]->v2.x > max_coordinates.x) {
			max_coordinates.x = objects[i]->v2.x;
		}
		if (objects[i]->v3.x > max_coordinates.x) {
			max_coordinates.x = objects[i]->v3.x;
		}
		

		if (objects[i]->v1.y < min_coordinates.y) {
			min_coordinates.y = objects[i]->v1.y;
		}
		if (objects[i]->v2.y < min_coordinates.y) {
			min_coordinates.y = objects[i]->v2.y;
		}
		if (objects[i]->v3.y < min_coordinates.y) {
			min_coordinates.y = objects[i]->v3.y;
		}

		if (objects[i]->v1.y > max_coordinates.y) {
			max_coordinates.y = objects[i]->v1.y;
		}
		if (objects[i]->v2.y > max_coordinates.y) {
			max_coordinates.x = objects[i]->v2.x;
		}
		if (objects[i]->v3.y > max_coordinates.y) {
			max_coordinates.y = objects[i]->v3.y;
		}

		if (objects[i]->v1.z < min_coordinates.z) {
			min_coordinates.z = objects[i]->v1.z;
		}
		if (objects[i]->v2.z < min_coordinates.z) {
			min_coordinates.z = objects[i]->v2.z;
		}
		if (objects[i]->v3.y < min_coordinates.y) {
			min_coordinates.z = objects[i]->v3.y;
		}

		if (objects[i]->v1.z > max_coordinates.z) {
			max_coordinates.z = objects[i]->v1.z;
		}
		if (objects[i]->v2.z > max_coordinates.z) {
			max_coordinates.z = objects[i]->v2.z;
		}
		if (objects[i]->v3.z > max_coordinates.z) {
			max_coordinates.z = objects[i]->v3.z;
		}
	}*/

	AABB *scene_bounding_box = new AABB(min_coordinates.x, min_coordinates.y, min_coordinates.z, max_coordinates.x, max_coordinates.y, max_coordinates.z);
	this->bounding_box = scene_bounding_box;
	
	for (int i = 0; i < 8; i++) {
		addChildNode(i, objects);
	}

}

bool OctreeNode::traverseOctree(ray *r, std::vector<object*> *objects) {
	int valid_children = 0;
	if (this->bounding_box->intersect(r)) {
		for (int i = 0; i < 8; i++) {
			if (child_valid[i]) {
				valid_children++;
				if (children[i]->traverseOctree(r, objects)) {
					break;
				}
			}
		}
		if (valid_children == 0) {
			objects = &this->node_objects;
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
			min_coordinates.x = this->bounding_box->max_coordinates.x;
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
			min_coordinates.x = this->bounding_box->max_coordinates.x;
			min_coordinates.y = (this->bounding_box->max_coordinates.y + this->bounding_box->min_coordinates.y) / 2;
			min_coordinates.z = (this->bounding_box->max_coordinates.z + this->bounding_box->min_coordinates.z) / 2;
			break;
		default:
			break;
	}
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
}

void OctreeNode::printOctree() {
	std::cout << "Level: " << this->level << " Position: " << this->position << " Objects: " << this->node_objects.size() << std::endl;
	for (int i = 0; i < 8; i++) {
		if (this->child_valid[i]) {
			this->children[i]->printOctree();
		}
	}
}

OctreeNode::~OctreeNode() {

}
