// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// ray.cpp

#include "ray.h"

bool ray::cast(std::vector<object*> &objects, std::vector<light*> &lights, double color[4], int bounces, int lastObject, std::default_random_engine &generator, int x, int y, int width, double *depthMap, int &primary_objID, vec3 &hit_normal, int &object_type, std::vector<int> &hit_list, int &shadowed) {
	double distance = std::numeric_limits<double>::max();
	int closestObject = -1;
	double diffuse_color[3] = {0,0,0};
	double reflection_color[4] = {0,0,0,0};
	double refraction_color[4] = {0,0,0,0};
	for (int j = 0; j < objects.size(); j++) {
		if (objects[j]->hit(this, objects, lights, color, distance)) {
			closestObject = j;	
		}
	}

	if (closestObject == -1) {
		// ray didn't hit any objects
		color[3] = 0;
		return false;
	}
	primary_objID = objects[closestObject]->objectID;
	hit_list.push_back(objects[closestObject]->objectID);
	if (x != -1 && y != -1) {
		depthMap[y * width + x] = distance;
	}

	object_type = objects[closestObject]->object_type;
  	std::normal_distribution<double> distribution(0, objects[closestObject]->mat->roughness);
	double rand_x = distribution(generator);
	double rand_y = distribution(generator);
	double rand_z = distribution(generator);

	double hitX = this->origin.x + distance * this->direction.x;
	double hitY = this->origin.y + distance * this->direction.y;
	double hitZ = this->origin.z + distance * this->direction.z;
	vec3 hitPoint(hitX, hitY, hitZ);
	vec3 normal = objects[closestObject]->getNormal(hitX, hitY, hitZ);
	normal.normalize();
	hit_normal.x = normal.x;
	hit_normal.y = normal.y;
	hit_normal.z = normal.z;
	double cosineNE = normal.dot(this->direction);
	if (cosineNE > 0) {
		normal.reverse();
	}

	if ((objects[closestObject]->mat->shininess.x > 0 || objects[closestObject]->mat->shininess.y > 0 || objects[closestObject]->mat->shininess.z > 0) && bounces > 0) {
		if (objects[closestObject]->mat->roughness > 0) {
			normal.x += rand_x;
			normal.y += rand_y;
			normal.z += rand_z;
			normal.normalize();
		}
		double n_dot_i = normal.dot(this->direction);
		double reflect_dir_x = this->direction.x - 2 * (n_dot_i) * normal.x;
		double reflect_dir_y = this->direction.y - 2 * (n_dot_i) * normal.y;
		double reflect_dir_z = this->direction.z - 2 * (n_dot_i) * normal.z;
		ray *reflect_ray = new ray(hitX, hitY, hitZ, reflect_dir_x, reflect_dir_y, reflect_dir_z);
		reflect_ray->cast(objects, lights, reflection_color, bounces - 1, -1, generator, x, y, width, depthMap, primary_objID, hit_normal, object_type, hit_list, shadowed);
		delete reflect_ray;
	}
	if ((objects[closestObject]->mat->transparency.x > 0 || objects[closestObject]->mat->transparency.y > 0 || objects[closestObject]->mat->transparency.z > 0) && bounces > 0) {
		double ior = objects[closestObject]->mat->ior;
		if (objects[closestObject]->mat->roughness > 0) {
			normal.x += rand_x;
			normal.y += rand_y;
			normal.z += rand_z;
			normal.normalize();
		}
		double n_dot_i = normal.dot(this->direction);
		if (n_dot_i < 0) {
			n_dot_i = -n_dot_i;
			ior = 1 / ior;
		}
		else {
			normal.reverse();
		}
		double k = 1.0 - (ior * ior) * (1.0 - n_dot_i * n_dot_i);
		double refract_dir_x;
		double refract_dir_y;
		double refract_dir_z;
		if (k < 0) {
			normal = objects[closestObject]->getNormal(hitX, hitY, hitZ);
			normal.normalize();
			if (objects[closestObject]->mat->roughness > 0) {
				normal.x += rand_x;
				normal.y += rand_y;
				normal.z += rand_z;
				normal.normalize();
			}
			n_dot_i = normal.dot(this->direction);
			refract_dir_x = this->direction.x - 2 * (n_dot_i) * normal.x;
			refract_dir_y = this->direction.y - 2 * (n_dot_i) * normal.y;
			refract_dir_z = this->direction.z - 2 * (n_dot_i) * normal.z;
			ray *refract_ray = new ray(hitX, hitY, hitZ, refract_dir_x, refract_dir_y, refract_dir_z);
			refract_ray->cast(objects, lights, refraction_color, bounces - 1, -1, generator, -1, -1, width, depthMap, primary_objID, hit_normal, object_type, hit_list, shadowed);
			delete refract_ray;
		}
		else {
			refract_dir_x = ior * this->direction.x + (ior * n_dot_i - sqrt(k)) * normal.x;
			refract_dir_y = ior * this->direction.y + (ior * n_dot_i - sqrt(k)) * normal.y;
			refract_dir_z = ior * this->direction.z + (ior * n_dot_i - sqrt(k)) * normal.z;
			ray *refract_ray = new ray(hitX + (.0001) * refract_dir_x, hitY + (.0001) * refract_dir_y, hitZ + (.0001) * refract_dir_z, refract_dir_x, refract_dir_y, refract_dir_z);
			refract_ray->cast(objects, lights, refraction_color, bounces - 1, closestObject, generator, -1, -1, width, depthMap, primary_objID, hit_normal, object_type, hit_list, shadowed);
			delete refract_ray;
		}
	}
	if (lastObject == closestObject) {
		// don't count diffuse color for refraction exit ray of object
	}
	else if (lights.size() > 0) {
		ray *light_ray;
		double distanceToLight;
		double distanceToBlock;
		for (int i = 0; i < lights.size(); i++) {
			if (lights[i]->type == 1) { // sun light

				distanceToLight = std::numeric_limits<double>::max();

				light_ray = new ray(hitX, hitY, hitZ, lights[i]->x, lights[i]->y, lights[i]->z);
				distanceToBlock = light_ray->castLight(objects, lights[i], distanceToLight);
				delete light_ray;
				if (distanceToBlock < distanceToLight) {
					diffuse_color[0] += 0;
					diffuse_color[1] += 0;
					diffuse_color[2] += 0;
					continue;
				}
				if (objects[closestObject]->mat->roughness > 0) {
					normal.x += rand_x;
					normal.y += rand_y;
					normal.z += rand_z;
					normal.normalize();
				}
				double cosineNL = normal.x * lights[i]->x + normal.y * lights[i]->y + normal.z * lights[i]->z;
				if (cosineNL > 0) {
					vec3 objectColor = objects[closestObject]->getColor();
					if (cosineNL < 0.000000000001) {
						color[0] += objectColor.x * lights[i]->c[0] / 6.0;
						color[1] += objectColor.y * lights[i]->c[1] / 6.0;
						color[2] += objectColor.z * lights[i]->c[2] / 6.0;
					}
					else if (cosineNL < 0.35) {
						color[0] += objectColor.x * lights[i]->c[0] / 2.0;
						color[1] += objectColor.y * lights[i]->c[1] / 2.0;
						color[2] += objectColor.z * lights[i]->c[2] / 2.0;
					}
					else{
						color[0] += objectColor.x * lights[i]->c[0];
						color[1] += objectColor.y * lights[i]->c[1];
						color[2] += objectColor.z * lights[i]->c[2];
					}
				}
				
			}
			else { // bulb light
				vec3 lightCoordinate(lights[i]->x, lights[i]->y, lights[i]->z);
				
				distanceToLight = hitPoint.distance(lightCoordinate);

				vec3 directionToLight(lightCoordinate.x - hitPoint.x, lightCoordinate.y - hitPoint.y, lightCoordinate.z - hitPoint.z);
				directionToLight.normalize();

				light_ray = new ray(hitX, hitY, hitZ, directionToLight.x, directionToLight.y, directionToLight.z);
				distanceToBlock = light_ray->castLight(objects, lights[i], distanceToLight);
				delete light_ray;
				vec3 objectColor = objects[closestObject]->getColor();
				if (distanceToBlock < distanceToLight) {
					if (objects[closestObject]->object_type == 1) {
						diffuse_color[0] += objectColor.x * lights[i]->c[0] / 6.0;
						diffuse_color[1] += objectColor.y * lights[i]->c[1] / 6.0;
						diffuse_color[2] += objectColor.z * lights[i]->c[2] / 6.0;
						continue;
					}
					else {
						shadowed = 1;
					}
				}
				if (objects[closestObject]->mat->roughness > 0) {
					normal.x += rand_x;
					normal.y += rand_y;
					normal.z += rand_z;
					normal.normalize();
				}
				double cosineNL = normal.x * directionToLight.x + normal.y * directionToLight.y + normal.z * directionToLight.z;
				if (objects[closestObject]->object_type == 1) {
					if (cosineNL < 0.000000000001) {
						diffuse_color[0] += objectColor.x * lights[i]->c[0] / 6.0;
						diffuse_color[1] += objectColor.y * lights[i]->c[1] / 6.0;
						diffuse_color[2] += objectColor.z * lights[i]->c[2] / 6.0;
					}
					else if (cosineNL < 0.35) {
						diffuse_color[0] += objectColor.x * lights[i]->c[0] / 2.0;
						diffuse_color[1] += objectColor.y * lights[i]->c[1] / 2.0;
						diffuse_color[2] += objectColor.z * lights[i]->c[2] / 2.0;
					}
					else{
						diffuse_color[0] += objectColor.x * lights[i]->c[0];
						diffuse_color[1] += objectColor.y * lights[i]->c[1];
						diffuse_color[2] += objectColor.z * lights[i]->c[2];
					}
				}
				else {
					diffuse_color[0] += (objectColor.x * lights[i]->c[0] * cosineNL) / (distanceToLight * distanceToLight);
					diffuse_color[1] += (objectColor.y * lights[i]->c[1] * cosineNL) / (distanceToLight * distanceToLight);
					diffuse_color[2] += (objectColor.z * lights[i]->c[2] * cosineNL) / (distanceToLight * distanceToLight);
				}
				
			}
		}
	}
	else { // default light (when there is no user lighting)
		vec3 objectColor = objects[closestObject]->getColor();
		diffuse_color[0] = objectColor.x;
		diffuse_color[1] = objectColor.y;
		diffuse_color[2] = objectColor.z;
	}
		color[0] += diffuse_color[0] * objects[closestObject]->mat->diffuse.x + reflection_color[0] * objects[closestObject]->mat->shininess.x + refraction_color[0] * objects[closestObject]->mat->transparency.x;
		color[1] += diffuse_color[1] * objects[closestObject]->mat->diffuse.y + reflection_color[1] * objects[closestObject]->mat->shininess.y + refraction_color[1] * objects[closestObject]->mat->transparency.y;
		color[2] += diffuse_color[2] * objects[closestObject]->mat->diffuse.z + reflection_color[2] * objects[closestObject]->mat->shininess.z + refraction_color[2] * objects[closestObject]->mat->transparency.z;
	color[3] += 1;
	return true;
}

double ray::castLight(std::vector<object*> &objects, light *target_light, double distance) {
	for (int i = 0; i < objects.size(); i++) {
		objects[i]->shadowHit(this, target_light, distance);
	}
	return distance;
}

double ray::detect_edge(std::vector<object*> &objects, std::vector<light*> &lights, double color[4], int bounces, int lastObject, std::default_random_engine &generator, int &stencil_objID, int &object_type, std::vector<int> &hit_list) {
	double distance = std::numeric_limits<double>::max();
	int closestObject = -1;
	double diffuse_color[3] = {0,0,0};
	double reflection_color[4] = {0,0,0,0};
	double refraction_color[4] = {0,0,0,0};
	for (int j = 0; j < objects.size(); j++) {
		if (objects[j]->hit(this, objects, lights, color, distance)) {
			closestObject = j;	
		}
	}
	if (closestObject != -1) {
		stencil_objID = objects[closestObject]->objectID;
	}
	object_type = objects[closestObject]->object_type;
	hit_list.push_back(objects[closestObject]->objectID);
	std::normal_distribution<double> distribution(0, objects[closestObject]->mat->roughness);
	double rand_x = distribution(generator);
	double rand_y = distribution(generator);
	double rand_z = distribution(generator);
	if ((objects[closestObject]->mat->shininess.x > 0 || objects[closestObject]->mat->shininess.y > 0 || objects[closestObject]->mat->shininess.z > 0) && bounces > 0) {
		double hitX = this->origin.x + distance * this->direction.x;
		double hitY = this->origin.y + distance * this->direction.y;
		double hitZ = this->origin.z + distance * this->direction.z;
		vec3 normal = objects[closestObject]->getNormal(hitX, hitY, hitZ);
		normal.normalize();
		if (objects[closestObject]->mat->roughness > 0) {
			normal.x += rand_x;
			normal.y += rand_y;
			normal.z += rand_z;
			normal.normalize();
		}
		double n_dot_i = normal.dot(this->direction);
		double reflect_dir_x = this->direction.x - 2 * (n_dot_i) * normal.x;
		double reflect_dir_y = this->direction.y - 2 * (n_dot_i) * normal.y;
		double reflect_dir_z = this->direction.z - 2 * (n_dot_i) * normal.z;
		ray *reflect_ray = new ray(hitX, hitY, hitZ, reflect_dir_x, reflect_dir_y, reflect_dir_z);
		vec3 hitnormal(0,0,0);
		distance = reflect_ray->detect_edge(objects, lights, color, bounces, -1, generator, stencil_objID, object_type, hit_list);
		delete reflect_ray;
	}
	return distance;
}

ray::~ray() {

}
