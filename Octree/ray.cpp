// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// ray.cpp

#include "ray.h"

#define PI 3.14159265

bool ray::cast(std::vector<light*> &lights, double color[4], double environmentColor[4], int bounces, int lastObject, std::default_random_engine &generator, int primary_ray, int x, int y, int width, double *depthMap, int &primary_objID, vec3 &hit_normal, int &object_type, std::vector<int> &hit_list, int &shadowed, int flat, int environment, int light_samples, int indirect_samples, int indirect_bounces, OctreeNode *octree) {
	double distance = std::numeric_limits<double>::max();
	int closestObject = -1;
	double direct_diffuse_color[3] = {0,0,0};
	double indirect_diffuse_color[4] = {0,0,0,0};
	double diffuse_color[3] = {0,0,0};
	double reflection_color[4] = {0,0,0,0};
	double refraction_color[4] = {0,0,0,0};
	std::vector<object*> objects;
	std::vector<std::vector<object*>> objs;
	octree->traverseOctree(this, objs);

	for (int r = 0; r < objs.size(); r++) {
		for (int e = 0; e < objs[r].size(); e++) {
			objects.push_back(objs[r][e]);
		}
	}
	int closestObjectID = -1;
	for (int j = 0; j < objects.size(); j++) {
		if (objects[j]->hit(this, objects, lights, color, distance)) {
			closestObject = j;
			closestObjectID = objects[j]->objectID;	
		}
	}
	if (closestObject == -1) {
		// ray didn't hit any objects
		if (environment) {
			if (environmentColor[3] == 1) { // gradient - color to white
				double t = 0.5 * (this->direction.y + 1.0);
				color[0] = (1.0 - t) * 1.0 + t * environmentColor[0];
				color[1] = (1.0 - t) * 1.0 + t * environmentColor[1];
				color[2] = (1.0 - t) * 1.0 + t * environmentColor[2];
				color[3] = 1;
			}
			else if (environmentColor[3] == 0) { // solid color
				color[0] = environmentColor[0];
				color[1] = environmentColor[1];
				color[2] = environmentColor[2];
				color[3] = 1;
			}	
		}
		else {
			color[3] = 0;
		}	
		return false;
	}
	primary_objID = objects[closestObject]->objectID;
	hit_list.push_back(objects[closestObject]->objectID);
	if (x != -1 && y != -1 && primary_ray) {
		depthMap[y * width + x] = distance;
	}
	int previous_obj_type = object_type;
	object_type = objects[closestObject]->object_type;
  	std::normal_distribution<double> distribution(0, objects[closestObject]->mat->roughness);
	double rand_x = distribution(generator);
	double rand_y = distribution(generator);
	double rand_z = distribution(generator);

	double hitX = this->origin.x + distance * this->direction.x;
	double hitY = this->origin.y + distance * this->direction.y;
	double hitZ = this->origin.z + distance * this->direction.z;
	vec3 hitPoint(hitX, hitY, hitZ);
	vec3 normal = objects[closestObject]->getNormal(hitX, hitY, hitZ, flat);
	normal.normalize();
	hit_normal.x = normal.x;
	hit_normal.y = normal.y;
	hit_normal.z = normal.z;
	double cosineNE = normal.dot(this->direction);
	if (cosineNE > 0) {
		normal.reverse();
	}

	vec3 objectEmission = objects[closestObject]->getEmission();

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
		reflect_ray->cast(lights, reflection_color, environmentColor, bounces - 1, -1, generator, primary_ray, x, y, width, depthMap, primary_objID, hit_normal, object_type, hit_list, shadowed, flat, environment, light_samples, indirect_samples, indirect_bounces, octree);
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
			normal = objects[closestObject]->getNormal(hitX, hitY, hitZ, flat);
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
			refract_ray->cast(lights, refraction_color, environmentColor, bounces - 1, -1, generator, primary_ray, -1, -1, width, depthMap, primary_objID, hit_normal, object_type, hit_list, shadowed, flat, environment, light_samples, indirect_samples, indirect_bounces, octree);
			delete refract_ray;
		}
		else {
			refract_dir_x = ior * this->direction.x + (ior * n_dot_i - sqrt(k)) * normal.x;
			refract_dir_y = ior * this->direction.y + (ior * n_dot_i - sqrt(k)) * normal.y;
			refract_dir_z = ior * this->direction.z + (ior * n_dot_i - sqrt(k)) * normal.z;
			ray *refract_ray = new ray(hitX + (.0001) * refract_dir_x, hitY + (.0001) * refract_dir_y, hitZ + (.0001) * refract_dir_z, refract_dir_x, refract_dir_y, refract_dir_z);
			refract_ray->cast(lights, refraction_color, environmentColor, bounces - 1, closestObject, generator, primary_ray, -1, -1, width, depthMap, primary_objID, hit_normal, object_type, hit_list, shadowed, flat, environment, light_samples, indirect_samples, indirect_bounces, octree);
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
					if (objects[closestObject]->object_type == 1) {
						if (cosineNL < 0.000000000001) {
							direct_diffuse_color[0] += objectColor.x * lights[i]->c[0] / 6.0;
							direct_diffuse_color[1] += objectColor.y * lights[i]->c[1] / 6.0;
							direct_diffuse_color[2] += objectColor.z * lights[i]->c[2] / 6.0;
						}
						else if (cosineNL < 0.35) {
							direct_diffuse_color[0] += objectColor.x * lights[i]->c[0] / 2.0;
							direct_diffuse_color[1] += objectColor.y * lights[i]->c[1] / 2.0;
							direct_diffuse_color[2] += objectColor.z * lights[i]->c[2] / 2.0;
						}
						else{
							direct_diffuse_color[0] += objectColor.x * lights[i]->c[0];
							direct_diffuse_color[1] += objectColor.y * lights[i]->c[1];
							direct_diffuse_color[2] += objectColor.z * lights[i]->c[2];
						}
					}
					else {
						direct_diffuse_color[0] += objectColor.x * lights[i]->c[0] * cosineNL;
						direct_diffuse_color[1] += objectColor.y * lights[i]->c[1] * cosineNL;
						direct_diffuse_color[2] += objectColor.z * lights[i]->c[2] * cosineNL;
					}
				}
				
			}
			else if (lights[i]->type == 2) { // bulb light
				vec3 lightCoordinate(lights[i]->x, lights[i]->y, lights[i]->z);
				
				distanceToLight = hitPoint.distance(lightCoordinate);

				vec3 directionToLight(lightCoordinate.x - hitPoint.x, lightCoordinate.y - hitPoint.y, lightCoordinate.z - hitPoint.z);
				directionToLight.normalize();

				light_ray = new ray(hitX + (.0001) * directionToLight.x, hitY + (.0001) * directionToLight.y, hitZ + (.0001) * directionToLight.z, directionToLight.x, directionToLight.y, directionToLight.z);

				std::vector<object*> lightobjects;
				std::vector<std::vector<object*>> lightobjs;
				octree->traverseOctree(light_ray, lightobjs);
				for (int r = 0; r < lightobjs.size(); r++) {
					for (int e = 0; e < lightobjs[r].size(); e++) {
						lightobjects.push_back(lightobjs[r][e]);
					}
				}
				double cosineLN = normal.dot(light_ray->direction);
				if (cosineLN <= 0) {
					distanceToBlock = 0;
				}
				else {
					distanceToBlock = light_ray->castLight(lightobjects, lights[i], distanceToLight);
				}
				delete light_ray;
				vec3 objectColor = objects[closestObject]->getColor();
				if (distanceToBlock < distanceToLight) {
					if (objects[closestObject]->object_type == 1) {
						direct_diffuse_color[0] += objectColor.x * lights[i]->c[0] / 6.0;
						direct_diffuse_color[1] += objectColor.y * lights[i]->c[1] / 6.0;
						direct_diffuse_color[2] += objectColor.z * lights[i]->c[2] / 6.0;
						continue;
					}
					else {
						shadowed = 1;
					}
					//continue;
				}
				if (objects[closestObject]->mat->roughness > 0) {
					normal.x += rand_x;
					normal.y += rand_y;
					normal.z += rand_z;
					normal.normalize();
				}
				double cosineNL = normal.x * directionToLight.x + normal.y * directionToLight.y + normal.z * directionToLight.z;
				if (objects[closestObject]->object_type == 1) {
					if (std::abs(cosineNL * lights[i]->c[0]) < 0.0000000001) {
						direct_diffuse_color[0] += objectColor.x / 6.0;
					}
					else if (std::abs(cosineNL * lights[i]->c[0]) < 0.45) {
						direct_diffuse_color[0] += objectColor.x / 2.0;
					}
					else{
						direct_diffuse_color[0] += objectColor.x;
					}

					if (std::abs(cosineNL * lights[i]->c[1]) < 0.0000000001) {
						direct_diffuse_color[1] += objectColor.y / 6.0;
					}
					else if (std::abs(cosineNL * lights[i]->c[1]) < 0.45) {
						direct_diffuse_color[1] += objectColor.y / 2.0;
					}
					else{
						direct_diffuse_color[1] += objectColor.y;
					}

					if (std::abs(cosineNL * lights[i]->c[2]) < 0.0000000001) {
						direct_diffuse_color[2] += objectColor.z / 6.0;
					}
					else if (std::abs(cosineNL * lights[i]->c[2]) < 0.45) {
						direct_diffuse_color[2] += objectColor.z / 2.0;
					}
					else{
						direct_diffuse_color[2] += objectColor.z;
					}
				}
				else {
					direct_diffuse_color[0] += (objectColor.x * lights[i]->c[0] * cosineNL) / (distanceToLight * distanceToLight);
					direct_diffuse_color[1] += (objectColor.y * lights[i]->c[1] * cosineNL) / (distanceToLight * distanceToLight);
					direct_diffuse_color[2] += (objectColor.z * lights[i]->c[2] * cosineNL) / (distanceToLight * distanceToLight);
				}
				if (objects[closestObject]->mat->eccentricity > 0) {
					vec3 halfway(0,0,0);
					halfway.x = -this->direction.x + directionToLight.x;
					halfway.y = -this->direction.y + directionToLight.y;
					halfway.z = -this->direction.z + directionToLight.z;
					halfway.normalize();
					double base_intensity = normal.dot(halfway);
					double specular_intensity = pow(base_intensity, objects[closestObject]->mat->eccentricity);
					double highlight_r = specular_intensity * lights[i]->c[0];
					double highlight_g = specular_intensity * lights[i]->c[1];
					double highlight_b = specular_intensity * lights[i]->c[2];
					if (highlight_r > 0.1 || highlight_g > 0.1 || highlight_b > 0.1) {
						direct_diffuse_color[0] += objectColor.x + 0.5;
						direct_diffuse_color[1] += objectColor.y + 0.5;
						direct_diffuse_color[2] += objectColor.z + 0.5;
						object_type = previous_obj_type;
					}
				}
				
			}
			else if (lights[i]->type == 3) { // area light
				int num_shadow_rays = 0;
				for (int n = 0; n < light_samples; n++) {
					
					std::uniform_real_distribution<double> distribution_x(-(lights[i]->light_center.x - lights[i]->p0.x), lights[i]->light_center.x - lights[i]->p0.x);
					double rand_x = distribution_x(generator);
					std::uniform_real_distribution<double> distribution_y(-(lights[i]->light_center.y - lights[i]->p0.y), lights[i]->light_center.y - lights[i]->p0.y);
					double rand_y = distribution_y(generator);
					std::uniform_real_distribution<double> distribution_z(-(lights[i]->light_center.z - lights[i]->p0.z), lights[i]->light_center.z - lights[i]->p0.z);
					double rand_z = distribution_z(generator);
					//std::cout << rand_x << " " << rand_y << " " << rand_z << std::endl;
					vec3 lightCoordinate(lights[i]->light_center.x + rand_x, lights[i]->light_center.y + rand_y, lights[i]->light_center.z + rand_z);
					
					distanceToLight = hitPoint.distance(lightCoordinate);

					vec3 directionToLight(lightCoordinate.x - hitPoint.x, lightCoordinate.y - hitPoint.y, lightCoordinate.z - hitPoint.z);
					directionToLight.normalize();

					light_ray = new ray(hitX, hitY, hitZ, directionToLight.x, directionToLight.y, directionToLight.z);
					distanceToBlock = light_ray->castLight(objects, lights[i], distanceToLight);
					delete light_ray;
					if (distanceToBlock < distanceToLight) {
						//direct_diffuse_color[0] += 0;
						//direct_diffuse_color[1] += 0;
						//direct_diffuse_color[2] += 0;
						
						//continue;
						if (objects[closestObject]->object_type == 0) {
							num_shadow_rays++;
						}
						else {
							continue;
						}
						
					}

					vec3 normal = objects[closestObject]->getNormal(hitX, hitY, hitZ, flat);
					normal.normalize();
					if (objects[closestObject]->mat->roughness > 0) {
						normal.x += rand_x;
						normal.y += rand_y;
						normal.z += rand_z;
						normal.normalize();
					}
					double cosineNE = normal.dot(this->direction);
					if (cosineNE > 0) {
						normal.reverse();
					}
					double cosineNL = normal.x * directionToLight.x + normal.y * directionToLight.y + normal.z * directionToLight.z;
					if (cosineNL > 0) {
						vec3 objectColor = objects[closestObject]->getColor();
						direct_diffuse_color[0] += (lights[i]->c[0] * cosineNL) / (distanceToLight * distanceToLight);
						direct_diffuse_color[1] += (lights[i]->c[1] * cosineNL) / (distanceToLight * distanceToLight);
						direct_diffuse_color[2] += (lights[i]->c[2] * cosineNL) / (distanceToLight * distanceToLight);
					}
				}
				direct_diffuse_color[0] /= light_samples;
				direct_diffuse_color[1] /= light_samples;
				direct_diffuse_color[2] /= light_samples;
				if ((double)num_shadow_rays / (double)light_samples > 0.6) {
					shadowed = 1;
				}

				
			}
		}
	}
	if (indirect_bounces > 0) {
		std::vector<int> indirect_hit_list;
		for (int s = 0; s < indirect_samples; s++) {
			double hitX = this->origin.x + distance * this->direction.x;
			double hitY = this->origin.y + distance * this->direction.y;
			double hitZ = this->origin.z + distance * this->direction.z;
			vec3 normal = objects[closestObject]->getNormal(hitX, hitY, hitZ, flat);
			normal.normalize();

			std::uniform_real_distribution<double> distribution(-1, 1);
			vec3 new_dir(distribution(generator), distribution(generator), distribution(generator));
			new_dir.normalize();
			if (normal.dot(new_dir)) {
				new_dir.reverse();
			}
			ray *indirect_diffuse_ray = new ray(hitX + (new_dir.x * .0001), hitY + (new_dir.y * .0001), hitZ + (new_dir.z * .0001), new_dir.x, new_dir.y, new_dir.z);
			double indirect_sample_color[4] {0,0,0,0};
			double cosineNE = normal.dot(this->direction);
			if (cosineNE > 0) {
				normal.reverse();
			}
			int indirect_object_type = 0;
			vec3 indirect_hit_normal(0,0,0);
			int indirect_shadowed = 0;
			int indirect_objID = 0;
			double indirect_depth[1];
			
			indirect_diffuse_ray->cast(lights, indirect_sample_color, environmentColor, bounces, -1, generator, 0, x, y, width, indirect_depth, indirect_objID, indirect_hit_normal, indirect_object_type, indirect_hit_list, indirect_shadowed, flat, environment, light_samples, indirect_samples, indirect_bounces - 1, octree);
			double cosineNL = normal.x * new_dir.x + normal.y * new_dir.y + normal.z * new_dir.z;
			cosineNL = std::abs(cosineNL);
			indirect_diffuse_color[0] += indirect_sample_color[0] * cosineNL * 2 * PI;
			indirect_diffuse_color[1] += indirect_sample_color[1] * cosineNL * 2 * PI;
			indirect_diffuse_color[2] += indirect_sample_color[2] * cosineNL * 2 * PI;
			delete indirect_diffuse_ray;
		}
		indirect_diffuse_color[0] /= indirect_samples;
		indirect_diffuse_color[1] /= indirect_samples;
		indirect_diffuse_color[2] /= indirect_samples;
	}
	vec3 objectColor = objects[closestObject]->getColor();
	double r_channel = direct_diffuse_color[0] + indirect_diffuse_color[0];
	double g_channel = direct_diffuse_color[1] + indirect_diffuse_color[1];
	double b_channel = direct_diffuse_color[2] + indirect_diffuse_color[2];	

	if (objects[closestObject]->object_type == 1) { 

		// 0.05, 0.7 for area lights
		if (r_channel < 0.1) {
			diffuse_color[0] = objectEmission.x + (0.35) * objects[closestObject]->mat->diffuse.x * objectColor.x;
		}
		else if (r_channel < 0.9) {
			diffuse_color[0] = objectEmission.x + (0.7) * objects[closestObject]->mat->diffuse.x * objectColor.x;
		}
		else {
			diffuse_color[0] = objectEmission.x + objects[closestObject]->mat->diffuse.x * objectColor.x;
		}

		if (g_channel < 0.1) {
			diffuse_color[1] = objectEmission.y + (0.35) * objects[closestObject]->mat->diffuse.y * objectColor.y;
		}
		else if (g_channel < 0.9) {
			diffuse_color[1] = objectEmission.y + (0.7) * objects[closestObject]->mat->diffuse.y * objectColor.y;
		}
		else {
			diffuse_color[1] = objectEmission.y + objects[closestObject]->mat->diffuse.y * objectColor.y;
		}

		if (b_channel < 0.1) {
			diffuse_color[2] = objectEmission.z + (0.35) * objects[closestObject]->mat->diffuse.z * objectColor.z;
		}
		else if (b_channel < 0.9) {
			diffuse_color[2] = objectEmission.z + (0.7) * objects[closestObject]->mat->diffuse.z * objectColor.z;
		}
		else {
			diffuse_color[2] = objectEmission.z + objects[closestObject]->mat->diffuse.z * objectColor.z;
		}
	}
	else { 
		diffuse_color[0] = objectEmission.x + (direct_diffuse_color[0] + indirect_diffuse_color[0]) * objects[closestObject]->mat->diffuse.x * objectColor.x;
		diffuse_color[1] = objectEmission.y + (direct_diffuse_color[1] + indirect_diffuse_color[1]) * objects[closestObject]->mat->diffuse.y * objectColor.y;
		diffuse_color[2] = objectEmission.z + (direct_diffuse_color[2] + indirect_diffuse_color[2]) * objects[closestObject]->mat->diffuse.z * objectColor.z;
	}

	color[0] += diffuse_color[0] + reflection_color[0] * objects[closestObject]->mat->shininess.x + refraction_color[0] * objects[closestObject]->mat->transparency.x;
	color[1] += diffuse_color[1] + reflection_color[1] * objects[closestObject]->mat->shininess.y + refraction_color[1] * objects[closestObject]->mat->transparency.y;
	color[2] += diffuse_color[2] + reflection_color[2] * objects[closestObject]->mat->shininess.z + refraction_color[2] * objects[closestObject]->mat->transparency.z;
	color[3] += 1;
	return true;
}

double ray::castLight(std::vector<object*> &objects, light *target_light, double distance) {
	for (int i = 0; i < objects.size(); i++) {
		objects[i]->shadowHit(this, target_light, distance);
	}
	return distance;
}

double ray::detect_edge(std::vector<light*> &lights, double color[4], int bounces, int lastObject, std::default_random_engine &generator, int &stencil_objID, int &object_type, std::vector<int> &hit_list, OctreeNode *octree) {
	double distance = std::numeric_limits<double>::max();
	int closestObject = -1;
	double diffuse_color[3] = {0,0,0};
	double reflection_color[4] = {0,0,0,0};
	double refraction_color[4] = {0,0,0,0};
	std::vector<object*> objects;
	std::vector<std::vector<object*>> objs;
	octree->traverseOctree(this, objs);
	for (int r = 0; r < objs.size(); r++) {
		for (int e = 0; e < objs[r].size(); e++) {
			objects.push_back(objs[r][e]);
		}
	}
	for (int j = 0; j < objects.size(); j++) {
		if (objects[j]->hit(this, objects, lights, color, distance)) {
			closestObject = j;	
		}
	}
	if (closestObject < 0) {
		stencil_objID = -1;
		return distance;
	}
	stencil_objID = objects[closestObject]->objectID;
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
		vec3 normal = objects[closestObject]->getNormal(hitX, hitY, hitZ, 0);
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
		distance = reflect_ray->detect_edge(lights, color, bounces, -1, generator, stencil_objID, object_type, hit_list, octree);
		delete reflect_ray;
	}
	return distance;
}

ray::~ray() {

}
