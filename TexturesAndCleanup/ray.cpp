// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// ray.cpp

#include "ray.h"

#define PI 3.14159265

// bool ray::cast(Scene *theScene, ImageBuffers *imageBuffers, double color[4], int bounces, int lastObject, int primary_ray, int x, int y, int &primary_objID, vec3 &hit_normal, int &object_type, std::vector<int> &hit_list, int &shadowed, int indirect_bounces) {
bool ray::cast(Scene *theScene, ImageBuffers *imageBuffers, HitRecord *hitRecord, double color[4], int bounces, int indirect_bounces) {
	double distance = std::numeric_limits<double>::max();
	int closestObject = -1;
	double direct_diffuse_color[3] = {0,0,0};
	double indirect_diffuse_color[4] = {0,0,0,0};
	double diffuse_color[3] = {0,0,0};
	double reflection_color[4] = {0,0,0,0};
	double refraction_color[4] = {0,0,0,0};
	std::vector<object*> objects;
	std::vector<std::vector<object*>> objs;
	theScene->octree->traverseOctree(this, objs);

	for (int r = 0; r < objs.size(); r++) {
		for (int e = 0; e < objs[r].size(); e++) {
			objects.push_back(objs[r][e]);
		}
	}
	int closestObjectID = -1;
	for (int j = 0; j < objects.size(); j++) {
		if (objects[j]->hit(this, objects, theScene->lights, color, distance)) { 
			closestObject = j;
			closestObjectID = objects[j]->objectID;	
		}
	}
	if (closestObject == -1) {
		// ray didn't hit any objects
		if (hitRecord->x != -1 && hitRecord->y != -1 && hitRecord->primary_ray) {
			imageBuffers->environmentMap[hitRecord->y * theScene->width + hitRecord->x] = 1;
		}
		if (theScene->environment) {
			if (theScene->environmentColor[3] == 1) { // gradient - color to white
				double t = 0.5 * (this->direction.y + 1.0);
				color[0] = (1.0 - t) * 1.0 + t * theScene->environmentColor[0];
				color[1] = (1.0 - t) * 1.0 + t * theScene->environmentColor[1];
				color[2] = (1.0 - t) * 1.0 + t * theScene->environmentColor[2];
				color[3] = 1;
			}
			else if (theScene->environmentColor[3] == 0) { // solid color
				color[0] = theScene->environmentColor[0];
				color[1] = theScene->environmentColor[1];
				color[2] = theScene->environmentColor[2];
				color[3] = 1;
			}	
		}
		else {
			color[3] = 0;
		}	
		return false;
	}
	
	hitRecord->primary_objID = objects[closestObject]->objectID;
	hitRecord->hit_list.push_back(objects[closestObject]->objectID);
	if (hitRecord->x != -1 && hitRecord->y != -1 && hitRecord->primary_ray) {
		imageBuffers->depthMap[hitRecord->y * theScene->width + hitRecord->x] = distance;
	}
	int previous_obj_type = hitRecord->object_type;
	hitRecord->object_type = objects[closestObject]->object_type;
  	std::normal_distribution<double> distribution(0, objects[closestObject]->mat->roughness);
	double rand_x = distribution(theScene->generator);
	double rand_y = distribution(theScene->generator);
	double rand_z = distribution(theScene->generator);

	double hitX = this->origin.x + distance * this->direction.x;
	double hitY = this->origin.y + distance * this->direction.y;
	double hitZ = this->origin.z + distance * this->direction.z;
	vec3 hitPoint(hitX, hitY, hitZ);
	vec3 normal = objects[closestObject]->getNormal(hitX, hitY, hitZ, theScene->flat);
	normal.normalize();
	double cosineNE = normal.dot(this->direction);
	if (cosineNE > 0) {
		normal.reverse();
	}
	hitRecord->hit_normal.x = normal.x;
	hitRecord->hit_normal.y = normal.y;
	hitRecord->hit_normal.z = normal.z;

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
		ray *reflect_ray = new ray(hitX + (.0001) * reflect_dir_x, hitY + (.0001) * reflect_dir_y, hitZ + (.0001) * reflect_dir_z, reflect_dir_x, reflect_dir_y, reflect_dir_z);
		hitRecord->lastObject = -1;
		reflect_ray->cast(theScene, imageBuffers, hitRecord, reflection_color, bounces - 1, indirect_bounces);
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
			normal = objects[closestObject]->getNormal(hitX, hitY, hitZ, theScene->flat);
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
			hitRecord->lastObject = -1;
			hitRecord->x = -1;
			hitRecord->y = -1;
			refract_ray->cast(theScene, imageBuffers, hitRecord, refraction_color, bounces - 1, indirect_bounces);
			delete refract_ray;
		}
		else {
			refract_dir_x = ior * this->direction.x + (ior * n_dot_i - sqrt(k)) * normal.x;
			refract_dir_y = ior * this->direction.y + (ior * n_dot_i - sqrt(k)) * normal.y;
			refract_dir_z = ior * this->direction.z + (ior * n_dot_i - sqrt(k)) * normal.z;
			ray *refract_ray = new ray(hitX + (.0001) * refract_dir_x, hitY + (.0001) * refract_dir_y, hitZ + (.0001) * refract_dir_z, refract_dir_x, refract_dir_y, refract_dir_z);
			hitRecord->lastObject = closestObject;
			hitRecord->x = -1;
			hitRecord->y = -1;
			refract_ray->cast(theScene, imageBuffers, hitRecord, refraction_color, bounces - 1, indirect_bounces);
			delete refract_ray;
		}
	}
	if (hitRecord->lastObject == closestObject) {
		// don't count diffuse color for refraction exit ray of object
	}
	else if (theScene->lights.size() > 0) {
		ray *light_ray;
		double distanceToLight;
		double distanceToBlock;
		for (int i = 0; i < theScene->lights.size(); i++) {
			if (theScene->lights[i]->type == 1) { // sun light

				distanceToLight = std::numeric_limits<double>::max();

				light_ray = new ray(hitX, hitY, hitZ, theScene->lights[i]->x, theScene->lights[i]->y, theScene->lights[i]->z);
				distanceToBlock = light_ray->castLight(objects, theScene->lights[i], distanceToLight);
				delete light_ray;
				if (distanceToBlock < distanceToLight) {
					if (objects[closestObject]->object_type == 1) {
						diffuse_color[0] += 0;
						diffuse_color[1] += 0;
						diffuse_color[2] += 0;
						continue;
					}
					else {
						hitRecord->shadowed = 1;
					}
				}
				if (objects[closestObject]->mat->roughness > 0) {
					normal.x += rand_x;
					normal.y += rand_y;
					normal.z += rand_z;
					normal.normalize();
				}
				// hitRecord->setAttributes(distance, hitPoint, normal);
				double cosineNL = normal.x * theScene->lights[i]->x + normal.y * theScene->lights[i]->y + normal.z * theScene->lights[i]->z;
				if (cosineNL > 0) {
					vec3 texCoords = objects[closestObject]->getTextureCoordinates(hitPoint);
					vec3 objectColor = objects[closestObject]->getColor(texCoords.x, texCoords.y, hitPoint);
					if (objects[closestObject]->object_type == 1) {
						if (cosineNL < 0.000000000001) {
							direct_diffuse_color[0] += objectColor.x * theScene->lights[i]->c[0] / 6.0;
							direct_diffuse_color[1] += objectColor.y * theScene->lights[i]->c[1] / 6.0;
							direct_diffuse_color[2] += objectColor.z * theScene->lights[i]->c[2] / 6.0;
						}
						else if (cosineNL < 0.35) {
							direct_diffuse_color[0] += objectColor.x * theScene->lights[i]->c[0] / 2.0;
							direct_diffuse_color[1] += objectColor.y * theScene->lights[i]->c[1] / 2.0;
							direct_diffuse_color[2] += objectColor.z * theScene->lights[i]->c[2] / 2.0;
						}
						else{
							direct_diffuse_color[0] += objectColor.x * theScene->lights[i]->c[0];
							direct_diffuse_color[1] += objectColor.y * theScene->lights[i]->c[1];
							direct_diffuse_color[2] += objectColor.z * theScene->lights[i]->c[2];
						}
					}
					else {
						direct_diffuse_color[0] += objectColor.x * theScene->lights[i]->c[0] * cosineNL;
						direct_diffuse_color[1] += objectColor.y * theScene->lights[i]->c[1] * cosineNL;
						direct_diffuse_color[2] += objectColor.z * theScene->lights[i]->c[2] * cosineNL;
					}
				}
				
			}
			else if (theScene->lights[i]->type == 2) { // bulb light
				vec3 lightCoordinate(theScene->lights[i]->x, theScene->lights[i]->y, theScene->lights[i]->z);
				
				distanceToLight = hitPoint.distance(lightCoordinate);

				vec3 directionToLight(lightCoordinate.x - hitPoint.x, lightCoordinate.y - hitPoint.y, lightCoordinate.z - hitPoint.z);
				directionToLight.normalize();

				light_ray = new ray(hitX + (.0001) * directionToLight.x, hitY + (.0001) * directionToLight.y, hitZ + (.0001) * directionToLight.z, directionToLight.x, directionToLight.y, directionToLight.z);

				std::vector<object*> lightobjects;
				std::vector<std::vector<object*>> lightobjs;
				theScene->octree->traverseOctree(light_ray, lightobjs);
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
					distanceToBlock = light_ray->castLight(lightobjects, theScene->lights[i], distanceToLight);
				}
				delete light_ray;
				vec3 texCoords = objects[closestObject]->getTextureCoordinates(hitPoint);
				vec3 objectColor = objects[closestObject]->getColor(texCoords.x, texCoords.y, hitPoint);
				if (distanceToBlock < distanceToLight) {
					if (objects[closestObject]->object_type == 1) {
						direct_diffuse_color[0] += objectColor.x * theScene->lights[i]->c[0] / 6.0;
						direct_diffuse_color[1] += objectColor.y * theScene->lights[i]->c[1] / 6.0;
						direct_diffuse_color[2] += objectColor.z * theScene->lights[i]->c[2] / 6.0;
						continue;
					}
					else {
						hitRecord->shadowed = 1;
					}
					//continue;
				}
				if (objects[closestObject]->mat->roughness > 0) {
					normal.x += rand_x;
					normal.y += rand_y;
					normal.z += rand_z;
					normal.normalize();
				}
				// hitRecord->setAttributes(distance, hitPoint, normal);
				double cosineNL = normal.x * directionToLight.x + normal.y * directionToLight.y + normal.z * directionToLight.z;
				if (objects[closestObject]->object_type == 1) {
					if (std::abs(cosineNL * theScene->lights[i]->c[0]) < 0.0000000001) {
						direct_diffuse_color[0] += objectColor.x / 6.0;
					}
					else if (std::abs(cosineNL * theScene->lights[i]->c[0]) < 0.45) {
						direct_diffuse_color[0] += objectColor.x / 2.0;
					}
					else{
						direct_diffuse_color[0] += objectColor.x;
					}

					if (std::abs(cosineNL * theScene->lights[i]->c[1]) < 0.0000000001) {
						direct_diffuse_color[1] += objectColor.y / 6.0;
					}
					else if (std::abs(cosineNL * theScene->lights[i]->c[1]) < 0.45) {
						direct_diffuse_color[1] += objectColor.y / 2.0;
					}
					else{
						direct_diffuse_color[1] += objectColor.y;
					}

					if (std::abs(cosineNL * theScene->lights[i]->c[2]) < 0.0000000001) {
						direct_diffuse_color[2] += objectColor.z / 6.0;
					}
					else if (std::abs(cosineNL * theScene->lights[i]->c[2]) < 0.45) {
						direct_diffuse_color[2] += objectColor.z / 2.0;
					}
					else{
						direct_diffuse_color[2] += objectColor.z;
					}
				}
				else {
					direct_diffuse_color[0] += (objectColor.x * theScene->lights[i]->c[0] * cosineNL) / (distanceToLight * distanceToLight);
					direct_diffuse_color[1] += (objectColor.y * theScene->lights[i]->c[1] * cosineNL) / (distanceToLight * distanceToLight);
					direct_diffuse_color[2] += (objectColor.z * theScene->lights[i]->c[2] * cosineNL) / (distanceToLight * distanceToLight);
				}
				if (objects[closestObject]->mat->eccentricity > 0) {
					vec3 halfway(0,0,0);
					halfway.x = -this->direction.x + directionToLight.x;
					halfway.y = -this->direction.y + directionToLight.y;
					halfway.z = -this->direction.z + directionToLight.z;
					halfway.normalize();
					double base_intensity = normal.dot(halfway);
					double specular_intensity = pow(base_intensity, objects[closestObject]->mat->eccentricity);
					double highlight_r = specular_intensity * theScene->lights[i]->c[0];
					double highlight_g = specular_intensity * theScene->lights[i]->c[1];
					double highlight_b = specular_intensity * theScene->lights[i]->c[2];
					if (highlight_r > 0.1 || highlight_g > 0.1 || highlight_b > 0.1) {
						direct_diffuse_color[0] += objectColor.x + 0.5;
						direct_diffuse_color[1] += objectColor.y + 0.5;
						direct_diffuse_color[2] += objectColor.z + 0.5;
						hitRecord->object_type = previous_obj_type;
					}
				}
				
			}
			else if (theScene->lights[i]->type == 3) { // area light
				int num_shadow_rays = 0;
				for (int n = 0; n < theScene->light_samples; n++) {
					
					std::uniform_real_distribution<double> distribution_x(-(theScene->lights[i]->light_center.x - theScene->lights[i]->p0.x), theScene->lights[i]->light_center.x - theScene->lights[i]->p0.x);
					double rand_x = distribution_x(theScene->generator);
					std::uniform_real_distribution<double> distribution_y(-(theScene->lights[i]->light_center.y - theScene->lights[i]->p0.y), theScene->lights[i]->light_center.y - theScene->lights[i]->p0.y);
					double rand_y = distribution_y(theScene->generator);
					std::uniform_real_distribution<double> distribution_z(-(theScene->lights[i]->light_center.z - theScene->lights[i]->p0.z), theScene->lights[i]->light_center.z - theScene->lights[i]->p0.z);
					double rand_z = distribution_z(theScene->generator);
					//std::cout << rand_x << " " << rand_y << " " << rand_z << std::endl;
					vec3 lightCoordinate(theScene->lights[i]->light_center.x + rand_x, theScene->lights[i]->light_center.y + rand_y, theScene->lights[i]->light_center.z + rand_z);
					
					distanceToLight = hitPoint.distance(lightCoordinate);

					vec3 directionToLight(lightCoordinate.x - hitPoint.x, lightCoordinate.y - hitPoint.y, lightCoordinate.z - hitPoint.z);
					directionToLight.normalize();

					light_ray = new ray(hitX, hitY, hitZ, directionToLight.x, directionToLight.y, directionToLight.z);
					distanceToBlock = light_ray->castLight(objects, theScene->lights[i], distanceToLight);
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

					vec3 normal = objects[closestObject]->getNormal(hitX, hitY, hitZ, theScene->flat);
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
					// hitRecord->setAttributes(distance, hitPoint, normal);
					double cosineNL = normal.x * directionToLight.x + normal.y * directionToLight.y + normal.z * directionToLight.z;
					if (cosineNL > 0) {
						vec3 texCoords = objects[closestObject]->getTextureCoordinates(hitPoint);
						vec3 objectColor = objects[closestObject]->getColor(texCoords.x, texCoords.y, hitPoint);
						direct_diffuse_color[0] += (theScene->lights[i]->c[0] * cosineNL) / (distanceToLight * distanceToLight);
						direct_diffuse_color[1] += (theScene->lights[i]->c[1] * cosineNL) / (distanceToLight * distanceToLight);
						direct_diffuse_color[2] += (theScene->lights[i]->c[2] * cosineNL) / (distanceToLight * distanceToLight);
					}
				}
				direct_diffuse_color[0] /= theScene->light_samples;
				direct_diffuse_color[1] /= theScene->light_samples;
				direct_diffuse_color[2] /= theScene->light_samples;
				if ((double)num_shadow_rays / (double)theScene->light_samples > 0.6) {
					hitRecord->shadowed = 1;
				}

				
			}
		}
	}
	if (indirect_bounces > 0) {
		std::vector<int> indirect_hit_list;
		for (int s = 0; s < theScene->indirect_samples; s++) {
			double hitX = this->origin.x + hitRecord->distance * this->direction.x;
			double hitY = this->origin.y + hitRecord->distance * this->direction.y;
			double hitZ = this->origin.z + hitRecord->distance * this->direction.z;
			vec3 normal = objects[closestObject]->getNormal(hitX, hitY, hitZ, theScene->flat);
			normal.normalize();

			std::uniform_real_distribution<double> distribution(-1, 1);
			vec3 new_dir(distribution(theScene->generator), distribution(theScene->generator), distribution(theScene->generator));
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
			// int indirect_object_type = 0;
			// vec3 indirect_hit_normal(0,0,0);
			// int indirect_shadowed = 0;
			// int indirect_objID = 0;
			// double indirect_depth[1];
			HitRecord *indirectHitRecord = new HitRecord(hitRecord->x, hitRecord->y);
			indirectHitRecord->primary_ray = 0;
			indirectHitRecord->primary_objID = 0;
			// indirect_diffuse_ray->cast(theScene, imageBuffers, indirect_sample_color, bounces, -1, 0, x, y, indirect_objID, indirect_hit_normal, indirect_object_type, indirect_hit_list, indirect_shadowed, indirect_bounces - 1);
			indirect_diffuse_ray->cast(theScene, imageBuffers, indirectHitRecord, indirect_sample_color, bounces, indirect_bounces - 1);
			double cosineNL = normal.x * new_dir.x + normal.y * new_dir.y + normal.z * new_dir.z;
			cosineNL = std::abs(cosineNL);
			indirect_diffuse_color[0] += indirect_sample_color[0] * cosineNL * 2 * PI;
			indirect_diffuse_color[1] += indirect_sample_color[1] * cosineNL * 2 * PI;
			indirect_diffuse_color[2] += indirect_sample_color[2] * cosineNL * 2 * PI;
			delete indirect_diffuse_ray;
			delete indirectHitRecord;
		}
		indirect_diffuse_color[0] /= theScene->indirect_samples;
		indirect_diffuse_color[1] /= theScene->indirect_samples;
		indirect_diffuse_color[2] /= theScene->indirect_samples;
	}
	vec3 texCoords = objects[closestObject]->getTextureCoordinates(hitPoint);
	vec3 objectColor = objects[closestObject]->getColor(texCoords.x, texCoords.y, hitPoint);
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

double ray::detect_edge(Scene *theScene, double color[4], int bounces, int lastObject, int &stencil_objID, int &object_type, std::vector<int> &hit_list) {
	double distance = std::numeric_limits<double>::max();
	int closestObject = -1;
	double diffuse_color[3] = {0,0,0};
	double reflection_color[4] = {0,0,0,0};
	double refraction_color[4] = {0,0,0,0};
	std::vector<object*> objects;
	std::vector<std::vector<object*>> objs;
	theScene->octree->traverseOctree(this, objs);
	for (int r = 0; r < objs.size(); r++) {
		for (int e = 0; e < objs[r].size(); e++) {
			objects.push_back(objs[r][e]);
		}
	}
	for (int j = 0; j < objects.size(); j++) {
		if (objects[j]->hit(this, objects, theScene->lights, color, distance)) {
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
	double rand_x = distribution(theScene->generator);
	double rand_y = distribution(theScene->generator);
	double rand_z = distribution(theScene->generator);
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
		ray *reflect_ray = new ray(hitX + (.0001) * reflect_dir_x, hitY + (.0001) * reflect_dir_y, hitZ + (.0001) * reflect_dir_z, reflect_dir_x, reflect_dir_y, reflect_dir_z);
		vec3 hitnormal(0,0,0);
		distance = reflect_ray->detect_edge(theScene, color, bounces, -1, stencil_objID, object_type, hit_list);
		delete reflect_ray;
	}
	return distance;
}

ray::~ray() {

}
