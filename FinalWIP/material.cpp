// Nick Moon
// nm9nz
// CS4810 HW3: Raytracer
// material.cpp

#include "material.h"

material::material() {
	vec3 diffuse_color(0,0,0);
	this->diffuse_color = diffuse_color;
	vec3 ambient(0,0,0);
	this->ambient = ambient;
	vec3 shininess(0,0,0);
	this->shininess = shininess;
	vec3 transparency(0,0,0);
	this->transparency = transparency;
	vec3 diffuse(1,1,1);
	this->diffuse = diffuse;
	this->ior = 0.0;
	this->roughness = 0.0;
	this->eccentricity = 0.0;
	this->diffuse_map = "";
	this->object_type = 1;
}

material::material(vec3 &shininess, vec3 &transparency, double ior, double roughness, double eccentricity) {
	vec3 diffuse_color(0,0,0);
	this->diffuse_color = diffuse_color;
	vec3 ambient(0,0,0);
	this->ambient = ambient;
	this->shininess = shininess;
	this->transparency = transparency;
	vec3 diffuse(0,0,0);
	this->diffuse = diffuse;
	this->diffuse.x = (1 - this->shininess.x) * (1 - transparency.x);
	this->diffuse.y = (1 - this->shininess.y) * (1 - transparency.y);
	this->diffuse.z = (1 - this->shininess.z) * (1 - transparency.z);
	this->transparency.x = (1 - this->shininess.x) * this->transparency.x;
	this->transparency.y = (1 - this->shininess.y) * this->transparency.y;
	this->transparency.z = (1 - this->shininess.z) * this->transparency.z;
	this->ior = ior;
	this->roughness = roughness;
	this->eccentricity = eccentricity;
	this->diffuse_map = "";
	this->object_type = 1;
}

material::material(vec3 &diffuse_color, vec3 &ambient, vec3 &diffuse, vec3 &shininess, vec3 &transparency, double ior, double roughness, double eccentricity, std::string diffuse_map) {
	this->diffuse_color = diffuse_color;
	this->ambient = ambient; 
	this->shininess = shininess;
	this->transparency = transparency;
	this->diffuse = diffuse;
	this->diffuse.x = (1 - this->shininess.x) * (1 - transparency.x);
	this->diffuse.y = (1 - this->shininess.y) * (1 - transparency.y);
	this->diffuse.z = (1 - this->shininess.z) * (1 - transparency.z);
	this->transparency.x = (1 - this->shininess.x) * this->transparency.x;
	this->transparency.y = (1 - this->shininess.y) * this->transparency.y;
	this->transparency.z = (1 - this->shininess.z) * this->transparency.z;
	this->ior = ior;
	this->roughness = roughness;
	this->eccentricity = eccentricity;
	this->diffuse_map = diffuse_map;
	this->object_type = 1;
}

material::~material() {

}
