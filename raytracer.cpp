// Nick Moon
// nm9nz
// CS4810 Raytracer

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>
#include <ctime>

#define cimg_use_png
#include "CImg.h"

#include "vec3.h"
#include "vertex.h"
#include "ray.h"
#include "light.h"
#include "material.h"
#include "object.h"
#include "sphere.h"
#include "plane.h"
#include "tri.h"
#include "particle.h"
#include "painter.h"
#include "style.h"


#define PI 3.14159265

double printProgress(int xi, int yi, int width, int height, double milestone, time_t start) {
	double progress = (double)yi / (double)height;
	if (progress - milestone >= 0.1) {
		std::cout << "Progress: " << (progress * 100) << "% done, time elapsed: " << difftime(time(NULL), start) << std::endl;
		return progress;
	}
	return milestone;
	
}



double clip(double color) {
	if (color < 0) {
		return 0;
	}
	else if (color > 255) {
		return 255;
	}
	else {
		return color;
	}
}

double maxItoD(int x, int y) {
	if (x >= y) {
		return x;
	}
	else {
		return y;
	}
}

vec3 cameraDirection(int camera_mode, double xi, double yi, double width, double height, vec3 &forward, vec3 &right, vec3 &up, int &shoot) {
	vec3 camera_direction(0,0,0);
	double sx = (2 * xi - width) / maxItoD(width,height);
	double sy = (height - 2 * yi) / maxItoD(width,height);
	if (camera_mode == 0) {
		camera_direction.x = forward.x + sx * right.x + sy * up.x;
		camera_direction.y = forward.y + sx * right.y + sy * up.y;
		camera_direction.z = forward.z + sx * right.z + sy * up.z;
	}
	else if (camera_mode == 1) {
		sx /= forward.magnitude();
		sy /= forward.magnitude();
		vec3 normalized_f = forward;
		normalized_f.normalize();
		double r_squared = sx * sx + sy * sy;
		if (r_squared > 1) {
			shoot = 0;
			return camera_direction;
		}
		double sz = sqrt(1 - r_squared);
		camera_direction.x = sz * normalized_f.x + sx * right.x + sy * up.x;
		camera_direction.y = sz * normalized_f.y + sx * right.y + sy * up.y;
		camera_direction.z = sz * normalized_f.z + sx * right.z + sy * up.z;
	}
	else if (camera_mode == 2) {
		double theta_x = xi * (360 / width);
		double theta_y = yi * (180 / height);
		
		double sx = cos(theta_x * PI/180);
		double sy = sin(theta_y * PI/180);
		
		camera_direction.x = sx;
		camera_direction.y = sy;
		camera_direction.z = forward.z;
	}
	camera_direction.normalize();
	return camera_direction;
}



int main(int argc, char** argv) {
	time_t start = time(0);
	if (argc != 2) {
		std::cout << "You must provide input file name as command line argument!" << std::endl;
		exit(1);
	}
	int width, height = 0;
	double milestone = 0.0;
	double x, y, z, r = 0;
	double A, B, C, D = 0;
	int objectID = 0;
	int object_type = 0;
	double ior = 1.458;
	int bounces = 4;
	double roughness = 0;
	int v1, v2, v3, v1n, v2n, v3n;
	vertex *vert1, *vert2, *vert3;
	double lastNormal[3] = {0,0,0};
	int spp = 1;
	int camera_mode = 0;
	double lastColor[3] = {1, 1, 1};
	double resultColor[4] = {0, 0, 0, 0};
	double sx, sy = 0;
	vec3 eye(0,0,0);
	vec3 forward(0,0,-1);
	vec3 right(1,0,0);
	vec3 up(0,1,0);
	vec3 direction(0,0,0);
	vec3 shininess(0,0,0);
	vec3 transparency(0,0,0);
	vec3 *vert1n;
	vec3 *vert2n;
	vec3 *vert3n;
	vec3 hit_normal(0,0,0);
	int stencil_ring_samples = 0;
	double stencil_rings = 0;
	double outline_cutoff = 1;
	int shoot = 1;
	int output_depth = 0;
	double e_slope = 1;
	int post_process_outline = 0;
	double stencil_radius = 0;
	std::string v1parts = "";
	std::string v2parts = "";
	std::string v3parts = "";
	std::string lineBuffer = "";
	std::string token = "";
	std::string outputFileName = "";
	std::ifstream inputFile;
	inputFile.open(argv[1]);

	std::getline(inputFile, lineBuffer);
	std::stringstream ss(lineBuffer);
	ss >> token;
	if (token.compare("png")) {
		std::cout << "You must specifiy 'png' at the beginning of the file!" << std::endl;
		inputFile.close();
		exit(1);	
	}
	ss >> width;
	ss >> height;
	ss >> outputFileName;
	
	cimg_library::CImg<float> myImage(width, height, 1, 4, 0);

    	unsigned char zero = 0;
    	myImage = zero;

	std::vector<object*> objects;
	std::vector<light*> lights;
	std::vector<vertex*> verts;
	std::vector<vec3*> normals;
	std::vector<std::shared_ptr<PaintParticle>> paint_particles;
        ray *primary_ray;
	
	while (std::getline(inputFile, lineBuffer)) {
		if (lineBuffer.empty()) {
			continue;
		}
		std::stringstream ss(lineBuffer);
		ss >> token;
		if (!token.compare("sphere")) {
			ss >> x;
			ss >> y;
			ss >> z;
			ss >> r;
			material *m = new material(shininess, transparency, ior, roughness);
			object *s = new sphere(x, y, z, r, lastColor, m, objectID, object_type);
			objects.push_back(s);
		}
		else if (!token.compare("sun")) {
			ss >> x;
			ss >> y;
			ss >> z;
			light *l = new light(1, x, y, z, lastColor);
			lights.push_back(l);
		}
		else if (!token.compare("color")) {
			ss >> lastColor[0];
			ss >> lastColor[1];
			ss >> lastColor[2];
		}
		else if (!token.compare("plane")) {
			ss >> A;
			ss >> B;
			ss >> C;
			ss >> D;
			material *m = new material(shininess, transparency, ior, roughness);
			object *p = new plane(A, B, C, D, lastColor, m, objectID);
			objects.push_back(p);
		}
		else if (!token.compare("xyz")) {
			ss >> x;
			ss >> y;
			ss >> z;
			vec3 coordinate(x,y,z);
			vec3 normal(lastNormal[0], lastNormal[1], lastNormal[2]);
			normal.normalize();
			vec3 color(lastColor[0], lastColor[1], lastColor[2]);
			vertex *v = new vertex(coordinate, normal, color);
			verts.push_back(v);
		}
		else if (!token.compare("trif")) {
			ss >> v1;
			ss >> v2;
			ss >> v3;
			if (v1 < 0) {
				vert1 = verts.at(verts.size() + v1);
			}
			else {
				vert1 = verts.at(v1 - 1);
			}
			if (v2 < 0) {
				vert2 = verts.at(verts.size() + v2);
			}
			else {
				vert2 = verts.at(v2 - 1);
			}
			if (v3 < 0) {
				vert3 = verts.at(verts.size() + v3);
			}
			else {
				vert3 = verts.at(v3 - 1);
			}
			material *m = new material(shininess, transparency, ior, roughness);
			tri *t = new tri(vert1, vert2, vert3, lastColor, m, objectID);
			objects.push_back(t);
		}
		else if (!token.compare("normal")) {
			ss >> lastNormal[0];
			ss >> lastNormal[1];
			ss >> lastNormal[2];
		}
		else if (!token.compare("bulb")) {
			ss >> x;
			ss >> y;
			ss >> z;
			light *l = new light(2, x, y, z, lastColor);
			lights.push_back(l);
		}
		else if (!token.compare("shininess")) {
			ss >> shininess.x;
			if (ss >> shininess.y) {
				ss >> shininess.z;
			}
			else {
				shininess.y = shininess.x;
				shininess.z = shininess.x;
			}
		}
		else if (!token.compare("transparency")) {
			ss >> transparency.x;
			if (ss >> transparency.y) {
				ss >> transparency.z;
			}
			else {
				transparency.y = transparency.x;
				transparency.z = transparency.x;
			}
		}
		else if (!token.compare("ior")) {
			ss >> ior;
		}
		else if (!token.compare("bounces")) {
			ss >> bounces;
		}
		else if (!token.compare("roughness")) {
			ss >> roughness;
		}
		else if (!token.compare("eye")) {
			ss >> eye.x;
			ss >> eye.y;
			ss >> eye.z;
		}
		else if (!token.compare("forward")) {
			ss >> forward.x;
			ss >> forward.y;
			ss >> forward.z;

			vec3 p = forward.cross(up);
			up = p.cross(forward);
			up.normalize();

			right = forward.cross(up);
			right.normalize();
		}
		else if (!token.compare("up")) {
			ss >> up.x;
			ss >> up.y;
			ss >> up.z;

			vec3 p = forward.cross(up);
			up = p.cross(forward);
			up.normalize();

			right = forward.cross(up);
			right.normalize();
		}
		else if (!token.compare("aa")) {
			ss >> spp;
		}
		else if (!token.compare("fisheye")) {
			camera_mode = 1;
		}
		else if (!token.compare("panorama")) {
			camera_mode = 2;
		}
		else if (!token.compare("dof")) {
			camera_mode = 2;
		}
		else if (!token.compare("v")) {
			ss >> x;
			ss >> y;
			ss >> z;
			vec3 coordinate(x,y,z);
			vec3 normal(lastNormal[0], lastNormal[1], lastNormal[2]);
			normal.normalize();
			vec3 color(lastColor[0], lastColor[1], lastColor[2]);
			vertex *v = new vertex(coordinate, normal, color);
			verts.push_back(v);
		}
		else if (!token.compare("vn")) {
			ss >> x;
			ss >> y;
			ss >> z;
			vec3 *v = new vec3(x, y, z);
			normals.push_back(v);
		}
		else if (!token.compare("f")) {
			ss >> v1parts;
			ss >> v2parts;
			ss >> v3parts;

			int v1Slash = v1parts.find("/", 0);
			int v1tStart = v1Slash + 1;
			int v1tSlash = v1parts.find("/", v1tStart);
			int v1nStart = v1tSlash + 1;
			v1 = std::stoi(v1parts.substr(0, v1Slash), nullptr);
			//v1t = std::stoi(v1parts.substr(v1tStart, v1tSlash - v1tStart), nullptr, 16);
			v1n = std::stoi(v1parts.substr(v1nStart), nullptr);

			int v2Slash = v2parts.find("/", 0);
			int v2tStart = v2Slash + 1;
			int v2tSlash = v2parts.find("/", v2tStart);
			int v2nStart = v2tSlash + 1;
			v2 = std::stoi(v2parts.substr(0, v2Slash), nullptr);
			//v2t = std::stoi(v2parts.substr(v2tStart, v2tSlash - v2tStart), nullptr, 16);
			v2n = std::stoi(v2parts.substr(v2nStart), nullptr);

			int v3Slash = v3parts.find("/", 0);
			int v3tStart = v3Slash + 1;
			int v3tSlash = v3parts.find("/", v3tStart);
			int v3nStart = v3tSlash + 1;
			v3 = std::stoi(v3parts.substr(0, v3Slash), nullptr);
			//v3t = std::stoi(v3parts.substr(v3tStart, v3tSlash - v3tStart), nullptr, 16);
			v3n = std::stoi(v3parts.substr(v3nStart), nullptr);

			if (v1 < 0) {
				vert1 = verts.at(verts.size() + v1);
			}
			else {
				vert1 = verts.at(v1 - 1);
			}
			if (v2 < 0) {
				vert2 = verts.at(verts.size() + v2);
			}
			else {
				vert2 = verts.at(v2 - 1);
			}
			if (v3 < 0) {
				vert3 = verts.at(verts.size() + v3);
			}
			else {
				vert3 = verts.at(v3 - 1);
			}
			if (v1n < 0) {
				vert1n = normals.at(normals.size() + v1n);
			}
			else {
				vert1n = normals.at(v1n - 1);
			}
			if (v2n < 0) {
				vert2n = normals.at(normals.size() + v2n);
			}
			else {
				vert2n = normals.at(v2n - 1);
			}
			if (v3n < 0) {
				vert3n = normals.at(normals.size() + v3n);
			}
			else {
				vert3n = normals.at(v3n - 1);
			}
			vert1 = new vertex(vert1->xyz, *vert1n, vert1->color);
			vert2 = new vertex(vert2->xyz, *vert2n, vert2->color);
			vert3 = new vertex(vert3->xyz, *vert3n, vert3->color);
			material *m = new material(shininess, transparency, ior, roughness);
			tri *t = new tri(vert1, vert2, vert3, lastColor, m, objectID, object_type);
			objects.push_back(t);
		}
		else if (!token.compare("output_depth")) {
			output_depth = 1;
		}
		else if (!token.compare("stencil_radius")) {
			ss >> stencil_radius;
		}
		else if (!token.compare("post_process_outline")) {
			post_process_outline = 1;
		}
		else if (!token.compare("stencil_ring_samples")) {
			ss >> stencil_ring_samples;
		}
		else if (!token.compare("stencil_rings")) {
			ss >> stencil_rings;
		}
		else if (!token.compare("outline_cutoff")) {
			ss >> outline_cutoff;
		}
		else if (!token.compare("objectID")) {
			ss >> objectID;
		}
		else if (!token.compare("object_type")) {
			ss >> object_type;
		}
		else {
				//skip line
		}
	}
	std::default_random_engine generator;
	generator.seed(time(NULL));
	std::uniform_real_distribution<double> distribution(-0.5, 0.5);
	std::default_random_engine generator2;
	generator2.seed(time(NULL));
	double depthMap[width * height];
	
	
	for (int yi = 0; yi < height; yi++) {
		for (int xi = 0; xi < width; xi++) {
			milestone = printProgress(xi, yi, width, height, milestone, start);
			resultColor[0] = 0;
			resultColor[1] = 0;
			resultColor[2] = 0;
			resultColor[3] = 0;
			if (spp > 1) {
				for (int sample = 0; sample < spp; sample++) {
					int primary_objID = -1;
					double rand_x = distribution(generator);
					double rand_y = distribution(generator);
					shoot = 1;
					direction = cameraDirection(camera_mode, xi + rand_x, yi + rand_y, width, height, forward, right, up, shoot);
					if (shoot) {
						primary_ray = new ray(eye.x, eye.y, eye.z, direction.x, direction.y, direction.z);
						primary_ray->cast(objects, lights, resultColor, bounces, -1, generator, xi, yi, width, depthMap, post_process_outline, primary_objID, hit_normal, object_type);
						delete primary_ray;
					}
				}
			}
			else {
				shoot = 1;
				direction = cameraDirection(camera_mode, xi, yi, width, height, forward, right, up, shoot);
				if (shoot) {
					primary_ray = new ray(eye.x, eye.y, eye.z, direction.x, direction.y, direction.z);
					int primary_objID = -1;
					int objectType = 0;
					if (primary_ray->cast(objects, lights, resultColor, bounces, -1, generator, xi, yi, width, depthMap, post_process_outline, primary_objID, hit_normal, objectType)) {
						if (objectType == 1) {
						    	std::uniform_real_distribution<double> distributionTheta(0, 360);
							int edge_rays = 0;
							double step_size = stencil_radius / stencil_rings;
							double radius = 0;
							int test_num = 0;
							for (double ring_radius = stencil_radius; ring_radius > 0; ring_radius -= stencil_rings) {
								for (test_num = 0; test_num < stencil_ring_samples; test_num++) {
									radius += step_size;
									double randTheta = distributionTheta(generator2);
									double randx = radius * cos(randTheta * PI / 180);
									double randy = radius * sin(randTheta * PI / 180);
									vec3 newdir = cameraDirection(camera_mode, xi + randx, yi + randy, width, height, forward, right, up, shoot);
									primary_ray->direction.x = newdir.x;
									primary_ray->direction.y = newdir.y;
									primary_ray->direction.z = newdir.z;
									int stencil_objID = -1;
									double edge_test = primary_ray->detect_edge(objects, lights, resultColor, bounces, -1, generator, stencil_objID);
									if (stencil_objID != primary_objID) {
										edge_rays++;
									}
									else if (std::abs(edge_test - depthMap[yi * width + xi])  > outline_cutoff) {
										edge_rays++;
									}
									primary_ray->direction.x = direction.x;
									primary_ray->direction.y = direction.y;
									primary_ray->direction.z = direction.z;
								}
							}
							double edge_strength = 1 - pow((std::abs(edge_rays - 0.5 * stencil_rings * stencil_ring_samples) / (0.5 * stencil_rings * stencil_ring_samples)),2);
							resultColor[0] = resultColor[0] * (1 - edge_strength);
							resultColor[1] = resultColor[1] * (1 - edge_strength);
							resultColor[2] = resultColor[2] * (1 - edge_strength);
						} else if (objectType == 0) {
							vec3 particle_color(resultColor[0], resultColor[1], resultColor[2]);
						    auto paint_particle = std::make_shared<PaintParticle>(xi, yi, particle_color, primary_ray->direction, hit_normal, depthMap[yi * width + xi], 1, 1); // none are edges for now
						    paint_particles.push_back(paint_particle);
						} 
					delete primary_ray;
				}
			}
			}
			resultColor[0] /= spp;
			resultColor[1] /= spp;
			resultColor[2] /= spp;
			resultColor[3] /= spp;
			if (output_depth) {
				myImage(xi, yi, 0, 0) = clip(depthMap[yi * width + xi] * 20);
				myImage(xi, yi, 0, 1) = clip(depthMap[yi * width + xi] * 20);
				myImage(xi, yi, 0, 2) = clip(depthMap[yi * width + xi] * 20);
				myImage(xi, yi, 0, 3) = clip(resultColor[3] * 255);
			}
			else {
				myImage(xi, yi, 0, 0) = clip(resultColor[0] * 255);
				myImage(xi, yi, 0, 1) = clip(resultColor[1] * 255);
				myImage(xi, yi, 0, 2) = clip(resultColor[2] * 255);
				myImage(xi, yi, 0, 3) = clip(resultColor[3] * 255);
			}
		}
	}
	if (post_process_outline) {
		for (int yi = 1; yi < height - 1; yi++) {
			for (int xi = 1; xi < width - 1; xi++) {
				for (int layer = 0; layer < stencil_ring_samples; layer++) {
					if (depthMap[yi * width + xi] - depthMap[yi * width + xi + (layer + 1)] > e_slope) {
						myImage(xi, yi, 0, 0) = 0;
						myImage(xi, yi, 0, 1) = 0;
						myImage(xi, yi, 0, 2) = 0;
					}
					else if (depthMap[yi * width + xi] - depthMap[yi * width + xi - (layer + 1)] > e_slope) {
						myImage(xi, yi, 0, 0) = 0;
						myImage(xi, yi, 0, 1) = 0;
						myImage(xi, yi, 0, 2) = 0;
					}
					else if (depthMap[yi * width + xi] - depthMap[(yi - (layer + 1)) * width + xi] > e_slope) {
						myImage(xi, yi, 0, 0) = 0;
						myImage(xi, yi, 0, 1) = 0;
						myImage(xi, yi, 0, 2) = 0;
					}
					else if (depthMap[yi * width + xi] - depthMap[(yi + (layer + 1)) * width + xi] > e_slope) {
						myImage(xi, yi, 0, 0) = 0;
						myImage(xi, yi, 0, 1) = 0;
						myImage(xi, yi, 0, 2) = 0;
					}
				}

			}
		}	
	}
	
	auto painter = Painter(Expressionist(), paint_particles, &myImage);
        painter.paint();


	inputFile.close();
	myImage.save_png(outputFileName.c_str());
	std::cout << "Progress: 100% done, time elapsed: " << difftime(time(NULL), start) << std::endl;
    	return 0;
}
