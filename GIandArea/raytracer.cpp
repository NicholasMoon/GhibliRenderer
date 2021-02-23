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
#include "octreenode.h"
#include "AABB.h"
#include "particle.h"
#include "painter.h"
#include "style.h"
#include "stroke.h"
#include "brush.h"



#define PI 3.14159265

double printProgress(int pixnum, int totalpixels, double milestone, time_t start) {
	double progress = (double)pixnum / (double)totalpixels;
	if (progress - milestone >= 0.1) {
		std::cout << "Progress: " << (progress * 100) << "% done, time elapsed: " << difftime(time(NULL), start) << std::endl;
		return progress;
	}
	return milestone;
	
}

void getWorldBoundaries(vec3 &min_coordinates, vec3 &max_coordinates, std::vector<object*> &objects) {
	for (int i = 0; i < objects.size(); i++) {
		objects[i]->updateWorldBoundaries(min_coordinates, max_coordinates);
	}
}

double Max_Double(double a, double b, double c) {
	if (a >= b && a >= c) {
		return a;
	}
	else if (b >= a && b >= c) {
		return b;
	}
	return c;
}

double Min_Double(double a, double b, double c) {
	if (a <= b && a <= c) {
		return a;
	}
	else if (b <= a && b <= c) {
		return b;
	}
	return c;
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

vec3 cameraDirection(double xi, double yi, double width, double height, vec3 &forward, vec3 &right, vec3 &up) {
	vec3 camera_direction(0,0,0);
	double sx = (2 * xi - width) / maxItoD(width,height);
	double sy = (height - 2 * yi) / maxItoD(width,height);
	camera_direction.x = forward.x + sx * right.x + sy * up.x;
	camera_direction.y = forward.y + sx * right.y + sy * up.y;
	camera_direction.z = forward.z + sx * right.z + sy * up.z;
	camera_direction.normalize();
	return camera_direction;
}

void getDrawingGradient(vec3 hit_normal, double *old_color, double *light_x_color, double *light_y_color, vec3 &final_gradient, int pixnum) {
    // See which is stronger - edge or light gradient and draw strokes perpendicular to this

    // Edge gradient
    double ex = hit_normal.x;
    double ey = hit_normal.y;
    vec3 edge_gradient(ex, ey, 0);

    // Light gradient 
    double old_light_total = old_color[0] + old_color[1] + old_color[2];
    double x_light_total = light_x_color[0] + light_x_color[1] + light_x_color[2];
    double y_light_total = light_y_color[0] + light_y_color[1] + light_y_color[2];

    vec3 light_gradient(x_light_total - old_light_total, y_light_total - old_light_total, 0);

    // Weighting
    double edge = edge_gradient.magnitude();
    double light = std::min(2.0 * light_gradient.magnitude(), pow(1 - edge, 2));

    // See if either vector is 0 length
    if (edge_gradient.magnitude() != 0) edge_gradient.normalize();
    if (light_gradient.magnitude() != 0) light_gradient.normalize();

    // final_gradient = (edge_gradient * edge + light_gradient * light).normalize();
	final_gradient = edge_gradient;
	// if (edge >= 0.6) {
    //     final_gradient = edge_gradient;
    // }
    // else {
    //     final_gradient = light_gradient;
    // }
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
	double eccentricity = 0.0;
	int v1, v2, v3, v4, v1n, v2n, v3n, v4n;
	vertex *vert1, *vert2, *vert3, *vert4;
	double lastNormal[3] = {0,0,0};
	double lastColor[3] = {1, 1, 1};
	double lastEmission[3] = {0, 0, 0};
	double resultColor[4] = {0, 0, 0, 0};
	double lightDxColor[4] = {0, 0, 0, 0};
	double lightDyColor[4] = {0, 0, 0, 0};
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
	vec3 light_dx_right(0,0,0);
	vec3 light_dy_up(0,0,0);
	vec3 stroke_gradient(0,0,0);
	int stencil_ring_samples = 0;
	double stencil_rings = 0;
	double outline_cutoff = 1;
	int shoot = 1;
	int spp = 1;
	int output_depth = 0;
	double e_slope = 1;
	double stencil_radius = 0;
	int light_samples = 1;
	int indirect_samples = 0;
	int indirect_bounces = 0;
	int flat = 0;
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
	ray *light_dx_ray, *light_dy_ray;
	stroke *paint_stroke;
	brush *paint_brush;
	brush *const_brush_small, *const_brush_medium, *const_brush_large;

	// Initialize stroke lengths (do this based on image resolution)
    int small_size = 1;
    int medium_size = 10;
    int large_size = 20;
    int strokeLengths[3] = {small_size, medium_size, large_size};

    // Initialize brush set
    const_brush_small = new constbrush(1);
    const_brush_small->create_mask();

    const_brush_medium = new constbrush(4);
    const_brush_medium->create_mask();

    const_brush_large = new constbrush(10);
    const_brush_large->create_mask();

    brush *brushSet[3] = {const_brush_small, const_brush_medium, const_brush_large};
	
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
			material *m = new material(shininess, transparency, ior, roughness, eccentricity);
			// object *s = new sphere(x, y, z, r, lastColor, m, objectID, object_type);
			object *s = new sphere(x, y, z, r, lastColor, lastEmission, m, objectID, object_type);
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
		else if (!token.compare("emission")) {
			ss >> lastEmission[0];
			ss >> lastEmission[1];
			ss >> lastEmission[2];
		}
		else if (!token.compare("plane")) {
			ss >> A;
			ss >> B;
			ss >> C;
			ss >> D;
			material *m = new material(shininess, transparency, ior, roughness, eccentricity);
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
			material *m = new material(shininess, transparency, ior, roughness, eccentricity);
			// tri *t = new tri(vert1, vert2, vert3, lastColor, m, objectID, object_type);
			tri *t = new tri(vert1, vert2, vert3, lastColor, lastEmission, m, objectID, object_type); 
			objects.push_back(t);
		}
		else if (!token.compare("stencil_radius")) {
			ss >> stencil_radius;
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
		else if (!token.compare("indirect_samples")) {
			ss >> indirect_samples;
		}
		else if (!token.compare("indirect_bounces")) {
			ss >> indirect_bounces;
		}
		else if (!token.compare("light_samples")) {
			ss >> light_samples;
		}
		else if (!token.compare("area_light")) {
			ss >> v1;
			ss >> v2;
			ss >> v3;
			ss >> v4;
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
			if (v4 < 0) {
				vert4 = verts.at(verts.size() + v4);
			}
			else {
				vert4 = verts.at(v4 - 1);
			}
			light *l = new light(3, vert1->xyz, vert2->xyz, vert3->xyz, vert4->xyz, lastColor);
			lights.push_back(l);
		}
		else if (!token.compare("flat")) {
			ss >> flat;
		}
		else if (!token.compare("eccentricity")) {
			ss >> eccentricity;
		}
		else if (!token.compare("spp")) {
			ss >> spp;
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
	double edgeMap[width * height][3];
	int objectTypeMap[width * height];
	int objectBoundaryMap[width * height];
	int objectIDMap[width * height];
	int shadowMap[width * height];
	float paintMap[width * height]; // paint density values

	for (int pm = 0; pm < width * height; pm++) {
		paintMap[pm] = 0;
	}

	// Randomize pixels
	std::vector<int> x_values;
    std::vector<int> y_values;

    for (int i=0; i < width; i++) {
        x_values.push_back(i);
    }
    for (int j=0; j < height; j++) {
        y_values.push_back(j);
    }
    auto rng = std::default_random_engine {};
    std::shuffle(x_values.begin(), x_values.end(), rng);
    std::shuffle(y_values.begin(), y_values.end(), rng);
    std::vector<std::tuple<int, int>> pixels;
    for (int x : x_values) {
        for (int y : y_values) {
            pixels.push_back(std::make_tuple (x, y));
        }
    }
    std::shuffle(pixels.begin(), pixels.end(), rng);
	int pixnum = 0;
	int totalpixels = width * height;
	for (std::tuple<int, int> curr_pixel : pixels) {
        int xi = std::get<0>(curr_pixel);
        int yi = std::get<1>(curr_pixel);
		pixnum++;
		milestone = printProgress(pixnum, totalpixels, milestone, start);
		resultColor[0] = 0;
		resultColor[1] = 0;
		resultColor[2] = 0;
		resultColor[3] = 0;
		lightDxColor[0] = 0;
		lightDxColor[1] = 0;
		lightDxColor[2] = 0;
		lightDxColor[3] = 0;
		lightDyColor[0] = 0;
		lightDyColor[1] = 0;
		lightDyColor[2] = 0;
		lightDyColor[3] = 0;
		double AAColor[4];
		AAColor[0] = 0;
		AAColor[1] = 0;
		AAColor[2] = 0;
		AAColor[3] = 0;
		int background_samples = 0;
		int foreground_samples = 0;
		for (int sample = 0; sample < spp; sample++) {
			if (sample > 0) {
				double rand_x = distribution(generator);
				double rand_y = distribution(generator);
				direction = cameraDirection(xi + rand_x, yi + rand_y, width, height, forward, right, up);
			}
			else {
				direction = cameraDirection(xi, yi, width, height, forward, right, up);
			}
			primary_ray = new ray(eye.x, eye.y, eye.z, direction.x, direction.y, direction.z);
			int primary_objID = -1;
			int objectType = 0;
			std::vector<int> hit_list;
			int shadowed = 0;
			if (primary_ray->cast(objects, lights, resultColor, bounces, -1, generator, 1, xi, yi, width, depthMap, primary_objID, hit_normal, objectType, hit_list, shadowed, flat, light_samples, indirect_samples, indirect_bounces)) {
				if (shadowed) {
					shadowMap[yi * width + xi] = 1;
				}
				if (objectType == 1) {
					foreground_samples++;
				}
				else {
					background_samples++;
				}
				objectTypeMap[yi * width + xi] = objectType;
				objectBoundaryMap[yi * width + xi] = hit_list.front();
				objectIDMap[yi * width + xi] = primary_objID;
				std::uniform_real_distribution<double> distributionTheta(0, 360);
				int edge_rays = 0;
				double step_size = stencil_radius / stencil_rings;
				double radius = 0;
				int test_num = 0;
				// resultColor[0] = 0;
				// resultColor[1] = 0;
				// resultColor[2] = 0;
				// resultColor[3] = 0;
				for (double ring_radius = stencil_radius; ring_radius > 0; ring_radius -= stencil_rings) {
					for (test_num = 0; test_num < stencil_ring_samples; test_num++) {
						int edgeObjectType = 0;
						radius += step_size;
						double randTheta = distributionTheta(generator2);
						double randx = radius * cos(randTheta * PI / 180);
						double randy = radius * sin(randTheta * PI / 180);
						vec3 newdir = cameraDirection(xi + randx, yi + randy, width, height, forward, right, up);
						primary_ray->direction.x = newdir.x;
						primary_ray->direction.y = newdir.y;
						primary_ray->direction.z = newdir.z;
						int stencil_objID = -1;
						std::vector<int> edge_hit_list;
						double edge_test = primary_ray->detect_edge(objects, lights, resultColor, bounces, -1, generator, stencil_objID, edgeObjectType, edge_hit_list);
						primary_ray->direction.x = direction.x;
						primary_ray->direction.y = direction.y;
						primary_ray->direction.z = direction.z;
						if (hit_list.size() == 0 && edge_hit_list.size() > 0) {
							edge_rays++;
						}
						else if  (hit_list.size() == 0 && edge_hit_list.size() == 0) {
							continue;
						}
						else if  (hit_list.size() > 0 && edge_hit_list.size() == 0) {
							edge_rays++;
						}
						else if (edge_hit_list[0] != hit_list[0] && edge_hit_list.size() != hit_list.size()) {
							edge_rays++;
						}
						else if (edgeObjectType == 0 && objectType == 0) {
							continue;
						}
						else if (stencil_objID != primary_objID) {
							edge_rays++;
						}
						else if (std::abs(edge_test - depthMap[yi * width + xi])  > outline_cutoff) {
							edge_rays++;
						}
					}
				}
				double edge_strength = 1 - pow((std::abs(edge_rays - 0.5 * stencil_rings * stencil_ring_samples) / (0.5 * stencil_rings * stencil_ring_samples)),2);
				//resultColor[0] = resultColor[0] * (1 - edge_strength);
				//resultColor[1] = resultColor[1] * (1 - edge_strength);
				//resultColor[2] = resultColor[2] * (1 - edge_strength);
				edgeMap[yi * width + xi][0] = 1 - edge_strength;
				edgeMap[yi * width + xi][1] = 1 - edge_strength;
				edgeMap[yi * width + xi][2] = 1 - edge_strength;

				if (objectType == 0) {
					vec3 particle_color(resultColor[0], resultColor[1], resultColor[2]);
					// auto paint_particle = std::make_shared<PaintParticle>(xi, yi, particle_color, primary_ray->direction, hit_normal, depthMap[yi * width + xi], 1, 1); // none are edges for now
					// paint_particles.push_back(paint_particle);

					// Get drawing gradient
					light_dx_right.x = right.x * 0.001;
					light_dx_right.y = right.y * 0.001;
					light_dx_right.z = right.z * 0.001;
					light_dy_up.x = up.x * 0.001;
					light_dy_up.y = up.y * 0.001;
					light_dy_up.z = up.z * 0.001;
					// vec3 light_dx_dir = cameraDirection(xi, yi, width, height, forward, light_dx_right, up); // adjust right vector
					// vec3 light_dy_dir = cameraDirection(xi, yi, width, height, forward, right, light_dy_up); // adjust up vector
					light_dx_ray = new ray(eye.x, eye.y, eye.z, direction.x + light_dx_right.x, direction.y + light_dx_right.y, direction.z + light_dx_right.z);
					light_dy_ray = new ray(eye.x, eye.y, eye.z, direction.x + light_dy_up.x, direction.y + light_dy_up.y, direction.z + light_dy_up.z);

					int lightObjectType = 0;
					int light_objID = -1;
					std::vector<int> light_hit_list;
					vec3 light_hit_normal(0,0,0);
					light_dx_ray->cast(objects, lights, lightDxColor, bounces, -1, generator, 1, xi, yi, width, depthMap, light_objID, light_hit_normal, lightObjectType, light_hit_list, shadowed, flat, light_samples, indirect_samples, indirect_bounces);
					light_dy_ray->cast(objects, lights, lightDyColor, bounces, -1, generator, 1, xi, yi, width, depthMap, light_objID, light_hit_normal, lightObjectType, light_hit_list, shadowed, flat, light_samples, indirect_samples, indirect_bounces);

					delete light_dx_ray;
					delete light_dy_ray;

					getDrawingGradient(hit_normal, resultColor, lightDxColor, lightDyColor, stroke_gradient, pixnum);
					double cosTheta = hit_normal.dot(primary_ray->direction);
					bool inside = false;
					if (cosTheta > 0) inside = true;
					
					// Make new stroke
					paint_stroke = new stroke(xi, yi, particle_color, primary_ray->direction, hit_normal, depthMap[yi * width + xi]);
					// Set these based on position and distance from camera
					paint_stroke->set_length(strokeLengths[1]);
					paint_stroke->set_curvature(0.8);
					paint_stroke->create(stroke_gradient, inside, 0, 0, myImage.width(), myImage.height());

					// Choose brush and paint
					paint_brush = brushSet[1];
					paint_brush->paint(paint_stroke, paintMap, objectTypeMap, objectBoundaryMap, objectIDMap, hit_list.front(), primary_objID, &myImage);
					delete paint_stroke;
				} 
				AAColor[0] += resultColor[0];
				AAColor[1] += resultColor[1];
				AAColor[2] += resultColor[2];
				AAColor[3] += resultColor[3];
				resultColor[0] = 0;
				resultColor[1] = 0;
				resultColor[2] = 0;
				resultColor[3] = 0;
				delete primary_ray;
			}
			if (foreground_samples > background_samples) {
				myImage(xi, yi, 0, 0) = clip((AAColor[0] / (double) spp) * 255);
				myImage(xi, yi, 0, 1) = clip((AAColor[1] / (double) spp) * 255);
				myImage(xi, yi, 0, 2) = clip((AAColor[2] / (double) spp) * 255);
				myImage(xi, yi, 0, 3) = clip(AAColor[3] * 255);
			}
		}
	}
	
	// auto painter = Painter(Expressionist(), paint_particles, &myImage);
    // painter.paint(objectTypeMap);

	for (int yk = 0; yk < height; yk++) { 
		for (int xk = 0; xk < width; xk++) {
			if (shadowMap[yk * width + xk] && objectTypeMap[yk * width + xk] == 0) {
				myImage(xk, yk, 0, 0) = clip(myImage(xk, yk, 0, 0) * .35);
				myImage(xk, yk, 0, 1) = clip(myImage(xk, yk, 0, 1) * .35);
				myImage(xk, yk, 0, 2) = clip(myImage(xk, yk, 0, 2) * .35);
			}
		}
	}

	for (int yj = 0; yj < height; yj++) { // this affects emissive objects
		for (int xj = 0; xj < width; xj++) {
			myImage(xj, yj, 0, 0) = clip(myImage(xj, yj, 0, 0) * edgeMap[yj * width + xj][0]);
			myImage(xj, yj, 0, 1) = clip(myImage(xj, yj, 0, 1) * edgeMap[yj * width + xj][1]);
			myImage(xj, yj, 0, 2) = clip(myImage(xj, yj, 0, 2) * edgeMap[yj * width + xj][2]);
		}
	}


	inputFile.close();
	myImage.save_png(outputFileName.c_str());
	std::cout << "Progress: 100% done, time elapsed: " << difftime(time(NULL), start) << std::endl;
    	return 0;
}
