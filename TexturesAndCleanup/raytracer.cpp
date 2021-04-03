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

// no includes
#include "vec3.h"
#include "image_buffers.h"
// includes only vec3
#include "vertex.h"
#include "light.h"
#include "material.h"
#include "style.h"
#include "particle.h"
// include multiple
#include "object.h"
#include "scene.h"
#include "sphere.h"
#include "plane.h"
#include "tri.h"
#include "texture.h"
#include "colortexture.h"
#include "imagetexture.h"
#include "ray.h"
#include "octreenode.h"
#include "AABB.h"
#include "painter.h"
#include "stroke.h"
#include "brush.h"
#include "utilities.h"




#define PI 3.14159265


int main(int argc, char** argv) {
	time_t start = time(0);
	if (argc != 2) {
		std::cout << "You must provide input file name as command line argument!" << std::endl;
		exit(1);
	}
	double milestone = 0.0;
	double resultColor[4] = {0,0,0,0};
	double lightDxColor[4] = {0, 0, 0, 0};
	double lightDyColor[4] = {0, 0, 0, 0};
	double sx, sy = 0;
	vec3 direction(0,0,0);
	vec3 hit_normal(0,0,0);
	vec3 light_dx_right(0,0,0);
	vec3 light_dy_up(0,0,0);
	vec3 stroke_gradient(0,0,0);
	
	Scene *theScene = new Scene(argv);

	cimg_library::CImg<float> myImage(theScene->width, theScene->height, 1, 4, 0);
    	unsigned char zero = 0;
    	myImage = zero;

	std::vector<stroke *> bottomLayer;
	std::vector<stroke *> middleLayer;
	std::vector<stroke *> topLayer;

	ImageBuffers *imageBuffers = new ImageBuffers(theScene->dimension);
	
	// std::vector<std::shared_ptr<PaintParticle>> paint_particles;
    	ray *primary_ray;
	ray *light_dx_ray, *light_dy_ray;
	stroke *paint_stroke;
	// brush *paint_brush;
	// brush *const_brush_small, *const_brush_medium, *const_brush_large;

	// Initialize stroke lengths (do this based on image resolution)
    int small_size = 1;
    int medium_size = 10;
    int large_size = 20;
    int strokeLengths[3] = {small_size, medium_size, large_size};

    // // Initialize brush set
    // const_brush_small = new brush(1, 2);
    // const_brush_small->create_mask();

    // const_brush_medium = new brush(4, 2);
    // const_brush_medium->create_mask();

    // const_brush_large = new brush(10, 2);
    // const_brush_large->create_mask();

    // brush *brushSet[3] = {const_brush_small, const_brush_medium, const_brush_large};
	
	std::uniform_real_distribution<double> distribution(-0.5, 0.5);
	std::default_random_engine generator2;
	generator2.seed(time(NULL));

	// // Randomize pixels
	// std::vector<int> x_values;
    // std::vector<int> y_values;

    // for (int i=0; i < theScene->width; i++) {
    //     x_values.push_back(i);
    // }
    // for (int j=0; j < theScene->height; j++) {
    //     y_values.push_back(j);
    // }
    // auto rng = std::default_random_engine {};
    // std::shuffle(x_values.begin(), x_values.end(), rng);
    // std::shuffle(y_values.begin(), y_values.end(), rng);
    // std::vector<std::tuple<int, int>> pixels;
    // for (int x : x_values) {
    //     for (int y : y_values) {
    //         pixels.push_back(std::make_tuple (x, y));
    //     }
    // }
    // std::shuffle(pixels.begin(), pixels.end(), rng);
	int pixnum = 0;
	for (int yi=0; yi < theScene->height; yi++) {
		for (int xi=0; xi < theScene->width; xi++) {
	// for (std::tuple<int, int> curr_pixel : pixels) {
    //     int xi = std::get<0>(curr_pixel);
    //     int yi = std::get<1>(curr_pixel);
			pixnum++;
			milestone = printProgress(pixnum, theScene->dimension, milestone, start);
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
			for (int sample = 0; sample < theScene->spp; sample++) {
				if (sample > 1) {
					double rand_x = distribution(theScene->generator);
					double rand_y = distribution(theScene->generator);
					direction = cameraDirection(xi + rand_x, yi + rand_y, theScene);
				}
				else {
					direction = cameraDirection(xi, yi, theScene);
				}
				primary_ray = new ray(theScene->eye.x, theScene->eye.y, theScene->eye.z, direction.x, direction.y, direction.z);
				int primary_objID = -1;
				int objectType = 0;
				std::vector<int> hit_list;
				int shadowed = 0;
				resultColor[0] = 0;
				resultColor[1] = 0;
				resultColor[2] = 0;
				resultColor[3] = 0;
				int index = yi * theScene->width + xi;
				int starting_bounces = theScene->bounces;
				int starting_indirect_bounces = theScene->indirect_bounces;
				if (primary_ray->cast(theScene, imageBuffers, resultColor, starting_bounces, -1, 1, xi, yi, primary_objID, hit_normal, objectType, hit_list, shadowed, starting_indirect_bounces)) {
					if (shadowed) {
						imageBuffers->shadowMap[index] = 1;
					}
					if (objectType == 1) {
						foreground_samples++;
					}
					else {
						background_samples++;
					}
					imageBuffers->objectTypeMap[index] = objectType;
					imageBuffers->objectBoundaryMap[index] = hit_list.front();
					imageBuffers->objectIDMap[index] = primary_objID;
					std::uniform_real_distribution<double> distributionTheta(0, 360);
					int edge_rays = 0;
					double step_size = theScene->stencil_radius / theScene->stencil_rings;
					double radius = 0;
					int test_num = 0;
					if (objectType == 1) {
						for (double ring_radius = theScene->stencil_radius; ring_radius > 0; ring_radius -= theScene->stencil_rings) {
							for (test_num = 0; test_num < theScene->stencil_ring_samples; test_num++) {
								int edgeObjectType = 0;
								radius += step_size;
								double randTheta = distributionTheta(generator2);
								double randx = radius * cos(randTheta * PI / 180);
								double randy = radius * sin(randTheta * PI / 180);
								vec3 newdir = cameraDirection(xi + randx, yi + randy, theScene);
								primary_ray->direction.x = newdir.x;
								primary_ray->direction.y = newdir.y;
								primary_ray->direction.z = newdir.z;
								int stencil_objID = -1;
								std::vector<int> edge_hit_list;
								double edge_test = primary_ray->detect_edge(theScene, resultColor, starting_bounces, -1, stencil_objID, edgeObjectType, edge_hit_list);
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
								else if (std::abs(edge_test - imageBuffers->depthMap[index]) > theScene->outline_cutoff) {
									edge_rays++;
								}
							}
						}
					}
					double edge_strength = 1 - pow((std::abs(edge_rays - 0.5 * theScene->stencil_rings * theScene->stencil_ring_samples) / (0.5 * theScene->stencil_rings * theScene->stencil_ring_samples)),2);
					imageBuffers->edgeMap[index] = 1 - edge_strength;
					imageBuffers->edgeMap[(yi + theScene->height) * theScene->width + xi] = 1 - edge_strength;
					imageBuffers->edgeMap[(yi + 2*theScene->height) * theScene->width + xi] = 1 - edge_strength;

					if (objectType == 0) {
						vec3 particle_color(resultColor[0], resultColor[1], resultColor[2]);

						// Get drawing gradient
						light_dx_right.x = theScene->right.x * 0.001;
						light_dx_right.y = theScene->right.y * 0.001;
						light_dx_right.z = theScene->right.z * 0.001;
						light_dy_up.x = theScene->up.x * 0.001;
						light_dy_up.y = theScene->up.y * 0.001;
						light_dy_up.z = theScene->up.z * 0.001;
						// vec3 light_dx_dir = cameraDirection(xi, yi, width, height, forward, light_dx_right, up); // adjust right vector
						// vec3 light_dy_dir = cameraDirection(xi, yi, width, height, forward, right, light_dy_up); // adjust up vector
						light_dx_ray = new ray(theScene->eye.x, theScene->eye.y, theScene->eye.z, direction.x + light_dx_right.x, direction.y + light_dx_right.y, direction.z + light_dx_right.z);
						light_dy_ray = new ray(theScene->eye.x, theScene->eye.y, theScene->eye.z, direction.x + light_dy_up.x, direction.y + light_dy_up.y, direction.z + light_dy_up.z);

						int lightObjectType = 0;
						int light_objID = -1;
						std::vector<int> light_hit_list;
						vec3 light_hit_normal(0,0,0);
						light_dx_ray->cast(theScene, imageBuffers, lightDxColor, starting_bounces, -1, 1, xi, yi, light_objID, light_hit_normal, lightObjectType, light_hit_list, shadowed, starting_indirect_bounces);
						light_dy_ray->cast(theScene, imageBuffers, lightDyColor, starting_bounces, -1, 1, xi, yi, light_objID, light_hit_normal, lightObjectType, light_hit_list, shadowed, starting_indirect_bounces);

						delete light_dx_ray;
						delete light_dy_ray;

						getDrawingGradient(hit_normal, resultColor, lightDxColor, lightDyColor, stroke_gradient, pixnum);
						double cosTheta = hit_normal.dot(primary_ray->direction);
						bool inside = false;
						if (cosTheta > 0) inside = true;
						
						// Make new stroke
						paint_stroke = new stroke(xi, yi, particle_color, primary_ray->direction, hit_normal, imageBuffers->depthMap[index], hit_list.front(), primary_objID);
						// Set these based on position and distance from camera
						paint_stroke->set_length(strokeLengths[1]);
						paint_stroke->set_curvature(0.8);
						paint_stroke->create(stroke_gradient, inside, 0, 0, myImage.width(), myImage.height());

						bottomLayer.push_back(paint_stroke);
						
						// // Choose brush and paint
						// paint_brush = brushSet[1];
						// paint_brush->paint(paint_stroke, imageBuffers, &myImage);
						// delete paint_stroke;
					} 
				}
				else {
					imageBuffers->edgeMap[index] = 1;
					imageBuffers->edgeMap[(yi + theScene->height) * theScene->width + xi] = 1;
					imageBuffers->edgeMap[(yi + 2*theScene->height) * theScene->width + xi] = 1;
				}
				AAColor[0] += resultColor[0];
				AAColor[1] += resultColor[1];
				AAColor[2] += resultColor[2];
				AAColor[3] += resultColor[3];
				delete primary_ray;
			}
			if (foreground_samples > background_samples || theScene->environment || imageBuffers->environmentMap[yi * theScene->width + xi]) {
				myImage(xi, yi, 0, 0) = clip((AAColor[0] / (double) theScene->spp) * 255);
				myImage(xi, yi, 0, 1) = clip((AAColor[1] / (double) theScene->spp) * 255);
				myImage(xi, yi, 0, 2) = clip((AAColor[2] / (double) theScene->spp) * 255);
				myImage(xi, yi, 0, 3) = clip(AAColor[3] * 255);
			}
		}
	}
	
	painter *background_painter = new painter(theScene->width, theScene->height, bottomLayer, middleLayer, topLayer);
    background_painter->paint(imageBuffers, &myImage);

	for (int yk = 0; yk < theScene->height; yk++) { 
		for (int xk = 0; xk < theScene->width; xk++) {
			if (imageBuffers->shadowMap[yk * theScene->width + xk] && imageBuffers->objectTypeMap[yk * theScene->width + xk] == 0) {
				myImage(xk, yk, 0, 0) = clip(myImage(xk, yk, 0, 0) * .35);
				myImage(xk, yk, 0, 1) = clip(myImage(xk, yk, 0, 1) * .35);
				myImage(xk, yk, 0, 2) = clip(myImage(xk, yk, 0, 2) * .35);
			}
		}
	}

	for (int yj = 0; yj < theScene->height; yj++) { 
		for (int xj = 0; xj < theScene->width; xj++) {
			myImage(xj, yj, 0, 0) = clip(myImage(xj, yj, 0, 0) * imageBuffers->edgeMap[yj * theScene->width + xj]);
			myImage(xj, yj, 0, 1) = clip(myImage(xj, yj, 0, 1) * imageBuffers->edgeMap[(yj + theScene->height) * theScene->width + xj]);
			myImage(xj, yj, 0, 2) = clip(myImage(xj, yj, 0, 2) * imageBuffers->edgeMap[(yj + 2*theScene->height) * theScene->width + xj]);
		}
	}


	theScene->inputFile.close();
	std::string image_directory = "images/";
    image_directory.append(theScene->outputFileName);
	myImage.save_png(image_directory.c_str());
	std::cout << "Progress: 100% done, time elapsed: " << difftime(time(NULL), start) << std::endl;
    	return 0;
}
