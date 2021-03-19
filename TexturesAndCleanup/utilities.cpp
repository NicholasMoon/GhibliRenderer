#include "utilities.h"

std::default_random_engine generator;
std::uniform_real_distribution<float> distribution(0.0, 1.0);

float random_number() {
    return distribution(generator);
}

double distance(double x, double x0, double y, double y0) {
    double x2 = (x - x0) * (x - x0);
    double y2 = (y - y0) * (y - y0);
    return sqrt(x2 + y2);
}

double sample_parabola(double a, double x) { 
    return a * (x * x);
}

// void alpha_composite(vec3 &c1, vec3 &c2, double a1, double a2, vec3 &final_color, double &final_alpha) {
//     final_alpha = a1 + a2 * (1-a1);
//     final_color = (c1 * a1) + c2 * (1.0 - a1);
//     // final_color.print();
//     // std::cout << final_alpha << std::endl;
// }

float clamp(float color, float min, float max) {
    float c = color;
    if (c >= max) {
        c = max;
    } else if (c <= min) {
        c = min;
    }

    return c;
}

void write_color(CImg<float> &img, int x, int y, vec3 color, int alpha) {
    img(x, y, 0, 0) = clamp(color.x, 0, 1) * 255;
    img(x, y, 0, 1) = clamp(color.y, 0, 1) * 255;
    img(x, y, 0, 2) = clamp(color.z, 0, 1) * 255;
    img(x, y, 0, 3) = alpha;
}

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

vec3 cameraDirection(double xi, double yi, Scene *theScene) {
	vec3 camera_direction(0,0,0);
	double sx = (2 * xi - theScene->width) / maxItoD(theScene->width,theScene->height);
	double sy = (theScene->height - 2 * yi) / maxItoD(theScene->width,theScene->height);
	camera_direction.x = theScene->forward.x + sx * theScene->right.x + sy * theScene->up.x;
	camera_direction.y = theScene->forward.y + sx * theScene->right.y + sy * theScene->up.y;
	camera_direction.z = theScene->forward.z + sx * theScene->right.z + sy * theScene->up.z;
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
