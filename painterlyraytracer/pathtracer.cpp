// To compile this file, type the following:
//   g++ pathtracer.cpp -lpng -lX11 -lpthread -o output
// To run, type the following:
//   ./output files/test_ground_objects_indirect.txt

#define cimg_use_png

#include <sstream>
#include <fstream>
#include <string>

#include "utilities.h"
#include "objects.h"
#include "lights.h"
#include "painter.h"

// IMAGE

int width = 0; 
int height = 0; 
int depth = 1;
int spectrum = 4;
std::string filename;
CImg<float> img;

// SCENE

Vec3 eye = Vec3(0, 0, 0);
Vec3 forward = Vec3(0, 0, -1);
Vec3 right = Vec3(1, 0, 0);
Vec3 up = Vec3(0, 1, 0);

Vec3 black = Vec3(0, 0, 0);
Vec3 gray = Vec3(0.5, 0.5, 0.5);
Vec3 white = Vec3(1, 1, 1);
float gamma_val = 1.0;
bool global_illumination = false;
int secondary_rays = 1;
int bounces = 2;

// PAINTING
bool background = false;
bool paint = false;

std::vector<Vec3> color_list;
std::vector<Vec3> emission_values;
std::vector<std::shared_ptr<Light>> scene_lights; 
std::vector<std::shared_ptr<SceneObj>> scene_objects;  
std::vector<std::shared_ptr<PaintParticle>> particles; 

void sphere(std::vector<std::string> &tokens) {
    float x = std::stof(tokens.at(1));
    float y = std::stof(tokens.at(2));
    float z = std::stof(tokens.at(3));
    Vec3 center = Vec3(x, y, z);
    float radius = std::stof(tokens.at(4));

    // Create new sphere and add to render list
    auto new_sphere = std::make_shared<Sphere>(center, radius);
    new_sphere->material = 1; // make a material class
    new_sphere->type = 1;
    if (background) {
        new_sphere->type = 2; 
    }
    new_sphere->color = color_list.back();
    new_sphere->emission = emission_values.back();
    scene_objects.push_back(new_sphere);
}

void color(std::vector<std::string> &tokens) {
    float r = std::stof(tokens.at(1));
    float g = std::stof(tokens.at(2));
    float b = std::stof(tokens.at(3));
    Vec3 new_color = Vec3(r, g, b);
    
    color_list.push_back(new_color);
}

void emission(std::vector<std::string> &tokens) {
    float e1 = std::stof(tokens.at(1));
    float e2 = std::stof(tokens.at(2));
    float e3 = std::stof(tokens.at(3));
    Vec3 new_emission = Vec3(e1, e2, e3);
    
    emission_values.push_back(new_emission);
}

void sun(std::vector<std::string> &tokens) {
    float x = std::stof(tokens.at(1));
    float y = std::stof(tokens.at(2));
    float z = std::stof(tokens.at(3));
    Vec3 direction = Vec3(x, y, z);
    Vec3 color = color_list.back();

    // Create new light and add to scene lights
    auto new_light = std::make_shared<Sun>(direction, color);
    scene_lights.push_back(new_light);
}

void bulb(std::vector<std::string> &tokens) {
    float x = std::stof(tokens.at(1));
    float y = std::stof(tokens.at(2));
    float z = std::stof(tokens.at(3));
    Vec3 position = Vec3(x, y, z);
    Vec3 color = color_list.back();

    // Create new light and add to scene lights
    auto new_light = std::make_shared<Bulb>(position, color);
    scene_lights.push_back(new_light);
}

void gi(std::vector<std::string> &tokens) {
    global_illumination = true;
}

void spp(std::vector<std::string> &tokens) {
    secondary_rays = std::stoi(tokens.at(1));
}

void change_eye(std::vector<std::string> &tokens) {
    float x = std::stof(tokens.at(1));
    float y = std::stof(tokens.at(2));
    float z = std::stof(tokens.at(3));

    eye = Vec3(x, y, z);
}

void add_to_background(std::vector<std::string> &tokens) {
    background = true;
}

void add_to_foreground(std::vector<std::string> &tokens) {
    background = false;
}

void paintstyle(std::vector<std::string> &tokens) {
    paint = true;
}

void create_png(std::vector<std::string> &tokens) {
    width = std::stoi(tokens.at(1));
    height = std::stoi(tokens.at(2));
    filename = tokens.at(3);

    // Create image where each pixel is a float and set all pixels to 0
    img.assign(width, height, depth, spectrum, 0); // make member function of struct

    // Set default color and emission values
    color_list.push_back(Vec3(1, 1, 1));
    emission_values.push_back(Vec3(0, 0, 0));
}

int process_tokens(std::vector<std::string> &tokens) {
    std::string keyword = tokens.at(0);

    if (keyword == "png") {
        create_png(tokens);
    } else if (keyword == "sphere") {
        sphere(tokens);
    } else if (keyword == "color") {
        color(tokens);
    } else if (keyword == "emission") {
        emission(tokens);
    } else if (keyword == "sun") {
        sun(tokens);
    } else if (keyword == "bulb") {
        bulb(tokens);
    } else if (keyword == "gi") {
        gi(tokens);
    } else if (keyword == "spp") {
        spp(tokens);
    } else if (keyword == "eye") {
        change_eye(tokens);
    } else if (keyword == "background") {
        add_to_background(tokens);
    } else if (keyword == "foreground") {
        add_to_foreground(tokens);
    } else if (keyword == "paint") {
        paintstyle(tokens);
    }
    return 0;
}

// RAYTRACING

Ray get_ray(int x, int y) {
    Ray ray;

    float denominator = std::max(width, height);
    float sx = (2 * x - width) / denominator;
    float sy = (height - 2 * y) / denominator;

    Vec3 sr_su = (right * sx) + (up * sy);
    ray.direction = (forward + sr_su).norm();
    ray.origin = eye;

    return ray;
}

bool get_intersection(Ray &ray, std::shared_ptr<SceneObj> &draw_object, Intersection &hit) {
    float min_distance = std::numeric_limits<float>::max(); 

    bool intersect = false;

    for (unsigned long int i=0; i < scene_objects.size(); i++) {
        Intersection intersection;
        bool found_intersection = scene_objects.at(i)->intersect(ray, intersection);
        if (found_intersection) {  
            intersect = true;
            float t = intersection.distance;
            if (t < min_distance) {
                draw_object = scene_objects.at(i);
                hit = intersection;
                min_distance = t;
            }
        }
    }
    return intersect;
}

int calc_shadows(Ray &shadow_ray, std::shared_ptr<Light> &light, Intersection &hit) {
    double light_distance = light->distance(hit);
    for (unsigned long int i=0; i < scene_objects.size(); i++) {
        auto shadow_object = scene_objects.at(i);
        Intersection shadow_hit;
        bool found_intersection = scene_objects.at(i)->intersect(shadow_ray, shadow_hit);
        if (found_intersection) {
            if (shadow_hit.distance < light_distance) {
                return 0;
            }
        }
    }
    return 1;
}

Vec3 sample_hemisphere(float r1, float r2) {
    float r = sqrt(1 - r1*r1);
    float phi = 2 * M_PI * r2;
    float x = r * cos(phi);
    float z = r * sin(phi);
    return Vec3(x, r1, z);
}

void create_orthonormal_system(Vec3 &n, Vec3 &rotX, Vec3 &rotY) {
    Vec3 closest_vec;
    if (n.x > n.y) {
        closest_vec = Vec3(n.z, 0, -n.x); // x-z plane 
    } else {
        closest_vec = Vec3(0, -n.z, n.y); // y-z plane
    }
    rotX = closest_vec.norm();
    rotY = n.cross(rotX);
}

void trace(Ray &ray, Vec3 &color, std::shared_ptr<SceneObj> &obj, Intersection &hit, int depth) {

    if (depth >= bounces) {
        color = black;
        return;
    }

    bool intersect = get_intersection(ray, obj, hit);
    if (!intersect) {
        color = gray;
        return;
    }

    int material = obj->material;
    Vec3 emission = obj->emission;


    if (material == 1) {
        Vec3 direct_lighting = Vec3(0, 0, 0);

        for (unsigned long int i=0; i < scene_lights.size(); i++) {
            auto light = scene_lights.at(i);
            Vec3 light_direction = light->direction(hit);
            double intensity = light->intensity(hit);

            // Shadow test
            Ray shadow_ray = Ray(hit.point, light_direction);
            int visibility = calc_shadows(shadow_ray, light, hit);

            // Lambert's law
            double cos_theta = light_direction.dot(hit.normal);
            if (cos_theta < 0) continue;

            Vec3 direct_color = light->color * (visibility * cos_theta * intensity);
            direct_lighting = direct_lighting + direct_color;
        }

        Vec3 indirect_lighting = Vec3(0, 0, 0);
        
        if (global_illumination) {
            Vec3 rotX, rotY;
            create_orthonormal_system(hit.normal, rotX, rotY);
            float pdf = 1 / (2 * M_PI);

            // Trace indirect rays and accumulate results (Monte Carlo integration)
            for (int i=0; i < secondary_rays; i++) {
                float r1 = random_number();
                float r2 = random_number();
                Vec3 random_direction = sample_hemisphere(r1, r2);
                // Transform sample to local coordinate system
                float sample_x = Vec3(rotY.x, hit.normal.x, rotX.x).dot(random_direction);
                float sample_y = Vec3(rotY.y, hit.normal.y, rotX.y).dot(random_direction);
                float sample_z = Vec3(rotY.z, hit.normal.z, rotX.z).dot(random_direction);
                Vec3 sample_direction = Vec3(sample_x, sample_y, sample_z);
                Ray diffuse_ray = Ray(hit.point, sample_direction);

                float cos_theta = sample_direction.dot(hit.normal);
                Vec3 trace_color;
                std::shared_ptr<SceneObj> trace_obj;
                Intersection trace_hit;
                // trace(diffuse_ray, trace_color, depth+1);
                trace(diffuse_ray, trace_color, trace_obj, trace_hit, depth+1);
                trace_color = trace_color * cos_theta;
                indirect_lighting = indirect_lighting + trace_color;
            }

            // Average by number of samples and divide by pdf of random variable (same for all rays)
            indirect_lighting = indirect_lighting / (secondary_rays * pdf);
        } 

        // Apply the rendering equation 
        Vec3 brdf = obj->color; 
        Vec3 point_color = (direct_lighting + indirect_lighting).multVbyV(brdf);
        color = emission + point_color;
    }
}

void render(int width, int height) {
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            Vec3 pixel = Vec3(0, 0, 0);
            Ray ray = get_ray(x, y);
            Vec3 color;
            std::shared_ptr<SceneObj> draw_object;
            Intersection hit;
            trace(ray, color, draw_object, hit, 0);
            pixel = pixel + color;
            if (hit.distance > 0) {
                if (draw_object->type == 1) {
                    write_color(img, x, y, color, 255);
                } else if (draw_object->type == 2) {
                    auto paint_particle = std::make_shared<PaintParticle>(x, y, color, ray.direction, hit.normal, hit.distance, 1, 1); // none are edges for now
                    particles.push_back(paint_particle);
                }   
            } else {
                write_color(img, x, y, color, 255);   
            }
        }
    }
}

int main(int argc, char** argv) {
    // Read file
    std::ifstream stream(argv[1]);
    std::string line;

    // Process each line
    while (std::getline(stream, line)) {
        std::vector<std::string> tokens;
        // Skip blank lines (fix later?)
        if (line.size() == 1) {
            continue;
        }
        std::string token;
        std::istringstream s(line);
        while (s >> token) {
            tokens.push_back(token);
        }
        process_tokens(tokens);
    }

    // Render image
    render(width, height);

    // Go to painterly pipeline
    if (paint) {
        auto painter = Painter(Expressionist(), particles, &img);
        painter.paint();
    }

    // Save the image
    img.save_png(filename.c_str());

    return 0;
}
