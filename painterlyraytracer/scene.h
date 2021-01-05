#ifndef SCENE_H
#define SCENE_H

#include <sstream>
#include <fstream>
#include <string>

#include "utilities.h"
#include "objects.h"
#include "lights.h"
#include "particle.h"

class Camera {
public:
    Vec3 eye, forward, right, up;
    Camera() { 
        eye = Vec3(0, 0, 0); 
        forward = Vec3(0, 0, -1); 
        right = Vec3(1, 0, 0); 
        up = Vec3(0, 1, 0); 
    }

    Ray get_ray(int x, int y, int width, int height) {
        Ray ray;

        float denominator = std::max(width, height);
        float sx = (2 * x - width) / denominator;
        float sy = (height - 2 * y) / denominator;

        Vec3 sr_su = (right * sx) + (up * sy);
        ray.direction = (forward + sr_su).norm();
        ray.origin = eye;

        return ray;
    }
    void set_eye(const Intersection &intersection) { return; }
    void set_forward(const Intersection &intersection) { return; }
    void set_up(const Intersection &intersection) { return; }
};

class Scene {
public:
    int width, height, depth, spectrum; 
    std::string filename;
    CImg<float> img;
    Vec3 black, gray, white;
    float gamma_val;
    bool global_illumination;
    int secondary_rays, bounces;
    bool background, paint;
    std::vector<Vec3> color_list;
    std::vector<Vec3> emission_values;
    std::vector<std::shared_ptr<Light>> scene_lights; 
    std::vector<std::shared_ptr<SceneObj>> scene_objects;
    std::vector<std::shared_ptr<PaintParticle>> particles;   
    Camera cam;
    Scene() { 
        width = 0; 
        height = 0; 
        depth = 1;
        spectrum = 4;
        black = Vec3(0, 0, 0); 
        gray = Vec3(0.5, 0.5, 0.5); 
        white = Vec3(1, 1, 1); 
        gamma_val = 1.0;
        global_illumination = false;
        secondary_rays = 1;
        bounces = 2;
        background = false;
        paint = false;
        cam = Camera();
    }
};

class SceneReader {
public:
    Scene scene;
    SceneReader() { scene = Scene(); }
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
        if (scene.background) {
            new_sphere->type = 2; 
        }
        new_sphere->color = scene.color_list.back();
        new_sphere->emission = scene.emission_values.back();
        scene.scene_objects.push_back(new_sphere);
    }

    void color(std::vector<std::string> &tokens) {
        float r = std::stof(tokens.at(1));
        float g = std::stof(tokens.at(2));
        float b = std::stof(tokens.at(3));
        Vec3 new_color = Vec3(r, g, b);
        
        scene.color_list.push_back(new_color);
    }

    void emission(std::vector<std::string> &tokens) {
        float e1 = std::stof(tokens.at(1));
        float e2 = std::stof(tokens.at(2));
        float e3 = std::stof(tokens.at(3));
        Vec3 new_emission = Vec3(e1, e2, e3);
        
        scene.emission_values.push_back(new_emission);
    }

    void sun(std::vector<std::string> &tokens) {
        float x = std::stof(tokens.at(1));
        float y = std::stof(tokens.at(2));
        float z = std::stof(tokens.at(3));
        Vec3 direction = Vec3(x, y, z);
        Vec3 color = scene.color_list.back();

        // Create new light and add to scene lights
        auto new_light = std::make_shared<Sun>(direction, color);
        scene.scene_lights.push_back(new_light);
    }

    void bulb(std::vector<std::string> &tokens) {
        float x = std::stof(tokens.at(1));
        float y = std::stof(tokens.at(2));
        float z = std::stof(tokens.at(3));
        Vec3 position = Vec3(x, y, z);
        Vec3 color = scene.color_list.back();

        // Create new light and add to scene lights
        auto new_light = std::make_shared<Bulb>(position, color);
        scene.scene_lights.push_back(new_light);
    }

    void gi(std::vector<std::string> &tokens) {
        scene.global_illumination = true;
    }

    void spp(std::vector<std::string> &tokens) {
        scene.secondary_rays = std::stoi(tokens.at(1));
    }

    void add_to_background(std::vector<std::string> &tokens) {
        scene.background = true;
    }

    void add_to_foreground(std::vector<std::string> &tokens) {
        scene.background = false;
    }

    void paintstyle(std::vector<std::string> &tokens) {
        scene.paint = true;
    }

    void create_png(std::vector<std::string> &tokens) {
        scene.width = std::stoi(tokens.at(1));
        scene.height = std::stoi(tokens.at(2));
        scene.filename = tokens.at(3);

        // Create image where each pixel is a float and set all pixels to 0
        scene.img.assign(scene.width, scene.height, scene.depth, scene.spectrum, 0); // make member function of struct
        // Set default color and emission values
        scene.color_list.push_back(Vec3(1, 1, 1));
        scene.emission_values.push_back(Vec3(0, 0, 0));
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
        } else if (keyword == "background") {
            add_to_background(tokens);
        } else if (keyword == "foreground") {
            add_to_foreground(tokens);
        } else if (keyword == "paint") {
            paintstyle(tokens);
        }
        return 0;
    }

    void process_scene_file(std::string scene_file) {
        std::ifstream stream(scene_file); 
        std::string line;

        // Process each line
        while (std::getline(stream, line)) {
            std::vector<std::string> tokens;
            // Skip blank lines 
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
    }
};

#endif