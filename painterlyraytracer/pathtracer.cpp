// To compile this file, type the following:
//   g++ pathtracer.cpp -lpng -lX11 -lpthread -o output
// To run, type the following:
//   ./output files/test_ground_objects_indirect.txt

#define cimg_use_png

#include "utilities.h"
#include "scene.h"
#include "objects.h"
#include "lights.h"
#include "painter.h"

// RAYTRACING

bool get_intersection(Ray &ray, std::shared_ptr<SceneObj> &draw_object, Intersection &hit, Scene &scene) {
    float min_distance = std::numeric_limits<float>::max(); 

    bool intersect = false;

    for (unsigned long int i=0; i < scene.scene_objects.size(); i++) {
        Intersection intersection;
        bool found_intersection = scene.scene_objects.at(i)->intersect(ray, intersection);
        if (found_intersection) {  
            intersect = true;
            float t = intersection.distance;
            if (t < min_distance) {
                draw_object = scene.scene_objects.at(i);
                hit = intersection;
                min_distance = t;
            }
        }
    }
    return intersect;
}

int calc_shadows(Ray &shadow_ray, std::shared_ptr<Light> &light, Intersection &hit, Scene &scene) {
    double light_distance = light->distance(hit);
    for (unsigned long int i=0; i < scene.scene_objects.size(); i++) {
        auto shadow_object = scene.scene_objects.at(i);
        Intersection shadow_hit;
        bool found_intersection = scene.scene_objects.at(i)->intersect(shadow_ray, shadow_hit);
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

void trace(Ray &ray, Vec3 &color, std::shared_ptr<SceneObj> &obj, Intersection &hit, Scene &scene, int depth) {

    if (depth >= scene.bounces) {
        color = scene.black;
        return;
    }

    bool intersect = get_intersection(ray, obj, hit, scene);
    if (!intersect) {
        color = scene.gray;
        return;
    }

    int material = obj->material;
    Vec3 emission = obj->emission;


    if (material == 1) {
        Vec3 direct_lighting = Vec3(0, 0, 0);

        for (unsigned long int i=0; i < scene.scene_lights.size(); i++) {
            auto light = scene.scene_lights.at(i);
            Vec3 light_direction = light->direction(hit);
            double intensity = light->intensity(hit);

            // Shadow test
            Ray shadow_ray = Ray(hit.point, light_direction);
            int visibility = calc_shadows(shadow_ray, light, hit, scene);

            // Lambert's law
            double cos_theta = light_direction.dot(hit.normal);
            if (cos_theta < 0) continue;

            Vec3 direct_color = light->color * (visibility * cos_theta * intensity);
            direct_lighting = direct_lighting + direct_color;
        }

        Vec3 indirect_lighting = Vec3(0, 0, 0);
        
        if (scene.global_illumination) {
            Vec3 rotX, rotY;
            create_orthonormal_system(hit.normal, rotX, rotY);
            float pdf = 1 / (2 * M_PI);

            // Trace indirect rays and accumulate results (Monte Carlo integration)
            for (int i=0; i < scene.secondary_rays; i++) {
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
                trace(diffuse_ray, trace_color, trace_obj, trace_hit, scene, depth+1);
                trace_color = trace_color * cos_theta;
                indirect_lighting = indirect_lighting + trace_color;
            }

            // Average by number of samples and divide by pdf of random variable (same for all rays)
            indirect_lighting = indirect_lighting / (scene.secondary_rays * pdf);
        } 

        // Apply the rendering equation 
        Vec3 brdf = obj->color; 
        Vec3 point_color = (direct_lighting + indirect_lighting).multVbyV(brdf);
        color = emission + point_color;
    }
}

void render(Scene &scene) {
    for (int x = 0; x < scene.width; x++) {
        for (int y = 0; y < scene.height; y++) {
            Vec3 pixel = Vec3(0, 0, 0);
            Ray ray = scene.cam.get_ray(x, y, scene.width, scene.height);
            Vec3 color;
            std::shared_ptr<SceneObj> draw_object;
            Intersection hit;
            trace(ray, color, draw_object, hit, scene, 0);
            pixel = pixel + color;
            if (hit.distance > 0) {
                if (draw_object->type == 1) {
                    write_color(scene.img, x, y, color, 255);
                } else if (draw_object->type == 2) {
                    auto paint_particle = std::make_shared<PaintParticle>(x, y, color, ray.direction, hit.normal, hit.distance, 1, 1); // none are edges for now
                    scene.particles.push_back(paint_particle);
                }   
            } else {
                write_color(scene.img, x, y, color, 255);   
            }
        }
    }
}

int main(int argc, char** argv) {
    // Read file
    auto reader = SceneReader();
    reader.process_scene_file(argv[1]);

    // Newly initialized scene
    auto scene = reader.scene;

    // Render scene
    render(scene);

    // Go to painterly pipeline
    if (scene.paint) {
        auto painter = Painter(Expressionist(), scene.particles, &scene.img);
        painter.paint();
    }

    // Save the image
    scene.img.save_png(scene.filename.c_str());

    return 0;
}
