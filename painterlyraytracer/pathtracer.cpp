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

// PAINTING

int get_mask_value(std::tuple<int, int> offset, int radius) {
        int x = (0 - std::get<0>(offset));
        int x2 = x * x;
        int y = (0 - std::get<1>(offset));
        int y2 = y * y;
        if (x2 + y2 <= radius * radius) return 1;
        return 0;
    }

void create_mask(int **brush, const int radius) {
    // brush.resize(radius, std::vector<int>(radius));
    for (int i=-radius; i < radius+1; i++) {
        for (int j=-radius; j < radius+1; j++) {
            auto offset = std::make_tuple (i, j);
            int mask_val = get_mask_value(offset, radius);
            // printf("(%d, %d): %d \n", i, j, mask_val);
            // brush.at(i+radius).at(j+radius) = mask_val; // fix this 
            brush[i+radius][j+radius] = mask_val;
        }
    }
}

std::vector<StrokeParticle> makeSplineStroke(Style style, std::shared_ptr<PaintParticle> &particle, int min_x, int max_x, int min_y, int max_y) {
    float x0 = particle->x;
    float y0 = particle->y;
    Vec3 color = particle->color;
    int size = particle->size;
    
    std::vector<StrokeParticle> stroke;
    StrokeParticle s0 = StrokeParticle(x0, y0, color, size);
    stroke.push_back(s0);

    if (particle->edge == 2) {
        return stroke;
    }

    float x = x0;
    float y = y0;
    float lastDx = 0;
    float lastDy = 0;

    for (int i=1; i < style.max_stroke_length; i++) {
        // Find tangent stroke
        // Vec3 up = Vec3(1, 0, 0);
        // float dot_nd = particle->normal.dot(up);
        // float mag_d = up.mag() * up.mag();
        // Vec3 o_prime = up * (dot_nd/mag_d);

        // float dot_nd = particle->normal.dot(particle->direction);
        // float mag_d = particle->direction.mag() * particle->direction.mag();
        // Vec3 o_prime = particle->direction * (dot_nd/mag_d);

        // Vec3 o_prime = particle->normal.cross(particle->direction); // orientation for edge particles
        Vec3 o_prime = particle->normal; 

        // Get unit vector of gradient
        float gx = o_prime.x;
        float gy = o_prime.y;
        // std::cout << "gx: " << gx << " gy: " << gy << std::endl;
        // if (gx == 0 && gy == 0) continue;
        float dx = -gy;
        float dy = -gx;
        // float dx = -1;
        // float dy = 1;

        if ((lastDx * dx + lastDy * dy) < 0) {
            dx = -dx;
            dy = -dy;
        }
        // std::cout << "bef dx: " << dx << " dy: " << dy << std::endl;
        dx = style.curvature_filter * dx + (1 - style.curvature_filter) * lastDx;
        dy = style.curvature_filter * dy + (1 - style.curvature_filter) * lastDy;
        // std::cout << "bet dx: " << dx << " dy: " << dy << std::endl;
        dx = dx / sqrt(dx * dx + dy * dy);
        dy = dy / sqrt(dx * dx + dy * dy);
        // std::cout << "after dx: " << dx << " dy: " << dy << std::endl;

        // float yDistort = Vec3(dx, dy, 0).mag() * style.curvature_filter * (x - 0.5) * (x - 0.5);

        // Filter stroke direction
        x = x + size * dx;
        y = y + size * dy;
        lastDx = dx;
        lastDy = dy;

        if (x < min_x || x >= max_x) continue;
        if (y < min_y || y >= max_y) continue;

        StrokeParticle s = StrokeParticle(x, y, color, size);
        stroke.push_back(s);
    }
    return stroke;
}

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

void render(Scene &scene, Style style, int radius, int **brush) {
    // start randomize pixels
    std::vector<int> x_values;
    std::vector<int> y_values;

    for (int i=0; i < scene.width; i++) {
        x_values.push_back(i);
    }
    for (int j=0; j < scene.height; j++) {
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
    // end randomize pixels

    // for (int x = 0; x < scene.width; x++) {
    //     for (int y = 0; y < scene.height; y++) {
    // for (int x : x_values) {
    //     for (int y : y_values) {
    for (std::tuple<int, int> curr_pixel : pixels) {
        int x = std::get<0>(curr_pixel);
        int y = std::get<1>(curr_pixel);
        Vec3 pixel = Vec3(0, 0, 0);
        Ray ray = scene.cam.get_ray(x, y, scene.width, scene.height);
        Vec3 color;
        std::shared_ptr<SceneObj> draw_object;
        Intersection hit;
        trace(ray, color, draw_object, hit, scene, 0);
        pixel = pixel + color;
        if (hit.distance) {
            if (draw_object->type == 1) {
                if (scene.img(x, y, 0, 3) == 0) write_color(scene.img, x, y, color, 255); // do if not already paint there
            } else if (draw_object->type == 2) {
                auto paint_particle = std::make_shared<PaintParticle>(x, y, color, ray.direction, hit.normal, hit.distance, 1, 1); // none are edges for now
                scene.particles.push_back(paint_particle);
                auto s = makeSplineStroke(style, paint_particle, 0, scene.width, 0, scene.height); // don't include stroke points that are behind objects
                for (unsigned long int m=0; m < s.size(); m++) {
                    auto p0 = s.at(m);  // get current point in stroke
                    Vec3 stroke_color = p0.color; 

                    int x = p0.x;
                    int y = p0.y;
                    // Paint surrounding pixels (maybe make function for this)
                    for (int j=-radius; j < radius+1; j++) {
                        for (int k=-radius; k < radius+1; k++) {
                            int curr_x = x + j;
                            int curr_y = y + k;
                            if (curr_x < 0 || curr_x >= scene.width) continue;
                            if (curr_y < 0 || curr_y >= scene.height) continue;
                            int paint_num = brush[j+radius][k+radius];
                            if (paint_num) {
                                write_color(scene.img, curr_x, curr_y, stroke_color, 255); // pick from shuffled list of x and y values
                            }
                        }
                    }
                }
            } 
        } else {
            // color.print();
            if (scene.img(x, y, 0, 3) == 0) write_color(scene.img, x, y, color, 255);  // do this only if there isn't color already there
        }
    }
    //     }
    // }   
}

int main(int argc, char** argv) {
    // Read file
    auto reader = SceneReader();
    reader.process_scene_file(argv[1]);

    // Newly initialized scene
    auto scene = reader.scene;

    // Initialize brush
    auto style = Expressionist();
    int mask_size = (style.brush_size * 2) + 1;
    int **brush = new int*[mask_size];
    for (int i=0; i < mask_size; i++) {
        brush[i] = new int[mask_size];
    }
    create_mask(brush, style.brush_size);

    // Render scene
    render(scene, style, style.brush_size, brush);

    // Go to painterly pipeline
    // if (scene.paint) {
    //     auto painter = Painter(Expressionist(), scene.particles, &scene.img);
    //     painter.paint();
    // }

    // Save the image
    scene.img.save_png(scene.filename.c_str());

    // Deallocate memory
    for (int i=0; i < mask_size; i++) {
        delete [] brush[i];
    }
    delete [] brush;

    return 0;
}
