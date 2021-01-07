#ifndef PAINTER_H
#define PAINTER_H

#include <algorithm>
#include "utilities.h"
#include "style.h"
#include "particle.h"

class Painter {
public:
    Style style;
    std::vector<std::shared_ptr<PaintParticle>> particles;
    CImg<float> *img;
    int mask_size;
    int **brush;
    Painter(Style s, std::vector<std::shared_ptr<PaintParticle>> p, CImg<float> *i) {
        style = s;
        particles = p;
        img = i;
        mask_size = (style.brush_size * 2) + 1;

        brush = new int*[mask_size];
        for (int i=0; i < mask_size; i++) {
            brush[i] = new int[mask_size];
        }
    }

    ~Painter() {
        for (int i=0; i < mask_size; i++) {
            delete [] brush[i];
        }
        delete [] brush;
    }

    int get_mask_value(std::tuple<int, int> offset, int radius) {
        int x = (0 - std::get<0>(offset));
        int x2 = x * x;
        int y = (0 - std::get<1>(offset));
        int y2 = y * y;
        if (x2 + y2 <= radius * radius) return 1;
        return 0;
    }

    void create_mask(const int radius) {
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

    std::vector<StrokeParticle> makeSplineStroke(std::shared_ptr<PaintParticle> &particle, int min_x, int max_x, int min_y, int max_y) {
        float x0 = particle->x;
        float y0 = particle->y;
        vec3 color = particle->color;
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
            vec3 o_prime = particle->normal; 

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

    void paint_layer(int radius) {
        std::shuffle(particles.begin(), particles.end(), std::random_device()); // add to utilities

        std::vector<vec3> edge_particles;
        // for (PaintParticle p : particles) {
        //     if (p.edge == 2) edge_particles.push_back(p.position);
        // }

        int min_x, max_x, min_y, max_y;
        if (edge_particles.size() == 0) {
            min_x = 0;
            max_x = img->width();
            min_y = 0;
            max_y = img->height();
        } else {
            // TODO: find best translation for this part
            std::cout << "still implementing" << std::endl;
        }

        // Introduce more randomness
        std::vector<std::vector<StrokeParticle>> strokes;
        for (unsigned long int i=0; i < particles.size(); i++) {
            auto s = makeSplineStroke(particles.at(i), min_x, max_x, min_y, max_y);
            strokes.push_back(s);
        }
        // std::shuffle(strokes.begin(), strokes.end(), std::random_device()); 

        for (unsigned long int i=0; i < strokes.size(); i++) {
            for (unsigned long int m=0; m < strokes.at(i).size(); m++) {
                auto p0 = strokes.at(i).at(m);  // get current point in stroke
                vec3 stroke_color = p0.color; 

                int x = p0.x;
                int y = p0.y;
                // Paint surrounding pixels (maybe make function for this)
                for (int j=-radius; j < radius+1; j++) {
                    for (int k=-radius; k < radius+1; k++) {
                        int curr_x = x + j;
                        int curr_y = y + k;
                        if (curr_x < 0 || curr_x >= img->width()) continue;
                        if (curr_y < 0 || curr_y >= img->height()) continue;
                        int paint_num = brush[j+radius][k+radius];
                        if (paint_num) {
                            write_color(*img, curr_x, curr_y, stroke_color, 255);
                        }
                    }
                }
            }
        }
    }

    void paint() {
        create_mask(style.brush_size);
        paint_layer(style.brush_size);
    }
};

#endif
