#include "imagetexture.h"

image_texture::image_texture(std::string file) {
    this->file = file;

    CImg<float> image(this->file.c_str());
    this->image = image;

    this->width = this->image.width();
    this->height = this->image.height();
}

vec3 image_texture::getColor(double u, double v, vec3 &hit_point) {
    // image texture mapping -> https://raytracing.github.io/books/RayTracingTheNextWeek.html#imagetexturemapping 
	vec3 texture_color(0,0,0);
    int width = this->width;
    int height = this->height;

    u = clamp(u, 0.0, 1.0);
    v = clamp(v, 0.0, 1.0);
    v = 1 - v;

    int x = (int)(u * width);
    int y = (int)(v * height);

    if (x >= width) x = width - 1;
    if (y >= height) y = height - 1; 

    texture_color.x = (this->image)(x, y, 0, 0) / 255;
    texture_color.y = (this->image)(x, y, 0, 1) / 255;
    texture_color.z = (this->image)(x, y, 0, 2) / 255;

	return texture_color;
}