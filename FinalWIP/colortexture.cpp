#include "colortexture.h"


vec3 color_texture::getColor(double u, double v, vec3 &hit_point) {
	return this->color;
}

color_texture::~color_texture() {

}
