// octreenode.cpp

#include "scene.h"

Scene::Scene(char** argv) {
	generator.seed(time(NULL));

	width = 0;
	height = 0;
	bounces = 4;
	environmentColor[0] = 0;
	environmentColor[1] = 0;
	environmentColor[2] = 0;
	environmentColor[3] = 0;
	stencil_ring_samples = 0;
	stencil_rings = 0.0;
	outline_cutoff = 1.0;
	shoot = 1;
	spp = 1;
	output_depth = 0;
	octree_box_limit_x8 = 2;
	max_octree_depth = 5;
	e_slope = 1.0;
	stencil_radius = 0.0;
	light_samples = 1;
	indirect_samples = 0;
	indirect_bounces = 0;
	flat = 0;


	vec3 _eye(0,0,0);
	eye = _eye;
	vec3 _forward(0,0,-1);
	forward = _forward;
	vec3 _right(1,0,0);
	right = _right;
	vec3 _up(0,1,0);
	up = _up;
	
	

	double x, y, z, r = 0;
	double A, B, C, D = 0;
	double u, v = 0;
	int objectID = 0;
	int object_type = 0;
	double ior = 1.458;
	double roughness = 0;
	double eccentricity = 0.0;
	int v1, v2, v3, v4, v1t, v2t, v3t, v4t, v1n, v2n, v3n, v4n;
	vertex *vert1, *vert2, *vert3, *vert4;
	double lastNormal[3] = {0,0,0};
	double lastUV[3] = {0,0,0};
	double lastColor[3] = {1, 1, 1};
	double lastEmission[3] = {0, 0, 0};
	vec3 shininess(0,0,0);
	vec3 transparency(0,0,0);
	vec3 *vert1t;
	vec3 *vert2t;
	vec3 *vert3t;
	vec3 *vert1n;
	vec3 *vert2n;
	vec3 *vert3n;
	std::string v1parts = "";
	std::string v2parts = "";
	std::string v3parts = "";
	std::string materialFile = "";
	std::string currMaterial = "";
	std::string lineBuffer = "";
	bool material_file_provided = false;
	bool color_provided = false;
	std::string token = "";
	outputFileName = "";
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
	dimension = width * height;
	ss >> outputFileName;
	
	std::vector<object*> objects;
	std::vector<vertex*> verts;
	std::vector<vec3*> normals;

	std::vector<vec3*> texcoords;
	std::map<std::string, material *> materials;
	std::vector<std::string> texture_files;
	std::map<std::string, image_texture *> texture_maps;

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
			color_provided = true;
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
			vec3 uv(lastUV[0], lastUV[0], lastUV[0]);
			vec3 normal(lastNormal[0], lastNormal[1], lastNormal[2]);
			normal.normalize();
			vec3 color(lastColor[0], lastColor[1], lastColor[2]);
			vertex *v = new vertex(coordinate, uv, normal, color);
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
		else if(!token.compare("mtllib")) {
			ss >> materialFile;
			std::string input_directory = "input_files/";
			input_directory.append(materialFile);
			material_file_provided = true;
			readMaterialFile(input_directory, materials, texture_files);
			initTextureMaps(texture_files, texture_maps);
		}
		else if (!token.compare("v")) {
			ss >> x;
			ss >> y;
			ss >> z;
			vec3 coordinate(x,y,z);
			vec3 uv(lastUV[0], lastUV[0], lastUV[0]);
			vec3 normal(lastNormal[0], lastNormal[1], lastNormal[2]);
			normal.normalize();
			vec3 color(lastColor[0], lastColor[1], lastColor[2]);
			vertex *v = new vertex(coordinate, uv, normal, color);
			verts.push_back(v);
		}
		else if (!token.compare("vt")) {
			ss >> u;
			ss >> v;
			vec3 *t = new vec3(u, v, 0);
			texcoords.push_back(t);
		}
		else if (!token.compare("vn")) {
			ss >> x;
			ss >> y;
			ss >> z;
			vec3 *v = new vec3(x, y, z);
			normals.push_back(v);
		}
		else if (!token.compare("usemtl")) {
			ss >> currMaterial;
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
			v1t = std::stoi(v1parts.substr(v1tStart, v1tSlash - v1tStart), nullptr);
			v1n = std::stoi(v1parts.substr(v1nStart), nullptr);

			int v2Slash = v2parts.find("/", 0);
			int v2tStart = v2Slash + 1;
			int v2tSlash = v2parts.find("/", v2tStart);
			int v2nStart = v2tSlash + 1;
			v2 = std::stoi(v2parts.substr(0, v2Slash), nullptr);
			v2t = std::stoi(v2parts.substr(v2tStart, v2tSlash - v2tStart), nullptr);
			v2n = std::stoi(v2parts.substr(v2nStart), nullptr);

			int v3Slash = v3parts.find("/", 0);
			int v3tStart = v3Slash + 1;
			int v3tSlash = v3parts.find("/", v3tStart);
			int v3nStart = v3tSlash + 1;
			v3 = std::stoi(v3parts.substr(0, v3Slash), nullptr);
			v3t = std::stoi(v3parts.substr(v3tStart, v3tSlash - v3tStart), nullptr);
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
			if (v1t < 0) {
				vert1t = texcoords.at(texcoords.size() + v1t);
			}
			else {
				vert1t = texcoords.at(v1t - 1);
			}
			if (v2t < 0) {
				vert2t = texcoords.at(texcoords.size() + v2t);
			}
			else {
				vert2t = texcoords.at(v2t - 1);
			}
			if (v3t < 0) {
				vert3t = texcoords.at(texcoords.size() + v3t);
			}
			else {
				vert3t = texcoords.at(v3t - 1);
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
			vert1 = new vertex(vert1->xyz, *vert1t, *vert1n, vert1->color);
			vert2 = new vertex(vert2->xyz, *vert2t, *vert2n, vert2->color);
			vert3 = new vertex(vert3->xyz, *vert3t, *vert3n, vert3->color);

			material *m;
			texture *tex;
			if (material_file_provided) {
				m = materials.at(currMaterial);
				if (m->diffuse_map != "") {
					tex = texture_maps.at(m->diffuse_map);
				}
				else {
					tex = new color_texture(m->diffuse_color);
				}
			}
			else if (color_provided) {
				m = new material(shininess,transparency,ior,roughness,eccentricity);
				vec3 last_color(lastColor[0], lastColor[1], lastColor[2]);
				tex = new color_texture(last_color);
			}
			tri *t = new tri(vert1, vert2, vert3, tex, lastEmission, m, objectID, object_type);
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
		else if (!token.compare("environment")) {
			ss >> environmentColor[0];
			ss >> environmentColor[1];
			ss >> environmentColor[2];
			ss >> environmentColor[3];
			environment = 1;
		}
		else if (!token.compare("octree_box_limit_x8")) {
			ss >> octree_box_limit_x8;
		}
		else if (!token.compare("max_octree_depth")) {
			ss >> max_octree_depth;
		}
		else {
				//skip line
		}
	}

	octree = new OctreeNode(octree_box_limit_x8, 0, -1, objects.size(), max_octree_depth);
	octree->buildOctree(objects);

}

void Scene::readMaterialFile(std::string material_file, std::map<std::string, material *> &materials, std::vector<std::string> &texture_files) {
	
	std::ifstream stream(material_file); 
    	std::string line;

	std::string curr_mat_name;
	material *curr_mat;

	while (std::getline(stream, line)) {
		if (line.empty()) {
			continue;
		}
		std::string token;
		std::stringstream ss(line);
		ss >> token;
		if (!token.compare("newmtl")) {
			// material name
			ss >> curr_mat_name;
			material *new_mat = new material();
			materials.insert(std::pair<std::string, material *>(curr_mat_name, new_mat)); // maybe check if key is already in there
			curr_mat = materials.at(curr_mat_name);
		}
		else if (!token.compare("illum")) {
			// illumination model
		}
		else if (!token.compare("Ka")) {
			// ambient color
			ss >> curr_mat->ambient.x;
			ss >> curr_mat->ambient.y;
			ss >> curr_mat->ambient.z;
		}
		else if (!token.compare("Kd")) {
			// diffuse color
			ss >> curr_mat->diffuse_color.x;
			ss >> curr_mat->diffuse_color.y;
			ss >> curr_mat->diffuse_color.z;
		}
		else if (!token.compare("Ks")) {
			// specular color
			//ss >> curr_mat->shininess.x;
			//ss >> curr_mat->shininess.y;
			//ss >> curr_mat->shininess.z;
		}
		else if (!token.compare("Ns")) {
			// eccentricity
			ss >> curr_mat->eccentricity;
		}
		else if (!token.compare("Tf")) {
			// transmission color
			// ss >> curr_mat->transparency.x;
			// ss >> curr_mat->transparency.y;
			// ss >> curr_mat->transparency.z;
		}
		else if (!token.compare("map_Kd")) {
			// diffuse texture map -> image name
			ss >> curr_mat->diffuse_map;
			std::string input_directory = "input_files/";
    			input_directory.append(curr_mat->diffuse_map);
			curr_mat->diffuse_map = input_directory;
			texture_files.push_back(curr_mat->diffuse_map);
		}
		else if (!token.compare("Ni")) {
			// index of refraction
			ss >> curr_mat->ior;
		}
		else {
			// ignore line 
		}
	}
}

void Scene::initTextureMaps(std::vector<std::string> &texture_files, std::map<std::string, image_texture *> &texture_maps) {
	for (int i=0; i < texture_files.size(); i++) {
		image_texture *new_texture = new image_texture(texture_files.at(i));
		texture_maps.insert(std::pair<std::string, image_texture *>(texture_files.at(i), new_texture));
	}
}

void Scene::PrintScene() {
    
}

Scene::~Scene() {

}
