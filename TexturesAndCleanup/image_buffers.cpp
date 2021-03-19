// octreenode.cpp

#include "image_buffers.h"

ImageBuffers::ImageBuffers(int dimension) {
	depthMap = new double[dimension];
	edgeMap = new double[dimension*3];
	objectTypeMap = new int[dimension];
	objectBoundaryMap = new int[dimension];
	objectIDMap = new int[dimension];
	shadowMap = new int[dimension];
	paintMap = new float[dimension];

	for (int pm = 0; pm < dimension; pm++) {
		paintMap[pm] = 0;
	}

}

ImageBuffers::~ImageBuffers() {

}
