#include <algorithm>
#include <iostream>
#include <thrust/detail/config/host_device.h>

#include "Grid.h"
#include "utils.h"
Grid::Grid() :cellSize() {

}

Grid::~Grid() {

}

void Grid::init(const PCISPH::Vec3 boxSize, float cellSize) {
	this->cellSize = cellSize;
	this->boxSize = boxSize;

	this->gridSize = PCISPH::iVec3(
		(int)floor(boxSize.x / cellSize) + 1,
		(int)floor(boxSize.y / cellSize) + 1,
		(int)floor(boxSize.z / cellSize) + 1
	);
	cellNumber = gridSize.x * gridSize.y * gridSize.z;

}

void Grid::update(const thrust::host_vector<PCISPH::Vec3>& positions, thrust::host_vector<size_t>& offset, std::function<void(size_t, size_t)> swap){
	// particle count in each grid cell
	std::vector<size_t> cellCount(cellNumber, 0);
	// pointer to particle array index
	std::vector<size_t> cellIndex(cellNumber, 0);

	size_t particleNumber = positions.size();
	// grid indices for each particle
	std::vector<size_t> indices(particleNumber);

#pragma omp parallel for
	for (size_t i = 0; i < particleNumber; i++) {
		size_t index = linearIndex(positions[i]);
		indices[i] = index;
		cellCount[index] += 1;
	}
	offset.resize(cellNumber);

	size_t index = 0;
#pragma omp parallel for
	for (size_t i = 0; i < cellNumber; i++) {
		offset[i] = index;
		cellIndex[i] = index;
		index += cellCount[i];
	}

	auto a = offset.data();
	
	offset.back() = index;

#pragma omp parallel for
	for (size_t i = 0; i < particleNumber; i++) {
		while (i < offset[indices[i]] || i >= cellIndex[indices[i]]) {
			size_t j = cellIndex[indices[i]]++;
			std::swap(indices[i], indices[j]);
			swap(i, j);
		}
	}
}

PCISPH::uVec3 Grid::getGridPos(const PCISPH::Vec3 pos) const {
	PCISPH::Vec3 p = pos;

	if (p.x < 0) p.x = 0;
	if (p.y < 0) p.y = 0;
	if (p.z < 0) p.z = 0;

	if (p.x > boxSize.x) p.x = boxSize.x;
	if (p.y > boxSize.y) p.y = boxSize.y;
	if (p.z > boxSize.z) p.z = boxSize.z;


	return PCISPH::uVec3(
		(PCISPH::uint)floor(p.x / cellSize),
		(PCISPH::uint)floor(p.y / cellSize),
		(PCISPH::uint)floor(p.z / cellSize)
	);
}

__device__ __host__ 
PCISPH::uVec3 getGridPos(const PCISPH::Vec3 &pos, PCISPH::Vec3 boxSize, float cellSize) {
	PCISPH::Vec3 p = pos;

	if (p.x < 0) p.x = 0;
	if (p.y < 0) p.y = 0;
	if (p.z < 0) p.z = 0;

	if (p.x > boxSize.x) p.x = boxSize.x;
	if (p.y > boxSize.y) p.y = boxSize.y;
	if (p.z > boxSize.z) p.z = boxSize.z;

	return PCISPH::uVec3(
		(PCISPH::uint)floor(p.x / cellSize),
		(PCISPH::uint)floor(p.y / cellSize),
		(PCISPH::uint)floor(p.z / cellSize)
	);
}


__device__ __host__ 
size_t Grid::linearIndex(const PCISPH::uVec3 gridPos) const {
	return gridPos.x + gridPos.y * gridSize.x + gridPos.z * gridSize.x * gridSize.y;
}

__device__ __host__ 
size_t Grid::linearIndex(const PCISPH::Vec3 pos) const {
	return linearIndex(getGridPos(pos));
}

/* Query the neighbor points
* Then do func(i, j)
* Expect Func format:
* void Func(size_t j);
*/
//template<typename Func>
