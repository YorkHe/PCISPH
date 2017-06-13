#include "Grid.h"
#include "utils.h"
#include <algorithm>
#include <iostream>

Grid::Grid() :cellSize(), offset() {

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

	offset.resize(cellNumber + 1);
}

void Grid::update(const std::vector<PCISPH::Vec3> &positions, std::function<void(size_t, size_t)> swap) {
	// particle count in each grid cell
	std::vector<size_t> cellCount(cellNumber, 0);
	// pointer to particle array index
	std::vector<size_t> cellIndex(cellNumber, 0);

	size_t particleNumber = positions.size();
	// grid indices for each particle
	std::vector<size_t> indices(particleNumber);

	for (size_t i = 0; i < particleNumber; i++) {
		size_t index = linearIndex(positions[i]);
		indices[i] = index;
		cellCount[index] += 1;
	}

	size_t index = 0;
	for (size_t i = 0; i < cellNumber; i++) {
		offset[i] = index;
		cellIndex[i] = index;
		index += cellCount[i];
	}
	offset.back() = index;

	for (size_t i = 0; i < particleNumber; i++) {
		while (i < offset[indices[i]] || i >= cellIndex[indices[i]]) {
			size_t j = cellIndex[indices[i]]++;
			std::swap(indices[i], indices[j]);
			swap(i, j);
		}
	}
}

/* Query the neighbor points
* Then do func(i, j)
* Expect Func format:
* void Func(size_t j);
*/
//template<typename Func>
void Grid::query(const PCISPH::Vec3 &pos, std::function<void(size_t)> func) const {
	glm::u64vec3 boundBoxMin = this->getGridPos(pos - PCISPH::Vec3(this->cellSize));
	glm::u64vec3 boundBoxMax = this->getGridPos(pos + PCISPH::Vec3(this->cellSize));

	for (size_t z = boundBoxMin.z; z <= boundBoxMax.z; z++) {
		for (size_t y = boundBoxMin.y; y <= boundBoxMax.y; y++) {
			for (size_t x = boundBoxMin.x; x <= boundBoxMax.x; x++) {
				size_t cellIndex = linearIndex(PCISPH::uVec3(x, y, z));
				for (size_t neighborIndex = offset[cellIndex]; neighborIndex < offset[cellIndex + 1]; neighborIndex++) {
					func(neighborIndex);
				}
			}
		}
	}
}


inline PCISPH::uVec3 Grid::getGridPos(const PCISPH::Vec3 &pos) const {
	PCISPH::Vec3 p(pos);
	if (PCISPH::minComponent(p) < 0) {
		p = PCISPH::Vec3(0);
	}
	if (PCISPH::maxComponent(p - boxSize) > 0) {
		p = boxSize;
	}
	return PCISPH::uVec3(
		(PCISPH::uint)floor(p.x / cellSize),
		(PCISPH::uint)floor(p.y / cellSize),
		(PCISPH::uint)floor(p.z / cellSize)
	);
}

inline size_t Grid::linearIndex(const PCISPH::uVec3 &gridPos) const {
	return gridPos.x + gridPos.y * gridSize.x + gridPos.z * gridSize.x * gridSize.y;
}

inline size_t Grid::linearIndex(const PCISPH::Vec3 &pos) const {
	return linearIndex(getGridPos(pos));
}